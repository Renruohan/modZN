"""
modZN common library
====================

The shared core of the two design tools ``modZN_linear.py`` and ``modZN_circular.py``.
It designs the overlapping zones between neighbouring fragments and the PCR primers
(phosphorothioate bonds are written as ``*``) of a modAssembly experiment.

Contents:

* thermodynamics (nearest-neighbour delta G, primer3 Tm / hairpin / dimer),
* the overlap search (``exact`` and ``greedy``) and the three overlap filters,
* primer design helpers, result containers, report and CSV output,
* the command-line and graphical (FreeSimpleGUI) interfaces.

Terms used throughout:

* fragment : one DNA piece to be assembled (a string over ``ACGT``).
* junction : the boundary between two neighbouring fragments. ``n`` fragments have
  ``n - 1`` junctions (linear) or ``n`` junctions (circular).
* overlap  : the overlapping zone of a junction: the last ``left`` bases of the upstream
  fragment plus the first ``overlap_length - left`` bases of the downstream fragment.
* layout   : the ``left`` value of every junction, as a tuple.
"""

from __future__ import annotations

import argparse
import csv
import multiprocessing
import os
import sys
import traceback
from bisect import bisect_left
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures.process import BrokenProcessPool
from dataclasses import dataclass, field, replace
from typing import Callable, Dict, List, Optional, Sequence, Set, Tuple

try:
    import primer3
except ImportError as exc:
    raise ImportError(
        "modZN requires the 'primer3-py' package. Install it with: pip install primer3-py"
    ) from exc


###############################################################
# CONSTANTS
###############################################################

MODIFICATION_MARK = "*"          # a phosphorothioate bond inside a primer sequence
VALID_BASES = "ACGT"
SEARCH_METHODS = ("exact", "greedy")
MIN_OVERLAP_LENGTH = 8
# primer3 cannot evaluate hairpins or dimers of longer sequences, and every overlap is
# evaluated in every run (the report states its structure), so this is a hard limit.
MAX_OVERLAP_LENGTH = 60
SIMILARITY_MISMATCH_FRACTION = 0.3   # overlaps differing at fewer positions than this are "too similar"
MAX_EXACT_SEARCH_NODES = 20_000_000  # safety limit of the exact search
# delta G values are rounded to 2 decimals, so differences carry float noise;
# the spread limit is compared with this tolerance.
DG_CAP_TOLERANCE = 1e-9
TM_MAX_NN_LENGTH = 200           # always use the nearest-neighbour model for Tm
MAX_AUTO_JOBS = 8                # most worker processes used when --jobs is left automatic
# Batches of fewer overlap pairs are evaluated in the calling process: starting workers
# would cost more than it saves.
MIN_PARALLEL_PAIRS = 24
# The exact search starts with the cross dimers evaluated so far. If one delta G window
# needs more search nodes than this, it evaluates all pairs inside the windows instead.
WINDOW_NODE_ALLOWANCE = 30_000
# Extra length (bases) a fragment with two neighbours needs so that the annealing region of its
# forward primer can reach the delta G of the reverse primer's.
MIN_ANNEALING_MARGIN = 30
SECONDARY_STRUCTURE_TEMP_C = 37.0    # temperature of the hairpin and dimer calculations
DEFAULT_HAIRPIN_DG_THRESHOLD = -5.0         # kcal/mol
DEFAULT_CROSS_DIMER_DG_THRESHOLD = -15.0   # kcal/mol
DEFAULT_SELF_DIMER_DG_THRESHOLD = -20.0    # kcal/mol
# Text shown in the report for each type of cross dimer.
CROSS_DIMER_MODE_LABEL = {
    "direct": "direct repeat",
    "inverted": "inverted repeat",
    "none": "no duplex found",
}


class DesignError(ValueError):
    """The input is invalid, or no design satisfies the constraints."""


class _WindowTooHard(Exception):
    """Internal signal: one window exceeded WINDOW_NODE_ALLOWANCE."""


class SearchBudgetExceeded(DesignError):
    """The exact search used up MAX_EXACT_SEARCH_NODES.

    This does not prove that no design exists, so callers may retry with fewer constraints.
    """


###############################################################
# PARAMETERS AND RESULT CONTAINERS
###############################################################

@dataclass(frozen=True)
class ThermoConditions:
    """Reaction conditions for the Tm, hairpin and dimer calculations."""
    na_mM: float = 67.0        # Na+ concentration, mM
    mg_mM: float = 2.5         # Mg2+ concentration, mM
    dntp_mM: float = 0.2       # dNTP concentration, mM
    oligo_uM: float = 0.02     # oligo concentration, uM

    def primer3_kwargs(self) -> Dict[str, float]:
        """The concentrations in the units primer3 expects (mM, mM, mM, nM)."""
        return dict(
            mv_conc=self.na_mM,
            dv_conc=self.mg_mM,
            dntp_conc=self.dntp_mM,
            dna_conc=self.oligo_uM * 1000.0,
        )


@dataclass(frozen=True)
class DesignParameters:
    """All settings of a design run."""
    overlap_length: int
    num_modifications: int
    search_method: str = "exact"               # "exact" or "greedy"
    hairpin_filter: bool = True
    hairpin_dg_threshold: float = DEFAULT_HAIRPIN_DG_THRESHOLD       # kcal/mol, 37 degC
    self_dimer_filter: bool = True
    self_dimer_dg_threshold: float = DEFAULT_SELF_DIMER_DG_THRESHOLD  # kcal/mol, 37 degC
    cross_dimer_filter: bool = True
    cross_dimer_dg_threshold: float = DEFAULT_CROSS_DIMER_DG_THRESHOLD  # kcal/mol, 37 degC
    max_dg_spread: float = 30.0                # largest allowed delta G difference between overlaps
    search_frame_length: int = 8               # how far an overlap may be shifted around the centre
    primer_extension_length: Optional[int] = None   # None: use the topology default
    # worker processes for the cross-dimer calculations:
    # None = automatic (all usable cores, at most MAX_AUTO_JOBS), 1 = no extra processes
    jobs: Optional[int] = None
    thermo: ThermoConditions = field(default_factory=ThermoConditions)


@dataclass
class OverlapInfo:
    """One designed overlapping zone."""
    junction: int              # 1-based junction number
    sequence: str
    left: int                  # bases taken from the upstream fragment
    right: int                 # bases taken from the downstream fragment
    delta_g: float             # stacking delta G used to compare overlaps (see OverlapSpace)
    tm: float
    hairpin_dg: float          # strongest hairpin of the two strands, kcal/mol at 37 degC
    self_dimer_dg: float = 0.0  # strongest self dimer of the two strands, kcal/mol at 37 degC


@dataclass
class Primer:
    """One designed primer, 5'->3'; ``*`` marks a phosphorothioate bond."""
    name: str
    sequence: str
    tm: float
    delta_g: float

    @property
    def length(self) -> int:
        """Length in bases (the ``*`` marks are not bases)."""
        return len(strip_modifications(self.sequence))


@dataclass
class DesignResult:
    """The complete result of a design run."""
    topology: str
    params: DesignParameters
    layout: Tuple[int, ...]
    dg_spread: float
    overlaps: List[OverlapInfo]
    primers: List[Primer]           # ordered F1-F, F1-R, F2-F, F2-R, ...
    warnings: List[str] = field(default_factory=list)
    # strongest cross dimer between two overlaps: (junction a, junction b, delta G, mode),
    # mode being "direct" or "inverted"
    worst_cross_dimer: Optional[Tuple[int, int, float, str]] = None


@dataclass(frozen=True)
class TopologySpec:
    """What differs between the linear and the circular tool."""
    name: str                                   # "linear" or "circular"
    circular: bool
    default_primer_extension_length: int
    design_primer_pairs: Callable[[Sequence[str], Tuple[int, ...], DesignParameters],
                                  List[Tuple[str, str]]]
    program_name: str                           # script name shown in --help


###############################################################
# SEQUENCE UTILITIES
###############################################################

_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C", MODIFICATION_MARK: MODIFICATION_MARK}


def reverse_complement(sequence: str) -> str:
    """Reverse complement; the ``*`` marks are mirrored with the sequence."""
    return "".join(_COMPLEMENT[base] for base in reversed(sequence))


def strip_modifications(sequence: str) -> str:
    """Remove all ``*`` marks."""
    return sequence.replace(MODIFICATION_MARK, "")


def spread_modifications(primer: str) -> str:
    """Spread the ``*`` marks over the bases that follow the first one.

    Input: ``OVERLAP***ANNEALING``. Output: ``OVERLAP A*N*N*...``, i.e. consecutive bonds
    starting one base into the annealing region (the bond at the overlap boundary stays
    unmodified, so the exonuclease stops one base inside the annealing region).
    """
    bases = list(primer)
    num_marks = bases.count(MODIFICATION_MARK)
    if num_marks == 0:
        return primer
    position = bases.index(MODIFICATION_MARK)
    bases = [base for base in bases if base != MODIFICATION_MARK]
    for _ in range(num_marks):
        bases.insert(position + 1, MODIFICATION_MARK)
        position += 2
    return "".join(bases)


def parse_fragment_sequences(text: str) -> List[str]:
    """Read fragment sequences, one per line, in assembly order.

    Blank lines and lines starting with ``>`` or ``#`` are skipped; spaces are removed
    and bases are upper-cased. Only A, C, G, T are accepted.
    """
    fragments: List[str] = []
    for line_number, raw_line in enumerate(text.splitlines(), start=1):
        line = "".join(raw_line.split())
        if not line or line.startswith((">", "#")):
            continue
        line = line.upper()
        invalid = sorted(set(line) - set(VALID_BASES))
        if invalid:
            raise DesignError(
                "Line %d contains characters other than A/C/G/T: %s" % (line_number, " ".join(invalid))
            )
        fragments.append(line)
    return fragments


###############################################################
# THERMODYNAMICS
###############################################################

# Nearest-neighbour stacking free energies (kcal/mol, 37 degC), by 5'->3' dinucleotide.
# Source: SantaLucia (1998) unified parameters.
NEAREST_NEIGHBOR_DG = {
    "AA": -1.00, "AT": -0.88, "AC": -1.44, "AG": -1.28,
    "TA": -0.58, "TT": -1.00, "TC": -1.30, "TG": -1.45,
    "CA": -1.45, "CT": -1.28, "CC": -1.84, "CG": -2.17,
    "GA": -1.30, "GT": -1.44, "GC": -2.24, "GG": -1.84,
}
TERMINAL_GC_DG = 0.98   # initiation term for a terminal G/C
TERMINAL_AT_DG = 1.03   # initiation term for a terminal A/T


def stacking_dg(sequence: str) -> float:
    """Sum of the stacking terms only (no initiation terms), rounded to 2 decimals.

    ``*`` marks are ignored.
    """
    sequence = strip_modifications(sequence)
    total = 0
    for i in range(len(sequence) - 1):
        value = NEAREST_NEIGHBOR_DG.get(sequence[i:i + 2])
        if value is not None:
            total = total + value
    return round(total, 2)


def primer_dg(primer: str) -> float:
    """Nearest-neighbour delta G of a primer: stacking plus both terminal terms.

    ``*`` marks are ignored.
    """
    sequence = strip_modifications(primer)
    if not sequence:
        raise DesignError("Cannot calculate delta G of an empty sequence.")
    total = 0
    for i in range(len(sequence) - 1):
        value = NEAREST_NEIGHBOR_DG.get(sequence[i:i + 2])
        if value is not None:
            total = total + value
    for terminal_base in (sequence[0], sequence[-1]):
        total = total + (TERMINAL_GC_DG if terminal_base in "GC" else TERMINAL_AT_DG)
    return round(total, 2)


def calculate_tm(sequence: str, conditions: ThermoConditions) -> float:
    """Melting temperature in C (primer3: SantaLucia nearest neighbour, Owczarzy salt
    correction). ``*`` marks are ignored."""
    return primer3.calc_tm(
        strip_modifications(sequence),
        tm_method="santalucia",
        salt_corrections_method="owczarzy",
        max_nn_length=TM_MAX_NN_LENGTH,
        **conditions.primer3_kwargs(),
    )


###############################################################
# OVERLAP FILTERS (hairpin, self dimer, cross dimer)
#
# After the exonuclease step every overlap exists in the reaction as two
# single-stranded overhangs: its top strand (on one fragment) and its bottom
# strand, the reverse complement (on the neighbouring fragment). They are
# different molecules, so the filters check both strands.
###############################################################

def _strongest(values) -> float:
    """The most negative delta G primer3 found; 0.0 if it found no structure at all."""
    found = [value for value in values if value is not None]
    return min(found) if found else 0.0


class HairpinFilter:
    """Rejects an overlap that folds into a stable hairpin.

    Both strands are checked and the stronger hairpin decides. An overlap is rejected only
    if its hairpin delta G (primer3, 37 degC) is more negative than ``dg_threshold``.
    """

    def __init__(self, conditions: ThermoConditions, dg_threshold: float):
        self.conditions = conditions
        self.dg_threshold = dg_threshold
        self._cache: Dict[str, float] = {}

    def hairpin_dg(self, sequence: str) -> float:
        """Strongest hairpin of the sequence or its reverse complement (0.0 if none)."""
        key = min(sequence, reverse_complement(sequence))      # both strands share one entry
        if key not in self._cache:
            kwargs = self.conditions.primer3_kwargs()
            self._cache[key] = _strongest(
                (result.dg / 1000.0) if result.structure_found else None
                for result in (primer3.calc_hairpin(strand, temp_c=SECONDARY_STRUCTURE_TEMP_C,
                                                    **kwargs)
                               for strand in (key, reverse_complement(key))))
        return self._cache[key]

    def passes(self, sequence: str) -> bool:
        return self.hairpin_dg(sequence) >= self.dg_threshold


class SelfDimerFilter:
    """Rejects an overlap that anneals to a copy of itself.

    Many copies of every fragment are in the reaction. If an overlap is partly
    self-complementary (e.g. a palindrome), its overhang can pair with the same overhang
    of another copy, joining a fragment to itself instead of to its neighbour. The hairpin
    filter does not see this (it needs a fold with a loop), nor does the cross-dimer
    filter (it compares different junctions).

    Both strands are checked. The pairing of the top strand with the bottom strand of the
    same overlap is the intended junction and is not counted.
    """

    def __init__(self, conditions: ThermoConditions, dg_threshold: float):
        self.conditions = conditions
        self.dg_threshold = dg_threshold
        self._cache: Dict[str, float] = {}

    def self_dimer_dg(self, sequence: str) -> float:
        """Strongest self dimer of the sequence or its reverse complement (0.0 if none)."""
        key = min(sequence, reverse_complement(sequence))      # both strands share one entry
        if key not in self._cache:
            kwargs = self.conditions.primer3_kwargs()
            self._cache[key] = _strongest(
                (result.dg / 1000.0) if result.structure_found else None
                for result in (primer3.calc_homodimer(strand, temp_c=SECONDARY_STRUCTURE_TEMP_C,
                                                      **kwargs)
                               for strand in (key, reverse_complement(key))))
        return self._cache[key]

    def passes(self, sequence: str) -> bool:
        return self.self_dimer_dg(sequence) >= self.dg_threshold


class OverlapFilter:
    """Hairpin and self-dimer filter combined: they judge one overlap on its own.

    Either may be missing. An overlap passes if it passes every filter that is present.
    Because they look at a single candidate, they are applied before the search by
    removing failing candidates.
    """

    def __init__(self, hairpin: Optional[HairpinFilter], self_dimer: Optional[SelfDimerFilter]):
        self.hairpin = hairpin
        self.self_dimer = self_dimer

    def passes(self, sequence: str) -> bool:
        return ((self.hairpin is None or self.hairpin.passes(sequence))
                and (self.self_dimer is None or self.self_dimer.passes(sequence)))

    def thresholds(self) -> str:
        """The active thresholds, for warning messages."""
        parts = []
        if self.hairpin is not None:
            parts.append("hairpin %.1f" % self.hairpin.dg_threshold)
        if self.self_dimer is not None:
            parts.append("self dimer %.1f" % self.self_dimer.dg_threshold)
        return ", ".join(parts) + " kcal/mol"

    def weakest_available(self, sequences: Sequence[str]) -> str:
        """The least problematic candidate per active check, for warning messages."""
        parts = []
        if self.hairpin is not None:
            parts.append("weakest hairpin %.2f kcal/mol"
                         % max(self.hairpin.hairpin_dg(s) for s in sequences))
        if self.self_dimer is not None:
            parts.append("weakest self dimer %.2f kcal/mol"
                         % max(self.self_dimer.self_dimer_dg(s) for s in sequences))
        return ", ".join(parts)


def resolve_jobs(jobs: Optional[int]) -> int:
    """Number of worker processes. ``None`` means automatic: all cores this process may
    use (this respects a cluster allocation or ``taskset``), at most ``MAX_AUTO_JOBS``."""
    if jobs is not None and jobs >= 1:
        return jobs
    try:
        cores = len(os.sched_getaffinity(0))
    except AttributeError:                      # platforms without sched_getaffinity
        cores = os.cpu_count() or 1
    return max(1, min(cores, MAX_AUTO_JOBS))


def _evaluate_pair(sequence_a: str, sequence_b: str, kwargs: Dict[str, float],
                   threshold: float, full: bool
                   ) -> Tuple[bool, Optional[Tuple[float, str]]]:
    """Evaluate the four duplexes of two overlaps (``sequence_a <= sequence_b``).

    Returns ``(passes, detail)``, where ``detail`` is ``(strongest delta G, mode)`` (see
    ``CrossDimerFilter.cross_dimer_detail``) or ``None`` if the evaluation stopped early.

    Unless ``full`` is set, it stops at the first duplex below ``threshold``: the pair is
    then known to fail, which saves up to three of the four primer3 calls. A pair that
    passes always has all four duplexes evaluated, so its detail is complete.

    ``passes`` is exactly "strongest delta G >= threshold" (0.0 if primer3 finds nothing).
    """
    worst, worst_mode = None, "none"
    for mode, first, second in CrossDimerFilter.duplexes(sequence_a, sequence_b):
        result = primer3.calc_heterodimer(first, second, temp_c=SECONDARY_STRUCTURE_TEMP_C,
                                          **kwargs)
        if result.structure_found:
            value = result.dg / 1000.0
            if worst is None or value < worst:
                worst, worst_mode = value, mode
            if not full and value < threshold:
                return False, None
    if worst is None:
        worst = 0.0                             # no duplex at all between these two overlaps
    return worst >= threshold, (worst, worst_mode)


def _evaluate_chunk(pairs: Sequence[Tuple[str, str]], kwargs: Dict[str, float],
                    threshold: float, full: bool):
    """Worker-process entry point: evaluate a list of overlap pairs."""
    return [_evaluate_pair(a, b, kwargs, threshold, full) for a, b in pairs]


class CrossDimerFilter:
    """Rejects two overlaps (of different junctions) that can anneal to each other.

    Two junctions can pair in four ways (see ``duplexes``). All four are evaluated with
    primer3 and the strongest decides.

    A pair is rejected only if its strongest cross dimer is more negative than
    ``dg_threshold``.
    """

    def __init__(self, conditions: ThermoConditions, dg_threshold: float, jobs: int = 1):
        self.conditions = conditions
        self.dg_threshold = dg_threshold
        self.jobs = max(1, jobs)
        # complete result of a pair: (strongest delta G, mode)
        self._cache: Dict[Tuple[str, str], Tuple[float, str]] = {}
        # verdict of a pair; may be known without a complete result (evaluation stopped early)
        self._verdicts: Dict[Tuple[str, str], bool] = {}

    @staticmethod
    def _key(sequence_a: str, sequence_b: str) -> Tuple[str, str]:
        """Cache key; the result does not depend on the order of the two overlaps."""
        return (sequence_a, sequence_b) if sequence_a <= sequence_b else (sequence_b, sequence_a)

    @staticmethod
    def duplexes(sequence_a: str, sequence_b: str) -> Tuple[Tuple[str, str, str], ...]:
        """The four duplexes two overlaps can form, as (mode, strand, strand).

            top(a) with bottom(b)       } direct repeat (a resembles b):
            bottom(a) with top(b)       }   fragments joined in the wrong order
            top(a) with top(b)          } inverted repeat (a resembles revcomp(b)):
            bottom(a) with bottom(b)    }   fragments joined head to head

        The two duplexes of a mode are different molecules (their strands are only
        reverse complements of each other); their unpaired ends differ, so their delta G
        differs too. All four really form in the reaction, so all four are evaluated.
        """
        rc_a, rc_b = reverse_complement(sequence_a), reverse_complement(sequence_b)
        return (
            ("direct", sequence_a, rc_b),      # top(a) : bottom(b)
            ("direct", rc_a, sequence_b),      # bottom(a) : top(b)
            ("inverted", sequence_a, sequence_b),   # top(a) : top(b)
            ("inverted", rc_a, rc_b),               # bottom(a) : bottom(b)
        )

    def _store(self, key: Tuple[str, str],
               result: Tuple[bool, Optional[Tuple[float, str]]]) -> None:
        passes, detail = result
        self._verdicts[key] = passes
        if detail is not None:
            self._cache[key] = detail

    def cross_dimer_detail(self, sequence_a: str, sequence_b: str) -> Tuple[float, str]:
        """Strongest cross dimer of two overlaps: (delta G in kcal/mol, mode).

        Mode is "direct" or "inverted"; (0.0, "none") if no duplex is found.
        """
        key = self._key(sequence_a, sequence_b)
        if key not in self._cache:
            self._store(key, _evaluate_pair(key[0], key[1], self.conditions.primer3_kwargs(),
                                            self.dg_threshold, True))
        return self._cache[key]

    def cross_dimer_dg(self, sequence_a: str, sequence_b: str) -> float:
        """Strongest cross dimer delta G in kcal/mol (0.0 if none)."""
        return self.cross_dimer_detail(sequence_a, sequence_b)[0]

    def passes(self, sequence_a: str, sequence_b: str) -> bool:
        """True if no duplex of the two overlaps is more negative than the threshold."""
        key = self._key(sequence_a, sequence_b)
        verdict = self._verdicts.get(key)
        if verdict is None:
            if key in self._cache:
                verdict = self._cache[key][0] >= self.dg_threshold
                self._verdicts[key] = verdict
            else:
                self._store(key, _evaluate_pair(key[0], key[1], self.conditions.primer3_kwargs(),
                                                self.dg_threshold, False))
                verdict = self._verdicts[key]
        return verdict

    def prefetch(self, pairs: Sequence[Tuple[str, str]], full: bool = False) -> None:
        """Evaluate many pairs at once (in worker processes if that pays off).

        The results go into the caches, so later ``passes`` / ``cross_dimer_detail`` calls
        are lookups. Each pair gets the same result as when evaluated alone: only the
        timing changes, never the design. ``full`` asks for complete results, not just
        verdicts.

        Processes are used instead of threads because primer3 keeps its working state in
        global variables and is not safe to run concurrently in one process.
        """
        todo = []
        seen = set()
        for sequence_a, sequence_b in pairs:
            key = self._key(sequence_a, sequence_b)
            if key in seen or key in self._cache or (not full and key in self._verdicts):
                continue
            seen.add(key)
            todo.append(key)
        if self.jobs <= 1 or len(todo) < MIN_PARALLEL_PAIRS:
            return                              # not worth starting workers: evaluate lazily
        kwargs = self.conditions.primer3_kwargs()
        size = max(1, -(-len(todo) // (self.jobs * 4)))          # about four chunks per worker
        chunks = [todo[i:i + size] for i in range(0, len(todo), size)]
        method = "fork" if "fork" in multiprocessing.get_all_start_methods() else None
        try:
            with ProcessPoolExecutor(max_workers=self.jobs,
                                     mp_context=multiprocessing.get_context(method)) as pool:
                futures = [pool.submit(_evaluate_chunk, chunk, kwargs, self.dg_threshold, full)
                           for chunk in chunks]
                results = [future.result() for future in futures]
        except (OSError, RuntimeError, BrokenProcessPool):
            return              # no workers available here: the lazy path gives the same values
        for chunk, chunk_results in zip(chunks, results):
            for key, result in zip(chunk, chunk_results):
                self._store(key, result)

    @property
    def evaluations(self) -> int:
        """Number of distinct pairs evaluated so far."""
        return len(self._cache.keys() | self._verdicts.keys())


###############################################################
# COMPATIBILITY BETWEEN OVERLAPS
###############################################################

Layout = Tuple[int, ...]

_BASE_CODE = {"A": 0, "C": 1, "G": 2, "T": 3}


def encode_sequence(sequence: str) -> int:
    """Pack a sequence into an integer, two bits per base, so that two sequences can be
    compared with one XOR (this makes the all-pairs similarity check fast)."""
    code = 0
    for index, base in enumerate(sequence):
        code |= _BASE_CODE[base] << (2 * index)
    return code


def _popcount(value: int) -> int:
    try:
        return value.bit_count()        # Python 3.10+
    except AttributeError:
        return bin(value).count("1")


def count_mismatches(code_a: int, code_b: int, low_bit_mask: int) -> int:
    """Number of differing bases between two packed sequences of equal length.

    ``low_bit_mask`` has the low bit of every two-bit field set (see ``low_bit_mask_for``).
    """
    difference = code_a ^ code_b
    return _popcount((difference | (difference >> 1)) & low_bit_mask)


def low_bit_mask_for(length: int) -> int:
    """Mask with the low bit of each of ``length`` two-bit fields set."""
    mask = 0
    for index in range(length):
        mask |= 1 << (2 * index)
    return mask


class ConflictTable:
    """Which pairs of candidate overlaps (at different junctions) cannot be used together.

    Two overlaps clash if they are positionally similar (fewer than
    ``SIMILARITY_MISMATCH_FRACTION`` of the positions differ) or, when a cross-dimer
    filter is given, if they can anneal to each other.

    The two checks have very different costs:

    * similarity is a few integer operations per pair, so all pairs are checked up front;
    * a cross dimer needs a primer3 call (about a millisecond), so pairs are evaluated on
      demand and cached. An ordinary search only looks at the pairs along its own path,
      far fewer than all of them.

    Only clashing pairs are stored.
    """

    def __init__(self, sequences: Sequence[Sequence[str]], overlap_length: int,
                 cross_dimer_filter: Optional[CrossDimerFilter]):
        self.sequences = sequences
        self.num_junctions = len(sequences)
        self.domain_size = len(sequences[0]) if sequences else 0
        self.cross_dimer_filter = cross_dimer_filter
        # clashes[(junction, index)][other junction] = indices at that junction that clash
        self.clashes: Dict[Tuple[int, int], Dict[int, Set[int]]] = {}
        self.num_conflicts = 0
        self.cross_dimer_queries = 0  # lookups, including cache hits (not primer3 calls)
        # junction pairs ("arcs") whose candidate combinations have all been evaluated
        self._materialised: Set[Tuple[int, int]] = set()
        # Arcs whose cross-dimer clashes are ignored because the design cannot satisfy them.
        # The similarity rule always stays in force.
        self.exempt_arcs: Set[Tuple[int, int]] = set()
        # the subset of exempt_arcs where *no* combination of candidates is compatible
        self.hopeless_arcs: Set[Tuple[int, int]] = set()
        self._dimer_entries: Set[Tuple[int, int, int, int]] = set()
        # The same clashes as bit masks over the candidate indices of the other junction
        # (bit i set = candidate i clashes), one set for the similarity rule and one for
        # the cross-dimer rule. They let the search test many candidates at once.
        self._hamming_masks: Dict[Tuple[int, int], Dict[int, int]] = {}
        self._dimer_masks: Dict[Tuple[int, int], Dict[int, int]] = {}

        max_mismatches = overlap_length * SIMILARITY_MISMATCH_FRACTION
        codes = [[encode_sequence(sequence) for sequence in per_junction]
                 for per_junction in sequences]
        mask = low_bit_mask_for(overlap_length)
        for ja in range(self.num_junctions - 1):
            for jb in range(ja + 1, self.num_junctions):
                for ia, code_a in enumerate(codes[ja]):
                    for ib, code_b in enumerate(codes[jb]):
                        if count_mismatches(code_a, code_b, mask) < max_mismatches:
                            self._add(ja, ia, jb, ib)

    def _add(self, ja: int, ia: int, jb: int, ib: int, from_cross_dimer: bool = False) -> None:
        self.clashes.setdefault((ja, ia), {}).setdefault(jb, set()).add(ib)
        self.clashes.setdefault((jb, ib), {}).setdefault(ja, set()).add(ia)
        masks = self._dimer_masks if from_cross_dimer else self._hamming_masks
        masks.setdefault((ja, ia), {})[jb] = masks.get((ja, ia), {}).get(jb, 0) | (1 << ib)
        masks.setdefault((jb, ib), {})[ja] = masks.get((jb, ib), {}).get(ja, 0) | (1 << ia)
        self.num_conflicts += 1
        if from_cross_dimer:
            key = (ja, ia, jb, ib) if (ja, ia) <= (jb, ib) else (jb, ib, ja, ia)
            self._dimer_entries.add(key)

    def _has(self, ja: int, ia: int, jb: int, ib: int) -> bool:
        return ib in self.clashes.get((ja, ia), {}).get(jb, ())

    def _is_dimer_entry(self, ja: int, ia: int, jb: int, ib: int) -> bool:
        key = (ja, ia, jb, ib) if (ja, ia) <= (jb, ib) else (jb, ib, ja, ia)
        return key in self._dimer_entries

    def incompatible(self, ja: int, ia: int, jb: int, ib: int) -> bool:
        """True if the two candidate overlaps must not be used together.

        Evaluates the cross dimer if it is not known yet.
        """
        exempt = (min(ja, jb), max(ja, jb)) in self.exempt_arcs
        if self._has(ja, ia, jb, ib):
            # an exempt arc drops the cross-dimer clash but keeps the similarity clash
            if not exempt or not self._is_dimer_entry(ja, ia, jb, ib):
                return True
            return False
        if self.cross_dimer_filter is None or exempt:
            return False
        if (min(ja, jb), max(ja, jb)) in self._materialised:
            return False            # the whole arc is known, so "no clash recorded" means compatible
        self.cross_dimer_queries += 1
        if self.cross_dimer_filter.passes(self.sequences[ja][ia], self.sequences[jb][ib]):
            return False
        self._add(ja, ia, jb, ib, from_cross_dimer=True)
        return True

    def arc_conflict_count(self, ja: int, jb: int) -> int:
        """How many candidate combinations between two junctions clash."""
        self.materialise_arc(ja, jb)
        return sum(1 for ia in range(self.domain_size) for ib in range(self.domain_size)
                   if self.incompatible(ja, ia, jb, ib))

    def hopeless_arc(self, ja: int, jb: int) -> bool:
        """True if *every* candidate combination between the two junctions clashes."""
        self.materialise_arc(ja, jb)
        return all(self.incompatible(ja, ia, jb, ib)
                   for ia in range(self.domain_size) for ib in range(self.domain_size))

    def exempt(self, ja: int, jb: int, hopeless: bool = False) -> None:
        """Stop applying the cross-dimer rule between these two junctions."""
        self.exempt_arcs.add((min(ja, jb), max(ja, jb)))
        if hopeless:
            self.hopeless_arcs.add((min(ja, jb), max(ja, jb)))

    def violates_when_exempt(self, ja: int, ia: int, jb: int, ib: int) -> bool:
        """True if this pair is only usable because its arc is exempt.

        The exemption is lifted for the check: while it is in force, the cross dimer
        is not looked at, so "no clash recorded" would prove nothing.
        """
        arc = (min(ja, jb), max(ja, jb))
        if arc not in self.exempt_arcs:
            return False
        self.exempt_arcs.discard(arc)
        try:
            return self.incompatible(ja, ia, jb, ib)
        finally:
            self.exempt_arcs.add(arc)

    def withdraw_exemption(self, ja: int, jb: int) -> None:
        """Apply the cross-dimer rule between these two junctions again."""
        arc = (min(ja, jb), max(ja, jb))
        self.exempt_arcs.discard(arc)
        self.hopeless_arcs.discard(arc)

    def clear_exemptions(self) -> None:
        """Apply the cross-dimer rule everywhere again.

        Exemptions belong to one search (they depend on which candidates it may use).
        Evaluated pairs stay cached.
        """
        self.exempt_arcs.clear()
        self.hopeless_arcs.clear()

    def materialise_arc(self, ja: int, jb: int) -> None:
        """Evaluate every candidate combination between two junctions."""
        key = (min(ja, jb), max(ja, jb))
        if key in self._materialised or self.cross_dimer_filter is None:
            self._materialised.add(key)
            return
        for ia in range(self.domain_size):
            for ib in range(self.domain_size):
                if self._has(ja, ia, jb, ib):
                    continue
                self.cross_dimer_queries += 1
                if not self.cross_dimer_filter.passes(self.sequences[ja][ia],
                                                      self.sequences[jb][ib]):
                    self._add(ja, ia, jb, ib, from_cross_dimer=True)
        self._materialised.add(key)

    def prefetch_pairs(self, pairs: Sequence[Tuple[int, int, int, int]]) -> None:
        """Evaluate many candidate pairs at once (see ``CrossDimerFilter.prefetch``).

        ``pairs`` are ``(junction a, index a, junction b, index b)``. Pairs that are already
        decided are skipped. This only changes when a pair is evaluated, not its outcome.
        """
        if self.cross_dimer_filter is None:
            return
        todo = []
        for ja, ia, jb, ib in pairs:
            arc = (min(ja, jb), max(ja, jb))
            if arc in self.exempt_arcs or arc in self._materialised or self._has(ja, ia, jb, ib):
                continue
            todo.append((self.sequences[ja][ia], self.sequences[jb][ib]))
        self.cross_dimer_filter.prefetch(todo)

    def materialise_arcs(self, arcs: Sequence[Tuple[int, int]]) -> None:
        """``materialise_arc`` for several arcs, evaluated together."""
        pending = sorted({(min(ja, jb), max(ja, jb)) for ja, jb in arcs} - self._materialised)
        self.prefetch_pairs([(ja, ia, jb, ib) for ja, jb in pending
                             for ia in range(self.domain_size) for ib in range(self.domain_size)])
        for ja, jb in pending:
            self.materialise_arc(ja, jb)

    def materialise_domain_pairs(self, domains: Sequence[Sequence[int]]) -> None:
        """Evaluate every pair of candidates from two different junctions' ``domains``.

        ``domains[j]`` lists the candidate indices of junction ``j``. Afterwards the
        clashes among them are completely known.
        """
        pairs = [(ja, ia, jb, ib)
                 for ja in range(len(domains) - 1) for jb in range(ja + 1, len(domains))
                 for ia in domains[ja] for ib in domains[jb]]
        self.prefetch_pairs(pairs)
        for ja, ia, jb, ib in pairs:
            self.incompatible(ja, ia, jb, ib)          # records the clashes just evaluated

    def clash_mask(self, junction: int, index: int, other_junction: int) -> int:
        """Bit mask of the candidates of ``other_junction`` *known* to clash with this one.

        Only evaluated pairs are known. The cross-dimer clashes of an exempt arc are left out.
        """
        entry = (junction, index)
        mask = self._hamming_masks.get(entry, {}).get(other_junction, 0)
        if (min(junction, other_junction), max(junction, other_junction)) in self.exempt_arcs:
            return mask
        return mask | self._dimer_masks.get(entry, {}).get(other_junction, 0)

    def clashing_indices(self, junction: int, index: int, other_junction: int) -> Set[int]:
        """Candidate indices of ``other_junction`` known to clash with this one."""
        stored = self.clashes.get((junction, index), {}).get(other_junction, frozenset())
        if not stored or (min(junction, other_junction),
                          max(junction, other_junction)) not in self.exempt_arcs:
            return stored
        return {other for other in stored
                if not self._is_dimer_entry(junction, index, other_junction, other)}

    @property
    def primer3_calls(self) -> int:
        """Number of distinct overlap pairs evaluated with primer3."""
        return self.cross_dimer_filter.evaluations if self.cross_dimer_filter else 0

    @property
    def is_empty(self) -> bool:
        return self.num_conflicts == 0


###############################################################
# OVERLAP SPACE
###############################################################


def search_window(overlap_length: int, search_frame_length: int) -> Tuple[int, ...]:
    """The allowed ``left`` values of a junction: a window around ``overlap_length // 2``."""
    first = int(overlap_length / 2 - search_frame_length / 2)
    stop = int(overlap_length / 2 + search_frame_length / 2 + 1)
    return tuple(range(max(0, first), min(overlap_length + 1, stop)))


class OverlapSpace:
    """All candidate overlaps of an assembly.

    At every junction, ``left`` (the number of bases taken from the upstream fragment) may
    vary in a small window around ``overlap_length // 2``; the remaining
    ``overlap_length - left`` bases come from the downstream fragment.

    Overlaps are compared by the stacking delta G of the overlap extended by one base on
    each side.
    """

    def __init__(self, fragments: Sequence[str], overlap_length: int, circular: bool,
                 search_frame_length: int, cross_dimer_filter: Optional["CrossDimerFilter"] = None):
        self.cross_dimer_filter = cross_dimer_filter
        self.fragments = list(fragments)
        self.num_fragments = len(self.fragments)
        self.overlap_length = overlap_length
        self.circular = circular
        self.num_junctions = self.num_fragments if circular else self.num_fragments - 1

        self.initial_left = overlap_length // 2
        self.left_choices = search_window(overlap_length, search_frame_length)
        self._index_of_left = {left: index for index, left in enumerate(self.left_choices)}

        self._sequences = [[self._build_overlap_sequence(j, left) for left in self.left_choices]
                           for j in range(self.num_junctions)]
        self._dg_table = [
            {left: stacking_dg(self.extended_sequence(j, left)) for left in self.left_choices}
            for j in range(self.num_junctions)
        ]
        self._conflicts: Optional[ConflictTable] = None

    @property
    def conflicts(self) -> ConflictTable:
        """The conflict table (built on first use)."""
        if self._conflicts is None:
            self._conflicts = ConflictTable(self._sequences, self.overlap_length,
                                            self.cross_dimer_filter)
        return self._conflicts

    def invalidate_conflicts(self) -> None:
        """Discard the conflict table (needed after the cross-dimer filter is detached)."""
        self._conflicts = None

    # -- sequences -------------------------------------------------------------------------

    def upstream_fragment(self, junction: int) -> str:
        return self.fragments[junction]

    def downstream_fragment(self, junction: int) -> str:
        return self.fragments[(junction + 1) % self.num_fragments]

    def _build_overlap_sequence(self, junction: int, left: int) -> str:
        upstream = self.upstream_fragment(junction)
        downstream = self.downstream_fragment(junction)
        return upstream[len(upstream) - left:] + downstream[:self.overlap_length - left]

    def index_of_left(self, left: int) -> int:
        """Position of a ``left`` value in ``left_choices``."""
        return self._index_of_left[left]

    def overlap_sequence(self, junction: int, left: int) -> str:
        """The overlap itself."""
        return self._sequences[junction][self._index_of_left[left]]

    def extended_sequence(self, junction: int, left: int) -> str:
        """The overlap extended by one base on each side (used for delta G)."""
        upstream = self.upstream_fragment(junction)
        downstream = self.downstream_fragment(junction)
        return upstream[len(upstream) - left - 1:] + downstream[:self.overlap_length - left + 1]

    # -- scoring ---------------------------------------------------------------------------

    def dg(self, junction: int, left: int) -> float:
        return self._dg_table[junction][left]

    def layout_dgs(self, layout: Layout) -> List[float]:
        return [self.dg(j, left) for j, left in enumerate(layout)]

    def dg_spread(self, layout: Layout) -> float:
        """Largest minus smallest overlap delta G of the layout."""
        dgs = self.layout_dgs(layout)
        return max(dgs) - min(dgs)

    # -- clashes between overlaps ---------------------------------------------------------------

    def cross_dimer_dg(self, junction_a: int, left_a: int, junction_b: int, left_b: int) -> float:
        """Strongest cross dimer of two overlaps (needs a cross-dimer filter)."""
        return self.cross_dimer_filter.cross_dimer_dg(self.overlap_sequence(junction_a, left_a),
                                                      self.overlap_sequence(junction_b, left_b))

    def too_similar(self, junction_a: int, left_a: int, junction_b: int, left_b: int) -> bool:
        """True if the two overlaps must not be used together (see ``ConflictTable``)."""
        return self.conflicts.incompatible(junction_a, self._index_of_left[left_a],
                                           junction_b, self._index_of_left[left_b])

    def layout_is_distinct(self, layout: Layout) -> bool:
        """True if no two overlaps of the layout clash."""
        for a in range(len(layout) - 1):
            for b in range(a + 1, len(layout)):
                if self.too_similar(a, layout[a], b, layout[b]):
                    return False
        return True

    def worst_cross_dimer(self, layout: Layout,
                          reporter: Optional["CrossDimerFilter"] = None
                          ) -> Optional[Tuple[int, int, float, str]]:
        """The strongest cross dimer between two overlaps of ``layout``.

        Returns ``(junction a, junction b, delta G, mode)`` (junctions numbered from 1),
        or ``None`` if there is only one junction. ``reporter`` lets the value be
        reported even when no cross-dimer filter was used.
        """
        measure = reporter or self.cross_dimer_filter
        if measure is None or len(layout) < 2:
            return None
        pairs = [(a, b) for a in range(len(layout) - 1) for b in range(a + 1, len(layout))]
        measure.prefetch([(self.overlap_sequence(a, layout[a]), self.overlap_sequence(b, layout[b]))
                          for a, b in pairs], full=True)
        candidates = []
        for a, b in pairs:
            dg, mode = measure.cross_dimer_detail(self.overlap_sequence(a, layout[a]),
                                                  self.overlap_sequence(b, layout[b]))
            candidates.append((a + 1, b + 1, dg, mode))
        return min(candidates, key=lambda item: item[2])


def build_overlap_space(fragments: Sequence[str], params: DesignParameters,
                        circular: bool, primer_extension_length: int,
                        cross_dimer_filter: Optional[CrossDimerFilter] = None) -> OverlapSpace:
    """Check the input and build the :class:`OverlapSpace`."""
    if len(fragments) < 2:
        raise DesignError("At least two fragments are required.")
    if params.overlap_length < MIN_OVERLAP_LENGTH:
        raise DesignError("The overlap length must be at least %d." % MIN_OVERLAP_LENGTH)
    if params.overlap_length > MAX_OVERLAP_LENGTH:
        raise DesignError(
            "The overlap length must not exceed %d: primer3, which evaluates the hairpins and "
            "dimers of every overlap, cannot handle longer sequences." % MAX_OVERLAP_LENGTH)
    if params.num_modifications < 0:
        raise DesignError("The number of modifications must not be negative.")
    if params.search_frame_length < 0:
        raise DesignError("The search frame length must not be negative.")
    if primer_extension_length < 1:
        raise DesignError("The primer extension length must be positive.")
    if params.search_method not in SEARCH_METHODS:
        raise DesignError("Unknown search method '%s' (use one of: %s)."
                          % (params.search_method, ", ".join(SEARCH_METHODS)))

    overlap_length = params.overlap_length
    window = search_window(overlap_length, params.search_frame_length)
    smallest_left, largest_left = window[0], window[-1]
    # shortest fragment for which every slice taken below stays inside the fragment
    min_length = max(
        largest_left + (overlap_length - smallest_left) + 2,      # both overlaps fit
        largest_left + primer_extension_length + 1,               # reverse-primer annealing region fits
        0 if circular else overlap_length + primer_extension_length,   # terminal primers (linear)
    )
    # A fragment with two neighbours gets a forward primer whose annealing region must be able to
    # grow as long as the reverse primer's; otherwise that region degenerates to a single base.
    length_with_neighbours = max(
        min_length,
        2 * (largest_left + primer_extension_length) + MIN_ANNEALING_MARGIN)
    for index, fragment in enumerate(fragments, start=1):
        has_two_neighbours = circular or 1 < index < len(fragments)
        required = length_with_neighbours if has_two_neighbours else min_length
        if len(fragment) < required:
            raise DesignError(
                "Fragment %d is too short (%d bp). With the current settings it must be at least "
                "%d bp long." % (index, len(fragment), required)
            )
    return OverlapSpace(fragments, overlap_length, circular, params.search_frame_length,
                        cross_dimer_filter)





###############################################################
# OVERLAP SEARCH
###############################################################

def greedy_search(space: OverlapSpace, max_dg_spread: float,
                  allowed: Optional[List[Set[int]]] = None) -> List[Layout]:
    """Greedy hill climbing.

    From the centred layout, repeatedly move only the weakest overlap (highest delta G) or the
    strongest one (lowest delta G) towards the others, depth first.  Every visited layout that
    meets the constraints is collected.  ``allowed`` optionally restricts the ``left`` values per
    junction (used by the overlap filters).  Returns the layouts in discovery order.
    """
    num_junctions = space.num_junctions
    choices = [[left for left in space.left_choices if allowed is None or left in allowed[j]]
               for j in range(num_junctions)]
    if any(not per_junction for per_junction in choices):
        return []
    visited: Set[Layout] = set()
    collected: List[Layout] = []
    collected_set: Set[Layout] = set()

    def moves(layout: Layout):
        """Lazily yield the neighbouring layouts that are worth visiting."""
        dgs = space.layout_dgs(layout)
        strongest = dgs.index(min(dgs))
        weakest = dgs.index(max(dgs))
        for left in choices[weakest]:        # improve the weakest overlap
            neighbour = layout[:weakest] + (left,) + layout[weakest + 1:]
            if neighbour not in visited and space.dg(weakest, left) < dgs[weakest]:
                yield neighbour
        for left in choices[strongest]:      # relax the strongest overlap
            neighbour = layout[:strongest] + (left,) + layout[strongest + 1:]
            if neighbour not in visited and space.dg(strongest, left) > dgs[strongest]:
                yield neighbour

    stack: List = []

    def visit(layout: Layout) -> None:
        visited.add(layout)
        if (space.dg_spread(layout) <= max_dg_spread + DG_CAP_TOLERANCE
                and layout not in collected_set
                and space.layout_is_distinct(layout)):
            collected.append(layout)
            collected_set.add(layout)
        stack.append(moves(layout))

    visit(tuple(min(per_junction, key=lambda left: (abs(left - space.initial_left), left))
                for per_junction in choices))
    while stack:
        try:
            neighbour = next(stack[-1])
        except StopIteration:
            stack.pop()
            continue
        visit(neighbour)
    return collected


def _mask_of(indices: Sequence[int]) -> int:
    """Bit mask with bit ``i`` set for every candidate index ``i``."""
    mask = 0
    for index in indices:
        mask |= 1 << index
    return mask


def arc_consistency(domains: List[List[int]], conflicts: ConflictTable,
                    junctions: Sequence[int]) -> Optional[List[List[int]]]:
    """AC-3 arc consistency: drop candidates that have no compatible partner at another junction.

    Deleting a candidate can remove the last partner of another one, so affected arcs are queued
    again until nothing changes.  Returns the reduced domains, or ``None`` if a domain becomes
    empty (the window has no layout).  ``domains`` holds candidate indices (positions in
    ``space.left_choices``).

    Only clashes that are already known are used; no new pair is evaluated here.  This is safe
    (a candidate is dropped only if every remaining partner is known to clash) and cheap.
    """
    reduced = [list(values) for values in domains]
    # Arcs between junctions already known to clash are processed first: they are the ones likely
    # to empty a domain, and the queue is abandoned as soon as that happens, so the remaining arcs
    # never have to be evaluated at all.
    suspect = {(min(ja, jb), max(ja, jb))
               for (ja, _), others in conflicts.clashes.items() for jb in others}
    def is_suspect(arc: Tuple[int, int]) -> bool:
        ja, jb = junctions[arc[0]], junctions[arc[1]]
        return (min(ja, jb), max(ja, jb)) in suspect

    # the queue is consumed from the end, so suspect arcs are placed last
    queue = sorted(((a, b) for a in range(len(junctions)) for b in range(len(junctions))
                    if a != b), key=is_suspect)
    queued = set(queue)
    alive = [_mask_of(values) for values in reduced]        # candidate sets as bit masks
    while queue:
        a, b = queue.pop()
        queued.discard((a, b))
        junction_a, junction_b = junctions[a], junctions[b]
        domain_b = alive[b]
        # a value of A has a partner unless every remaining value of B is known to clash with it
        survivors = [value for value in reduced[a]
                     if domain_b & ~conflicts.clash_mask(junction_a, value, junction_b)]
        if len(survivors) == len(reduced[a]):
            continue
        if not survivors:
            return None
        reduced[a] = survivors
        alive[a] = _mask_of(survivors)
        for c in range(len(junctions)):
            if c != a and c != b and (c, a) not in queued:
                queue.append((c, a))
                queued.add((c, a))
    return reduced


def exact_search(space: OverlapSpace, max_dg_spread: float,
                 allowed: Optional[List[Set[int]]] = None) -> Optional[Layout]:
    """Exact search for the layout with the smallest delta G spread.

    Each junction picks one ``left`` from its search window.  A layout is valid if all overlap
    delta G values fit in a window ``[low, high]`` and no two overlaps clash.  The answer is the
    valid layout with the smallest ``high - low``; ties go to the lowest ``low``, then to the most
    centred candidates, so the result is fully determined by the input.

    A two-pointer sweep over the distinct delta G values finds the window (both pointers only
    move forward).  Each window is tested by :func:`arc_consistency` plus a depth-first search
    with forward checking.  Cross dimers are evaluated lazily; if one window is too hard
    (more than ``WINDOW_NODE_ALLOWANCE`` nodes), they are evaluated up front, in parallel if
    possible.  This only changes when pairs are evaluated, never the result.

    ``allowed`` optionally restricts the ``left`` values per junction.  Returns ``None`` if no
    layout satisfies the constraints.
    """
    num_junctions = space.num_junctions
    options: List[List[Tuple[float, int]]] = []     # per junction: (delta G, index), by delta G
    for j in range(num_junctions):
        choices = [(space.dg(j, left), index) for index, left in enumerate(space.left_choices)
                   if allowed is None or left in allowed[j]]
        if not choices:
            return None
        options.append(sorted(choices))
    sorted_dgs = [[dg for dg, _ in per_junction] for per_junction in options]

    pool = sorted({dg for per_junction in options for dg, _ in per_junction})
    conflicts = space.conflicts
    # exemptions are per search: they depend on which candidates ``allowed`` leaves
    conflicts.clear_exemptions()
    center = space.initial_left
    left_of = space.left_choices
    budget = [MAX_EXACT_SEARCH_NODES]
    junctions = range(num_junctions)

    def covers_all_junctions(low: float, high: float) -> bool:
        for dgs in sorted_dgs:
            position = bisect_left(dgs, low)
            if position == len(dgs) or dgs[position] > high:
                return False
        return True

    # per junction: (delta G, index), most centred candidate first
    by_centredness = [sorted(per_junction, key=lambda item: (abs(left_of[item[1]] - center),
                                                             item[1]))
                      for per_junction in options]

    def window_domains(low: float, high: float) -> List[List[int]]:
        """Candidate indices per junction, best (most centred) first."""
        return [[index for dg, index in per_junction if low <= dg <= high]
                for per_junction in by_centredness]

    eager = [False]          # True once a window was too hard to settle lazily

    def backtrack(domains: List[List[int]], allowance: Optional[int] = None) -> Optional[Layout]:
        """Depth-first search inside one window, limited by the global node budget.

        Junctions are assigned in order, most centred candidate first.  Forward checking removes,
        at each choice, the later candidates known to clash with it and abandons the branch when a
        later junction has none left.  It only prunes branches without a solution, so the layout
        found is the same as in a plain depth-first search.  Pairs not yet evaluated are checked
        when their second junction is reached.
        """
        chosen: List[int] = []
        nodes = [0]

        def extend(level: int, alive: List[int]) -> Optional[Layout]:
            if level == num_junctions:
                return tuple(left_of[index] for index in chosen)
            live = alive[level]
            for index in domains[level]:
                if not (live >> index) & 1:
                    continue                   # already known to clash with an earlier choice
                nodes[0] += 1
                if allowance is not None and nodes[0] > allowance:
                    raise _WindowTooHard()
                budget[0] -= 1
                if budget[0] < 0:
                    raise SearchBudgetExceeded(
                        "The exact search exceeded its time budget. This happens when the "
                        "constraints are very hard to satisfy, for example with highly repetitive "
                        "sequences or when many pairs of overlaps are rejected by the cross-dimer "
                        "filter. Try a smaller --search-frame-length or a different overlap "
                        "length."
                    )
                if any(conflicts.incompatible(earlier, chosen[earlier], level, index)
                       for earlier in range(level)):
                    continue                   # a pair that had not been evaluated yet clashes
                narrowed = alive[:]
                for later in range(level + 1, num_junctions):
                    remaining = alive[later] & ~conflicts.clash_mask(level, index, later)
                    if not remaining:
                        break                  # a later junction would be left without a candidate
                    narrowed[later] = remaining
                else:
                    chosen.append(index)
                    found = extend(level + 1, narrowed)
                    if found is not None:
                        return found
                    chosen.pop()
            return None

        return extend(0, [_mask_of(per_junction) for per_junction in domains])

    def centred_layout(low: float, high: float) -> Optional[Layout]:
        """Most centred layout inside the window, ignoring pairwise constraints."""
        chosen = []
        for per_junction in by_centredness:
            for dg, index in per_junction:
                if low <= dg <= high:
                    chosen.append(left_of[index])
                    break
            else:
                return None            # this junction has no candidate inside the window
        return tuple(chosen)

    def propagate_with_knowledge(domains: List[List[int]]) -> Optional[List[List[int]]]:
        """Arc consistency that evaluates the cross dimers it needs.

        Alternate propagation and evaluation of the pairs among still-alive candidates until
        nothing changes.  Removed candidates are never evaluated; results stay cached for later
        windows.  Returns ``None`` if a junction runs out of candidates (window refuted).
        """
        current = arc_consistency(domains, conflicts, junctions)
        while current is not None:
            conflicts.materialise_domain_pairs(current)
            narrowed = arc_consistency(current, conflicts, junctions)
            if narrowed is None:
                return None
            if [len(values) for values in narrowed] == [len(values) for values in current]:
                return narrowed        # no candidate disappeared: the knowledge is complete here
            current = narrowed
        return None

    def solve_window(low: float, high: float) -> Optional[Layout]:
        """The best layout inside this window, or ``None`` if the window admits none."""
        if not covers_all_junctions(low, high):
            return None
        domains = window_domains(low, high)
        if conflicts.is_empty and conflicts.cross_dimer_filter is None:
            # no binary constraints at all: every junction simply takes its most centred candidate
            return tuple(left_of[per_junction[0]] for per_junction in domains)
        # Propagation first: it can refute a window outright.  Lazily it only uses cross dimers
        # evaluated so far; once a window proved too hard (``eager``), the conflicts among all
        # candidates are evaluated up front, which makes propagation plus a short search enough.
        if eager[0]:
            propagated = propagate_with_knowledge(domains)
        else:
            propagated = arc_consistency(domains, conflicts, junctions)
        if propagated is None:
            return None                # propagation alone refuted the window
        try:
            return backtrack(propagated, None if eager[0] else WINDOW_NODE_ALLOWANCE)
        except _WindowTooHard:
            eager[0] = True
            propagated = propagate_with_knowledge(domains)
            if propagated is None:
                return None
            return backtrack(propagated)

    def sweep(window_solver) -> Tuple[Optional[float], Optional[Layout]]:
        """Two-pointer sweep: for every ``low``, find the smallest feasible ``high``.

        Widening a window never makes it infeasible and raising ``low`` never lowers the smallest
        feasible ``high``, so both pointers only move forward.  Windows that cannot beat the best
        width so far are skipped.
        """
        best_width: Optional[float] = None
        best_layout: Optional[Layout] = None
        limit = max_dg_spread          # no window wider than this is worth testing
        high_index = 0
        for low_index, low in enumerate(pool):
            if high_index < low_index:
                high_index = low_index
            layout = None
            while high_index < len(pool):
                if pool[high_index] - low > limit + DG_CAP_TOLERANCE:
                    break              # no window starting at this low can still be useful
                layout = window_solver(low, pool[high_index])
                if layout is not None:
                    break
                layout = None
                high_index += 1
            if high_index == len(pool):
                break                  # infeasible here, hence for every larger low as well
            if layout is None:
                continue               # stopped on the width limit, not on a solution
            width = pool[high_index] - low
            if best_width is None or width < best_width:
                best_width, best_layout = width, layout
                limit = min(max_dg_spread, best_width)
        return best_width, best_layout

    # Optimistic pass: find the optimum ignoring pairwise constraints (no primer3 calls).  Such
    # constraints can only worsen it, so if this layout satisfies them it is the answer.
    _optimistic_width, optimistic = sweep(centred_layout)
    if optimistic is None:
        return None
    indices = [space.index_of_left(left) for left in optimistic]
    # evaluate all pairs together for the check below
    conflicts.prefetch_pairs([(a, indices[a], b, indices[b]) for a in range(num_junctions - 1)
                              for b in range(a + 1, num_junctions)])
    clashing_pairs = [(a, b) for a in range(num_junctions - 1)
                      for b in range(a + 1, num_junctions)
                      if conflicts.incompatible(a, indices[a], b, indices[b])]
    if not clashing_pairs:
        return optimistic

    if conflicts.cross_dimer_filter is None:
        # Only the similarity rule is active (never exempted): nothing to give up.
        return sweep(solve_window)[1]

    # A junction pair where no combination of overlaps is compatible can never be satisfied.
    # It must also clash in the optimistic layout, so only those pairs are tested.  The caller
    # reports exempted pairs via :attr:`ConflictTable.exempt_arcs`.
    # evaluate all candidate combinations of the clashing pairs together
    conflicts.materialise_arcs(clashing_pairs)
    for a, b in clashing_pairs:
        if conflicts.hopeless_arc(a, b):
            conflicts.exempt(a, b, hopeless=True)

    def keep_only_violated(layout: Layout) -> Layout:
        """Undo exemptions the layout did not need, so only real losses are reported."""
        indices = [space.index_of_left(left) for left in layout]
        for arc in sorted(conflicts.exempt_arcs):
            x, y = arc
            if not conflicts.violates_when_exempt(x, indices[x], y, indices[y]):
                conflicts.withdraw_exemption(x, y)
        return layout

    layout = sweep(solve_window)[1]
    if layout is not None:
        return keep_only_violated(layout)

    # Still nothing: a pair may be unsatisfiable only within the usable windows.  Give up one
    # pair at a time, most constraining first, until a layout exists.
    remaining = sorted((pair for pair in clashing_pairs
                        if tuple(sorted(pair)) not in conflicts.exempt_arcs),
                       key=lambda pair: -conflicts.arc_conflict_count(*pair))
    for a, b in remaining:
        conflicts.exempt(a, b)
        layout = sweep(solve_window)[1]
        if layout is not None:
            return keep_only_violated(layout)
    return None


@dataclass
class LayoutSelection:
    layout: Layout
    warnings: List[str]


def search_layout(space: OverlapSpace, params: DesignParameters,
                  allowed: Optional[List[Set[int]]]) -> Optional[Layout]:
    """Run the configured search method; ``None`` if ``allowed`` admits no valid layout."""
    if params.search_method == "greedy":
        ranked = sorted(greedy_search(space, params.max_dg_spread, allowed), key=space.dg_spread)
        return ranked[0] if ranked else None
    return exact_search(space, params.max_dg_spread, allowed)


def _select_with_overlap_filter_relaxation(space: OverlapSpace, params: DesignParameters,
                                           overlap_filter: Optional[OverlapFilter]
                                           ) -> Optional[LayoutSelection]:
    """Best overlap layout under the constraints attached to ``space``, or ``None``.

    Overlaps failing the hairpin/self-dimer filter are removed before the search.  If that leaves
    no layout, the filter is relaxed junction by junction (not globally):

    1. junctions where no candidate passes are released first;
    2. then the remaining junctions are released one at a time, most constrained first.

    Each relaxation produces a warning naming the junctions and the least problematic candidate.
    """
    warnings: List[str] = []

    if overlap_filter is None:
        layout = search_layout(space, params, None)
        return LayoutSelection(layout, warnings) if layout is not None else None

    num_junctions = space.num_junctions
    all_choices = set(space.left_choices)
    passing = [{left for left in space.left_choices
                if overlap_filter.passes(space.overlap_sequence(j, left))}
               for j in range(num_junctions)]

    def describe(junctions: Sequence[int]) -> str:
        return ", ".join(
            "%d (%s)" % (j + 1, overlap_filter.weakest_available(
                [space.overlap_sequence(j, left) for left in space.left_choices]))
            for j in sorted(junctions))

    released: Set[int] = set()

    def allowed_sets() -> List[Set[int]]:
        return [all_choices if j in released else passing[j] for j in range(num_junctions)]

    def relaxation_warning(impossible: Sequence[int]) -> str:
        """Warning text; built after the search so it matches the returned design."""
        extra = sorted(released.difference(impossible))
        parts = []
        if impossible:
            parts.append("no overlap passes it at junction%s %s"
                         % ("" if len(impossible) == 1 else "s", describe(impossible)))
        if extra:
            parts.append("no layout exists with the filter on at every remaining junction, so it "
                         "was also relaxed at junction%s %s"
                         % ("" if len(extra) == 1 else "s", describe(extra)))
        still_filtered = [j for j in range(num_junctions) if j not in released]
        if still_filtered:
            closing = ("Junction%s %s %s still filtered."
                       % ("" if len(still_filtered) == 1 else "s",
                          ", ".join(str(j + 1) for j in still_filtered),
                          "is" if len(still_filtered) == 1 else "are"))
        else:
            closing = ("The filter is now off at every junction, i.e. the reported layout is the "
                       "best layout without hairpin and self-dimer filtering.")
        return ("Overlap filter (%s): %s. %s"
                % (overlap_filter.thresholds(), "; ".join(parts), closing))

    # Step 1: junctions with no passing overlap can never be satisfied.
    impossible = [j for j in range(num_junctions) if not passing[j]]
    released.update(impossible)

    layout = search_layout(space, params, allowed_sets())
    if layout is not None:
        if impossible:
            warnings.append(relaxation_warning(impossible))
        return LayoutSelection(layout, warnings)

    # Step 2: release the others one by one, most constrained first.
    for junction in sorted((j for j in range(num_junctions) if j not in released),
                           key=lambda j: (len(passing[j]), j)):
        released.add(junction)
        layout = search_layout(space, params, allowed_sets())
        if layout is None:
            continue
        warnings.append(relaxation_warning(impossible))
        return LayoutSelection(layout, warnings)

    return None


def select_overlap_layout(space: OverlapSpace, params: DesignParameters,
                          overlap_filter: Optional[OverlapFilter]) -> LayoutSelection:
    """Select the best overlap layout, relaxing filters only as far as necessary.

    Hairpin/self-dimer filters are relaxed per junction first.  If that fails and a cross-dimer
    filter is attached, that filter is switched off as a last resort and the search repeated,
    with a warning giving the strongest cross dimer in the result.  Exhausting the exact
    search's node budget counts as "no layout found" and is handled the same way.
    """
    cross_dimer_filter = space.cross_dimer_filter
    exhausted_budget = False
    try:
        selection = _select_with_overlap_filter_relaxation(space, params, overlap_filter)
    except SearchBudgetExceeded:
        if cross_dimer_filter is None:
            raise                      # nothing left to relax: the sequences themselves are hard
        selection, exhausted_budget = None, True

    if selection is not None:
        table = space.conflicts
        hopeless = sorted(table.hopeless_arcs)
        relaxed = sorted(table.exempt_arcs - table.hopeless_arcs)
        if hopeless or relaxed:
            names = lambda arcs: ", ".join("%d/%d" % (a + 1, b + 1) for a, b in arcs)
            parts = []
            if hopeless:
                # may stem from the cross-dimer or similarity rule, so no attribution
                parts.append("no shift of the overlaps inside the search frame makes junction%s %s "
                             "compatible, so the cross-dimer part of that constraint was dropped "
                             "there; separating them needs different fragment boundaries, a "
                             "different overlap length or a larger search frame"
                             % ("" if len(hopeless) == 1 else "s", names(hopeless)))
            if relaxed:
                parts.append("no layout satisfies the filter at every remaining junction pair, so "
                             "it was additionally given up between junction%s %s"
                             % ("" if len(relaxed) == 1 else "s", names(relaxed)))
            # one message, so the closing statement holds for the whole relaxation
            selection.warnings.append(
                "Cross-dimer filter (threshold %.1f kcal/mol): %s. Every other pair of overlaps is "
                "still filtered."
                % (cross_dimer_filter.dg_threshold, "; ".join(parts)))
        return selection

    if cross_dimer_filter is not None:
        space.cross_dimer_filter = None
        space.invalidate_conflicts()        # the cached table was built with the filter attached
        try:
            selection = _select_with_overlap_filter_relaxation(space, params, overlap_filter)
        finally:
            space.cross_dimer_filter = cross_dimer_filter
            space.invalidate_conflicts()
        if selection is not None:
            worst = space.worst_cross_dimer(selection.layout, reporter=cross_dimer_filter)
            if worst is None:        # single junction: there is no pair to report
                return selection
            selection.warnings.insert(0, (
                "Cross-dimer filter (threshold %.1f kcal/mol): %s, so the cross-dimer filter was "
                "switched off. The reported layout contains a cross dimer of %.2f kcal/mol between "
                "junctions %d and %d (%s), which may cause misassembly; consider different "
                "fragment boundaries, a different overlap length or a larger search frame.%s"
                % (cross_dimer_filter.dg_threshold,
                   "the search could not satisfy it within its time budget" if exhausted_budget
                   else "no overlap layout keeps every pair of overlaps below it",
                   worst[2], worst[0], worst[1], CROSS_DIMER_MODE_LABEL[worst[3]],
                   # exact gives the filter up per junction pair, keeping most protection
                   " Use --search-method exact to give the filter up only between the junction "
                   "pairs that cannot be satisfied, instead of for the whole assembly."
                   if params.search_method == "greedy" else "")))
            return selection

    raise DesignError(
        "No overlap layout satisfies the constraints (maximum delta G spread %.1f kcal/mol). "
        "Try a different overlap length or a larger search frame." % params.max_dg_spread
    )


###############################################################
# PRIMER HELPERS
###############################################################

def find_annealing_end(fragment: str, first: int, stop: int, target_dg: float) -> int:
    """Length of the 5' part of ``fragment`` whose delta G is closest to ``target_dg``.

    Scans lengths from ``first`` to ``stop - 1`` for the last crossing of the target and returns
    the closer side.  If the target is never crossed, returns ``first``.
    """
    end = first
    for j in range(first, stop):
        dg_short = primer_dg(fragment[:j])
        dg_long = primer_dg(fragment[:j + 1])
        if dg_short - target_dg >= 0 and dg_long - target_dg <= 0:
            end = j if -dg_short + target_dg > dg_long - target_dg else j + 1
    return end


def design_reverse_primer(fragments: Sequence[str], layout: Layout, index: int,
                          params: DesignParameters) -> str:
    """Reverse primer of fragment ``index``: 5' end = overlap with the next fragment."""
    fragment = fragments[index]
    length = len(fragment)
    extension = params.primer_extension_length
    left = layout[index]
    right = params.overlap_length - left
    following = fragments[(index + 1) % len(fragments)]
    marks = MODIFICATION_MARK * params.num_modifications
    template = (fragment[length - left - extension:length - left] + marks
                + fragment[length - left:] + following[:right])
    return spread_modifications(reverse_complement(template))


def design_forward_primer(fragments: Sequence[str], layout: Layout, index: int,
                          params: DesignParameters) -> str:
    """Forward primer of fragment ``index``: 5' end = overlap with the previous fragment.

    The annealing region is sized so its delta G matches the reverse primer's.
    """
    num_fragments = len(fragments)
    previous = (index - 1) % num_fragments
    fragment = fragments[index]
    previous_fragment = fragments[previous]
    length = len(fragment)
    extension = params.primer_extension_length
    left = layout[index]
    previous_left = layout[previous]
    previous_right = params.overlap_length - previous_left
    marks = MODIFICATION_MARK * params.num_modifications

    target_dg = primer_dg(fragment[length - left - extension:])
    annealing_end = find_annealing_end(fragment, previous_right + 1, length - left - extension, target_dg)
    return spread_modifications(
        previous_fragment[len(previous_fragment) - previous_left:] + fragment[:previous_right] + marks
        + fragment[previous_right:annealing_end]
    )


def design_standard_primer_pair(fragments: Sequence[str], layout: Layout, index: int,
                                params: DesignParameters) -> Tuple[str, str]:
    """(forward, reverse) primers of a fragment that overlaps with both of its neighbours."""
    return (design_forward_primer(fragments, layout, index, params),
            design_reverse_primer(fragments, layout, index, params))


###############################################################
# DESIGN PIPELINE
###############################################################

def run_design(fragments: Sequence[str], params: DesignParameters, spec: TopologySpec) -> DesignResult:
    """Run the complete design (overlaps, primers, Tm) for the given topology."""
    if params.primer_extension_length is None:
        params = replace(params, primer_extension_length=spec.default_primer_extension_length)

    cross_dimer_reporter = CrossDimerFilter(params.thermo, params.cross_dimer_dg_threshold,
                                            jobs=resolve_jobs(params.jobs))
    space = build_overlap_space(fragments, params, spec.circular, params.primer_extension_length,
                                cross_dimer_reporter if params.cross_dimer_filter else None)
    hairpin_filter = (HairpinFilter(params.thermo, params.hairpin_dg_threshold)
                      if params.hairpin_filter else None)
    self_dimer_filter = (SelfDimerFilter(params.thermo, params.self_dimer_dg_threshold)
                         if params.self_dimer_filter else None)
    overlap_filter = (OverlapFilter(hairpin_filter, self_dimer_filter)
                      if hairpin_filter is not None or self_dimer_filter is not None else None)
    selection = select_overlap_layout(space, params, overlap_filter)
    layout = selection.layout

    # the report lists hairpin/self-dimer values even when the filters are off
    reporter = hairpin_filter or HairpinFilter(params.thermo, params.hairpin_dg_threshold)
    self_dimer_reporter = (self_dimer_filter
                           or SelfDimerFilter(params.thermo, params.self_dimer_dg_threshold))
    overlaps = []
    for junction, left in enumerate(layout):
        sequence = space.overlap_sequence(junction, left)
        overlaps.append(OverlapInfo(
            junction=junction + 1,
            sequence=sequence,
            left=left,
            right=params.overlap_length - left,
            delta_g=space.dg(junction, left),
            tm=calculate_tm(sequence, params.thermo),
            hairpin_dg=reporter.hairpin_dg(sequence),
            self_dimer_dg=self_dimer_reporter.self_dimer_dg(sequence),
        ))

    primers: List[Primer] = []
    for index, (forward, reverse) in enumerate(spec.design_primer_pairs(fragments, layout, params), start=1):
        primers.append(Primer("F%d-F" % index, forward, calculate_tm(forward, params.thermo),
                              primer_dg(forward)))
        # reverse-primer delta G is taken on the template strand it anneals to
        primers.append(Primer("F%d-R" % index, reverse, calculate_tm(reverse, params.thermo),
                              primer_dg(reverse_complement(reverse))))

    return DesignResult(
        topology=spec.name, params=params, layout=layout, dg_spread=space.dg_spread(layout),
        overlaps=overlaps, primers=primers, warnings=selection.warnings,
        worst_cross_dimer=space.worst_cross_dimer(layout, reporter=cross_dimer_reporter),
    )


###############################################################
# OUTPUT
###############################################################

def format_report(result: DesignResult) -> str:
    """Human-readable report of a design result."""
    params = result.params
    lines = [
        "modZN design report (%s assembly)" % result.topology,
        "Search method: %s" % params.search_method,
        "Hairpin filter: %s" % ("on (threshold %.1f kcal/mol)" % params.hairpin_dg_threshold
                                if params.hairpin_filter else "off"),
        "Self-dimer filter: %s" % ("on (threshold %.1f kcal/mol)" % params.self_dimer_dg_threshold
                                   if params.self_dimer_filter else "off"),
        "Cross-dimer filter: %s" % ("on (threshold %.1f kcal/mol)" % params.cross_dimer_dg_threshold
                                    if params.cross_dimer_filter else "off"),
        "Distribution Information: %s" % [[left, params.overlap_length - left] for left in result.layout],
        "Delta G spread between overlaps: %.2f kcal/mol" % result.dg_spread,
    ]
    if result.worst_cross_dimer is not None:
        junction_a, junction_b, worst_dg, mode = result.worst_cross_dimer
        lines.append("Strongest cross dimer between overlaps: %.2f kcal/mol "
                     "(junctions %d and %d, %s)"
                     % (worst_dg, junction_a, junction_b, CROSS_DIMER_MODE_LABEL[mode]))
    for warning in result.warnings:
        lines.append("WARNING: " + warning)
    lines += ["", "Overlapping Zone:"]
    for overlap in result.overlaps:
        lines += [
            overlap.sequence,
            "Delta G (kcal/mol): %s" % overlap.delta_g,
            "Tm Value (degC): %.1f" % overlap.tm,
            "Hairpin Delta G (kcal/mol): %.2f" % overlap.hairpin_dg,
            "Self-dimer Delta G (kcal/mol): %.2f" % overlap.self_dimer_dg,
            "",
        ]
    for fragment_index in range(len(result.primers) // 2):
        forward = result.primers[2 * fragment_index]
        reverse = result.primers[2 * fragment_index + 1]
        lines.append("Fragment %d:" % (fragment_index + 1))
        for label, primer in (("Forward Primer", forward), ("Reverse Primer", reverse)):
            lines += [
                "%s: %s" % (label, primer.sequence),
                "Tm Value (degC): %.1f" % primer.tm,
                "Delta G (kcal/mol): %s" % primer.delta_g,
            ]
        lines.append("")
    return "\n".join(lines)


def write_primer_csv(result: DesignResult, path: str) -> None:
    """Write the primer table (name, sequence, length, Tm, delta G) as CSV."""
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Primer Name", "Sequence", "Length", "Tm Value (degC)", "Delta G Value (kcal/mol)"])
        for primer in result.primers:
            writer.writerow([primer.name, primer.sequence, primer.length,
                             "%.1f" % primer.tm, primer.delta_g])


###############################################################
# COMMAND-LINE INTERFACE
###############################################################

def build_argument_parser(spec: TopologySpec) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=spec.program_name,
        description="modZN: overlap and primer design for %s modAssembly. "
                    "Run without arguments (or with --gui) to open the graphical interface."
                    % spec.name,
    )
    parser.add_argument("--gui", action="store_true", help="open the graphical interface")
    parser.add_argument("-i", "--input", help="text file with one fragment sequence per line, in assembly order")
    parser.add_argument("-l", "--overlap-length", type=int, help="length of the overlapping zones")
    parser.add_argument("-m", "--num-modifications", type=int,
                        help="number of phosphorothioate modifications per primer")
    parser.add_argument("-o", "--output-csv", default="data.csv",
                        help="output CSV file with the primers (default: %(default)s)")
    parser.add_argument("--report", help="optionally write the full text report to this file")

    search = parser.add_argument_group("overlap search")
    search.add_argument("--search-method", choices=SEARCH_METHODS, default="exact",
                        help="'exact' finds the global optimum, 'greedy' is a hill-climbing "
                             "search (default: %(default)s)")
    search.add_argument("--max-dg-spread", type=float, default=30.0,
                        help="maximum delta G difference between overlaps, kcal/mol (default: %(default)s)")
    search.add_argument("--search-frame-length", type=int, default=8,
                        help="width of the search window around the centred overlap (default: %(default)s)")
    search.add_argument("--no-hairpin-filter", action="store_true",
                        help="do not filter overlaps by hairpin strength")
    search.add_argument("--hairpin-dg-threshold", type=float, default=DEFAULT_HAIRPIN_DG_THRESHOLD,
                        help="overlaps whose hairpin delta G (37 degC, kcal/mol, strongest of the two "
                             "strands) is more negative than this are rejected "
                             "(default: %(default)s)")
    search.add_argument("--no-self-dimer-filter", action="store_true",
                        help="do not reject overlaps that pair with a copy of themselves")
    search.add_argument("--self-dimer-dg-threshold", type=float,
                        default=DEFAULT_SELF_DIMER_DG_THRESHOLD,
                        help="overlaps whose self dimer delta G (37 degC, kcal/mol, strongest of the "
                             "two strands) is more negative than this are rejected "
                             "(default: %(default)s)")
    search.add_argument("--no-cross-dimer-filter", action="store_true",
                        help="do not reject pairs of overlaps that can anneal to each other")
    search.add_argument("--cross-dimer-dg-threshold", type=float,
                        default=DEFAULT_CROSS_DIMER_DG_THRESHOLD,
                        help="pairs of overlaps whose cross dimer delta G (37 degC, kcal/mol) is more "
                             "negative than this are not used together (default: %(default)s)")

    primers = parser.add_argument_group("primers and conditions")
    primers.add_argument("--primer-extension-length", type=int,
                         default=spec.default_primer_extension_length,
                         help="length of the primer extension region (default: %(default)s)")
    primers.add_argument("--jobs", type=int, default=None,
                         help="worker processes for the cross-dimer thermodynamics; 1 runs "
                              "everything in one process (default: all cores this process may use, "
                              "at most %d)" % MAX_AUTO_JOBS)
    primers.add_argument("--na-conc", type=float, default=ThermoConditions.na_mM,
                         help="Na+ concentration in mM for Tm calculation (default: %(default)s)")
    primers.add_argument("--mg-conc", type=float, default=ThermoConditions.mg_mM,
                         help="Mg2+ concentration in mM (default: %(default)s)")
    primers.add_argument("--dntp-conc", type=float, default=ThermoConditions.dntp_mM,
                         help="dNTP concentration in mM (default: %(default)s)")
    primers.add_argument("--oligo-conc", type=float, default=ThermoConditions.oligo_uM,
                         help="oligo concentration in uM (default: %(default)s)")
    return parser


def run_cli(spec: TopologySpec, argv: Optional[Sequence[str]] = None) -> int:
    """Command-line entry point. Returns the process exit code."""
    parser = build_argument_parser(spec)
    args = parser.parse_args(argv)
    missing = [flag for flag, value in (("--input", args.input), ("--overlap-length", args.overlap_length),
                                         ("--num-modifications", args.num_modifications)) if value is None]
    if missing:
        parser.error("the following arguments are required: " + ", ".join(missing))
    if args.jobs is not None and args.jobs < 1:
        parser.error("--jobs must be at least 1")

    params = DesignParameters(
        overlap_length=args.overlap_length,
        num_modifications=args.num_modifications,
        search_method=args.search_method,
        hairpin_filter=not args.no_hairpin_filter,
        hairpin_dg_threshold=args.hairpin_dg_threshold,
        self_dimer_filter=not args.no_self_dimer_filter,
        self_dimer_dg_threshold=args.self_dimer_dg_threshold,
        cross_dimer_filter=not args.no_cross_dimer_filter,
        cross_dimer_dg_threshold=args.cross_dimer_dg_threshold,
        max_dg_spread=args.max_dg_spread,
        search_frame_length=args.search_frame_length,
        primer_extension_length=args.primer_extension_length,
        jobs=args.jobs,
        thermo=ThermoConditions(args.na_conc, args.mg_conc, args.dntp_conc, args.oligo_conc),
    )
    try:
        with open(args.input, "r") as handle:
            fragments = parse_fragment_sequences(handle.read())
        result = run_design(fragments, params, spec)
    except (DesignError, OSError) as error:
        print("Error: %s" % error, file=sys.stderr)
        return 1

    report = format_report(result)
    print(report)
    try:
        write_primer_csv(result, args.output_csv)
        print("Primer table written to %s" % args.output_csv)
        if args.report:
            with open(args.report, "w") as handle:
                handle.write(report + "\n")
            print("Report written to %s" % args.report)
    except OSError as error:
        print("Error: %s" % error, file=sys.stderr)
        return 1
    return 0


###############################################################
# GRAPHICAL USER INTERFACE (FreeSimpleGUI)
###############################################################

def _read_form(values: Dict[str, object]) -> Tuple[List[str], DesignParameters]:
    """Convert the values of the input window into fragments and design parameters."""
    try:
        expected_fragments = int(str(values["-NUM-"]).strip())
        overlap_length = int(str(values["-OVERLAP-"]).strip())
        num_modifications = int(str(values["-MODS-"]).strip())
    except ValueError:
        raise DesignError("The number of fragments, the overlap length and the number of "
                          "modifications must be integers.")
    fragments = parse_fragment_sequences(str(values["-SEQUENCES-"]))
    if len(fragments) != expected_fragments:
        raise DesignError("%d fragment sequences were entered, but %d fragments were expected."
                          % (len(fragments), expected_fragments))
    params = DesignParameters(
        overlap_length=overlap_length,
        num_modifications=num_modifications,
        search_method=str(values["-METHOD-"]),
        hairpin_filter=bool(values["-HAIRPIN-"]),
        self_dimer_filter=bool(values["-SELFDIMER-"]),
        cross_dimer_filter=bool(values["-CROSSDIMER-"]),
    )
    return fragments, params


def run_gui(spec: TopologySpec) -> None:
    """Open the graphical interface for the given topology."""
    try:
        import FreeSimpleGUI as sg
    except ImportError:
        print("FreeSimpleGUI is not installed (pip install FreeSimpleGUI); use the command-line "
              "interface instead (see --help).", file=sys.stderr)
        raise SystemExit(1)

    sg.theme("LightBlue")
    sg.theme_button_color(("blue", "#6D9F85"))
    font = ("Helvetica", 14)
    try:
        from Title_png import img as title_image
    except ImportError:
        title_image = None

    def banner():
        """A *fresh* banner row each time: an element must not be placed in two windows."""
        return [[sg.Image(data=title_image)]] if title_image else []

    layout = banner() + [
        [sg.Text("Number of fragments:", font=font, size=(36, 1)), sg.Input(key="-NUM-", font=font, size=(6, 1))],
        [sg.Text("Overlap length:", font=font, size=(36, 1)), sg.Input(key="-OVERLAP-", font=font, size=(6, 1))],
        [sg.Text("Number of modifications per primer:", font=font, size=(36, 1)),
         sg.Input(key="-MODS-", font=font, size=(6, 1))],
        [sg.Text("Overlap search method:", font=font, size=(36, 1)),
         sg.Combo(list(SEARCH_METHODS), default_value="exact", readonly=True, key="-METHOD-", font=font)],
        [sg.Checkbox("Filter overlaps by hairpin", default=True, key="-HAIRPIN-", font=font),
         sg.Checkbox("Filter overlaps by self dimer", default=True, key="-SELFDIMER-", font=font),
         sg.Checkbox("Filter overlap pairs by cross dimer", default=True, key="-CROSSDIMER-",
                     font=font)],
        [sg.Text("Fragment sequences, one per line, in assembly order:", font=font)],
        [sg.Multiline(key="-SEQUENCES-", font=("Courier", 11), size=(90, 10))],
        [sg.Button("Primer Design", font=font), sg.Button("Exit", font=font)],
    ]
    window = sg.Window("modZN for primer design (%s assembly)" % spec.name, layout,
                       element_justification="center")

    while True:
        event, values = window.read()
        if event in (sg.WIN_CLOSED, "Exit"):
            break
        if event != "Primer Design":
            continue
        try:
            fragments, params = _read_form(values)
            result = run_design(fragments, params, spec)
        except DesignError as error:
            sg.popup_error(str(error), title="modZN")
            continue
        except Exception as error:  # keep the window alive on unexpected errors
            traceback.print_exc()
            sg.popup_error("Unexpected error: %s" % error, title="modZN")
            continue
        _show_result_window(sg, banner, result, font)

    window.close()


def _show_result_window(sg, banner, result: DesignResult, font) -> None:
    """Show the report and let the user export the primer table.

    ``banner`` is a callable returning a fresh banner row, because a FreeSimpleGUI element may
    only ever be placed in one window; reusing one makes window creation fail.
    """
    layout = banner() + [
        [sg.Multiline(format_report(result), font=("Courier", 11), size=(90, 25), disabled=True)],
        [sg.Button("Export CSV", font=font), sg.Button("Close", font=font)],
    ]
    window = sg.Window("modZN results", layout, element_justification="center", modal=True)
    while True:
        event, _values = window.read()
        if event in (sg.WIN_CLOSED, "Close"):
            break
        if event == "Export CSV":
            path = sg.popup_get_file("Save primer table", save_as=True, default_path="data.csv",
                                     default_extension=".csv", file_types=(("CSV files", "*.csv"),))
            if path:
                try:
                    write_primer_csv(result, path)
                except OSError as error:       # unwritable location, full disk, ...
                    sg.popup_error("Could not save the primer table: %s" % error, title="modZN")
                else:
                    sg.popup("Saved to %s" % path, title="modZN")
    window.close()


def run_application(spec: TopologySpec, argv: Optional[Sequence[str]] = None) -> int:
    """Start the GUI when no arguments are given (or --gui is used), otherwise run the command line."""
    arguments = list(sys.argv[1:] if argv is None else argv)
    if not arguments or "--gui" in arguments:
        run_gui(spec)
        return 0
    return run_cli(spec, arguments)

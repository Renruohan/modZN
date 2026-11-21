"""
modZN for LINEAR assembly (sequence design with GUI and command-line interface)

Designs the overlapping zones and the PCR primers of a linear modAssembly.

    python modZN_linear.py            # graphical interface
    python modZN_linear.py --help     # command-line interface

Shared algorithms live in modZN_common.py; this file only adds what is specific
to linear assemblies (primers of the two terminal fragments).
"""

import sys
from typing import List, Sequence, Tuple

from modZN_common import (
    MODIFICATION_MARK,
    DesignParameters,
    TopologySpec,
    design_reverse_primer,
    design_standard_primer_pair,
    reverse_complement,
    run_application,
    run_design,
    spread_modifications,
)

# Primer extension length (bases); override with --primer-extension-length.
DEFAULT_PRIMER_EXTENSION_LENGTH = 25


def design_primer_pairs(fragments: Sequence[str], layout: Tuple[int, ...],
                        params: DesignParameters) -> List[Tuple[str, str]]:
    """(forward, reverse) primers for every fragment of a linear assembly.

    The first forward primer and the last reverse primer have no neighbour: they are plain
    ``overlap_length + extension`` bases of the fragment end, without modifications.
    """
    overlap_length = params.overlap_length
    extension = params.primer_extension_length
    marks = MODIFICATION_MARK * params.num_modifications
    last = len(fragments) - 1

    pairs = []
    for index, fragment in enumerate(fragments):
        if index == 0:
            forward = fragment[:overlap_length + extension]
            reverse = design_reverse_primer(fragments, layout, index, params)
        elif index == last:
            previous = fragments[index - 1]
            previous_left = layout[index - 1]
            previous_right = overlap_length - previous_left
            forward = spread_modifications(
                previous[len(previous) - previous_left:] + fragment[:previous_right] + marks
                + fragment[previous_right:previous_right + extension]
            )
            reverse = reverse_complement(fragment[len(fragment) - (overlap_length + extension):])
        else:
            forward, reverse = design_standard_primer_pair(fragments, layout, index, params)
        pairs.append((forward, reverse))
    return pairs


SPEC = TopologySpec(
    name="linear",
    circular=False,
    default_primer_extension_length=DEFAULT_PRIMER_EXTENSION_LENGTH,
    design_primer_pairs=design_primer_pairs,
    program_name="modZN_linear.py",
)


def design(fragments: Sequence[str], params: DesignParameters):
    """Design overlaps and primers for a linear assembly (see modZN_common.run_design)."""
    return run_design(fragments, params, SPEC)


def run_gui() -> None:
    """Open the graphical interface for linear assembly."""
    from modZN_common import run_gui as _run_gui
    _run_gui(SPEC)


def main(argv=None) -> int:
    """GUI when called without arguments, command-line interface otherwise."""
    return run_application(SPEC, argv)


if __name__ == "__main__":
    sys.exit(main())

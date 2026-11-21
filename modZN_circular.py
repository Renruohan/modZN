"""
modZN for CIRCULAR assembly (sequence design with GUI and command-line interface)

Designs the overlapping zones and the PCR primers of a circular modAssembly
(for example a plasmid).

    python modZN_circular.py            # graphical interface
    python modZN_circular.py --help     # command-line interface

Shared algorithms live in modZN_common.py.  In a circular assembly every fragment
has two neighbours (the last is followed by the first), so all are treated alike.
"""

import sys
from typing import List, Sequence, Tuple

from modZN_common import (
    DesignParameters,
    TopologySpec,
    design_standard_primer_pair,
    run_application,
    run_design,
)

# Primer extension length (bases); override with --primer-extension-length.
DEFAULT_PRIMER_EXTENSION_LENGTH = 25


def design_primer_pairs(fragments: Sequence[str], layout: Tuple[int, ...],
                        params: DesignParameters) -> List[Tuple[str, str]]:
    """(forward, reverse) primers for every fragment of a circular assembly."""
    return [design_standard_primer_pair(fragments, layout, index, params)
            for index in range(len(fragments))]


SPEC = TopologySpec(
    name="circular",
    circular=True,
    default_primer_extension_length=DEFAULT_PRIMER_EXTENSION_LENGTH,
    design_primer_pairs=design_primer_pairs,
    program_name="modZN_circular.py",
)


def design(fragments: Sequence[str], params: DesignParameters):
    """Design overlaps and primers for a circular assembly (see modZN_common.run_design)."""
    return run_design(fragments, params, SPEC)


def run_gui() -> None:
    """Open the graphical interface for circular assembly."""
    from modZN_common import run_gui as _run_gui
    _run_gui(SPEC)


def main(argv=None) -> int:
    """GUI when called without arguments, command-line interface otherwise."""
    return run_application(SPEC, argv)


if __name__ == "__main__":
    sys.exit(main())

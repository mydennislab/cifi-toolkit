"""
cifi - toolkit for downstream processing of CiFi long reads.

https://voles.dennislab.org
"""

__version__ = "1.0.0"

from ._core import (
    FilterResult,
    ProcessingResult,
    SingleEnzymeQCResult,
    Statistics,
    filter_bam,
    find_all_degenerate,
    get_enzyme_info,
    has_degenerate_bases,
    is_bam_file,
    list_enzymes,
    process_reads,
    process_reads_custom,
    revcomp,
    revcomp_degenerate,
    run_qc_analysis_custom,
)

__all__ = [
    "process_reads",
    "process_reads_custom",
    "filter_bam",
    "run_qc_analysis_custom",
    "list_enzymes",
    "get_enzyme_info",
    "is_bam_file",
    "find_all_degenerate",
    "has_degenerate_bases",
    "revcomp",
    "revcomp_degenerate",
    "ProcessingResult",
    "FilterResult",
    "Statistics",
    "SingleEnzymeQCResult",
    "__version__",
]

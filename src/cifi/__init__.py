"""
cifi - toolkit for downstream processing of CiFi long reads.

https://voles.dennislab.org
"""

__version__ = "1.0.0"

from ._core import (
    ContactsResult,
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
    pa5_position,
    parse_segment_name,
    process_reads,
    process_reads_custom,
    reconstruct_contacts,
    revcomp,
    revcomp_degenerate,
    run_qc_analysis_custom,
    segment_name,
)

__all__ = [
    "process_reads",
    "process_reads_custom",
    "filter_bam",
    "reconstruct_contacts",
    "run_qc_analysis_custom",
    "list_enzymes",
    "get_enzyme_info",
    "is_bam_file",
    "find_all_degenerate",
    "has_degenerate_bases",
    "revcomp",
    "revcomp_degenerate",
    "segment_name",
    "parse_segment_name",
    "pa5_position",
    "ProcessingResult",
    "FilterResult",
    "ContactsResult",
    "Statistics",
    "SingleEnzymeQCResult",
    "__version__",
]

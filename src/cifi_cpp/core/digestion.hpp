#pragma once

#include "enzyme.hpp"
#include "../stats/statistics.hpp"
#include "../io/writer.hpp"
#include "../io/table_writer.hpp"
#include <cstdint>
#include <string>
#include <vector>
#include <memory>

namespace cifi {

struct ProcessingConfig {
    EnzymeInfo enzyme;
    int min_segments = 3;
    // Minimum length of an *emitted* read, guaranteed for both R1 and R2.
    int min_segment_len = 60;
    // Drop the 5' site remnant from segments that begin at a cut.
    bool strip_overhang = true;
    // Reverse complement R2 (after any strip). Off by default: segments keep
    // the orientation they were sequenced in, as in the Pore-C convention.
    bool revcomp_r2 = false;
    bool fast_mode = false;
};

/**
 * What extract_segments kept, and what it discarded on the way.
 */
struct SegmentExtraction {
    std::vector<std::pair<size_t, size_t>> segments;
    // 1-based index of the cut-delimited span each kept segment came from,
    // parallel to segments. Dropped spans leave gaps; see segment_name.hpp.
    std::vector<uint32_t> span_index;
    // The cut span each kept segment came from, before the 5' trim; parallel
    // to segments. segments[i] lies inside spans[i].
    std::vector<std::pair<size_t, size_t>> spans;
    // Cut spans in the read, sites + 1, empty ones included: the range of
    // span_index.
    uint32_t spans_total = 0;
    uint64_t dropped_short = 0;   // non-empty spans below min_emit_len after trimming
    uint64_t bases_dropped = 0;   // bases in those spans, post-trim
    uint64_t bases_trimmed = 0;   // bases removed by the 5' site-remnant trim
};

struct ProcessingResult {
    uint64_t reads_in = 0;
    uint64_t reads_out = 0;
    uint64_t reads_skipped = 0;
    uint64_t pairs_written = 0;
    uint64_t total_segments = 0;

    // Input profile, over every read including ones later skipped
    uint64_t total_bases_in = 0;
    uint64_t gc_bases_in = 0;
    uint64_t total_sites = 0;

    // Yield, over reads that passed
    uint64_t bases_out_r1 = 0;
    uint64_t bases_out_r2 = 0;
    // Unique-segment output; both stay 0 when no segments writer is given
    uint64_t segments_written = 0;
    uint64_t bases_out_segments = 0;
    uint64_t segments_dropped_short = 0;
    uint64_t bases_dropped_short = 0;
    uint64_t bases_trimmed_overhang = 0;
    // Bases in reads that no segment survived from, so that
    // kept + trimmed + dropped + this == total_bases_in
    uint64_t bases_in_filtered_reads = 0;
    // Rows of the segments table; 0 when no table writer is given
    uint64_t segments_table_rows = 0;

    // Filtering reason counters
    uint64_t filtered_few_sites = 0;    // Reads with < min_segments sites
    uint64_t filtered_short_segments = 0;  // Reads where all segments were too short after length filtering

    Statistics segment_length_stats;
    Statistics sites_per_read_stats;
    // Input-read length, over every read; the others cover passing reads only.
    Statistics read_length_stats;
    Statistics segments_per_read_stats;
    Statistics pairs_per_read_stats;

    ProcessingResult(bool fast_mode = false)
        : segment_length_stats(fast_mode, 100)
        , sites_per_read_stats(fast_mode, 1)
        , read_length_stats(fast_mode, 100)
        , segments_per_read_stats(fast_mode, 1)
        , pairs_per_read_stats(fast_mode, 1) {}
};

/**
 * Header lines of the segments table (see segments_table_row): the
 * ##key=value metadata, then the #columns line.
 */
std::vector<std::string> segments_table_header(const std::string& command,
                                               const std::string& version);

/**
 * Process a single read: digest and write all pairwise contacts.
 * Returns true if read passed filters and was processed.
 *
 * When out_segments is given, each retained segment of a passing read is
 * also written once, in read order and native orientation, named per
 * segment_name.hpp. When out_table is given, the same segments get one row
 * each in the segments table, in the same order. The R1/R2 output is the
 * same whatever is given.
 */
bool process_single_read(
    const std::string& name,
    const std::string& sequence,
    const std::string& quality,
    const ProcessingConfig& config,
    FastqWriter& out_r1,
    FastqWriter& out_r2,
    ProcessingResult& result,
    FastqWriter* out_segments = nullptr,
    TableWriter* out_table = nullptr
);

/**
 * Extract the emitted span of each segment.
 *
 * lead_trim is removed from segments that begin at a cut site; the read's
 * leading segment does not begin at one and so is never trimmed. Spans shorter
 * than min_emit_len are dropped, which makes min_emit_len a guarantee about the
 * emitted read rather than about the untrimmed source segment.
 */
SegmentExtraction extract_segments(
    const std::string& sequence,
    const EnzymeInfo& enzyme,
    int min_emit_len,
    int lead_trim = 0
);

} // namespace cifi

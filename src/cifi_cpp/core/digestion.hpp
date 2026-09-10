#pragma once

#include "enzyme.hpp"
#include "../stats/statistics.hpp"
#include "../io/writer.hpp"
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
    uint64_t segments_dropped_short = 0;
    uint64_t bases_dropped_short = 0;
    uint64_t bases_trimmed_overhang = 0;
    // Bases in reads that no segment survived from, so that
    // kept + trimmed + dropped + this == total_bases_in
    uint64_t bases_in_filtered_reads = 0;

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
 * Process a single read: digest and write all pairwise contacts.
 * Returns true if read passed filters and was processed.
 */
bool process_single_read(
    const std::string& name,
    const std::string& sequence,
    const std::string& quality,
    const ProcessingConfig& config,
    FastqWriter& out_r1,
    FastqWriter& out_r2,
    ProcessingResult& result
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

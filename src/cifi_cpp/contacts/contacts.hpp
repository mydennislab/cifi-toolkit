#pragma once

#include "../stats/statistics.hpp"
#include <cstdint>
#include <string>

namespace cifi {

struct ContactsConfig {
    // Segments under this MAPQ take no part in any contact. 1 follows the
    // CiFi paper and the previous assembly pipeline; cifi filter's default of
    // 30 is tuned for contact maps rather than scaffolding.
    int min_mapq = 1;
    // BAM decompression threads
    int threads = 4;
};

struct ContactsResult {
    uint64_t records_seen = 0;           // alignment records in the input
    uint64_t segments_seen = 0;          // distinct segment QNAMEs
    uint64_t reads_seen = 0;             // distinct original reads (QNAME groups)
    uint64_t primary_mapped = 0;         // segments whose primary record is mapped
    uint64_t unmapped = 0;               // segments whose primary record is unmapped
    uint64_t secondary_ignored = 0;      // records
    uint64_t supplementary_ignored = 0;  // records
    uint64_t duplicate_primary = 0;      // further primary records of a segment already seen
    uint64_t below_mapq = 0;             // mapped primaries under min_mapq
    uint64_t usable_segments = 0;        // mapped primaries at or above min_mapq
    uint64_t reads_with_contacts = 0;    // reads with >= 2 usable segments
    uint64_t contacts_written = 0;
    // Records the R1/R2 route would have mapped for the same segments, i.e.
    // n(n-1) per read: both mates of every pair. Against records_seen this is
    // the mapping work saved by mapping unique segments.
    uint64_t pair_mates_equivalent = 0;
    uint64_t max_usable_in_read = 0;
    uint64_t max_contacts_in_read = 0;
    std::string sort_order;              // @HD SO of the input, empty when absent

    // Over contributing reads only. Bin size 1 keeps these exact while
    // holding one counter per distinct value instead of one per read.
    Statistics usable_per_read_stats{true, 1};
    Statistics contacts_per_read_stats{true, 1};
};

/**
 * PA5 position of one alignment.
 *
 * yahs reads a PA5 position verbatim (link.c, dump_links_from_pa5_file:
 * sscanf then a clamp to len-1), so the value has to be the same one yahs
 * would derive itself from a name-sorted BAM: the 0-based midpoint
 * s/2 + e/2 + (s&1 && e&1) with s = pos0 and e = the exclusive 0-based end,
 * which equals floor((pos0 + end0) / 2). Despite the "zero-based to
 * one-based" comment next to that formula, no +1 is applied there, and none
 * is applied here. end0 comes from bam_endpos (M/=/X/D/N consume reference).
 */
inline uint32_t pa5_position(int64_t pos0, int64_t end0) {
    return static_cast<uint32_t>(pos0 + (end0 - pos0) / 2);
}

/**
 * Stream a name-grouped BAM of uniquely emitted segments and write every
 * pairwise contact between the usable segments of each read as a PA5 row.
 *
 * Memory is bounded by the segments of the read currently streaming past;
 * the quadratic expansion happens once per read, at the group boundary.
 */
ContactsResult reconstruct_contacts(
    const std::string& input_path,
    const std::string& output_path,
    const ContactsConfig& config
);

} // namespace cifi

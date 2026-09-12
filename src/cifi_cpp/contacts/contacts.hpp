#pragma once

#include "../stats/statistics.hpp"
#include <cstdint>
#include <string>

namespace cifi {

/**
 * How a contact is written; both forms carry the same pairs and MAPQs.
 *
 * PA5 has one row per contact with a point position per segment. yahs
 * rebuilds an interval around that point from its global --read-length
 * (link.c, dump_links_from_pa5_file: rl >>= 1 at line 1726, then
 * [p - rl, p + rl] at 1762-1765), which suits fixed-length Hi-C reads but
 * not CiFi segments, whose lengths vary by orders of magnitude.
 *
 * BED has two consecutive rows per contact, each with the segment's own
 * aligned span. yahs takes those spans as they are for coverage
 * (link.c 1614-1621) and derives from them the midpoint PA5 would have
 * carried (link.c 1633-1634). No read length enters into it, so BED is the
 * form meant for scaffolding CiFi data.
 */
enum class ContactsFormat { PA5, BED };

struct ContactsConfig {
    // Segments under this MAPQ take no part in any contact. 1 follows the
    // CiFi paper and the previous assembly pipeline; cifi filter's default of
    // 30 is tuned for contact maps rather than scaffolding.
    int min_mapq = 1;
    // BAM decompression threads
    int threads = 4;
    ContactsFormat format = ContactsFormat::PA5;
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
    std::string group_order;             // @HD GO, "query" when grouped but unsorted

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
 * pairwise contact between the usable segments of each read, as one PA5 row
 * or two BED rows per contact according to config.format.
 *
 * Memory is bounded by the segments of the read currently streaming past;
 * the quadratic expansion happens once per read, at the group boundary.
 * The output takes its final name only once it is complete; a run that
 * fails leaves no partial file behind.
 */
ContactsResult reconstruct_contacts(
    const std::string& input_path,
    const std::string& output_path,
    const ContactsConfig& config
);

} // namespace cifi

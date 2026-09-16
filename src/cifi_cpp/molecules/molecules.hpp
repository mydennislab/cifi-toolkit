#pragma once

#include <cstdint>
#include <string>

namespace cifi {

/**
 * The long molecule table: `cifi molecules`.
 *
 * A CiFi read is one molecule of ordered segments. The digest names each
 * retained segment <read>__CIFI_SEG__<span_index> and, with
 * --segments-table, records where in the read it sat. Once the segments
 * are aligned, this table puts the molecule back together: one row per
 * (segment, alignment record), molecules in BAM order, segments in read
 * order, records of a segment by rank. Anything derived from it is never a
 * second source of truth.
 *
 * Two indices, never overloaded: span_index is the digest's cut-span
 * index (gaps where spans were dropped, the k of the segment name);
 * retained_index is the dense 1..N order of the molecule's segments.
 */
struct MoleculesConfig {
    // The segments table of the digest, for spans_total, read_length and
    // the read coordinates; empty writes '.' in those columns. When given,
    // it is held in memory (about 55 bytes per segment plus a per-molecule
    // index: 86 MB resident for 1.1 M segments in 200,000 molecules) and
    // the molecule and segment sets of the two inputs must agree, see
    // extract_molecules.
    std::string table;
    // Drop secondary and supplementary rows (--candidates primary)
    bool primary_only = false;
    // BAM decompression threads
    int threads = 4;
    // Recorded in the header
    std::string command;
    std::string version;
};

struct MoleculesResult {
    uint64_t records_seen = 0;
    uint64_t molecules = 0;                 // distinct reads (QNAME groups)
    uint64_t segments = 0;                  // distinct segments, table ones included
    uint64_t rows_written = 0;
    uint64_t primary_mapped = 0;            // records: mapped, neither secondary nor supplementary
    uint64_t unmapped = 0;                  // records
    uint64_t secondary = 0;                 // records
    uint64_t supplementary = 0;             // records
    uint64_t duplicate_primary = 0;         // a further primary record of a segment already seen
    uint64_t segments_with_candidates = 0;  // segments with more than one record
    uint64_t molecules_with_candidates = 0; // molecules with such a segment
    uint64_t rows_dropped_candidates = 0;   // secondary/supplementary rows left out (primary_only)
    uint64_t table_molecules = 0;           // loaded from the segments table
    uint64_t table_segments = 0;
    uint64_t table_segments_missing_from_bam = 0;   // written as U rows
    uint64_t table_molecules_missing_from_bam = 0;  // no row at all
    std::string sort_order;                 // @HD SO of the input, empty when absent
    std::string group_order;                // @HD GO, "query" when grouped but unsorted
};

/**
 * Stream a name-grouped BAM of segment alignments (secondary, supplementary
 * and unmapped records included) and write the long molecule table as a
 * bgzip TSV. One molecule is in memory at a time.
 *
 * With a segments table: a BAM segment absent from the table, a BAM
 * molecule absent from the table, or a segment whose length differs
 * between the two is an error naming the segment. A table segment absent
 * from the BAM while its molecule is present is written as an unmapped
 * row so the molecule stays complete, and counted; a table molecule with
 * no record at all is only counted.
 *
 * The output takes its final name only once complete.
 */
MoleculesResult extract_molecules(
    const std::string& input_path,
    const std::string& output_path,
    const MoleculesConfig& config
);

} // namespace cifi

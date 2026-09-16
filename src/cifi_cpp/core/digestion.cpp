#include "digestion.hpp"
#include "segment_name.hpp"
#include <algorithm>
#include <stdexcept>

namespace cifi {

SegmentExtraction extract_segments(
    const std::string& sequence,
    const EnzymeInfo& enzyme,
    int min_emit_len,
    int lead_trim
) {
    auto sites = find_all_degenerate(sequence, enzyme.site);

    // Build cut positions
    std::vector<size_t> cuts;
    cuts.push_back(0);
    for (size_t pos : sites) {
        cuts.push_back(pos + enzyme.cut_offset);
    }
    cuts.push_back(sequence.length());

    // Take each segment's emitted span, then filter on that length. Only a
    // segment that begins at a cut carries the site remnant; the read's
    // leading segment starts at position 0 and so keeps its full length.
    SegmentExtraction out;
    out.spans_total = static_cast<uint32_t>(cuts.size() - 1);
    for (size_t i = 0; i < cuts.size() - 1; i++) {
        size_t start = cuts[i];
        size_t end = cuts[i + 1];
        size_t span_start = start;

        // Span 0 runs from the read start and carries no remnant; every later
        // span begins at a cut, including one that falls at offset 0.
        if (i > 0) {
            // clamp so a span shorter than the remnant cannot underflow
            size_t trimmed = std::min(start + static_cast<size_t>(lead_trim), end);
            out.bases_trimmed += trimmed - start;
            start = trimmed;
        }

        if (end <= start) {
            continue;  // empty span (adjacent cuts, or a site at either end)
        }
        if (static_cast<int>(end - start) >= min_emit_len) {
            out.segments.push_back({start, end});
            out.span_index.push_back(static_cast<uint32_t>(i + 1));
            out.spans.push_back({span_start, end});
        } else {
            out.dropped_short++;
            out.bases_dropped += end - start;
        }
    }

    return out;
}

static const char SEGMENTS_TABLE_COLUMNS[] =
    "molecule_id\tspan_index\tretained_index\tsegments_kept\tspans_total\t"
    "read_length\tread_start\tread_end\tspan_start\tspan_end\toriginal_len\t"
    "processed_len\ttrimmed_5p\tterminal\tenzyme\tcut_offset";

// The header is ##key=value metadata, one item per line, then the #columns
// line; the vocabulary is fixed (README, "Molecule tables"). Nothing in it
// varies between two runs of the same command.
std::vector<std::string> segments_table_header(const std::string& command,
                                               const std::string& version) {
    return {
        "##cifi_format=segments",
        "##format_version=1",
        "##tool_version=" + version,
        "##command=" + command,
        "##read_coordinates=0-based-half-open",
        "##span_index=1-based-original-digest-span",
        "##retained_index=1-based-dense-retained-order",
        "##terminal=T5,T3,I,T5T3",
        std::string("#columns: ") + SEGMENTS_TABLE_COLUMNS,
    };
}

// One table row per retained segment. Two indices, never confused: span_index
// is the cut span the segment came from (gaps where spans were dropped, the
// same k as in the FASTQ name), retained_index its rank among the segments
// kept (dense 1..N). read_length is the whole read, so the read coordinates
// can be turned into PAF query coordinates of the molecule. terminal marks
// the read ends by span: T5 for the first span, T3 for the last, I in
// between; a read of one span is T5T3.
static void write_segments_table_row(TableWriter& out, const std::string& name,
                                     size_t read_length,
                                     const SegmentExtraction& extraction, size_t i,
                                     const EnzymeInfo& enzyme, std::string& line) {
    const auto& [start, end] = extraction.segments[i];
    const auto& [span_start, span_end] = extraction.spans[i];
    uint32_t span_index = extraction.span_index[i];
    const char* terminal = "I";
    if (span_index == 1 && span_index == extraction.spans_total) terminal = "T5T3";
    else if (span_index == 1) terminal = "T5";
    else if (span_index == extraction.spans_total) terminal = "T3";

    line.clear();
    line += name;
    line += '\t'; line += std::to_string(span_index);
    line += '\t'; line += std::to_string(i + 1);
    line += '\t'; line += std::to_string(extraction.segments.size());
    line += '\t'; line += std::to_string(extraction.spans_total);
    line += '\t'; line += std::to_string(read_length);
    line += '\t'; line += std::to_string(start);
    line += '\t'; line += std::to_string(end);
    line += '\t'; line += std::to_string(span_start);
    line += '\t'; line += std::to_string(span_end);
    line += '\t'; line += std::to_string(span_end - span_start);
    line += '\t'; line += std::to_string(end - start);
    line += '\t'; line += std::to_string(start - span_start);
    line += '\t'; line += terminal;
    // a named enzyme by its name, a custom site by the site itself
    line += '\t'; line += enzyme.name == "Custom" ? enzyme.site : enzyme.name;
    line += '\t'; line += std::to_string(enzyme.cut_offset);
    line += '\n';
    out.write(line);
}

bool process_single_read(
    const std::string& name,
    const std::string& sequence,
    const std::string& quality,
    const ProcessingConfig& config,
    FastqWriter& out_r1,
    FastqWriter& out_r2,
    ProcessingResult& result,
    FastqWriter* out_segments,
    TableWriter* out_table
) {
    // Input profile: recorded for every read, including ones skipped below,
    // so the report describes what was fed in rather than what survived.
    result.read_length_stats.add(static_cast<int>(sequence.length()));
    result.total_bases_in += sequence.length();
    for (char c : sequence) {
        if (c == 'G' || c == 'C' || c == 'g' || c == 'c') result.gc_bases_in++;
    }

    // Find sites (supports degenerate IUPAC bases)
    auto sites = find_all_degenerate(sequence, config.enzyme.site);
    result.total_sites += sites.size();

    // Early exit: not enough sites
    if (static_cast<int>(sites.size()) < config.min_segments - 1) {
        result.reads_skipped++;
        result.filtered_few_sites++;
        result.bases_in_filtered_reads += sequence.length();
        return false;
    }

    // The site remnant belongs to the segment, not to whichever slot the
    // segment lands in, so trim once here and let both mates read the same
    // spans. min_segment_len then bounds the emitted read directly.
    int lead_trim = config.strip_overhang ? config.enzyme.overhang_length() : 0;

    auto extraction = extract_segments(sequence, config.enzyme,
                                       config.min_segment_len, lead_trim);
    const auto& segments = extraction.segments;

    // Record what extraction discarded before deciding the read's fate, so the
    // totals still add up for reads that are filtered out below.
    result.segments_dropped_short += extraction.dropped_short;
    result.bases_dropped_short += extraction.bases_dropped;
    result.bases_trimmed_overhang += extraction.bases_trimmed;

    // Check segment count after length filtering
    if (static_cast<int>(segments.size()) < config.min_segments) {
        result.reads_skipped++;
        result.filtered_short_segments++;
        for (const auto& [s0, e0] : segments) {
            result.bases_in_filtered_reads += e0 - s0;
        }
        return false;
    }

    // Record stats for passing reads only
    result.sites_per_read_stats.add(static_cast<int>(sites.size()));
    for (const auto& [start, end] : segments) {
        result.segment_length_stats.add(static_cast<int>(end - start));
    }

    result.reads_out++;
    result.total_segments += segments.size();
    result.segments_per_read_stats.add(static_cast<int>(segments.size()));
    result.pairs_per_read_stats.add(
        static_cast<int>(segments.size() * (segments.size() - 1) / 2));

    // The parser finds the marker from the right, so a name that already
    // carries it would still split, but its segments could no longer be
    // told apart from those of a read named after them. The table keys on
    // the same names, so it refuses such reads as well.
    if ((out_segments || out_table) &&
        name.find(SEGMENT_NAME_SEP) != std::string::npos) {
        throw std::runtime_error("read name already contains " +
                                 std::string(SEGMENT_NAME_SEP) + ": " + name);
    }

    // Each retained segment once, as it lies in the read: this is what gets
    // mapped when contacts are reconstructed after alignment, so the R2
    // reverse complement (a paired-FASTQ concern) does not apply here.
    if (out_segments) {
        for (size_t i = 0; i < segments.size(); i++) {
            const auto& [start, end] = segments[i];
            out_segments->write(segment_name(name, extraction.span_index[i]),
                                sequence.substr(start, end - start),
                                quality.substr(start, end - start));
            result.bases_out_segments += end - start;
        }
        result.segments_written += segments.size();
    }

    // The same segments, one table row each, in the same order.
    if (out_table) {
        std::string line;
        for (size_t i = 0; i < segments.size(); i++) {
            write_segments_table_row(*out_table, name, sequence.length(), extraction, i,
                                     config.enzyme, line);
        }
        result.segments_table_rows += segments.size();
    }

    // Generate ALL pairs (n choose 2)
    for (size_t i = 0; i < segments.size(); i++) {
        for (size_t j = i + 1; j < segments.size(); j++) {
            const auto& f1 = segments[i];
            const auto& f2 = segments[j];

            std::string seq1 = sequence.substr(f1.first, f1.second - f1.first);
            std::string qual1 = quality.substr(f1.first, f1.second - f1.first);
            std::string seq2 = sequence.substr(f2.first, f2.second - f2.first);
            std::string qual2 = quality.substr(f2.first, f2.second - f2.first);

            // Segments arrive already trimmed, so R2 differs from R1 only by
            // the optional reverse complement.
            std::string r2_seq, r2_qual;
            if (config.revcomp_r2) {
                r2_seq = revcomp(seq2);
                r2_qual.assign(qual2.rbegin(), qual2.rend());
            } else {
                r2_seq = std::move(seq2);
                r2_qual = std::move(qual2);
            }

            // Both mates carry one name: R1/R2 FASTQ has no flags, so the name
            // is the only thing that identifies a pair. A "/1" or "/2" suffix
            // would make the two files disagree on every pair. These indices
            // count kept segments (unlike the span indices in the unique
            // segment names) and are kept as they were for compatibility.
            std::string pair_name = name + "_" + std::to_string(i) + "_" + std::to_string(j - i - 1);

            out_r1.write(pair_name, seq1, qual1);
            out_r2.write(pair_name, r2_seq, r2_qual);
            result.bases_out_r1 += seq1.length();
            result.bases_out_r2 += r2_seq.length();

            result.pairs_written++;
        }
    }

    return true;
}

} // namespace cifi

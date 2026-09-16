#include "contacts.hpp"
#include "../core/segment_name.hpp"
#include "../io/grouped_bam.hpp"
#include "../io/hts_handles.hpp"
#include "../io/writer.hpp"

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <parallel_hashmap/phmap.h>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

namespace cifi {

namespace {

struct UsableSegment {
    uint32_t index;   // span index from the QNAME
    int32_t tid;
    uint32_t start;   // 0-based reference start of the primary alignment
    uint32_t end;     // exclusive end from the CIGAR (bam_endpos)
    uint8_t mapq;
};

// Everything known about the read whose segments are currently streaming
// past. It is reset at each group boundary, which is what bounds memory by
// the segments of one read rather than by the BAM.
struct ReadGroup {
    std::string name;
    bool active = false;
    std::vector<UsableSegment> usable;
    phmap::flat_hash_set<uint32_t> seen;       // any record of the segment
    phmap::flat_hash_set<uint32_t> primaries;  // a primary record was taken

    void start(const std::string& read) {
        name = read;
        active = true;
        usable.clear();
        seen.clear();
        primaries.clear();
    }
};

void append_uint(std::string& s, uint64_t v) {
    char tmp[24];
    int n = 0;
    do {
        tmp[n++] = static_cast<char>('0' + v % 10);
        v /= 10;
    } while (v);
    while (n) s += tmp[--n];
}

// PA5 row: pair_name contig1 pos1 contig2 pos2 mapq1 mapq2, one per contact.
void append_pa5(std::string& line, const std::string& name,
                const char* contig_a, const UsableSegment& a,
                const char* contig_b, const UsableSegment& b) {
    line += name;
    line += '\t';
    line += contig_a;
    line += '\t';
    append_uint(line, pa5_position(a.start, a.end));
    line += '\t';
    line += contig_b;
    line += '\t';
    append_uint(line, pa5_position(b.start, b.end));
    line += '\t';
    append_uint(line, a.mapq);
    line += '\t';
    append_uint(line, b.mapq);
    line += '\n';
}

// BED record: contig start end pair_name mapq, two per contact, back to back.
//
// This is the shape yahs's BED reader expects (link.c, dump_links_from_bed_file):
// it holds one record and pairs it with the next only if the names match
// (lines 1581-1592; is_read_pair in asset.c 175-188 accepts identical names,
// so no /1 /2 suffix is needed), otherwise the earlier record is dropped
// and the later one becomes the held record (1652-1660). Columns are read as
// "%s %u %u %s %hhu" (1585, 1591). The link position is the midpoint
// s/2 + e/2 + (s&1 && e&1) of the two integers (1633-1634), the same value
// the PA5 row carries, and [start, end] goes into the coverage track as
// given (1614-1621); a PA5 input has only the point, so there the interval
// is manufactured from --read-length (1726, 1762-1765).
void append_bed(std::string& line, const char* contig, const UsableSegment& s,
                const std::string& name) {
    line += contig;
    line += '\t';
    append_uint(line, s.start);
    line += '\t';
    append_uint(line, s.end);
    line += '\t';
    line += name;
    line += '\t';
    append_uint(line, s.mapq);
    line += '\n';
}

void flush_group(ReadGroup& group, const sam_hdr_t* hdr, TextWriter& out,
                 ContactsFormat format, ContactsResult& result) {
    if (!group.active) return;
    uint64_t n = group.seen.size();
    result.pair_mates_equivalent += n * (n - 1);

    size_t k = group.usable.size();
    if (k < 2) return;

    // Segment order in the input depends on how it was sorted (natural or
    // lexicographic names); ordering by span index makes the output the same
    // either way.
    std::sort(group.usable.begin(), group.usable.end(),
              [](const UsableSegment& a, const UsableSegment& b) { return a.index < b.index; });

    std::string name, line;
    for (size_t i = 0; i < k; i++) {
        const auto& a = group.usable[i];
        const char* contig_a = sam_hdr_tid2name(hdr, a.tid);
        for (size_t j = i + 1; j < k; j++) {
            const auto& b = group.usable[j];
            const char* contig_b = sam_hdr_tid2name(hdr, b.tid);
            name = contact_name(group.name, a.index, b.index);
            line.clear();
            if (format == ContactsFormat::BED) {
                append_bed(line, contig_a, a, name);
                append_bed(line, contig_b, b, name);
            } else {
                append_pa5(line, name, contig_a, a, contig_b, b);
            }
            out.write(line);
        }
    }

    uint64_t contacts = static_cast<uint64_t>(k) * (k - 1) / 2;
    result.reads_with_contacts++;
    result.contacts_written += contacts;
    result.max_usable_in_read = std::max<uint64_t>(result.max_usable_in_read, k);
    result.max_contacts_in_read = std::max(result.max_contacts_in_read, contacts);
    result.usable_per_read_stats.add(static_cast<int>(k));
    result.contacts_per_read_stats.add(static_cast<int>(contacts));
}

}  // namespace

ContactsResult reconstruct_contacts(
    const std::string& input_path,
    const std::string& output_path,
    const ContactsConfig& config
) {
    HtsFilePtr fp(hts_open(input_path.c_str(), "r"));
    if (!fp) throw std::runtime_error("Cannot open: " + input_path);
    if (config.threads > 1) hts_set_threads(fp.get(), config.threads);

    SamHeaderPtr hdr(sam_hdr_read(fp.get()));
    if (!hdr) throw std::runtime_error("Cannot read header: " + input_path);

    ContactsResult result;

    // Only a name-grouped input keeps a read's segments together; refuse the
    // one order that is known to scatter them rather than stream through it
    // and emit a fraction of the contacts. The caller needs both tags to
    // tell a grouped input from an undeclared one.
    read_group_order(hdr.get(), input_path, "cifi contacts", result.sort_order,
                     result.group_order);

    TextWriter out(output_path);
    ReadGroup group;
    RecentReads recent(256);
    BamRecordPtr b(bam_init1());

    int ret;
    while ((ret = sam_read1(fp.get(), hdr.get(), b.get())) >= 0) {
        result.records_seen++;
        auto seg = parse_segment_name(bam_get_qname(b.get()));

        if (!group.active || seg.read != group.name) {
            flush_group(group, hdr.get(), out, config.format, result);
            if (group.active) recent.push(group.name);
            if (recent.contains(seg.read)) {
                throw std::runtime_error(
                    "segments of read " + seg.read + " are not grouped together in " +
                    input_path + "; sort by name first: samtools sort -n");
            }
            group.start(seg.read);
            result.reads_seen++;
        }
        if (group.seen.insert(seg.index).second) result.segments_seen++;

        uint16_t flag = b->core.flag;
        if (flag & BAM_FSECONDARY) {
            result.secondary_ignored++;
            continue;
        }
        if (flag & BAM_FSUPPLEMENTARY) {
            result.supplementary_ignored++;
            continue;
        }
        // One primary record per segment. A second one (concatenated BAMs,
        // say) is counted as the anomaly it is; the first one stands.
        if (!group.primaries.insert(seg.index).second) {
            result.duplicate_primary++;
            continue;
        }
        if ((flag & BAM_FUNMAP) || b->core.tid < 0) {
            result.unmapped++;
            continue;
        }
        result.primary_mapped++;
        if (static_cast<int>(b->core.qual) < config.min_mapq) {
            result.below_mapq++;
            continue;
        }
        result.usable_segments++;
        group.usable.push_back({seg.index, b->core.tid,
                                static_cast<uint32_t>(b->core.pos),
                                static_cast<uint32_t>(bam_endpos(b.get())),
                                b->core.qual});
    }
    if (ret < -1) throw std::runtime_error("Truncated or corrupt input: " + input_path);

    flush_group(group, hdr.get(), out, config.format, result);
    out.close();
    return result;
}

}  // namespace cifi

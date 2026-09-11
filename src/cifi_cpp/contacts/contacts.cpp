#include "contacts.hpp"
#include "../core/segment_name.hpp"
#include "../io/writer.hpp"

#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/sam.h>
#include <parallel_hashmap/phmap.h>

#include <algorithm>
#include <deque>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace cifi {

namespace {

struct HtsFileCloser {
    void operator()(htsFile* fp) const { if (fp) hts_close(fp); }
};
struct HeaderDestroyer {
    void operator()(sam_hdr_t* h) const { if (h) sam_hdr_destroy(h); }
};
struct RecordDestroyer {
    void operator()(bam1_t* b) const { if (b) bam_destroy1(b); }
};

struct UsableSegment {
    uint32_t index;   // span index from the QNAME
    int32_t tid;
    uint32_t pos;     // PA5 position
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

// The last few finished reads. A BAM that is grouped only in stretches, or a
// coordinate-sorted one missing its SO tag, brings a read back after its
// group was already flushed; that would silently lose contacts, so it is an
// error. Fixed size: a sanity net, not an index of the file.
class RecentReads {
public:
    explicit RecentReads(size_t capacity) : capacity_(capacity) {}

    bool contains(const std::string& name) const { return set_.count(name) > 0; }

    void push(const std::string& name) {
        order_.push_back(name);
        set_.insert(name);
        if (order_.size() > capacity_) {
            set_.erase(order_.front());
            order_.pop_front();
        }
    }

private:
    size_t capacity_;
    std::deque<std::string> order_;
    phmap::flat_hash_set<std::string> set_;
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

void flush_group(ReadGroup& group, const sam_hdr_t* hdr, TextWriter& out,
                 ContactsResult& result) {
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

    std::string line;
    for (size_t i = 0; i < k; i++) {
        const auto& a = group.usable[i];
        const char* contig_a = sam_hdr_tid2name(hdr, a.tid);
        for (size_t j = i + 1; j < k; j++) {
            const auto& b = group.usable[j];
            line = contact_name(group.name, a.index, b.index);
            line += '\t';
            line += contig_a;
            line += '\t';
            append_uint(line, a.pos);
            line += '\t';
            line += sam_hdr_tid2name(hdr, b.tid);
            line += '\t';
            append_uint(line, b.pos);
            line += '\t';
            append_uint(line, a.mapq);
            line += '\t';
            append_uint(line, b.mapq);
            line += '\n';
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
    std::unique_ptr<htsFile, HtsFileCloser> fp(hts_open(input_path.c_str(), "r"));
    if (!fp) throw std::runtime_error("Cannot open: " + input_path);
    if (config.threads > 1) hts_set_threads(fp.get(), config.threads);

    std::unique_ptr<sam_hdr_t, HeaderDestroyer> hdr(sam_hdr_read(fp.get()));
    if (!hdr) throw std::runtime_error("Cannot read header: " + input_path);

    ContactsResult result;

    // Only a name-grouped input keeps a read's segments together; refuse the
    // one order that is known to scatter them rather than stream through it
    // and emit a fraction of the contacts.
    kstring_t so = KS_INITIALIZE;
    if (sam_hdr_find_tag_hd(hdr.get(), "SO", &so) == 0 && ks_len(&so) > 0) {
        result.sort_order.assign(ks_str(&so), ks_len(&so));
    }
    ks_free(&so);
    if (result.sort_order == "coordinate") {
        throw std::runtime_error(
            input_path + " is sorted by coordinate; cifi contacts needs the segments "
            "of each read grouped together. Sort by name first: "
            "samtools sort -n -o segments.ns.bam " + input_path);
    }

    TextWriter out(output_path);
    ReadGroup group;
    RecentReads recent(256);
    std::unique_ptr<bam1_t, RecordDestroyer> b(bam_init1());

    int ret;
    while ((ret = sam_read1(fp.get(), hdr.get(), b.get())) >= 0) {
        result.records_seen++;
        auto seg = parse_segment_name(bam_get_qname(b.get()));

        if (!group.active || seg.read != group.name) {
            flush_group(group, hdr.get(), out, result);
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
                                pa5_position(b->core.pos, bam_endpos(b.get())),
                                b->core.qual});
    }
    if (ret < -1) throw std::runtime_error("Truncated or corrupt input: " + input_path);

    flush_group(group, hdr.get(), out, result);
    out.close();
    return result;
}

}  // namespace cifi

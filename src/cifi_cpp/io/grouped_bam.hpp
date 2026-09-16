#pragma once

#include <htslib/kstring.h>
#include <htslib/sam.h>
#include <parallel_hashmap/phmap.h>

#include <deque>
#include <stdexcept>
#include <string>

namespace cifi {

/**
 * What `cifi contacts` and `cifi molecules` need from a BAM grouped by
 * read name: the header's declaration of its order, and a net that catches
 * an input that is not grouped after all.
 */

/**
 * The @HD SO and GO tags, empty when absent. minimap2 declares its own
 * grouping as SO:unsorted GO:query; samtools sort -n as SO:queryname. Only
 * coordinate order is known to scatter a read's segments, so that one is
 * refused here, naming the tool; anything else is streamed.
 */
inline void read_group_order(sam_hdr_t* hdr, const std::string& input_path,
                             const char* tool, std::string& sort_order,
                             std::string& group_order) {
    kstring_t tag = KS_INITIALIZE;
    if (sam_hdr_find_tag_hd(hdr, "SO", &tag) == 0 && ks_len(&tag) > 0) {
        sort_order.assign(ks_str(&tag), ks_len(&tag));
    }
    if (sam_hdr_find_tag_hd(hdr, "GO", &tag) == 0 && ks_len(&tag) > 0) {
        group_order.assign(ks_str(&tag), ks_len(&tag));
    }
    ks_free(&tag);
    if (sort_order == "coordinate") {
        throw std::runtime_error(
            input_path + " is sorted by coordinate; " + tool + " needs the segments "
            "of each read grouped together. Sort by name first: "
            "samtools sort -n -o segments.ns.bam " + input_path);
    }
}

/**
 * The last few finished reads. A BAM that is grouped only in stretches, or
 * a coordinate-sorted one missing its SO tag, brings a read back after its
 * group was already flushed; that would silently lose part of it, so it is
 * an error. Fixed size: a sanity net, not an index of the file.
 */
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

} // namespace cifi

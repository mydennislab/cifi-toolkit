// References:
// - T-Digest: https://arxiv.org/abs/1902.04023
// - GK sketch: "Space-Efficient Online Computation of Quantile Summaries" (SIGMOD 2001)
#pragma once

#include <vector>
#include <utility>
#include <cstdint>
#include <climits>
#include <parallel_hashmap/phmap.h>

namespace cifi {

/**
 * Unified statistics collector with exact and fast modes.
 *
 * - Exact mode: stores all values, computes precise median/percentiles
 * - Fast mode: streaming stats + binned histogram, O(1) memory
 */
class Statistics {
public:
    explicit Statistics(bool fast_mode = false, int bin_size = 100);

    void add(int value);
    void clear();

    // Accessors (always available)
    uint64_t count() const { return count_; }
    uint64_t sum() const { return sum_; }
    int min() const { return min_; }
    int max() const { return max_; }
    double mean() const;

    // May be approximate in fast mode
    double median() const;
    double percentile(double p) const;

    // Histogram access (for reports)
    std::vector<std::pair<int, uint64_t>> get_histogram() const;

    /**
     * Equal-width histogram over [min, max] as (bin edges, counts).
     *
     * Built here rather than from values() so a report does not have to copy
     * every recorded value across the language boundary. Integer-valued data
     * should pass num_bins <= max-min+1 to avoid structurally empty bins.
     */
    std::pair<std::vector<double>, std::vector<uint64_t>> binned(
        int num_bins, bool integer_bins = false) const;

    // Raw values (exact mode only, empty in fast mode)
    const std::vector<int>& values() const { return values_; }

    bool is_fast_mode() const { return fast_mode_; }

private:
    bool fast_mode_;
    int bin_size_;

    // Running stats (O(1) memory)
    uint64_t count_ = 0;
    uint64_t sum_ = 0;
    int min_ = INT32_MAX;
    int max_ = INT32_MIN;

    // Exact mode: raw values
    std::vector<int> values_;

    // Fast mode: binned histogram
    phmap::flat_hash_map<int, uint64_t> histogram_;

    // Helper for approximate percentile from histogram
    double percentile_from_histogram(double p) const;
};

} // namespace cifi

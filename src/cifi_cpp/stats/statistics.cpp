#include "statistics.hpp"
#include <algorithm>
#include <cmath>

namespace cifi {

Statistics::Statistics(bool fast_mode, int bin_size)
    : fast_mode_(fast_mode), bin_size_(bin_size) {}

void Statistics::add(int value) {
    count_++;
    sum_ += value;
    if (value < min_) min_ = value;
    if (value > max_) max_ = value;

    if (fast_mode_) {
        histogram_[value / bin_size_]++;
    } else {
        values_.push_back(value);
    }
}

void Statistics::clear() {
    count_ = 0;
    sum_ = 0;
    min_ = INT32_MAX;
    max_ = INT32_MIN;
    values_.clear();
    histogram_.clear();
}

double Statistics::mean() const {
    return count_ > 0 ? static_cast<double>(sum_) / count_ : 0.0;
}

double Statistics::median() const {
    return percentile(0.5);
}

double Statistics::percentile(double p) const {
    if (count_ == 0) return 0.0;

    // Normalize: accept both 0-1 and 0-100 ranges
    if (p > 1.0) p /= 100.0;
    if (p < 0.0) p = 0.0;
    if (p > 1.0) p = 1.0;

    if (fast_mode_) {
        return percentile_from_histogram(p);
    }

    // Exact mode: sort and index
    std::vector<int> sorted = values_;
    std::sort(sorted.begin(), sorted.end());

    size_t idx = static_cast<size_t>(p * (sorted.size() - 1));
    return sorted[idx];
}

double Statistics::percentile_from_histogram(double p) const {
    if (histogram_.empty()) return 0.0;

    // Get sorted bins
    std::vector<std::pair<int, uint64_t>> bins(histogram_.begin(), histogram_.end());
    std::sort(bins.begin(), bins.end());

    // Find target count
    uint64_t target = static_cast<uint64_t>(p * count_);
    uint64_t cumulative = 0;

    for (const auto& [bin, cnt] : bins) {
        cumulative += cnt;
        if (cumulative >= target) {
            double bin_start = bin * bin_size_;
            // With unit bins the bin IS the value; interpolating would report a
            // median of value + 0.5, which can land above max().
            if (bin_size_ == 1) return bin_start;
            return bin_start + bin_size_ / 2.0;
        }
    }

    if (bin_size_ == 1) return bins.back().first;
    return bins.back().first * bin_size_ + bin_size_ / 2.0;
}

std::vector<std::pair<int, uint64_t>> Statistics::get_histogram() const {
    if (fast_mode_) {
        std::vector<std::pair<int, uint64_t>> result(histogram_.begin(), histogram_.end());
        std::sort(result.begin(), result.end());
        return result;
    }

    // Build histogram from raw values
    phmap::flat_hash_map<int, uint64_t> hist;
    for (int v : values_) {
        hist[v / bin_size_]++;
    }

    std::vector<std::pair<int, uint64_t>> result(hist.begin(), hist.end());
    std::sort(result.begin(), result.end());
    return result;
}


std::pair<std::vector<double>, std::vector<uint64_t>> Statistics::binned(
    int num_bins, bool integer_bins) const {
    std::vector<double> edges;
    std::vector<uint64_t> counts;
    if (count_ == 0 || num_bins < 1) return {edges, counts};

    if (min_ == max_) {
        edges.push_back(min_);
        counts.push_back(count_);
        return {edges, counts};
    }

    // Integer data is binned on half-integer edges so each bar is centred on
    // the value it counts, instead of drifting up to a unit away from it.
    const double lo = integer_bins ? min_ - 0.5 : min_;
    const double hi = integer_bins ? max_ + 0.5 : max_;
    const double width = (hi - lo) / num_bins;
    edges.reserve(num_bins + 1);
    for (int i = 0; i <= num_bins; ++i) edges.push_back(lo + i * width);
    counts.assign(num_bins, 0);

    auto place = [&](double value, uint64_t weight) {
        int idx = static_cast<int>((value - lo) / width);
        if (idx < 0) idx = 0;
        if (idx >= num_bins) idx = num_bins - 1;
        counts[idx] += weight;
    };

    if (fast_mode_) {
        // Bin midpoints stand in for the values they represent.
        for (const auto& [bin, cnt] : histogram_) {
            double centre = bin_size_ == 1 ? static_cast<double>(bin)
                                           : bin * bin_size_ + bin_size_ / 2.0;
            place(centre, cnt);
        }
    } else {
        for (int v : values_) place(v, 1);
    }

    return {edges, counts};
}

} // namespace cifi

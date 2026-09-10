// src/cifi_cpp/core/enzyme.cpp
#include "enzyme.hpp"
#include <unordered_map>
#include <algorithm>

namespace cifi {

static const std::unordered_map<std::string, std::pair<std::string, int>> ENZYMES = {
    // 4-cutters
    {"NlaIII", {"CATG", 4}},
    {"DpnII",  {"GATC", 0}},
    {"Sau3AI", {"GATC", 0}},
    {"MboI",   {"GATC", 0}},
    // 6-cutters
    {"HindIII", {"AAGCTT", 1}},
};

int EnzymeInfo::overhang_length() const {
    int site_len = static_cast<int>(site.length());
    if (cut_offset == 0) return site_len;
    // A cut_offset past the site would give a negative remnant, which becomes a
    // huge value once widened to size_t at the call site.
    if (cut_offset >= site_len) return 0;
    return site_len - cut_offset;
}

std::optional<EnzymeInfo> get_enzyme(const std::string& name) {
    auto it = ENZYMES.find(name);
    if (it == ENZYMES.end()) return std::nullopt;

    return EnzymeInfo{name, it->second.first, it->second.second};
}

std::vector<std::string> list_enzymes() {
    std::vector<std::string> names;
    names.reserve(ENZYMES.size());
    for (const auto& [name, _] : ENZYMES) {
        names.push_back(name);
    }
    std::sort(names.begin(), names.end());
    return names;
}

std::vector<size_t> find_all(const std::string& text, const std::string& pattern) {
    std::vector<size_t> positions;
    size_t pos = 0;
    while ((pos = text.find(pattern, pos)) != std::string::npos) {
        positions.push_back(pos);
        pos++;
    }
    return positions;
}

// Check if a base matches an IUPAC degenerate code
static bool base_matches(char base, char pattern_char) {
    // Convert to uppercase for comparison
    char b = (base >= 'a' && base <= 'z') ? base - 32 : base;
    char p = (pattern_char >= 'a' && pattern_char <= 'z') ? pattern_char - 32 : pattern_char;

    // Exact match
    if (b == p) return true;

    // IUPAC degenerate base codes
    switch (p) {
        case 'N': return (b == 'A' || b == 'C' || b == 'G' || b == 'T');
        case 'R': return (b == 'A' || b == 'G');  // puRine
        case 'Y': return (b == 'C' || b == 'T');  // pYrimidine
        case 'W': return (b == 'A' || b == 'T');  // Weak
        case 'S': return (b == 'C' || b == 'G');  // Strong
        case 'M': return (b == 'A' || b == 'C');  // aMino
        case 'K': return (b == 'G' || b == 'T');  // Keto
        case 'B': return (b == 'C' || b == 'G' || b == 'T');  // not A
        case 'D': return (b == 'A' || b == 'G' || b == 'T');  // not C
        case 'H': return (b == 'A' || b == 'C' || b == 'T');  // not G
        case 'V': return (b == 'A' || b == 'C' || b == 'G');  // not T
        default: return false;
    }
}

bool has_degenerate_bases(const std::string& pattern) {
    for (char c : pattern) {
        char upper = (c >= 'a' && c <= 'z') ? c - 32 : c;
        if (upper == 'N' || upper == 'R' || upper == 'Y' || upper == 'W' ||
            upper == 'S' || upper == 'M' || upper == 'K' || upper == 'B' ||
            upper == 'D' || upper == 'H' || upper == 'V') {
            return true;
        }
    }
    return false;
}

std::vector<size_t> find_all_degenerate(const std::string& text, const std::string& pattern) {
    std::vector<size_t> positions;

    if (pattern.empty() || text.size() < pattern.size()) {
        return positions;
    }

    // If no degenerate bases, use faster exact matching
    if (!has_degenerate_bases(pattern)) {
        return find_all(text, pattern);
    }

    // Degenerate pattern matching
    size_t text_len = text.size();
    size_t pattern_len = pattern.size();

    for (size_t i = 0; i <= text_len - pattern_len; ++i) {
        bool match = true;
        for (size_t j = 0; j < pattern_len; ++j) {
            if (!base_matches(text[i + j], pattern[j])) {
                match = false;
                break;
            }
        }
        if (match) {
            positions.push_back(i);
        }
    }

    return positions;
}

std::string revcomp(const std::string& seq) {
    std::string rc;
    rc.reserve(seq.size());
    for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
        switch (*it) {
            case 'A': case 'a': rc += 'T'; break;
            case 'T': case 't': rc += 'A'; break;
            case 'G': case 'g': rc += 'C'; break;
            case 'C': case 'c': rc += 'G'; break;
            default: rc += 'N'; break;
        }
    }
    return rc;
}

std::string revcomp_degenerate(const std::string& seq) {
    std::string rc;
    rc.reserve(seq.size());
    for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
        char c = *it;
        char upper = (c >= 'a' && c <= 'z') ? c - 32 : c;
        switch (upper) {
            case 'A': rc += 'T'; break;
            case 'T': rc += 'A'; break;
            case 'G': rc += 'C'; break;
            case 'C': rc += 'G'; break;
            case 'N': rc += 'N'; break;
            case 'R': rc += 'Y'; break;  // A/G -> T/C
            case 'Y': rc += 'R'; break;  // C/T -> G/A
            case 'W': rc += 'W'; break;  // A/T -> T/A (self-complement)
            case 'S': rc += 'S'; break;  // C/G -> G/C (self-complement)
            case 'M': rc += 'K'; break;  // A/C -> T/G
            case 'K': rc += 'M'; break;  // G/T -> C/A
            case 'B': rc += 'V'; break;  // C/G/T -> G/C/A
            case 'D': rc += 'H'; break;  // A/G/T -> T/C/A
            case 'H': rc += 'D'; break;  // A/C/T -> T/G/A
            case 'V': rc += 'B'; break;  // A/C/G -> T/G/C
            default: rc += 'N'; break;
        }
    }
    return rc;
}

} // namespace cifi

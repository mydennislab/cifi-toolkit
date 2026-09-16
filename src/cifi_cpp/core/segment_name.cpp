#include "segment_name.hpp"

#include <cstring>
#include <stdexcept>

namespace cifi {

static const size_t SEP_LEN = std::strlen(SEGMENT_NAME_SEP);

std::string segment_name(const std::string& read_name, uint32_t span_index) {
    std::string out;
    out.reserve(read_name.size() + SEP_LEN + 10);
    out += read_name;
    out += SEGMENT_NAME_SEP;
    out += std::to_string(span_index);
    return out;
}

ParsedSegmentName parse_segment_name(const std::string& qname) {
    size_t at = qname.rfind(SEGMENT_NAME_SEP);
    if (at == std::string::npos || at == 0) {
        throw std::invalid_argument(
            "not a cifi segment name (expected <read>" + std::string(SEGMENT_NAME_SEP) +
            "<index>, as written by cifi digest --segments-out): " + qname);
    }

    size_t digits = at + SEP_LEN;
    if (digits == qname.size()) {
        throw std::invalid_argument("segment name has no index: " + qname);
    }
    uint64_t index = 0;
    for (size_t i = digits; i < qname.size(); i++) {
        char c = qname[i];
        if (c < '0' || c > '9') {
            throw std::invalid_argument("segment index is not a number: " + qname);
        }
        index = index * 10 + static_cast<uint64_t>(c - '0');
        if (index > UINT32_MAX) {
            throw std::invalid_argument("segment index out of range: " + qname);
        }
    }
    // Spans count from 1; a 0 means the producer disagrees with the contract.
    if (index == 0) {
        throw std::invalid_argument("segment index must be >= 1: " + qname);
    }

    return {qname.substr(0, at), static_cast<uint32_t>(index)};
}

std::string contact_name(const std::string& read_name, uint32_t i, uint32_t j) {
    std::string out = segment_name(read_name, i);
    out += SEGMENT_NAME_SEP;
    out += std::to_string(j);
    return out;
}

} // namespace cifi

#pragma once

#include <cstdint>
#include <string>

namespace cifi {

/**
 * Naming contract for uniquely emitted segments.
 *
 * `cifi digest --segments-out` names each segment
 *
 *     <original_read_name>__CIFI_SEG__<k>
 *
 * where k is the 1-based index of the cut-delimited span the segment came
 * from. Spans that were dropped (too short, or empty) leave gaps rather than
 * renumbering their neighbours, so k is stable whatever the filters did.
 *
 * The suffix is found from the right, so nothing is assumed about the
 * original name: PacBio names carry '/', ':' and '_' in varying numbers. The
 * marker must not occur in the original name; digest refuses such reads.
 *
 * This scheme is separate from the R1/R2 pair names, which keep their
 * kept-order form "<read>_<i>_<j-i-1>" for backward compatibility. Only the
 * segment scheme is parsed downstream (by `cifi contacts`).
 */
constexpr const char SEGMENT_NAME_SEP[] = "__CIFI_SEG__";

std::string segment_name(const std::string& read_name, uint32_t span_index);

struct ParsedSegmentName {
    std::string read;
    uint32_t index;
};

/**
 * Split a segment QNAME back into read name and span index.
 *
 * Throws std::invalid_argument for anything that does not follow the
 * contract, so a BAM of something other than digested segments fails
 * loudly instead of collapsing into one enormous group.
 */
ParsedSegmentName parse_segment_name(const std::string& qname);

/**
 * Contact pair name (PA5 column 1, BED column 4): the two segment names
 * joined, sharing the read prefix, e.g. "<read>__CIFI_SEG__1__CIFI_SEG__3".
 * Unique per unordered pair since i < j is enforced by the caller.
 */
std::string contact_name(const std::string& read_name, uint32_t i, uint32_t j);

} // namespace cifi

#pragma once

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <memory>

namespace cifi {

/**
 * Owning pointers for the htslib objects a BAM reader holds.
 *
 * Reading loops throw on a bad record name or a failed write; with raw
 * handles the open file, header and record would leak on that path. These
 * release them during unwinding, so the readers need no cleanup code before
 * a throw.
 */
struct HtsFileCloser {
    void operator()(htsFile* fp) const { if (fp) hts_close(fp); }
};
struct HeaderDestroyer {
    void operator()(sam_hdr_t* h) const { if (h) sam_hdr_destroy(h); }
};
struct RecordDestroyer {
    void operator()(bam1_t* b) const { if (b) bam_destroy1(b); }
};

using HtsFilePtr = std::unique_ptr<htsFile, HtsFileCloser>;
using SamHeaderPtr = std::unique_ptr<sam_hdr_t, HeaderDestroyer>;
using BamRecordPtr = std::unique_ptr<bam1_t, RecordDestroyer>;

} // namespace cifi

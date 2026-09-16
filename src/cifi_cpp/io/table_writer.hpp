#pragma once

#include <htslib/bgzf.h>
#include <string>
#include <vector>

namespace cifi {

/**
 * Writer for the bgzip TSV tables (segments, modifications, molecules).
 *
 * Same contract as TextWriter: lines are collected into a block before they
 * reach the file, the file is written under a temporary name beside the
 * requested one and takes that name only in close(), and a run that throws
 * leaves nothing behind. Compression is BGZF whatever the name, so the
 * tables read with any gzip reader and can also be indexed and read
 * block-wise by htslib.
 *
 * The header (the `#` lines) is written first through header(), or handed
 * to close() when it holds something only the end of the run knows (the
 * modification modes seen). In that case the data went to the temporary
 * without a header, and close() puts the header in front of it: BGZF is a
 * run of independent gzip members, so the header becomes its own block and
 * the data blocks are copied behind it byte for byte, without
 * recompression. The EOF marker of the data file is kept, as it must end
 * the file.
 */
class TableWriter {
public:
    explicit TableWriter(const std::string& path);
    ~TableWriter();

    // The header lines, without newlines, before any write().
    void header(const std::vector<std::string>& lines);
    void write(const std::string& line);  // line carries its own newline
    // Flush, then move into place. The overload takes the header of a table
    // whose header() was never called.
    void close();
    void close(const std::vector<std::string>& lines);

private:
    void flush();
    void discard();
    void finish_data();
    void prepend_header(const std::vector<std::string>& lines);

    std::string path_;
    std::string tmp_path_;   // empty once the file is in place or discarded
    BGZF* bgzf_ = nullptr;
    std::string buf_;
    bool header_written_ = false;
};

} // namespace cifi

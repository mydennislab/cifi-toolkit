#pragma once

#include <string>
#include <memory>
#include <fstream>
#include <zlib.h>

namespace cifi {

/**
 * Abstract FASTQ writer interface.
 */
class FastqWriter {
public:
    virtual ~FastqWriter() = default;
    virtual void write(const std::string& name,
                       const std::string& seq,
                       const std::string& qual) = 0;
    virtual void close() = 0;
};

/**
 * Plain text FASTQ writer.
 */
class PlainFastqWriter : public FastqWriter {
public:
    explicit PlainFastqWriter(const std::string& path);
    ~PlainFastqWriter() override;

    void write(const std::string& name,
               const std::string& seq,
               const std::string& qual) override;
    void close() override;

private:
    std::ofstream out_;
    std::string path_;
    std::string buf_;  // reused record buffer
};

/**
 * Gzip-compressed FASTQ writer.
 */
class GzipFastqWriter : public FastqWriter {
public:
    explicit GzipFastqWriter(const std::string& path);
    ~GzipFastqWriter() override;

    void write(const std::string& name,
               const std::string& seq,
               const std::string& qual) override;
    void close() override;

private:
    gzFile gz_;
    std::string buf_;  // reused record buffer, avoids a malloc per read
};

/**
 * Factory: creates appropriate writer based on path and flags.
 */
std::unique_ptr<FastqWriter> make_writer(const std::string& path, bool force_gzip);

/**
 * Helper: check if path ends with .gz
 */
bool ends_with_gz(const std::string& path);

} // namespace cifi

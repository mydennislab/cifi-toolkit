#include "writer.hpp"
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <algorithm>
#include <limits>
#include <random>
#include <fcntl.h>
#include <unistd.h>

namespace cifi {

// PlainFastqWriter

PlainFastqWriter::PlainFastqWriter(const std::string& path)
    : out_(path), path_(path) {
    if (!out_) {
        throw std::runtime_error("Cannot open for writing: " + path);
    }
}

PlainFastqWriter::~PlainFastqWriter() {
    // close() reports flush failures by throwing, which would terminate if it
    // escaped a destructor. An explicit close() is the path that surfaces them.
    try {
        close();
    } catch (...) {
    }
}

void PlainFastqWriter::write(const std::string& name,
                              const std::string& seq,
                              const std::string& qual) {
    buf_.clear();
    buf_.reserve(name.size() + seq.size() + qual.size() + 6);
    buf_ += '@';
    buf_ += name;
    buf_ += '\n';
    buf_ += seq;
    buf_ += "\n+\n";
    buf_ += qual;
    buf_ += '\n';

    out_.write(buf_.data(), static_cast<std::streamsize>(buf_.size()));
    // A failed write here would otherwise leave R1 and R2 with different record
    // counts and still exit successfully.
    if (!out_) {
        throw std::runtime_error("Failed writing FASTQ record " + name +
                                 " to " + path_);
    }
}

void PlainFastqWriter::close() {
    if (out_.is_open()) {
        // The final flush happens here; a failure now truncates the file.
        out_.close();
        if (out_.fail()) {
            throw std::runtime_error("Failed closing " + path_);
        }
    }
}

// GzipFastqWriter

GzipFastqWriter::GzipFastqWriter(const std::string& path) {
    gz_ = gzopen(path.c_str(), "wb");
    if (!gz_) {
        throw std::runtime_error("Cannot open for gzip writing: " + path);
    }
}

GzipFastqWriter::~GzipFastqWriter() {
    // close() reports flush failures by throwing, which would terminate if it
    // escaped a destructor. An explicit close() is the path that surfaces them.
    try {
        close();
    } catch (...) {
    }
}

void GzipFastqWriter::write(const std::string& name,
                             const std::string& seq,
                             const std::string& qual) {
    // gzprintf() formats into a fixed internal buffer (8KB by default) and
    // silently writes nothing when the record does not fit. CiFi segments run
    // to tens of kb, so build the record ourselves and hand it to gzwrite,
    // which has no length limit.
    buf_.clear();
    buf_.reserve(name.size() + seq.size() + qual.size() + 6);
    buf_ += '@';
    buf_ += name;
    buf_ += '\n';
    buf_ += seq;
    buf_ += "\n+\n";
    buf_ += qual;
    buf_ += '\n';

    if (buf_.size() > static_cast<size_t>(std::numeric_limits<int>::max())) {
        throw std::runtime_error("FASTQ record too large to write: " + name);
    }

    int written = gzwrite(gz_, buf_.data(), static_cast<unsigned>(buf_.size()));
    if (written != static_cast<int>(buf_.size())) {
        int err = 0;
        const char* msg = gzerror(gz_, &err);
        throw std::runtime_error("Failed writing FASTQ record " + name + ": " +
                                 (msg ? msg : "short write"));
    }
}

void GzipFastqWriter::close() {
    if (gz_) {
        // gzclose performs the final deflate flush, so a failure here can
        // truncate an otherwise complete file.
        int rc = gzclose(gz_);
        gz_ = nullptr;
        if (rc != Z_OK) {
            throw std::runtime_error("Failed closing gzip output (zlib code " +
                                     std::to_string(rc) + ")");
        }
    }
}

// Factory

bool ends_with_gz(const std::string& path) {
    if (path.size() < 3) return false;
    std::string suffix = path.substr(path.size() - 3);
    std::transform(suffix.begin(), suffix.end(), suffix.begin(), ::tolower);
    return suffix == ".gz";
}

std::unique_ptr<FastqWriter> make_writer(const std::string& path, bool force_gzip) {
    bool use_gzip = force_gzip || ends_with_gz(path);

    if (use_gzip) {
        std::string gz_path = ends_with_gz(path) ? path : path + ".gz";
        return std::make_unique<GzipFastqWriter>(gz_path);
    }
    return std::make_unique<PlainFastqWriter>(path);
}

// TextWriter

static const size_t TEXT_WRITER_BLOCK = 1 << 20;

// A temporary name beside path that this run alone holds. A fixed name such
// as path + ".tmp" would let two runs asked for the same output truncate one
// file and interleave in it, and the survivor's rename would promote the
// mixture; O_EXCL cannot hand the same name to both. The descriptor is
// closed again because the writer opens the name through the stream or
// zlib; the file stays claimed either way, with the usual permissions.
static std::string claim_temp_name(const std::string& path) {
    static const char alnum[] =
        "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
    std::random_device rd;
    std::uniform_int_distribution<int> pick(0, sizeof(alnum) - 2);
    for (int attempt = 0; attempt < 100; attempt++) {
        std::string candidate = path + ".tmp.";
        for (int i = 0; i < 6; i++) candidate += alnum[pick(rd)];
        int fd = open(candidate.c_str(), O_WRONLY | O_CREAT | O_EXCL, 0666);
        if (fd >= 0) {
            close(fd);
            return candidate;
        }
        if (errno != EEXIST) {
            throw std::runtime_error("Cannot open for writing: " + path + ": " +
                                     std::strerror(errno));
        }
    }
    throw std::runtime_error("No free temporary name beside " + path);
}

TextWriter::TextWriter(const std::string& path)
    : path_(path), tmp_path_(claim_temp_name(path)) {
    // A sibling of the final path: rename() only replaces atomically within
    // one filesystem. Compression follows the final name, as before. The
    // name is already claimed, so a failure to open it must drop the claim:
    // no destructor runs for a constructor that throws.
    if (ends_with_gz(path)) {
        gz_ = gzopen(tmp_path_.c_str(), "wb");
        if (!gz_) {
            discard();
            throw std::runtime_error("Cannot open for gzip writing: " + path);
        }
    } else {
        out_.open(tmp_path_);
        if (!out_) {
            discard();
            throw std::runtime_error("Cannot open for writing: " + path);
        }
    }
    buf_.reserve(TEXT_WRITER_BLOCK + 4096);
}

TextWriter::~TextWriter() {
    // Reaching the destructor with the temporary still present means the
    // caller never got to close(): an exception is unwinding through it.
    // The partial file must not be promoted, so it is dropped here rather
    // than completed. Nothing here throws.
    discard();
}

void TextWriter::write(const std::string& line) {
    buf_ += line;
    if (buf_.size() >= TEXT_WRITER_BLOCK) {
        flush();
    }
}

void TextWriter::flush() {
    if (buf_.empty()) return;
    if (gz_) {
        int written = gzwrite(gz_, buf_.data(), static_cast<unsigned>(buf_.size()));
        if (written != static_cast<int>(buf_.size())) {
            int err = 0;
            const char* msg = gzerror(gz_, &err);
            throw std::runtime_error("Failed writing " + path_ + ": " +
                                     (msg ? msg : "short write"));
        }
    } else {
        out_.write(buf_.data(), static_cast<std::streamsize>(buf_.size()));
        if (!out_) {
            throw std::runtime_error("Failed writing " + path_);
        }
    }
    buf_.clear();
}

void TextWriter::discard() {
    if (gz_) {
        gzclose(gz_);
        gz_ = nullptr;
    }
    if (out_.is_open()) out_.close();
    if (!tmp_path_.empty()) {
        std::remove(tmp_path_.c_str());
        tmp_path_.clear();
    }
}

void TextWriter::close() {
    if (tmp_path_.empty()) return;  // already in place, or discarded
    try {
        if (gz_) {
            flush();
            int rc = gzclose(gz_);
            gz_ = nullptr;
            if (rc != Z_OK) {
                throw std::runtime_error("Failed closing " + path_ + " (zlib code " +
                                         std::to_string(rc) + ")");
            }
        } else if (out_.is_open()) {
            flush();
            out_.close();
            if (out_.fail()) {
                throw std::runtime_error("Failed closing " + path_);
            }
        }
        // The final flush succeeded; only now does the file take its name.
        if (std::rename(tmp_path_.c_str(), path_.c_str()) != 0) {
            throw std::runtime_error("Cannot move " + tmp_path_ + " to " + path_ + ": " +
                                     std::strerror(errno));
        }
        tmp_path_.clear();
    } catch (...) {
        discard();
        throw;
    }
}

} // namespace cifi

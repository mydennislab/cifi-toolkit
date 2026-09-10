#include "writer.hpp"
#include <stdexcept>
#include <algorithm>
#include <limits>

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

} // namespace cifi

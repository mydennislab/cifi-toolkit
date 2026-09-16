#include "table_writer.hpp"
#include "writer.hpp"

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <memory>
#include <stdexcept>

namespace cifi {

static const size_t TABLE_WRITER_BLOCK = 1 << 20;

TableWriter::TableWriter(const std::string& path)
    : path_(path), tmp_path_(claim_temp_name(path)) {
    // The name is already claimed, so a failure to open it must drop the
    // claim: no destructor runs for a constructor that throws.
    bgzf_ = bgzf_open(tmp_path_.c_str(), "w");
    if (!bgzf_) {
        discard();
        throw std::runtime_error("Cannot open for bgzip writing: " + path);
    }
    buf_.reserve(TABLE_WRITER_BLOCK + 4096);
}

TableWriter::~TableWriter() {
    // Reaching the destructor with the temporary still present means the
    // caller never got to close(): an exception is unwinding through it.
    // The partial file must not be promoted. Nothing here throws.
    discard();
}

void TableWriter::header(const std::vector<std::string>& lines) {
    for (const auto& line : lines) {
        buf_ += line;
        buf_ += '\n';
    }
    header_written_ = true;
}

void TableWriter::write(const std::string& line) {
    buf_ += line;
    if (buf_.size() >= TABLE_WRITER_BLOCK) {
        flush();
    }
}

void TableWriter::flush() {
    if (buf_.empty()) return;
    ssize_t written = bgzf_write(bgzf_, buf_.data(), buf_.size());
    if (written < 0 || static_cast<size_t>(written) != buf_.size()) {
        throw std::runtime_error("Failed writing " + path_);
    }
    buf_.clear();
}

void TableWriter::discard() {
    if (bgzf_) {
        bgzf_close(bgzf_);
        bgzf_ = nullptr;
    }
    if (!tmp_path_.empty()) {
        std::remove(tmp_path_.c_str());
        tmp_path_.clear();
    }
}

// The last flush and the EOF block; a failure here truncates the file.
void TableWriter::finish_data() {
    flush();
    int rc = bgzf_close(bgzf_);
    bgzf_ = nullptr;
    if (rc != 0) {
        throw std::runtime_error("Failed closing " + path_);
    }
}

void TableWriter::close() {
    if (tmp_path_.empty()) return;  // already in place, or discarded
    if (!header_written_) {
        throw std::logic_error("TableWriter::close() before header(): " + path_);
    }
    try {
        finish_data();
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

void TableWriter::close(const std::vector<std::string>& lines) {
    if (tmp_path_.empty()) return;
    if (header_written_) {
        throw std::logic_error("TableWriter: header given twice for " + path_);
    }
    try {
        finish_data();
        prepend_header(lines);
        tmp_path_.clear();
    } catch (...) {
        discard();
        throw;
    }
}

// The header as BGZF blocks of its own, then the data file behind them.
// Written to a second temporary, so the requested name still appears only
// once the whole file is there.
void TableWriter::prepend_header(const std::vector<std::string>& lines) {
    std::string text;
    for (const auto& line : lines) {
        text += line;
        text += '\n';
    }

    std::string out_path = claim_temp_name(path_);
    try {
        std::ofstream out(out_path, std::ios::binary);
        if (!out) throw std::runtime_error("Cannot open for writing: " + path_);

        std::unique_ptr<char[]> block(new char[BGZF_MAX_BLOCK_SIZE]);
        for (size_t at = 0; at < text.size(); at += BGZF_BLOCK_SIZE) {
            size_t take = std::min<size_t>(BGZF_BLOCK_SIZE, text.size() - at);
            size_t dlen = BGZF_MAX_BLOCK_SIZE;
            if (bgzf_compress(block.get(), &dlen, text.data() + at, take, -1) != 0) {
                throw std::runtime_error("Failed compressing the header of " + path_);
            }
            out.write(block.get(), static_cast<std::streamsize>(dlen));
        }

        std::ifstream in(tmp_path_, std::ios::binary);
        if (!in) throw std::runtime_error("Cannot reopen " + tmp_path_);
        std::vector<char> copy(1 << 20);
        while (in) {
            in.read(copy.data(), static_cast<std::streamsize>(copy.size()));
            out.write(copy.data(), in.gcount());
        }
        if (in.bad() || !out) throw std::runtime_error("Failed writing " + path_);
        out.close();
        if (out.fail()) throw std::runtime_error("Failed closing " + path_);

        if (std::rename(out_path.c_str(), path_.c_str()) != 0) {
            throw std::runtime_error("Cannot move " + out_path + " to " + path_ + ": " +
                                     std::strerror(errno));
        }
    } catch (...) {
        std::remove(out_path.c_str());
        throw;
    }
    std::remove(tmp_path_.c_str());
}

} // namespace cifi

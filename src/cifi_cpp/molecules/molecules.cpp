#include "molecules.hpp"
#include "../core/segment_name.hpp"
#include "../io/grouped_bam.hpp"
#include "../io/hts_handles.hpp"
#include "../io/table_writer.hpp"

#include <htslib/bgzf.h>
#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/sam.h>
#include <parallel_hashmap/phmap.h>

#include <algorithm>
#include <climits>
#include <stdexcept>
#include <string>
#include <vector>

namespace cifi {

namespace {

const char MOLECULES_COLUMNS[] =
    "molecule_id\tspan_index\tretained_index\tsegment_count\tspans_total\t"
    "read_length\tread_start\tread_end\tseg_len\tref\tref_start\tref_end\tstrand\t"
    "mapq\taln\trank\tas\tnm\tqstart\tqend\tmlen\tblen\tcigar";

// One alignment record of one segment, reduced to what the row carries.
struct Candidate {
    uint32_t span_index;
    char aln;              // P primary, S secondary, L supplementary, U unmapped
    bool first_primary;    // the record that takes rank 0
    int32_t tid;
    int64_t ref_start;
    int64_t ref_end;
    bool reverse;
    int mapq;
    bool has_as, has_nm;
    int64_t as, nm;
    int64_t ambiguous;     // minimap2's nn tag: ambiguous bases in the alignment, else 0
    uint64_t seg_len;      // whole segment, clips included (one value per segment)
    uint64_t qstart, qend; // aligned query interval, PAF convention
    uint64_t mlen, blen;   // matching bases, alignment block length
    std::string cigar;
    uint64_t order;        // arrival order, the last tie-break
};

// The read whose records are currently streaming past; reset at each
// group boundary, which is what bounds memory by one molecule.
struct Molecule {
    std::string name;
    bool active = false;
    std::vector<Candidate> records;
    phmap::flat_hash_set<uint32_t> primaries;   // spans whose rank-0 record was taken

    void start(const std::string& read) {
        name = read;
        active = true;
        records.clear();
        primaries.clear();
    }
};

// --- the segments table, held in memory ------------------------------------

struct TableSegment {
    uint32_t span_index;
    uint32_t spans_total;
    uint64_t read_start;
    uint64_t read_end;
};

struct TableMolecule {
    uint32_t first;   // into segments
    uint32_t count;
    uint64_t read_length;
    bool seen;
};

class SegmentsTable {
public:
    bool loaded() const { return loaded_; }

    void load(const std::string& path, MoleculesResult& result) {
        BGZF* fp = bgzf_open(path.c_str(), "r");
        if (!fp) throw std::runtime_error("Cannot open segments table: " + path);
        kstring_t line = KS_INITIALIZE;
        std::vector<std::string> fields;
        int col_molecule = -1, col_span = -1, col_total = -1, col_start = -1, col_end = -1;
        int col_read_length = -1;
        try {
            while (bgzf_getline(fp, '\n', &line) >= 0) {
                if (line.l == 0) continue;
                if (line.s[0] == '#') {
                    const char prefix[] = "#columns:";
                    if (line.l > sizeof(prefix) - 1 &&
                        std::equal(prefix, prefix + sizeof(prefix) - 1, line.s)) {
                        size_t at = sizeof(prefix) - 1;
                        while (at < line.l && line.s[at] == ' ') at++;
                        split(std::string(line.s + at, line.l - at), fields);
                        for (size_t i = 0; i < fields.size(); i++) {
                            if (fields[i] == "molecule_id") col_molecule = static_cast<int>(i);
                            else if (fields[i] == "span_index") col_span = static_cast<int>(i);
                            else if (fields[i] == "spans_total") col_total = static_cast<int>(i);
                            else if (fields[i] == "read_length") col_read_length = static_cast<int>(i);
                            else if (fields[i] == "read_start") col_start = static_cast<int>(i);
                            else if (fields[i] == "read_end") col_end = static_cast<int>(i);
                        }
                    }
                    continue;
                }
                if (col_molecule < 0 || col_span < 0 || col_total < 0 || col_read_length < 0 ||
                    col_start < 0 || col_end < 0) {
                    throw std::runtime_error(
                        path + " has no #columns line naming molecule_id, span_index, "
                        "spans_total, read_length, read_start and read_end; is it a cifi "
                        "segments table?");
                }
                split(std::string(line.s, line.l), fields);
                if (static_cast<int>(fields.size()) <= std::max({col_molecule, col_span, col_total,
                                                                  col_read_length, col_start,
                                                                  col_end})) {
                    throw std::runtime_error("Short row in segments table " + path + ": " +
                                             std::string(line.s, line.l));
                }
                add(fields[col_molecule], parse_uint(fields[col_span], path),
                    parse_uint(fields[col_total], path),
                    parse_uint(fields[col_read_length], path),
                    parse_uint(fields[col_start], path), parse_uint(fields[col_end], path), path);
            }
        } catch (...) {
            ks_free(&line);
            bgzf_close(fp);
            throw;
        }
        ks_free(&line);
        bgzf_close(fp);
        for (auto& molecule : molecules_) {
            std::sort(segments_.begin() + molecule.first,
                      segments_.begin() + molecule.first + molecule.count,
                      [](const TableSegment& a, const TableSegment& b) {
                          return a.span_index < b.span_index;
                      });
        }
        result.table_molecules = molecules_.size();
        result.table_segments = segments_.size();
        loaded_ = true;
    }

    TableMolecule* find(const std::string& name) {
        auto it = index_.find(name);
        return it == index_.end() ? nullptr : &molecules_[it->second];
    }

    const TableSegment* segment(const TableMolecule& molecule, uint32_t span_index) const {
        auto begin = segments_.begin() + molecule.first;
        auto end = begin + molecule.count;
        auto it = std::lower_bound(begin, end, span_index,
                                   [](const TableSegment& s, uint32_t k) { return s.span_index < k; });
        return (it != end && it->span_index == span_index) ? &*it : nullptr;
    }

    const TableSegment* segments_of(const TableMolecule& molecule) const {
        return segments_.data() + molecule.first;
    }

    uint64_t unseen_molecules() const {
        uint64_t n = 0;
        for (const auto& molecule : molecules_) n += !molecule.seen;
        return n;
    }

private:
    static void split(const std::string& text, std::vector<std::string>& out) {
        out.clear();
        size_t start = 0;
        while (true) {
            size_t tab = text.find('\t', start);
            out.push_back(text.substr(start, tab == std::string::npos ? std::string::npos
                                                                      : tab - start));
            if (tab == std::string::npos) break;
            start = tab + 1;
        }
    }

    static uint64_t parse_uint(const std::string& text, const std::string& path) {
        if (text.empty()) throw std::runtime_error("Empty number in segments table " + path);
        uint64_t value = 0;
        for (char c : text) {
            if (c < '0' || c > '9') {
                throw std::runtime_error("Not a number in segments table " + path + ": " + text);
            }
            value = value * 10 + static_cast<uint64_t>(c - '0');
        }
        return value;
    }

    // Rows of one molecule are adjacent in the table (read order); a
    // molecule coming back later would leave two entries, so it is refused.
    // The read length is a property of the molecule, held once; every row
    // of the molecule repeats it and must agree.
    void add(const std::string& name, uint64_t span, uint64_t total, uint64_t read_length,
             uint64_t start, uint64_t end, const std::string& path) {
        if (molecules_.empty() || name != last_name_) {
            if (index_.count(name)) {
                throw std::runtime_error("Segments table " + path +
                                         " is not grouped by molecule: " + name + " recurs");
            }
            index_.emplace(name, static_cast<uint32_t>(molecules_.size()));
            molecules_.push_back({static_cast<uint32_t>(segments_.size()), 0, read_length, false});
            last_name_ = name;
        } else if (molecules_.back().read_length != read_length) {
            throw std::runtime_error("Segments table " + path + ": the rows of " + name +
                                     " disagree on read_length");
        }
        if (end > read_length) {
            throw std::runtime_error("Segments table " + path + ": segment " +
                                     segment_name(name, static_cast<uint32_t>(span)) +
                                     " ends past read_length");
        }
        segments_.push_back({static_cast<uint32_t>(span), static_cast<uint32_t>(total),
                             start, end});
        molecules_.back().count++;
    }

    bool loaded_ = false;
    std::string last_name_;
    phmap::flat_hash_map<std::string, uint32_t> index_;
    std::vector<TableMolecule> molecules_;
    std::vector<TableSegment> segments_;
};

// --- alignment geometry ----------------------------------------------------

// What the row derives from the CIGAR, following paftools.js sam2paf: the
// query length is taken from the CIGAR (M/=/X, I, S and H bases; a
// secondary record carries no SEQ, so l_qseq cannot serve), the query
// interval swaps its clips on the reverse strand, matching bases are the
// M/=/X bases minus the mismatches (NM minus gap bases, or the X bases of
// an extended CIGAR), and the block length is M/=/X plus I plus D. One
// step beyond sam2paf: minimap2 counts an ambiguous base (N) as a mismatch
// in NM and leaves it out of its PAF block length, reporting the count as
// nn; subtracting nn here makes blen the value minimap2 -c writes.
void alignment_geometry(const bam1_t* b, Candidate& c) {
    const uint32_t* cigar = bam_get_cigar(b);
    uint32_t n = b->core.n_cigar;
    uint64_t M = 0, I = 0, D = 0, X = 0, clipped = 0;
    uint64_t clip_first = 0, clip_last = 0;
    bool have_M = false, have_ext = false;
    c.cigar.clear();
    for (uint32_t i = 0; i < n; i++) {
        uint32_t len = bam_cigar_oplen(cigar[i]);
        char op = bam_cigar_opchr(cigar[i]);
        c.cigar += std::to_string(len);
        c.cigar += op;
        switch (op) {
            case 'M': M += len; have_M = true; break;
            case '=': M += len; have_ext = true; break;
            case 'X': M += len; X += len; have_ext = true; break;
            case 'I': I += len; break;
            case 'D': D += len; break;
            case 'S':
            case 'H':
                clipped += len;
                if (i == 0) clip_first = len;
                else if (i + 1 == n) clip_last = len;
                break;
            default: break;  // N, P: no query bases
        }
    }
    uint64_t qlen = M + I + clipped;
    c.seg_len = qlen;
    if (c.reverse) {
        c.qstart = clip_last;
        c.qend = qlen - clip_first;
    } else {
        c.qstart = clip_first;
        c.qend = qlen - clip_last;
    }
    uint64_t mismatches;
    if (have_ext && !have_M) {
        mismatches = X;
    } else if (c.has_nm) {
        uint64_t gaps = I + D;
        mismatches = c.nm > static_cast<int64_t>(gaps) ? static_cast<uint64_t>(c.nm) - gaps : 0;
    } else {
        mismatches = 0;
    }
    c.mlen = M > mismatches ? M - mismatches : 0;
    uint64_t blen = M + I + D;
    uint64_t ambiguous = c.ambiguous > 0 ? static_cast<uint64_t>(c.ambiguous) : 0;
    c.blen = blen > ambiguous ? blen - ambiguous : 0;
}

bool aux_int(const bam1_t* b, const char tag[2], int64_t& value) {
    uint8_t* aux = bam_aux_get(b, tag);
    if (!aux) return false;
    char type = *aux;
    if (type != 'c' && type != 'C' && type != 's' && type != 'S' && type != 'i' &&
        type != 'I') {
        return false;
    }
    value = bam_aux2i(aux);
    return true;
}

// --- rows -------------------------------------------------------------------

void append_uint(std::string& s, uint64_t v) {
    char tmp[24];
    int n = 0;
    do {
        tmp[n++] = static_cast<char>('0' + v % 10);
        v /= 10;
    } while (v);
    while (n) s += tmp[--n];
}

void append_row(std::string& line, const std::string& name, const Candidate& c,
                uint32_t retained_index, size_t segment_count, const TableMolecule* tm,
                const TableSegment* seg, uint32_t rank, const sam_hdr_t* hdr) {
    line.clear();
    line += name;
    line += '\t'; append_uint(line, c.span_index);
    line += '\t'; append_uint(line, retained_index);
    line += '\t'; append_uint(line, segment_count);
    if (seg) {
        line += '\t'; append_uint(line, seg->spans_total);
        line += '\t'; append_uint(line, tm->read_length);
        line += '\t'; append_uint(line, seg->read_start);
        line += '\t'; append_uint(line, seg->read_end);
    } else {
        line += "\t.\t.\t.\t.";
    }
    line += '\t'; append_uint(line, c.seg_len);
    if (c.aln == 'U') {
        line += "\t.\t.\t.\t.\t.\tU\t";
        append_uint(line, rank);
        line += "\t.\t.\t.\t.\t.\t.\t*\n";
        return;
    }
    line += '\t'; line += sam_hdr_tid2name(hdr, c.tid);
    line += '\t'; append_uint(line, static_cast<uint64_t>(c.ref_start));
    line += '\t'; append_uint(line, static_cast<uint64_t>(c.ref_end));
    line += '\t'; line += c.reverse ? '-' : '+';
    line += '\t'; append_uint(line, static_cast<uint64_t>(c.mapq));
    line += '\t'; line += c.aln;
    line += '\t'; append_uint(line, rank);
    line += '\t'; if (c.has_as) line += std::to_string(c.as); else line += '.';
    line += '\t'; if (c.has_nm) line += std::to_string(c.nm); else line += '.';
    line += '\t'; append_uint(line, c.qstart);
    line += '\t'; append_uint(line, c.qend);
    line += '\t'; append_uint(line, c.mlen);
    line += '\t'; append_uint(line, c.blen);
    line += '\t'; line += c.cigar;
    line += '\n';
}

// Rank order within a segment: the first primary-class record, then the
// rest by score (missing scores last), then reference interval (contig,
// start, end), forward strand before reverse, CIGAR text, and only then
// arrival, so that records that differ in anything the row shows rank the
// same whatever order the BAM had them in. The whole molecule sorts by
// span first.
bool before(const Candidate& a, const Candidate& b) {
    if (a.span_index != b.span_index) return a.span_index < b.span_index;
    if (a.first_primary != b.first_primary) return a.first_primary;
    int64_t sa = a.has_as ? a.as : INT64_MIN;
    int64_t sb = b.has_as ? b.as : INT64_MIN;
    if (sa != sb) return sa > sb;
    if (a.tid != b.tid) return a.tid < b.tid;
    if (a.ref_start != b.ref_start) return a.ref_start < b.ref_start;
    if (a.ref_end != b.ref_end) return a.ref_end < b.ref_end;
    if (a.reverse != b.reverse) return !a.reverse;
    if (a.cigar != b.cigar) return a.cigar < b.cigar;
    return a.order < b.order;
}

void flush_molecule(Molecule& m, SegmentsTable& table, const sam_hdr_t* hdr,
                    TableWriter& out, const MoleculesConfig& config,
                    MoleculesResult& result) {
    if (!m.active) return;
    auto& records = m.records;

    // The molecule in the table, when there is one: every BAM segment must
    // be there with the same length, and table segments the BAM lacks are
    // completed as unmapped rows.
    TableMolecule* tm = nullptr;
    if (table.loaded()) {
        tm = table.find(m.name);
        if (!tm) {
            uint32_t first = UINT32_MAX;
            for (const auto& c : records) first = std::min(first, c.span_index);
            throw std::runtime_error("molecule " + m.name + " is in the BAM but not in the "
                                     "segments table (" + segment_name(m.name, first) + ")");
        }
        tm->seen = true;
        phmap::flat_hash_map<uint32_t, uint64_t> longest;   // span -> longest record
        for (const auto& c : records) {
            auto& len = longest[c.span_index];
            len = std::max(len, c.seg_len);
        }
        // Checked in span order, so the segment an error names is the first
        std::vector<std::pair<uint32_t, uint64_t>> in_bam(longest.begin(), longest.end());
        std::sort(in_bam.begin(), in_bam.end());
        for (const auto& [span, seg_len] : in_bam) {
            const TableSegment* seg = table.segment(*tm, span);
            std::string qname = segment_name(m.name, span);
            if (!seg) {
                throw std::runtime_error("segment " + qname +
                                         " is in the BAM but not in the segments table");
            }
            uint64_t table_len = seg->read_end - seg->read_start;
            if (seg_len != table_len) {
                throw std::runtime_error(
                    "segment " + qname + " is " + std::to_string(seg_len) +
                    " bp in the BAM and " + std::to_string(table_len) +
                    " bp in the segments table; is the table from the same digest?");
            }
        }
        const TableSegment* segs = table.segments_of(*tm);
        for (uint32_t i = 0; i < tm->count; i++) {
            if (longest.count(segs[i].span_index)) continue;
            Candidate u{};
            u.span_index = segs[i].span_index;
            u.aln = 'U';
            u.first_primary = true;
            u.tid = -1;
            u.seg_len = segs[i].read_end - segs[i].read_start;
            u.cigar = "*";
            u.order = records.size();
            records.push_back(std::move(u));
            result.table_segments_missing_from_bam++;
        }
    }

    std::sort(records.begin(), records.end(), before);

    // Dense order of the distinct spans, which of them have several records,
    // and one length per segment: the records of a segment agree on it
    // when every CIGAR carries its clips, and the longest stands otherwise.
    std::vector<uint32_t> spans;
    bool molecule_has_candidates = false;
    for (size_t i = 0, first = 0; i < records.size(); i++) {
        bool last = i + 1 == records.size() || records[i + 1].span_index != records[i].span_index;
        if (i == 0 || records[i].span_index != records[i - 1].span_index) {
            spans.push_back(records[i].span_index);
            first = i;
        }
        if (last) {
            if (i > first) {
                result.segments_with_candidates++;
                molecule_has_candidates = true;
            }
            uint64_t seg_len = 0;
            for (size_t j = first; j <= i; j++) seg_len = std::max(seg_len, records[j].seg_len);
            for (size_t j = first; j <= i; j++) records[j].seg_len = seg_len;
        }
    }
    result.molecules++;
    result.segments += spans.size();
    if (molecule_has_candidates) result.molecules_with_candidates++;

    std::string line;
    uint32_t retained_index = 0, rank = 0;
    for (size_t i = 0; i < records.size(); i++) {
        const Candidate& c = records[i];
        if (i == 0 || c.span_index != records[i - 1].span_index) {
            retained_index++;
            rank = 0;
        } else {
            rank++;
        }
        if (config.primary_only && (c.aln == 'S' || c.aln == 'L')) {
            result.rows_dropped_candidates++;
            continue;
        }
        const TableSegment* seg = tm ? table.segment(*tm, c.span_index) : nullptr;
        append_row(line, m.name, c, retained_index, spans.size(), tm, seg, rank, hdr);
        out.write(line);
        result.rows_written++;
    }
}

}  // namespace

MoleculesResult extract_molecules(
    const std::string& input_path,
    const std::string& output_path,
    const MoleculesConfig& config
) {
    MoleculesResult result;
    SegmentsTable table;
    if (!config.table.empty()) table.load(config.table, result);

    HtsFilePtr fp(hts_open(input_path.c_str(), "r"));
    if (!fp) throw std::runtime_error("Cannot open: " + input_path);
    if (config.threads > 1) hts_set_threads(fp.get(), config.threads);

    SamHeaderPtr hdr(sam_hdr_read(fp.get()));
    if (!hdr) throw std::runtime_error("Cannot read header: " + input_path);
    read_group_order(hdr.get(), input_path, "cifi molecules", result.sort_order,
                     result.group_order);

    // ##key=value metadata, then the #columns line. This table is the
    // canonical one; a view derived from it says canonical=false.
    TableWriter out(output_path);
    out.header({
        "##cifi_format=molecules",
        "##format_version=1",
        "##tool_version=" + config.version,
        "##command=" + config.command,
        "##canonical=true",
        "##read_coordinates=0-based-half-open",
        "##reference_coordinates=0-based-half-open",
        "##span_index=1-based-original-digest-span",
        "##retained_index=1-based-dense-retained-order",
        "##aln=P-primary,S-secondary,L-supplementary,U-unmapped",
        "##rank=0-primary,then-AS-descending",
        "##segments_table=" + (config.table.empty() ? std::string("none") : config.table),
        std::string("#columns: ") + MOLECULES_COLUMNS,
    });

    Molecule molecule;
    RecentReads recent(256);
    BamRecordPtr rec(bam_init1());
    bam1_t* b = rec.get();

    int ret;
    while ((ret = sam_read1(fp.get(), hdr.get(), b)) >= 0) {
        result.records_seen++;
        auto seg = parse_segment_name(bam_get_qname(b));

        if (!molecule.active || seg.read != molecule.name) {
            flush_molecule(molecule, table, hdr.get(), out, config, result);
            if (molecule.active) recent.push(molecule.name);
            if (recent.contains(seg.read)) {
                throw std::runtime_error(
                    "segments of read " + seg.read + " are not grouped together in " +
                    input_path + "; sort by name first: samtools sort -n");
            }
            molecule.start(seg.read);
        }

        Candidate c{};
        c.span_index = seg.index;
        c.order = molecule.records.size();
        uint16_t flag = b->core.flag;
        bool primary_class = !(flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY));
        if (primary_class) {
            c.first_primary = molecule.primaries.insert(seg.index).second;
            if (!c.first_primary) result.duplicate_primary++;
        }
        if ((flag & BAM_FUNMAP) || b->core.tid < 0) {
            c.aln = 'U';
            c.tid = -1;
            c.seg_len = static_cast<uint64_t>(b->core.l_qseq);
            c.cigar = "*";
            result.unmapped++;
        } else {
            c.aln = (flag & BAM_FSECONDARY) ? 'S' : (flag & BAM_FSUPPLEMENTARY) ? 'L' : 'P';
            if (c.aln == 'S') result.secondary++;
            else if (c.aln == 'L') result.supplementary++;
            else result.primary_mapped++;
            c.tid = b->core.tid;
            c.ref_start = b->core.pos;
            c.ref_end = bam_endpos(b);
            c.reverse = flag & BAM_FREVERSE;
            c.mapq = b->core.qual;
            c.has_as = aux_int(b, "AS", c.as);
            c.has_nm = aux_int(b, "NM", c.nm);
            if (!aux_int(b, "nn", c.ambiguous)) c.ambiguous = 0;
            alignment_geometry(b, c);
        }
        molecule.records.push_back(std::move(c));
    }
    if (ret < -1) throw std::runtime_error("Truncated or corrupt input: " + input_path);

    flush_molecule(molecule, table, hdr.get(), out, config, result);
    if (table.loaded()) result.table_molecules_missing_from_bam = table.unseen_molecules();
    out.close();
    return result;
}

}  // namespace cifi

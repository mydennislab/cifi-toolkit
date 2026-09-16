"""Unique-segment output of `cifi digest --segments-out`.

Each retained segment of a passing read is written exactly once, named
<read>__CIFI_SEG__<span>, so it can be mapped a single time and the pairwise
contacts reconstructed afterwards. The R1/R2 output must not change when the
option is used, and R2 reverse complementation must not reach this file.
"""

import gzip
import os
import subprocess
import sys

import pytest

import cifi
from cifi import parse_segment_name, process_reads, segment_name

HINDIII_SITE = "AAGCTT"
HINDIII_CUT_OFFSET = 1
OVERHANG = len(HINDIII_SITE) - HINDIII_CUT_OFFSET  # 5, i.e. "AGCTT"
SEP = "__CIFI_SEG__"

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCATGCA")


def revcomp(seq):
    return seq.translate(_COMPLEMENT)[::-1]


def block(unit, length):
    seq = (unit * (length // len(unit) + 1))[:length]
    assert HINDIII_SITE not in seq
    return seq


def make_read(*source_lengths):
    """A read whose HindIII segments have exactly `source_lengths` bp (pre-trim)."""
    units = ["ACGT", "ACCG", "AGGC", "ATTG", "ACTA", "AGTC"]
    parts = []
    for i, slen in enumerate(source_lengths):
        last = i == len(source_lengths) - 1
        filler = slen - (0 if i == 0 else OVERHANG) - (0 if last else 1)
        assert filler >= 0, f"segment {i} too short to construct"
        parts.append(block(units[i % len(units)], filler))
    return HINDIII_SITE.join(parts)


def ramped_quality(length):
    return "".join(chr(33 + ((i * 7) % 41)) for i in range(length))


def expected_spans(sequence, min_emit_len, lead_trim=OVERHANG):
    """(span_index, start, end) of every kept segment, per the digest model."""
    sites, start = [], 0
    while (idx := sequence.find(HINDIII_SITE, start)) != -1:
        sites.append(idx)
        start = idx + 1
    cuts = [0] + [p + HINDIII_CUT_OFFSET for p in sites] + [len(sequence)]
    out = []
    for k, (a, b) in enumerate(zip(cuts, cuts[1:]), start=1):
        s = min(a + lead_trim, b) if a > 0 else a
        if b > s and b - s >= min_emit_len:
            out.append((k, s, b))
    return out


def write_fastq(path, records):
    with open(path, "w") as fh:
        for name, seq, qual in records:
            fh.write(f"@{name}\n{seq}\n+\n{qual}\n")


def read_fastq(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        lines = [ln.rstrip("\n") for ln in fh]
    return [(lines[i][1:], lines[i + 1], lines[i + 3]) for i in range(0, len(lines), 4)]


def run_digest(tmp_path, records, *, min_segments=2, min_segment_len=60,
               strip_overhang=True, revcomp_r2=False, segments_name="segments.fastq"):
    """records: list of (name, seq) or (name, seq, qual). Returns (r1, r2, segs, result)."""
    tmp_path.mkdir(parents=True, exist_ok=True)
    recs = [r if len(r) == 3 else (r[0], r[1], "I" * len(r[1])) for r in records]
    fq = tmp_path / "in.fastq"
    write_fastq(fq, recs)
    r1, r2 = tmp_path / "out_R1.fastq", tmp_path / "out_R2.fastq"
    segs = tmp_path / segments_name if segments_name else None
    result = process_reads(
        str(fq), str(r1), str(r2), "HindIII",
        min_segments, min_segment_len, strip_overhang, False, False, revcomp_r2,
        str(segs) if segs else "",
    )
    seg_records = read_fastq(segs) if segs else []
    return read_fastq(r1), read_fastq(r2), seg_records, result


# --- A: every retained segment once, in order, with its quality ------------

def test_four_segments_give_four_unique_records_and_six_pairs(tmp_path):
    read = make_read(400, 300, 250, 200)
    qual = ramped_quality(len(read))
    r1, r2, segs, result = run_digest(tmp_path, [("read1", read, qual)])

    spans = expected_spans(read, 60)
    assert len(spans) == 4
    assert [name for name, _, _ in segs] == [f"read1{SEP}{k}" for k, _, _ in spans]
    assert [seq for _, seq, _ in segs] == [read[s:e] for _, s, e in spans]
    assert [q for _, _, q in segs] == [qual[s:e] for _, s, e in spans]
    assert len(r1) == len(r2) == result.pairs_written == 6
    assert result.segments_written == 4
    assert result.bases_out_segments == sum(e - s for _, s, e in spans)


def test_segment_indices_follow_read_order(tmp_path):
    _, _, segs, _ = run_digest(tmp_path, [("r", make_read(400, 300, 250, 200, 150))])

    indices = [parse_segment_name(name)[1] for name, _, _ in segs]
    assert indices == [1, 2, 3, 4, 5]


# --- B: a dropped segment leaves a gap in the numbering ------------------

def test_short_segment_is_absent_from_segments_and_pairs(tmp_path):
    """Source 64bp emits 59bp against a 60bp cutoff; it must vanish everywhere."""
    read = make_read(400, 60 + OVERHANG - 1, 300, 200)
    r1, r2, segs, result = run_digest(tmp_path, [("read1", read)])

    assert [parse_segment_name(n)[1] for n, _, _ in segs] == [1, 3, 4]
    assert [seq for _, seq, _ in segs] == [read[s:e] for _, s, e in expected_spans(read, 60)]
    assert result.segments_written == result.total_segments == 3
    assert len(r1) == len(r2) == result.pairs_written == 3
    # the dropped bases appear in neither output
    dropped = expected_spans(read, 1)[1]
    assert not any(read[dropped[1]:dropped[2]] in seq for _, seq, _ in segs + r1 + r2)


def test_naming_is_deterministic_across_runs(tmp_path):
    read = make_read(400, 60 + OVERHANG - 1, 300, 200)
    _, _, a, _ = run_digest(tmp_path / "a", [("read1", read)])
    _, _, b, _ = run_digest(tmp_path / "b", [("read1", read)])

    assert a == b


# --- C: reads failing min-segments contribute nothing ---------------------

def test_read_below_min_segments_emits_no_segments(tmp_path):
    r1, r2, segs, result = run_digest(tmp_path, [("read1", make_read(400, 300))],
                                      min_segments=3)

    assert segs == [] and r1 == [] and r2 == []
    assert result.segments_written == 0
    assert result.reads_out == 0


def test_read_left_with_one_segment_after_length_filter_emits_nothing(tmp_path):
    read = make_read(400, 60 + OVERHANG - 1)  # second segment drops, one remains
    _, _, segs, result = run_digest(tmp_path, [("read1", read)], min_segments=2)

    assert segs == []
    assert result.segments_written == 0
    assert result.filtered_short_segments == 1


def test_only_passing_reads_reach_the_segments_file(tmp_path):
    reads = [
        ("pass", make_read(400, 300, 200)),
        ("few_sites", block("ACGT", 500)),
        ("too_few", make_read(400, 300)),
    ]
    _, _, segs, result = run_digest(tmp_path, reads, min_segments=3)

    assert {parse_segment_name(n)[0] for n, _, _ in segs} == {"pass"}
    assert result.segments_written == 3


# --- D: overhang stripping trims sequence and quality alike --------------

def test_strip_overhang_trims_quality_with_the_sequence(tmp_path):
    read = make_read(400, 300, 250)
    qual = ramped_quality(len(read))
    _, _, segs, _ = run_digest(tmp_path, [("read1", read, qual)])

    for (name, seq, q), (_, s, e) in zip(segs, expected_spans(read, 60)):
        assert (seq, q) == (read[s:e], qual[s:e]), name
        assert len(seq) == len(q)
    # every segment after the first begins at a cut, so its remnant is gone
    for _, seq, _ in segs[1:]:
        assert not seq.startswith("AGCTT")


def test_no_strip_keeps_the_remnant_and_its_quality(tmp_path):
    read = make_read(400, 300, 250)
    qual = ramped_quality(len(read))
    _, _, segs, _ = run_digest(tmp_path, [("read1", read, qual)], strip_overhang=False)

    spans = expected_spans(read, 60, lead_trim=0)
    assert [(seq, q) for _, seq, q in segs] == [(read[s:e], qual[s:e]) for _, s, e in spans]
    for _, seq, _ in segs[1:]:
        assert seq.startswith("AGCTT")


# --- E: R2 reverse complementation stays a paired-FASTQ operation ---------

def test_revcomp_r2_changes_r2_but_not_the_segments(tmp_path):
    read = make_read(400, 300, 250)
    qual = ramped_quality(len(read))
    _, r2_plain, segs_plain, _ = run_digest(tmp_path / "plain", [("read1", read, qual)])
    _, r2_rc, segs_rc, _ = run_digest(tmp_path / "rc", [("read1", read, qual)], revcomp_r2=True)

    assert segs_rc == segs_plain
    assert r2_rc != r2_plain
    assert [seq for _, seq, _ in r2_rc] == [revcomp(seq) for _, seq, _ in r2_plain]
    # native orientation: the segment reads as it lies in the source read
    for _, seq, _ in segs_rc:
        assert seq in read


# --- the pairs output is untouched by the option --------------------------

def test_pairs_output_is_identical_with_and_without_segments_out(tmp_path):
    reads = [("a", make_read(400, 300, 250, 200)), ("b", make_read(400, 60 + OVERHANG - 1, 300))]
    r1_without, r2_without, _, res_without = run_digest(tmp_path / "w", reads, segments_name=None)
    r1_with, r2_with, segs, res_with = run_digest(tmp_path / "s", reads)

    assert (r1_with, r2_with) == (r1_without, r2_without)
    assert res_with.pairs_written == res_without.pairs_written
    assert res_without.segments_written == 0
    assert res_with.segments_written == len(segs) == 6


def test_pair_names_keep_their_kept_order_scheme(tmp_path):
    """The two naming schemes coexist: pairs by kept order, segments by span."""
    read = make_read(400, 60 + OVERHANG - 1, 300, 200)
    r1, _, segs, _ = run_digest(tmp_path, [("read1", read)])

    assert [n for n, _, _ in r1] == ["read1_0_0", "read1_0_1", "read1_1_0"]
    assert [n for n, _, _ in segs] == [f"read1{SEP}1", f"read1{SEP}3", f"read1{SEP}4"]


# --- naming contract -----------------------------------------------------

@pytest.mark.parametrize("read_name", [
    "m84039_240101_000000_s1/44242450/ccs",
    "m64411e_240417_155217/1/ccs:0026:0332",
    "read:with:colons_and__underscores",
    "a",
    "1__2__3",
])
def test_segment_name_round_trips_awkward_read_names(read_name):
    for k in (1, 7, 123456):
        assert parse_segment_name(segment_name(read_name, k)) == (read_name, k)


@pytest.mark.parametrize("bad", [
    "plain_read_name",
    "read" + SEP,            # no index
    "read" + SEP + "x",      # not a number
    "read" + SEP + "0",      # spans count from 1
    "read" + SEP + "-1",
    "read" + SEP + "1.5",
    SEP + "3",               # no read name
    "",
])
def test_parse_rejects_names_outside_the_contract(bad):
    with pytest.raises(ValueError):
        parse_segment_name(bad)


def test_digest_refuses_read_names_that_already_carry_the_marker(tmp_path):
    with pytest.raises(RuntimeError, match=SEP):
        run_digest(tmp_path, [(f"read{SEP}1", make_read(400, 300, 250))])


def test_segment_names_survive_a_sam_round_trip(tmp_path):
    """The name written to FASTQ is what an aligner puts in QNAME, verbatim."""
    _, _, segs, _ = run_digest(tmp_path, [("m84039/1/ccs", make_read(400, 300, 250))])

    for name, _, _ in segs:
        assert " " not in name and "\t" not in name
        read, k = parse_segment_name(name)
        assert read == "m84039/1/ccs" and k >= 1


# --- compression by extension, input formats ------------------------------

def test_gz_extension_writes_gzip_and_plain_writes_plain(tmp_path):
    read = make_read(400, 300, 250)
    _, _, segs_gz, _ = run_digest(tmp_path / "gz", [("r", read)], segments_name="s.fastq.gz")
    _, _, segs_plain, _ = run_digest(tmp_path / "plain", [("r", read)], segments_name="s.fastq")

    assert (tmp_path / "gz" / "s.fastq.gz").read_bytes()[:2] == b"\x1f\x8b"
    assert (tmp_path / "plain" / "s.fastq").read_bytes()[:1] == b"@"
    assert segs_gz == segs_plain


def test_sam_input_yields_the_same_segments_as_fastq(tmp_path):
    read = make_read(400, 300, 250)
    qual = ramped_quality(len(read))
    _, _, from_fastq, _ = run_digest(tmp_path / "fq", [("r1", read, qual)])

    sam = tmp_path / "in.sam"
    sam.write_text("@HD\tVN:1.6\tSO:unknown\n"
                   f"r1\t4\t*\t0\t0\t*\t*\t0\t0\t{read}\t{qual}\n")
    r1, r2 = tmp_path / "o_R1.fastq", tmp_path / "o_R2.fastq"
    process_reads(str(sam), str(r1), str(r2), "HindIII", 2, 60, True, False, False, False,
                  str(tmp_path / "segs.fastq"))

    assert read_fastq(tmp_path / "segs.fastq") == from_fastq


def test_cli_reports_the_segments_file(tmp_path):
    fq = tmp_path / "in.fastq"
    write_fastq(fq, [("r1", make_read(400, 300, 250, 200), "I" * len(make_read(400, 300, 250, 200)))])
    segs = tmp_path / "out.segments.fastq.gz"
    proc = subprocess.run(
        [sys.executable, "-m", "cifi.cli", "digest", str(fq), "-e", "HindIII",
         "-o", str(tmp_path / "out"), "--segments-out", str(segs), "--no-report"],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stderr
    assert str(segs) in proc.stdout
    assert len(read_fastq(segs)) == 4
    import json
    stats = json.loads((tmp_path / "out_stats.json").read_text())
    assert stats["results"]["segments_written"] == 4
    assert stats["output"]["segments"] == str(segs)


@pytest.mark.skipif(not os.path.isdir("/proc/self/fd"), reason="needs /proc to count open files")
def test_a_refused_bam_read_does_not_leak_the_open_file(tmp_path):
    """The reader throws mid-file; the htslib handle must go with the exception."""
    sam = tmp_path / "in.sam"
    sam.write_text("@HD\tVN:1.6\tSO:unknown\n"
                   f"r1\t4\t*\t0\t0\t*\t*\t0\t0\t{make_read(400, 300, 250)}\t*\n"
                   f"read{SEP}1\t4\t*\t0\t0\t*\t*\t0\t0\t{make_read(400, 300, 250)}\t*\n")
    r1, r2 = tmp_path / "o_R1.fastq", tmp_path / "o_R2.fastq"

    def failing_digest():
        with pytest.raises(RuntimeError, match=SEP):
            process_reads(str(sam), str(r1), str(r2), "HindIII", 2, 60, True, False, False, False,
                          str(tmp_path / "segs.fastq"))

    failing_digest()   # first call: whatever lazy state htslib keeps is set up now
    before = len(os.listdir("/proc/self/fd"))
    for _ in range(20):
        failing_digest()
    assert len(os.listdir("/proc/self/fd")) == before


# --- the segments table: one row per retained segment ----------------------
#
# `--segments-table` writes the digest's view of each molecule: where every
# retained segment sat in the read, which cut span it came from, and what
# the trim removed. Rows correspond one to one with --segments-out records.

TABLE_COLUMNS = [
    "molecule_id", "span_index", "retained_index", "segments_kept", "spans_total",
    "read_length", "read_start", "read_end", "span_start", "span_end", "original_len",
    "processed_len", "trimmed_5p", "terminal", "enzyme", "cut_offset",
]


def read_table(path):
    """(header lines, rows as dicts) of a bgzip TSV with a #columns line."""
    with gzip.open(path, "rt") as fh:
        lines = [ln.rstrip("\n") for ln in fh]
    header = [ln for ln in lines if ln.startswith("#")]
    columns = next(ln for ln in header if ln.startswith("#columns:"))
    names = columns[len("#columns:"):].split()
    rows = [dict(zip(names, ln.split("\t"))) for ln in lines if not ln.startswith("#")]
    for row in rows:
        assert len(row) == len(names), row
    return header, rows


def header_metadata(header):
    """The ##key=value lines of a native table header as a dict, checking
    the shape: every line but the last is ##key=value, keys are unique,
    and the last line is #columns:."""
    assert header[-1].startswith("#columns: "), header[-1]
    meta = {}
    for line in header[:-1]:
        assert line.startswith("##") and "=" in line, line
        key, value = line[2:].split("=", 1)
        assert key and key not in meta, line
        meta[key] = value
    return meta


def run_digest_table(tmp_path, records, *, min_segments=2, min_segment_len=60,
                     strip_overhang=True, segments_name="segments.fastq", enzyme="HindIII",
                     site=None, cut_offset=None, command="cifi digest test"):
    """Like run_digest, with --segments-table; returns (segs, table rows, header, result)."""
    from cifi import process_reads_custom
    tmp_path.mkdir(parents=True, exist_ok=True)
    recs = [r if len(r) == 3 else (r[0], r[1], "I" * len(r[1])) for r in records]
    fq = tmp_path / "in.fastq"
    write_fastq(fq, recs)
    r1, r2 = tmp_path / "out_R1.fastq", tmp_path / "out_R2.fastq"
    segs = tmp_path / segments_name if segments_name else None
    table = tmp_path / "segments.tsv.gz"
    args = (str(fq), str(r1), str(r2))
    common = (min_segments, min_segment_len, strip_overhang, False, False, False,
              str(segs) if segs else "", str(table), command, "1.1.0-test")
    if site is not None:
        result = process_reads_custom(*args, site, cut_offset, *common)
    else:
        result = process_reads(*args, enzyme, *common)
    header, rows = read_table(table)
    return (read_fastq(segs) if segs else []), rows, header, result


def test_table_rows_correspond_to_the_segment_records(tmp_path):
    reads = [("a", make_read(400, 300, 250, 200)), ("b", make_read(400, 60 + OVERHANG - 1, 300))]
    segs, rows, header, result = run_digest_table(tmp_path, reads)

    assert [f"{r['molecule_id']}{SEP}{r['span_index']}" for r in rows] == [n for n, _, _ in segs]
    assert result.segments_table_rows == len(rows) == len(segs) == 6
    assert list(rows[0]) == TABLE_COLUMNS


def test_table_indices_and_totals(tmp_path):
    """span_index has the gap of the dropped span; retained_index is dense."""
    read = make_read(400, 60 + OVERHANG - 1, 300, 200)
    _, rows, _, _ = run_digest_table(tmp_path, [("r", read)])

    assert [int(r["span_index"]) for r in rows] == [1, 3, 4]
    assert [int(r["retained_index"]) for r in rows] == [1, 2, 3]
    assert {r["segments_kept"] for r in rows} == {"3"}
    assert {r["spans_total"] for r in rows} == {"4"}
    assert {r["molecule_id"] for r in rows} == {"r"}


def test_table_coordinates_slice_the_read_to_the_segment_sequence(tmp_path):
    read = make_read(400, 300, 60 + OVERHANG - 1, 250, 200)
    qual = ramped_quality(len(read))
    segs, rows, _, _ = run_digest_table(tmp_path, [("r", read, qual)])

    spans = expected_spans(read, 60)
    assert [(int(r["read_start"]), int(r["read_end"])) for r in rows] == [(s, e) for _, s, e in spans]
    assert {r["read_length"] for r in rows} == {str(len(read))}
    assert all(int(r["read_end"]) <= int(r["read_length"]) for r in rows)
    for row, (_, seq, _) in zip(rows, segs):
        s, e = int(row["read_start"]), int(row["read_end"])
        assert read[s:e] == seq
        assert int(row["processed_len"]) == e - s == len(seq)
        ss, se = int(row["span_start"]), int(row["span_end"])
        assert int(row["original_len"]) == se - ss
        assert int(row["trimmed_5p"]) == s - ss
        assert se == e
    # the leading span starts at the read start and carries no remnant
    assert rows[0]["span_start"] == "0" and rows[0]["trimmed_5p"] == "0"
    # every later retained span begins at a cut and lost the 5 bp remnant
    assert {r["trimmed_5p"] for r in rows[1:]} == {str(OVERHANG)}
    # the untrimmed spans tile the read between the cuts
    untrimmed = expected_spans(read, 1, lead_trim=0)
    by_index = {k: (s, e) for k, s, e in untrimmed}
    for row in rows:
        assert (int(row["span_start"]), int(row["span_end"])) == by_index[int(row["span_index"])]


def test_table_without_strip_keeps_span_and_read_coordinates_equal(tmp_path):
    read = make_read(400, 300, 250)
    _, rows, _, _ = run_digest_table(tmp_path, [("r", read)], strip_overhang=False)

    for row in rows:
        assert row["trimmed_5p"] == "0"
        assert row["read_start"] == row["span_start"] and row["read_end"] == row["span_end"]
        assert row["original_len"] == row["processed_len"]


def test_table_terminal_flags(tmp_path):
    reads = [
        ("four", make_read(400, 300, 250, 200)),
        ("gap_at_end", make_read(400, 300, 60 + OVERHANG - 1)),   # last span dropped
        ("single", block("ACGT", 500)),                            # no site at all
    ]
    _, rows, _, _ = run_digest_table(tmp_path, reads, min_segments=1)

    flags = {(r["molecule_id"], int(r["span_index"])): r["terminal"] for r in rows}
    assert flags[("four", 1)] == "T5" and flags[("four", 4)] == "T3"
    assert flags[("four", 2)] == flags[("four", 3)] == "I"
    # the retained last segment of gap_at_end is span 2 of 3: interior
    assert flags[("gap_at_end", 1)] == "T5" and flags[("gap_at_end", 2)] == "I"
    assert ("gap_at_end", 3) not in flags
    assert flags[("single", 1)] == "T5T3"
    single = next(r for r in rows if r["molecule_id"] == "single")
    assert (single["spans_total"], single["segments_kept"]) == ("1", "1")


def test_table_records_the_enzyme_and_cut(tmp_path):
    _, rows_named, _, _ = run_digest_table(tmp_path / "named", [("r", make_read(400, 300, 250))])
    read_gatc = "GATC".join([block("ACGT", 300), block("ACCG", 200), block("AGGC", 250)])
    _, rows_custom, _, _ = run_digest_table(tmp_path / "custom", [("r", read_gatc)],
                                            site="GATC", cut_offset=0)

    assert {(r["enzyme"], r["cut_offset"]) for r in rows_named} == {("HindIII", "1")}
    assert {(r["enzyme"], r["cut_offset"]) for r in rows_custom} == {("GATC", "0")}
    assert len(rows_custom) == 3


def test_table_header_lines_and_bgzip(tmp_path):
    _, rows, header, _ = run_digest_table(tmp_path, [("r", make_read(400, 300, 250))],
                                          command="cifi digest in.fastq -e HindIII -o out")

    assert header == [
        "##cifi_format=segments",
        "##format_version=1",
        "##tool_version=1.1.0-test",
        "##command=cifi digest in.fastq -e HindIII -o out",
        "##read_coordinates=0-based-half-open",
        "##span_index=1-based-original-digest-span",
        "##retained_index=1-based-dense-retained-order",
        "##terminal=T5,T3,I,T5T3",
        "#columns: " + "\t".join(TABLE_COLUMNS),
    ]
    assert header_metadata(header)["cifi_format"] == "segments"
    raw = (tmp_path / "segments.tsv.gz").read_bytes()
    # gzip magic, then the BGZF extra subfield 'BC' that marks a block
    assert raw[:2] == b"\x1f\x8b" and raw[12:14] == b"BC"
    # the file is a run of complete blocks, ending in the 28-byte EOF block
    assert raw.endswith(bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000"))


def test_table_is_written_without_segments_out(tmp_path):
    segs, rows, _, result = run_digest_table(tmp_path, [("r", make_read(400, 300, 250))],
                                             segments_name=None)

    assert segs == [] and len(rows) == 3
    assert result.segments_written == 0 and result.segments_table_rows == 3


def test_table_lists_only_passing_reads(tmp_path):
    reads = [
        ("pass", make_read(400, 300, 200)),
        ("few_sites", block("ACGT", 500)),
        ("too_few", make_read(400, 300)),
    ]
    _, rows, _, _ = run_digest_table(tmp_path, reads, min_segments=3)

    assert {r["molecule_id"] for r in rows} == {"pass"}


def test_fastq_outputs_are_byte_identical_with_and_without_the_table(tmp_path):
    import hashlib

    def md5(path):
        return hashlib.md5(path.read_bytes()).hexdigest()

    reads = [("a", make_read(400, 300, 250, 200), ramped_quality(len(make_read(400, 300, 250, 200)))),
             ("b", make_read(400, 60 + OVERHANG - 1, 300))]
    run_digest(tmp_path / "w", reads)
    run_digest_table(tmp_path / "t", reads)

    for name in ("out_R1.fastq", "out_R2.fastq", "segments.fastq"):
        assert md5(tmp_path / "w" / name) == md5(tmp_path / "t" / name), name
    assert not (tmp_path / "w" / "segments.tsv.gz").exists()


def test_cli_writes_the_table_and_records_it_in_the_stats(tmp_path):
    import json
    fq = tmp_path / "in.fastq"
    read = make_read(400, 300, 250, 200)
    write_fastq(fq, [("r1", read, "I" * len(read))])
    table = tmp_path / "out.segments.tsv.gz"
    proc = subprocess.run(
        [sys.executable, "-m", "cifi.cli", "digest", str(fq), "-e", "HindIII",
         "-o", str(tmp_path / "out"), "--segments-table", str(table), "--no-report"],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stderr
    assert str(table) in proc.stdout
    header, rows = read_table(table)
    assert len(rows) == 4
    meta = header_metadata(header)
    assert meta["command"].startswith("cifi digest ") and "--segments-table" in meta["command"]
    assert meta["tool_version"] == cifi.__version__
    stats = json.loads((tmp_path / "out_stats.json").read_text())
    assert stats["output"]["segments_table"] == str(table)
    assert stats["results"]["segments_table_rows"] == 4
    assert stats["parameters"]["segments_table"] == str(table)
    # the statistics repeat the table's provenance lines
    assert stats["command"] == meta["command"]
    assert stats["cifi_version"] == meta["tool_version"]
    assert stats["formats"] == {"segments_table": {"name": "segments", "version": 1}}

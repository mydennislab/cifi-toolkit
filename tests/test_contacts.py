"""`cifi contacts`: pairwise contacts from mapped unique segments.

The input is a name-grouped BAM (SAM here; htslib reads both) of segments
named by `cifi digest --segments-out`. Every read with k usable segments
must give exactly C(k,2) contacts, never a self pair, never a duplicate, and
never a pair across reads. PA5 positions follow yahs's own name-sorted-BAM
convention: the 0-based alignment midpoint. BED carries the aligned span of
each segment instead, two consecutive records per contact.
"""

import gzip
import json
import random
import re
import shutil
import subprocess
import sys

import pytest

from cifi import pa5_position, reconstruct_contacts

SEP = "__CIFI_SEG__"
CONTIGS = (("ctg1", 100000), ("ctg2", 50000))


def query_length(cigar):
    return sum(int(n) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar) if op in "MIS=X")


def sam_record(qname, flag, rname, pos, mapq, cigar):
    """One SAM line; SEQ is sized to the CIGAR so htslib accepts it."""
    if flag & 4:
        return f"{qname}\t{flag}\t*\t0\t0\t*\t*\t0\t0\tACGT\t*"
    seq = "A" * query_length(cigar)
    return f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{seq}\t*"


def write_sam(path, records, sort_order="queryname", contigs=CONTIGS, group_order=None):
    hd = "@HD\tVN:1.6" + (f"\tSO:{sort_order}" if sort_order else "")
    header = [hd + (f"\tGO:{group_order}" if group_order else "")]
    header += [f"@SQ\tSN:{name}\tLN:{length}" for name, length in contigs]
    with open(path, "w") as fh:
        fh.write("\n".join(header + [sam_record(*r) for r in records]) + "\n")


def read_rows(path, columns):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        rows = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]
    for row in rows:
        assert len(row) == columns, row
    return rows


def read_pa5(path):
    return read_rows(path, 7)


def read_bed(path):
    return read_rows(path, 5)


def bed_pairs(rows):
    """(name, record, record) per two consecutive rows, as yahs reads them."""
    assert len(rows) % 2 == 0, "BED rows come in pairs"
    pairs = []
    for a, b in zip(rows[0::2], rows[1::2]):
        assert a[3] == b[3], (a, b)
        pairs.append((a[3], a, b))
    return pairs


def bed_as_pa5(rows):
    """Reduce BED pairs to the PA5 rows yahs would derive: floor((start + end) / 2)."""
    def mid(rec):
        return str((int(rec[1]) + int(rec[2])) // 2)
    return [[name, a[0], mid(a), b[0], mid(b), a[4], b[4]] for name, a, b in bed_pairs(rows)]


def split_pair_name(name):
    read, i, j = name.rsplit(SEP, 2)
    return read, int(i), int(j)


def run_contacts(tmp_path, records, *, mapq=1, sort_order="queryname", output=None,
                 fmt="pa5", contigs=CONTIGS, group_order=None):
    tmp_path.mkdir(parents=True, exist_ok=True)
    sam = tmp_path / "in.sam"
    write_sam(sam, records, sort_order, contigs, group_order)
    out = tmp_path / (output or f"out.{fmt}")
    result = reconstruct_contacts(str(sam), str(out), mapq, 1, fmt)
    rows = read_bed(out) if fmt == "bed" else read_pa5(out)
    return rows, result


def seg(read, k):
    return f"{read}{SEP}{k}"


# --- F: four mapped segments -> six contacts, each pair once -------------

def test_four_mapped_segments_give_six_contacts(tmp_path):
    records = [(seg("r1", k), 0, "ctg1", 1000 * k, 60, "500M") for k in (1, 2, 3, 4)]
    rows, result = run_contacts(tmp_path, records)

    pairs = [split_pair_name(row[0]) for row in rows]
    assert len(rows) == result.contacts_written == 6
    assert all(read == "r1" and i < j for read, i, j in pairs)
    assert sorted((i, j) for _, i, j in pairs) == [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
    assert len({row[0] for row in rows}) == 6
    assert result.reads_with_contacts == 1 and result.usable_segments == 4


@pytest.mark.parametrize("k,expected", [(2, 1), (3, 3), (5, 10), (17, 136)])
def test_contact_count_is_k_choose_2(tmp_path, k, expected):
    records = [(seg("r", i), 0, "ctg1", 100 * i, 60, "50M") for i in range(1, k + 1)]
    rows, result = run_contacts(tmp_path, records)

    assert len(rows) == result.contacts_written == expected
    assert result.max_usable_in_read == k and result.max_contacts_in_read == expected


def test_single_usable_segment_gives_no_contact(tmp_path):
    rows, result = run_contacts(tmp_path, [(seg("r", 1), 0, "ctg1", 100, 60, "50M")])

    assert rows == [] and result.contacts_written == 0
    assert result.reads_seen == 1 and result.reads_with_contacts == 0


# --- G: MAPQ threshold ----------------------------------------------------

def test_mapq_threshold_selects_the_usable_segments(tmp_path):
    records = [
        (seg("r", 1), 0, "ctg1", 100, 60, "50M"),
        (seg("r", 2), 0, "ctg1", 5000, 30, "50M"),
        (seg("r", 3), 0, "ctg2", 100, 0, "50M"),
        (seg("r", 4), 4, "*", 0, 0, "*"),
    ]
    rows, result = run_contacts(tmp_path, records, mapq=1)

    assert [split_pair_name(row[0])[1:] for row in rows] == [(1, 2)]
    assert rows[0][5:] == ["60", "30"]
    assert result.records_seen == 4 and result.segments_seen == 4 and result.reads_seen == 1
    assert result.primary_mapped == 3 and result.unmapped == 1
    assert result.below_mapq == 1 and result.usable_segments == 2
    assert result.reads_with_contacts == 1 and result.contacts_written == 1


def test_raising_the_threshold_drops_the_lower_segment(tmp_path):
    records = [
        (seg("r", 1), 0, "ctg1", 100, 60, "50M"),
        (seg("r", 2), 0, "ctg1", 5000, 30, "50M"),
        (seg("r", 3), 0, "ctg2", 100, 0, "50M"),
    ]
    rows, result = run_contacts(tmp_path, records, mapq=31)

    assert rows == []
    assert result.below_mapq == 2 and result.usable_segments == 1


def test_mapq_zero_threshold_keeps_everything_mapped(tmp_path):
    records = [
        (seg("r", 1), 0, "ctg1", 100, 0, "50M"),
        (seg("r", 2), 0, "ctg1", 5000, 0, "50M"),
    ]
    rows, _ = run_contacts(tmp_path, records, mapq=0)

    assert len(rows) == 1


# --- H: primary alignments only ------------------------------------------

def test_secondary_and_supplementary_records_are_ignored(tmp_path):
    records = [
        (seg("r", 1), 0, "ctg1", 1001, 60, "100M"),
        (seg("r", 1), 256, "ctg2", 7001, 0, "100M"),       # secondary
        (seg("r", 2), 0, "ctg1", 5001, 60, "60M40S"),
        (seg("r", 2), 2048, "ctg2", 9001, 60, "60H40M"),   # supplementary
        (seg("r", 3), 0, "ctg2", 2001, 60, "100M"),
    ]
    rows, result = run_contacts(tmp_path, records)

    assert len(rows) == 3
    assert result.secondary_ignored == 1 and result.supplementary_ignored == 1
    assert result.segments_seen == 3 and result.usable_segments == 3
    # the supplementary placement must not leak into any row
    assert not any("ctg2" in row and row[2] == "9020" for row in rows)
    by_pair = {split_pair_name(row[0])[1:]: row for row in rows}
    assert by_pair[(1, 2)][1:5] == ["ctg1", "1050", "ctg1", "5030"]


def test_secondary_only_segment_is_seen_but_unusable(tmp_path):
    records = [
        (seg("r", 1), 0, "ctg1", 100, 60, "50M"),
        (seg("r", 2), 256, "ctg1", 900, 0, "50M"),
        (seg("r", 3), 0, "ctg2", 100, 60, "50M"),
    ]
    rows, result = run_contacts(tmp_path, records)

    assert len(rows) == 1 and result.segments_seen == 3 and result.usable_segments == 2


def test_duplicate_primary_records_are_counted_once(tmp_path):
    records = [
        (seg("r", 1), 0, "ctg1", 100, 60, "50M"),
        (seg("r", 1), 0, "ctg2", 300, 60, "50M"),   # second primary for the same segment
        (seg("r", 2), 0, "ctg1", 900, 60, "50M"),
    ]
    rows, result = run_contacts(tmp_path, records)

    assert len(rows) == 1
    assert result.duplicate_primary == 1 and result.usable_segments == 2
    assert rows[0][1:3] == ["ctg1", "124"], "the first record stands"


# --- I: never across reads -------------------------------------------------

def test_segments_of_different_reads_are_never_paired(tmp_path):
    reads = ["m84039/1/ccs", "m84039/10/ccs", "m84039/1/ccs_extra", "other"]
    records = []
    for r in reads:
        records += [(seg(r, k), 0, "ctg1", 100 * k, 60, "50M") for k in (1, 2, 3)]
    rows, result = run_contacts(tmp_path, records)

    assert len(rows) == 3 * len(reads)
    assert result.reads_seen == len(reads) and result.reads_with_contacts == len(reads)
    for row in rows:
        read, i, j = split_pair_name(row[0])
        assert read in reads and i < j
        assert row[0] == seg(read, i) + SEP + str(j)
    assert {split_pair_name(row[0])[0] for row in rows} == set(reads)


def test_reads_with_shared_name_prefixes_stay_separate(tmp_path):
    """Grouping is by the exact read name, not by a prefix of it."""
    records = [
        (seg("r1", 1), 0, "ctg1", 100, 60, "50M"),
        (seg("r1", 2), 0, "ctg1", 900, 60, "50M"),
        (seg("r10", 1), 0, "ctg1", 100, 60, "50M"),
        (seg("r10", 2), 0, "ctg1", 900, 60, "50M"),
    ]
    rows, result = run_contacts(tmp_path, records)

    assert sorted(split_pair_name(row[0])[0] for row in rows) == ["r1", "r10"]
    assert result.reads_seen == 2


def test_regrouped_read_is_an_error_not_a_partial_result(tmp_path):
    """A read coming back after its group was flushed means the input is not grouped."""
    records = [
        (seg("a", 1), 0, "ctg1", 100, 60, "50M"),
        (seg("b", 1), 0, "ctg1", 100, 60, "50M"),
        (seg("a", 2), 0, "ctg1", 900, 60, "50M"),
    ]
    with pytest.raises(RuntimeError, match="samtools sort -n"):
        run_contacts(tmp_path, records, sort_order=None)


# --- J: coordinate semantics ----------------------------------------------

@pytest.mark.parametrize("pos,cigar,expected", [
    # 5S100M2I3D50M4S: reference length 153 -> 1000 + 153 // 2
    (1001, "5S100M2I3D50M4S", 1076),
    # odd/odd carry in yahs's s/2 + e/2 + (s&1 && e&1): 7 + 10 // 2
    (8, "10M", 12),
    (1, "11M", 5),
    # =, X and N consume reference like M and D; I, S, H do not
    (101, "10=5X10N5M3I2S", 115),
    (1, "1M", 0),
])
def test_pa5_position_is_the_zero_based_alignment_midpoint(tmp_path, pos, cigar, expected):
    records = [
        (seg("r", 1), 0, "ctg1", pos, 60, cigar),
        (seg("r", 2), 0, "ctg2", 1, 60, "1M"),
    ]
    rows, _ = run_contacts(tmp_path, records)

    assert rows == [[seg("r", 1) + SEP + "2", "ctg1", str(expected), "ctg2", "0", "60", "60"]]


def test_strand_does_not_change_the_position(tmp_path):
    fwd = [(seg("r", 1), 0, "ctg1", 1001, 60, "153M"), (seg("r", 2), 0, "ctg2", 1, 60, "1M")]
    rev = [(seg("r", 1), 16, "ctg1", 1001, 60, "153M"), (seg("r", 2), 0, "ctg2", 1, 60, "1M")]
    rows_fwd, _ = run_contacts(tmp_path / "f", fwd)
    rows_rev, _ = run_contacts(tmp_path / "r", rev)

    assert rows_fwd == rows_rev
    assert rows_fwd[0][2] == "1076"


def test_pa5_position_matches_yahs_formula():
    def yahs(s, e):
        return s // 2 + e // 2 + (1 if (s & 1 and e & 1) else 0)

    for s, e in [(0, 10), (0, 11), (1, 11), (1, 12), (99, 199), (99, 200), (5, 6), (5, 7),
                 (1000, 1153), (7, 17)]:
        assert pa5_position(s, e) == yahs(s, e) == (s + e) // 2


def test_rows_carry_both_mapqs_and_the_named_segment_order(tmp_path):
    records = [
        (seg("r", 3), 0, "ctg2", 1, 17, "20M"),
        (seg("r", 1), 0, "ctg1", 1, 42, "20M"),
    ]
    rows, _ = run_contacts(tmp_path, records)

    # sorted by span index whatever the input order, mapq columns follow suit
    assert rows == [[seg("r", 1) + SEP + "3", "ctg1", "10", "ctg2", "10", "42", "17"]]


# --- K: malformed names -----------------------------------------------------

@pytest.mark.parametrize("qname", [
    "plain_read",
    "read" + SEP,
    "read" + SEP + "abc",
    "read" + SEP + "0",
    SEP + "1",
])
def test_names_outside_the_contract_fail_loudly(tmp_path, qname):
    records = [(qname, 0, "ctg1", 100, 60, "50M")]
    with pytest.raises((RuntimeError, ValueError), match=SEP):
        run_contacts(tmp_path, records)


def test_a_bad_name_in_the_middle_of_a_valid_file_still_fails(tmp_path):
    records = [
        (seg("r", 1), 0, "ctg1", 100, 60, "50M"),
        ("m84039/1/ccs", 0, "ctg1", 100, 60, "50M"),
        (seg("r", 2), 0, "ctg1", 100, 60, "50M"),
    ]
    with pytest.raises((RuntimeError, ValueError)):
        run_contacts(tmp_path, records)


# --- input order and output format ----------------------------------------

def test_coordinate_sorted_input_is_refused(tmp_path):
    records = [(seg("r", k), 0, "ctg1", 100 * k, 60, "50M") for k in (1, 2)]
    with pytest.raises(RuntimeError, match="samtools sort -n"):
        run_contacts(tmp_path, records, sort_order="coordinate")


def test_missing_sort_order_streams_grouped_input(tmp_path):
    records = [(seg("r", k), 0, "ctg1", 100 * k, 60, "50M") for k in (1, 2, 3)]
    rows, result = run_contacts(tmp_path, records, sort_order=None)

    assert len(rows) == 3 and result.sort_order == ""


def test_gz_output_is_gzip(tmp_path):
    records = [(seg("r", k), 0, "ctg1", 100 * k, 60, "50M") for k in (1, 2, 3)]
    rows, _ = run_contacts(tmp_path, records, output="out.pa5.gz")

    assert (tmp_path / "out.pa5.gz").read_bytes()[:2] == b"\x1f\x8b"
    assert len(rows) == 3


def test_output_has_no_header_or_comment_lines(tmp_path):
    """yahs parses every non-blank line as a record, so none may be a comment."""
    records = [(seg("r", k), 0, "ctg1", 100 * k, 60, "50M") for k in (1, 2)]
    run_contacts(tmp_path, records)

    text = (tmp_path / "out.pa5").read_text()
    assert not text.startswith("#") and "\n#" not in text
    assert text.count("\n") == 1


def test_cli_writes_pa5_and_stats_json(tmp_path):
    sam = tmp_path / "in.sam"
    write_sam(sam, [
        (seg("r", 1), 0, "ctg1", 100, 60, "50M"),
        (seg("r", 2), 0, "ctg1", 900, 60, "50M"),
        (seg("r", 3), 0, "ctg2", 900, 0, "50M"),
        (seg("r", 4), 4, "*", 0, 0, "*"),
    ])
    out = tmp_path / "sample.pa5"
    proc = subprocess.run(
        [sys.executable, "-m", "cifi.cli", "contacts", str(sam), "-o", str(out), "-q", "1"],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stderr
    assert "Warning" not in proc.stderr, proc.stderr
    assert len(read_pa5(out)) == 1

    stats = json.loads((tmp_path / "sample_contacts_stats.json").read_text())
    res = stats["results"]
    assert res["records_seen"] == 4 and res["segments_seen"] == 4 and res["reads_seen"] == 1
    assert res["primary_mapped"] == 3 and res["unmapped"] == 1 and res["below_mapq"] == 1
    assert res["usable_segments"] == 2 and res["reads_with_contacts"] == 1
    assert res["contacts_written"] == 1
    assert res["pair_mates_equivalent"] == 4 * 3
    assert stats["parameters"]["mapq_threshold"] == 1
    assert (tmp_path / "sample_contacts_report.html").exists()


def test_cli_fails_clearly_on_coordinate_sorted_input(tmp_path):
    sam = tmp_path / "coord.sam"
    write_sam(sam, [(seg("r", 1), 0, "ctg1", 100, 60, "50M")], sort_order="coordinate")
    proc = subprocess.run(
        [sys.executable, "-m", "cifi.cli", "contacts", str(sam), "-o", str(tmp_path / "o.pa5")],
        capture_output=True, text=True,
    )
    assert proc.returncode == 1
    assert "samtools sort -n" in proc.stderr


def test_query_grouped_header_is_accepted_without_warning(tmp_path):
    """minimap2's own output declares SO:unsorted GO:query; that is grouped input."""
    sam = tmp_path / "grouped.sam"
    write_sam(sam, [(seg("r", k), 0, "ctg1", 100 * k, 60, "50M") for k in (1, 2)],
              sort_order="unsorted", group_order="query")
    proc = subprocess.run(
        [sys.executable, "-m", "cifi.cli", "contacts", str(sam), "-o", str(tmp_path / "o.pa5"),
         "--no-report", "--no-json"],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stderr
    assert "Warning" not in proc.stderr, proc.stderr


# --- L: BED output ----------------------------------------------------------
#
# yahs's BED reader (link.c, dump_links_from_bed_file) pairs each record with
# the next one when the names match, takes columns as contig, start, end,
# name, MAPQ, links the two midpoints and feeds the spans to its coverage
# normalisation. These tests pin the file to that reader.

def test_bed_writes_the_alignment_span_and_mapq_of_each_segment(tmp_path):
    records = [
        (seg("r", 1), 0, "chr1", 101, 60, "100M"),   # 0-based [100, 200)
        (seg("r", 2), 0, "chr2", 501, 35, "150M"),   # 0-based [500, 650)
    ]
    rows, result = run_contacts(tmp_path, records, fmt="bed",
                                contigs=(("chr1", 1000), ("chr2", 1000)))

    name = seg("r", 1) + SEP + "2"
    assert rows == [["chr1", "100", "200", name, "60"], ["chr2", "500", "650", name, "35"]]
    assert (tmp_path / "out.bed").read_text() == (
        f"chr1\t100\t200\t{name}\t60\nchr2\t500\t650\t{name}\t35\n")
    assert result.contacts_written == 1


def test_bed_spans_follow_each_alignment_not_a_read_length(tmp_path):
    """Segment lengths vary by orders of magnitude and every span is its own."""
    cigars = {1: "50M", 2: "1200M", 3: "7M", 4: "300M2I40D100M5S", 5: "10S2000M3I"}
    ref_len = {k: sum(int(n) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", c) if op in "MDN=X")
               for k, c in cigars.items()}
    records = [(seg("r", k), 0, "ctg1", 1 + 5000 * k, 60, c) for k, c in cigars.items()]
    rows, result = run_contacts(tmp_path, records, fmt="bed")

    assert result.contacts_written == 10 and len(rows) == 20
    spans = {}
    for name, a, b in bed_pairs(rows):
        _, i, j = split_pair_name(name)
        for k, rec in ((i, a), (j, b)):
            spans.setdefault(k, set()).add((int(rec[1]), int(rec[2])))
    assert spans == {k: {(5000 * k, 5000 * k + ref_len[k])} for k in cigars}
    assert sorted(e - s for k in cigars for s, e in spans[k]) == [7, 50, 440, 1200, 2000]


def test_bed_records_of_a_contact_are_adjacent(tmp_path):
    records = [(seg("r1", k), 0, "ctg1", 1000 * k, 60, "500M") for k in (1, 2, 3, 4)]
    rows, result = run_contacts(tmp_path, records, fmt="bed")

    assert result.contacts_written == 6 and len(rows) == 12
    pairs = bed_pairs(rows)
    assert [split_pair_name(name)[1:] for name, _, _ in pairs] == [
        (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
    assert len({name for name, _, _ in pairs}) == 6
    # a record never sits next to one of another contact
    for k in range(0, 12, 2):
        assert rows[k][3] == rows[k + 1][3]
        if k + 2 < 12:
            assert rows[k + 1][3] != rows[k + 2][3]


MIXED_RECORDS = [
    (seg("a", 1), 0, "ctg1", 1001, 60, "5S100M2I3D50M4S"),
    (seg("a", 2), 16, "ctg2", 8, 42, "11M"),
    (seg("a", 3), 0, "ctg1", 4001, 0, "50M"),            # below MAPQ 1
    (seg("a", 4), 4, "*", 0, 0, "*"),                     # unmapped
    (seg("a", 5), 0, "ctg2", 20001, 17, "10=5X10N5M3I2S"),
    (seg("b", 2), 0, "ctg1", 7, 60, "10M"),
    (seg("b", 2), 256, "ctg2", 7001, 0, "10M"),           # secondary
    (seg("b", 3), 0, "ctg1", 9001, 30, "60M40S"),
    (seg("b", 3), 2048, "ctg2", 9001, 60, "60H40M"),      # supplementary
    (seg("b", 7), 0, "ctg2", 1, 60, "1M"),
    (seg("c", 1), 0, "ctg1", 100, 60, "50M"),             # alone: no contact
]


def test_bed_midpoints_reproduce_the_pa5_rows(tmp_path):
    """Same BAM, both formats: contacts, contigs, midpoints and MAPQs agree."""
    pa5_rows, pa5_result = run_contacts(tmp_path / "pa5", MIXED_RECORDS, fmt="pa5")
    bed_rows, bed_result = run_contacts(tmp_path / "bed", MIXED_RECORDS, fmt="bed")

    assert pa5_rows and bed_as_pa5(bed_rows) == pa5_rows
    assert len(bed_rows) == 2 * len(pa5_rows) == 2 * 6
    # the midpoint yahs computes from the BED span is the PA5 value
    for name, a, b in bed_pairs(bed_rows):
        row = next(r for r in pa5_rows if r[0] == name)
        assert pa5_position(int(a[1]), int(a[2])) == int(row[2])
        assert pa5_position(int(b[1]), int(b[2])) == int(row[4])
    for field in ("records_seen", "segments_seen", "reads_seen", "primary_mapped", "unmapped",
                  "secondary_ignored", "supplementary_ignored", "duplicate_primary",
                  "below_mapq", "usable_segments", "reads_with_contacts", "contacts_written",
                  "pair_mates_equivalent", "max_usable_in_read", "max_contacts_in_read"):
        assert getattr(bed_result, field) == getattr(pa5_result, field), field


@pytest.mark.parametrize("mapq", [1, 31])
def test_filtering_is_the_same_in_both_formats(tmp_path, mapq):
    pa5_rows, pa5_result = run_contacts(tmp_path / "pa5", MIXED_RECORDS, fmt="pa5", mapq=mapq)
    bed_rows, bed_result = run_contacts(tmp_path / "bed", MIXED_RECORDS, fmt="bed", mapq=mapq)

    assert {name for name, _, _ in bed_pairs(bed_rows)} == {row[0] for row in pa5_rows}
    assert bed_result.below_mapq == pa5_result.below_mapq == (1 if mapq == 1 else 3)
    assert bed_result.unmapped == pa5_result.unmapped == 1
    assert bed_result.secondary_ignored == pa5_result.secondary_ignored == 1
    assert bed_result.supplementary_ignored == pa5_result.supplementary_ignored == 1
    assert bed_result.contacts_written == pa5_result.contacts_written == len(pa5_rows)
    # placements of the ignored records reach neither output
    excluded = {("ctg1", "4000"), ("ctg2", "7000"), ("ctg2", "9000")}
    assert not any((rec[0], rec[1]) in excluded for rec in bed_rows)
    assert not any((row[1], row[2]) in excluded or (row[3], row[4]) in excluded
                   for row in pa5_rows)
    assert all(int(rec[4]) >= mapq for rec in bed_rows)


def test_bed_gz_output_is_gzip(tmp_path):
    records = [(seg("r", k), 0, "ctg1", 100 * k, 60, "50M") for k in (1, 2, 3)]
    rows, _ = run_contacts(tmp_path, records, fmt="bed", output="out.bed.gz")

    assert (tmp_path / "out.bed.gz").read_bytes()[:2] == b"\x1f\x8b"
    assert len(rows) == 6


def test_unknown_format_is_rejected(tmp_path):
    sam = tmp_path / "in.sam"
    write_sam(sam, [(seg("r", 1), 0, "ctg1", 100, 60, "50M")])
    with pytest.raises(ValueError, match="pa5 or bed"):
        reconstruct_contacts(str(sam), str(tmp_path / "out.vcf"), 1, 1, "vcf")
    assert not (tmp_path / "out.vcf").exists()


def test_cli_infers_bed_from_the_output_name(tmp_path):
    sam = tmp_path / "in.sam"
    write_sam(sam, [(seg("r", k), 0, "ctg1", 100 * k, 60, "50M") for k in (1, 2)])
    out = tmp_path / "sample.bed"
    proc = subprocess.run(
        [sys.executable, "-m", "cifi.cli", "contacts", str(sam), "-o", str(out), "-q", "1"],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stderr
    assert "Warning" not in proc.stderr, proc.stderr
    assert len(read_bed(out)) == 2

    stats = json.loads((tmp_path / "sample_contacts_stats.json").read_text())
    assert stats["parameters"]["output_format"] == "bed"
    assert stats["results"]["contacts_written"] == 1
    assert (tmp_path / "sample_contacts_report.html").exists()


@pytest.mark.parametrize("output,flag,expected", [
    ("x.pa5", None, "pa5"),
    ("x.pa5.gz", None, "pa5"),
    ("x.bed.gz", None, "bed"),
    ("x.txt", "bed", "bed"),         # any name goes with an explicit format
    ("x.pa5", "bed", "bed"),         # explicit wins, with a warning
    ("x.bed", "PA5", "pa5"),
])
def test_cli_format_option_and_inference(tmp_path, output, flag, expected):
    sam = tmp_path / "in.sam"
    write_sam(sam, [(seg("r", k), 0, "ctg1", 100 * k, 60, "50M") for k in (1, 2)])
    out = tmp_path / output
    cmd = [sys.executable, "-m", "cifi.cli", "contacts", str(sam), "-o", str(out),
           "--no-report", "--quiet"]
    if flag:
        cmd += ["--format", flag]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr

    rows = read_bed(out) if expected == "bed" else read_pa5(out)
    assert len(rows) == (2 if expected == "bed" else 1)
    stats = json.loads((tmp_path / "x_contacts_stats.json").read_text())
    assert stats["parameters"]["output_format"] == expected
    named = "bed" if ".bed" in output else "pa5" if ".pa5" in output else None
    if named and named != expected:
        assert "yahs will read it as" in proc.stderr
    else:
        assert "Warning" not in proc.stderr, proc.stderr


@pytest.mark.parametrize("output", ["contacts.txt", "contacts.tsv", "contacts", "contacts.gz"])
def test_cli_does_not_guess_the_format_of_an_unknown_name(tmp_path, output):
    """No --format and no .bed/.pa5: a usage error, not a silently chosen PA5.

    yahs stops on such a name as well (yahs.c: "unknown link file format")
    unless given --file-type, so a guessed format would fail one step later.
    """
    sam = tmp_path / "in.sam"
    write_sam(sam, [(seg("r", k), 0, "ctg1", 100 * k, 60, "50M") for k in (1, 2)])
    out = tmp_path / output
    proc = subprocess.run(
        [sys.executable, "-m", "cifi.cli", "contacts", str(sam), "-o", str(out),
         "--no-report", "--quiet"],
        capture_output=True, text=True,
    )
    assert proc.returncode == 2, proc.stderr
    assert "--format" in proc.stderr and ".bed" in proc.stderr and ".pa5" in proc.stderr
    assert [p.name for p in tmp_path.iterdir()] == ["in.sam"], "nothing was written"


# --- atomic output ----------------------------------------------------------

# Read r is complete and flushed (its contact written) when the error comes.
BAD_NAME = [
    (seg("r", 1), 0, "ctg1", 100, 60, "50M"),
    (seg("r", 2), 0, "ctg1", 900, 60, "50M"),
    (seg("s", 1), 0, "ctg1", 100, 60, "50M"),
    ("m84039/1/ccs", 0, "ctg1", 100, 60, "50M"),   # outside the naming contract
]
REGROUPED = [
    (seg("a", 1), 0, "ctg1", 100, 60, "50M"),
    (seg("a", 2), 0, "ctg1", 900, 60, "50M"),
    (seg("b", 1), 0, "ctg1", 100, 60, "50M"),
    (seg("a", 3), 0, "ctg1", 1700, 60, "50M"),    # a comes back after its group
]
# Three reads of 150 segments: 33,525 contacts, more than the writer's 1 MB
# block, so records have reached the disk before the bad name arrives.
MANY_THEN_BAD = [(seg(f"big{r}", k), 0, "ctg1", 1 + 10 * k, 60, "50M")
                 for r in range(3) for k in range(1, 151)]
MANY_THEN_BAD.append(("m84039/1/ccs", 0, "ctg1", 100, 60, "50M"))


def test_the_large_case_flushes_before_the_error(tmp_path):
    rows, _ = run_contacts(tmp_path, MANY_THEN_BAD[:-1])
    assert len(rows) == 3 * 150 * 149 // 2
    assert (tmp_path / "out.pa5").stat().st_size > 1 << 20


@pytest.mark.parametrize("records", [BAD_NAME, REGROUPED, MANY_THEN_BAD],
                         ids=["bad_name", "regrouped", "flushed_then_bad"])
@pytest.mark.parametrize("output", ["out.pa5", "out.pa5.gz", "out.bed", "out.bed.gz"])
def test_a_failed_run_leaves_no_output_file(tmp_path, records, output):
    """Contacts were written before the error; none of them may surface."""
    fmt = "bed" if ".bed" in output else "pa5"
    with pytest.raises((RuntimeError, ValueError)):
        run_contacts(tmp_path, records, fmt=fmt, output=output, sort_order=None)

    assert not (tmp_path / output).exists()
    assert [p.name for p in tmp_path.iterdir()] == ["in.sam"], "no temporary left behind"


def test_a_completed_run_leaves_only_the_output(tmp_path):
    records = [(seg("r", k), 0, "ctg1", 100 * k, 60, "50M") for k in (1, 2)]
    run_contacts(tmp_path, records, fmt="bed")

    assert sorted(p.name for p in tmp_path.iterdir()) == ["in.sam", "out.bed"]


def test_a_failed_run_keeps_the_previous_output(tmp_path):
    out = tmp_path / "out.pa5"
    out.write_text("previous\n")
    with pytest.raises((RuntimeError, ValueError)):
        run_contacts(tmp_path, BAD_NAME)

    assert out.read_text() == "previous\n"


def test_concurrent_runs_on_one_output_name_do_not_share_a_temporary(tmp_path):
    """Each run writes its own temporary; whichever lands last is complete.

    Three reads of 150 segments with MAPQs cycling 0/10/40/60, so different
    thresholds give different files; a temporary shared between the runs
    would have them truncate and overwrite each other's bytes.
    """
    records = [(seg(f"big{r}", k), 0, "ctg1", 1 + 10 * k, (0, 10, 40, 60)[k % 4], "50M")
               for r in range(3) for k in range(1, 151)]
    sam = tmp_path / "in.sam"
    write_sam(sam, records)
    other = tmp_path / "out.pa5.tmp"
    other.write_text("temporary of another run\n")

    def cmd(mapq, out):
        return [sys.executable, "-m", "cifi.cli", "contacts", str(sam), "-o", str(out),
                "-q", str(mapq), "--no-report", "--no-json", "--quiet"]

    thresholds = (1, 11, 41)
    procs = [subprocess.Popen(cmd(q, tmp_path / "out.pa5")) for q in thresholds]
    assert [p.wait() for p in procs] == [0, 0, 0]
    for q in thresholds:
        subprocess.run(cmd(q, tmp_path / f"q{q}.pa5"), check=True)

    expected = {q: (tmp_path / f"q{q}.pa5").read_bytes() for q in thresholds}
    assert len(set(expected.values())) == 3
    assert (tmp_path / "out.pa5").read_bytes() in expected.values()
    assert other.read_text() == "temporary of another run\n"
    assert sorted(p.name for p in tmp_path.iterdir()) == [
        "in.sam", "out.pa5", "out.pa5.tmp", "q1.pa5", "q11.pa5", "q41.pa5"]


# --- yahs end to end --------------------------------------------------------

def write_reference(path, contigs, rng, width=60):
    """FASTA plus the .fai yahs opens beside it, written by hand."""
    fai, offset = [], 0
    with open(path, "w") as fa:
        for name, length in contigs:
            header = f">{name}\n"
            fa.write(header)
            offset += len(header)
            fai.append(f"{name}\t{length}\t{offset}\t{width}\t{width + 1}")
            seq = "".join(rng.choice("ACGT") for _ in range(length))
            for i in range(0, length, width):
                line = seq[i:i + width] + "\n"
                fa.write(line)
                offset += len(line)
    with open(f"{path}.fai", "w") as fh:
        fh.write("\n".join(fai) + "\n")


def synthetic_segments(contigs, rng, n_reads):
    """Segments of varying length, mostly near each other on one contig."""
    records = []
    for r in range(n_reads):
        cname, clen = contigs[rng.randrange(len(contigs))]
        anchor = rng.randrange(clen)
        for k in range(1, rng.choice([2, 3, 3, 4, 5, 6]) + 1):
            length = rng.randint(100, 3000)
            if rng.random() < 0.15:
                cname2, clen2 = contigs[rng.randrange(len(contigs))]
                pos = rng.randrange(1, clen2 - length)
                records.append((seg(f"read{r}", k), 0, cname2, pos, 60, f"{length}M"))
            else:
                pos = min(max(1, anchor + int(rng.gauss(0, 20000))), clen - length)
                flag = 16 if rng.random() < 0.5 else 0
                records.append((seg(f"read{r}", k), flag, cname, pos, 60, f"{length}M"))
    return records


@pytest.mark.skipif(shutil.which("yahs") is None, reason="yahs is not on PATH")
@pytest.mark.parametrize("fmt", ["bed", "pa5"])
def test_yahs_reads_the_contacts_and_runs_to_completion(tmp_path, fmt):
    """yahs pairs every BED record with its neighbour; no --read-length is passed."""
    rng = random.Random(11)
    contigs = (("ctg1", 300000), ("ctg2", 200000), ("ctg3", 150000))
    ref = tmp_path / "ref.fa"
    write_reference(ref, contigs, rng)
    rows, result = run_contacts(tmp_path, synthetic_segments(contigs, rng, 600),
                                fmt=fmt, contigs=contigs)

    proc = subprocess.run(
        [shutil.which("yahs"), "-o", str(tmp_path / f"yahs_{fmt}"), str(ref),
         str(tmp_path / f"out.{fmt}")],
        capture_output=True, text=True, cwd=tmp_path,
    )
    assert proc.returncode == 0, proc.stderr
    dumped = re.search(r"dumped (\d+) read pairs from (\d+) records", proc.stderr)
    assert dumped, proc.stderr
    pairs, records = int(dumped.group(1)), int(dumped.group(2))
    assert pairs == result.contacts_written > 500
    assert records == len(rows) == (2 * pairs if fmt == "bed" else pairs)
    assert (tmp_path / f"yahs_{fmt}_scaffolds_final.agp").exists()
    assert (tmp_path / f"yahs_{fmt}_scaffolds_final.fa").exists()

"""`cifi contacts`: pairwise contacts from mapped unique segments.

The input is a name-grouped BAM (SAM here; htslib reads both) of segments
named by `cifi digest --segments-out`. Every read with k usable segments
must give exactly C(k,2) PA5 rows, never a self pair, never a duplicate, and
never a pair across reads. Positions follow yahs's own name-sorted-BAM
convention: the 0-based alignment midpoint.
"""

import gzip
import json
import re
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


def write_sam(path, records, sort_order="queryname"):
    header = ["@HD\tVN:1.6" + (f"\tSO:{sort_order}" if sort_order else "")]
    header += [f"@SQ\tSN:{name}\tLN:{length}" for name, length in CONTIGS]
    with open(path, "w") as fh:
        fh.write("\n".join(header + [sam_record(*r) for r in records]) + "\n")


def read_pa5(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        rows = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]
    for row in rows:
        assert len(row) == 7, row
    return rows


def split_pair_name(name):
    read, i, j = name.rsplit(SEP, 2)
    return read, int(i), int(j)


def run_contacts(tmp_path, records, *, mapq=1, sort_order="queryname", output="out.pa5"):
    tmp_path.mkdir(parents=True, exist_ok=True)
    sam = tmp_path / "in.sam"
    write_sam(sam, records, sort_order)
    out = tmp_path / output
    result = reconstruct_contacts(str(sam), str(out), mapq, 1)
    return read_pa5(out), result


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

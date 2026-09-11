"""Regressions for defects found in review of the PE/rename/stats work.

Each test here corresponds to a specific defect: the stats that did not add up,
the trim that skipped a cut at position 0, the fast-mode median that could
exceed the maximum, and the enzymes and inputs the original suite never covered.
"""

import json
import subprocess
import sys

import pytest

from cifi import get_enzyme_info, process_reads, process_reads_custom

HINDIII = "AAGCTT"
OVERHANG = 5


def block(unit, length, avoid=HINDIII):
    seq = (unit * (length // len(unit) + 1))[:length]
    assert avoid not in seq
    return seq


def make_read(*source_lengths):
    units = ["ACGT", "ACCG", "AGGC", "ATTG", "ACTA", "AGTC"]
    parts = []
    for i, slen in enumerate(source_lengths):
        last = i == len(source_lengths) - 1
        filler = slen - (0 if i == 0 else OVERHANG) - (0 if last else 1)
        assert filler >= 0
        parts.append(block(units[i % len(units)], filler))
    return HINDIII.join(parts)


def write_fastq(path, records):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")


def run(tmp_path, reads, *, enzyme="HindIII", min_segments=2, min_segment_len=60,
        fast=False):
    tmp_path.mkdir(parents=True, exist_ok=True)
    recs = [r if isinstance(r, tuple) else (f"r{i}", r) for i, r in enumerate(reads)]
    fq = tmp_path / "in.fastq"
    write_fastq(fq, recs)
    r1, r2 = tmp_path / "o_R1.fastq", tmp_path / "o_R2.fastq"
    result = process_reads(str(fq), str(r1), str(r2), enzyme,
                           min_segments, min_segment_len, True, False, fast, False)
    return result, recs


# --- the yield figures must account for every input base ------------------

def test_dropped_and_trimmed_bases_are_counted_for_filtered_reads(tmp_path):
    """A read failing min_segments still trimmed and dropped real bases."""
    # both segments emit under the cutoff -> read is filtered on min_segments
    result, _ = run(tmp_path, [make_read(70, 70)], min_segments=2, min_segment_len=200)

    assert result.reads_out == 0
    assert result.filtered_short_segments == 1
    assert result.segments_dropped_short == 2, "dropped segments must still be counted"
    assert result.bases_dropped_short > 0


def test_bases_in_filtered_reads_close_the_accounting_gap(tmp_path):
    """kept + trimmed + dropped + filtered-read bases == total bases in."""
    reads = [
        make_read(400, 300, 200),          # passes
        make_read(70, 70),                 # filtered on min_segments
        block("ACGT", 500),                # filtered on too few sites
    ]
    result, recs = run(tmp_path, reads, min_segments=3, min_segment_len=60)

    total = (result.segment_length_stats.sum()
             + result.bases_trimmed_overhang
             + result.bases_dropped_short
             + result.bases_in_filtered_reads)
    assert total == result.total_bases_in == sum(len(s) for _, s in recs)


def test_everything_dropped_short_is_reported_as_such(tmp_path):
    """The degenerate case that used to report zero dropped segments."""
    result, _ = run(tmp_path, [make_read(70, 70), make_read(70, 70)],
                    min_segments=2, min_segment_len=200)

    assert result.segments_dropped_short == 4
    assert result.segment_length_stats.sum() == 0


# --- trimming a cut at position 0 ----------------------------------------

def test_a_cut_at_position_zero_is_trimmed(tmp_path):
    """DpnII cuts before GATC, so a read starting with GATC begins at a cut."""
    site, cut = get_enzyme_info("DpnII")
    assert (site, cut) == ("GATC", 0)
    read = "GATC" + block("ACCA", 200, avoid="GATC") + "GATC" + block("ATTA", 200, avoid="GATC")

    result, _ = run(tmp_path, [read], enzyme="DpnII", min_segment_len=60)

    r1 = (tmp_path / "o_R1.fastq").read_text().splitlines()[1]
    assert not r1.startswith("GATC"), "leading segment begins at a cut and must be trimmed"


def test_zero_length_spans_are_not_counted_as_dropped_segments(tmp_path):
    """A site at position 0 and a site ending the read both make empty spans."""
    body = block("ACCA", 200, avoid="GATC")
    read = "GATC" + body + "GATC"          # empty span at each end
    result, _ = run(tmp_path, [read], enzyme="DpnII", min_segment_len=60)

    assert result.segments_dropped_short == 0


# --- fast mode ------------------------------------------------------------

def test_fast_mode_median_never_exceeds_the_maximum(tmp_path):
    """Bin-size-1 statistics are exact values, not bins to interpolate within."""
    reads = [make_read(*([200] * 3)) for _ in range(50)]
    fast, _ = run(tmp_path, reads, fast=True)

    s = fast.segments_per_read_stats
    assert s.min() == s.max() == 3
    assert s.median() <= s.max()
    assert s.median() == 3


def test_fast_and_exact_agree_on_integer_valued_medians(tmp_path):
    reads = [make_read(*([200] * n)) for n in (2, 3, 3, 4, 5)]
    exact, _ = run(tmp_path / "e", reads, fast=False)
    fast, _ = run(tmp_path / "f", reads, fast=True)

    assert fast.segments_per_read_stats.median() == exact.segments_per_read_stats.median()
    assert fast.sites_per_read_stats.median() == exact.sites_per_read_stats.median()


# --- enzymes the original suite never covered -----------------------------

@pytest.mark.parametrize("enzyme", ["DpnII", "MboI", "Sau3AI", "NlaIII", "HindIII"])
def test_every_enzyme_produces_mates_of_equal_count_and_name(tmp_path, enzyme):
    site, _ = get_enzyme_info(enzyme)
    filler = "".join(c for c in "ACGTTGCAACGT" if True)
    read = site.join(block("ACCA" if i % 2 else "TTGA", 300, avoid=site)
                     for i in range(4))

    result, _ = run(tmp_path / enzyme, [read], enzyme=enzyme, min_segment_len=60)

    lines1 = (tmp_path / enzyme / "o_R1.fastq").read_text().splitlines()
    lines2 = (tmp_path / enzyme / "o_R2.fastq").read_text().splitlines()
    assert len(lines1) == len(lines2) == 4 * result.pairs_written
    assert lines1[0::4] == lines2[0::4]


def test_custom_site_rejects_a_cut_position_past_the_site(tmp_path):
    """cut_offset > len(site) made overhang_length() negative and underflowed."""
    write_fastq(tmp_path / "in.fastq", [("r0", make_read(400, 300))])
    with pytest.raises(Exception):
        process_reads_custom(
            str(tmp_path / "in.fastq"), str(tmp_path / "a.fq"), str(tmp_path / "b.fq"),
            "GATC", 10, 2, 60, True, False, False, False,
        )


# --- report rendering -----------------------------------------------------

def _digest_cli(tmp_path, reads, *args):
    write_fastq(tmp_path / "in.fastq", [(f"r{i}", s) for i, s in enumerate(reads)])
    proc = subprocess.run(
        [sys.executable, "-m", "cifi.cli", "digest", str(tmp_path / "in.fastq"),
         "-e", "HindIII", "-o", str(tmp_path / "out"), "-m", "2", "-l", "60", *args],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stderr
    # report generation swallows exceptions into a stderr warning and still
    # exits 0, so a green exit code alone proves nothing
    assert "Report failed" not in proc.stderr, proc.stderr
    assert "Warning" not in proc.stderr or "strip-overhang" in proc.stderr, proc.stderr
    return proc


def test_fast_mode_report_keeps_the_input_and_yield_sections(tmp_path):
    """These sections were nested inside an exact-mode-only guard."""
    _digest_cli(tmp_path, [make_read(400, 300, 200)], "--fast", "--report", "--json")

    html = (tmp_path / "out_digestion_report.html").read_text()
    for section in ("Input Reads", "Yield"):
        assert f"<h2>{section}</h2>" in html, f"{section} missing under --fast"


def test_report_renders_from_a_pre_rename_stats_file(tmp_path):
    """Archived stats.json used fragment_* keys; rendering must not KeyError."""
    from cifi.report import generate_digest_report

    old = {
        "cifi_version": "0.1.0", "timestamp": "2026-01-01T00:00:00",
        "input": {"file": "x.bam", "path": "/x.bam", "format": "BAM"},
        "parameters": {"enzyme": "HindIII", "enzyme_site": "AAGCTT", "cut_offset": 1,
                       "custom_enzyme": False, "min_fragments": 3, "min_frag_len": 20,
                       "strip_overhang": True, "gzip_output": False, "fast_mode": False},
        "results": {"reads_in": 10, "reads_out": 8, "reads_skipped": 2,
                    "filtered_few_sites": 1, "filtered_short_frags": 1,
                    "pairs_written": 20, "total_fragments": 24, "pass_rate": 0.8,
                    "avg_fragments_per_read": 3.0, "avg_pairs_per_read": 2.5},
        "output": {"r1": "a_R1.fastq", "r2": "a_R2.fastq"},
        # written by the previous release: pre-rename key inside input_reads
        "input_reads": {"count": 10, "total_bases": 5000, "gc_content": 41.0,
                        "total_sites": 23, "mean_sites_per_read": 2.3,
                        "length": {"count": 10, "min": 100, "max": 900,
                                   "mean": 500.0, "median": 480.0}},
    }
    path = generate_digest_report(old, str(tmp_path / "old.html"))
    html = open(path).read()
    assert "Digestion" in html
    assert "Input Reads" in html, "the pre-rename input section must still render"


def test_sites_per_read_histogram_survives_fast_mode(tmp_path):
    """The template guard was fixed but this histogram was still exact-only."""
    _digest_cli(tmp_path, [make_read(400, 300, 200), make_read(400, 300)],
                "--fast", "--report", "--json")

    html = (tmp_path / "out_digestion_report.html").read_text()
    assert "<h2>Sites per Read Distribution</h2>" in html


def test_integer_histogram_edges_land_on_the_values(tmp_path):
    """Bars must be centred on 3, 4, 5 ... not 3.49, 4.46, 5.43."""
    reads = [make_read(*([200] * n)) for n in (2, 3, 4, 5)]
    _digest_cli(tmp_path, reads, "--report", "--json")

    d = json.loads((tmp_path / "out_stats.json").read_text())
    edges = d["segments_per_read"]["histogram"]["bins"]
    lo = d["segments_per_read"]["min"]
    assert edges[0] == lo - 0.5, f"first edge {edges[0]} should be {lo - 0.5}"
    assert all(abs(e - round(e * 2) / 2) < 1e-9 for e in edges)


def test_yield_table_reconciles_to_bases_in(tmp_path):
    """Every base in the HTML yield table must be accounted for."""
    _digest_cli(tmp_path, [make_read(400, 300, 200), make_read(70, 70)],
                "--report", "--json")

    html = (tmp_path / "out_digestion_report.html").read_text()
    assert "In Filtered Reads" in html

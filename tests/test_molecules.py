"""`cifi molecules` (the long molecule table).

The long table has one row per (segment, alignment record) of a name-grouped
BAM of digested segments, with the digest's read coordinates joined in from
the segments table.

Most tests run on a synthetic dataset built here: a HindIII-free reference
of four contigs (the last carries a copy of a piece of the first, so some
segments have secondary placements), reads assembled from known pieces on
either strand with a few substitutions, pieces too short to keep, a piece
of junk that maps nowhere and a chimeric piece that aligns in two parts.
The reads go through `cifi digest` as an unaligned BAM with MM/ML tags,
then minimap2 and samtools sort -n. Needs minimap2 and samtools on PATH.
"""

import gzip
import json
import random
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from test_segments_out import header_metadata, read_table

from cifi import __version__ as CIFI_VERSION
from cifi import extract_molecules

SEP = "__CIFI_SEG__"
HINDIII = "AAGCTT"
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
_COMPLEMENT = str.maketrans("ACGT", "TGCA")

MINIMAP2, SAMTOOLS = shutil.which("minimap2"), shutil.which("samtools")
needs_tools = pytest.mark.skipif(MINIMAP2 is None or SAMTOOLS is None,
                                 reason="minimap2 and samtools are needed")

MOLECULE_COLUMNS = [
    "molecule_id", "span_index", "retained_index", "segment_count", "spans_total",
    "read_length", "read_start", "read_end", "seg_len", "ref", "ref_start", "ref_end",
    "strand", "mapq", "aln", "rank", "as", "nm", "qstart", "qend", "mlen", "blen", "cigar",
]
READ_COLUMNS = ("spans_total", "read_length", "read_start", "read_end")


def revcomp(seq):
    return seq.translate(_COMPLEMENT)[::-1]


def cifi(*args, check=True):
    proc = subprocess.run([sys.executable, "-m", "cifi.cli", *map(str, args)],
                          capture_output=True, text=True)
    if check:
        assert proc.returncode == 0, proc.stderr
    return proc


def run(cmd, **kw):
    return subprocess.run(cmd, check=True, capture_output=True, text=True, **kw)


# --- synthetic dataset ------------------------------------------------------

def without_site(seq, rng):
    seq = list(seq)
    text = "".join(seq)
    idx = text.find(HINDIII)
    while idx != -1:
        seq[idx + 2] = rng.choice("AT")
        text = "".join(seq)
        idx = text.find(HINDIII, idx + 1)
    return text


def make_reference(rng, length=40000):
    ref = {f"ctg{i}": without_site("".join(rng.choice("ACGT") for _ in range(length)), rng)
           for i in (1, 2, 3, 4)}
    # a 3 kb copy of ctg1 inside ctg4: the source of secondary placements
    piece = ref["ctg1"][10000:13000]
    ref["ctg4"] = ref["ctg4"][:20000] + piece + ref["ctg4"][23000:]
    return ref


def mutate(seq, rng, n_subs):
    out = list(seq)
    for pos in rng.sample(range(len(out)), n_subs):
        out[pos] = rng.choice([b for b in "ACGT" if b != out[pos]])
    return "".join(out)


def piece(ref, rng, contig, start, end, strand="+", subs=0):
    seq = ref[contig][start:end]
    seq = mutate(seq, rng, subs) if subs else seq
    return revcomp(seq) if strand == "-" else seq


# name -> list of pieces; each piece is a sequence with a label for the checks
def make_reads(ref, rng):
    junk = without_site("".join(rng.choice("ACGT") for _ in range(900)), rng)
    chimera = ref["ctg2"][30000:31500] + "GG" + ref["ctg3"][5000:6200]
    return {
        # four mapped pieces on both strands, a few substitutions
        "m84039_240101_000000_s1/1/ccs": [
            piece(ref, rng, "ctg1", 1000, 2500, "+", subs=3),
            piece(ref, rng, "ctg2", 8000, 8700, "-"),
            piece(ref, rng, "ctg1", 30000, 32200, "-", subs=5),
            piece(ref, rng, "ctg3", 100, 1300, "+", subs=1),
        ],
        # a 30 bp piece is dropped by --min-segment-len: span_index gap
        "m84039_240101_000000_s1/2/ccs": [
            piece(ref, rng, "ctg2", 20000, 21000, "+"),
            piece(ref, rng, "ctg2", 20500, 20530, "+"),
            piece(ref, rng, "ctg3", 20000, 22500, "-", subs=4),
            piece(ref, rng, "ctg4", 1000, 1800, "+"),
        ],
        # junk in the middle: an unmapped segment
        "m84039_240101_000000_s1/3/ccs": [
            piece(ref, rng, "ctg3", 30000, 31000, "+"),
            junk,
            piece(ref, rng, "ctg3", 35000, 36000, "-", subs=2),
        ],
        # pieces of the duplicated region: MAPQ 0 and, with -N 5, secondaries
        "m84039_240101_000000_s1/4/ccs": [
            piece(ref, rng, "ctg1", 10500, 12500, "+"),
            piece(ref, rng, "ctg2", 1000, 2000, "+"),
            piece(ref, rng, "ctg1", 11000, 12000, "-", subs=1),
        ],
        # a chimeric piece aligns as primary plus supplementary
        "m84039_240101_000000_s1/5/ccs": [
            piece(ref, rng, "ctg4", 30000, 31000, "+"),
            chimera,
            piece(ref, rng, "ctg1", 20000, 21000, "-"),
        ],
        # only two pieces: fails --min-segments 3, absent everywhere
        "m84039_240101_000000_s1/6/ccs": [
            piece(ref, rng, "ctg1", 5000, 6000, "+"),
            piece(ref, rng, "ctg2", 5000, 6000, "+"),
        ],
        # every piece too short except one: also filtered
        "m84039_240101_000000_s1/7/ccs": [
            piece(ref, rng, "ctg1", 5000, 5040, "+"),
            piece(ref, rng, "ctg2", 5000, 6000, "+"),
            piece(ref, rng, "ctg2", 7000, 7030, "+"),
        ],
    }


def mm_tags(read, rng):
    """MM/ML on a sample of the A, C and T bases: A+a., C+m?, T-a."""
    calls = {"A+a.": [], "C+m?": [], "T-a.": []}
    for pos, base in enumerate(read):
        if base == "A" and rng.random() < 0.05:
            calls["A+a."].append(pos)
        elif base == "C" and rng.random() < 0.3:
            calls["C+m?"].append(pos)
        elif base == "T" and rng.random() < 0.05:
            calls["T-a."].append(pos)
    items, mls = [], []
    for key, positions in calls.items():
        base = key[0]
        deltas, seen = [], 0
        for pos in positions:
            n = read.count(base, 0, pos)
            deltas.append(n - seen)
            seen = n + 1
            mls.append(rng.randrange(256))
        items.append(key + "," + ",".join(map(str, deltas)))
    return "MM:Z:" + ";".join(items) + ";", "ML:B:C," + ",".join(map(str, mls))


def write_fasta(path, ref):
    with open(path, "w") as fh:
        for name, seq in ref.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")


class Dataset:
    def __init__(self, work):
        rng = random.Random(2024)
        self.work = work
        self.ref = make_reference(rng)
        self.reads = {name: HINDIII.join(pieces) for name, pieces in make_reads(self.ref, rng).items()}
        self.fa = work / "ref.fa"
        write_fasta(self.fa, self.ref)
        run([SAMTOOLS, "faidx", str(self.fa)])
        sam = work / "reads.sam"
        with open(sam, "w") as fh:
            fh.write("@HD\tVN:1.6\n")
            for name, read in self.reads.items():
                mm, ml = mm_tags(read, rng)
                qual = "".join(chr(33 + ((i * 7) % 41)) for i in range(len(read)))
                fh.write(f"{name}\t4\t*\t0\t0\t*\t*\t0\t0\t{read}\t{qual}\t{mm}\t{ml}\n")
        self.bam = work / "reads.bam"
        run([SAMTOOLS, "view", "-b", "-o", str(self.bam), str(sam)])

        self.segments_fq = work / "dg.segments.fastq"
        self.table = work / "dg.segments.tsv.gz"
        cifi("digest", self.bam, "-e", "HindIII", "-o", work / "dg", "--no-report",
             "--segments-out", self.segments_fq, "--segments-table", self.table)
        self.digest_stats = json.loads((work / "dg_stats.json").read_text())

        common = ["-x", "map-hifi", "--no-hash-name", "-t", "2", str(self.fa), str(self.segments_fq)]
        self.ns_bam = self.map(work / "segs.ns.bam", ["--secondary=no"], common)
        self.sec_bam = self.map(work / "segs.sec.ns.bam", ["-N", "5", "--secondary=yes"], common)
        self.paf = work / "segs.paf"
        with open(self.paf, "w") as fh:
            subprocess.run([MINIMAP2, "-c", "--secondary=no", *common], check=True, stdout=fh,
                           stderr=subprocess.DEVNULL)

        self.long = work / "mol.tsv.gz"
        self.long_stats = self.molecules(self.ns_bam, self.long, "--table", self.table)
        self.long_notable = work / "mol.notable.tsv.gz"
        self.notable_stats = self.molecules(self.ns_bam, self.long_notable)

    def map(self, out, flags, common):
        sam = out.with_suffix(".sam")
        with open(sam, "w") as fh:
            subprocess.run([MINIMAP2, "-a", *flags, *common], check=True, stdout=fh,
                           stderr=subprocess.DEVNULL)
        run([SAMTOOLS, "sort", "-n", "-o", str(out), str(sam)])
        return out

    def molecules(self, bam, out, *extra):
        cifi("molecules", bam, "-o", out, *extra)
        return json.loads(Path(str(out).replace(".tsv.gz", "") + "_molecules_stats.json").read_text())


@pytest.fixture(scope="module")
def data(tmp_path_factory):
    if MINIMAP2 is None or SAMTOOLS is None:
        pytest.skip("minimap2 and samtools are needed")
    return Dataset(tmp_path_factory.mktemp("molecules"))


def read_fastq(path):
    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh]
    return {lines[i][1:]: lines[i + 1] for i in range(0, len(lines), 4)}


def read_sam(path):
    with open(path) as fh:
        return [ln.rstrip("\n").split("\t") for ln in fh if not ln.startswith("@")]


def by_molecule(rows):
    out = {}
    for row in rows:
        out.setdefault(row["molecule_id"], []).append(row)
    return out


# --- the long table on the synthetic dataset -------------------------------

@needs_tools
def test_the_dataset_has_what_the_tests_need(data):
    stats = data.digest_stats["results"]
    assert stats["reads_in"] == 7 and stats["reads_out"] == 5
    _, table = read_table(data.table)
    assert {r["molecule_id"].split("/")[1] for r in table} == {"1", "2", "3", "4", "5"}
    gaps = [r for r in table if r["molecule_id"].endswith("/2/ccs")]
    assert [int(r["span_index"]) for r in gaps] == [1, 3, 4]
    flags = [int(f[1]) for f in read_sam(data.ns_bam.with_suffix(".sam"))]
    assert any(f & 4 for f in flags), "the junk piece must be unmapped"
    assert any(f & 2048 for f in flags), "the chimeric piece must yield a supplementary"
    sec = [int(f[1]) for f in read_sam(data.sec_bam.with_suffix(".sam"))]
    assert any(f & 256 for f in sec), "the duplicated region must yield secondaries"


@needs_tools
def test_columns_and_header(data):
    header, rows = read_table(data.long)
    command = header_metadata(header)["command"]
    assert command.startswith("cifi molecules ") and f"--table {data.table}" in command
    assert header == [
        "##cifi_format=molecules",
        "##format_version=1",
        f"##tool_version={CIFI_VERSION}",
        f"##command={command}",
        "##canonical=true",
        "##read_coordinates=0-based-half-open",
        "##reference_coordinates=0-based-half-open",
        "##span_index=1-based-original-digest-span",
        "##retained_index=1-based-dense-retained-order",
        "##aln=P-primary,S-secondary,L-supplementary,U-unmapped",
        "##rank=0-primary,then-AS-descending",
        f"##segments_table={data.table}",
        "#columns: " + "\t".join(MOLECULE_COLUMNS),
    ]
    assert list(rows[0]) == MOLECULE_COLUMNS
    raw = data.long.read_bytes()
    assert raw[:2] == b"\x1f\x8b" and raw[12:14] == b"BC"


@needs_tools
def test_indices_and_totals_follow_the_segments_table(data):
    _, rows = read_table(data.long)
    _, table = read_table(data.table)
    segments = {(r["molecule_id"], r["span_index"]): r for r in table}
    fastq = read_fastq(data.segments_fq)

    assert {(r["molecule_id"], r["span_index"]) for r in rows} == set(segments)
    for mol, mrows in by_molecule(rows).items():
        spans = sorted({int(r["span_index"]) for r in mrows})
        for row in mrows:
            seg = segments[(mol, row["span_index"])]
            assert int(row["retained_index"]) == spans.index(int(row["span_index"])) + 1
            assert row["retained_index"] == seg["retained_index"]
            assert row["segment_count"] == seg["segments_kept"] == str(len(spans))
            assert row["spans_total"] == seg["spans_total"]
            assert (row["read_start"], row["read_end"]) == (seg["read_start"], seg["read_end"])
            assert row["read_length"] == seg["read_length"] == str(len(data.reads[mol]))
            assert int(row["seg_len"]) == int(seg["read_end"]) - int(seg["read_start"])
            assert int(row["seg_len"]) == len(fastq[f"{mol}{SEP}{row['span_index']}"])


@needs_tools
def test_rows_are_ordered_by_molecule_retained_index_and_rank(data):
    _, rows = read_table(data.long)
    # molecules in BAM order: the name-sorted BAM's order of first appearance
    names = run([SAMTOOLS, "view", str(data.ns_bam)]).stdout.splitlines()
    bam_molecules = []
    for ln in names:
        mol = ln.split("\t")[0].rsplit(SEP, 1)[0]
        if mol not in bam_molecules:
            bam_molecules.append(mol)
    seen = []
    for row in rows:
        if row["molecule_id"] not in seen:
            seen.append(row["molecule_id"])
    assert seen == bam_molecules
    for mol, mrows in by_molecule(rows).items():
        keys = [(int(r["retained_index"]), int(r["rank"])) for r in mrows]
        assert keys == sorted(keys), mol
        assert len(set(keys)) == len(keys)


@needs_tools
def test_unmapped_segment_is_a_u_row(data):
    _, rows = read_table(data.long)
    junk = [r for r in rows if r["molecule_id"].endswith("/3/ccs") and r["span_index"] == "2"]
    assert len(junk) == 1
    row = junk[0]
    assert row["aln"] == "U" and row["rank"] == "0"
    for col in ("ref", "ref_start", "ref_end", "strand", "mapq", "as", "nm", "qstart", "qend",
                "mlen", "blen"):
        assert row[col] == ".", col
    assert row["cigar"] == "*"
    assert row["seg_len"] != "."
    others = [r for r in rows if r["molecule_id"].endswith("/3/ccs") and r["span_index"] != "2"]
    assert {r["aln"] for r in others} == {"P"}


@needs_tools
def test_primary_and_supplementary_rows_reproduce_minimap2s_paf(data):
    """The 12 PAF columns regenerated from the table equal minimap2 -c's own."""
    _, rows = read_table(data.long)
    tlen = {name: int(length) for name, length, *_ in
            (ln.split("\t") for ln in open(f"{data.fa}.fai"))}
    regenerated = set()
    for r in rows:
        if r["aln"] not in ("P", "L"):
            continue
        regenerated.add((f"{r['molecule_id']}{SEP}{r['span_index']}", int(r["seg_len"]),
                         int(r["qstart"]), int(r["qend"]), r["strand"], r["ref"], tlen[r["ref"]],
                         int(r["ref_start"]), int(r["ref_end"]), int(r["mlen"]), int(r["blen"]),
                         int(r["mapq"])))
    minimap = set()
    for ln in open(data.paf):
        f = ln.rstrip("\n").split("\t")
        assert "tp:A:P" in f[12:]
        minimap.add((f[0], *map(int, f[1:4]), f[4], f[5], *map(int, f[6:12])))
    assert regenerated == minimap
    assert len(minimap) > 10
    # the chimeric piece: one P and one L row, both in the PAF
    chimera = [r for r in rows if r["molecule_id"].endswith("/5/ccs") and r["span_index"] == "2"]
    assert sorted(r["aln"] for r in chimera) == ["L", "P"]
    assert [r["rank"] for r in sorted(chimera, key=lambda r: r["aln"] != "P")] == ["0", "1"]


@needs_tools
def test_as_nm_and_cigar_are_the_records_own(data):
    _, rows = read_table(data.long)
    records = {}
    for f in read_sam(data.ns_bam.with_suffix(".sam")):
        tags = dict(t.split(":", 1) for t in f[11:])
        key = (f[0], int(f[1]) & 0x900, f[2], int(f[3]) - 1)
        records[key] = (tags.get("AS", ".").split(":")[-1], tags.get("NM", ".").split(":")[-1], f[5])
    for r in rows:
        if r["aln"] == "U":
            continue
        kind = {"P": 0, "S": 0x100, "L": 0x800}[r["aln"]]
        key = (f"{r['molecule_id']}{SEP}{r['span_index']}", kind, r["ref"], int(r["ref_start"]))
        assert (r["as"], r["nm"], r["cigar"]) == records[key], key


@needs_tools
def test_without_the_table_read_columns_are_dots_and_the_rest_agree(data):
    _, with_table = read_table(data.long)
    header, without = read_table(data.long_notable)
    assert header_metadata(header)["segments_table"] == "none"
    assert len(without) == len(with_table)
    for a, b in zip(with_table, without):
        assert tuple(b[col] for col in READ_COLUMNS) == (".",) * len(READ_COLUMNS)
        for col in MOLECULE_COLUMNS:
            if col not in READ_COLUMNS:
                assert a[col] == b[col], col


@needs_tools
def test_stats_json(data):
    stats = data.long_stats
    res = stats["results"]
    _, rows = read_table(data.long)
    assert res["molecules"] == 5
    assert res["segments"] == len({(r["molecule_id"], r["span_index"]) for r in rows})
    assert res["rows_written"] == len(rows)
    by_aln = {k: sum(1 for r in rows if r["aln"] == k) for k in "PSLU"}
    assert res["records"] == {"primary": by_aln["P"], "secondary": by_aln["S"],
                              "supplementary": by_aln["L"], "unmapped": by_aln["U"]}
    assert res["molecules_with_candidates"] == 1 and res["segments_with_candidates"] == 1
    assert res["table_segments_missing_from_bam"] == 0
    assert res["table_molecules_missing_from_bam"] == 0
    assert stats["parameters"]["candidates"] == "all"
    assert stats["parameters"]["table"] == str(data.table)
    assert stats["input"]["sort_order"] == "queryname"


# --- secondary alignments: rank and --candidates ---------------------------

@needs_tools
def test_secondaries_are_ranked_by_alignment_score(data, tmp_path):
    out = tmp_path / "sec.tsv.gz"
    stats = data.molecules(data.sec_bam, out, "--table", data.table)
    _, rows = read_table(out)
    dup = [r for r in rows if r["molecule_id"].endswith("/4/ccs")]
    with_secondary = {}
    for r in dup:
        with_secondary.setdefault(r["span_index"], []).append(r)
    ranked = {k: v for k, v in with_secondary.items() if len(v) > 1}
    assert ranked, "the duplicated-region segments must carry secondaries"
    for recs in ranked.values():
        assert recs[0]["aln"] == "P" and recs[0]["rank"] == "0"
        assert [int(r["rank"]) for r in recs] == list(range(len(recs)))
        assert {r["aln"] for r in recs[1:]} == {"S"}
        scores = [int(r["as"]) for r in recs[1:]]
        assert scores == sorted(scores, reverse=True)
        assert {r["ref"] for r in recs} == {"ctg1", "ctg4"}
    assert stats["results"]["records"]["secondary"] > 0
    assert stats["results"]["molecules_with_candidates"] >= 1

    primary = tmp_path / "primary.tsv.gz"
    pstats = data.molecules(data.sec_bam, primary, "--table", data.table, "--candidates", "primary")
    _, prows = read_table(primary)
    assert {r["aln"] for r in prows} <= {"P", "U"}
    assert [(r["molecule_id"], r["span_index"]) for r in prows] == \
        [(r["molecule_id"], r["span_index"]) for r in rows if r["aln"] in ("P", "U")]
    assert pstats["parameters"]["candidates"] == "primary"
    assert pstats["results"]["rows_dropped_candidates"] == len(rows) - len(prows) > 0


# --- hand-built inputs: rank tie-breaks, errors -----------------------------

CONTIGS = (("ctg1", 100000), ("ctg2", 50000))


def query_length(cigar):
    return sum(int(n) for n, op in CIGAR_RE.findall(cigar) if op in "MIS=X")


def sam_record(qname, flag, rname, pos, mapq, cigar, tags=()):
    if flag & 4:
        return "\t".join([qname, str(flag), "*", "0", "0", "*", "*", "0", "0", "ACGTACGT", "*"])
    seq = "A" * query_length(cigar)
    return "\t".join([qname, str(flag), rname, str(pos), str(mapq), cigar, "*", "0", "0",
                      seq, "*", *tags])


def write_sam(path, records, sort_order="queryname"):
    hd = "@HD\tVN:1.6" + (f"\tSO:{sort_order}" if sort_order else "")
    header = [hd] + [f"@SQ\tSN:{n}\tLN:{ln}" for n, ln in CONTIGS]
    path.write_text("\n".join(header + [sam_record(*r) for r in records]) + "\n")


def seg(read, k):
    return f"{read}{SEP}{k}"


def run_molecules(tmp_path, records, *, table=None, candidates="all", sort_order="queryname"):
    tmp_path.mkdir(parents=True, exist_ok=True)
    sam = tmp_path / "in.sam"
    write_sam(sam, records, sort_order)
    out = tmp_path / "out.tsv.gz"
    result = extract_molecules(str(sam), str(out), table or "", candidates, 1,
                               "cifi molecules test", "1.1.0-test")
    _, rows = read_table(out)
    return rows, result


def write_segments_table(path, rows):
    """rows: (molecule_id, span_index, retained_index, segments_kept, spans_total,
    read_start, read_end). The read length is the molecule's last read_end
    plus a tail of 7 bp."""
    cols = ["molecule_id", "span_index", "retained_index", "segments_kept", "spans_total",
            "read_length", "read_start", "read_end", "span_start", "span_end", "original_len",
            "processed_len", "trimmed_5p", "terminal", "enzyme", "cut_offset"]
    read_length = {}
    for mol, *_, e in rows:
        read_length[mol] = max(read_length.get(mol, 0), e + 7)
    with gzip.open(path, "wt") as fh:
        fh.write("##cifi_format=segments\n##format_version=1\n##tool_version=test\n"
                 "##command=test\n")
        fh.write("#columns: " + "\t".join(cols) + "\n")
        for mol, k, ri, kept, total, s, e in rows:
            fh.write("\t".join(map(str, [mol, k, ri, kept, total, read_length[mol], s, e, s, e,
                                         e - s, e - s, 0, "I", "HindIII", 1])) + "\n")


def test_rank_orders_by_score_then_position(tmp_path):
    records = [
        (seg("r", 1), 0, "ctg1", 1001, 60, "100M", ("AS:i:180", "NM:i:2")),
        (seg("r", 1), 256, "ctg2", 501, 0, "100M", ("AS:i:150",)),
        (seg("r", 1), 256, "ctg1", 9001, 0, "100M", ("AS:i:170",)),
        (seg("r", 1), 256, "ctg1", 5001, 0, "100M", ("AS:i:170",)),
        (seg("r", 1), 2048, "ctg2", 7001, 60, "60H40M", ("AS:i:80",)),
        (seg("r", 2), 0, "ctg1", 3001, 60, "50M", ("AS:i:100", "NM:i:0")),
    ]
    rows, result = run_molecules(tmp_path, records)

    first = [r for r in rows if r["span_index"] == "1"]
    assert [(r["aln"], r["rank"], r["ref"], r["ref_start"], r["as"]) for r in first] == [
        ("P", "0", "ctg1", "1000", "180"),
        ("S", "1", "ctg1", "5000", "170"),   # equal scores: lower position first
        ("S", "2", "ctg1", "9000", "170"),
        ("S", "3", "ctg2", "500", "150"),
        ("L", "4", "ctg2", "7000", "80"),
    ]
    assert [r["nm"] for r in first] == ["2", ".", ".", ".", "."]
    assert [r["retained_index"] for r in rows] == ["1"] * 5 + ["2"]
    assert {r["segment_count"] for r in rows} == {"2"}
    assert {r["spans_total"] for r in rows} == {"."}
    assert result.molecules == 1 and result.segments == 2 and result.rows_written == 6
    assert result.secondary == 3 and result.supplementary == 1 and result.primary_mapped == 2
    assert result.segments_with_candidates == 1 and result.molecules_with_candidates == 1


def test_rank_ties_break_on_end_strand_and_cigar_whatever_the_input_order(tmp_path):
    """Three secondaries with the same score, contig and start: the rank
    follows ref_end, then strand, then the CIGAR, so both BAM orders give
    the same table."""
    tied = [
        (seg("r", 1), 256, "ctg1", 5001, 0, "100M", ("AS:i:170",)),
        (seg("r", 1), 272, "ctg1", 5001, 0, "100M", ("AS:i:170",)),
        (seg("r", 1), 256, "ctg1", 5001, 0, "50M50S", ("AS:i:170",)),
        (seg("r", 1), 256, "ctg1", 5001, 0, "40M1I59M", ("AS:i:170",)),
    ]
    primary = (seg("r", 1), 0, "ctg1", 1001, 60, "100M", ("AS:i:180",))
    expected = [
        ("P", "0", "1000", "1100", "+", "100M"),
        ("S", "1", "5000", "5050", "+", "50M50S"),      # shortest reference interval
        ("S", "2", "5000", "5099", "+", "40M1I59M"),
        ("S", "3", "5000", "5100", "+", "100M"),        # same interval: forward first
        ("S", "4", "5000", "5100", "-", "100M"),
    ]
    tables = []
    for i, order in enumerate((tied, tied[::-1], tied[2:] + tied[:2])):
        rows, _ = run_molecules(tmp_path / str(i), [primary, *order])
        got = [(r["aln"], r["rank"], r["ref_start"], r["ref_end"], r["strand"], r["cigar"])
               for r in rows]
        assert got == expected, i
        tables.append((tmp_path / str(i) / "out.tsv.gz").read_bytes())
    assert tables[0] == tables[1] == tables[2]


def test_query_interval_and_lengths_follow_paf_conventions(tmp_path):
    records = [
        (seg("r", 1), 0, "ctg1", 1001, 60, "5S100M2I3D50M4S", ("NM:i:9",)),
        (seg("r", 2), 16, "ctg1", 2001, 60, "10S80M7S", ("NM:i:1",)),
        (seg("r", 3), 2048, "ctg2", 3001, 60, "20H30M5S"),
        (seg("r", 4), 0, "ctg2", 4001, 60, "10=2X8=1I9=", ("NM:i:3",)),
        # minimap2 counts ambiguous bases in NM and leaves them out of blen (nn)
        (seg("r", 5), 0, "ctg2", 5001, 60, "100M", ("NM:i:7", "nn:i:5")),
    ]
    rows, _ = run_molecules(tmp_path, records)

    got = {r["span_index"]: (r["seg_len"], r["qstart"], r["qend"], r["ref_start"], r["ref_end"],
                             r["mlen"], r["blen"], r["strand"]) for r in rows}
    # mlen = M - (NM - I - D): 150 - (9 - 2 - 3) = 146; blen = M + I + D = 155
    assert got["1"] == ("161", "5", "157", "1000", "1153", "146", "155", "+")
    # reverse strand: clips swap ends
    assert got["2"] == ("97", "7", "87", "2000", "2080", "79", "80", "-")
    # hard clips count in the length; no NM: mismatches taken as zero
    assert got["3"] == ("55", "20", "50", "3000", "3030", "30", "30", "+")
    # extended CIGAR: matches are the = bases
    assert got["4"] == ("30", "0", "30", "4000", "4029", "27", "30", "+")
    # nn: 5 of the 7 NM are ambiguous bases; mlen = 100 - 7, blen = 100 - 5
    assert got["5"] == ("100", "0", "100", "5000", "5100", "93", "95", "+")


def test_unmapped_only_segment_and_duplicate_primary(tmp_path):
    records = [
        (seg("r", 1), 4, "*", 0, 0, "*"),
        (seg("r", 3), 0, "ctg1", 1001, 60, "100M"),
        (seg("r", 3), 0, "ctg2", 1001, 60, "100M"),   # a second primary record
        (seg("r", 7), 0, "ctg1", 5001, 60, "100M"),
    ]
    rows, result = run_molecules(tmp_path, records)

    assert [(r["span_index"], r["retained_index"], r["aln"], r["rank"]) for r in rows] == [
        ("1", "1", "U", "0"), ("3", "2", "P", "0"), ("3", "2", "P", "1"), ("7", "3", "P", "0")]
    assert result.duplicate_primary == 1 and result.unmapped == 1
    assert {r["segment_count"] for r in rows} == {"3"}


def test_candidates_primary_keeps_unmapped_rows(tmp_path):
    records = [
        (seg("r", 1), 4, "*", 0, 0, "*"),
        (seg("r", 2), 0, "ctg1", 1001, 60, "100M"),
        (seg("r", 2), 256, "ctg2", 1001, 0, "100M"),
        (seg("r", 2), 2048, "ctg2", 3001, 60, "50H50M"),
    ]
    rows, result = run_molecules(tmp_path, records, candidates="primary")

    assert [(r["span_index"], r["aln"]) for r in rows] == [("1", "U"), ("2", "P")]
    assert result.rows_dropped_candidates == 2 and result.rows_written == 2


def test_the_table_fills_in_read_coordinates_and_missing_segments(tmp_path):
    table = tmp_path / "segments.tsv.gz"
    tmp_path.mkdir(exist_ok=True)
    write_segments_table(table, [
        ("r", 1, 1, 3, 4, 0, 100), ("r", 3, 2, 3, 4, 205, 305), ("r", 4, 3, 3, 4, 310, 410),
        ("s", 1, 1, 2, 2, 0, 100), ("s", 2, 2, 2, 2, 105, 205),
    ])
    records = [
        (seg("r", 1), 0, "ctg1", 1001, 60, "100M"),
        (seg("r", 4), 0, "ctg1", 5001, 60, "100M"),     # span 3 is not in the BAM
        (seg("s", 2), 0, "ctg2", 1001, 60, "100M"),
        (seg("s", 1), 0, "ctg2", 3001, 60, "100M"),
    ]
    rows, result = run_molecules(tmp_path, records, table=str(table))

    assert [(r["molecule_id"], r["span_index"], r["retained_index"], r["aln"]) for r in rows] == [
        ("r", "1", "1", "P"), ("r", "3", "2", "U"), ("r", "4", "3", "P"),
        ("s", "1", "1", "P"), ("s", "2", "2", "P")]
    r3 = rows[1]
    assert (r3["read_start"], r3["read_end"], r3["seg_len"], r3["spans_total"]) == \
        ("205", "305", "100", "4")
    assert [r["read_length"] for r in rows] == ["417"] * 3 + ["212"] * 2
    assert {r["segment_count"] for r in rows[:3]} == {"3"}
    assert result.table_segments_missing_from_bam == 1
    assert result.table_molecules_missing_from_bam == 0


def test_a_table_molecule_absent_from_the_bam_is_counted(tmp_path):
    table = tmp_path / "segments.tsv.gz"
    tmp_path.mkdir(exist_ok=True)
    write_segments_table(table, [("r", 1, 1, 2, 2, 0, 100), ("r", 2, 2, 2, 2, 105, 205),
                                 ("gone", 1, 1, 2, 2, 0, 100), ("gone", 2, 2, 2, 2, 105, 205)])
    records = [(seg("r", k), 0, "ctg1", 1001 * k, 60, "100M") for k in (1, 2)]
    rows, result = run_molecules(tmp_path, records, table=str(table))

    assert {r["molecule_id"] for r in rows} == {"r"}
    assert result.table_molecules_missing_from_bam == 1


@pytest.mark.parametrize("records,message", [
    ([(seg("r", 1), 0, "ctg1", 1001, 60, "100M"), (seg("r", 2), 0, "ctg1", 2001, 60, "100M"),
      (seg("r", 5), 0, "ctg1", 3001, 60, "100M")], seg("r", 5)),
    ([(seg("x", 1), 0, "ctg1", 1001, 60, "100M"), (seg("x", 2), 0, "ctg1", 2001, 60, "100M")],
     seg("x", 1)),
    ([(seg("r", 1), 0, "ctg1", 1001, 60, "100M"), (seg("r", 2), 0, "ctg1", 2001, 60, "90M")],
     seg("r", 2)),
], ids=["segment_not_in_table", "molecule_not_in_table", "length_differs"])
def test_table_mismatches_name_the_offending_segment(tmp_path, records, message):
    table = tmp_path / "segments.tsv.gz"
    tmp_path.mkdir(exist_ok=True)
    write_segments_table(table, [("r", 1, 1, 2, 2, 0, 100), ("r", 2, 2, 2, 2, 105, 205)])
    with pytest.raises(RuntimeError, match=re.escape(message)):
        run_molecules(tmp_path, records, table=str(table))
    assert not (tmp_path / "out.tsv.gz").exists()


def test_regrouped_and_coordinate_sorted_inputs_are_refused(tmp_path):
    regrouped = [
        (seg("a", 1), 0, "ctg1", 1001, 60, "100M"),
        (seg("b", 1), 0, "ctg1", 1001, 60, "100M"),
        (seg("a", 2), 0, "ctg1", 2001, 60, "100M"),
    ]
    with pytest.raises(RuntimeError, match="samtools sort -n"):
        run_molecules(tmp_path / "regrouped", regrouped, sort_order=None)
    with pytest.raises(RuntimeError, match="samtools sort -n"):
        run_molecules(tmp_path / "coordinate", regrouped[:1], sort_order="coordinate")
    assert not (tmp_path / "regrouped" / "out.tsv.gz").exists()


@pytest.mark.parametrize("qname", ["plain_read", "read" + SEP, "read" + SEP + "0"])
def test_names_outside_the_contract_fail(tmp_path, qname):
    with pytest.raises((RuntimeError, ValueError), match=SEP):
        run_molecules(tmp_path, [(qname, 0, "ctg1", 1001, 60, "100M")])


def test_cli_writes_the_table_and_stats(tmp_path):
    sam = tmp_path / "in.sam"
    write_sam(sam, [
        (seg("r", 1), 0, "ctg1", 1001, 60, "100M", ("AS:i:200", "NM:i:0")),
        (seg("r", 2), 4, "*", 0, 0, "*"),
        (seg("r", 3), 0, "ctg2", 1001, 60, "100M", ("AS:i:200", "NM:i:0")),
    ])
    out = tmp_path / "x.cifi-molecules.tsv.gz"
    proc = cifi("molecules", sam, "-o", out)
    assert "Warning" not in proc.stderr, proc.stderr
    header, rows = read_table(out)
    assert len(rows) == 3
    meta = header_metadata(header)
    assert meta["command"].startswith("cifi molecules ") and meta["segments_table"] == "none"
    stats = json.loads((tmp_path / "x.cifi-molecules_molecules_stats.json").read_text())
    # the statistics repeat the table's provenance lines
    assert stats["command"] == meta["command"]
    assert stats["cifi_version"] == meta["tool_version"]
    assert stats["format"] == {"name": "molecules", "version": 1}
    assert stats["input"]["path"] == str(sam) and stats["input"]["table"] is None
    assert stats["results"]["molecules"] == 1 and stats["results"]["segments"] == 3
    assert stats["results"]["records"] == {"primary": 2, "secondary": 0, "supplementary": 0,
                                           "unmapped": 1}
    assert stats["output"]["file"] == str(out)
    assert str(out) in proc.stdout

    proc = cifi("molecules", sam, "-o", tmp_path / "y.tsv.gz", "--no-json", "--quiet")
    assert proc.stdout == ""
    assert not (tmp_path / "y_molecules_stats.json").exists()


def test_a_failed_run_leaves_no_output(tmp_path):
    sam = tmp_path / "in.sam"
    write_sam(sam, [(seg("r", 1), 0, "ctg1", 1001, 60, "100M"), ("bad", 0, "ctg1", 1, 60, "10M")])
    proc = cifi("molecules", sam, "-o", tmp_path / "out.tsv.gz", check=False)
    assert proc.returncode == 1 and SEP in proc.stderr
    assert sorted(p.name for p in tmp_path.iterdir()) == ["in.sam"]


# --- the native table headers ----------------------------------------------

def read_table_plain(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        lines = [ln.rstrip("\n") for ln in fh]
    header = [ln for ln in lines if ln.startswith("#")]
    names = next(ln for ln in header if ln.startswith("#columns:"))[len("#columns:"):].split()
    rows = [dict(zip(names, ln.split("\t"))) for ln in lines if not ln.startswith("#")]
    return header, rows


NATIVE_HEADER_KEYS = {
    "segments": ["cifi_format", "format_version", "tool_version", "command", "read_coordinates",
                 "span_index", "retained_index", "terminal"],
    "molecules": ["cifi_format", "format_version", "tool_version", "command", "canonical",
                  "read_coordinates", "reference_coordinates", "span_index", "retained_index",
                  "aln", "rank", "segments_table"],
}


@needs_tools
def test_native_headers_parse_to_the_documented_keys(data):
    """Each native table's header is ##key=value lines in a fixed order,
    parsed here into a dict, then the #columns line."""
    tables = {"segments": data.table, "molecules": data.long}
    metas = {}
    for name, path in tables.items():
        header, rows = read_table_plain(path)
        assert rows, name
        meta = metas[name] = header_metadata(header)
        assert list(meta) == NATIVE_HEADER_KEYS[name], name
        assert meta["cifi_format"] == name and meta["format_version"] == "1"
        assert meta["tool_version"] == CIFI_VERSION
        assert meta["command"].startswith("cifi ")
        assert all(v == v.strip() and "\t" not in v for v in meta.values()), name
    seg, mol = metas["segments"], metas["molecules"]
    assert seg["span_index"] == mol["span_index"] == "1-based-original-digest-span"
    assert seg["retained_index"] == mol["retained_index"] == "1-based-dense-retained-order"
    assert mol["canonical"] == "true"
    assert mol["segments_table"] == str(data.table)

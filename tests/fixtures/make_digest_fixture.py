"""Build the digest regression input and its expected outputs.

The expected files under tests/fixtures/digest/ were produced by cifi 1.0.0,
before unique-segment output existed. test_backward_compat.py re-runs the same
option sets and expects byte-identical R1/R2 and matching statistics, so a
change to pairing, naming, orientation or trimming shows up here first.

Regenerate the input (and, with --expected, the outputs from the currently
installed cifi) with:

    python tests/fixtures/make_digest_fixture.py [--expected]
"""

import gzip
import random
import json
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "digest"

HINDIII = "AAGCTT"
GATC = "GATC"

# Each entry: (name, option list). The name is the fixture prefix.
OPTION_SETS = [
    ("hindiii_default", ["-e", "HindIII"]),
    ("hindiii_m2_revcomp_nostrip", ["-e", "HindIII", "-m", "2", "--revcomp-r2", "--no-strip-overhang"]),
    ("gatc_custom_l30", ["--site", "GATC", "--cut-pos", "0", "-m", "2", "-l", "30"]),
    ("hindiii_m2_gzip_fast", ["-e", "HindIII", "-m", "2", "--gzip", "--fast"]),
]



def _relative(path):
    """Express a path under the fixture directory relative to the repo root."""
    path = str(path)
    root = str(OUT)
    if path.startswith(root):
        return "tests/fixtures/digest" + path[len(root):]
    return path

def random_block(rng, length, avoid=(HINDIII, GATC)):
    """Random sequence free of the recognition sites used by the option sets."""
    while True:
        seq = "".join(rng.choice("ACGT") for _ in range(length))
        if not any(site in seq for site in avoid):
            return seq


def join_with(site, rng, lengths):
    return site.join(random_block(rng, n) for n in lengths)


def make_reads():
    rng = random.Random(20240911)
    reads = [
        ("m84039_240101_000000_s1/1/ccs", join_with(HINDIII, rng, [400, 300, 250, 200])),
        # the 59bp span (64 - 5 overhang) falls under the default cutoff
        ("m84039_240101_000000_s1/2/ccs", join_with(HINDIII, rng, [400, 59, 300])),
        ("m84039_240101_000000_s1/3/ccs", random_block(rng, 900)),
        ("m84039_240101_000000_s1/4/ccs", join_with(HINDIII, rng, [500, 450])),
        ("m84039_240101_000000_s1/5/ccs", join_with(HINDIII, rng, [3000, 120, 800, 61, 2200, 150])),
        # site at the very start and at the very end of the read
        ("m84039_240101_000000_s1/6/ccs", HINDIII + join_with(HINDIII, rng, [300, 300, 300]) + HINDIII),
        ("read:with:colons_and__underscores", join_with(GATC, rng, [200, 80, 35, 150, 400])),
        ("m84039_240101_000000_s1/8/ccs", join_with(GATC, rng, [1500, 1200]) + HINDIII + random_block(rng, 700)),
    ]
    out = []
    for name, seq in reads:
        # ramped quality so a shifted or reversed slice is detectable
        qual = "".join(chr(33 + ((i * 7) % 41)) for i in range(len(seq)))
        out.append((name, seq, qual))
    return out


def write_input(path):
    with open(path, "w") as fh:
        for name, seq, qual in make_reads():
            fh.write(f"@{name}\n{seq}\n+\n{qual}\n")


def write_expected(input_fastq):
    for prefix, opts in OPTION_SETS:
        work = OUT / "work"
        work.mkdir(exist_ok=True)
        subprocess.run(
            [sys.executable, "-m", "cifi.cli", "digest", str(input_fastq),
             "-o", str(work / prefix), "--no-report", *opts],
            check=True, capture_output=True, text=True,
        )
        for mate in ("R1", "R2"):
            gz = work / f"{prefix}_{mate}.fastq.gz"
            plain = work / f"{prefix}_{mate}.fastq"
            if gz.exists():
                with gzip.open(gz, "rt") as src, open(OUT / f"{prefix}_{mate}.fastq", "w") as dst:
                    shutil.copyfileobj(src, dst)
            else:
                shutil.copy(plain, OUT / f"{prefix}_{mate}.fastq")
        # The stats file records where it was generated; keep the fixture
        # portable by rewriting those paths relative to the repository.
        stats = json.loads((work / f"{prefix}_stats.json").read_text())
        stats["input"]["path"] = _relative(stats["input"]["path"])
        for key, value in stats.get("output", {}).items():
            stats["output"][key] = _relative(value)
        (OUT / f"{prefix}_stats.json").write_text(json.dumps(stats, indent=2) + "\n")
        shutil.rmtree(work)


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    input_fastq = OUT / "input.fastq"
    write_input(input_fastq)
    if "--expected" in sys.argv:
        write_expected(input_fastq)

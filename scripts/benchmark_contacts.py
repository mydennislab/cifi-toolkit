#!/usr/bin/env python
"""Time the segments route, step by step.

Every step runs as a child process and reports wall time and peak RSS, so
the cost of `cifi contacts` can be read next to the mapping it replaces.
Nothing large is bundled: point it at a real sample, or let it synthesise
one (tests/e2e/synth.py).

    python scripts/benchmark_contacts.py --workdir build/bench                 # synthetic
    python scripts/benchmark_contacts.py --workdir build/bench --synthetic-reads 20000
    python scripts/benchmark_contacts.py --workdir build/bench \\
        --reference asm.fa --reads sample.cifi.bam --threads 32
    python scripts/benchmark_contacts.py --workdir build/bench --bam sample.segments.ns.bam
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent


def timed(cmd, stdout=None):
    """Run cmd; return (wall seconds, peak RSS in MB) of the child."""
    start = time.perf_counter()
    proc = subprocess.Popen(cmd, stdout=stdout, stderr=subprocess.DEVNULL)
    _, status, usage = os.wait4(proc.pid, 0)
    wall = time.perf_counter() - start
    if os.waitstatus_to_exitcode(status) != 0:
        raise RuntimeError(f"failed: {' '.join(map(str, cmd))}")
    return wall, usage.ru_maxrss / 1024


def count_fastq_records(path):
    import gzip
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        return sum(1 for i, _ in enumerate(fh) if i % 4 == 0)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--bam", help="name-grouped segments BAM; benchmarks cifi contacts only")
    ap.add_argument("--reference", help="assembly FASTA (with --reads)")
    ap.add_argument("--reads", help="CiFi reads, BAM or FASTQ (with --reference)")
    ap.add_argument("--enzyme", default="HindIII")
    ap.add_argument("--min-segments", type=int, default=3)
    ap.add_argument("--mapq", type=int, default=1)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--synthetic-reads", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args(argv)

    workdir = Path(args.workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    rows = []
    cifi = [sys.executable, "-m", "cifi.cli"]

    bam = args.bam
    if bam is None:
        minimap2, samtools = shutil.which("minimap2"), shutil.which("samtools")
        if not (minimap2 and samtools):
            print("minimap2 and samtools are needed to build the segments BAM", file=sys.stderr)
            return 1
        reference, reads = args.reference, args.reads
        if not (reference and reads):
            sys.path.insert(0, str(REPO / "tests" / "e2e"))
            from synth import write_dataset
            reference, reads = write_dataset(workdir / "synthetic", args.synthetic_reads, args.seed)

        segments = workdir / "sample.segments.fastq.gz"
        wall, rss = timed(cifi + ["digest", str(reads), "-e", args.enzyme, "-m",
                                  str(args.min_segments), "-o", str(workdir / "sample"),
                                  "--gzip", "--segments-out", str(segments),
                                  "--no-report", "--no-json"], stdout=subprocess.DEVNULL)
        n_segments = count_fastq_records(segments)
        n_pairs = count_fastq_records(workdir / "sample_R1.fastq.gz")
        rows.append(("cifi digest --segments-out", wall, rss,
                     f"{n_segments:,} segments, {n_pairs:,} pairs"))

        sam = workdir / "sample.segments.sam"
        with open(sam, "w") as out:
            wall, rss = timed([minimap2, "-t", str(args.threads), "-ax", "map-hifi",
                               str(reference), str(segments)], stdout=out)
        rows.append(("minimap2 -ax map-hifi (segments)", wall, rss,
                     f"{n_segments:,} records; the pairs route would map {2 * n_pairs:,}"))

        bam = workdir / "sample.segments.ns.bam"
        wall, rss = timed([samtools, "sort", "-n", "-@", str(args.threads), "-o", str(bam),
                           str(sam)])
        rows.append(("samtools sort -n", wall, rss, ""))
        sam.unlink()

    pa5 = workdir / "sample.pa5"
    wall, rss = timed(cifi + ["contacts", str(bam), "-o", str(pa5), "-q", str(args.mapq),
                              "-t", str(args.threads), "--no-report", "--quiet"])
    stats = json.loads((workdir / "sample_contacts_stats.json").read_text())["results"]
    rows.append(("cifi contacts", wall, rss,
                 f"{stats['records_seen']:,} records in, {stats['contacts_written']:,} contacts out"))

    print(f"{'step':<36} {'wall (s)':>9} {'max RSS (MB)':>13}  notes")
    for step, wall, rss, note in rows:
        print(f"{step:<36} {wall:>9.2f} {rss:>13.0f}  {note}")
    print()
    print(f"segments seen:        {stats['segments_seen']:,}")
    print(f"usable segments:      {stats['usable_segments']:,}")
    print(f"contacts written:     {stats['contacts_written']:,}")
    print(f"R1/R2 mates avoided:  {stats['pair_mates_equivalent']:,} "
          f"({stats['mapping_work_reduction']:.1f}x fewer sequences to map)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

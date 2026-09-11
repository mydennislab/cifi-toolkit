"""Old route versus new route: the same contacts, from far fewer alignments.

Old route (what the R1/R2 pairs give when mapped):
    cifi digest -> R1/R2 -> minimap2 (each mate on its own) -> both mates
    primary, mapped and at or above the MAPQ threshold -> alignment midpoints

New route:
    cifi digest --segments-out -> minimap2 (each segment once)
    -> samtools sort -n -> cifi contacts -> PA5

Both are reduced to {(read, span_i, span_j): ((contig, pos, mapq) x 2)} and
compared exactly: identities, contigs, positions and MAPQs. Positions are
the 0-based alignment midpoints on both sides, computed here from the SAM
CIGAR independently of the C++ implementation.

One caveat lies in the aligner, not the toolkit: minimap2 breaks ties between
equally scoring placements with a hash of the query name, and the two routes
name the same sequence differently. On real data a small share of ambiguous
(low-MAPQ) segments therefore lands elsewhere, and the pairs route even
disagrees with itself, since it maps each segment once per pair under a new
name each time. With `--minimap2-args --no-hash-name` the tie-break no longer
depends on the name and the two routes agree exactly.

    python tests/e2e/equivalence.py --workdir build/e2e            # synthetic
    python tests/e2e/equivalence.py --reference asm.fa --reads reads.bam --workdir build/e2e
"""

import argparse
import gzip
import re
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SEP = "__CIFI_SEG__"
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def find_tools():
    """Paths of minimap2 and samtools, or None if either is missing."""
    minimap2, samtools = shutil.which("minimap2"), shutil.which("samtools")
    if minimap2 and samtools:
        return minimap2, samtools
    return None


def run(cmd, **kw):
    subprocess.run(cmd, check=True, **kw)


def reference_length(cigar):
    return sum(int(n) for n, op in CIGAR_RE.findall(cigar) if op in "MDN=X")


def midpoint(pos1, cigar):
    """PA5 position: 0-based start plus half the reference span, floored."""
    pos0 = pos1 - 1
    return pos0 + reference_length(cigar) // 2


def primary_alignments(sam_path):
    """{qname: (contig, midpoint, mapq)} for mapped primary records; unmapped -> None."""
    out = {}
    with open(sam_path) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            f = line.rstrip("\n").split("\t")
            qname, flag = f[0], int(f[1])
            if flag & (0x100 | 0x800):
                continue
            if flag & 0x4:
                out[qname] = None
            else:
                out[qname] = (f[2], midpoint(int(f[3]), f[5]), int(f[4]))
    return out


def fastq_names(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        return [ln[1:].rstrip("\n").split()[0] for i, ln in enumerate(fh) if i % 4 == 0]


def split_pair_name(name):
    """R1/R2 pair name <read>_<i>_<j-i-1> in kept-order indices -> (read, i, j)."""
    read, i, d = name.rsplit("_", 2)
    return read, int(i), int(i) + int(d) + 1


def split_contact_name(name):
    read, i, j = name.rsplit(SEP, 2)
    return read, int(i), int(j)


def canonical(a, b):
    return tuple(sorted([a, b]))


def minimap2_cmd(minimap2, threads, reference, query, extra=()):
    return [minimap2, "-t", str(threads), "-ax", "map-hifi", *extra, str(reference), str(query)]


def digest(reads, prefix, enzyme, min_segments, segments_out=None):
    cmd = [sys.executable, "-m", "cifi.cli", "digest", str(reads), "-e", enzyme,
           "-m", str(min_segments), "-o", str(prefix), "--no-report", "--no-json"]
    if segments_out:
        cmd += ["--segments-out", str(segments_out)]
    run(cmd, stdout=subprocess.DEVNULL)


def old_route(workdir, tools, reference, reads, enzyme, min_segments, mapq, threads,
              aligner_args=()):
    """Contacts as the R1/R2 route yields them, keyed by kept-order indices."""
    minimap2, _ = tools
    prefix = workdir / "old"
    digest(reads, prefix, enzyme, min_segments)
    r1, r2 = f"{prefix}_R1.fastq", f"{prefix}_R2.fastq"
    aln = {}
    for mate, fastq in (("R1", r1), ("R2", r2)):
        sam = workdir / f"old_{mate}.sam"
        with open(sam, "w") as out:
            run(minimap2_cmd(minimap2, threads, reference, fastq, aligner_args),
                stdout=out, stderr=subprocess.DEVNULL)
        aln[mate] = primary_alignments(sam)

    contacts = {}
    names = fastq_names(r1)
    for name in names:
        a, b = aln["R1"].get(name), aln["R2"].get(name)
        if a is None or b is None or min(a[2], b[2]) < mapq:
            continue
        contacts[split_pair_name(name)] = canonical(a, b)

    # The pairs route maps a segment once per pair it is in, each time under
    # another name. Placements that disagree between those copies are the
    # aligner breaking a tie by query name; they cannot agree with any single
    # placement, whichever route made it.
    placements = {}
    for name, a in aln["R1"].items():
        read, i, _ = split_pair_name(name)
        placements.setdefault((read, i), set()).add(a)
    for name, a in aln["R2"].items():
        read, _, j = split_pair_name(name)
        placements.setdefault((read, j), set()).add(a)
    inconsistent = sum(1 for p in placements.values() if len(p) > 1)
    return contacts, {"pairs": len(names), "mapped_records": 2 * len(names),
                      "segments": len(placements), "inconsistent_segments": inconsistent}


def new_route(workdir, tools, reference, reads, enzyme, min_segments, mapq, threads,
              aligner_args=()):
    """Contacts from cifi contacts, keyed by span indices, plus the span order per read."""
    minimap2, samtools = tools
    prefix = workdir / "new"
    segments = workdir / "new.segments.fastq.gz"
    digest(reads, prefix, enzyme, min_segments, segments_out=segments)

    bam = workdir / "new.segments.ns.bam"
    mm = subprocess.Popen(minimap2_cmd(minimap2, threads, reference, segments, aligner_args),
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    run([samtools, "sort", "-n", "-@", "2", "-o", str(bam)], stdin=mm.stdout)
    mm.stdout.close()
    if mm.wait() != 0:
        raise RuntimeError("minimap2 failed")

    pa5 = workdir / "new.pa5"
    run([sys.executable, "-m", "cifi.cli", "contacts", str(bam), "-o", str(pa5),
         "-q", str(mapq), "--no-report", "--no-json", "--quiet"])

    contacts = {}
    with open(pa5) as fh:
        for line in fh:
            name, c1, p1, c2, p2, q1, q2 = line.rstrip("\n").split("\t")
            contacts[split_contact_name(name)] = canonical((c1, int(p1), int(q1)),
                                                           (c2, int(p2), int(q2)))

    # kept-order index -> span index, per read, from the order segments were written
    spans = {}
    for qname in fastq_names(segments):
        read, k = split_contact_name(qname + SEP + "0")[:2]
        spans.setdefault(read, []).append(k)
    return contacts, spans, {"segments": sum(len(v) for v in spans.values()),
                             "mapped_records": sum(len(v) for v in spans.values())}


def compare(old, new, spans):
    """Translate the old keys to span indices and diff the two contact sets."""
    translated = {}
    for (read, i, j), value in old.items():
        order = spans[read]
        translated[(read, order[i], order[j])] = value
    only_old = sorted(set(translated) - set(new))
    only_new = sorted(set(new) - set(translated))
    differing = sorted(k for k in set(translated) & set(new) if translated[k] != new[k])
    return {
        "old_contacts": len(translated),
        "new_contacts": len(new),
        "shared_identical": len(set(translated) & set(new)) - len(differing),
        "only_old": only_old,
        "only_new": only_new,
        "differing": [(k, translated[k], new[k]) for k in differing],
    }


def run_equivalence(workdir, reference, reads, enzyme="HindIII", min_segments=3,
                    mapq=1, threads=4, aligner_args=()):
    tools = find_tools()
    if tools is None:
        raise RuntimeError("minimap2 and samtools are required")
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    old, old_work = old_route(workdir, tools, reference, reads, enzyme, min_segments, mapq,
                              threads, aligner_args)
    new, spans, new_work = new_route(workdir, tools, reference, reads, enzyme, min_segments,
                                     mapq, threads, aligner_args)
    report = compare(old, new, spans)
    report["old_mapped_records"] = old_work["mapped_records"]
    report["new_mapped_records"] = new_work["mapped_records"]
    report["old_segments"] = old_work["segments"]
    report["old_inconsistent_segments"] = old_work["inconsistent_segments"]
    report["identical"] = not (report["only_old"] or report["only_new"] or report["differing"])
    return report


def print_report(report):
    print(f"old route contacts:      {report['old_contacts']:,}")
    print(f"new route contacts:      {report['new_contacts']:,}")
    print(f"identical contacts:      {report['shared_identical']:,}")
    print(f"only in old route:       {len(report['only_old']):,}")
    print(f"only in new route:       {len(report['only_new']):,}")
    print(f"same pair, other values: {len(report['differing']):,}")
    old_n, new_n = report["old_mapped_records"], report["new_mapped_records"]
    ratio = old_n / new_n if new_n else 0
    print(f"records mapped, old:     {old_n:,} (two mates per pair)")
    print(f"records mapped, new:     {new_n:,} (each segment once)")
    print(f"mapping work reduction:  {ratio:.1f}x")
    print(f"segments the old route placed inconsistently between its own pairs: "
          f"{report['old_inconsistent_segments']:,} of {report['old_segments']:,}")
    for k in report["only_old"][:5]:
        print(f"  only old: {k}")
    for k in report["only_new"][:5]:
        print(f"  only new: {k}")
    for k, a, b in report["differing"][:5]:
        print(f"  differs:  {k} old={a} new={b}")
    print("RESULT: identical" if report["identical"] else "RESULT: DIFFERENT")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--reference", help="FASTA; synthetic when omitted")
    ap.add_argument("--reads", help="CiFi reads (BAM/FASTQ); synthetic when omitted")
    ap.add_argument("--enzyme", default="HindIII")
    ap.add_argument("--min-segments", type=int, default=3)
    ap.add_argument("--mapq", type=int, default=1)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--synthetic-reads", type=int, default=500)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--minimap2-args", default="",
                    help="extra minimap2 options for both routes, e.g. --no-hash-name")
    args = ap.parse_args(argv)

    if find_tools() is None:
        print("minimap2 and samtools not found on PATH; nothing to compare", file=sys.stderr)
        return 0

    workdir = Path(args.workdir)
    reference, reads = args.reference, args.reads
    if not (reference and reads):
        sys.path.insert(0, str(HERE))
        from synth import write_dataset
        reference, reads = write_dataset(workdir / "synthetic", args.synthetic_reads, args.seed)

    report = run_equivalence(workdir, reference, reads, args.enzyme, args.min_segments,
                             args.mapq, args.threads, args.minimap2_args.split())
    print_report(report)
    return 0 if report["identical"] else 1


if __name__ == "__main__":
    sys.exit(main())

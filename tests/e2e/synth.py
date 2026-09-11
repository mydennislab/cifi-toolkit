"""Synthetic reference and CiFi reads for end-to-end checks.

Reads are concatemers of genomic pieces joined by HindIII sites, in the
orientation and with the imperfections the real thing has: either strand,
a low substitution rate, pieces too short to keep, pieces from a duplicated
region (MAPQ 0) and pieces of junk that map nowhere. The reference is kept
free of HindIII sites so every cut the digest makes is a ligation junction.

    python tests/e2e/synth.py --out build/synth --reads 2000 --seed 7
"""

import argparse
import random
from pathlib import Path

HINDIII = "AAGCTT"
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def revcomp(seq):
    return seq.translate(_COMPLEMENT)[::-1]


def random_sequence(rng, length):
    return "".join(rng.choice("ACGT") for _ in range(length))


def without_site(seq, rng):
    """Break every HindIII site (on either strand it reads the same) by one base."""
    seq = list(seq)
    text = "".join(seq)
    idx = text.find(HINDIII)
    while idx != -1:
        seq[idx + 2] = rng.choice("AT")   # AAGCTT -> AAACTT / AATCTT
        text = "".join(seq)
        idx = text.find(HINDIII, idx + 1)
    return text


def make_reference(rng, n_contigs=4, length=60000, duplicated=3000):
    """Contigs of random sequence; the last one carries a copy of a piece of the first."""
    ref = {f"ctg{i + 1}": without_site(random_sequence(rng, length), rng)
           for i in range(n_contigs)}
    if n_contigs >= 2 and duplicated:
        src, dst = f"ctg1", f"ctg{n_contigs}"
        piece = ref[src][10000:10000 + duplicated]
        at = length // 2
        ref[dst] = ref[dst][:at] + piece + ref[dst][at + duplicated:]
    return ref


def mutate(seq, rng, rate):
    out = list(seq)
    for i in range(len(out)):
        if rng.random() < rate:
            out[i] = rng.choice([b for b in "ACGT" if b != out[i]])
    return "".join(out)


def make_reads(rng, ref, n_reads, sub_rate=0.005):
    """Yield (name, sequence); the sequence is pieces joined by HindIII sites."""
    names = list(ref)
    dup_contig = names[-1]
    dup_start = len(ref[dup_contig]) // 2
    for i in range(n_reads):
        k = rng.choice([1, 2, 3, 3, 4, 4, 5, 6, 8, 12, 17])
        pieces = []
        for _ in range(k):
            roll = rng.random()
            if roll < 0.05:
                piece = random_sequence(rng, rng.randint(200, 1500))      # junk, unmapped
            elif roll < 0.15:
                start = dup_start + rng.randint(0, 2000)                   # duplicated, MAPQ 0
                piece = ref[dup_contig][start:start + rng.randint(300, 900)]
            elif roll < 0.25:
                contig = rng.choice(names)                                 # too short to keep
                start = rng.randint(0, len(ref[contig]) - 100)
                piece = ref[contig][start:start + rng.randint(20, 59)]
            else:
                contig = rng.choice(names)
                length = rng.randint(80, 2500)
                start = rng.randint(0, len(ref[contig]) - length)
                piece = ref[contig][start:start + length]
            piece = mutate(piece, rng, sub_rate)
            if rng.random() < 0.5:
                piece = revcomp(piece)
            pieces.append(piece)
        if i % 50 == 7:
            name = f"read:{i}:with:colons_and__underscores"
        else:
            name = f"m84039_240101_000000_s1/{1000 + i}/ccs"
        yield name, HINDIII.join(pieces)


def write_fasta(path, ref):
    with open(path, "w") as fh:
        for name, seq in ref.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")


def write_fastq(path, reads):
    with open(path, "w") as fh:
        for name, seq in reads:
            fh.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")


def write_dataset(out_dir, n_reads=500, seed=7):
    """Write ref.fa and reads.fastq under out_dir; returns their paths."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(seed)
    ref = make_reference(rng)
    ref_path, reads_path = out_dir / "ref.fa", out_dir / "reads.fastq"
    write_fasta(ref_path, ref)
    write_fastq(reads_path, make_reads(rng, ref, n_reads))
    return ref_path, reads_path


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--reads", type=int, default=500)
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()
    ref_path, reads_path = write_dataset(args.out, args.reads, args.seed)
    print(ref_path)
    print(reads_path)

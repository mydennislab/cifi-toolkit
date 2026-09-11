"""Segment trimming, the length cutoff, and R2 orientation.

The overhang is a property of the *segment* - only segments that begin at a
cut site carry one - not of whichever slot (R1 or R2) the segment lands in. So
each segment is trimmed once, at extraction, and `-l` is then a guarantee about
the length of an emitted read rather than of the untrimmed source segment.
"""

import pytest

from cifi import process_reads

HINDIII_SITE = "AAGCTT"
HINDIII_CUT_OFFSET = 1
OVERHANG = len(HINDIII_SITE) - HINDIII_CUT_OFFSET  # 5, i.e. "AGCTT"

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCATGCA")


def revcomp(seq):
    return seq.translate(_COMPLEMENT)[::-1]


def block(unit, length):
    seq = (unit * (length // len(unit) + 1))[:length]
    assert HINDIII_SITE not in seq
    return seq


def make_read(*source_lengths):
    """Build a read whose HindIII segments have exactly `source_lengths` bp.

    These are *source* lengths, before any overhang trim. Every segment but the
    first opens with the 5bp overhang, so its emitted length is 5 bp shorter.
    """
    units = ["ACGT", "ACCG", "AGGC", "ATTG", "ACTA", "AGTC"]
    parts = []
    for i, flen in enumerate(source_lengths):
        last = i == len(source_lengths) - 1
        filler = flen - (0 if i == 0 else OVERHANG) - (0 if last else 1)
        assert filler >= 0, f"segment {i} too short to construct"
        parts.append(block(units[i % len(units)], filler))
    return HINDIII_SITE.join(parts)


def expected_emitted(sequence, min_emit_len, lead_trim=OVERHANG):
    """Replicate the intended model: trim at extraction, then filter on emitted length."""
    sites, start = [], 0
    while (idx := sequence.find(HINDIII_SITE, start)) != -1:
        sites.append(idx)
        start = idx + 1
    cuts = [0] + [p + HINDIII_CUT_OFFSET for p in sites] + [len(sequence)]
    out = []
    for a, b in zip(cuts, cuts[1:]):
        s = a + lead_trim if a > 0 else a   # only segments starting at a cut are trimmed
        if b > s and b - s >= min_emit_len:
            out.append(sequence[s:b])
    return out


def write_fastq(path, records):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")


def read_fastq(path):
    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh]
    return [(lines[i][1:], lines[i + 1]) for i in range(0, len(lines), 4)]


def run_digest(tmp_path, sequence, *, min_segments=2, min_segment_len=60,
               strip_overhang=True, revcomp_r2=False):
    tmp_path.mkdir(parents=True, exist_ok=True)
    fq = tmp_path / "in.fastq"
    write_fastq(fq, [("read1", sequence)])
    r1, r2 = tmp_path / "out_R1.fastq", tmp_path / "out_R2.fastq"
    result = process_reads(
        str(fq), str(r1), str(r2), "HindIII",
        min_segments, min_segment_len, strip_overhang, False, False, revcomp_r2,
    )
    return read_fastq(r1), read_fastq(r2), result


# --- the length guarantee -------------------------------------------------

@pytest.mark.parametrize("cutoff", [60, 70])
def test_no_emitted_read_falls_below_the_cutoff(tmp_path, cutoff):
    read = make_read(cutoff, cutoff + OVERHANG, cutoff + 1, 400, cutoff + OVERHANG - 1)
    r1, r2, _ = run_digest(tmp_path, read, min_segment_len=cutoff)

    assert r1, "expected surviving pairs"
    for name, seq in r1 + r2:
        assert len(seq) >= cutoff, f"{name} is {len(seq)}bp, below the {cutoff}bp cutoff"


def test_a_segment_at_a_cut_survives_at_exactly_the_cutoff(tmp_path):
    """source == cutoff + overhang emits exactly `cutoff` bp, so it must be kept."""
    r1, r2, result = run_digest(tmp_path, make_read(400, 60 + OVERHANG), min_segment_len=60)

    assert result.pairs_written == 1
    assert len(r2[0][1]) == 60


def test_a_segment_one_bp_short_of_the_cutoff_is_dropped(tmp_path):
    """Emitting 59bp against a 60bp cutoff: the segment must not survive."""
    read = make_read(400, 60 + OVERHANG - 1, 400)
    r1, r2, result = run_digest(tmp_path, read, min_segment_len=60)

    assert result.total_segments == 2, "the just-too-short segment should be gone"
    assert result.pairs_written == 1
    assert [len(seq) for _, seq in r1] == [400]


def test_a_segment_barely_longer_than_the_overhang_cannot_emit_a_1bp_read(tmp_path):
    """The reported 1bp read: a 6bp segment trimmed of its 5bp overhang."""
    r1, r2, result = run_digest(tmp_path, make_read(400, OVERHANG + 1, 400), min_segment_len=6)

    assert result.total_segments == 2, "the 6bp segment emits 1bp and must be dropped"
    for name, seq in r1 + r2:
        assert len(seq) >= 6, f"{name} is {len(seq)}bp"


# --- trimming is a property of the segment, not of the slot ---------------

def test_the_read_leading_segment_is_not_trimmed(tmp_path):
    """It does not begin at a cut, so it carries no overhang and keeps full length."""
    r1, _, result = run_digest(tmp_path, make_read(60, 400), min_segment_len=60)

    assert result.pairs_written == 1, "a 60bp leading segment must clear a 60bp cutoff"
    assert len(r1[0][1]) == 60


def test_a_segment_emits_the_same_bases_as_r1_and_as_r2(tmp_path):
    """Same segment, either slot, same sequence - that is what 'not slot-based' means."""
    read = make_read(400, 300, 200)
    r1, r2, _ = run_digest(tmp_path, read, min_segment_len=60)
    segs = expected_emitted(read, 60)
    assert len(segs) == 3

    # pairs are (0,1), (0,2), (1,2)
    assert [seq for _, seq in r1] == [segs[0], segs[0], segs[1]]
    assert [seq for _, seq in r2] == [segs[1], segs[2], segs[2]]
    # segment 1 appears as R2 of pair 0 and as R1 of pair 2, identically
    assert r2[0][1] == r1[2][1]


def test_no_strip_keeps_the_overhang_on_every_segment(tmp_path):
    read = make_read(400, 300)
    _, r2, _ = run_digest(tmp_path, read, min_segment_len=60, strip_overhang=False)
    segs = expected_emitted(read, 60, lead_trim=0)

    assert [seq for _, seq in r2] == [segs[1]]
    assert r2[0][1].startswith("AGCTT")


# --- R2 orientation -------------------------------------------------------

def test_r2_is_not_reverse_complemented_by_default(tmp_path):
    """Native concatemer orientation is the default, matching the Pore-C convention."""
    read = make_read(400, 300)
    _, r2, _ = run_digest(tmp_path, read, min_segment_len=60)

    assert r2[0][1] == expected_emitted(read, 60)[1]


def test_revcomp_r2_reverse_complements_after_stripping(tmp_path):
    """Opt-in revcomp composes with stripping rather than replacing it."""
    read = make_read(400, 300)
    _, r2, _ = run_digest(tmp_path, read, min_segment_len=60, revcomp_r2=True)

    assert r2[0][1] == revcomp(expected_emitted(read, 60)[1])


def test_revcomp_r2_composes_with_no_strip(tmp_path):
    """The two options are independent; both off, both on, or either alone."""
    read = make_read(400, 300)
    _, r2, _ = run_digest(tmp_path, read, min_segment_len=60,
                          strip_overhang=False, revcomp_r2=True)

    assert r2[0][1] == revcomp(expected_emitted(read, 60, lead_trim=0)[1])
    assert r2[0][1].endswith(revcomp("AGCTT"))


def test_revcomp_r2_leaves_r1_alone(tmp_path):
    read = make_read(400, 300)
    plain_r1, _, _ = run_digest(tmp_path / "a", read, min_segment_len=60)
    rc_r1, _, _ = run_digest(tmp_path / "b", read, min_segment_len=60, revcomp_r2=True)

    assert plain_r1 == rc_r1


def test_quality_is_reversed_alongside_the_sequence(tmp_path):
    """A revcomp'd read must carry its quality string reversed too."""
    read = make_read(400, 300)
    fq = tmp_path / "in.fastq"
    # ramped quality so a reversal is detectable
    quals = "".join(chr(33 + (i % 40)) for i in range(len(read)))
    with open(fq, "w") as fh:
        fh.write(f"@read1\n{read}\n+\n{quals}\n")
    r1p, r2p = tmp_path / "o_R1.fastq", tmp_path / "o_R2.fastq"
    process_reads(str(fq), str(r1p), str(r2p), "HindIII", 2, 60, True, False, False, True)

    with open(r2p) as fh:
        lines = [ln.rstrip("\n") for ln in fh]
    seq, qual = lines[1], lines[3]
    assert len(seq) == len(qual)

    # the emitted R2 quality must be the source slice reversed
    start = read.find(HINDIII_SITE) + HINDIII_CUT_OFFSET + OVERHANG
    assert qual == quals[start:start + len(seq)][::-1]

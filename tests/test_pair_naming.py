"""Mates must be pairable by name.

R1/R2 FASTQ carries no SAM flags, so the read name is the only thing that
identifies a pair. Appending "/1" and "/2" makes the two files disagree on every
pair: tools that group mates by name (samtools fastq, pairtools, most QC) see
two unrelated singletons instead of one pair.
"""

import pytest

from cifi import process_reads

HINDIII_SITE = "AAGCTT"
HINDIII_CUT_OFFSET = 1
HINDIII_OVERHANG = len(HINDIII_SITE) - HINDIII_CUT_OFFSET


def block(unit, length):
    seq = (unit * (length // len(unit) + 1))[:length]
    assert HINDIII_SITE not in seq
    return seq


def make_read(*segment_lengths):
    """Build a read whose HindIII segments have exactly `segment_lengths` bp."""
    units = ["ACGT", "ACCG", "AGGC", "ATTG", "ACTA", "AGTC"]
    parts = []
    for i, flen in enumerate(segment_lengths):
        last = i == len(segment_lengths) - 1
        filler = flen - (0 if i == 0 else HINDIII_OVERHANG) - (0 if last else 1)
        assert filler >= 0, f"segment {i} too short to construct"
        parts.append(block(units[i % len(units)], filler))
    return HINDIII_SITE.join(parts)


def write_fastq(path, records):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")


def read_names(path):
    with open(path) as fh:
        return [ln[1:].rstrip("\n") for i, ln in enumerate(fh) if i % 4 == 0]


def run_digest(tmp_path, sequence, *, read_name="read1", min_segments=2, min_segment_len=20):
    fq = tmp_path / "in.fastq"
    write_fastq(fq, [(read_name, sequence)])
    r1, r2 = tmp_path / "out_R1.fastq", tmp_path / "out_R2.fastq"
    result = process_reads(
        str(fq), str(r1), str(r2), "HindIII",
        min_segments, min_segment_len, True, False, False,
    )
    return read_names(r1), read_names(r2), result


def test_mates_share_an_identical_name(tmp_path):
    """The whole point: R1[i] and R2[i] must be groupable by name."""
    n1, n2, _ = run_digest(tmp_path, make_read(200, 200, 200))

    assert n1 == n2 != []


def test_names_carry_no_mate_suffix(tmp_path):
    """A trailing /1 or /2 is what made the two files disagree."""
    n1, n2, _ = run_digest(tmp_path, make_read(200, 200, 200))

    for name in n1 + n2:
        assert not name.endswith("/1")
        assert not name.endswith("/2")


def test_each_pair_gets_a_unique_name(tmp_path):
    """Segments repeat across pairs, so the pair - not the segment - is the identity."""
    n1, _, result = run_digest(tmp_path, make_read(*([200] * 6)))

    assert len(n1) == result.pairs_written == 15  # 6 choose 2
    assert len(set(n1)) == len(n1)


def test_names_stay_unique_across_multiple_source_reads(tmp_path):
    """Two reads in one file must not collide in the pair namespace."""
    fq = tmp_path / "in.fastq"
    write_fastq(fq, [
        ("readA", make_read(200, 200, 200)),
        ("readB", make_read(200, 200, 200)),
    ])
    r1, r2 = tmp_path / "out_R1.fastq", tmp_path / "out_R2.fastq"
    process_reads(str(fq), str(r1), str(r2), "HindIII", 2, 20, True, False, False)

    n1, n2 = read_names(r1), read_names(r2)
    assert n1 == n2
    assert len(n1) == 6  # 3 pairs from each read
    assert len(set(n1)) == 6


def test_name_keeps_the_source_read_as_a_prefix(tmp_path):
    """Provenance back to the originating concatemer must survive."""
    n1, _, _ = run_digest(tmp_path, make_read(200, 200), read_name="m84039/44242450/ccs")

    assert len(n1) == 1
    assert n1[0].startswith("m84039/44242450/ccs")


@pytest.mark.parametrize("n_segments,expected_pairs", [(2, 1), (3, 3), (4, 6), (5, 10)])
def test_name_count_tracks_pairs_written(tmp_path, n_segments, expected_pairs):
    """Neither file may gain or lose a name relative to the pair counter."""
    n1, n2, result = run_digest(tmp_path, make_read(*([200] * n_segments)))

    assert result.pairs_written == expected_pairs
    assert len(n1) == len(n2) == expected_pairs

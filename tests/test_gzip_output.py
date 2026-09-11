"""Gzip output must contain exactly the same records as plain output.

zlib's gzprintf() formats into a fixed internal buffer (8KB by default) and
returns 0 without writing anything when the result does not fit. CiFi segments
routinely run to several kb, so a writer built on gzprintf drops long records on
the floor - and because R1 and R2 have different length distributions, it drops
a different number from each file, leaving the two mates out of step.
"""

import gzip

import pytest

from cifi import process_reads

HINDIII_SITE = "AAGCTT"
HINDIII_CUT_OFFSET = 1
HINDIII_OVERHANG = len(HINDIII_SITE) - HINDIII_CUT_OFFSET

# A FASTQ record is name + seq + "+" + qual, so a sequence a little over half
# the 8KB buffer is already enough to overflow it.
GZPRINTF_BUFFER = 8192


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


def read_fastq(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        lines = [ln.rstrip("\n") for ln in fh]
    return [(lines[i][1:], lines[i + 1]) for i in range(0, len(lines), 4)]


def run_digest(tmp_path, sequence, *, gzip_output, min_segments=2, min_segment_len=20):
    tmp_path.mkdir(parents=True, exist_ok=True)
    fq = tmp_path / "in.fastq"
    write_fastq(fq, [("read1", sequence)])
    suffix = ".fastq.gz" if gzip_output else ".fastq"
    r1 = tmp_path / f"out_R1{suffix}"
    r2 = tmp_path / f"out_R2{suffix}"
    result = process_reads(
        str(fq), str(r1), str(r2), "HindIII",
        min_segments, min_segment_len, True, gzip_output, False,
    )
    return read_fastq(r1), read_fastq(r2), result


# segments comfortably either side of the gzprintf buffer
LONG_READ_SEGMENTS = (300, 6000, 400, 9000, 500)


def test_gzip_output_keeps_records_larger_than_the_gzprintf_buffer(tmp_path):
    """A single long record must survive; gzprintf silently drops it."""
    read = make_read(300, 6000)
    r1, r2, result = run_digest(tmp_path, read, gzip_output=True)

    assert result.pairs_written == 1
    assert len(r1) == 1, "the long R2 record was dropped by the gzip writer"
    assert len(r2) == 1
    assert len(r2[0][1]) == 6000 - HINDIII_OVERHANG


def test_gzip_r1_and_r2_have_equal_record_counts(tmp_path):
    """The reported symptom: the two mate files must not drift apart."""
    r1, r2, result = run_digest(tmp_path, make_read(*LONG_READ_SEGMENTS), gzip_output=True)

    assert len(r1) == len(r2), f"R1 has {len(r1)} records, R2 has {len(r2)}"
    assert len(r1) == result.pairs_written


def test_every_written_pair_reaches_both_files(tmp_path):
    """pairs_written is the contract: neither file may be short of it."""
    r1, r2, result = run_digest(tmp_path, make_read(*LONG_READ_SEGMENTS), gzip_output=True)

    assert result.pairs_written == 10  # 5 choose 2
    assert len(r1) == result.pairs_written
    assert len(r2) == result.pairs_written


def test_gzip_output_is_byte_identical_to_plain_output(tmp_path):
    """Compression must not change which records are emitted, or their content."""
    read = make_read(*LONG_READ_SEGMENTS)
    plain_r1, plain_r2, _ = run_digest(tmp_path / "plain", read, gzip_output=False)
    gz_r1, gz_r2, _ = run_digest(tmp_path / "gz", read, gzip_output=True)

    assert gz_r1 == plain_r1
    assert gz_r2 == plain_r2


@pytest.mark.parametrize("segment_len", [
    GZPRINTF_BUFFER // 2 - 100,   # record just under the buffer
    GZPRINTF_BUFFER // 2 + 100,   # record just over the buffer
    GZPRINTF_BUFFER,              # comfortably over
    GZPRINTF_BUFFER * 3,          # far over
])
def test_records_survive_across_the_buffer_boundary(tmp_path, segment_len):
    """Sweep the boundary so an off-by-one in the writer cannot hide."""
    r1, r2, result = run_digest(tmp_path, make_read(300, segment_len), gzip_output=True)

    assert len(r1) == len(r2) == result.pairs_written == 1
    assert len(r2[0][1]) == segment_len - HINDIII_OVERHANG


def test_plain_output_also_keeps_long_records(tmp_path):
    """Regression guard for the writer that was already correct."""
    r1, r2, result = run_digest(tmp_path, make_read(300, GZPRINTF_BUFFER * 2), gzip_output=False)

    assert len(r1) == len(r2) == result.pairs_written == 1

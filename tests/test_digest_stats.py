"""Statistics collected during digestion.

The digest report used to describe only what came out: segment lengths and pair
counts. These cover the input profile (what was fed in), the yield (how much
sequence survived, and what it was lost to), and the per-read distributions
behind the averages.

Populations differ by design and the tests pin that down: input stats cover
every read including ones that were later skipped, while segment and pair stats
cover only reads that passed the filters.
"""

import pytest

from cifi import process_reads

HINDIII_SITE = "AAGCTT"
HINDIII_CUT_OFFSET = 1
OVERHANG = len(HINDIII_SITE) - HINDIII_CUT_OFFSET  # 5


def block(unit, length):
    seq = (unit * (length // len(unit) + 1))[:length]
    assert HINDIII_SITE not in seq
    return seq


def make_read(*source_lengths):
    """A read whose HindIII segments have exactly `source_lengths` bp (pre-trim)."""
    units = ["ACGT", "ACCG", "AGGC", "ATTG", "ACTA", "AGTC"]
    parts = []
    for i, slen in enumerate(source_lengths):
        last = i == len(source_lengths) - 1
        filler = slen - (0 if i == 0 else OVERHANG) - (0 if last else 1)
        assert filler >= 0, f"segment {i} too short to construct"
        parts.append(block(units[i % len(units)], filler))
    return HINDIII_SITE.join(parts)


def write_fastq(path, records):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")


def read_lengths(path):
    with open(path) as fh:
        return [len(ln.rstrip("\n")) for i, ln in enumerate(fh) if i % 4 == 1]


def run(tmp_path, reads, *, min_segments=2, min_segment_len=60):
    """reads: list of sequences, or (name, sequence) tuples."""
    tmp_path.mkdir(parents=True, exist_ok=True)
    recs = [r if isinstance(r, tuple) else (f"read{i}", r) for i, r in enumerate(reads)]
    fq = tmp_path / "in.fastq"
    write_fastq(fq, recs)
    r1, r2 = tmp_path / "out_R1.fastq", tmp_path / "out_R2.fastq"
    result = process_reads(
        str(fq), str(r1), str(r2), "HindIII",
        min_segments, min_segment_len, True, False, False, False,
    )
    return result, read_lengths(r1), read_lengths(r2), recs


# --- input profile: covers every read, including skipped ones -------------

def test_total_bases_in_counts_every_input_read(tmp_path):
    """A read skipped for too few sites still consumed input sequence."""
    passing = make_read(400, 300, 200)
    skipped = block("ACGT", 500)          # no HindIII site at all -> skipped
    result, _, _, recs = run(tmp_path, [passing, skipped], min_segments=3)

    assert result.reads_in == 2
    assert result.reads_out == 1
    assert result.total_bases_in == len(passing) + len(skipped)


def test_gc_bases_in_matches_the_input_sequence(tmp_path):
    reads = [make_read(400, 300), block("ACGT", 200), block("GGCC", 120)]
    result, _, _, recs = run(tmp_path, reads)

    expected = sum(s.count("G") + s.count("C") for _, s in recs)
    assert result.gc_bases_in == expected


def test_read_length_stats_cover_all_reads_not_just_passing(tmp_path):
    passing = make_read(400, 300)
    skipped = block("ACGT", 90)
    result, _, _, _ = run(tmp_path, [passing, skipped])

    stats = result.read_length_stats
    assert stats.count() == 2
    assert stats.min() == min(len(passing), len(skipped))
    assert stats.max() == max(len(passing), len(skipped))


def test_total_sites_counts_sites_in_every_read(tmp_path):
    reads = [make_read(400, 300, 200), make_read(400, 300), block("ACGT", 200)]
    result, _, _, recs = run(tmp_path, reads)

    expected = sum(s.count(HINDIII_SITE) for _, s in recs)
    assert expected == 2 + 1 + 0
    assert result.total_sites == expected


# --- yield ----------------------------------------------------------------

def test_bases_out_match_the_written_files(tmp_path):
    result, r1_lens, r2_lens, _ = run(tmp_path, [make_read(400, 300, 200)])

    assert result.bases_out_r1 == sum(r1_lens)
    assert result.bases_out_r2 == sum(r2_lens)


def test_bases_trimmed_counts_only_segments_that_begin_at_a_cut(tmp_path):
    """3 segments -> 2 begin at a cut, so 2 x overhang bases come off."""
    result, _, _, _ = run(tmp_path, [make_read(400, 300, 200)])

    assert result.total_segments == 3
    assert result.bases_trimmed_overhang == 2 * OVERHANG


def test_bases_dropped_accounts_for_segments_below_the_cutoff(tmp_path):
    """The dropped segment emits 59bp against a 60bp cutoff."""
    result, _, _, _ = run(tmp_path, [make_read(400, 60 + OVERHANG - 1, 400)])

    assert result.segments_dropped_short == 1
    assert result.total_segments == 2
    assert result.bases_dropped_short == 60 - 1


def test_every_input_base_is_accounted_for(tmp_path):
    """Cuts tile the read end to end, so kept + trimmed + dropped == the read.

    An exact invariant: any base unaccounted for is a bug in the extraction
    bookkeeping, not a rounding difference.
    """
    passing = make_read(400, 60 + OVERHANG - 1, 300)
    result, _, _, _ = run(tmp_path, [passing])

    kept = result.segment_length_stats.sum()
    assert kept + result.bases_trimmed_overhang + result.bases_dropped_short == len(passing)
    assert result.total_bases_in == len(passing)


# --- per-read distributions ----------------------------------------------

def test_segments_per_read_distribution(tmp_path):
    reads = [make_read(*([200] * 3)), make_read(*([200] * 5))]
    result, _, _, _ = run(tmp_path, reads)

    stats = result.segments_per_read_stats
    assert stats.count() == 2
    assert stats.min() == 3
    assert stats.max() == 5
    assert result.total_segments == 8


@pytest.mark.parametrize("n,expected_pairs", [(2, 1), (3, 3), (4, 6), (5, 10)])
def test_pairs_per_read_distribution(tmp_path, n, expected_pairs):
    result, _, _, _ = run(tmp_path, [make_read(*([200] * n))])

    stats = result.pairs_per_read_stats
    assert stats.count() == 1
    assert stats.min() == stats.max() == expected_pairs
    assert result.pairs_written == expected_pairs


def test_skipped_reads_are_absent_from_segment_and_pair_distributions(tmp_path):
    """Only passing reads contribute; that is what makes the populations differ."""
    result, _, _, _ = run(tmp_path, [make_read(200, 200), block("ACGT", 300)])

    assert result.reads_in == 2
    assert result.reads_out == 1
    assert result.segments_per_read_stats.count() == 1
    assert result.pairs_per_read_stats.count() == 1
    assert result.read_length_stats.count() == 2


# --- filtering breakdown --------------------------------------------------

def test_filter_causes_are_counted_separately(tmp_path):
    no_sites = block("ACGT", 300)                    # too few sites
    too_short = make_read(70, 70)                    # segments below a 200bp cutoff
    ok = make_read(400, 400)
    result, _, _, _ = run(tmp_path, [no_sites, too_short, ok], min_segment_len=200)

    assert result.reads_in == 3
    assert result.reads_out == 1
    assert result.filtered_few_sites == 1
    assert result.filtered_short_segments == 1
    assert result.reads_skipped == 2

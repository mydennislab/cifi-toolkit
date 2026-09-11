"""End to end: the segments route reproduces the pairs route exactly.

Needs minimap2 and samtools on PATH; skipped otherwise. The reference and
reads are synthetic (tests/e2e/synth.py), so the check runs anywhere.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

from equivalence import find_tools, run_equivalence  # noqa: E402
from synth import write_dataset  # noqa: E402

pytestmark = pytest.mark.skipif(find_tools() is None,
                                reason="minimap2 and samtools are needed for the e2e check")


@pytest.fixture(scope="module")
def report(tmp_path_factory):
    work = tmp_path_factory.mktemp("e2e")
    reference, reads = write_dataset(work / "synthetic", n_reads=400, seed=7)
    return run_equivalence(work, reference, reads, mapq=1, threads=2)


def test_both_routes_yield_the_same_contacts(report):
    assert report["only_old"] == [], report["only_old"][:5]
    assert report["only_new"] == [], report["only_new"][:5]
    assert report["differing"] == [], report["differing"][:5]
    assert report["identical"]


def test_the_comparison_is_not_vacuous(report):
    assert report["old_contacts"] == report["new_contacts"] > 100
    assert report["shared_identical"] == report["new_contacts"]


def test_the_segments_route_maps_far_fewer_records(report):
    assert report["new_mapped_records"] < report["old_mapped_records"]
    # every pair costs two mates, every segment one record: n(n-1) vs n per read
    assert report["old_mapped_records"] / report["new_mapped_records"] > 2

"""Digest without --segments-out must match what cifi 1.0.0 produced.

The fixtures under tests/fixtures/digest/ were written by the release before
unique-segment output existed (see make_digest_fixture.py). Pair names,
pairing order, orientation, trimming and the statistics must all still agree.
"""

import gzip
import importlib.util
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

FIXTURES = Path(__file__).resolve().parent / "fixtures"
DIGEST = FIXTURES / "digest"


def _option_sets():
    spec = importlib.util.spec_from_file_location("make_digest_fixture",
                                                  FIXTURES / "make_digest_fixture.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.OPTION_SETS


def _read_text(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        return fh.read()


def _comparable(stats):
    """Everything in the stats file except what legitimately varies per run."""
    stats = json.loads(json.dumps(stats))
    for key in ("cifi_version", "timestamp", "output"):
        stats.pop(key, None)
    stats["input"].pop("path", None)
    return stats


@pytest.mark.parametrize("prefix,opts", _option_sets())
def test_digest_matches_pre_change_fixture(tmp_path, prefix, opts):
    shutil.copy(DIGEST / "input.fastq", tmp_path / "input.fastq")
    proc = subprocess.run(
        [sys.executable, "-m", "cifi.cli", "digest", str(tmp_path / "input.fastq"),
         "-o", str(tmp_path / prefix), "--no-report", *opts],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stderr

    for mate in ("R1", "R2"):
        produced = tmp_path / f"{prefix}_{mate}.fastq"
        if not produced.exists():
            produced = tmp_path / f"{prefix}_{mate}.fastq.gz"
        assert _read_text(produced) == (DIGEST / f"{prefix}_{mate}.fastq").read_text(), mate

    produced = json.loads((tmp_path / f"{prefix}_stats.json").read_text())
    expected = json.loads((DIGEST / f"{prefix}_stats.json").read_text())
    assert _comparable(produced) == _comparable(expected)
    assert "segments_written" not in produced["results"]
    assert "segments" not in produced["output"]

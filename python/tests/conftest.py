"""Shared pytest configuration for the somatic_workflow python suite.

Responsibilities:
  * put the python/ directory on sys.path so the CLI scripts import as
    top-level modules (``import genotype``, ``import merge_pileups``, ...);
  * expose the ``--update-golden`` flag used by the regression tier;
  * provide common fixtures (output dir, reference .dict, golden-file helper).
"""

import pathlib
import sys

import pytest

# python/ — the directory holding the scripts under test (parent of tests/).
PYDIR = pathlib.Path(__file__).resolve().parent.parent
if str(PYDIR) not in sys.path:
    sys.path.insert(0, str(PYDIR))

# tests/ — so test modules can ``from helpers import ...`` under importlib mode.
TESTDIR = pathlib.Path(__file__).resolve().parent
if str(TESTDIR) not in sys.path:
    sys.path.insert(0, str(TESTDIR))

FIXTURES = pathlib.Path(__file__).resolve().parent / "fixtures"
GOLDEN = pathlib.Path(__file__).resolve().parent / "golden"


def pytest_addoption(parser):
    parser.addoption(
        "--update-golden",
        action="store_true",
        default=False,
        help="Rewrite golden baselines instead of asserting against them.",
    )


@pytest.fixture
def update_golden(request):
    return request.config.getoption("--update-golden")


@pytest.fixture
def tmp_out(tmp_path):
    """A fresh output directory for a script run."""
    out = tmp_path / "out"
    out.mkdir()
    return out


@pytest.fixture(scope="session")
def fixtures_dir():
    return FIXTURES


def _build_synth(out, chr_prefixed):
    import json
    import types

    import pandas as pd

    from tests.synthetic import build_cohort, generate_cohort, write_truth

    cohort = build_cohort()
    write_truth(cohort, out)
    files = generate_cohort(cohort, out, seed=0, chr_prefixed=chr_prefixed)
    truth = {
        name: pd.read_csv(out / "truth" / f"{name}.truth.tsv", sep="\t")
        for name in ["patient", "clones", "segments", "germline_sites",
                     "somatic_variants"]
    }
    manifest = json.loads((out / "manifest.json").read_text())
    return types.SimpleNamespace(dir=out, files=files, truth=truth,
                                 manifest=manifest, cohort=cohort,
                                 chr_prefixed=chr_prefixed)


@pytest.fixture(scope="session")
def synth(tmp_path_factory):
    """The synthetic golden cohort (canonical, no chr prefix), generated once.

    Namespace: dir (Path), files (index), truth (dict of DataFrames), manifest, cohort.
    """
    return _build_synth(tmp_path_factory.mktemp("synth"), chr_prefixed=False)


@pytest.fixture(scope="session")
def synth_chr(tmp_path_factory):
    """The synthetic cohort emitted with chr-prefixed contig naming."""
    return _build_synth(tmp_path_factory.mktemp("synth_chr"), chr_prefixed=True)


@pytest.fixture(scope="session")
def ref_dict(tmp_path_factory):
    """A tiny 3-contig reference .dict (GRCh37-style, no ``chr`` prefix)."""
    d = tmp_path_factory.mktemp("ref") / "tiny.dict"
    d.write_text(
        "@HD\tVN:1.6\n"
        "@SQ\tSN:1\tLN:1000\n"
        "@SQ\tSN:2\tLN:1000\n"
        "@SQ\tSN:X\tLN:1000\n"
    )
    return str(d)

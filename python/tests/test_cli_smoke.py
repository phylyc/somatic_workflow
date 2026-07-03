"""Smoke tests: every script imports and its argparse wiring works.

These run the scripts as subprocesses (exercising ``if __name__ == '__main__'``)
and assert ``--help`` exits 0. Cheap guard against import/argparse breakage.
"""

import pytest

from helpers import run_script

SCRIPTS = [
    "merge_pileups.py",
    "pileup_to_allelic_counts.py",
    "acs_conversion.py",
    "validate_inputs.py",
    "harmonize_copy_ratios.py",
    "genotype.py",
    "calculate_cancer_cell_fraction.py",
    "map_to_absolute_copy_number.py",
]


@pytest.mark.parametrize("script", SCRIPTS)
def test_help_exits_zero(script):
    proc = run_script(script, ["--help"], expect_returncode=0)
    assert proc.stdout, f"{script} --help produced no output"


@pytest.mark.parametrize("script", SCRIPTS)
def test_module_imports(script):
    """Importing the module must not crash or have side effects."""
    import importlib

    mod = importlib.import_module(script[:-3])
    assert hasattr(mod, "main"), f"{script} has no main()"

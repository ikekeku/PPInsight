"""Tests for the ppinsight umbrella CLI."""

import subprocess
import sys


def _run_cli(*args):
    """Run the ppinsight CLI in a subprocess and return the result."""
    return subprocess.run(
        [sys.executable, "-m", "ppinsight.cli", *args],
        capture_output=True, text=True, timeout=10,
    )


class TestUmbrellaCLI:
    def test_help(self):
        result = _run_cli("--help")
        assert result.returncode == 0
        assert "ppinsight" in result.stdout
        assert "Commands:" in result.stdout

    def test_no_args_shows_help(self):
        result = _run_cli()
        assert result.returncode == 0
        assert "Commands:" in result.stdout

    def test_unknown_command(self):
        result = _run_cli("nonexistent")
        assert result.returncode == 2

    def test_subcommand_help_fetch(self):
        result = _run_cli("fetch", "--help")
        assert result.returncode == 0
        assert "UniProt" in result.stdout or "accession" in result.stdout.lower()

    def test_subcommand_help_collect(self):
        result = _run_cli("collect", "--help")
        assert result.returncode == 0
        assert (
            "directories" in result.stdout.lower()
            or "output" in result.stdout.lower()
        )

    def test_subcommand_help_compare(self):
        result = _run_cli("compare", "--help")
        assert result.returncode == 0

    def test_subcommand_help_lightdock(self):
        result = _run_cli("lightdock", "--help")
        assert result.returncode == 0
        assert "receptor" in result.stdout.lower()

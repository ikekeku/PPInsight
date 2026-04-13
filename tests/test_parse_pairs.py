"""Tests for the parse_pairs module."""

import textwrap

import pandas as pd
import pytest

from ppinsight.parse_pairs import (
    _extract_refs,
    _strip_refs,
    main,
    parse_interaction_table,
)

# ── Helpers ──────────────────────────────────────────────────────

def _write_tsv(path, text):
    """Write *text* (dedented) to *path* as a TSV file."""
    with open(path, "w") as f:
        f.write(textwrap.dedent(text))


def _rtk_table(tmp_path):
    """Create a minimal RTK-interactome-style table and return its path."""
    tsv = tmp_path / "table.tsv"
    # Columns: (blank) | Family | Member | Int1-5 | NonInt1-4
    lines = [
        "\t\tProtein A\t\t\t\t\tProtein B\t\t\t\t",
        "\tFamily\tMember\tInteractions\t\t\t\t\tNon-Interactions\t\t\t",
        "\tEphrin\tEPHA1\tERBB2 [6]\tFGFR2 [7]\t\t\t\tFGFR3\tMET\t\t",
        "\t\tEPHA2\tMET [10][12]\t\t\t\t\tFGFR1\t\t\t",
        "\tFGFR\tFGFR1\tERBB2 [3]\t\t\t\t\t\t\t\t",
    ]
    tsv.write_text("\n".join(lines) + "\n")
    return tsv


# ── Unit tests: helper functions ─────────────────────────────────

class TestStripRefs:
    def test_single_ref(self):
        assert _strip_refs("ERBB2 [6]") == "ERBB2"

    def test_multiple_refs(self):
        assert _strip_refs("MET [10][12]") == "MET"

    def test_no_refs(self):
        assert _strip_refs("FGFR3") == "FGFR3"

    def test_empty(self):
        assert _strip_refs("") == ""

    def test_whitespace(self):
        assert _strip_refs("  ERBB2  [6]  ") == "ERBB2"


class TestExtractRefs:
    def test_single(self):
        assert _extract_refs("ERBB2 [6]") == "6"

    def test_multiple(self):
        assert _extract_refs("MET [10][12]") == "10,12"

    def test_none(self):
        assert _extract_refs("FGFR3") == ""

    def test_empty(self):
        assert _extract_refs("") == ""


# ── Integration tests: parse_interaction_table ───────────────────

class TestParseInteractionTable:
    def test_basic_parse(self, tmp_path):
        tsv = _rtk_table(tmp_path)
        df = parse_interaction_table(tsv)

        assert isinstance(df, pd.DataFrame)
        assert set(df.columns) == {
            "proteinA", "proteinB", "label", "family", "references",
        }
        assert len(df) > 0

    def test_interactions_extracted(self, tmp_path):
        tsv = _rtk_table(tmp_path)
        df = parse_interaction_table(tsv)

        interactions = df[df["label"] == "interaction"]
        pairs = set(zip(
            interactions["proteinA"],
            interactions["proteinB"],
            strict=False,
        ))
        assert ("EPHA1", "ERBB2") in pairs
        assert ("EPHA1", "FGFR2") in pairs
        assert ("EPHA2", "MET") in pairs
        assert ("FGFR1", "ERBB2") in pairs

    def test_non_interactions_extracted(self, tmp_path):
        tsv = _rtk_table(tmp_path)
        df = parse_interaction_table(tsv)

        non_int = df[df["label"] == "non-interaction"]
        pairs = set(zip(
            non_int["proteinA"],
            non_int["proteinB"],
            strict=False,
        ))
        assert ("EPHA1", "FGFR3") in pairs
        assert ("EPHA1", "MET") in pairs
        assert ("EPHA2", "FGFR1") in pairs

    def test_family_propagation(self, tmp_path):
        """Continuation rows inherit family from the previous non-empty row."""
        tsv = _rtk_table(tmp_path)
        df = parse_interaction_table(tsv)

        epha2 = df[df["proteinA"] == "EPHA2"]
        assert all(epha2["family"] == "Ephrin")

        fgfr1 = df[df["proteinA"] == "FGFR1"]
        assert all(fgfr1["family"] == "FGFR")

    def test_references_extracted(self, tmp_path):
        tsv = _rtk_table(tmp_path)
        df = parse_interaction_table(tsv)

        row = df[(df["proteinA"] == "EPHA1") & (df["proteinB"] == "ERBB2")]
        assert row.iloc[0]["references"] == "6"

        row = df[(df["proteinA"] == "EPHA2") & (df["proteinB"] == "MET")]
        assert row.iloc[0]["references"] == "10,12"

    def test_deduplication(self, tmp_path):
        """Duplicate (proteinA, proteinB, label) rows are dropped."""
        tsv = _rtk_table(tmp_path)
        df = parse_interaction_table(tsv)
        dups = df.duplicated(subset=["proteinA", "proteinB", "label"])
        assert not dups.any()

    def test_empty_table_raises(self, tmp_path):
        tsv = tmp_path / "empty.tsv"
        lines = [
            "\t\tProtein A\t\t\tProtein B\t\t",
            "\tFamily\tMember\tInteractions\t\tNon-Interactions\t",
        ]
        tsv.write_text("\n".join(lines) + "\n")
        with pytest.raises(ValueError, match="No protein pairs"):
            parse_interaction_table(tsv)

    def test_no_header_raises(self, tmp_path):
        """Table without recognizable headers should raise."""
        tsv = tmp_path / "bad.tsv"
        tsv.write_text("col1\tcol2\tcol3\nA\tB\tC\n")
        with pytest.raises(ValueError):
            parse_interaction_table(tsv)


# ── CLI tests ────────────────────────────────────────────────────

class TestCLI:
    def test_csv_output(self, tmp_path):
        tsv = _rtk_table(tmp_path)
        out = tmp_path / "pairs.csv"
        main([str(tsv), "-o", str(out)])
        assert out.exists()
        df = pd.read_csv(out)
        assert "proteinA" in df.columns
        assert len(df) > 0

    def test_tsv_output(self, tmp_path):
        tsv = _rtk_table(tmp_path)
        out = tmp_path / "pairs.tsv"
        main([str(tsv), "-o", str(out)])
        assert out.exists()
        df = pd.read_csv(out, sep="\t")
        assert "proteinA" in df.columns

    def test_stats_flag(self, tmp_path, capsys):
        tsv = _rtk_table(tmp_path)
        out = tmp_path / "pairs.csv"
        main([str(tsv), "-o", str(out), "--stats"])
        captured = capsys.readouterr()
        assert "Interactions:" in captured.out
        assert "Non-interactions:" in captured.out
        assert "Unique proteins:" in captured.out

    def test_default_output_name(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        tsv = _rtk_table(tmp_path)
        main([str(tsv)])
        assert (tmp_path / "pairs.csv").exists()

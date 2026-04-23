"""Path-level tests for pdb_to_rosetta output directory behavior."""

from ppinsight import pdb_to_rosetta


def test_make_output_dir_auto_increments_when_existing(tmp_path):
    rec = tmp_path / "rec.pdb"
    lig = tmp_path / "lig.pdb"
    rec.write_text("ATOM\n", encoding="utf-8")
    lig.write_text("ATOM\n", encoding="utf-8")

    run1 = pdb_to_rosetta._make_output_dir(
        str(rec),
        str(lig),
        base_root=str(tmp_path / "out"),
        method="rosetta_runs",
    )
    run2 = pdb_to_rosetta._make_output_dir(
        str(rec),
        str(lig),
        base_root=str(tmp_path / "out"),
        method="rosetta_runs",
    )

    assert run1 != run2
    assert run2.endswith("_1")

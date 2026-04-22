from ppinsight.utils import resolve_input_path


def test_resolve_input_path_accepts_accession_stem(tmp_path):
    pdb_file = tmp_path / "P69905.pdb"
    pdb_file.write_text("ATOM\n", encoding="utf-8")

    resolved = resolve_input_path("P69905", search_root=str(tmp_path))

    assert resolved == str(pdb_file.resolve())


def test_resolve_input_path_accepts_legacy_ent_files(tmp_path):
    ent_file = tmp_path / "pdb1a00.ent"
    ent_file.write_text("ATOM\n", encoding="utf-8")

    resolved = resolve_input_path("pdb1a00", search_root=str(tmp_path))

    assert resolved == str(ent_file.resolve())

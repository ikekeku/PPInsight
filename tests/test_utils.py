from ppinsight.utils import (
    copy_pdb_selected_chains,
    dbref_chains_for_accession,
    find_column,
    normalize_column_name,
    resolve_input_path,
    resolve_traceability_path,
)


def test_normalize_column_name_lowercases_and_replaces_spaces():
    assert normalize_column_name(" Protein A ") == "protein_a"


def test_find_column_matches_case_insensitively():
    columns = ["Model", "proteinA", "Output Path"]

    assert find_column(columns, "model") == "Model"
    assert find_column(columns, "proteina") == "proteinA"
    assert find_column(columns, "output_path") == "Output Path"
    assert find_column(columns, "missing") is None


def test_resolve_traceability_path_handles_relative_paths(tmp_path):
    resolved = resolve_traceability_path("models/pose_1.pdb", base_dir=tmp_path)

    assert resolved == str((tmp_path / "models" / "pose_1.pdb").resolve())


def test_resolve_traceability_path_treats_placeholder_values_as_empty(tmp_path):
    assert resolve_traceability_path("nan", base_dir=tmp_path) == ""
    assert resolve_traceability_path(" - ", base_dir=tmp_path) == ""


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


def test_dbref_chains_for_accession_returns_matching_chains(tmp_path):
    pdb_file = tmp_path / "P12345.pdb"
    pdb_file.write_text(
        "DBREF  1ABC A    1     1  UNP    P12345   TEST_HUMAN      1      1\n"
        "DBREF2 1ABC B P12345\n"
        "DBREF  1ABC X    1     1  UNP    Q99999   OTHER_HUMAN     1      1\n",
        encoding="utf-8",
    )

    assert dbref_chains_for_accession(pdb_file, "P12345") == ["A", "B"]


def test_copy_pdb_selected_chains_drops_unselected_coordinates_and_ter(tmp_path):
    src = tmp_path / "source.pdb"
    dst = tmp_path / "filtered.pdb"
    src_contents = (
        "HEADER    TEST\n"
        "ATOM      1  N   GLY A   1      11.000  11.000  11.000"
        "  1.00 20.00           N\n"
        "TER       2      GLY A   1\n"
        "ATOM      3  N   SER X   2      12.000  12.000  12.000"
        "  1.00 20.00           N\n"
        "TER       4      SER X   2\n"
        "ATOM      5  N   TYR B   3      13.000  13.000  13.000"
        "  1.00 20.00           N\n"
        "TER       6      TYR B   3\n"
        "END\n"
    )
    src.write_text(src_contents, encoding="utf-8")

    copy_pdb_selected_chains(src, dst, {"A", "B"})

    assert dst.read_text(encoding="utf-8").splitlines() == [
        "HEADER    TEST",
        (
            "ATOM      1  N   GLY A   1      11.000  11.000  11.000"
            "  1.00 20.00           N"
        ),
        "TER       2      GLY A   1",
        (
            "ATOM      5  N   TYR B   3      13.000  13.000  13.000"
            "  1.00 20.00           N"
        ),
        "TER       6      TYR B   3",
        "END",
    ]

from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

from src.iatlascbioportalexport import lens_files


@pytest.fixture
def mock_syn(tmp_path):
    """
    Create a mock Synapse client whose .get() writes pre-baked files
    into tmp_path and returns objects with .path.
    """
    syn = MagicMock()

    # Prepare "downloaded" files
    downloaded_data = tmp_path / "orig_data.txt"
    downloaded_data.write_text("DATA CONTENT\n", encoding="utf-8")

    downloaded_meta = tmp_path / "orig_meta.txt"
    downloaded_meta.write_text(
        "\n".join(
            [
                "cancer_study_identifier: OLD",
                "genetic_alteration_type: MRNA_EXPRESSION",
                "datatype: CONTINUOUS",
                "stable_id: rna_seq_mrna",
                "profile_name: mRNA expression",
                "profile_description: Expression levels",
                "data_filename: old_name.txt",
                "show_profile_in_analysis_tab: false",
                "",
            ]
        ),
        encoding="utf-8",
    )

    def get_side_effect(synid, downloadLocation=None):
        # Simulate syn.get returning an object with .path
        if synid == "synDATA":
            return SimpleNamespace(path=str(downloaded_data))
        if synid == "synMETA":
            return SimpleNamespace(path=str(downloaded_meta))
        raise ValueError(f"Unexpected synid: {synid}")

    syn.get.side_effect = get_side_effect
    return syn


def test_parse_metadata_kv_basic():
    text = (
        "cancer_study_identifier: jkdjakf\n"
        "genetic_alteration_type: MRNA_EXPRESSION\n"
        "datatype: CONTINUOUS\n"
    )
    assert lens_files.parse_metadata_kv(text) == {
        "cancer_study_identifier": "jkdjakf",
        "genetic_alteration_type": "MRNA_EXPRESSION",
        "datatype": "CONTINUOUS",
    }


def test_parse_metadata_kv_ignores_comments_and_blank_lines():
    text = (
        "# a comment line\n"
        "\n"
        "key1: val1\n"
        "   # indented comment\n"
        "key2: val2\n"
        "\n"
    )
    assert lens_files.parse_metadata_kv(text) == {"key1": "val1", "key2": "val2"}


def test_parse_metadata_kv_strips_whitespace_around_key_and_value():
    text = "  key1   :   val1  \n" "\tkey2:\tval2\t\n"
    assert lens_files.parse_metadata_kv(text) == {"key1": "val1", "key2": "val2"}


def test_parse_metadata_kv_skips_lines_without_colon():
    text = "key1: val1\n" "this line is malformed\n" "key2: val2\n" "also_malformed\n"
    assert lens_files.parse_metadata_kv(text) == {"key1": "val1", "key2": "val2"}


def test_parse_metadata_kv_allows_colon_in_value():
    text = "url: https://example.org/a:b:c\n"
    assert lens_files.parse_metadata_kv(text) == {"url": "https://example.org/a:b:c"}


def test_parse_metadata_kv_last_value_wins_on_duplicate_keys():
    text = "key1: val1\n" "key1: val2\n"
    assert lens_files.parse_metadata_kv(text) == {"key1": "val2"}


def test_write_metadata_kv_basic():
    meta = {"a": "1", "b": "2"}
    out = lens_files.write_metadata_kv(meta)

    # order may vary; test as a set of lines
    lines = [line for line in out.strip().splitlines()]
    assert set(lines) == {"a: 1", "b: 2"}
    assert out.endswith("\n")


def test_write_metadata_kv_empty_dict():
    assert lens_files.write_metadata_kv({}) == "\n"


def test_roundtrip_parse_then_write_then_parse_is_idempotent():
    text = (
        "# comment\n"
        "key1: val1\n"
        "key2: val2 with spaces\n"
        "url: https://example.org/a:b\n"
        "\n"
    )
    parsed = lens_files.parse_metadata_kv(text)
    written = lens_files.write_metadata_kv(parsed)
    reparsed = lens_files.parse_metadata_kv(written)
    assert reparsed == parsed


@pytest.mark.parametrize(
    "text, expected",
    [
        ("k:v\n", {"k": "v"}),
        ("k: v\nk2:v2\n", {"k": "v", "k2": "v2"}),
        ("   \n#x\nk: v \n", {"k": "v"}),
        ("no_colon\nk: v\n", {"k": "v"}),
    ],
)
def test_parse_metadata_kv_parametrized(text, expected):
    assert lens_files.parse_metadata_kv(text) == expected


def test_download_and_patch_files_happy_path(mock_syn, tmp_path, capsys):
    out_dir = tmp_path

    lens_files.download_and_patch_files(
        syn=mock_syn,
        data_synid="synDATA",
        meta_synid="synMETA",
        out_dir=str(out_dir),
        out_data_filename="data_rna_seq_mrna.txt",
        out_meta_filename="meta_rna_seq_mrna.txt",
        cancer_study_identifier="something",
        data_filename_in_meta="data_rna_seq_mrna.txt",
    )

    # syn.get called correctly
    mock_syn.get.assert_any_call("synDATA", downloadLocation=str(out_dir))
    mock_syn.get.assert_any_call("synMETA", downloadLocation=str(out_dir))
    assert mock_syn.get.call_count == 2

    # data file copied to standardized name
    final_data = out_dir / "data_rna_seq_mrna.txt"
    assert final_data.exists()
    assert final_data.read_text(encoding="utf-8") == "DATA CONTENT\n"

    # metadata file written to standardized name and patched
    final_meta = out_dir / "meta_rna_seq_mrna.txt"
    assert final_meta.exists()
    meta_text = final_meta.read_text(encoding="utf-8")

    # updated keys present
    assert "cancer_study_identifier: something" in meta_text
    assert "data_filename: data_rna_seq_mrna.txt" in meta_text

    # an unrelated key preserved
    assert "datatype: CONTINUOUS" in meta_text

    # prints something helpful
    captured = capsys.readouterr()
    assert "Wrote:" in captured.out
    assert "data_rna_seq_mrna.txt" in captured.out
    assert "meta_rna_seq_mrna.txt" in captured.out


def test_download_and_patch_files_overwrites_prior_outputs(mock_syn, tmp_path):
    out_dir = tmp_path

    # Pre-create outputs with different content to ensure overwrite
    (out_dir / "data_rna_seq_mrna.txt").write_text("OLD DATA", encoding="utf-8")
    (out_dir / "meta_rna_seq_mrna.txt").write_text("OLD META", encoding="utf-8")

    lens_files.download_and_patch_files(
        syn=mock_syn,
        data_synid="synDATA",
        meta_synid="synMETA",
        out_dir=str(out_dir),
        out_data_filename="data_rna_seq_mrna.txt",
        out_meta_filename="meta_rna_seq_mrna.txt",
        cancer_study_identifier="something",
        data_filename_in_meta="data_rna_seq_mrna.txt",
    )

    assert (out_dir / "data_rna_seq_mrna.txt").read_text(
        encoding="utf-8"
    ) == "DATA CONTENT\n"
    assert "cancer_study_identifier: something" in (
        out_dir / "meta_rna_seq_mrna.txt"
    ).read_text(encoding="utf-8")

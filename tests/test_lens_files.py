import pytest

from src.iatlascbioportalexport import lens_files


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
    text = (
        "  key1   :   val1  \n"
        "\tkey2:\tval2\t\n"
    )
    assert lens_files.parse_metadata_kv(text) == {"key1": "val1", "key2": "val2"}


def test_parse_metadata_kv_skips_lines_without_colon():
    text = (
        "key1: val1\n"
        "this line is malformed\n"
        "key2: val2\n"
        "also_malformed\n"
    )
    assert lens_files.parse_metadata_kv(text) == {"key1": "val1", "key2": "val2"}


def test_parse_metadata_kv_allows_colon_in_value():
    text = "url: https://example.org/a:b:c\n"
    assert lens_files.parse_metadata_kv(text) == {"url": "https://example.org/a:b:c"}


def test_parse_metadata_kv_last_value_wins_on_duplicate_keys():
    text = (
        "key1: val1\n"
        "key1: val2\n"
    )
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

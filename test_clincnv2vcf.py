import tempfile

import pandas as pd
import pytest

from clincnv2vcf import clean_info, get_DELDUP_by_CN, parse_metadata


def test_parse_metadata() -> None:
    # Test case 1: file with no metadata
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as f:
        f.write("This is a test file with no metadata.\n")
        file_path = f.name
    assert parse_metadata(file_path) == ([], {})

    # Test case 2: file with metadata but no key-value pairs
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as f:
        f.write("## This is a test file with metadata but no key-value pairs.\n")
        f.write("##\n")
        f.write("## More metadata with no key-value pairs.\n")
        file_path = f.name
    assert parse_metadata(file_path) == ([0, 1, 2], {})

    # Test case 3: file with metadata and key-value pairs
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as f:
        f.write("## This is a test file with metadata and key-value pairs.\n")
        f.write("## Key1: Value1\n")
        f.write("## Key2: Value2\n")
        file_path = f.name
    assert parse_metadata(file_path) == (
        [0, 1, 2],
        {"Key1": "Value1", "Key2": "Value2"},
    )


def test_get_DELDUP_by_CN_normal() -> None:
    assert get_DELDUP_by_CN(0) == "DEL"
    assert get_DELDUP_by_CN(0, ltgt=True) == "<DEL>"
    assert get_DELDUP_by_CN(4) == "DUP"
    assert get_DELDUP_by_CN(4, ltgt=True) == "<DUP>"
    assert get_DELDUP_by_CN(100) == "DUP"
    assert get_DELDUP_by_CN(100, ltgt=True) == "<DUP>"
    assert get_DELDUP_by_CN(3, male_x=True) == "DUP"
    assert get_DELDUP_by_CN(3, ltgt=True, male_x=True) == "<DUP>"
    assert get_DELDUP_by_CN(100, male_x=True) == "DUP"
    assert get_DELDUP_by_CN(100, ltgt=True, male_x=True) == "<DUP>"


def test_get_DELDUP_by_CN_edge_cases() -> None:
    assert get_DELDUP_by_CN(0, male_x=True) == "DEL"
    assert get_DELDUP_by_CN(0, ltgt=True, male_x=True) == "<DEL>"
    assert get_DELDUP_by_CN(1, male_x=True) == "."
    assert get_DELDUP_by_CN(1, ltgt=True, male_x=True) == "."
    assert get_DELDUP_by_CN(2, male_x=True) == "DUP"
    assert get_DELDUP_by_CN(2, ltgt=True, male_x=True) == "<DUP>"
    assert get_DELDUP_by_CN(1) == "DEL"
    assert get_DELDUP_by_CN(1, ltgt=True) == "<DEL>"
    assert get_DELDUP_by_CN(2) == "."
    assert get_DELDUP_by_CN(2, ltgt=True) == "."
    assert get_DELDUP_by_CN(3) == "DUP"
    assert get_DELDUP_by_CN(3, ltgt=True) == "<DUP>"


def test_get_DELDUP_by_CN_error() -> None:
    with pytest.raises(TypeError):
        get_DELDUP_by_CN(None)  # type: ignore


def test_clean_info_with_none() -> None:
    assert clean_info(None) == "."  # type: ignore


def test_clean_info_with_nan() -> None:
    assert clean_info(pd.NA) == "."  # type: ignore


def test_clean_info_with_alphanumeric_chars() -> None:
    assert clean_info("abc123") == "abc123"


def test_clean_info_with_special_chars() -> None:
    assert clean_info("a#b$c%d^") == "a_b_c_d_"


def test_clean_info_with_mixed_chars() -> None:
    assert clean_info("a#b123$c45%d^") == "a_b123_c45_d_"

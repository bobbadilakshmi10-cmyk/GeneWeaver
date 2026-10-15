
from pathlib import Path

import pytest

from member1.fasta_pipeline import (
    clean_dna,
    encode_dna,
    prepare_data,
)


DATA_FILE = Path("data/sample.fasta")


def test_lowercase_conversion():
    """Lowercase DNA should become uppercase."""
    assert clean_dna("atgc") == "ATGC"


def test_invalid_base():
    """Unsupported DNA characters should raise an error."""
    with pytest.raises(ValueError):
        clean_dna("ATGX")


def test_dna_encoding():
    """Check numerical values for DNA bases."""
    assert encode_dna("ACGTN").tolist() == [
        0, 1, 2, 3, -1
    ]


def test_overlapping_chunks():
    """Check the correct overlapping chunks."""
    records = prepare_data(
        DATA_FILE,
        chunk_size=5,
        step=1,
    )

    assert len(records) > 0

    # Check the first two chunks of sequence_1.
    assert records[0]["sequence"] == "ATGCG"
    assert records[1]["sequence"] == "TGCGT"

    # Consecutive chunks should overlap by four bases.
    assert (
        records[0]["sequence"][1:]
        == records[1]["sequence"][:-1]
    )


def test_chunk_positions():
    """Check that chunk positions use 1-based indexing."""
    records = prepare_data(
        DATA_FILE,
        chunk_size=5,
        step=1,
    )

    assert records[0]["position"] == 1
    assert records[1]["position"] == 2


def test_invalid_chunk_size():
    """A zero or negative chunk size must be rejected."""
    with pytest.raises(ValueError):
        prepare_data(DATA_FILE, chunk_size=0)


def test_missing_file():
    """A missing FASTA file should raise an error."""
    with pytest.raises(FileNotFoundError):
        prepare_data("data/missing.fasta")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])


# Import Path to handle file paths
from pathlib import Path

# Import pytest to write and run test cases
import pytest

# Import the functions to be tested from the DNA pipeline
from member1.fasta_pipeline import (
    clean_dna,
    encode_dna,
    prepare_data,
)

# Define the path of the sample FASTA file
DATA_FILE = Path("data/sample.fasta")


def test_lowercase_conversion():
    """Test whether lowercase DNA becomes uppercase."""

    # Convert lowercase DNA into uppercase and check the result
    assert clean_dna("atgc") == "ATGC"


def test_invalid_base():
    """Test whether invalid DNA characters raise an error."""

    # Check that an unsupported character (X) raises ValueError
    with pytest.raises(ValueError):
        clean_dna("ATGX")


def test_dna_encoding():
    """Test whether DNA bases are converted into numbers."""

    # Encode A=0, C=1, G=2, T=3, and N=-1
    # Convert the NumPy array into a list for comparison
    assert encode_dna("ACGTN").tolist() == [
        0, 1, 2, 3, -1
    ]


def test_overlapping_chunks():
    """Test whether consecutive DNA chunks overlap correctly."""

    # Prepare DNA data using chunks of 5 bases
    # step=1 means each new chunk starts one base later
    records = prepare_data(
        DATA_FILE,
        chunk_size=5,
        step=1,
    )

    # Make sure at least one DNA chunk is available
    assert len(records) > 0

    # Check the first chunk of sequence_1
    assert records[0]["sequence"] == "ATGCG"

    # Check the second chunk, which starts one base later
    assert records[1]["sequence"] == "TGCGT"

    # Verify that the last four bases of the first chunk
    # match the first four bases of the second chunk
    assert (
        records[0]["sequence"][1:]
        == records[1]["sequence"][:-1]
    )


def test_chunk_positions():
    """Test whether DNA chunk positions start from 1."""

    # Prepare DNA data with overlapping chunks
    records = prepare_data(
        DATA_FILE,
        chunk_size=5,
        step=1,
    )

    # The first chunk should start at position 1
    assert records[0]["position"] == 1

    # The second chunk should start at position 2
    assert records[1]["position"] == 2


def test_invalid_chunk_size():
    """Test whether zero or negative chunk sizes are rejected."""

    # A chunk size of zero should raise ValueError
    with pytest.raises(ValueError):
        prepare_data(DATA_FILE, chunk_size=0)


def test_missing_file():
    """Test whether a missing FASTA file raises an error."""

    # Try to read a file that does not exist
    # FileNotFoundError should be raised
    with pytest.raises(FileNotFoundError):
        prepare_data("data/missing.fasta")


# Run the tests when this file is executed directly
if __name__ == "__main__":

    # Run pytest on this file and display detailed test results
    pytest.main([__file__, "-v"])

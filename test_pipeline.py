
# Import Path to work with file and folder paths
from pathlib import Path

# Import pytest to create and execute test cases
import pytest

# Import DNA processing functions from the Member 1 FASTA pipeline
from member1.fasta_pipeline import (
    clean_dna,
    encode_dna,
    prepare_data,
)

# Store the location of the sample FASTA file used for testing
DATA_FILE = Path("data/sample.fasta")


# Test Case 1: Verify that lowercase DNA bases are converted to uppercase
def test_lowercase_conversion():
    """Check whether lowercase DNA sequences are converted to uppercase."""

    # Pass a lowercase DNA sequence to the cleaning function
    # Expected result should contain only uppercase DNA bases
    assert clean_dna("atgc") == "ATGC"


# Test Case 2: Verify that invalid DNA characters are rejected
def test_invalid_base():
    """Check whether invalid characters raise a ValueError."""

    # X is not a valid DNA base
    # The function should raise ValueError when an invalid base is found
    with pytest.raises(ValueError):
        clean_dna("ATGX")


# Test Case 3: Verify that DNA bases are encoded into numerical values
def test_dna_encoding():
    """Check whether DNA bases are converted into their correct numeric codes."""

    # Expected DNA encoding:
    # A = 0, C = 1, G = 2, T = 3, N = -1

    # Encode the DNA sequence and convert the NumPy array into a Python list
    # Compare the actual encoded values with the expected values
    assert encode_dna("ACGTN").tolist() == [
        0, 1, 2, 3, -1
    ]


# Test Case 4: Verify that DNA sequences are divided into overlapping chunks
def test_overlapping_chunks():
    """Check whether consecutive DNA chunks overlap correctly."""

    # Read the sample FASTA file and divide sequences into chunks
    # chunk_size=5 means each chunk contains 5 DNA bases
    # step=1 means the starting position moves forward by one base
    records = prepare_data(
        DATA_FILE,
        chunk_size=5,
        step=1,
    )

    # Verify that the FASTA file produces at least one DNA chunk
    assert len(records) > 0

    # Check whether the first DNA chunk is generated correctly
    assert records[0]["sequence"] == "ATGCG"

    # Check whether the second DNA chunk starts one base later
    assert records[1]["sequence"] == "TGCGT"

    # Verify the overlapping region between the first and second chunks
    # The last four bases of the first chunk must match
    # the first four bases of the second chunk
    assert (
        records[0]["sequence"][1:]
        == records[1]["sequence"][:-1]
    )


# Test Case 5: Verify that DNA chunk positions start from 1
def test_chunk_positions():
    """Check whether DNA chunks have the correct starting positions."""

    # Prepare DNA chunks from the sample FASTA file
    # Each chunk contains 5 bases with an overlap of 4 bases
    records = prepare_data(
        DATA_FILE,
        chunk_size=5,
        step=1,
    )

    # Verify that the first DNA chunk starts at position 1
    assert records[0]["position"] == 1

    # Verify that the second DNA chunk starts at position 2
    assert records[1]["position"] == 2


# Test Case 6: Verify that invalid chunk sizes are rejected
def test_invalid_chunk_size():
    """Check whether zero or negative chunk sizes raise ValueError."""

    # A chunk size must be a positive integer
    # Passing zero should raise a ValueError
    with pytest.raises(ValueError):
        prepare_data(DATA_FILE, chunk_size=0)


# Test Case 7: Verify that missing FASTA files are handled correctly
def test_missing_file():
    """Check whether a missing FASTA file raises FileNotFoundError."""

    # Provide the path of a file that does not exist
    # The function should raise FileNotFoundError
    with pytest.raises(FileNotFoundError):
        prepare_data("data/missing.fasta")


# Main Entry Point: Execute test cases when this file runs directly
if __name__ == "__main__":

    # Run all test functions in this Python file
    # -v enables verbose output to display individual test results
    pytest.main([__file__, "-v"])

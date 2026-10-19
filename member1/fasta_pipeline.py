
# Import Path to check and handle file paths.
from pathlib import Path

# Import StringIO to treat a string like a file.
from io import StringIO

# Import NumPy to store DNA bases as numerical arrays.
import numpy as np

# Import SeqIO from BioPython to read FASTA files.
from Bio import SeqIO


# Map each DNA base to a numerical value.
DNA_ENCODING = {
    "A": 0,    # Adenine
    "C": 1,    # Cytosine
    "G": 2,    # Guanine
    "T": 3,    # Thymine
    "N": -1,   # Unknown or undetermined DNA base
}


# Function to clean and validate a DNA sequence.
def clean_dna(sequence):

    # Remove whitespace and convert all letters to uppercase.
    sequence = "".join(sequence.split()).upper()

    # Find characters that are not valid DNA bases.
    invalid = set(sequence) - set(DNA_ENCODING)

    # Raise an error if invalid characters are found.
    if invalid:
        raise ValueError(f"Invalid DNA bases: {sorted(invalid)}")

    # Return the cleaned DNA sequence.
    return sequence


# Function to convert DNA letters into numerical values.
def encode_dna(sequence):

    # Clean and validate the DNA sequence first.
    sequence = clean_dna(sequence)

    # Convert each DNA base into its corresponding number.
    # Store the numbers in a NumPy array using 8-bit integers.
    return np.array(
        [DNA_ENCODING[base] for base in sequence],
        dtype=np.int8,
    )


# Function to read a FASTA file and prepare DNA chunks.
# file_path: path of the FASTA file.
# chunk_size: maximum number of bases in each chunk.
# step: number of bases to move before creating the next chunk.
def prepare_data(file_path, chunk_size=5, step=1):

    # Check that chunk_size and step are positive.
    if chunk_size <= 0 or step <= 0:
        raise ValueError("chunk_size and step must be positive")

    # Convert the file path into a Path object.
    path = Path(file_path)

    # Check whether the FASTA file exists.
    if not path.is_file():
        raise FileNotFoundError(f"FASTA file not found: {path}")

    # Open the FASTA file for reading.
    # utf-8-sig also handles a possible UTF-8 byte-order mark.
    with path.open("r", encoding="utf-8-sig") as file:

        # Remove leading/trailing whitespace and skip blank lines.
        lines = [
            line.strip()
            for line in file
            if line.strip()
        ]

    # Find the index of the first FASTA header.
    # A FASTA header starts with the character ">".
    # Return None if no header is found.
    first_header = next(
        (i for i, line in enumerate(lines)
         if line.startswith(">")),
        None,
    )

    # Stop if the file does not contain a FASTA header.
    if first_header is None:
        raise ValueError(
            "No FASTA header found. Each sequence must start with >."
        )

    # Ignore any accidental text before the first FASTA header.
    lines = lines[first_header:]

    # Combine the cleaned lines into one text string.
    fasta_text = "\n".join(lines) + "\n"

    # Parse the FASTA text using BioPython.
    # StringIO allows SeqIO to read the text like a file.
    # list() collects all FASTA records into a list.
    fasta_records = list(
        SeqIO.parse(StringIO(fasta_text), "fasta")
    )

    # Stop if no FASTA records were parsed.
    if not fasta_records:
        raise ValueError("No valid FASTA sequences found.")

    # Create an empty list to store all DNA chunks.
    records = []

    # Process each FASTA record one by one.
    for fasta_record in fasta_records:

        # Get the DNA sequence, convert it to a string,
        # and clean and validate the sequence.
        dna = clean_dna(str(fasta_record.seq))

        # Raise an error if the DNA sequence is empty.
        if not dna:
            raise ValueError(f"Empty DNA sequence: {fasta_record.id}")

        # Convert the complete DNA sequence into numbers.
        encoded = encode_dna(dna)

        # Calculate the starting index of each chunk.
        # If DNA is shorter than chunk_size, create one short chunk.
        # Otherwise, create chunks using the specified step.
        starts = (
            [0]
            if len(dna) < chunk_size
            else range(0, len(dna) - chunk_size + 1, step)
        )

        # Process each starting index to create a DNA chunk.
        for start in starts:

            # Calculate the ending index without exceeding DNA length.
            end = min(start + chunk_size, len(dna))

            # Store the chunk details in a dictionary.
            records.append({
                # Store the FASTA sequence identifier.
                "id": fasta_record.id,

                # Convert the zero-based index to a one-based position.
                "position": start + 1,

                # Extract the DNA letters for this chunk.
                "sequence": dna[start:end],

                # Extract the corresponding encoded numbers.
                # copy() creates a separate NumPy array for this chunk.
                "encoded": encoded[start:end].copy(),

                # Calculate the number of bases in this chunk.
                "length": end - start,
            })

    # Return the list containing all prepared DNA chunks.
    return records


from pathlib import Path
from io import StringIO

import numpy as np
from Bio import SeqIO

DNA_ENCODING = {
    "A": 0,
    "C": 1,
    "G": 2,
    "T": 3,
    "N": -1,
}


def clean_dna(sequence):
    sequence = "".join(sequence.split()).upper()

    invalid = set(sequence) - set(DNA_ENCODING)
    if invalid:
        raise ValueError(f"Invalid DNA bases: {sorted(invalid)}")

    return sequence


def encode_dna(sequence):
    sequence = clean_dna(sequence)
    return np.array(
        [DNA_ENCODING[base] for base in sequence],
        dtype=np.int8,
    )


def prepare_data(file_path, chunk_size=5, step=1):
    if chunk_size <= 0 or step <= 0:
        raise ValueError("chunk_size and step must be positive")

    path = Path(file_path)

    if not path.is_file():
        raise FileNotFoundError(f"FASTA file not found: {path}")

    # Read the file and remove blank lines.
    with path.open("r", encoding="utf-8-sig") as file:
        lines = [
            line.strip()
            for line in file
            if line.strip()
        ]

    # Remove comments and accidental text before the first header.
    first_header = next(
        (i for i, line in enumerate(lines)
         if line.startswith(">")),
        None,
    )

    if first_header is None:
        raise ValueError(
            "No FASTA header found. Each sequence must start with >."
        )

    lines = lines[first_header:]

    # Parse the cleaned text using BioPython.
    fasta_text = "\n".join(lines) + "\n"
    fasta_records = list(
        SeqIO.parse(StringIO(fasta_text), "fasta")
    )

    if not fasta_records:
        raise ValueError("No valid FASTA sequences found.")

    records = []

    for fasta_record in fasta_records:
        dna = clean_dna(str(fasta_record.seq))

        if not dna:
            raise ValueError(f"Empty DNA sequence: {fasta_record.id}")

        encoded = encode_dna(dna)

        # Create overlapping chunks.
        starts = (
            [0]
            if len(dna) < chunk_size
            else range(0, len(dna) - chunk_size + 1, step)
        )

        for start in starts:
            end = min(start + chunk_size, len(dna))

            records.append({
                "id": fasta_record.id,
                "position": start + 1,
                "sequence": dna[start:end],
                "encoded": encoded[start:end].copy(),
                "length": end - start,
            })

    return records

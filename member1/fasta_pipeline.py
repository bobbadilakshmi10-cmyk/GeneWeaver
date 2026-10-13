
# Define a function to read a FASTA file and split DNA sequences into chunks.
# file_path: path of the input FASTA file.
# chunk_size=5: each DNA chunk contains up to 5 characters by default.
def prepare_data(file_path, chunk_size=5):

    # Create an empty list to store all the DNA chunks.
    records = []

    # Store the ID of the current DNA sequence.
    # Initially, no DNA sequence has been started.
    current_id = None

    # Create an empty list to collect lines of the current DNA sequence.
    sequence_parts = []

    # Open the FASTA file in read mode.
    # The 'with' statement automatically closes the file after reading.
    with open(file_path, "r") as file:

        # Read the file one line at a time.
        for line in file:

            # Remove spaces and newline characters from both ends.
            line = line.strip()

            # If the line is empty, skip it and continue to the next line.
            if not line:
                continue

            # FASTA header lines begin with the '>' symbol.
            if line.startswith(">"):

                # Check whether a previous DNA sequence has been collected.
                if current_id is not None:

                    # Join all collected DNA parts into one string.
                    # Convert all DNA letters to uppercase.
                    dna = "".join(sequence_parts).upper()

                    # Move through the DNA sequence in steps of chunk_size.
                    for i in range(0, len(dna), chunk_size):

                        # Extract a DNA chunk starting at position i.
                        # The chunk contains at most chunk_size characters.
                        records.append({
                            # Store the ID of the current DNA sequence.
                            "id": current_id,

                            # Store the extracted DNA chunk.
                            "sequence": dna[i:i + chunk_size]
                        })

                # Remove '>' and get the first word as the DNA sequence ID.
                # Example: '>sequence_1 sample' becomes 'sequence_1'.
                current_id = line[1:].split()[0]

                # Reset the list to collect the new DNA sequence.
                sequence_parts = []

            # If the line is not a header, it contains DNA sequence data.
            else:

                # Add this DNA line to the list of sequence parts.
                sequence_parts.append(line)

        # After the loop, process the final DNA sequence.
        # It does not have another header after it to trigger processing.
        if current_id is not None:

            # Combine all parts of the final DNA sequence.
            # Convert the complete sequence to uppercase.
            dna = "".join(sequence_parts).upper()

            # Split the final DNA sequence into chunks.
            for i in range(0, len(dna), chunk_size):

                # Add each final DNA chunk and its sequence ID to records.
                records.append({
                    # Store the final DNA sequence ID.
                    "id": current_id,

                    # Store the DNA chunk from position i.
                    "sequence": dna[i:i + chunk_size]
                })

    # Return the list containing all DNA chunks.
    return records

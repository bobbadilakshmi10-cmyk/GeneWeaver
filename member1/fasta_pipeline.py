
def prepare_data(file_path, chunk_size=5):
    records = []
    current_id = None
    sequence_parts = []

    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if current_id is not None:
                    dna = "".join(sequence_parts).upper()

                    for i in range(0, len(dna), chunk_size):
                        records.append({
                            "id": current_id,
                            "sequence": dna[i:i + chunk_size]
                        })

                current_id = line[1:].split()[0]
                sequence_parts = []
            else:
                sequence_parts.append(line)

        if current_id is not None:
            dna = "".join(sequence_parts).upper()

            for i in range(0, len(dna), chunk_size):
                records.append({
                    "id": current_id,
                    "sequence": dna[i:i + chunk_size]
                })

    return records

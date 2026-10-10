
from member1.fasta_pipeline import prepare_data

records = prepare_data("data/sample.fasta", chunk_size=5)

assert len(records) > 0, "No DNA records found!"

for record in records:
    assert "id" in record, "DNA ID is missing"
    assert "sequence" in record, "DNA sequence is missing"
    assert len(record["sequence"]) <= 5, "Chunk is too long"

print("All pipeline tests passed!")
print("Total DNA chunks:", len(records))


from member1.fasta_pipeline import prepare_data

filename = "data/sample.fasta"
records = prepare_data(filename, chunk_size=5)

print("GeneWeaver started successfully")

for record in records:
    print("DNA ID:", record["id"])
    print("DNA Sequence:", record["sequence"])

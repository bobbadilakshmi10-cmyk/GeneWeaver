

# Import the prepare_data function from the Member 1 DNA pipeline
from member1.fasta_pipeline import prepare_data

# Store the path of the input FASTA file
filename = "data/sample.fasta"

# Read the FASTA file, process the DNA sequences,
# and split each sequence into chunks of 5 bases
records = prepare_data(filename, chunk_size=5)

# Display a message to confirm that GeneWeaver has started successfully
print("GeneWeaver started successfully")

# Loop through each DNA record returned by the pipeline
for record in records:

    # Print the unique ID or name of the DNA record
    print("DNA ID:", record["id"])

    # Print the DNA sequence chunk stored in the record
    print("DNA Sequence:", record["sequence"])


# Import Path to check whether the FASTA file exists
from pathlib import Path

# Import Rich modules to display colorful output and progress bars
from rich.console import Console
from rich.progress import (
    Progress,
    BarColumn,
    TextColumn,
    TimeElapsedColumn,
)
from rich.table import Table

# Import the function that reads and prepares DNA data
from member1.fasta_pipeline import prepare_data


# Create a Console object to print colorful messages
console = Console()

# Specify the location of the input FASTA file
DATA_FILE = Path("data/sample.fasta")


# Function to display the final DNA analysis results
def display_results(results):
    """Display analysis results supplied by the team."""

    # Create a table for the analysis results
    table = Table(title="GeneWeaver Analysis Results")

    # Add column headings to the table
    table.add_column("DNA ID")
    table.add_column("Position", justify="right")
    table.add_column("DNA Chunk")
    table.add_column("Mismatches", justify="right")
    table.add_column("Score", justify="right")
    table.add_column("Highlighted Bases")

    # Loop through each analysis result
    for result in results:
        # Add the result values as a row in the table
        table.add_row(
            str(result["id"]),
            str(result["position"]),
            str(result["sequence"]),
            str(result["mismatches"]),
            str(result["score"]),
            str(result["highlighted_bases"]),
        )

    # Print the completed analysis table
    console.print(table)


# Main function to run the DNA data preparation pipeline
def main():
    """Run the Member 1 data preparation pipeline."""

    # Display a green success message when the program starts
    console.print(
        "\n[bold green]GeneWeaver Started Successfully[/bold green]"
    )

    # Check whether the input FASTA file exists
    if not DATA_FILE.is_file():
        # Display an error message if the file is missing
        console.print(
            f"[red]FASTA file not found: {DATA_FILE}[/red]"
        )
        return  # Stop the program if the file does not exist

    try:
        # Create a progress bar to show the preparation status
        with Progress(
            TextColumn("{task.description}"),  # Display the task name
            BarColumn(),                       # Display the progress bar
            TextColumn("{task.percentage:>3.0f}%"),  # Show percentage
            TimeElapsedColumn(),               # Show elapsed time
            console=console,                   # Use the Rich console
        ) as progress:

            # Create a progress task with a total value of 1
            task = progress.add_task(
                "Preparing DNA data",
                total=1,
            )

            # Read the FASTA file and split DNA into overlapping chunks
            # chunk_size=5 means each chunk contains up to 5 DNA bases
            # step=1 means the next chunk starts one base later
            records = prepare_data(
                DATA_FILE,
                chunk_size=5,
                step=1,
            )

            # Mark the preparation task as completed
            progress.update(task, advance=1)

        # Display the total number of prepared DNA chunks
        console.print(
            f"\n[green]Prepared {len(records)} DNA chunks.[/green]"
        )

        # Create a table to display the prepared DNA data
        table = Table(title="Prepared DNA Chunks")

        # Add column headings
        table.add_column("DNA ID")
        table.add_column("Position", justify="right")
        table.add_column("DNA Chunk")
        table.add_column("Encoded Values")
        table.add_column("Length", justify="right")

        # Loop through every prepared DNA chunk
        for record in records:
            # Add the DNA details as a row in the table
            table.add_row(
                record["id"],  # DNA record identifier
                str(record["position"]),  # Starting position of the chunk
                record["sequence"],  # DNA sequence, such as ATGCG
                str(record["encoded"].tolist()),  # Numeric DNA encoding
                str(record["length"]),  # Number of bases in the chunk
            )

        # Display the prepared DNA chunks table
        console.print(table)

        # Inform the user that analysis results are not integrated yet
        console.print(
            "\n[yellow]Waiting for the analysis module to provide "
            "mismatch counts, scores and highlighted bases.[/yellow]"
        )

        # After integrating the analysis module, pass its results here:
        # display_results(analysis_results)

    # Handle invalid DNA data or file input/output errors
    except (ValueError, OSError) as error:
        # Display the error in red
        console.print(f"[bold red]Error: {error}[/bold red]")


# Run main() only when this file is executed directly
if __name__ == "__main__":
    main()


from pathlib import Path

from rich.console import Console
from rich.progress import (
    Progress,
    BarColumn,
    TextColumn,
    TimeElapsedColumn,
)
from rich.table import Table

from member1.fasta_pipeline import prepare_data


console = Console()
DATA_FILE = Path("data/sample.fasta")


def display_results(results):
    """Display analysis results supplied by the team."""

    table = Table(title="GeneWeaver Analysis Results")

    table.add_column("DNA ID")
    table.add_column("Position", justify="right")
    table.add_column("DNA Chunk")
    table.add_column("Mismatches", justify="right")
    table.add_column("Score", justify="right")
    table.add_column("Highlighted Bases")

    for result in results:
        table.add_row(
            str(result["id"]),
            str(result["position"]),
            str(result["sequence"]),
            str(result["mismatches"]),
            str(result["score"]),
            str(result["highlighted_bases"]),
        )

    console.print(table)


def main():
    """Run the Member 1 data preparation pipeline."""

    console.print(
        "\n[bold green]GeneWeaver Started Successfully[/bold green]"
    )

    if not DATA_FILE.is_file():
        console.print(
            f"[red]FASTA file not found: {DATA_FILE}[/red]"
        )
        return

    try:
        # The first progress task represents pipeline preparation.
        with Progress(
            TextColumn("{task.description}"),
            BarColumn(),
            TextColumn("{task.percentage:>3.0f}%"),
            TimeElapsedColumn(),
            console=console,
        ) as progress:

            task = progress.add_task(
                "Preparing DNA data",
                total=1,
            )

            records = prepare_data(
                DATA_FILE,
                chunk_size=5,
                step=1,
            )

            progress.update(task, advance=1)

        console.print(
            f"\n[green]Prepared {len(records)} DNA chunks.[/green]"
        )

        # Display the prepared data before analysis integration.
        table = Table(title="Prepared DNA Chunks")

        table.add_column("DNA ID")
        table.add_column("Position", justify="right")
        table.add_column("DNA Chunk")
        table.add_column("Encoded Values")
        table.add_column("Length", justify="right")

        for record in records:
            table.add_row(
                record["id"],
                str(record["position"]),
                record["sequence"],
                str(record["encoded"].tolist()),
                str(record["length"]),
            )

        console.print(table)

        console.print(
            "\n[yellow]Waiting for the analysis module to provide "
            "mismatch counts, scores and highlighted bases.[/yellow]"
        )

        # Later, after analysis integration, call:
        # display_results(analysis_results)

    except (ValueError, OSError) as error:
        console.print(f"[bold red]Error: {error}[/bold red]")


if __name__ == "__main__":
    main()

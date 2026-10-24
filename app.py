
from pathlib import Path
from io import StringIO

from flask import Flask, jsonify, render_template, request
from Bio import SeqIO

from member1.fasta_pipeline import clean_dna
from analysis.scoring import scan_sequences


# Initialize Flask application.
app = Flask(__name__)

# Project folders.
BASE_DIR = Path(__file__).resolve().parent
UPLOAD_FOLDER = BASE_DIR / "data"

# Maximum uploaded file size: 5 MB.
app.config["MAX_CONTENT_LENGTH"] = 5 * 1024 * 1024


@app.route("/")
def home():
    """Open the GeneWeaver frontend."""
    return render_template("index.html")


@app.route("/api/health")
def health():
    """Check whether the application is running."""
    return jsonify({
        "status": "success",
        "message": "GeneWeaver is running"
    })


@app.route("/api/analyze", methods=["POST"])
def analyze():
    """Read a FASTA file and scan for similar DNA sequences."""

    try:
        # Get the uploaded file and form values.
        uploaded_file = request.files.get("fasta_file")
        target = request.form.get("target", "").strip().upper()

        max_mismatches = int(
            request.form.get("max_mismatches", 3)
        )

        step = int(request.form.get("step", 1))

        # Validate uploaded file.
        if uploaded_file is None or not uploaded_file.filename:
            return jsonify({
                "error": "Please upload a FASTA file."
            }), 400

        allowed_extensions = {".fa", ".fasta", ".fna"}
        extension = Path(uploaded_file.filename).suffix.lower()

        if extension not in allowed_extensions:
            return jsonify({
                "error": "Upload a .fa, .fasta, or .fna file."
            }), 400

        # Validate target DNA.
        if not target:
            return jsonify({
                "error": "Enter a target DNA sequence."
            }), 400

        target = clean_dna(target)

        if not target or set(target) - set("ACGT"):
            return jsonify({
                "error": "The target must contain only A, C, G, and T."
            }), 400

        # Validate mismatch count and scan step.
        if max_mismatches < 0 or max_mismatches > len(target):
            return jsonify({
                "error": "Mismatch count must be between 0 and target length."
            }), 400

        if step < 1:
            return jsonify({
                "error": "Step must be at least 1."
            }), 400

        # Read the uploaded file.
        file_text = uploaded_file.read().decode("utf-8-sig")

        # Remove blank lines and semicolon-prefixed comment lines.
        cleaned_lines = []

        for line in file_text.splitlines():
            line = line.strip()

            if not line:
                continue

            if line.startswith(";"):
                continue

            cleaned_lines.append(line)

        cleaned_text = "\n".join(cleaned_lines)

        # Parse the cleaned FASTA content.
        records = list(
            SeqIO.parse(StringIO(cleaned_text), "fasta")
        )

        if not records:
            return jsonify({
                "error": (
                    "No valid FASTA sequences were found. "
                    "Check that the file contains headers beginning with >."
                )
            }), 400

        # Prepare DNA sequences for scanning.
        sequences = []

        for record in records:
            sequence = clean_dna(str(record.seq))

            if sequence:
                sequences.append({
                    "id": record.id,
                    "sequence": sequence
                })

        if not sequences:
            return jsonify({
                "error": "The FASTA file contains no usable sequences."
            }), 400

        # Run the existing GeneWeaver scoring function.
        results = scan_sequences(
            sequences=sequences,
            target=target,
            max_mismatches=max_mismatches,
            step=step
        )

        # Limit results returned to the frontend.
        displayed_results = results[:5000]

        return jsonify({
            "status": "success",
            "sequence_count": len(sequences),
            "results_count": len(results),
            "results": displayed_results,
            "message": "DNA scan completed."
        })

    except (ValueError, UnicodeDecodeError) as error:
        return jsonify({
            "error": str(error)
        }), 400

    except Exception:
        app.logger.exception("GeneWeaver analysis failed")

        return jsonify({
            "error": (
                "Analysis failed. Check the terminal for "
                "the detailed error message."
            )
        }), 500


@app.errorhandler(413)
def file_too_large(error):
    """Handle files larger than 5 MB."""
    return jsonify({
        "error": "File is too large. Maximum size is 5 MB."
    }), 413


if __name__ == "__main__":
    app.run(debug=True)

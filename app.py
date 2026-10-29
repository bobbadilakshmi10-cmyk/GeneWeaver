
# Import Path to work with file and folder paths.
from pathlib import Path

# Import StringIO to treat a string as a file for FASTA parsing.
from io import StringIO

# Import Flask functions to create the application, handle requests,
# return JSON responses, and display HTML templates.
from flask import Flask, jsonify, render_template, request

# Import BioPython's SeqIO module to read FASTA DNA sequence files.
from Bio import SeqIO

# Import the function that cleans and validates DNA sequences.
from member1.fasta_pipeline import clean_dna

# Import the function that scans DNA sequences for target matches.
from analysis.scoring import scan_sequences


# --------------------------------------------------
# 1. INITIALIZE THE FLASK APPLICATION
# --------------------------------------------------

# Create the Flask application instance.
app = Flask(__name__)

# Get the absolute path of the current project directory.
BASE_DIR = Path(__file__).resolve().parent

# Define the folder used to store project data.
UPLOAD_FOLDER = BASE_DIR / "data"

# Set the maximum allowed request size to 5 MB.
app.config["MAX_CONTENT_LENGTH"] = 5 * 1024 * 1024


# --------------------------------------------------
# 2. HOME PAGE ROUTE
# --------------------------------------------------

@app.route("/")
def home():
    """Display the GeneWeaver frontend."""

    # Load index.html from the templates folder.
    return render_template("index.html")


# --------------------------------------------------
# 3. HEALTH CHECK API
# --------------------------------------------------

@app.route("/api/health")
def health():
    """Check whether the Flask application is running."""

    # Return the application status as a JSON response.
    return jsonify({
        "status": "success",
        "message": "GeneWeaver is running"
    })


# --------------------------------------------------
# 4. DNA ANALYSIS API
# --------------------------------------------------

# Accept POST requests containing a FASTA file and form data.
@app.route("/api/analyze", methods=["POST"])
def analyze():
    """Read a FASTA file and scan for similar DNA sequences."""

    try:
        # Retrieve the uploaded FASTA file from the request.
        uploaded_file = request.files.get("fasta_file")

        # Read the target DNA sequence, remove surrounding spaces,
        # and convert all letters to uppercase.
        target = request.form.get("target", "").strip().upper()

        # Get the maximum allowed number of mismatches.
        # Use 3 as the default value if the field is missing.
        max_mismatches = int(
            request.form.get("max_mismatches", 3)
        )

        # Get the scanning step size.
        # A step of 1 checks every possible starting position.
        step = int(request.form.get("step", 1))

        # ------------------------------------------
        # 4.1 VALIDATE THE UPLOADED FILE
        # ------------------------------------------

        # Check whether the user selected a file.
        if uploaded_file is None or not uploaded_file.filename:
            return jsonify({
                "error": "Please upload a FASTA file."
            }), 400

        # Define the supported FASTA file extensions.
        allowed_extensions = {".fa", ".fasta", ".fna"}

        # Extract the uploaded file's extension.
        extension = Path(uploaded_file.filename).suffix.lower()

        # Reject files with unsupported extensions.
        if extension not in allowed_extensions:
            return jsonify({
                "error": "Upload a .fa, .fasta, or .fna file."
            }), 400

        # ------------------------------------------
        # 4.2 VALIDATE THE TARGET DNA SEQUENCE
        # ------------------------------------------

        # Ensure the user entered a target sequence.
        if not target:
            return jsonify({
                "error": "Enter a target DNA sequence."
            }), 400

        # Clean the target sequence using the pipeline function.
        target = clean_dna(target)

        # Allow only A, C, G, and T in the target sequence.
        # Reject empty sequences and sequences containing
        # characters other than these four DNA bases.
        if not target or set(target) - set("ACGT"):
            return jsonify({
                "error": "The target must contain only A, C, G, and T."
            }), 400

        # ------------------------------------------
        # 4.3 VALIDATE SCANNING PARAMETERS
        # ------------------------------------------

        # The mismatch count cannot be negative or greater
        # than the length of the target sequence.
        if max_mismatches < 0 or max_mismatches > len(target):
            return jsonify({
                "error": "Mismatch count must be between 0 and target length."
            }), 400

        # The step size must be at least 1.
        if step < 1:
            return jsonify({
                "error": "Step must be at least 1."
            }), 400

        # ------------------------------------------
        # 4.4 READ THE FASTA FILE
        # ------------------------------------------

        # Read the uploaded file and decode it as UTF-8.
        # utf-8-sig also removes a UTF-8 byte-order mark if present.
        file_text = uploaded_file.read().decode("utf-8-sig")

        # Create an empty list to store cleaned FASTA lines.
        cleaned_lines = []

        # Process the uploaded file one line at a time.
        for line in file_text.splitlines():

            # Remove spaces from the beginning and end.
            line = line.strip()

            # Ignore empty lines.
            if not line:
                continue

            # Ignore comment lines beginning with a semicolon.
            if line.startswith(";"):
                continue

            # Keep all other lines for FASTA parsing.
            cleaned_lines.append(line)

        # Combine the cleaned lines into one string.
        cleaned_text = "\n".join(cleaned_lines)

        # ------------------------------------------
        # 4.5 PARSE THE FASTA RECORDS
        # ------------------------------------------

        # Convert the cleaned text into a file-like object.
        # SeqIO.parse reads the DNA records in FASTA format.
        records = list(
            SeqIO.parse(StringIO(cleaned_text), "fasta")
        )

        # Check whether any FASTA records were detected.
        if not records:
            return jsonify({
                "error": (
                    "No valid FASTA sequences were found. "
                    "Check that the file contains headers beginning with >."
                )
            }), 400

        # ------------------------------------------
        # 4.6 PREPARE THE DNA SEQUENCES
        # ------------------------------------------

        # Create a list to store cleaned DNA records.
        sequences = []

        # Process each record obtained from the FASTA file.
        for record in records:

            # Convert the sequence to a string and clean it.
            sequence = clean_dna(str(record.seq))

            # Add non-empty sequences to the scanning list.
            if sequence:
                sequences.append({
                    "id": record.id,       # FASTA record identifier.
                    "sequence": sequence   # Cleaned DNA sequence.
                })

        # Check whether any usable DNA sequences remain.
        if not sequences:
            return jsonify({
                "error": "The FASTA file contains no usable sequences."
            }), 400

        # ------------------------------------------
        # 4.7 SCAN DNA SEQUENCES
        # ------------------------------------------

        # Call the existing scoring module to find target matches.
        # The function receives the sequences, target DNA,
        # allowed mismatches, and scanning step.
        results = scan_sequences(
            sequences=sequences,
            target=target,
            max_mismatches=max_mismatches,
            step=step
        )

        # ------------------------------------------
        # 4.8 PREPARE RESULTS FOR THE FRONTEND
        # ------------------------------------------

        # Limit the number of results sent to the browser
        # to avoid returning an excessively large response.
        displayed_results = results[:5000]

        # Return the scan status, sequence count, total matches,
        # displayed results, and a success message as JSON.
        return jsonify({
            "status": "success",
            "sequence_count": len(sequences),
            "results_count": len(results),
            "results": displayed_results,
            "message": "DNA scan completed."
        })

    # ----------------------------------------------
    # 4.9 HANDLE EXPECTED INPUT ERRORS
    # ----------------------------------------------

    # Handle invalid integer values or files that cannot
    # be decoded using UTF-8.
    except (ValueError, UnicodeDecodeError) as error:
        return jsonify({
            "error": str(error)
        }), 400

    # ----------------------------------------------
    # 4.10 HANDLE UNEXPECTED ERRORS
    # ----------------------------------------------

    except Exception:
        # Write the full exception traceback to the terminal
        # or configured application logs for debugging.
        app.logger.exception("GeneWeaver analysis failed")

        # Return a general error message to the frontend
        # without exposing internal technical details.
        return jsonify({
            "error": (
                "Analysis failed. Check the terminal for "
                "the detailed error message."
            )
        }), 500


# --------------------------------------------------
# 5. HANDLE FILES LARGER THAN 5 MB
# --------------------------------------------------

@app.errorhandler(413)
def file_too_large(error):
    """Handle requests that exceed the configured size limit."""

    # Return an error response when an uploaded file is too large.
    return jsonify({
        "error": "File is too large. Maximum size is 5 MB."
    }), 413


# --------------------------------------------------
# 6. RUN THE FLASK APPLICATION
# --------------------------------------------------

# Run the development server only when this file is executed directly.
if __name__ == "__main__":

    # Enable debug mode to help identify errors during development.
    app.run(debug=True)

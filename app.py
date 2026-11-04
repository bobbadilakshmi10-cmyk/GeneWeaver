# ============================================================
# GENEWEAVER - FLASK APPLICATION
# This file handles the web interface, FASTA file processing,
# DNA sequence validation, and DNA target sequence analysis.
# ============================================================


# Import Path to work with project file and folder locations.
from pathlib import Path

# Import StringIO to convert text into a file-like object.
# This allows BioPython to read uploaded FASTA text.
from io import StringIO

# Import Flask components:
# Flask          -> Creates the web application.
# jsonify        -> Sends data as JSON responses.
# render_template -> Displays HTML pages.
# request        -> Receives files and form data from the user.
from flask import Flask, jsonify, render_template, request

# Import BioPython's SeqIO module to parse FASTA files.
from Bio import SeqIO

# Import the custom function that cleans DNA sequences.
from member1.fasta_pipeline import clean_dna

# Import the function that scans sequences and finds
# matches based on the target and allowed mismatches.
from analysis.scoring import scan_sequences


# ============================================================
# 1. FLASK APPLICATION CONFIGURATION
# ============================================================

# Create the main Flask application object.
app = Flask(__name__)

# Find the absolute path of the folder containing this app.py file.
BASE_DIR = Path(__file__).resolve().parent

# Define the project's data folder.
# This can be used to store DNA or project-related data files.
UPLOAD_FOLDER = BASE_DIR / "data"

# Limit the maximum request/upload size to 5 MB.
# This prevents very large files from being uploaded.
app.config["MAX_CONTENT_LENGTH"] = 5 * 1024 * 1024


# ============================================================
# 2. HOME PAGE ROUTE
# ============================================================

# When the user opens the main website URL (/),
# this function displays the GeneWeaver web page.
@app.route("/")
def home():
    """Display the GeneWeaver home page."""

    # Load index.html from the templates folder.
    return render_template("index.html")


# ============================================================
# 3. APPLICATION HEALTH CHECK API
# ============================================================

# This API checks whether the Flask application is running properly.
@app.route("/api/health")
def health():
    """Return the current application status."""

    # Send a success message in JSON format.
    return jsonify({
        "status": "success",
        "message": "GeneWeaver is running"
    })


# ============================================================
# 4. DNA ANALYSIS API
# ============================================================

# This route accepts a POST request from the frontend.
# The request contains:
# 1. A FASTA DNA file.
# 2. A target DNA sequence.
# 3. Maximum allowed mismatches.
# 4. Scanning step size.
@app.route("/api/analyze", methods=["POST"])
def analyze():
    """Read a FASTA file and search for target DNA matches."""

    try:
        # ====================================================
        # 4.1 GET DATA FROM THE USER REQUEST
        # ====================================================

        # Get the uploaded FASTA file using its form field name.
        uploaded_file = request.files.get("fasta_file")

        # Get the target DNA sequence entered by the user.
        # strip() removes extra spaces.
        # upper() converts lowercase DNA letters to uppercase.
        target = request.form.get("target", "").strip().upper()

        # Get the maximum number of mismatches allowed.
        # If the user does not provide a value, use 3 by default.
        max_mismatches = int(
            request.form.get("max_mismatches", 3)
        )

        # Get the scanning step size.
        # Step = 1 means every possible position is checked.
        step = int(request.form.get("step", 1))

        # ====================================================
        # 4.2 VALIDATE THE UPLOADED FASTA FILE
        # ====================================================

        # Check whether the user selected and uploaded a file.
        if uploaded_file is None or not uploaded_file.filename:
            return jsonify({
                "error": "Please upload a FASTA file."
            }), 400

        # Define the file extensions supported by the application.
        allowed_extensions = {".fa", ".fasta", ".fna"}

        # Extract the file extension and convert it to lowercase.
        extension = Path(uploaded_file.filename).suffix.lower()

        # Reject the file if it is not a supported FASTA format.
        if extension not in allowed_extensions:
            return jsonify({
                "error": "Upload a .fa, .fasta, or .fna file."
            }), 400

        # ====================================================
        # 4.3 VALIDATE THE TARGET DNA SEQUENCE
        # ====================================================

        # Check whether the user entered a target DNA sequence.
        if not target:
            return jsonify({
                "error": "Enter a target DNA sequence."
            }), 400

        # Clean the target DNA using the custom pipeline function.
        # This helps standardize the sequence before analysis.
        target = clean_dna(target)

        # Check that the target contains only valid DNA bases.
        # A -> Adenine
        # C -> Cytosine
        # G -> Guanine
        # T -> Thymine
        if not target or set(target) - set("ACGT"):
            return jsonify({
                "error": "The target must contain only A, C, G, and T."
            }), 400

        # ====================================================
        # 4.4 VALIDATE SCANNING SETTINGS
        # ====================================================

        # The mismatch value cannot be negative.
        # It also cannot be greater than the target sequence length.
        if max_mismatches < 0 or max_mismatches > len(target):
            return jsonify({
                "error": "Mismatch count must be between 0 and target length."
            }), 400

        # The scanning step must be 1 or greater.
        if step < 1:
            return jsonify({
                "error": "Step must be at least 1."
            }), 400

        # ====================================================
        # 4.5 READ AND CLEAN THE UPLOADED FASTA FILE
        # ====================================================

        # Read the uploaded file and decode it into text.
        # utf-8-sig also handles files containing a UTF-8
        # byte-order mark (BOM).
        file_text = uploaded_file.read().decode("utf-8-sig")

        # Create an empty list to store valid FASTA lines.
        cleaned_lines = []

        # Process the uploaded file one line at a time.
        for line in file_text.splitlines():

            # Remove unnecessary spaces from the line.
            line = line.strip()

            # Ignore blank lines in the uploaded file.
            if not line:
                continue

            # Ignore FASTA comment lines beginning with a semicolon.
            if line.startswith(";"):
                continue

            # Keep valid headers and DNA sequence lines.
            cleaned_lines.append(line)

        # Combine all cleaned lines into one formatted text string.
        cleaned_text = "\n".join(cleaned_lines)

        # ====================================================
        # 4.6 PARSE FASTA DNA RECORDS
        # ====================================================

        # Convert the FASTA text into a file-like object.
        # BioPython's SeqIO.parse() reads each FASTA record.
        records = list(
            SeqIO.parse(StringIO(cleaned_text), "fasta")
        )

        # Check whether any valid FASTA records were found.
        # FASTA records normally begin with a header using ">".
        if not records:
            return jsonify({
                "error": (
                    "No valid FASTA sequences were found. "
                    "Check that the file contains headers beginning with >."
                )
            }), 400

        # ====================================================
        # 4.7 PREPARE DNA SEQUENCES FOR ANALYSIS
        # ====================================================

        # Create an empty list to store cleaned DNA records.
        sequences = []

        # Process every DNA record from the FASTA file.
        for record in records:

            # Convert the BioPython sequence object to a string
            # and clean the DNA sequence.
            sequence = clean_dna(str(record.seq))

            # Add only valid, non-empty sequences to the list.
            if sequence:
                sequences.append({
                    # Store the FASTA sequence identifier.
                    "id": record.id,

                    # Store the cleaned DNA sequence.
                    "sequence": sequence
                })

        # Check whether usable DNA sequences are available.
        if not sequences:
            return jsonify({
                "error": "The FASTA file contains no usable sequences."
            }), 400

        # ====================================================
        # 4.8 SCAN DNA SEQUENCES
        # ====================================================

        # Call the scoring module to search for DNA sequences
        # that match the target sequence.
        #
        # sequences      -> DNA sequences from the FASTA file.
        # target         -> DNA pattern entered by the user.
        # max_mismatches -> Maximum differences allowed.
        # step           -> Number of positions to move after each scan.
        results = scan_sequences(
            sequences=sequences,
            target=target,
            max_mismatches=max_mismatches,
            step=step
        )

        # ====================================================
        # 4.9 PREPARE RESULTS FOR THE FRONTEND
        # ====================================================

        # Limit the results sent to the browser to 5000.
        # This helps prevent very large API responses.
        displayed_results = results[:5000]

        # Return the DNA analysis results in JSON format.
        return jsonify({
            "status": "success",

            # Number of valid DNA sequences processed.
            "sequence_count": len(sequences),

            # Total number of matches found.
            "results_count": len(results),

            # Results displayed to the user.
            "results": displayed_results,

            # Success message.
            "message": "DNA scan completed."
        })

    # ========================================================
    # 4.10 HANDLE EXPECTED USER INPUT ERRORS
    # ========================================================

    # Handle invalid number input and file decoding problems.
    except (ValueError, UnicodeDecodeError) as error:
        return jsonify({
            "error": str(error)
        }), 400

    # ========================================================
    # 4.11 HANDLE UNEXPECTED APPLICATION ERRORS
    # ========================================================

    except Exception:
        # Save the complete error information in the terminal/log.
        # This helps developers identify and fix problems.
        app.logger.exception("GeneWeaver analysis failed")

        # Send a safe error message to the user without exposing
        # internal application details.
        return jsonify({
            "error": (
                "Analysis failed. Check the terminal for "
                "the detailed error message."
            )
        }), 500


# ============================================================
# 5. HANDLE FILES LARGER THAN THE ALLOWED SIZE
# ============================================================

# This function is automatically called when a request
# exceeds the configured 5 MB size limit.
@app.errorhandler(413)
def file_too_large(error):
    """Return an error when the uploaded file is too large."""

    # Send a clear error message to the frontend.
    return jsonify({
        "error": "File is too large. Maximum size is 5 MB."
    }), 413


# ============================================================
# 6. START THE FLASK APPLICATION
# ============================================================

# Run the application only when this file is executed directly.
# This prevents the development server from starting when
# app.py is imported into another Python file.
if __name__ == "__main__":

    # Start the Flask development server.
    # Debug mode helps display detailed errors during development.
    app.run(debug=True)
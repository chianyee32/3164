# File: connection.py

import os
import sys
import subprocess
from flask import Flask, request, jsonify, send_from_directory
from flask_cors import CORS

# ─── 1. Configuration ──────────────────────────────────────────────────────────
BASE_DIR = os.path.dirname(__file__)

# Create folders for uploads and predictions
UPLOAD_TEMP_DIR = os.path.join(BASE_DIR, "uploads_temp")
os.makedirs(UPLOAD_TEMP_DIR, exist_ok=True)
PREDICTION_FOLDER = os.path.join(BASE_DIR, "uploads_prediction")
os.makedirs(PREDICTION_FOLDER, exist_ok=True)

# ─── 2. Flask App Initialization ───────────────────────────────────────────────
app = Flask(
    __name__,
    static_folder=BASE_DIR,    
    static_url_path=""         
)
app.config["UPLOAD_FOLDER"] = UPLOAD_TEMP_DIR
app.config["PREDICTION_FOLDER"] = PREDICTION_FOLDER
CORS(app, resources={r"/*": {"origins": "*"}})

# ─── 3. Static File Routes ─────────────────────────────────────────────────────
@app.route('/')
def index():
    return app.send_static_file('introduction.html')

# ─── 4. API Routes ──────────────────────────────────────────────────────────────
@app.route('/upload_prediction', methods=['POST'])
def upload_prediction():
    if "file" not in request.files:
        return jsonify({"error": "No file uploaded"}), 400
    f = request.files["file"]
    if f.filename == "":
        return jsonify({"error": "No file selected"}), 400

    # save incoming CSV
    in_path = os.path.join(app.config["UPLOAD_FOLDER"], f.filename)
    f.save(in_path)

    try:
        # Run your predictor using the same Python interpreter
        python_exe = sys.executable
        script = os.path.join(BASE_DIR, "predict_chemoresistance.py")
        proc = subprocess.run(
            [python_exe, script, in_path],
            capture_output=True,
            text=True
        )

        print("Prediction stdout:", proc.stdout)
        print("Prediction stderr:", proc.stderr)
        if proc.returncode != 0:
            return jsonify({
                "error": "Prediction failed.",
                "details": proc.stderr.strip()
            }), 400

        # move output CSV into uploads_prediction
        src = os.path.join(BASE_DIR, "predictions_output.csv")
        dst = os.path.join(app.config["PREDICTION_FOLDER"], "predictions_output.csv")
        if not os.path.exists(src):
            return jsonify({"error": "Output file not found"}), 500
        os.replace(src, dst)

        return jsonify({
            "message": "Prediction succeeded",
            "download_url": f"/download_prediction/predictions_output.csv"
        })

    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route('/download_prediction/<filename>')
def download_prediction_file(filename):
    # serves files out of uploads_prediction
    return send_from_directory(
        app.config["PREDICTION_FOLDER"],
        filename,
        as_attachment=True
    )


# ─── 5. Run Server ─────────────────────────────────────────────────────────────
if __name__ == '__main__':
    app.run(debug=True, port=5000)

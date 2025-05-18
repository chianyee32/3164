# File: connection.py

import os
import sys
import subprocess
import requests
from flask import Flask, request, jsonify, send_from_directory
from flask_cors import CORS

# ─── 1. Configuration ──────────────────────────────────────────────────────────
BASE_DIR = os.path.dirname(__file__)

# S3-hosted assets; env vars point to URLs, filename keys map to local filenames
ASSETS = {
    'MODEL_URL': '4_cancers_chemoresistance_model.h5',
    'SCALER_URL': 'scaler.pkl',
    'CCLE_TRANSCRIPT_URL': 'CCLE_Transcriptomics_cleaned.csv',
    'CCLE_PROTEIN_URL': 'CCLE_Proteomics_cleaned.csv',
    'DRUG_INFO_URL': 'cleaned_drug_info.csv',
    'PUBCHEM_URL': 'cleaned_pubchem.csv',
}

# Download any missing assets at startup
for env_var, filename in ASSETS.items():
    url = os.environ.get(env_var)
    dest = os.path.join(BASE_DIR, filename)
    if url and not os.path.exists(dest):
        print(f"Downloading {env_var} from {url} to {dest}...")
        resp = requests.get(url, stream=True)
        resp.raise_for_status()
        with open(dest, 'wb') as f:
            for chunk in resp.iter_content(8192):
                f.write(chunk)
        print(f"Downloaded {filename}.")

# Create folders for user uploads and predictions
UPLOAD_TEMP_DIR = os.path.join(BASE_DIR, "uploads_temp")
os.makedirs(UPLOAD_TEMP_DIR, exist_ok=True)
PREDICTION_FOLDER = os.path.join(BASE_DIR, "uploads_prediction")
os.makedirs(PREDICTION_FOLDER, exist_ok=True)

# ─── 2. Flask App Initialization ───────────────────────────────────────────────
app = Flask(__name__, static_folder=BASE_DIR, static_url_path="")
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
    if 'file' not in request.files:
        return jsonify({'error': 'No file uploaded'}), 400
    f = request.files['file']
    if f.filename == '':
        return jsonify({'error': 'No file selected'}), 400

    # Save incoming CSV
    in_path = os.path.join(app.config['UPLOAD_FOLDER'], f.filename)
    f.save(in_path)

    try:
        python_exe = sys.executable
        script = os.path.join(BASE_DIR, 'predict_chemoresistance.py')
        proc = subprocess.run([python_exe, script, in_path], capture_output=True, text=True)

        print('Prediction stdout:', proc.stdout)
        print('Prediction stderr:', proc.stderr)
        if proc.returncode != 0:
            return jsonify({'error': 'Prediction failed.', 'details': proc.stderr.strip()}), 400

        # Move output CSV
        src = os.path.join(BASE_DIR, 'predictions_output.csv')
        dst = os.path.join(app.config['PREDICTION_FOLDER'], 'predictions_output.csv')
        if not os.path.exists(src):
            return jsonify({'error': 'Output file not found'}), 500
        os.replace(src, dst)

        return jsonify({'message': 'Prediction succeeded', 'download_url': '/download_prediction/predictions_output.csv'})
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/download_prediction/<filename>')
def download_prediction_file(filename):
    return send_from_directory(app.config['PREDICTION_FOLDER'], filename, as_attachment=True)

# ─── 5. Run Server ─────────────────────────────────────────────────────────────
if __name__ == '__main__':
    port = int(os.environ.get('PORT', 5000))
    app.run(host='0.0.0.0', port=port)

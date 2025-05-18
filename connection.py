# File: connection.py

import os, sys, subprocess
from flask import Flask, request, jsonify, send_from_directory
from flask_cors import CORS

BASE_DIR         = os.path.dirname(__file__)
UPLOAD_FOLDER    = os.path.join(BASE_DIR, "uploads_temp")
PREDICTION_FOLDER= os.path.join(BASE_DIR, "uploads_prediction")
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(PREDICTION_FOLDER, exist_ok=True)

app = Flask(__name__, static_folder=BASE_DIR, static_url_path="")
CORS(app)

@app.route('/')
def index():
    return app.send_static_file('introduction.html')

@app.route('/upload_prediction', methods=['POST'])
def upload_prediction():
    f = request.files.get('file')
    if not f or f.filename == '':
        return jsonify(error="No file uploaded"), 400

    in_path = os.path.join(UPLOAD_FOLDER, f.filename)
    f.save(in_path)

    # invoke our predictor
    proc = subprocess.run(
        [sys.executable, os.path.join(BASE_DIR, "predict_chemoresistance.py"), in_path],
        capture_output=True, text=True
    )
    print(proc.stdout, proc.stderr, file=sys.stderr)

    if proc.returncode != 0:
        return jsonify(error="Prediction failed", details=proc.stderr), 400

    # move out
    out_src = os.path.join(BASE_DIR, "predictions_output.csv")
    out_dst = os.path.join(PREDICTION_FOLDER, "predictions_output.csv")
    if not os.path.exists(out_src):
        return jsonify(error="Output missing"), 500
    os.replace(out_src, out_dst)

    return jsonify(
        message="Prediction succeeded",
        download_url=f"/download_prediction/predictions_output.csv"
    )

@app.route('/download_prediction/<name>')
def download_prediction(name):
    return send_from_directory(PREDICTION_FOLDER, name, as_attachment=True)

if __name__ == '__main__':
    # bind to $PORT on Heroku
    port = int(os.environ.get("PORT", 5000))
    app.run(host="0.0.0.0", port=port)

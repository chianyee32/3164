# File: predict_chemoresistance.py

import os
import sys
import traceback
import pandas as pd
import joblib
from tensorflow.keras.models import load_model
from preprocessing_model import preprocess_user_dataset, preprocess_new_data

# Ensure console uses UTF-8 encoding so any emojis or special chars from imports don't error
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
    sys.stderr.reconfigure(encoding='utf-8')

# Base directory of this script
BASE_DIR = os.path.dirname(__file__)

# Absolute paths for model, scaler, and data files
MODEL_PATH       = os.path.join(BASE_DIR, '4_cancers_chemoresistance_model.h5')
SCALER_PATH      = os.path.join(BASE_DIR, 'scaler.pkl')
PREPROCESSED_CSV = os.path.join(BASE_DIR, 'user_preprocessed_output.csv')
OUTPUT_CSV       = os.path.join(BASE_DIR, 'predictions_output.csv')

# Suppress TensorFlow INFO and oneDNN messages
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
os.environ['TF_ENABLE_ONEDNN_OPTS']  = '0'


def load_pipeline():
    """
    Load the trained Keras model and the scaler object.
    """
    model  = load_model(MODEL_PATH)
    scaler = joblib.load(SCALER_PATH)
    return model, scaler


def make_predictions(model, X_scaled):
    """
    Generate raw predictions (LN_IC50 values) and return a flattened array.
    """
    preds = model.predict(X_scaled)
    return preds.flatten()


def save_predictions(original_data: pd.DataFrame, predictions, output_file: str = None):
    """
    Combine original metadata with predictions, classify sensitivity,
    and write to CSV.
    """
    if output_file is None:
        output_file = OUTPUT_CSV

    # Build results DataFrame
    results_df = original_data[
        ['DRUG_ID', 'DRUG_NAME', 'COSMIC_ID', 'CCLE_Name', 'CANCER_TYPE']
    ].copy()
    results_df['Predicted_LN_IC50'] = predictions

    # Label sensitivity
    def classify(ln):
        if ln < 2.36:
            return 'High'
        elif ln <= 5.26:
            return 'Intermediate'
        else:
            return 'Low'

    results_df['Sensitivity'] = results_df['Predicted_LN_IC50'].apply(classify)

    # Remove duplicates and sort
    before = len(results_df)
    results_df = (
        results_df.drop_duplicates(
            subset=['DRUG_ID','COSMIC_ID','CCLE_Name','DRUG_NAME','CANCER_TYPE'],
            keep='first'
        )
        .sort_values(by='Predicted_LN_IC50', ascending=True)
    )
    after = len(results_df)
    if before != after:
        print(f"Dropped {before - after} duplicate rows.")

    # Write output
    results_df.to_csv(output_file, index=False)
    print(f"Predictions saved to {output_file}")


def main(input_file_path: str):
    """
    End-to-end pipeline: preprocess input, run model, and save output.
    """
    try:
        print(f"Starting prediction for: {input_file_path}")

        # 1. Preprocess raw input
        preprocessed_df = preprocess_user_dataset(
            input_file_path,
            output_csv_path=PREPROCESSED_CSV
        )

        # 2. Load model + scaler
        model, scaler = load_pipeline()

        # 3. Scale preprocessed data
        _, X_scaled = preprocess_new_data(PREPROCESSED_CSV, scaler)

        # 4. Predict
        preds = make_predictions(model, X_scaled)

        # 5. Save predictions
        save_predictions(preprocessed_df, preds, output_file=OUTPUT_CSV)

        print("Prediction completed successfully.")

    except Exception:
        # Print full traceback for debugging
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python predict_chemoresistance.py <input_file_path>", file=sys.stderr)
        sys.exit(1)
    main(sys.argv[1])

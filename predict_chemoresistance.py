# File: predict_chemoresistance.py

import os
import sys
import joblib
import requests
import pandas as pd
from tensorflow.keras.models import load_model
from preprocessing_model import preprocess_user_dataset, preprocess_new_data
import traceback

BASE_DIR      = os.path.dirname(__file__)
MODEL_PATH    = os.path.join(BASE_DIR, 'model.h5')
SCALER_PATH   = os.path.join(BASE_DIR, 'scaler.pkl')
PREP_CSV      = os.path.join(BASE_DIR, 'user_preprocessed_output.csv')
OUTPUT_CSV    = os.path.join(BASE_DIR, 'predictions_output.csv')

ENV_MODEL_URL  = os.environ.get('MODEL_URL')
ENV_SCALER_URL = os.environ.get('SCALER_URL')


def _download_if_missing(url: str, dst: str):
    if not os.path.exists(dst):
        resp = requests.get(url, stream=True)
        resp.raise_for_status()
        with open(dst, 'wb') as f:
            for chunk in resp.iter_content(1024*1024):
                f.write(chunk)
        print(f"Downloaded {os.path.basename(dst)}")


def load_pipeline():
    if ENV_MODEL_URL:
        _download_if_missing(ENV_MODEL_URL, MODEL_PATH)
    if ENV_SCALER_URL:
        _download_if_missing(ENV_SCALER_URL, SCALER_PATH)

    model  = load_model(MODEL_PATH)
    scaler = joblib.load(SCALER_PATH)
    return model, scaler


def make_predictions(model, X):
    return model.predict(X).flatten()


def save_predictions(df: pd.DataFrame, preds: pd.Series):
    res = df[['DRUG_ID','DRUG_NAME','COSMIC_ID','CCLE_Name','CANCER_TYPE']].copy()
    res['Predicted_LN_IC50'] = preds
    res['Sensitivity'] = pd.cut(
        res['Predicted_LN_IC50'],
        bins=[-float('inf'),2.36,5.26,float('inf')],
        labels=['High','Intermediate','Low']
    )
    res = res.drop_duplicates(subset=res.columns.tolist()).sort_values('Predicted_LN_IC50')
    res.to_csv(OUTPUT_CSV, index=False)
    print(f"✅ Saved predictions to {OUTPUT_CSV}")


def main(input_path):
    try:
        print(f"📂 Input: {input_path}")
        # 1) preprocess
        prepped = preprocess_user_dataset(input_path, PREP_CSV)

        # 2) load model/scaler (downloads if needed)
        model, scaler = load_pipeline()

        # 3) align & scale
        df_raw, Xs    = preprocess_new_data(PREP_CSV, scaler)

        # 4) predict & save
        preds = make_predictions(model, Xs)
        save_predictions(df_raw, preds)

        print("✅ Done")
        return 0

    except Exception:
        traceback.print_exc(file=sys.stderr)
        return 1


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python predict_chemoresistance.py <input_csv>", file=sys.stderr)
        sys.exit(1)
    sys.exit(main(sys.argv[1]))

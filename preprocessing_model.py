# File: preprocessing_model.py
"""
Preprocessing pipeline that pulls its reference CSVs either from disk
(Preprocess_files/) or, if missing, downloads them from environment URLs.
"""

import os
import pandas as pd
import numpy as np
import requests
from typing import Tuple, Optional
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.DataStructs import ConvertToNumpyArray

BASE_DIR       = os.path.dirname(__file__)
PREPROCESS_DIR = os.path.join(BASE_DIR, 'Preprocess_files')
os.makedirs(PREPROCESS_DIR, exist_ok=True)

ENV_URLS = {
    'CCLE_Transcriptomics_cleaned.csv': os.environ.get('CCLE_TRANSCRIPT_URL'),
    'CCLE_Proteomics_cleaned.csv':     os.environ.get('CCLE_PROTEIN_URL'),
    'cleaned_drug_info.csv':           os.environ.get('DRUG_INFO_URL'),
    'cleaned_pubchem.csv':             os.environ.get('PUBCHEM_URL'),
}


def _ensure_file(name: str) -> str:
    """
    Make sure `name` exists in PREPROCESS_DIR. If missing and an env URL
    is set, download it.
    Returns the full local path.
    """
    local = os.path.join(PREPROCESS_DIR, name)
    if not os.path.exists(local):
        url = ENV_URLS.get(name)
        if not url:
            raise FileNotFoundError(
                f"{name!r} not found locally and no environment URL provided."
            )
        resp = requests.get(url, stream=True)
        resp.raise_for_status()
        with open(local, 'wb') as f:
            for chunk in resp.iter_content(1024*1024):
                f.write(chunk)
        print(f"Downloaded {name} from {url}")
    return local


def preprocess_user_dataset(
    user_csv_path: str,
    output_csv_path: Optional[str] = None
) -> pd.DataFrame:
    """
    1) Validate that user CSV has required transcriptomics & proteomics features.
    2) Merge drug_info & pubchem to get ISOSMILES.
    3) Generate 256‐bit Morgan fingerprints.
    4) Save a flattened CSV for prediction.
    """
    if output_csv_path is None:
        output_csv_path = os.path.join(BASE_DIR, 'user_preprocessed_output.csv')

    df = pd.read_csv(user_csv_path)
    print(f"Loaded user dataset: {user_csv_path}")

    # load reference columns
    trans_ref = pd.read_csv(
        _ensure_file('CCLE_Transcriptomics_cleaned.csv'),
        nrows=0
    ).columns
    prot_ref  = pd.read_csv(
        _ensure_file('CCLE_Proteomics_cleaned.csv'),
        nrows=0
    ).columns

    meta = {'DRUG_ID','DRUG_NAME','CCLE_Name','COSMIC_ID','CANCER_TYPE','LN_IC50'}
    user_feats = set(df.columns) - meta

    missing_t = set(trans_ref) - user_feats
    missing_p = set(prot_ref)  - user_feats
    if missing_t or missing_p:
        msg = "Missing omics features."
        if missing_t: msg += f" Transcriptomics: {list(missing_t)[:5]}…"
        if missing_p: msg += f" Proteomics: {list(missing_p)[:5]}…"
        raise ValueError(msg)

    # load auxiliary tables
    drug_info    = pd.read_csv(_ensure_file('cleaned_drug_info.csv'))
    pubchem_data = pd.read_csv(_ensure_file('cleaned_pubchem.csv'))

    # get ISOSMILES
    if 'ISOSMILES' in df:
        merged = df.copy()
    elif 'PubCHEM' in df:
        merged = df.merge(pubchem_data, on='PubCHEM', how='left')
    else:
        tmp    = df.merge(drug_info, on='DRUG_ID', how='left')
        tmp    = tmp.dropna(subset=['PubCHEM'])
        merged = tmp.merge(pubchem_data, on='PubCHEM', how='left')

    merged = merged.dropna(subset=['ISOSMILES']).drop_duplicates().reset_index(drop=True)
    print(f"Chemical merge complete: {merged.shape[0]} rows")

    # Morgan fingerprints
    fps = []
    for smi in merged['ISOSMILES']:
        m = Chem.MolFromSmiles(smi)
        if m:
            bv = rdMolDescriptors.GetMorganFingerprintAsBitVect(m, 2, nBits=256)
            arr = np.zeros((256,), dtype=int)
            ConvertToNumpyArray(bv, arr)
            fps.append(arr)
        else:
            fps.append(np.zeros((256,), dtype=int))
    fps_df = pd.DataFrame(fps)

    out = pd.concat([merged.reset_index(drop=True), fps_df], axis=1)
    drop_cols = [c for c in ('PubCHEM','ISOSMILES') if c in out]
    out = out.drop(columns=drop_cols)

    out.to_csv(output_csv_path, index=False)
    print(f"Preprocessed CSV written to: {output_csv_path}")
    return out


def preprocess_new_data(
    prepped_csv: str,
    scaler
) -> Tuple[pd.DataFrame, np.ndarray]:
    """
    Aligns the 256+omics features to the scaler’s expected ordering,
    then returns both the raw DataFrame and the scaled numpy array.
    """
    df = pd.read_csv(prepped_csv).loc[:, ~pd.read_csv(prepped_csv, nrows=1).columns.duplicated()]
    X = df.drop(columns=['DRUG_ID','DRUG_NAME','COSMIC_ID','CCLE_Name','CANCER_TYPE'], errors='ignore')

    if hasattr(scaler, 'feature_names_in_'):
        # align columns
        exp = list(scaler.feature_names_in_)
        missing = set(exp) - set(X.columns)
        for m in missing:
            X[m] = 0
        X = X[exp]

    Xs = scaler.transform(X)
    return df, Xs

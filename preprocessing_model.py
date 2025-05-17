# File: preprocessing_model.py
# Refined preprocessing module for chemoresistance prediction pipeline

import os
import pandas as pd
import numpy as np
from typing import Tuple, Optional
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.DataStructs import ConvertToNumpyArray

# Directory paths
BASE_DIR = os.path.dirname(__file__)
PREPROCESS_DIR = os.path.join(BASE_DIR, 'Preprocess_files')


def align_features(
    X_new: pd.DataFrame,
    expected_features: list[str]
) -> pd.DataFrame:
    """
    Aligns columns of X_new to match expected_features:
    - Adds missing columns (filled with zeros).
    - Drops extra columns.
    - Reorders columns to expected order.
    """
    # Add missing features
    missing = set(expected_features) - set(X_new.columns)
    if missing:
        missing_df = pd.DataFrame(
            0.0,
            index=X_new.index,
            columns=list(missing)
        )
        X_new = pd.concat([X_new, missing_df], axis=1)

    # Drop extra features
    extra = set(X_new.columns) - set(expected_features)
    if extra:
        X_new = X_new.drop(columns=list(extra))

    # Reorder to match training
    X_new = X_new.reindex(columns=expected_features)
    return X_new


def preprocess_user_dataset(
    user_csv_path: str,
    output_csv_path: Optional[str] = None
) -> pd.DataFrame:
    """
    Loads a raw user CSV, validates required omics features,
    merges chemical information, generates Morgan fingerprints,
    and saves a preprocessed CSV for downstream prediction.

    Returns the preprocessed DataFrame.
    """
    # Determine output path
    if output_csv_path is None:
        output_csv_path = os.path.join(BASE_DIR, 'user_preprocessed_output.csv')

    # Load and validate
    user_df = pd.read_csv(user_csv_path)
    print(f"Loaded user dataset: {user_csv_path}")

    # Reference feature sets
    trans_ref = pd.read_csv(
        os.path.join(PREPROCESS_DIR, 'CCLE_Transcriptomics_cleaned.csv'),
        nrows=0
    ).columns
    prot_ref  = pd.read_csv(
        os.path.join(PREPROCESS_DIR, 'CCLE_Proteomics_cleaned.csv'),
        nrows=0
    ).columns

    metadata_cols = {
        'DRUG_ID', 'DRUG_NAME', 'CCLE_Name',
        'COSMIC_ID', 'CANCER_TYPE', 'LN_IC50'
    }
    features = set(user_df.columns) - metadata_cols

    missing_trans = set(trans_ref) - features
    missing_prot  = set(prot_ref) - features
    if missing_trans or missing_prot:
        msg = "Both transcriptomics and proteomics features are required."
        if missing_trans:
            msg += f" Missing transcriptomics: {list(missing_trans)[:5]}..."
        if missing_prot:
            msg += f" Missing proteomics: {list(missing_prot)[:5]}..."
        raise ValueError(msg)
    print("Required omics features present.")

    # Merge chemical information to get ISOSMILES
    drug_info    = pd.read_csv(os.path.join(PREPROCESS_DIR, 'cleaned_drug_info.csv'))
    pubchem_data = pd.read_csv(os.path.join(PREPROCESS_DIR, 'cleaned_pubchem.csv'))

    if 'ISOSMILES' in user_df.columns:
        merged_df = user_df.copy()
    elif 'PubCHEM' in user_df.columns:
        merged_df = pd.merge(user_df, pubchem_data, on='PubCHEM', how='left')
        merged_df = merged_df.dropna(subset=['ISOSMILES'])
    else:
        tmp = pd.merge(user_df, drug_info, on='DRUG_ID', how='left')
        tmp = tmp.dropna(subset=['PubCHEM'])
        merged_df = pd.merge(tmp, pubchem_data, on='PubCHEM', how='left')
        merged_df = merged_df.dropna(subset=['ISOSMILES'])

    merged_df = merged_df.drop_duplicates().reset_index(drop=True)
    print(f"Chemical merge complete ({merged_df.shape[0]} rows).")

    # Generate Morgan fingerprints
    print("Generating Morgan fingerprints...")
    fps = []
    for smi in merged_df['ISOSMILES']:
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=256)
            arr = np.zeros((256,), dtype=int)
            ConvertToNumpyArray(fp, arr)
            fps.append(arr)
        else:
            print(f"Invalid SMILES skipped: {smi}")
            fps.append(np.zeros((256,), dtype=int))

    fps_df = pd.DataFrame(fps)
    enriched = merged_df.join(fps_df)
    # Drop columns no longer needed
    drop_cols = [c for c in ['PubCHEM', 'ISOSMILES'] if c in enriched.columns]
    preprocessed = enriched.drop(columns=drop_cols)

    # Save out
    preprocessed.to_csv(output_csv_path, index=False)
    print(f"Preprocessed data saved to: {output_csv_path}")
    return preprocessed


def preprocess_new_data(
    file_path: str,
    scaler
) -> Tuple[pd.DataFrame, np.ndarray]:
    """
    Loads the user-preprocessed CSV, aligns its feature columns
    to the scaler's expected inputs, scales them, and returns
    both the raw DataFrame and scaled numpy array.
    """
    data = pd.read_csv(file_path)
    # remove duplicate columns if any
    data = data.loc[:, ~data.columns.duplicated()]

    # drop metadata columns before scaling
    X_raw = data.drop(
        columns=[ 'DRUG_ID', 'DRUG_NAME', 'COSMIC_ID', 'CCLE_Name', 'CANCER_TYPE' ],
        errors='ignore'
    )

    if hasattr(scaler, 'feature_names_in_'):
        expected = list(scaler.feature_names_in_)
        X_aligned = align_features(X_raw, expected)
    else:
        X_aligned = X_raw

    X_scaled = scaler.transform(X_aligned)
    return data, X_scaled

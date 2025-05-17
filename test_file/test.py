import pandas as pd

# Load the dataset
df_test = pd.read_csv("random_rows_for_website_cleaned.csv")

# Check if the column exists
if "TSPAN6 (7105)" in df_test.columns:
    # Duplicate the column manually by creating a new DataFrame with duplicated column
    df_duplicate = df_test.copy()
    df_duplicate = pd.concat([df_duplicate, df_test[["TSPAN6 (7105)"]]], axis=1)
    
    # Rename the last column back to the same name to force a duplicate
    df_duplicate.columns.values[-1] = "TSPAN6 (7105)"

    print("⚠️ Successfully added duplicate 'TSPAN6 (7105)' column.")
else:
    print("❌ Column 'TSPAN6 (7105)' not found.")

# Save the DataFrame (even with duplicate columns)
df_duplicate.to_csv("user_input_dataset_with_duplicate_transcript.csv", index=False)

print("✅ File saved with duplicate column.")

#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np

def main():
    # Parent directory containing replicate subfolders
    parent_dir = os.getcwd()

    # Identify and sort replicate directories
    replicates = sorted([
        d for d in os.listdir(parent_dir)
        if os.path.isdir(os.path.join(parent_dir, d))
    ])

    # Read density_data.csv from each replicate
    dfs = []
    for rep in replicates:
        csv_path = os.path.join(parent_dir, rep, 'density_data.csv')
        if os.path.isfile(csv_path):
            dfs.append(pd.read_csv(csv_path))
        else:
            print(f"Warning: '{csv_path}' not found; skipping replicate '{rep}'.")

    if not dfs:
        print("No density_data.csv files found in any replicate folders.")
        return

    # Verify all DataFrames share the same shape and columns
    ref_cols = dfs[0].columns
    ref_rows = len(dfs[0])
    for idx, df in enumerate(dfs[1:], start=1):
        if not df.columns.equals(ref_cols):
            raise ValueError(
                f"Column mismatch in replicate '{replicates[idx]}'."
            )
        if len(df) != ref_rows:
            raise ValueError(
                f"Row count mismatch in replicate '{replicates[idx]}'."
            )

    # Stack data and compute per-row means across replicates
    stacked = np.stack([df.values for df in dfs], axis=2)
    row_means = np.mean(stacked, axis=2)

    # Build summary DataFrame
    summary_df = pd.DataFrame(row_means, columns=ref_cols)

    # Remove the third column (index 2)
    if summary_df.shape[1] >= 3:
        summary_df.drop(summary_df.columns[2], axis=1, inplace=True)

    # Save to CSV
    summary_df.to_csv('density_summary.csv', index=False)
    print("Created 'density_summary.csv' with per-row means across replicates (excluding third column).")

if __name__ == '__main__':
    main()

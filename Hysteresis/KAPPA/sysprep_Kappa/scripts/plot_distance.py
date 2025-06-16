#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('agg')   # for non‑interactive backends
import matplotlib.pyplot as plt
import sys
import os

def main(csv_path='density_summary.csv'):
    if not os.path.exists(csv_path):
        print(f"ERROR: '{csv_path}' not found.")
        sys.exit(1)

    # --- 1) Load & coerce ---
    df = pd.read_csv(csv_path, header=None, usecols=[0,1,2],
                     names=['ions','z','norm_density'])
    # force numeric, drop failures
    df['ions']         = pd.to_numeric(df['ions'],         errors='coerce').astype('Int64')
    df['z']            = pd.to_numeric(df['z'],            errors='coerce')
    df['norm_density'] = pd.to_numeric(df['norm_density'], errors='coerce')
    df = df.dropna(subset=['ions','z','norm_density'])
    df['ions'] = df['ions'].astype(int)  # now safe

    # clamp alpha to [0,1]
    df['norm_density'] = df['norm_density'].clip(0.0, 1.0)

    if df.empty:
        print("ERROR: No valid data found in CSV.")
        sys.exit(1)

    # --- 2) Set up figure ---
    fig, ax = plt.subplots(
        figsize=(3.0 * 8.25 / 2.54,   # reproduce your original ~3" height
                 3.0 * 8.25 / 2.54)
    )
    plt.tick_params(axis="both", labelsize=8)

    ax.set_xlabel('Number of ions', fontsize=9)
    ax.set_ylabel('Distance / Å',     fontsize=9)

    # x‑ticks every 5 (or adjust to your spacing)
    xticks = np.arange(df['ions'].min(), df['ions'].max()+1, 5)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks, rotation=90, fontsize=8)

    # grid only on x
    ax.grid(axis="x", color='lightgray', linewidth=0.5)

    # limits
    ax.set_xlim(df['ions'].min() - 1, df['ions'].max() + 1)
    ax.set_ylim(df['z'].min()-1, df['z'].max())

    # --- 3) Scatter plot ---
    # one call: alpha accepts an array
    ax.scatter(
        df['ions'], df['z'],
        s=40,
        c='blue',
        alpha=df['norm_density'],
        marker='s',
        linewidths=0
    )

    # --- 4) Save ---
    out_png = 'distance_graph.png'
    plt.savefig(out_png, format='png', dpi=500)
    print(f"Plot saved to '{out_png}'")

if __name__ == '__main__':
    # optionally take CSV path as first arg
    if len(sys.argv) > 1:
        main(sys.argv[1])
    else:
        main()


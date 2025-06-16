#!/usr/bin/env python3
import csv
import glob
import matplotlib
import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
import os
import sys

# Use non-interactive backend
matplotlib.use('agg')

# read the arguments
if len(sys.argv) != 4:
    print("Usage: script.py <atom_name> <sign> <surface_area>")
    sys.exit(1)

atom = sys.argv[1]
sign = int(sys.argv[2])
area = int(float(sys.argv[3]))
print(f"Processing atom: {atom}, sign: {sign}, area: {area}")

# find all potential files
files = glob.glob('*/potential*xvg')
print(f"Found {len(files)} potential files: {files}")

if not files:
    print("ERROR: No potential files found. Check directory structure and file pattern.")
    sys.exit(1)

# get number of ions in the model from the file name
particle_numbers = []
for fname in files:
    try:
        part_num = fname.split('_')[1].replace('potential','').replace('.xvg','')
        particle_numbers.append(int(part_num))
        print(f"Processed file {fname}: particle number = {part_num}")
    except Exception as e:
        print(f"Error processing file {fname}: {e}")

if not particle_numbers:
    print("ERROR: Could not extract any particle numbers from filenames.")
    sys.exit(1)

path = os.getcwd()
print(f"Current working directory: {path}")

# set up the plot
plt.figure(figsize=(3.0 * 8.25 / 2.54, 3.0 * 8.25 / 2.54))
plt.tick_params(axis="both", labelsize=8)
plt.xlabel('Number of ions', fontsize=9)
plt.ylabel('Distance / Å', fontsize=9)
plt.xticks(np.arange(0, max(particle_numbers)+1, 5.0), rotation=90)
plt.grid(axis="x", color='lightgray', linewidth=0.5)
plt.xlim(min(particle_numbers)-1, max(particle_numbers)+1)
plt.ylim(2, 10)

dens = []
valid_data = False

# process each directory
for i in particle_numbers:
    try:
        target_dir = os.path.join(path, str(i))
        print(f"Changing to directory: {target_dir}")
        os.chdir(target_dir)

        if not os.path.exists('NVT.gro') or not os.path.exists('NVT.xtc'):
            print(f"ERROR: Required files missing in {target_dir}")
            dens.append(np.zeros(200))
            os.chdir(path)
            continue

        u = mda.Universe('NVT.gro', 'NVT.xtc')
        coords = []
        for frame in u.trajectory:
            center = u.select_atoms(f"name {atom}")
            if len(center) == 0:
                continue
            coords.append(center.positions[:, 2])

        os.chdir(path)

        if not coords:
            print(f"WARNING: No coordinate data collected for particle number {i}")
            dens.append(np.zeros(200))
            continue

        flat = np.concatenate(coords)
        hist, bins = np.histogram(flat, bins=200, density=False)
        axis = (bins[:-1] + bins[1:]) / 2.0
        # shift so zero is surface at 10 Å
        axis = axis - 10

        if np.any(hist > 0):
            valid_data = True
        else:
            print(f"WARNING: Histogram for {i} is all zeros")

        dens.append(hist)

    except Exception as e:
        print(f"ERROR processing particle {i}: {e}")
        dens.append(np.zeros(200))
        os.chdir(path)

if not valid_data:
    print("ERROR: No valid data found in any directory.")
    plt.text(0.5, 0.5, "NO DATA FOUND", 
             ha='center', va='center', transform=plt.gca().transAxes, fontsize=14)
    plt.savefig('distance_graph_no_data.png', bbox_inches='tight')
    sys.exit(1)

# normalize
max_vals = [np.max(d) for d in dens if np.any(d > 0)]
norm = max(max_vals) if max_vals else 1.0
print(f"Using normalization factor: {norm}")

# --- CSV output ---
csv_output_path = os.path.join(path, 'density_data.csv')
with open(csv_output_path, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['particle_number', 'z_position', 'density', 'normalized_density'])
    for idx, pn in enumerate(particle_numbers):
        hist = dens[idx]
        if hist.size == 0:
            continue
        z_positions = axis  # same axis used for plotting
        normed = hist / norm if norm else hist
        for z, raw, normval in zip(z_positions, hist, normed):
            writer.writerow([pn, z, raw, normval])

print(f"Wrote CSV → {csv_output_path}")
# -------------------

# plotting
for idx, pn in enumerate(particle_numbers):
    h = dens[idx]
    if not np.any(h > 0):
        print(f"Skipping plot for {pn} (all zeros)")
        continue
    alpha = np.clip(h / norm, 0.0, 1.0)
    plt.scatter([pn]*len(axis), axis, s=40, c='blue', alpha=alpha, marker='s', linewidths=0)

plt.savefig('distance_graph.png', bbox_inches='tight')
print("Saved plot → distance_graph.png")


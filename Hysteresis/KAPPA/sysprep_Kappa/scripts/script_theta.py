#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script analyzes subfolders whose names are the total number of ions.
It reads the file density_<ion-number>.xvg from each subfolder and finds the highest density peak,
integrates the area from (highest peak - 0.2) to (highest peak + 0.2), then converts that integrated
density into the contact layer charge density by multiplying it with the elementary charge (≈1.6e-19 C)
and 10^20, yielding units of µC·cm⁻².
For the x-axis, the total number of ions (obtained from each folder name) is converted into the surface
charge density by dividing it by the old conversion factor and then multiplying it with the elementary
charge (≈1.6e-19 C) and 10^20, so that the x-axis represents the surface charge density in µC·cm⁻².

Additionally, the script performs an incremental zero-intercept linear regression test:
  - It starts with the first 3 points.
  - It then adds one additional point at a time; for each candidate group, it computes Pearson’s
    correlation coefficient (with means subtracted) and squares it to obtain R².
  - If adding the new point yields R² > 0.99, the candidate point is added to the LG group; otherwise,
    the process stops (excluding that candidate point and any subsequent ones).
  - A dashed grey line is drawn for the LG model extended over the entire x-axis.
  - For the remaining points, a horizontal dashed grey line is drawn at the y value of the first point excluded
    from the LG model.
  - The intersection of these two lines is computed (as excluded_y / m_candidate), and its rounded
    value is displayed in a red label ("Mono ≈ i µC·cm⁻²") that is placed in the upper right corner of the graph.

The x-axis limits are set dynamically from (min(surface charge density) – 10, but not below 0) to (max(surface charge density) + 10),
and similarly for the y-axis.
"""

import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import glob
import sys

# Conversion factor used in the original conversion from integrated density to number of ions.
CONVERSION_FACTOR = 4.262577 * 3.9376
# Elementary charge (in coulombs)
ELEMENTARY_CHARGE = 1.60217663e-19

def read_density_file(filepath):
    """
    Reads a density_xvg file and returns two numpy arrays: x and density.
    Lines starting with '#' or '@', or empty lines, are skipped.
    """
    x_vals = []
    y_vals = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or line.startswith('@') or not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                x_vals.append(float(parts[0]))
                y_vals.append(float(parts[1]))
            except ValueError:
                continue
    return np.array(x_vals), np.array(y_vals)

def integrate_peak(x, y, window=0.2):
    """
    Finds the highest density peak and integrates the density from 
    (x_peak - window) to (x_peak + window).
    """
    max_idx = np.argmax(y)
    x_peak = x[max_idx]
    left_bound = x_peak - window
    right_bound = x_peak + window
    integration_mask = (x >= left_bound) & (x <= right_bound)
    x_integ = x[integration_mask]
    y_integ = y[integration_mask]
    if len(x_integ) < 2:
        return 0.0
    area = np.trapezoid(y_integ, x_integ)
    return area

def main():
    base_dir = os.getcwd()
    subfolders = [d for d in os.listdir(base_dir) if os.path.isdir(d) and d.isdigit()]
    subfolders = sorted(subfolders, key=lambda s: int(s))
    
    total_ions = []  # Total ions, taken from folder names.
    contact_layer_charge_density_list = []  # Converted, non-integer contact layer charge density in µC/cm².
    
    for folder in subfolders:
        file_name = f"density_{folder}.xvg"
        file_path = os.path.join(base_dir, folder, file_name)
        if not os.path.isfile(file_path):
            print(f"Warning: {file_path} not found. Skipping folder {folder}.")
            continue
        x, y = read_density_file(file_path)
        if x.size == 0 or y.size == 0:
            print(f"Warning: No valid data in {file_path}.")
            continue
        # Use a window of 0.2 surrounding the highest peak.
        area = integrate_peak(x, y, window=0.2)
        # Convert the integrated density into contact layer charge density (µC/cm²):
        cl_charge_density = area * ELEMENTARY_CHARGE * 1e20
        print(f"Folder {folder}: integrated area = {area}, converted to {cl_charge_density:.4g} µC/cm² in contact layer.")
        total_ions.append(int(folder))
        contact_layer_charge_density_list.append(cl_charge_density)
    
    # Convert total ions (folder names) into surface charge density (µC/cm²) by dividing by the conversion factor,
    # then multiplying by the elementary charge and 1e20.
    total_ions_array = np.array(total_ions)
    surface_charge_density = (total_ions_array / CONVERSION_FACTOR) * (ELEMENTARY_CHARGE * 1e20)
    contact_layer_charge_density = np.array(contact_layer_charge_density_list)
    
    plt.figure(figsize=(8, 6))
    plt.scatter(surface_charge_density, contact_layer_charge_density, marker='o', color='blue')
    plt.xlabel("|Surface charge density| [µC·cm⁻²]", fontsize=12)
    plt.ylabel("|Contact layer charge density| [µC·cm⁻²]", fontsize=12)
    plt.title("Contact Layer vs Surface Charge Density Value", fontsize=14)
    plt.grid(True)
    
    # Dynamic axis limits.
    xmin = max(surface_charge_density.min() - 10, 0)
    xmax = surface_charge_density.max() + 10
    plt.xlim(xmin, xmax)
    ymin = max(contact_layer_charge_density.min() - 10, 0)
    ymax = contact_layer_charge_density.max() + 10
    plt.ylim(ymin, ymax)
    
    # Incremental zero-intercept linear regression (LG model) based on Pearson's r².
    n_points = len(surface_charge_density)
    excluded_y = None  # This will store the y value of the first excluded candidate.
    if n_points >= 3:
        candidate_indices = [0, 1, 2]
        # Try adding one point at a time.
        for i in range(3, n_points):
            temp_indices = candidate_indices + [i]
            x_temp = surface_charge_density[temp_indices]
            y_temp = contact_layer_charge_density[temp_indices]
            mean_x = np.mean(x_temp)
            mean_y = np.mean(y_temp)
            r_num = np.sum((x_temp - mean_x) * (y_temp - mean_y))
            r_den = np.sqrt(np.sum((x_temp - mean_x)**2) * np.sum((y_temp - mean_y)**2))
            r_temp = r_num / r_den if r_den != 0 else 0.0
            r2_candidate = r_temp**2
            print(f"Testing candidate indices {temp_indices}: r² = {r2_candidate:.4f}")
            if r2_candidate > 0.999:
                candidate_indices.append(i)
            else:
                # Save the y-value of the first excluded point and exit the loop.
                excluded_y = contact_layer_charge_density[i]
                break
        
        candidate_x = surface_charge_density[candidate_indices]
        candidate_y = contact_layer_charge_density[candidate_indices]
        m_candidate = np.sum(candidate_x * candidate_y) / np.sum(candidate_x**2)
        # Draw the LG regression line extended over the entire x-axis.
        x_line = np.linspace(xmin, xmax, 100)
        y_line = m_candidate * x_line
        plt.plot(x_line, y_line, color='grey', linestyle='dashed')
        
        # Draw the horizontal dashed line at the excluded point's y value.
        if excluded_y is not None:
            plt.hlines(excluded_y, xmin, xmax, colors="grey", linestyles="dashed")
            # Compute the intersection x-coordinate of the LG line and the horizontal line.
            if m_candidate != 0:
                x_int = excluded_y / m_candidate
                rounded_val = round(x_int,1)
                label_text = f"Mono ≈ {rounded_val} µC·cm⁻²"
                # Place the label in the upper right corner.
                plt.text(xmax - 5, ymax - 5, label_text, horizontalalignment='right', 
                         verticalalignment='top', fontsize=12, color='red')
        else:
            # If all points were included, no excluded point exists to display.
            print("All points were included in the LG model; no excluded candidate to draw a horizontal line.")
    
    plt.tight_layout()
    plt.savefig("contact_layer_density.png", format="png", bbox_inches='tight', dpi = 300)
    plt.close()
    print("Plot saved as contact_layer_density.png")

if __name__ == '__main__':
    main()


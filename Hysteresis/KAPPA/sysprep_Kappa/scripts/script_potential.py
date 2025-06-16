#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script performs two analyses:

1. Contact Layer Analysis:
   - Analyzes subfolders whose names are the total number of ions
   - Reads density_<ion-number>.xvg from each subfolder and finds the highest density peak
   - Integrates the area around the peak and converts to contact layer charge density
   - Performs incremental zero-intercept linear regression to find the monolayer transition point

2. Surface Charge vs Potential Analysis:
   - Uses the monolayer transition point (x_int) from the first analysis
   - Analyzes potential_<ion-number>.xvg files from the same subfolders
   - Creates a plot of surface charge density vs potential
   - Performs linear regression to predict potential at the monolayer transition point
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import glob
import csv
from matplotlib.ticker import FuncFormatter
from matplotlib.ticker import MultipleLocator

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

def analyze_contact_layer():
    """
    Performs the contact layer analysis and returns the monolayer transition point (x_int).
    Shows actual signs of charge densities on the plot.
    """
    base_dir = os.getcwd()
    subfolders = [d for d in os.listdir(base_dir) if os.path.isdir(d) and d.isdigit()]
    subfolders = sorted(subfolders, key=lambda s: int(s))
    
    # First, gather all necessary data from all folders
    folders_data = []
    for folder in subfolders:
        folder_data = {}
        folder_data['folder'] = folder
        folder_data['ions'] = int(folder)
        
        # Get density file data
        density_file = os.path.join(base_dir, folder, f"density_{folder}.xvg")
        if not os.path.isfile(density_file):
            print(f"Warning: {density_file} not found. Skipping folder {folder}.")
            continue
            
        # Get potential file data to determine charge sign
        potential_file = os.path.join(base_dir, folder, f"potential_{folder}.xvg")
        if not os.path.isfile(potential_file):
            print(f"Warning: {potential_file} not found. Skipping folder {folder}.")
            continue
            
        # Extract potential value from the last line
        last_line = ""
        with open(potential_file, 'r') as f:
            for line in f:
                if not line.startswith(('#', '@')) and line.strip():
                    last_line = line
        
        if not last_line:
            print(f"Warning: No valid data found in {potential_file}")
            continue
            
        parts = last_line.split()
        if len(parts) >= 2:
            potential_value = float(parts[1])
            folder_data['potential'] = potential_value
        else:
            print(f"Warning: Invalid format in last line of {potential_file}")
            continue
        
        # Process density data
        x, y = read_density_file(density_file)
        if x.size == 0 or y.size == 0:
            print(f"Warning: No valid data in {density_file}.")
            continue
            
        # Use a window of 0.2 surrounding the highest peak
        area = integrate_peak(x, y, window=0.2)
        
        # Convert the integrated density into contact layer charge density (µC/cm²)
        cl_charge_density = area * ELEMENTARY_CHARGE * 1e20
        
        # Determine signs based on potential polarity
        # If potential is positive, surface charge is negative, contact layer charge is positive
        # If potential is negative, surface charge is positive, contact layer charge is negative
        surface_charge = (folder_data['ions'] / CONVERSION_FACTOR) * (ELEMENTARY_CHARGE * 1e20)
        
        if potential_value > 0:
            folder_data['surface_charge'] = -abs(surface_charge)
            folder_data['contact_layer_charge'] = abs(cl_charge_density)
        else:
            folder_data['surface_charge'] = abs(surface_charge)
            folder_data['contact_layer_charge'] = -abs(cl_charge_density)
            
        print(f"Folder {folder}: potential = {potential_value:.3f} V, " 
              f"surface charge = {folder_data['surface_charge']:.4g} µC/cm², "
              f"contact layer charge = {folder_data['contact_layer_charge']:.4g} µC/cm²")
              
        folders_data.append(folder_data)
    
    # Convert to arrays for processing
    if not folders_data:
        print("No valid data found in any folder")
        return None
        
    # Extract surface charge and contact layer charge arrays with proper signs
    surface_charge_density = np.array([d['surface_charge'] for d in folders_data])
    contact_layer_charge_density = np.array([d['contact_layer_charge'] for d in folders_data])
    
    # For plotting in first quadrant, take absolute values but remember the signs
    abs_surface_charge = np.abs(surface_charge_density)
    abs_contact_layer = np.abs(contact_layer_charge_density)
    
    # This will store the monolayer transition point (with sign)
    x_int = None
    
    # Create figure with absolute values for first quadrant plotting
    plt.figure(figsize=(8, 6))
    plt.scatter(abs_surface_charge, abs_contact_layer, marker='o', color='blue', zorder=3)
    plt.xlabel("Surface Charge Density [µC·cm⁻²]", fontsize=12)
    plt.ylabel("Contact Layer Charge Density [µC·cm⁻²]", fontsize=12)
    plt.title("Contact Layer vs. Surface Charge Density", fontsize=14)
    #plt.grid(True, color='grey', alpha=0.5, linestyle='-')  # Continuous grey lines with same style as second plot
    
    # Dynamic axis limits with rounding to nearest multiple of 5, never below 0
    if len(abs_surface_charge) > 0:
        xmin_raw = abs_surface_charge.min() - 10
        xmax_raw = abs_surface_charge.max() + 10
        # Round xmin down to nearest multiple of 5, but never below 0
        xmin = max(5 * np.floor(xmin_raw / 5), 0)
        # Round xmax up to nearest multiple of 5
        xmax = 5 * np.ceil(xmax_raw / 5)
    else:
        xmin, xmax = 0, 100
    plt.xlim(xmin, xmax)
    
    if len(abs_contact_layer) > 0:
        ymin_raw = abs_contact_layer.min() - 10
        ymax_raw = abs_contact_layer.max() + 10
        # Round ymin down to nearest multiple of 5, but never below 0
        ymin = max(5 * np.floor(ymin_raw / 5), 0)
        # Round ymax up to nearest multiple of 5
        ymax = 5 * np.ceil(ymax_raw / 5)
    else:
        ymin, ymax = 0, 100
    plt.ylim(ymin, ymax)
    
    # Use FuncFormatter to show minus signs on the axis ticks and format to 0 decimal places
    ax = plt.gca()
      #Set grid lines interval
    #ax.xaxis.set_major_locator(MultipleLocator(10))
    #ax.yaxis.set_major_locator(MultipleLocator(10))
    # For x-axis, determine if we need to show negative signs
    if np.any(surface_charge_density < 0):
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"-{int(x)}" if np.any(surface_charge_density < 0) and x > 0 else f"{int(x)}"))
    else:
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{int(x)}"))
    
    # For y-axis, determine if we need to show negative signs
    if np.any(contact_layer_charge_density < 0):
        ax.yaxis.set_major_formatter(FuncFormatter(lambda y, pos: f"-{int(y)}" if np.any(contact_layer_charge_density < 0) and y > 0 else f"{int(y)}"))
    else:
        ax.yaxis.set_major_formatter(FuncFormatter(lambda y, pos: f"{int(y)}"))
    
    # Incremental zero-intercept linear regression (LG model) based on Pearson's r²
    n_points = len(abs_surface_charge)
    excluded_y_idx = None  # This will store the index of the first excluded candidate
    
    if n_points >= 3:
        candidate_indices = [0, 1, 2]
        # Try adding one point at a time
        for i in range(3, n_points):
            temp_indices = candidate_indices + [i]
            x_temp = abs_surface_charge[temp_indices]
            y_temp = abs_contact_layer[temp_indices]
            mean_x = np.mean(x_temp)
            mean_y = np.mean(y_temp)
            r_num = np.sum((x_temp - mean_x) * (y_temp - mean_y))
            r_den = np.sqrt(np.sum((x_temp - mean_x)**2) * np.sum((y_temp - mean_y)**2))
            r_temp = r_num / r_den if r_den != 0 else 0.0
            r2_candidate = r_temp**2
            print(f"Testing candidate indices {temp_indices}: r² = {r2_candidate:.4f}")
            if r2_candidate > 0.9995:
                candidate_indices.append(i)
            else:
                # Save the index of the first excluded point and exit the loop
                excluded_y_idx = i
                break
        
        candidate_x = abs_surface_charge[candidate_indices]
        candidate_y = abs_contact_layer[candidate_indices]
        
        # Calculate slope using the absolute values
        m_candidate = np.sum(candidate_x * candidate_y) / np.sum(candidate_x**2)
        
        # Draw the LG regression line extended
        x_line = np.linspace(xmin, xmax, 100)
        y_line = m_candidate * x_line
        plt.plot(x_line, y_line, color='grey', linestyle='dashed')
        
        # If we have an excluded point, handle the monolayer transition calculation
        if excluded_y_idx is not None:
            excluded_y_abs = abs_contact_layer[excluded_y_idx]
            plt.hlines(excluded_y_abs, xmin, xmax, colors="grey", linestyles="dashed")
            
            # Compute the intersection x-coordinate of the LG line and the horizontal line
            if m_candidate != 0:
                x_int_abs = excluded_y_abs / m_candidate
                
                # The sign of x_int should match the sign of the surface charge density
                x_int_sign = np.sign(surface_charge_density[0])
                x_int = x_int_abs * x_int_sign
                
                rounded_val = round(x_int_abs, 1)
                # Format the label with sign, but plot in first quadrant
                if x_int_sign < 0:
                    label_text = f"$\\sigma_M = -{rounded_val}$ µC·cm⁻²"
                else:
                    label_text = f"$\\sigma_M = {rounded_val}$ µC·cm⁻²"
                
                # Place the label at [xmax-3, ymax-3]
                plt.annotate(label_text, xy=(xmax-1, ymax-3), 
                         bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8),
                         horizontalalignment='right',
                         verticalalignment='top')
                
                # Print both the absolute and signed values
                print(f"Monolayer transition point: |x_int| = {x_int_abs:.4f}, x_int = {x_int:.4f} µC/cm²")
                
                # No red point on the first plot as requested
            else:
                print("Warning: m_candidate is zero, cannot compute x_int")
        else:
            # If all points were included, no excluded point exists
            print("All points were included in the LG model; no excluded candidate to draw a horizontal line.")
    
    # Save the data used for the plot to theta_sigma.txt - sorted by magnitude
    # Create sorted indices based on absolute values of surface charge density
    sorted_indices = np.argsort(np.abs(surface_charge_density))
    
    #with open("theta_sigma.txt", "w") as f:
     #   f.write("# Surface Charge Density [µC·cm⁻²]\tContact Layer Charge Density [µC·cm⁻²]\n")
      #  for idx in sorted_indices:
       #     f.write(f"{surface_charge_density[idx]:.7f}\t{contact_layer_charge_density[idx]:.7f}\n")
    #print("Data saved to theta_sigma.txt (sorted by magnitude of surface charge density)")
    
        # ——— write out CSV instead of TXT ———
    with open("theta_sigma.csv", "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([
            "surface_charge_density_µC/cm2",
            "contact_layer_charge_density_µC/cm2"
        ])
        for idx in sorted_indices:
            writer.writerow([
                f"{surface_charge_density[idx]:.7f}",
                f"{contact_layer_charge_density[idx]:.7f}"
            ])
    print("Data saved to theta_sigma.csv (sorted by magnitude)")




    #plt.tight_layout()
    plt.savefig("contact_layer_density.png", format="png", dpi=300)
    #bbox_inches='tight'Add to the plt.savefig if needed
    plt.close()
    print("Plot saved as contact_layer_density.png")
    
    return x_int

def analyze_surface_charge_potential(x_int):
    """
    Performs the surface charge vs potential analysis using the x_int value from
    the contact layer analysis. Plots in the first quadrant while showing the sign in labels.
    """
    if x_int is None:
        print("Warning: x_int not available. Using default value of 0.0")
        x_int = 0.0
    
    # Store original x_int (with sign) for later use in prediction
    x_int_original = x_int
    print(f"Original x_int with sign: {x_int_original}")
    
    # For threshold comparison, we need the absolute value
    x_int_abs = abs(x_int)
    threshold = x_int_abs + 0.1
    print(f"Using threshold |sigma| < {threshold}")

    # Dictionary to hold folder data
    folders_data = []

    # Loop over subfolders in the current directory that are named as numbers (number of ions)
    for folder in glob.glob("*"):
        if os.path.isdir(folder) and folder.isnumeric():
            folder_data = {}
            folder_data['folder'] = folder
            folder_data['ions'] = int(folder)

            # Construct file paths for density and potential data
            density_file = os.path.join(folder, f"density_{folder}.xvg")
            potential_file = os.path.join(folder, f"potential_{folder}.xvg")
            
            # Check if necessary files exist
            if not os.path.isfile(density_file) or not os.path.isfile(potential_file):
                print(f"Missing required files in folder {folder}. Skipping.")
                continue
            
            # Extract the potential value from the last line of the potential file
            last_line = ""
            with open(potential_file, 'r') as f:
                for line in f:
                    if not line.startswith(('#', '@')) and line.strip():
                        last_line = line
            
            if not last_line:
                print(f"Warning: No valid data found in {potential_file}")
                continue
                
            # Parse the last line to get the potential value
            parts = last_line.split()
            if len(parts) >= 2:
                potential_value = float(parts[1])
                folder_data['potential'] = potential_value
            else:
                print(f"Warning: Invalid format in last line of {potential_file}")
                continue
            
            # Calculate the surface charge density (sigma) in µC·cm⁻²
            sigma = (folder_data['ions'] / CONVERSION_FACTOR) * ELEMENTARY_CHARGE * 1e20

            # Apply the sign based on potential polarity
            if potential_value > 0:
                sigma = -abs(sigma)
            else:
                sigma = abs(sigma)
                
            folder_data['sigma'] = sigma
            
            # Store the original values (with signs)
            folder_data['original_potential'] = potential_value
            folder_data['original_sigma'] = sigma
            
            # For plotting in first quadrant:
            # Take absolute values but remember the signs
            folder_data['abs_potential'] = abs(potential_value)
            folder_data['abs_sigma'] = abs(sigma)
            
            print(f"Folder {folder}: potential = {potential_value:.3f} V, "
                  f"sigma = {sigma:.3f} µC/cm²")
                  
            folders_data.append(folder_data)

    # Extract plotting data
    if not folders_data:
        print("No valid data found in any folder")
        return
        
    # For plotting in first quadrant with absolute values
    abs_potential = np.array([d['abs_potential'] for d in folders_data])
    abs_sigma = np.array([d['abs_sigma'] for d in folders_data])
    
    # Keep original values with signs for calculations
    original_potential = np.array([d['original_potential'] for d in folders_data])
    original_sigma = np.array([d['original_sigma'] for d in folders_data])
    
    # Save the data used for the plot to sigma_potential.txt - sorted by magnitude
    # Create sorted indices based on absolute values of surface charge density
    sorted_indices = np.argsort(np.abs(original_sigma))
    
    #with open("sigma_potential.txt", "w") as f:
     #   f.write("# Surface Charge Density [µC·cm⁻²]\tPotential [V]\n")
      #  for idx in sorted_indices:
       #     f.write(f"{original_sigma[idx]:.7f}\t{original_potential[idx]:.7f}\n")
    #print("Data saved to sigma_potential.txt (sorted by magnitude of surface charge density)")
    
    
        # ——— write out CSV instead of TXT ———
    with open("sigma_potential.csv", "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([
            "surface_charge_density_µC/cm2",
            "potential_V"
        ])
        for idx in sorted_indices:
            writer.writerow([
                f"{original_sigma[idx]:.7f}",
                f"{original_potential[idx]:.7f}"
            ])
    print("Data saved to sigma_potential.csv (sorted by magnitude)")


    
    # Create the scatter plot with absolute values (first quadrant)
    plt.figure(figsize=(8, 6))
    plt.scatter(abs_potential, abs_sigma, marker='o', color='blue', zorder=3)
    plt.xlabel("Potential [V]", fontsize=12)
    plt.ylabel("Surface Charge Density [µC·cm⁻²]", fontsize=12)
    plt.title("Surface Charge Density vs. Potential", fontsize=14)
    
    # Dynamic axis limits with rounding to nearest multiple of 5, never below 0
    if abs_potential.size:
        xmin_raw = abs_potential.min() - 5
        xmax_raw = abs_potential.max() + 5
        # Round xmin down to nearest multiple of 5, but never below 0
        xmin = max(5 * np.floor(xmin_raw / 5), 0)
        # Round xmax up to nearest multiple of 5
        xmax = 5 * np.ceil(xmax_raw / 5)
    else:
        xmin, xmax = 0, 5
    plt.xlim(xmin, xmax)
    
    if abs_sigma.size:
        ymin_raw = abs_sigma.min() - 10
        ymax_raw = abs_sigma.max() + 10
        # Round ymin down to nearest multiple of 5, but never below 0
        ymin = max(5 * np.floor(ymin_raw / 5), 0)
        # Round ymax up to nearest multiple of 5
        ymax = 5 * np.ceil(ymax_raw / 5)
    else:
        ymin, ymax = 0, 10
    plt.ylim(ymin, ymax)
    
    # Add grid with continuous grey lines
    #plt.grid(True, color='grey', alpha=0.5, linestyle='-')
    
    # Use FuncFormatter to show minus signs on the axis ticks with 0 decimal places
    ax = plt.gca()
    # For x-axis, determine if we need to show negative signs
    if np.any(original_potential < 0):
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"-{int(x)}" if np.any(original_potential < 0) and x > 0 else f"{int(x)}"))
    else:
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{int(x)}"))
    
    # For y-axis, determine if we need to show negative signs
    if np.any(original_sigma < 0):
        ax.yaxis.set_major_formatter(FuncFormatter(lambda y, pos: f"-{int(y)}" if np.any(original_sigma < 0) and y > 0 else f"{int(y)}"))
    else:
        ax.yaxis.set_major_formatter(FuncFormatter(lambda y, pos: f"{int(y)}"))
    
    # Add a horizontal dashed line at y = sigma_M (absolute value for first quadrant plotting)
    if x_int is not None:
        plt.hlines(abs(x_int), xmin, xmax, colors="grey", linestyles="dashed", zorder=2)
        print(f"Added horizontal dashed line at y = |sigma_M| = {abs(x_int):.4f} µC/cm²")
    
    # Variable to store the predicted potential for writing to file
    predicted_potential = None
    
    # Perform linear regression on points with |sigma| values around |x_int|
    if original_potential.size > 0 and original_sigma.size > 0:
        # Create a mask for points with absolute values below the threshold
        mask = abs_sigma < threshold
        
        print(f"Number of points matching criteria: {np.sum(mask)}")
        
        # Only perform regression if we have enough points
        if np.sum(mask) >= 2:  # Need at least 2 points for regression
            # We need to use the original values with signs for accurate regression
            x_regression = original_potential[mask]
            y_regression = original_sigma[mask]
            
            # Calculate slope and intercept using NumPy's polyfit
            slope, intercept = np.polyfit(x_regression, y_regression, 1)
            
            # Find where the regression line equals x_int (with sign)
            # y = mx + b, solving for x when y = x_int: x = (x_int - b) / m
            if slope != 0:
                predicted_potential = (x_int_original - intercept) / slope
            else:
                predicted_potential = 0
                print("Warning: Slope is zero, cannot predict potential")
            
            # Print information about the regression
            print(f"Linear regression: y = {slope:.4f}x + {intercept:.4f}")
            print(f"Predicted potential at x_int = {x_int_original}: {predicted_potential:.4f} V")
            
            # Add annotation with the predicted potential using consistent format, but showing sign
            rounded_val = round(abs(predicted_potential), 2)
            if predicted_potential < 0:
                legend_text = f"$\\phi_M = -{rounded_val}$ V"
            else:
                legend_text = f"$\\phi_M = {rounded_val}$ V"
                
            # Position the label at [xmax-3, ymax-3]
            plt.annotate(legend_text, xy=(xmax-1, ymax-3), 
                       bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8),
                       horizontalalignment='right',
                       verticalalignment='top')
            
            # Write the predicted potential to a text file with 5 decimal places
            with open("predicted_potential.txt", "w") as f:
                f.write(f"$\\phi_M$ = {predicted_potential:.5f} V\n")
                print("Predicted potential written to predicted_potential.txt with 5 decimal places")
                
        else:
            print(f"Not enough points matching criteria for regression")
    else:
        print("No data points available for regression")

    # Save the plot to the file 'surfacecharge_potential.png'
    plt.savefig("surfacecharge_potential.png", dpi=300)
    print("Plot saved as surfacecharge_potential.png")

def main():
    print("Starting contact layer analysis...")
    x_int = analyze_contact_layer()
    
    print("\nStarting surface charge vs potential analysis...")
    analyze_surface_charge_potential(x_int)
    
    print("\nAnalysis complete. Generated files:")
    print("1. contact_layer_density.png")
    print("2. surfacecharge_potential.png")
    print("3. predicted_potential.txt")
    print("4. theta_sigma.csv")
    print("5. sigma_potential.csv")

if __name__ == '__main__':
    main()

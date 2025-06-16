#!/usr/bin/env python3
import os
import glob
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, NullLocator
from matplotlib.ticker import MultipleLocator

def main():
    # Working in current directory containing numeric subfolders
    parent_dir = os.getcwd()
    subdirs = sorted(
        d for d in glob.glob(os.path.join(parent_dir, '*'))
        if os.path.isdir(d) and os.path.basename(d).isdigit()
    )
    if not subdirs:
        print('No numeric subdirectories found in', parent_dir)
        sys.exit(1)

    # --- Theta–sigma processing ---
    theta_frames = []
    for folder in subdirs:
        path = os.path.join(folder, 'theta_sigma.csv')
        if os.path.isfile(path):
            df = pd.read_csv(path)
            if df.shape[1] >= 2:
                theta_frames.append(df.iloc[:, :2].values)

    if not theta_frames:
        print('No valid theta_sigma.csv files found; exiting.')
        sys.exit(1)

    theta_all = np.stack(theta_frames, axis=0)
    x_vals = theta_all[:, :, 0]
    y_vals = theta_all[:, :, 1]
    x_mean = x_vals.mean(axis=0)
    y_mean = y_vals.mean(axis=0)
    y_std = y_vals.std(axis=0)

    abs_x = np.abs(x_mean)
    abs_y = np.abs(y_mean)

    pd.DataFrame({
        'sigma_mean': abs_x,
        'contact_layer_mean': abs_y,
        'contact_layer_std': y_std
    }).to_csv('theta_sigma_summary.csv', index=False)
    print('Wrote summary data to theta_sigma_summary.csv')

    # Axis‐limit helper (clamped to first quadrant)
    def padded_limits(arr, pad=10):
        lo = max(5 * np.floor((arr.min() - pad) / 5), 0)
        hi = 5 * np.ceil((arr.max() + pad) / 5)
        return lo, hi

    x_lo, x_hi = padded_limits(abs_x)
    y_lo, y_hi = padded_limits(abs_y)

    # Zero‐intercept fit on abs‐values
    n_total = len(abs_x)
    n = min(3, n_total)
    exclude_i = None
    while True:
        xi = abs_x[:n]
        yi = abs_y[:n]
        denom = np.dot(xi, xi)
        slope = np.dot(xi, yi) / denom if denom else 0
        num = np.dot(xi, yi)
        den = np.sqrt(denom * np.dot(yi, yi))
        r2 = (num/den)**2 if den else 0
        if r2 >= 0.9999 and n < n_total:
            n += 1
            continue
        if r2 < 0.9999:
            exclude_i = n - 1
            n -= 1
            xi = abs_x[:n]
            yi = abs_y[:n]
            denom = np.dot(xi, xi)
            slope = np.dot(xi, yi) / denom if denom else 0
        break

    if exclude_i is not None:
        y_excl = abs_y[exclude_i]
        x_int_abs = y_excl / slope if slope else 0
    else:
        x_int_abs = abs_x[n-1]
    sign0 = np.sign(x_mean[0]) if x_mean.size else 1
    x_int = x_int_abs * sign0

    # --- Plot 1: Contact layer vs. Surface charge ---
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.errorbar(
        abs_x, abs_y, yerr=y_std,
        fmt='o', color='blue', ecolor='black',
        capsize=4, markersize=6, linestyle='none', zorder=3
    )
    ax.plot(
        [x_lo, x_hi],
        [slope * x_lo, slope * x_hi],
        '--', linewidth=1.5, color='grey', alpha=0.7, zorder=2
    )
    if exclude_i is not None:
        ax.hlines(
            abs_y[exclude_i],
            x_lo, x_hi,
            linestyle='--', linewidth=1.5,
            color='grey', alpha=0.7, zorder=2
        )

    ax.set_xlabel('Surface charge density / μC·cm⁻²', fontsize=14, labelpad=8)
    ax.set_ylabel('Contact layer charge density / μC·cm⁻²', fontsize=14, labelpad=8)
    ax.set_title('Contact layer vs. Surface charge', fontsize=16, pad=12)
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(y_lo, y_hi)

    # show ticks with sign from means
    sign_x = np.sign(x_mean[0]) if x_mean.size else 1
    sign_y = np.sign(y_mean[0]) if y_mean.size else 1
    if sign_x < 0:
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, pos: f'-{int(v)}' if v > 0 else f'{int(v)}'))
    else:
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, pos: f'{int(v)}'))
    if sign_y < 0:
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, pos: f'-{int(v)}' if v > 0 else f'{int(v)}'))
    else:
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, pos: f'{int(v)}'))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_minor_locator(NullLocator())

    ax.legend([], [], title=f'$σ_M$ = {x_int:.2f} μC·cm⁻²',
              fontsize=12, title_fontsize=12, loc='upper left')
    fig.savefig('theta_sigma_plot.png', dpi=500)
    print('Saved plot to theta_sigma_plot.png')

    # --- Plot 2: Surface charge vs. Potential ---
    pot_frames = []
    for folder in subdirs:
        path = os.path.join(folder, 'sigma_potential.csv')
        if os.path.isfile(path):
            df = pd.read_csv(path)
            if df.shape[1] >= 2:
                pot_frames.append(df.values)

    if not pot_frames:
        print('No valid sigma_potential.csv files found; skipping potential plot.')
        return

    pot_all = np.stack(pot_frames, axis=0)
    sigma_vals2 = pot_all[:, :, 0]
    pot_vals2 = pot_all[:, :, 1]
    pot_std = pot_vals2.std(axis=0)

    pot_mean = pot_vals2.mean(axis=0)
    sigma_mean2 = sigma_vals2.mean(axis=0)

    abs_pot2 = np.abs(pot_mean)
    abs_sigma2 = np.abs(sigma_mean2)

    pd.DataFrame({
        'potential_mean': abs_pot2,
        'potential_std': pot_std,
        'sigma_mean': abs_sigma2
    }).to_csv('sigma_potential_summary.csv', index=False)
    print('Wrote summary data to sigma_potential_summary.csv')

    x2_lo, x2_hi = padded_limits(abs_pot2)
    y2_lo, y2_hi = padded_limits(abs_sigma2)

    fig2, ax2 = plt.subplots(figsize=(7, 5))
    ax2.errorbar(
        abs_pot2, abs_sigma2, xerr=pot_std,
        fmt='o', color='blue', ecolor='black',
        capsize=4, markersize=6, linestyle='none', zorder=3
    )
    ax2.hlines(
        x_int_abs, x2_lo, x2_hi,
        linestyle='--', linewidth=1.5, color='grey', alpha=0.7, zorder=2
    )

    ax2.set_xlabel('Potential / V', fontsize=14, labelpad=8)
    ax2.set_ylabel('Surface charge density / μC·cm⁻²', fontsize=14, labelpad=8)
    ax2.set_title('Surface charge vs. Potential', fontsize=16, pad=12)
    ax2.set_xlim(x2_lo, x2_hi)
    ax2.set_ylim(y2_lo, y2_hi)
    
    xi2 = abs_pot2[:n]
    yi2 = abs_sigma2[:n]
    slope2, intercept2 = np.polyfit(xi2, yi2, 1)
    x_cross = (x_int_abs - intercept2) / slope2 if slope2 else float('nan')
    
    
    #xi2 = abs_pot2[:n]
    #yi2 = abs_sigma2[:n]
    #denom2 = np.dot(xi2, xi2)
    #slope2 = np.dot(xi2, yi2) / denom2 if denom2 else 0
    #x_cross = x_int_abs / slope2 if slope2 else float('nan')

    x_line = [x2_lo, x2_hi]
    y_line = [slope2*x + intercept2 for x in x_line]
    # Draw regression line across the plot
    ax2.plot(
    x_line, y_line,
    '--', linewidth=1.5, color='grey', alpha=0.7, zorder=2
)
    #ax2.plot(
    #    [x2_lo, x2_hi],
    #    [slope2 * x2_lo, slope2 * x2_hi],
    #    '--', linewidth=1.5, color='grey', alpha=0.7, zorder=2
    #)
    
    
    # show ticks with sign from means
    sign_x2 = np.sign(pot_mean[0]) if pot_mean.size else 1
    sign_y2 = np.sign(sigma_mean2[0]) if sigma_mean2.size else 1
    if sign_x2 < 0:
        ax2.xaxis.set_major_formatter(FuncFormatter(lambda v, pos: f'-{int(v)}' if v > 0 else f'{int(v)}'))
    else:
        ax2.xaxis.set_major_formatter(FuncFormatter(lambda v, pos: f'{int(v)}'))
    if sign_y2 < 0:
        ax2.yaxis.set_major_formatter(FuncFormatter(lambda v, pos: f'-{int(v)}' if v > 0 else f'{int(v)}'))
    else:
        ax2.yaxis.set_major_formatter(FuncFormatter(lambda v, pos: f'{int(v)}'))
    ax2.xaxis.set_minor_locator(NullLocator())
    ax2.yaxis.set_minor_locator(NullLocator())

    ax2.tick_params(axis='both', which='major', labelsize=12)
    fig2.savefig('sigma_potential_plot.png', dpi=500)
    print('Saved plot to sigma_potential_plot.png')

    with open('uM.txt', 'w') as f:
        f.write(f'{x_cross:.6f}\n')
    print('Wrote crossing point to uM.txt')


if __name__ == '__main__':
    main()


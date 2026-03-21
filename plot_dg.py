#!/usr/bin/env python3
"""
Plot DG solution from one or more VTK files produced by DG1DIO::writeToVTK.

Usage:
    # Single snapshot
    python plot_dg.py output_000000.vtk

    # Multiple snapshots overlaid on one axes
    python plot_dg.py output_000000.vtk output_000010.vtk output_000020.vtk

    # Animate all frames (sorted automatically)
    python plot_dg.py --animate output_*.vtk

    # Compare with exact solution sin(x - a*t) for advection speed a=1
    python plot_dg.py --exact 1.0 output_000050.vtk
"""

import argparse
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


# ─────────────────────────────────────────────────────────────────────────────
#  Parser for our specific VTK format
# ─────────────────────────────────────────────────────────────────────────────

def read_vtk(filename):
    """
    Read a legacy ASCII VTK file written by DG1DIO::writeToVTK.
    Returns:
        x      – 1D numpy array of point x-coordinates
        u      – 1D numpy array of POINT_DATA scalar values
        K      – number of DG elements (= number of POLY_LINE cells)
        Np     – nodes per element
        t      – time extracted from the header comment (or None)
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    # strip blank lines and normalise
    lines = [l.strip() for l in lines if l.strip()]

    # --- extract time from header comment line (line 2) ---
    t = None
    header = lines[1]               # e.g. "DG1D solution: u  t=1.2345"
    if 't=' in header:
        try:
            t = float(header.split('t=')[-1].split()[0])
        except ValueError:
            pass

    # --- locate sections ---
    def find(keyword):
        for i, l in enumerate(lines):
            if l.upper().startswith(keyword.upper()):
                return i
        raise ValueError(f"Keyword '{keyword}' not found in {filename}")

    # POINTS
    pt_idx  = find('POINTS')
    n_pts   = int(lines[pt_idx].split()[1])
    coords  = np.array([
        float(lines[pt_idx + 1 + i].split()[0])   # only x; y and z are 0
        for i in range(n_pts)
    ])

    # CELLS – deduce Np from the first cell entry
    cell_idx = find('CELLS')
    K        = int(lines[cell_idx].split()[1])
    first_cell_line = lines[cell_idx + 1].split()
    Np       = int(first_cell_line[0])             # first number is vertex count

    # POINT_DATA scalar values
    scalar_idx = find('SCALARS')
    # data starts two lines after SCALARS (skip LOOKUP_TABLE line)
    data_start = scalar_idx + 2
    u = np.array([float(lines[data_start + i]) for i in range(n_pts)])

    return coords, u, K, Np, t


def split_elements(x, u, K, Np):
    x2 = x.reshape(K, Np)
    u2 = u.reshape(K, Np)
    # sort nodes within each element by physical x position
    # (GLL nodes are not stored in spatial order)
    idx = np.argsort(x2, axis=1)
    x2 = np.take_along_axis(x2, idx, axis=1)
    u2 = np.take_along_axis(u2, idx, axis=1)
    return x2, u2

# ─────────────────────────────────────────────────────────────────────────────
#  Plotting helpers
# ─────────────────────────────────────────────────────────────────────────────

def plot_snapshot(ax, filename, label=None, color=None, alpha=1.0):
    """Plot one VTK file onto ax. Returns the time value (or None)."""
    x, u, K, Np, t = read_vtk(filename)
    x2, u2 = split_elements(x, u, K, Np)

    lbl = label or (f"t = {t:.4f}" if t is not None else filename)
    first = True
    for k in range(K):
        ax.plot(x2[k], u2[k],
                color=color,
                alpha=alpha,
                label=lbl if first else None)
        first = False
        # vertical dashed line at each element boundary (except last)
        if k < K - 1:
            ax.axvline(x2[k, -1], color='lightgray', linewidth=0.5,
                       linestyle='--', zorder=0)
    return t


def add_exact(ax, x_range, t, a, n_pts=500, label='exact'):
    """Overlay sin(x - a*t) for linear advection."""
    if t is None:
        return
    xc = np.linspace(x_range[0], x_range[1], n_pts)
    ax.plot(xc, np.sin(xc - a * t), 'k--', linewidth=1, label=label)


# ─────────────────────────────────────────────────────────────────────────────
#  Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description='Plot DG1D VTK output')
    parser.add_argument('files', nargs='+', help='VTK file(s)')
    parser.add_argument('--animate', action='store_true',
                        help='Animate files as a time series')
    parser.add_argument('--exact', type=float, default=None, metavar='A',
                        help='Overlay exact solution sin(x-At) with wave speed A')
    parser.add_argument('--save', type=str, default=None,
                        help='Save figure to file (e.g. out.png or out.gif)')
    parser.add_argument('--interval', type=int, default=100,
                        help='Animation frame interval in ms (default 100)')
    args = parser.parse_args()

    # expand any globs that the shell didn't expand (Windows / quoted patterns)
    files = []
    for pattern in args.files:
        expanded = sorted(glob.glob(pattern))
        files.extend(expanded if expanded else [pattern])

    if not files:
        print("No files found.", file=sys.stderr)
        sys.exit(1)

    # ── Single / multi-overlay plot ─────────────────────────────────────────
    if not args.animate or len(files) == 1:
        fig, ax = plt.subplots(figsize=(9, 4))
        cmap = plt.cm.viridis
        colors = [cmap(i / max(len(files) - 1, 1)) for i in range(len(files))]

        t_last = None
        for i, fname in enumerate(files):
            t_last = plot_snapshot(ax, fname, color=colors[i],
                                   alpha=0.8 if len(files) > 1 else 1.0)

        if args.exact is not None:
            # use the last file's time for the exact solution overlay
            x0, _, _, _, _ = read_vtk(files[-1])
            add_exact(ax, (x0.min(), x0.max()), t_last, args.exact)

        ax.set_xlabel('x')
        ax.set_ylabel('u')
        ax.set_title('DG1D solution')
        if len(files) <= 8:
            ax.legend(fontsize=8)
        plt.tight_layout()

        if args.save:
            fig.savefig(args.save, dpi=150)
            print(f"Saved to {args.save}")
        else:
            plt.show()

    # ── Animation ───────────────────────────────────────────────────────────
    else:
        # read first file to set axis limits
        x0, u0, K, Np, _ = read_vtk(files[0])
        u_all = np.concatenate([read_vtk(f)[1] for f in files])
        u_min, u_max = u_all.min(), u_all.max()
        pad = (u_max - u_min) * 0.1 or 0.1

        fig, ax = plt.subplots(figsize=(9, 4))
        ax.set_xlim(x0.min(), x0.max())
        ax.set_ylim(u_min - pad, u_max + pad)
        ax.set_xlabel('x')
        ax.set_ylabel('u')

        x2, u2 = split_elements(x0, u0, K, Np)
        lines_plot = [ax.plot(x2[k], u2[k], 'steelblue')[0] for k in range(K)]
        exact_line, = ax.plot([], [], 'k--', linewidth=1, label='exact')
        title = ax.set_title('')
        for k in range(K - 1):
            ax.axvline(x2[k, -1], color='lightgray', linewidth=0.5,
                       linestyle='--', zorder=0)

        def update(frame_idx):
            fname = files[frame_idx]
            x, u, K2, Np2, t = read_vtk(fname)
            x2, u2 = split_elements(x, u, K2, Np2)
            for k, line in enumerate(lines_plot):
                line.set_data(x2[k], u2[k])
            title.set_text(f"t = {t:.4f}" if t is not None else fname)
            if args.exact is not None and t is not None:
                xc = np.linspace(x.min(), x.max(), 500)
                exact_line.set_data(xc, np.sin(xc - args.exact * t))
            return lines_plot + [title, exact_line]

        ani = animation.FuncAnimation(
            fig, update, frames=len(files),
            interval=args.interval, blit=False)

        if args.exact is not None:
            ax.legend(fontsize=8)

        if args.save:
            ani.save(args.save, writer='pillow', fps=1000 // args.interval)
            print(f"Saved animation to {args.save}")
        else:
            plt.show()


if __name__ == '__main__':
    main()

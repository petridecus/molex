#!/usr/bin/env python3
"""Cross-validate molex xtal pipeline against gemmi reference implementation.

Compares intermediate values at each pipeline stage:
  1. Fc amplitudes (same grid density -> FFT -> V/N normalization)
  2. Scaling parameters (Wilson k estimate)
  3. Sigma-A (D and sigma_sq per resolution bin)
  4. R-factors (R-work, R-free)

Uses gemmi's DensityCalculatorX for Fc (IT92 form factors, same as molex),
followed by numpy FFT, to produce bit-comparable structure factors.

Usage:
  python3 examples/cross_validate.py [--structures 1AKI:19 2H5C:154 ...]
  python3 examples/cross_validate.py --all
"""

import sys
import os
import argparse

import gemmi
import numpy as np


EXAMPLES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)))

ALL_STRUCTURES = [
    ("3LZT", 1),
    ("3NIR", 4),
    ("6E6O", 5),
    ("1AKI", 19),
    ("7KOM", 23),
    ("1X6Z", 92),
    ("1G6X", 96),
    ("4XDX", 152),
    ("2H5C", 154),
    ("3ZUC", 178),
]


def load_reflections(sf_path: str):
    """Parse SF-CIF -> (h, k, l, fobs, sigma, free) arrays."""
    doc = gemmi.cif.read(sf_path)
    block = doc[0]

    h = np.array([int(x) for x in block.find_values("_refln.index_h")])
    k = np.array([int(x) for x in block.find_values("_refln.index_k")])
    l = np.array([int(x) for x in block.find_values("_refln.index_l")])

    fobs = None
    for col in ["_refln.F_meas_au", "_refln.F_obs"]:
        try:
            vals = list(block.find_values(col))
            fobs = np.array([float(x) if x not in ("?", ".") else 0.0 for x in vals])
            break
        except Exception:
            continue
    if fobs is None:
        for col in ["_refln.intensity_meas", "_refln.I_obs"]:
            try:
                vals = list(block.find_values(col))
                ints = np.array([float(x) if x not in ("?", ".") else 0.0 for x in vals])
                fobs = np.sqrt(np.maximum(ints, 0.0))
                break
            except Exception:
                continue
    if fobs is None:
        raise ValueError("No Fobs column found")

    sigma = np.ones_like(fobs)
    for col in ["_refln.F_meas_sigma_au", "_refln.F_obs_sigma"]:
        try:
            vals = list(block.find_values(col))
            sigma = np.array([float(x) if x not in ("?", ".") else 1.0 for x in vals])
            break
        except Exception:
            continue

    free = np.zeros(len(h), dtype=bool)
    try:
        flags = list(block.find_values("_refln.pdbx_r_free_flag"))
        free = np.array([str(x) == "1" for x in flags])
    except Exception:
        try:
            status = list(block.find_values("_refln.status"))
            free = np.array([str(x).lower() == "f" for x in status])
        except Exception:
            pass

    valid = fobs > 0
    n_valid = int(np.sum(valid))
    free_v = free[valid] if len(free) == len(h) else np.zeros(n_valid, dtype=bool)

    # If no free flags, assign ~5% randomly
    if np.sum(free_v) == 0:
        rng = np.random.default_rng(42)
        free_v = rng.random(n_valid) < 0.05

    return h[valid], k[valid], l[valid], fobs[valid], sigma[valid], free_v


def compute_fc_grid(st, sg_number, d_min, grid_dims=None):
    """Compute Fc via density splatting + FFT (same method as molex).

    Returns the full complex FFT array (un-half, V/N normalized).
    """
    dc = gemmi.DensityCalculatorX()
    dc.d_min = d_min
    dc.rate = 1.5  # same as molex oversampling
    dc.set_grid_cell_and_spacegroup(st)

    if grid_dims is not None:
        dc.grid.set_size(grid_dims[0], grid_dims[1], grid_dims[2])

    dc.put_model_density_on_grid(st[0])
    grid = dc.grid

    arr = np.array(grid, copy=False)
    fc_full = np.fft.fftn(arr)

    V = st.cell.volume
    N = arr.size
    fc_full *= V / N  # crystallographic normalization

    return fc_full, (grid.nu, grid.nv, grid.nw)


def extract_fc_at_hkl(fc_grid, h, k, l, dims):
    """Extract Fc values at specific Miller indices from the full FFT grid."""
    nu, nv, nw = dims
    fc_vals = np.zeros(len(h), dtype=complex)
    for i in range(len(h)):
        u = int(h[i]) % nu
        v = int(k[i]) % nv
        w = int(l[i]) % nw
        fc_vals[i] = fc_grid[u, v, w]
    return fc_vals


def wilson_k(fobs, fc_amp):
    """Estimate k_overall from Wilson plot (median ratio)."""
    valid = (fobs > 1) & (fc_amp > 0.01)
    if np.sum(valid) < 20:
        return 1.0
    return float(np.median(fobs[valid] / fc_amp[valid]))


def r_factor(fobs, fc_amp, k, mask):
    """R = sum|Fo - k*Fc| / sum|Fo| over mask."""
    fo = fobs[mask]
    fc = fc_amp[mask]
    denom = np.sum(np.abs(fo))
    if denom < 1e-10:
        return 0.0
    return float(np.sum(np.abs(fo - k * fc)) / denom)


def sigma_a_bins(fobs, fc_amp, cell, h, k, l, free, n_bins=20):
    """Compute sigma-A D and sigma_sq per resolution bin."""
    s2 = np.array([cell.calculate_1_d2([int(h[i]), int(k[i]), int(l[i])])
                    for i in range(len(h))])
    work = ~free
    s2_w, fo_w, fc_w = s2[work], fobs[work], fc_amp[work]

    s2_min, s2_max = s2_w.min(), s2_w.max()
    bw = (s2_max - s2_min) / n_bins

    d_bins = np.full(n_bins, np.nan)
    sq_bins = np.full(n_bins, np.nan)

    for b in range(n_bins):
        lo = s2_min + b * bw
        mask = (s2_w >= lo) & (s2_w < lo + bw)
        if np.sum(mask) < 3:
            continue
        fo_b, fc_b = fo_w[mask], fc_w[mask]
        d = np.clip(np.sum(fo_b * fc_b) / np.sum(fc_b ** 2), 0, 1)
        d_bins[b] = d
        sq_bins[b] = np.mean((fo_b - d * fc_b) ** 2)

    return d_bins, sq_bins


def validate_one(pdb: str, sg_number: int, verbose=True):
    """Run full cross-validation for one structure."""
    coord_path = os.path.join(EXAMPLES_DIR, f"{pdb}.cif")
    sf_path = os.path.join(EXAMPLES_DIR, f"{pdb}-sf.cif")

    if not os.path.exists(coord_path) or not os.path.exists(sf_path):
        return {"pdb": pdb, "error": "files not found"}

    st = gemmi.read_structure(coord_path)
    st.setup_entities()
    cell = st.cell

    h, k, l, fobs, sigma, free = load_reflections(sf_path)

    s2_all = np.array([cell.calculate_1_d2([int(h[i]), int(k[i]), int(l[i])])
                        for i in range(len(h))])
    d_min = 1.0 / np.sqrt(s2_all.max())

    if verbose:
        print(f"\n{'='*60}")
        print(f"  {pdb} | SG {sg_number} | {len(h)} reflections | d_min={d_min:.2f} A")
        print(f"{'='*60}")

    # ── Stage 1: Fc calculation ──────────────────────────────────────

    fc_grid, dims = compute_fc_grid(st, sg_number, d_min)
    fc_complex = extract_fc_at_hkl(fc_grid, h, k, l, dims)
    fc_amp = np.abs(fc_complex)

    if verbose:
        print(f"\n  Stage 1: Fc Calculation (density splat + FFT, V/N normalized)")
        print(f"    Grid: {dims}")
        print(f"    Fc mean:   {fc_amp.mean():.2f}")
        print(f"    Fc median: {np.median(fc_amp):.2f}")
        print(f"    Fc range:  [{fc_amp.min():.2f}, {fc_amp.max():.2f}]")
        print(f"    Fobs mean: {fobs.mean():.2f}")

    # ── Stage 2: Scaling ─────────────────────────────────────────────

    k_ov = wilson_k(fobs, fc_amp)
    work, test = ~free, free
    rw = r_factor(fobs, fc_amp, k_ov, work)
    rf = r_factor(fobs, fc_amp, k_ov, test)

    if verbose:
        print(f"\n  Stage 2: Scaling")
        print(f"    k_overall (Wilson): {k_ov:.4f}")
        print(f"    R-work:  {rw:.4f}")
        print(f"    R-free:  {rf:.4f}")
        print(f"    Gap:     {rf - rw:.4f}")

    # ── Stage 3: Sigma-A ─────────────────────────────────────────────

    fc_scaled = fc_amp * k_ov
    d_vals, sq_vals = sigma_a_bins(fobs, fc_scaled, cell, h, k, l, free)

    if verbose:
        print(f"\n  Stage 3: Sigma-A (on k-scaled Fc)")
        valid_d = d_vals[~np.isnan(d_vals)]
        valid_sq = sq_vals[~np.isnan(sq_vals)]
        if len(valid_d) > 0:
            print(f"    D range:       [{valid_d.min():.4f}, {valid_d.max():.4f}]")
            print(f"    D mean:        {valid_d.mean():.4f}")
            print(f"    sigma_sq range: [{valid_sq.min():.4f}, {valid_sq.max():.4f}]")
            print(f"    D bins[0:5]:   {d_vals[:5]}")
            print(f"    sq bins[0:5]:  {sq_vals[:5]}")

    # ── Stage 4: Map coefficients (simplified 2Fo-DFc) ───────────────

    # Use bin-averaged D for a simplified map coefficient computation
    d_mean = np.nanmean(d_vals) if np.any(~np.isnan(d_vals)) else 0.5
    map_amp = 2.0 * fobs - d_mean * fc_scaled
    frac_neg = np.mean(map_amp < 0)

    if verbose:
        print(f"\n  Stage 4: Map Coefficients (simplified 2Fo-DFc)")
        print(f"    D_mean used: {d_mean:.4f}")
        print(f"    2Fo-DFc mean:      {map_amp.mean():.2f}")
        print(f"    2Fo-DFc neg frac:  {frac_neg:.3f}")

    return {
        "pdb": pdb,
        "sg": sg_number,
        "n_refl": len(h),
        "d_min": d_min,
        "grid": dims,
        "fc_mean": float(fc_amp.mean()),
        "fc_median": float(np.median(fc_amp)),
        "fobs_mean": float(fobs.mean()),
        "k_overall": k_ov,
        "r_work": rw,
        "r_free": rf,
        "d_mean": float(d_mean),
        "sigma_sq_mean": float(np.nanmean(sq_vals)) if np.any(~np.isnan(sq_vals)) else None,
    }


def main():
    parser = argparse.ArgumentParser(description="Cross-validate molex vs gemmi")
    parser.add_argument("--all", action="store_true", help="Run all 10 structures")
    parser.add_argument("structures", nargs="*",
                        help="PDB:SG pairs (e.g. 1AKI:19)")
    args = parser.parse_args()

    if args.all:
        pairs = ALL_STRUCTURES
    elif args.structures:
        pairs = []
        for s in args.structures:
            pdb, sg = s.split(":")
            pairs.append((pdb, int(sg)))
    else:
        pairs = [("1AKI", 19)]

    print("=" * 60)
    print("  GEMMI REFERENCE CROSS-VALIDATION")
    print("  Method: DensityCalculatorX (IT92) + numpy FFT + V/N norm")
    print("=" * 60)

    results = []
    for pdb, sg in pairs:
        try:
            r = validate_one(pdb, sg, verbose=True)
            results.append(r)
        except Exception as e:
            print(f"\n  {pdb}: ERROR - {e}")
            results.append({"pdb": pdb, "error": str(e)})

    # Summary table
    print(f"\n\n{'='*76}")
    print("  GEMMI REFERENCE SUMMARY")
    print(f"{'='*76}")
    print(f"{'PDB':<6} {'SG':>4} {'Refls':>7} {'Fc_mean':>8} {'k_ov':>7} "
          f"{'Rwork':>7} {'Rfree':>7} {'D_mean':>7}")
    print("-" * 76)
    for r in results:
        if "error" in r:
            print(f"{r['pdb']:<6} ERROR: {r['error']}")
        else:
            print(f"{r['pdb']:<6} {r['sg']:>4} {r['n_refl']:>7} "
                  f"{r['fc_mean']:>8.2f} {r['k_overall']:>7.4f} "
                  f"{r['r_work']:>7.4f} {r['r_free']:>7.4f} "
                  f"{r['d_mean']:>7.4f}")


if __name__ == "__main__":
    main()

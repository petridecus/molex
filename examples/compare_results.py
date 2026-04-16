#!/usr/bin/env python3
"""Compare molex xtal pipeline output against gemmi reference, stage by stage.

Runs both molex (via dump_intermediates) and gemmi on each structure,
then prints a side-by-side comparison of intermediate values.

Usage:
  python3 examples/compare_results.py
"""

import subprocess
import os
import sys

import gemmi
import numpy as np


EXAMPLES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)))
MOLEX_BIN = os.path.join(
    EXAMPLES_DIR, "..", "target", "release", "examples", "dump_intermediates"
)

STRUCTURES = [
    ("6E6O", 5, "C 1 2 1"),
    ("1AKI", 19, "P 21 21 21"),
    ("7KOM", 23, "I 2 2 2"),
    ("1X6Z", 92, "P 41 21 2"),
    ("3ZUC", 178, "P 61 2 2"),
]


def run_molex(pdb, sg):
    """Run molex dump_intermediates and parse key=value output."""
    coord = os.path.join(EXAMPLES_DIR, f"{pdb}.cif")
    sf = os.path.join(EXAMPLES_DIR, f"{pdb}-sf.cif")
    result = subprocess.run(
        [MOLEX_BIN, coord, sf, str(sg)],
        capture_output=True, text=True, timeout=120,
    )
    if result.returncode != 0:
        return {"error": result.stderr.strip()[:200]}

    vals = {}
    for line in result.stdout.strip().split("\n"):
        if "=" in line:
            key, val = line.split("=", 1)
            vals[key] = val
    return vals


def run_gemmi(pdb, sg):
    """Compute reference values using gemmi."""
    coord_path = os.path.join(EXAMPLES_DIR, f"{pdb}.cif")
    sf_path = os.path.join(EXAMPLES_DIR, f"{pdb}-sf.cif")

    st = gemmi.read_structure(coord_path)
    st.setup_entities()
    cell = st.cell

    # Load reflections
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
        return {"error": "no Fobs"}

    valid = fobs > 0
    h, k, l, fobs = h[valid], k[valid], l[valid], fobs[valid]

    # Free flags
    free = np.zeros(len(h), dtype=bool)
    try:
        flags = list(block.find_values("_refln.pdbx_r_free_flag"))
        free_raw = np.array([str(x) == "1" for x in flags])
        free = free_raw[valid] if len(free_raw) == len(valid) else free
    except Exception:
        pass
    try:
        status = list(block.find_values("_refln.status"))
        free_raw = np.array([str(x).lower() == "f" for x in status])
        if np.sum(free_raw) > 0:
            free = free_raw[valid] if len(free_raw) == len(valid) else free
    except Exception:
        pass
    if np.sum(free) == 0:
        rng = np.random.default_rng(42)
        free = rng.random(len(h)) < 0.05

    # Compute Fc via density grid + FFT
    dc = gemmi.DensityCalculatorX()
    s2_all = np.array([cell.calculate_1_d2([int(h[i]), int(k[i]), int(l[i])])
                        for i in range(len(h))])
    d_min = 1.0 / np.sqrt(s2_all.max())
    dc.d_min = d_min
    dc.rate = 1.5
    dc.set_grid_cell_and_spacegroup(st)
    dc.put_model_density_on_grid(st[0])

    grid = dc.grid
    dims = (grid.nu, grid.nv, grid.nw)
    arr = np.array(grid, copy=False)
    fc_full = np.fft.fftn(arr)
    V = cell.volume
    N = arr.size
    fc_full *= V / N

    # Extract Fc at reflection positions
    fc_vals = np.zeros(len(h), dtype=complex)
    for i in range(len(h)):
        u = int(h[i]) % dims[0]
        v = int(k[i]) % dims[1]
        w = int(l[i]) % dims[2]
        fc_vals[i] = fc_full[u, v, w]
    fc_amp = np.abs(fc_vals)

    # Wilson k
    ok = (fobs > 1) & (fc_amp > 0.01)
    k_ov = float(np.median(fobs[ok] / fc_amp[ok])) if np.sum(ok) > 20 else 1.0

    # R-factors
    work = ~free
    rw_num = np.sum(np.abs(fobs[work] - k_ov * fc_amp[work]))
    rw_den = np.sum(np.abs(fobs[work]))
    r_work = float(rw_num / rw_den) if rw_den > 0 else 0.0

    rf_num = np.sum(np.abs(fobs[free] - k_ov * fc_amp[free]))
    rf_den = np.sum(np.abs(fobs[free]))
    r_free = float(rf_num / rf_den) if rf_den > 0 else 0.0

    # Sigma-A on k-scaled Fc
    fc_sc = fc_amp * k_ov
    s2_w = s2_all[work]
    fo_w = fobs[work]
    fc_w = fc_sc[work]
    s2_min_val = s2_w.min()
    s2_max_val = s2_w.max()
    bw = (s2_max_val - s2_min_val) / 20
    d_bins = []
    sq_bins = []
    for b in range(20):
        lo = s2_min_val + b * bw
        mask = (s2_w >= lo) & (s2_w < lo + bw)
        if np.sum(mask) < 3:
            d_bins.append(np.nan)
            sq_bins.append(np.nan)
            continue
        fo_b, fc_b = fo_w[mask], fc_w[mask]
        d = np.clip(np.sum(fo_b * fc_b) / np.sum(fc_b**2), 0, 1)
        d_bins.append(d)
        sq_bins.append(np.mean((fo_b - d * fc_b) ** 2))

    return {
        "grid": f"{dims[0]},{dims[1]},{dims[2]}",
        "n_refl": len(h),
        "fc_mean": float(fc_amp.mean()),
        "k_overall": k_ov,
        "r_work": r_work,
        "r_free": r_free,
        "d_bins": d_bins,
        "sigma_sq_bins": sq_bins,
    }


def fmt(val, width=10, prec=4):
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return f"{'N/A':>{width}}"
    if isinstance(val, float):
        return f"{val:>{width}.{prec}f}"
    return f"{str(val):>{width}}"


def compare(pdb, sg, sg_name):
    print(f"\n{'='*72}")
    print(f"  {pdb}  SG {sg} ({sg_name})")
    print(f"{'='*72}")

    molex = run_molex(pdb, sg)
    if "error" in molex:
        print(f"  molex ERROR: {molex['error']}")
        return False

    ref = run_gemmi(pdb, sg)
    if "error" in ref:
        print(f"  gemmi ERROR: {ref['error']}")
        return False

    # Grid
    m_grid = molex.get("MOLEX_GRID", "?")
    g_grid = ref.get("grid", "?")
    grid_match = m_grid == g_grid
    print(f"\n  Grid dimensions:    molex={m_grid}  gemmi={g_grid}  {'MATCH' if grid_match else 'DIFFER'}")

    # Scaling
    m_k = float(molex.get("MOLEX_K_OVERALL", "0"))
    g_k = ref["k_overall"]
    # molex uses LM + bulk solvent; gemmi uses Wilson median. Compare ratio.
    print(f"\n  {'Parameter':<20} {'molex':>12} {'gemmi':>12} {'ratio':>8} {'notes':}")
    print(f"  {'-'*68}")
    print(f"  {'k_overall':<20} {m_k:>12.4f} {g_k:>12.4f} {m_k/g_k if g_k else 0:>8.3f}  "
          f"{'(molex: LM+solvent, gemmi: Wilson median)'}")

    m_rw = float(molex.get("MOLEX_R_WORK", "0"))
    g_rw = ref["r_work"]
    m_rf = float(molex.get("MOLEX_R_FREE", "0"))
    g_rf = ref["r_free"]
    rw_diff = abs(m_rw - g_rw)
    rf_diff = abs(m_rf - g_rf)
    print(f"  {'R-work':<20} {m_rw:>12.4f} {g_rw:>12.4f} {'':>8}  diff={rw_diff:.4f}")
    print(f"  {'R-free':<20} {m_rf:>12.4f} {g_rf:>12.4f} {'':>8}  diff={rf_diff:.4f}")

    m_ksol = float(molex.get("MOLEX_K_SOL", "0"))
    m_bsol = float(molex.get("MOLEX_B_SOL", "0"))
    print(f"  {'k_sol':<20} {m_ksol:>12.4f} {'N/A':>12}  (gemmi: no bulk solvent model)")
    print(f"  {'B_sol':<20} {m_bsol:>12.1f} {'N/A':>12}")

    # Sigma-A D bins
    m_d_str = molex.get("MOLEX_D_BINS", "")
    m_d = [float(x) for x in m_d_str.split(",") if x] if m_d_str else []
    g_d = ref.get("d_bins", [])

    print(f"\n  Sigma-A D (first 5 bins):")
    print(f"  {'bin':<6} {'molex':>10} {'gemmi':>10}")
    for i in range(min(5, len(m_d), len(g_d))):
        g_val = g_d[i] if not np.isnan(g_d[i]) else "N/A"
        print(f"  {i:<6} {m_d[i]:>10.4f} {fmt(g_val, 10)}")

    # Sigma-sq bins
    m_sq_str = molex.get("MOLEX_SIGMA_SQ_BINS", "")
    m_sq = [float(x) for x in m_sq_str.split(",") if x] if m_sq_str else []
    g_sq = ref.get("sigma_sq_bins", [])

    print(f"\n  Sigma-A sigma_sq (first 5 bins):")
    print(f"  {'bin':<6} {'molex':>12} {'gemmi':>12} {'ratio':>8}")
    for i in range(min(5, len(m_sq), len(g_sq))):
        g_val = g_sq[i]
        if not np.isnan(g_val) and g_val > 0:
            ratio = m_sq[i] / g_val
            print(f"  {i:<6} {m_sq[i]:>12.1f} {g_val:>12.1f} {ratio:>8.3f}")
        else:
            print(f"  {i:<6} {m_sq[i]:>12.1f} {'N/A':>12}")

    # Verdict
    r_match = rw_diff < 0.15  # within 15 percentage points
    print(f"\n  R-work within 0.15 of gemmi: {'YES' if r_match else 'NO'}")
    return r_match


def main():
    if not os.path.exists(MOLEX_BIN):
        print(f"Build first: cargo build --example dump_intermediates --features xtal --release")
        sys.exit(1)

    print("=" * 72)
    print("  MOLEX vs GEMMI CROSS-VALIDATION")
    print("  molex: density splat + FFT + V/N norm + LM scaling + bulk solvent")
    print("  gemmi: DensityCalculatorX + numpy FFT + V/N norm + Wilson scaling")
    print("=" * 72)

    n_pass = 0
    n_total = 0
    for pdb, sg, sg_name in STRUCTURES:
        n_total += 1
        if compare(pdb, sg, sg_name):
            n_pass += 1

    print(f"\n\n{'='*72}")
    print(f"  RESULT: {n_pass}/{n_total} structures within tolerance")
    print(f"{'='*72}")

    if n_pass < n_total:
        sys.exit(1)


if __name__ == "__main__":
    main()

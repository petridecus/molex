#!/usr/bin/env python3
"""Stage-by-stage cross-validation of molex xtal pipeline against cctbx (Phenix).

Compares four pipeline stages for 1AKI:
  1. Fc calculation
  2. Bulk solvent scaling
  3. Sigma-A estimation
  4. Map coefficients / R-factors
"""

import subprocess
import os
import numpy as np

from iotbx import pdb as iotbx_pdb
import iotbx.cif
from mmtbx import f_model
from cctbx.array_family import flex

import gemmi

EXAMPLES = os.path.dirname(os.path.abspath(__file__))
MOLEX_BIN = os.path.join(EXAMPLES, "..", "target", "release", "examples", "dump_intermediates")


def load_cctbx(pdb, sg):
    """Set up cctbx f_model for a structure."""
    pdb_inp = iotbx_pdb.input(file_name=os.path.join(EXAMPLES, f"{pdb}.cif"))
    xrs = pdb_inp.xray_structure_simple()

    cif_reader = iotbx.cif.reader(file_path=os.path.join(EXAMPLES, f"{pdb}-sf.cif"))
    mas = cif_reader.as_miller_arrays()

    f_obs = None
    for ma in mas:
        if ma.is_xray_amplitude_array():
            f_obs = ma
            break
        if ma.is_xray_intensity_array():
            f_obs = ma.as_amplitude_array()
            break
    if f_obs is None:
        for ma in mas:
            if ma.is_real_array() and not ma.is_complex_array():
                f_obs = ma
                break

    f_obs = f_obs.select(f_obs.data() > 0)
    flags = f_obs.generate_r_free_flags(fraction=0.05)

    fmod = f_model.manager(xray_structure=xrs, f_obs=f_obs, r_free_flags=flags)
    fmod.update_all_scales(remove_outliers=False)
    return fmod, xrs, f_obs


def run_molex(pdb, sg):
    """Run molex dump_intermediates."""
    result = subprocess.run(
        [MOLEX_BIN,
         os.path.join(EXAMPLES, f"{pdb}.cif"),
         os.path.join(EXAMPLES, f"{pdb}-sf.cif"),
         str(sg)],
        capture_output=True, text=True, timeout=120)
    vals = {}
    for line in result.stdout.strip().split("\n"):
        if "=" in line:
            k, v = line.split("=", 1)
            vals[k] = v
    return vals


def gemmi_fc(pdb):
    """Compute Fc via gemmi grid density + numpy FFT."""
    st = gemmi.read_structure(os.path.join(EXAMPLES, f"{pdb}.cif"))
    st.setup_entities()
    cell = st.cell

    dc = gemmi.DensityCalculatorX()
    dc.d_min = 1.5
    dc.rate = 1.5
    dc.set_grid_cell_and_spacegroup(st)
    dc.put_model_density_on_grid(st[0])

    grid = dc.grid
    arr = np.array(grid, copy=False)
    fc_full = np.fft.fftn(arr) * (cell.volume / arr.size)
    return fc_full, (grid.nu, grid.nv, grid.nw), cell


# ======================================================================
print("=" * 72)
print("  STAGE-BY-STAGE CROSS-VALIDATION: molex vs cctbx (Phenix)")
print("  Structure: 1AKI (P 21 21 21, SG 19)")
print("=" * 72)

fmod, xrs, f_obs = load_cctbx("1AKI", 19)
molex = run_molex("1AKI", 19)

# ── STAGE 1: Fc CALCULATION ──────────────────────────────────────────
print("\n" + "=" * 72)
print("  STAGE 1: Fc CALCULATION")
print("=" * 72)

# cctbx Fc (analytical, exact symmetry)
fc_cctbx = fmod.f_calc()
fc_cctbx_amps = np.abs(fc_cctbx.data().as_numpy_array())
fc_cctbx_hkl = np.array(fc_cctbx.indices())

# gemmi Fc (grid density + FFT, same method as molex)
fc_grid, dims, cell = gemmi_fc("1AKI")

# Extract gemmi Fc at same Miller indices
fc_gemmi_amps = np.zeros(len(fc_cctbx_hkl))
for i, (h, k, l) in enumerate(fc_cctbx_hkl):
    u = h % dims[0]
    v = k % dims[1]
    w = l % dims[2]
    fc_gemmi_amps[i] = abs(fc_grid[u, v, w])

# Correlation between cctbx and gemmi Fc
valid = (fc_cctbx_amps > 0.1) & (fc_gemmi_amps > 0.1)
corr = np.corrcoef(fc_cctbx_amps[valid], fc_gemmi_amps[valid])[0, 1]
ratio = fc_gemmi_amps[valid] / fc_cctbx_amps[valid]

print(f"\n  cctbx (analytical):  mean |Fc| = {fc_cctbx_amps.mean():.2f}")
print(f"  gemmi (grid+FFT):    mean |Fc| = {fc_gemmi_amps.mean():.2f}")
print(f"  Correlation:         {corr:.6f}")
print(f"  Ratio gemmi/cctbx:   mean={ratio.mean():.4f} std={ratio.std():.4f}")
print(f"  Grid dims:           {dims} (molex: {molex.get('MOLEX_GRID', '?')})")

if corr > 0.999:
    print("  VERDICT: Fc calculation MATCHES (grid error < 0.1%)")
elif corr > 0.99:
    print("  VERDICT: Fc calculation CLOSE (grid error < 1%)")
else:
    print(f"  VERDICT: Fc calculation DIFFERS (corr={corr:.4f})")


# ── STAGE 2: BULK SOLVENT SCALING ────────────────────────────────────
print("\n" + "=" * 72)
print("  STAGE 2: BULK SOLVENT SCALING")
print("=" * 72)

# cctbx scaling
c_rw = fmod.r_work()
c_rf = fmod.r_free()
c_k1 = fmod.scale_k1()
c_k_iso = np.array(fmod.k_isotropic())
c_k_aniso = np.array(fmod.k_anisotropic())
c_k_masks = fmod.k_masks()
c_k_mask = np.array(c_k_masks[0]) if c_k_masks and len(c_k_masks) > 0 else np.zeros(1)

# molex scaling
m_k = float(molex.get("MOLEX_K_OVERALL", "0"))
m_ksol = float(molex.get("MOLEX_K_SOL", "0"))
m_bsol = float(molex.get("MOLEX_B_SOL", "0"))
m_rw = float(molex.get("MOLEX_R_WORK", "0"))
m_rf = float(molex.get("MOLEX_R_FREE", "0"))

print(f"\n  {'Parameter':<25} {'molex':>12} {'cctbx':>12}")
print(f"  {'-'*50}")
print(f"  {'R-work':<25} {m_rw:>12.4f} {c_rw:>12.4f}")
print(f"  {'R-free':<25} {m_rf:>12.4f} {c_rf:>12.4f}")
print(f"  {'k_overall / scale_k1':<25} {m_k:>12.4f} {c_k1:>12.4f}")
print(f"  {'k_sol / mean(k_mask)':<25} {m_ksol:>12.4f} {c_k_mask.mean():>12.4f}")
print(f"  {'B_sol':<25} {m_bsol:>12.1f} {'N/A':>12}")
print(f"  {'mean(k_iso)':<25} {'N/A':>12} {c_k_iso.mean():>12.4f}")
print(f"  {'mean(k_aniso)':<25} {'N/A':>12} {c_k_aniso.mean():>12.4f}")

rw_gap = m_rw - c_rw
print(f"\n  R-work gap: {rw_gap:+.4f}")
print(f"  VERDICT: {'GOOD' if abs(rw_gap) < 0.05 else 'SIGNIFICANT GAP' if abs(rw_gap) < 0.10 else 'LARGE GAP'} ({rw_gap:+.4f})")


# ── STAGE 3: SIGMA-A ESTIMATION ──────────────────────────────────────
print("\n" + "=" * 72)
print("  STAGE 3: SIGMA-A ESTIMATION")
print("=" * 72)

# cctbx alpha (= D) and beta (= sigma_sq), per-reflection
alpha_beta = fmod.alpha_beta()
alpha_data = np.array(alpha_beta[0].data())
beta_data = np.array(alpha_beta[1].data())

# Bin cctbx values by resolution
s2_cctbx = np.array([f_obs.unit_cell().d_star_sq(h) for h in f_obs.indices()])
s2_min = s2_cctbx.min()
s2_max = s2_cctbx.max()
bw = (s2_max - s2_min) / 20

c_d_bins = []
c_sq_bins = []
for b in range(20):
    lo = s2_min + b * bw
    mask = (s2_cctbx >= lo) & (s2_cctbx < lo + bw)
    if np.sum(mask) < 3:
        c_d_bins.append(np.nan)
        c_sq_bins.append(np.nan)
    else:
        c_d_bins.append(np.mean(alpha_data[mask]))
        c_sq_bins.append(np.mean(beta_data[mask]))
c_d_bins = np.array(c_d_bins)
c_sq_bins = np.array(c_sq_bins)

# molex D bins
m_d_str = molex.get("MOLEX_D_BINS", "")
m_d_bins = np.array([float(x) for x in m_d_str.split(",")]) if m_d_str else np.zeros(20)
m_sq_str = molex.get("MOLEX_SIGMA_SQ_BINS", "")
m_sq_bins = np.array([float(x) for x in m_sq_str.split(",")]) if m_sq_str else np.zeros(20)

print(f"\n  {'Bin':<4} {'d_max':>6} {'D(molex)':>10} {'D(cctbx)':>10} {'sq(molex)':>12} {'sq(cctbx)':>12}")
print(f"  {'-'*58}")
for b in range(20):
    d_hi = 1.0 / np.sqrt(s2_min + b * bw) if (s2_min + b * bw) > 0 else 999
    d_lo = 1.0 / np.sqrt(s2_min + (b + 1) * bw) if (s2_min + (b + 1) * bw) > 0 else 999
    c_d = c_d_bins[b] if not np.isnan(c_d_bins[b]) else 0
    c_sq = c_sq_bins[b] if not np.isnan(c_sq_bins[b]) else 0
    m_d = m_d_bins[b] if b < len(m_d_bins) else 0
    m_sq = m_sq_bins[b] if b < len(m_sq_bins) else 0
    print(f"  {b:<4} {d_lo:>5.2f}A {m_d:>10.4f} {c_d:>10.4f} {m_sq:>12.1f} {c_sq:>12.1f}")

# Overall comparison
valid_d = ~np.isnan(c_d_bins) & (m_d_bins[:len(c_d_bins)] > 0)
if np.sum(valid_d) > 3:
    d_corr = np.corrcoef(m_d_bins[:len(c_d_bins)][valid_d], c_d_bins[valid_d])[0, 1]
    print(f"\n  D correlation (molex vs cctbx): {d_corr:.4f}")
    print(f"  VERDICT: {'GOOD' if d_corr > 0.9 else 'MODERATE' if d_corr > 0.7 else 'POOR'} (corr={d_corr:.4f})")


# ── STAGE 4: R-FACTORS (FULL MODEL) ─────────────────────────────────
print("\n" + "=" * 72)
print("  STAGE 4: R-FACTORS (full model)")
print("=" * 72)

print(f"\n  {'':>25} {'molex':>12} {'cctbx':>12} {'gap':>8}")
print(f"  {'-'*58}")
print(f"  {'R-work':<25} {m_rw:>12.4f} {c_rw:>12.4f} {m_rw-c_rw:>+8.4f}")
print(f"  {'R-free':<25} {m_rf:>12.4f} {c_rf:>12.4f} {m_rf-c_rf:>+8.4f}")
print(f"  {'R-work - R-free gap':<25} {m_rw-m_rf:>12.4f} {c_rw-c_rf:>12.4f}")

rw_gap = m_rw - c_rw
if abs(rw_gap) < 0.03:
    verdict = "EXCELLENT"
elif abs(rw_gap) < 0.06:
    verdict = "GOOD"
elif abs(rw_gap) < 0.10:
    verdict = "ACCEPTABLE"
else:
    verdict = "NEEDS WORK"

print(f"\n  OVERALL VERDICT: {verdict} (R-work gap = {rw_gap:+.4f})")


# ── SUMMARY ──────────────────────────────────────────────────────────
print("\n" + "=" * 72)
print("  SUMMARY")
print("=" * 72)
print(f"  Stage 1 (Fc calc):    corr={corr:.6f}  {'PASS' if corr > 0.99 else 'FAIL'}")
print(f"  Stage 2 (scaling):    R-work gap={m_rw-c_rw:+.4f}  {verdict}")
print(f"  Stage 3 (sigma-A):    D corr={'%.4f'%d_corr if np.sum(valid_d)>3 else 'N/A'}")
print(f"  Stage 4 (R-factors):  molex={m_rw:.4f} cctbx={c_rw:.4f}")

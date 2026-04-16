#!/usr/bin/env python3
"""Cross-validate Fc (structure factor) calculations: cctbx vs gemmi vs molex."""

import os
import sys
import subprocess
import math

import numpy as np
import gemmi

from cctbx import crystal, miller, sgtbx, xray
from cctbx.array_family import flex
import iotbx.cif
import iotbx.pdb

EXAMPLES_DIR = os.path.dirname(os.path.abspath(__file__))
MOLEX_BIN = os.path.join(
    EXAMPLES_DIR, "..", "target", "release", "examples", "dump_intermediates"
)

PDB_ID = "1AKI"
SG_NUMBER = 19
SG_NAME = "P 21 21 21"
COORD_PATH = os.path.join(EXAMPLES_DIR, f"{PDB_ID}.cif")
SF_PATH = os.path.join(EXAMPLES_DIR, f"{PDB_ID}-sf.cif")


def load_reflections_gemmi(sf_path):
    doc = gemmi.cif.read(sf_path)
    block = doc[0]
    h = np.array([int(x) for x in block.find_values("_refln.index_h")])
    k = np.array([int(x) for x in block.find_values("_refln.index_k")])
    l_arr = np.array([int(x) for x in block.find_values("_refln.index_l")])
    fobs = None
    for col in ["_refln.F_meas_au", "_refln.F_obs"]:
        try:
            vals = list(block.find_values(col))
            fobs = np.array([float(x) if x not in ("?", ".") else 0.0 for x in vals])
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
        if len(flags) == len(h):
            free = np.array([str(x) == "1" for x in flags])
    except Exception:
        pass
    if np.sum(free) == 0:
        try:
            status = list(block.find_values("_refln.status"))
            if len(status) == len(h):
                free = np.array([str(x).strip().lower() == "f" for x in status])
        except Exception:
            pass
    valid = fobs > 0
    if np.sum(free[valid]) == 0:
        rng = np.random.default_rng(42)
        free_v = rng.random(int(np.sum(valid))) < 0.05
    else:
        free_v = free[valid]
    return h[valid], k[valid], l_arr[valid], fobs[valid], sigma[valid], free_v


def compute_fc_cctbx(coord_path, h, k, l):
    pdb_inp = iotbx.pdb.input(file_name=coord_path)
    xrs = pdb_inp.xray_structure_simple()
    sel = xrs.element_selection("H")
    xrs = xrs.select(~sel)
    print(f"  cctbx: {xrs.scatterers().size()} atoms (H removed)")
    print(f"  cctbx: unit cell = {xrs.unit_cell()}")
    print(f"  cctbx: space group = {xrs.space_group_info()}")
    sym = xrs.crystal_symmetry()
    indices = flex.miller_index([(int(h[i]), int(k[i]), int(l[i])) for i in range(len(h))])
    ms = miller.set(crystal_symmetry=sym, indices=indices, anomalous_flag=False)
    # Use the xray.structure_factors.from_scatterers_direct API
    from cctbx.xray import structure_factors as sf_module
    fc_result = sf_module.from_scatterers_direct(
        xray_structure=xrs,
        miller_set=ms,
    )
    fc_data = fc_result.f_calc().data()
    fc_complex = np.array([(complex(c).real, complex(c).imag) for c in fc_data])
    return fc_complex[:, 0] + 1j * fc_complex[:, 1]


def round_up_to_smooth(n):
    n = int(math.ceil(n))
    while True:
        m = n
        for p in [2, 3, 5, 7]:
            while m % p == 0:
                m //= p
        if m == 1:
            return n
        n += 1


def molex_grid_dims(cell, d_min, sg_number):
    grid_factors_map = {
        1: (1,1,1), 4: (1,2,1), 5: (2,2,1), 19: (2,2,2),
        23: (2,2,2), 92: (2,4,2), 96: (2,4,2),
        152: (2,2,6), 154: (2,2,6), 178: (2,2,6),
    }
    gf = grid_factors_map.get(sg_number, (1,1,1))
    grid_spacing = d_min / 3.0
    nu = round_up_to_smooth(math.ceil(cell.a / grid_spacing))
    nv = round_up_to_smooth(math.ceil(cell.b / grid_spacing))
    nw = round_up_to_smooth(math.ceil(cell.c / grid_spacing))
    while nu % gf[0] != 0: nu = round_up_to_smooth(nu + 1)
    while nv % gf[1] != 0: nv = round_up_to_smooth(nv + 1)
    while nw % gf[2] != 0: nw = round_up_to_smooth(nw + 1)
    return nu, nv, nw


def compute_fc_gemmi(coord_path, h, k, l, cell, sg_number, grid_dims=None):
    st = gemmi.read_structure(coord_path)
    st.setup_entities()
    st.remove_hydrogens()
    s2_all = np.array([cell.calculate_1_d2([int(h[i]), int(k[i]), int(l[i])])
                       for i in range(len(h))])
    d_min = 1.0 / np.sqrt(s2_all.max())
    if grid_dims is None:
        nu, nv, nw = molex_grid_dims(cell, d_min, sg_number)
    else:
        nu, nv, nw = grid_dims
    print(f"  gemmi: grid dims = ({nu}, {nv}, {nw})")
    print(f"  gemmi: d_min = {d_min:.4f} A")
    dc = gemmi.DensityCalculatorX()
    dc.d_min = d_min
    dc.rate = 1.5
    dc.set_grid_cell_and_spacegroup(st)
    dc.grid.set_size(nu, nv, nw)
    dc.put_model_density_on_grid(st[0])
    arr = np.array(dc.grid, copy=False)
    fc_full = np.fft.fftn(arr)
    V = cell.volume
    N = arr.size
    fc_full *= V / N
    fc_vals = np.zeros(len(h), dtype=complex)
    for i in range(len(h)):
        u = int(h[i]) % nu
        v = int(k[i]) % nv
        w = int(l[i]) % nw
        fc_vals[i] = fc_full[u, v, w]
    return fc_vals, (nu, nv, nw)


def run_molex(pdb, sg):
    coord = os.path.join(EXAMPLES_DIR, f"{pdb}.cif")
    sf = os.path.join(EXAMPLES_DIR, f"{pdb}-sf.cif")
    if not os.path.exists(MOLEX_BIN):
        return {"error": "binary not found at " + MOLEX_BIN}
    try:
        result = subprocess.run(
            [MOLEX_BIN, coord, sf, str(sg)],
            capture_output=True, text=True, timeout=120,
        )
    except Exception as e:
        return {"error": str(e)}
    if result.returncode != 0:
        return {"error": result.stderr.strip()[:200]}
    vals = {}
    for line in result.stdout.strip().split("\n"):
        if "=" in line:
            key, val = line.split("=", 1)
            vals[key] = val
    return vals


def main():
    print("=" * 78)
    print("  Fc CROSS-VALIDATION: cctbx (analytical) vs gemmi/molex (grid FFT)")
    print(f"  Structure: {PDB_ID}  |  SG: {SG_NUMBER} ({SG_NAME})")
    print("=" * 78)

    h, k, l, fobs, sigma, free = load_reflections_gemmi(SF_PATH)
    n_refl = len(h)
    print(f"\n  Loaded {n_refl} reflections from SF file")

    st = gemmi.read_structure(COORD_PATH)
    cell = st.cell
    print(f"  Unit cell: a={cell.a:.3f} b={cell.b:.3f} c={cell.c:.3f}")
    print(f"             alpha={cell.alpha:.2f} beta={cell.beta:.2f} gamma={cell.gamma:.2f}")
    print(f"  Volume: {cell.volume:.2f} A^3")

    s2_all = np.array([cell.calculate_1_d2([int(h[i]), int(k[i]), int(l[i])])
                       for i in range(len(h))])
    d_min = 1.0 / np.sqrt(s2_all.max())
    d_array = 1.0 / np.sqrt(s2_all)
    print(f"  d_min = {d_min:.4f} A  |  d_max = {d_array.max():.2f} A")

    # --- cctbx ---
    print(f"\n{'~'*78}")
    print("  Computing Fc with cctbx (analytical direct summation)...")
    fc_cctbx = compute_fc_cctbx(COORD_PATH, h, k, l)
    fc_cctbx_amp = np.abs(fc_cctbx)
    fc_cctbx_phase = np.angle(fc_cctbx, deg=True)
    print(f"  cctbx: |Fc| mean={fc_cctbx_amp.mean():.2f}  median={np.median(fc_cctbx_amp):.2f}  "
          f"range=[{fc_cctbx_amp.min():.2f}, {fc_cctbx_amp.max():.2f}]")

    # --- gemmi ---
    print(f"\n{'~'*78}")
    print("  Computing Fc with gemmi (DensityCalculatorX + FFT + V/N)...")
    molex_dims = molex_grid_dims(cell, d_min, SG_NUMBER)
    fc_gemmi, gemmi_dims = compute_fc_gemmi(COORD_PATH, h, k, l, cell, SG_NUMBER, molex_dims)
    fc_gemmi_amp = np.abs(fc_gemmi)
    fc_gemmi_phase = np.angle(fc_gemmi, deg=True)
    print(f"  gemmi: |Fc| mean={fc_gemmi_amp.mean():.2f}  median={np.median(fc_gemmi_amp):.2f}  "
          f"range=[{fc_gemmi_amp.min():.2f}, {fc_gemmi_amp.max():.2f}]")

    # --- molex ---
    print(f"\n{'~'*78}")
    print("  Running molex dump_intermediates...")
    molex_out = run_molex(PDB_ID, SG_NUMBER)
    molex_grid = None
    fc_molex_proxy = None
    if "error" in molex_out:
        print(f"  molex ERROR: {molex_out['error']}")
    else:
        molex_grid_str = molex_out.get("MOLEX_GRID", "")
        dims_parts = molex_grid_str.split(",")
        molex_grid = tuple(int(x) for x in dims_parts)
        print(f"  molex: grid = {molex_grid}")
        print(f"  molex: n_atoms = {molex_out.get('MOLEX_N_ATOMS','?')}")
        print(f"  molex: n_reflections = {molex_out.get('MOLEX_N_REFLECTIONS','?')}")
        print(f"  molex: k_overall = {molex_out.get('MOLEX_K_OVERALL','?')}")
        print(f"  molex: R-work = {molex_out.get('MOLEX_R_WORK','?')}")
        print(f"  molex: R-free = {molex_out.get('MOLEX_R_FREE','?')}")

        if molex_grid != gemmi_dims:
            print(f"\n  Grid mismatch: gemmi computed={gemmi_dims} vs molex={molex_grid}")
            print("  Recomputing gemmi Fc with molex's exact grid...")
            fc_molex_proxy, _ = compute_fc_gemmi(COORD_PATH, h, k, l, cell, SG_NUMBER, molex_grid)
        else:
            print(f"  Grid match: both use {gemmi_dims}")
            fc_molex_proxy = fc_gemmi.copy()

    # ===================================================================
    print(f"\n{'='*78}")
    print("  COMPARISON RESULTS")
    print(f"{'='*78}")

    valid = (fc_cctbx_amp > 0.1) & (fc_gemmi_amp > 0.1)
    n_valid = int(np.sum(valid))
    print(f"\n  Using {n_valid}/{n_refl} reflections with |Fc| > 0.1")

    # 1. Correlation
    print(f"\n  {'~'*70}")
    print("  1. Fc AMPLITUDE CORRELATION (Pearson r)")
    print(f"  {'~'*70}")
    corr_cctbx_gemmi = np.corrcoef(fc_cctbx_amp[valid], fc_gemmi_amp[valid])[0, 1]
    print(f"  cctbx vs gemmi:       r = {corr_cctbx_gemmi:.8f}")
    corr_cctbx_molex = None
    corr_gemmi_molex = None
    if fc_molex_proxy is not None:
        fc_molex_amp = np.abs(fc_molex_proxy)
        valid_m = valid & (fc_molex_amp > 0.1)
        corr_cctbx_molex = np.corrcoef(fc_cctbx_amp[valid_m], fc_molex_amp[valid_m])[0, 1]
        corr_gemmi_molex = np.corrcoef(fc_gemmi_amp[valid_m], fc_molex_amp[valid_m])[0, 1]
        print(f"  cctbx vs molex-proxy: r = {corr_cctbx_molex:.8f}")
        print(f"  gemmi vs molex-proxy: r = {corr_gemmi_molex:.8f}")

    # 2. Ratio statistics
    print(f"\n  {'~'*70}")
    print("  2. Fc AMPLITUDE RATIO STATISTICS")
    print(f"  {'~'*70}")
    ratio_cg = fc_cctbx_amp[valid] / fc_gemmi_amp[valid]
    print(f"  |Fc_cctbx| / |Fc_gemmi|:")
    print(f"    mean   = {ratio_cg.mean():.6f}")
    print(f"    median = {np.median(ratio_cg):.6f}")
    print(f"    std    = {ratio_cg.std():.6f}")
    print(f"    min    = {ratio_cg.min():.6f}")
    print(f"    max    = {ratio_cg.max():.6f}")
    if fc_molex_proxy is not None:
        ratio_cm = fc_cctbx_amp[valid_m] / fc_molex_amp[valid_m]
        print(f"\n  |Fc_cctbx| / |Fc_molex_proxy|:")
        print(f"    mean   = {ratio_cm.mean():.6f}")
        print(f"    median = {np.median(ratio_cm):.6f}")
        print(f"    std    = {ratio_cm.std():.6f}")

    # 3. Systematic scale
    print(f"\n  {'~'*70}")
    print("  3. SYSTEMATIC SCALE ANALYSIS")
    print(f"  {'~'*70}")
    k_ls = np.sum(fc_cctbx_amp[valid] * fc_gemmi_amp[valid]) / np.sum(fc_gemmi_amp[valid]**2)
    print(f"  LS scale factor (cctbx/gemmi): k = {k_ls:.6f}")
    scaled_gemmi = k_ls * fc_gemmi_amp[valid]
    resid = np.abs(fc_cctbx_amp[valid] - scaled_gemmi) / fc_cctbx_amp[valid]
    print(f"  After LS scale correction:")
    print(f"    Mean relative error:   {resid.mean()*100:.4f}%")
    print(f"    Median relative error: {np.median(resid)*100:.4f}%")
    print(f"    Max relative error:    {resid.max()*100:.4f}%")

    # Resolution-dependent
    print(f"\n  Resolution-dependent grid sampling error (after LS scaling):")
    print(f"  {'Res bin (A)':<18} {'N':>6} {'Mean %err':>10} {'Med %err':>10} {'Max %err':>10}")
    print(f"  {'-'*58}")
    d_vals = d_array[valid]
    resid_pct = resid * 100
    d_edges = [999.0, 5.0, 3.0, 2.5, 2.0, 1.8, 1.6, d_min - 0.01]
    for i in range(len(d_edges) - 1):
        hi, lo = d_edges[i], d_edges[i + 1]
        mask = (d_vals <= hi) & (d_vals > lo)
        n_bin = int(np.sum(mask))
        if n_bin > 0:
            me = resid_pct[mask].mean()
            md = np.median(resid_pct[mask])
            mx = resid_pct[mask].max()
            print(f"  {lo:5.2f} - {hi:5.2f}       {n_bin:>6}   {me:>8.4f}%  {md:>8.4f}%  {mx:>8.4f}%")

    # 4. Phase agreement
    # NOTE: cctbx uses exp(+2pi i h.r) convention while gemmi/numpy FFT uses
    # exp(-2pi i h.r). This means Fc_cctbx = conj(Fc_gemmi). We conjugate
    # the gemmi phases to match cctbx convention before comparing.
    print(f"\n  {'~'*70}")
    print("  4. PHASE AGREEMENT (cctbx vs gemmi)")
    print(f"     [gemmi phases conjugated to match cctbx exp(+2pi i h.r) convention]")
    print(f"  {'~'*70}")
    fc_gemmi_phase_conj = -fc_gemmi_phase  # conjugate = negate imaginary = negate phase
    phase_diff = fc_cctbx_phase[valid] - fc_gemmi_phase_conj[valid]
    phase_diff = (phase_diff + 180) % 360 - 180
    weights = fc_cctbx_amp[valid]
    weighted_phase_err = np.average(np.abs(phase_diff), weights=weights)
    print(f"  Phase difference (degrees):")
    print(f"    Mean |delta_phi|:            {np.mean(np.abs(phase_diff)):.4f}")
    print(f"    Amplitude-weighted |dphi|:   {weighted_phase_err:.4f}")
    print(f"    Median |delta_phi|:          {np.median(np.abs(phase_diff)):.4f}")
    print(f"    Std delta_phi:               {phase_diff.std():.4f}")
    print(f"    Max |delta_phi|:             {np.max(np.abs(phase_diff)):.4f}")

    print(f"\n  Phase error by resolution:")
    print(f"  {'Res bin (A)':<18} {'N':>6} {'Mean |dphi|':>12} {'Wtd |dphi|':>12}")
    print(f"  {'-'*52}")
    for i in range(len(d_edges) - 1):
        hi, lo = d_edges[i], d_edges[i + 1]
        mask = (d_vals <= hi) & (d_vals > lo)
        n_bin = int(np.sum(mask))
        if n_bin > 0:
            me = np.mean(np.abs(phase_diff[mask]))
            wtd = np.average(np.abs(phase_diff[mask]), weights=weights[mask])
            print(f"  {lo:5.2f} - {hi:5.2f}       {n_bin:>6}   {me:>10.4f}    {wtd:>10.4f}")

    # 5. R-factors
    print(f"\n  {'~'*70}")
    print("  5. R-FACTOR FROM RAW Fc (no bulk solvent, simple Wilson k scaling)")
    print(f"  {'~'*70}")
    work = ~free
    test = free
    for label, fc_amp in [("cctbx", fc_cctbx_amp), ("gemmi", fc_gemmi_amp)]:
        good = work & (fobs > 1) & (fc_amp > 0.01)
        if np.sum(good) > 20:
            k_wilson = float(np.median(fobs[good] / fc_amp[good]))
        else:
            k_wilson = 1.0
        r_work = float(np.sum(np.abs(fobs[work] - k_wilson * fc_amp[work])) / np.sum(fobs[work]))
        r_free_val = float(np.sum(np.abs(fobs[test] - k_wilson * fc_amp[test])) / np.sum(fobs[test]))
        print(f"  {label:<7}: k_wilson={k_wilson:.4f}  R-work={r_work:.4f}  R-free={r_free_val:.4f}")

    # ===================================================================
    print(f"\n{'='*78}")
    print("  SUMMARY TABLE")
    print(f"{'='*78}")
    print()
    print(f"  {'Metric':<48} {'Value':>15}")
    print(f"  {'-'*65}")
    print(f"  {'Structure':<48} {PDB_ID + ' (SG ' + str(SG_NUMBER) + ')':>15}")
    print(f"  {'Number of reflections':<48} {n_refl:>15}")
    print(f"  {'d_min (A)':<48} {d_min:>15.4f}")
    print(f"  {'Grid dimensions':<48} {str(gemmi_dims):>15}")
    print()
    print(f"  {'Pearson r: cctbx vs gemmi |Fc|':<48} {corr_cctbx_gemmi:>15.8f}")
    if corr_cctbx_molex is not None:
        print(f"  {'Pearson r: cctbx vs molex-proxy |Fc|':<48} {corr_cctbx_molex:>15.8f}")
        print(f"  {'Pearson r: gemmi vs molex-proxy |Fc|':<48} {corr_gemmi_molex:>15.8f}")
    print()
    print(f"  {'|Fc_cctbx|/|Fc_gemmi| mean':<48} {ratio_cg.mean():>15.6f}")
    print(f"  {'|Fc_cctbx|/|Fc_gemmi| std':<48} {ratio_cg.std():>15.6f}")
    print(f"  {'LS scale factor k (cctbx = k * gemmi)':<48} {k_ls:>15.6f}")
    print()
    print(f"  {'Mean relative |Fc| error (LS-scaled)':<48} {resid.mean()*100:>14.4f}%")
    print(f"  {'Median relative |Fc| error (LS-scaled)':<48} {np.median(resid)*100:>14.4f}%")
    print(f"  {'Max relative |Fc| error (LS-scaled)':<48} {resid.max()*100:>14.4f}%")
    print()
    print(f"  {'Mean |phase diff| cctbx vs gemmi (deg)':<48} {np.mean(np.abs(phase_diff)):>15.4f}")
    print(f"  {'Amplitude-weighted |phase diff| (deg)':<48} {weighted_phase_err:>15.4f}")
    print(f"  {'Max |phase diff| (deg)':<48} {np.max(np.abs(phase_diff)):>15.4f}")
    print()
    if molex_grid is not None and "error" not in molex_out:
        print(f"  {'molex grid dims':<48} {molex_out.get('MOLEX_GRID',''):>15}")
        print(f"  {'molex R-work (with bulk solvent + scaling)':<48} {molex_out.get('MOLEX_R_WORK',''):>15}")
        print(f"  {'molex R-free (with bulk solvent + scaling)':<48} {molex_out.get('MOLEX_R_FREE',''):>15}")
    print()
    print(f"  {'-'*65}")
    print(f"  The grid FFT method (gemmi/molex) introduces ~{resid.mean()*100:.2f}% mean")
    print(f"  relative |Fc| error vs cctbx analytical. Phase agreement is")
    print(f"  ~{np.mean(np.abs(phase_diff)):.2f} deg mean. Error grows at high resolution")
    print(f"  (near d_min) due to grid sampling limitations.")
    print(f"{'='*78}")


if __name__ == "__main__":
    main()

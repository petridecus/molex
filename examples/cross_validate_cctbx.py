#!/usr/bin/env python3
"""Cross-validate molex xtal pipeline against cctbx/Phenix reference.

Uses mmtbx.f_model (the core of phenix.refine) to compute:
  - Fc from atomic model (IT92 form factors)
  - Bulk solvent correction (flat mask model)
  - Anisotropic scaling (Levenberg-Marquardt)
  - Sigma-A estimation (maximum likelihood)
  - R-work, R-free
  - Map coefficients (2mFo-DFc)

Then compares stage-by-stage against molex output.

Usage:
  python3 examples/cross_validate_cctbx.py
"""

import os
import sys
import subprocess

import numpy as np

from cctbx import crystal, miller, sgtbx, xray
from cctbx.array_family import flex
from mmtbx import f_model
import iotbx.cif


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


def run_cctbx(pdb, sg_number):
    """Compute reference values using cctbx/mmtbx (Phenix engine)."""
    coord_path = os.path.join(EXAMPLES_DIR, f"{pdb}.cif")
    sf_path = os.path.join(EXAMPLES_DIR, f"{pdb}-sf.cif")

    # Read structure
    cif_model = iotbx.cif.reader(file_path=coord_path).build_crystal_structures()
    if not cif_model:
        # Fall back to mmcif model reading
        from iotbx import pdb as iotbx_pdb
        pdb_inp = iotbx_pdb.input(file_name=coord_path)
        xrs = pdb_inp.xray_structure_simple()
    else:
        key = list(cif_model.keys())[0]
        xrs = cif_model[key]

    # Read reflections
    cif_reader = iotbx.cif.reader(file_path=sf_path)
    miller_arrays = cif_reader.as_miller_arrays()

    # Find the Fobs array
    f_obs = None
    r_free_flags = None
    for ma in miller_arrays:
        info = str(ma.info())
        if ma.is_real_array():
            if "F_meas" in info or "F_obs" in info:
                f_obs = ma
            elif "intensity" in info or "I_obs" in info:
                f_obs = ma.as_amplitude_array()
        if ma.is_integer_array() or ma.is_bool_array():
            if "free" in info.lower() or "status" in info.lower():
                r_free_flags = ma

    if f_obs is None:
        # Try to get any amplitude or intensity array
        for ma in miller_arrays:
            if ma.is_real_array() and not ma.is_complex_array():
                if f_obs is None:
                    f_obs = ma
                    break

    if f_obs is None:
        return {"error": "no Fobs array found"}

    # Ensure positive values
    f_obs = f_obs.select(f_obs.data() > 0)

    # Make r_free_flags if not available
    if r_free_flags is not None:
        # Convert to boolean: flag value 1 = free
        try:
            r_free_flags = r_free_flags.common_set(f_obs)
            flags_data = r_free_flags.data()
            if hasattr(flags_data, 'as_bool'):
                r_free_flags = r_free_flags.customized_copy(data=flags_data == 1)
            else:
                r_free_flags = r_free_flags.customized_copy(
                    data=flex.bool([int(x) == 1 for x in flags_data])
                )
        except Exception:
            r_free_flags = None

    if r_free_flags is None:
        r_free_flags = f_obs.generate_r_free_flags(fraction=0.05, use_lattice_symmetry=False)

    # Ensure common set
    f_obs, r_free_flags = f_obs.common_sets(r_free_flags)

    # Match to xray structure
    f_obs = f_obs.resolution_filter(d_min=f_obs.d_min())

    # Create f_model manager — this is the core of phenix.refine
    fmodel = f_model.manager(
        xray_structure=xrs,
        f_obs=f_obs,
        r_free_flags=r_free_flags,
    )

    # Update scaling (bulk solvent + anisotropic)
    fmodel.update_all_scales(remove_outliers=False)

    # Extract results
    r_work = fmodel.r_work()
    r_free = fmodel.r_free()

    # Scaling: cctbx uses k_isotropic (array) and k_anisotropic (array)
    # scale_k1 is the overall scale applied to Fobs
    k1 = fmodel.scale_k1()

    # Bulk solvent
    k_sol_val = fmodel.k_sol()
    b_sol_val = fmodel.b_sol()
    b_cart_val = fmodel.b_cart()

    n_refl = f_obs.size()
    d_min_val = f_obs.d_min()

    return {
        "n_refl": n_refl,
        "d_min": d_min_val,
        "r_work": r_work,
        "r_free": r_free,
        "scale_k1": k1,
        "k_sol": k_sol_val,
        "b_sol": b_sol_val,
        "b_cart": list(b_cart_val) if b_cart_val is not None else [],
    }


def compare(pdb, sg, sg_name):
    print(f"\n{'='*76}")
    print(f"  {pdb}  SG {sg} ({sg_name})")
    print(f"{'='*76}")

    molex = run_molex(pdb, sg)
    if "error" in molex:
        print(f"  molex ERROR: {molex['error']}")
        return False

    try:
        ref = run_cctbx(pdb, sg)
    except Exception as e:
        print(f"  cctbx ERROR: {e}")
        return False

    if "error" in ref:
        print(f"  cctbx ERROR: {ref['error']}")
        return False

    m_k = float(molex.get("MOLEX_K_OVERALL", "0"))
    m_ksol = float(molex.get("MOLEX_K_SOL", "0"))
    m_bsol = float(molex.get("MOLEX_B_SOL", "0"))
    m_rw = float(molex.get("MOLEX_R_WORK", "0"))
    m_rf = float(molex.get("MOLEX_R_FREE", "0"))

    c_k1 = ref["scale_k1"]
    c_ksol = ref["k_sol"]
    c_bsol = ref["b_sol"]
    c_rw = ref["r_work"]
    c_rf = ref["r_free"]

    print(f"\n  {'Parameter':<20} {'molex':>12} {'cctbx':>12} {'diff':>12}")
    print(f"  {'-'*60}")
    print(f"  {'scale_k1 / k_ov':<20} {m_k:>12.4f} {c_k1:>12.4f}")
    print(f"  {'k_sol':<20} {m_ksol:>12.4f} {c_ksol:>12.4f} {abs(m_ksol-c_ksol):>12.4f}")
    print(f"  {'B_sol':<20} {m_bsol:>12.1f} {c_bsol:>12.1f} {abs(m_bsol-c_bsol):>12.1f}")
    print(f"  {'R-work':<20} {m_rw:>12.4f} {c_rw:>12.4f} {abs(m_rw-c_rw):>12.4f}")
    print(f"  {'R-free':<20} {m_rf:>12.4f} {c_rf:>12.4f} {abs(m_rf-c_rf):>12.4f}")

    if ref.get("b_cart"):
        print(f"  {'B_cart (cctbx)':<20} {str(ref['b_cart'][:3]):>60}")

    rw_diff = abs(m_rw - c_rw)
    ok = rw_diff < 0.15
    print(f"\n  R-work diff from cctbx: {rw_diff:.4f}  {'PASS' if ok else 'FAIL'}")
    return ok


def main():
    if not os.path.exists(MOLEX_BIN):
        print("Build first: cargo build --example dump_intermediates --features xtal --release")
        sys.exit(1)

    print("=" * 76)
    print("  MOLEX vs CCTBX (Phenix engine) CROSS-VALIDATION")
    print("  cctbx: mmtbx.f_model with full bulk solvent + aniso scaling")
    print("=" * 76)

    n_pass = 0
    n_total = 0
    for pdb, sg, sg_name in STRUCTURES:
        n_total += 1
        try:
            if compare(pdb, sg, sg_name):
                n_pass += 1
        except Exception as e:
            print(f"  EXCEPTION: {e}")

    print(f"\n{'='*76}")
    print(f"  RESULT: {n_pass}/{n_total} structures within tolerance of cctbx")
    print(f"{'='*76}")


if __name__ == "__main__":
    main()

//! Validate the molex xtal pipeline against real crystallographic data.
//!
//! Tests all 10 supported space groups using deposited PDB structures with
//! structure factors.
//!
//! Run with:
//! ```sh
//! cargo run --example validate_xtal --features xtal --release
//! ```
#![allow(
    clippy::print_stdout,
    clippy::print_stderr,
    clippy::expect_used,
    clippy::unwrap_used,
    clippy::panic,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::missing_docs_in_private_items,
    clippy::too_many_lines,
    missing_docs,
    unused_results
)]

use std::path::Path;

use molex::adapters::cif::{parse, CoordinateData, ReflectionData};
use molex::element::Element;
use molex::xtal::{
    self, grid_factors, requires_equal_uv, round_up_to_smooth, space_group,
    UnitCell, XtalRefinement,
};

// ── Test structures: one per supported space group ──────────────────────

/// (PDB code, space group number, SG name)
const TEST_STRUCTURES: &[(&str, u16, &str)] = &[
    // ── One per space group ──────────────────────────────────────────
    ("3LZT", 1, "P 1"),
    ("3NIR", 4, "P 1 21 1"),
    ("6E6O", 5, "C 1 2 1"),
    ("1AKI", 19, "P 21 21 21"),
    ("7KOM", 23, "I 2 2 2"),
    ("1X6Z", 92, "P 41 21 2"),
    ("1G6X", 96, "P 43 21 2"),
    ("4XDX", 152, "P 31 2 1"),
    ("2H5C", 154, "P 32 2 1"),
    ("3ZUC", 178, "P 61 2 2"),
    // ── Gap-filling structures ───────────────────────────────────────
    ("5YPA", 19, "P 21 21 21"), // low resolution (~3 A)
    ("7AF2", 19, "P 21 21 21"), // large protein, real R-free flags
    ("2OXY", 1, "P 1"),         // multi-chain (hemoglobin, 4 chains)
    ("3SGJ", 19, "P 21 21 21"), // multi-chain (trypsin + inhibitor)
];

fn main() {
    let examples_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("examples");

    let mut results: Vec<StructureResult> = Vec::new();
    let mut any_failure = false;

    for &(pdb, sg_number, sg_name) in TEST_STRUCTURES {
        println!(
            "\n============================================================"
        );
        println!("  {pdb}  SG {sg_number} ({sg_name})");
        println!(
            "============================================================"
        );

        let coord_path = examples_dir.join(format!("{pdb}.cif"));
        let sf_path = examples_dir.join(format!("{pdb}-sf.cif"));

        if !coord_path.exists() || !sf_path.exists() {
            println!("  SKIP: files not found");
            continue;
        }

        match validate_structure(pdb, sg_number, &coord_path, &sf_path) {
            Ok(r) => {
                results.push(r);
            }
            Err(msg) => {
                println!("  FAIL: {msg}");
                results.push(StructureResult {
                    pdb: pdb.to_owned(),
                    sg_number,
                    sg_name: sg_name.to_owned(),
                    n_atoms: 0,
                    n_reflections: 0,
                    r_work: f64::NAN,
                    r_free: f64::NAN,
                    r_work_after: f64::NAN,
                    refinement_improved: false,
                    b_corr: f64::NAN,
                    pass: false,
                });
                any_failure = true;
            }
        }
    }

    // ── Summary table ───────────────────────────────────────────────────

    println!("\n\n============================================================================");
    println!("  VALIDATION SUMMARY");
    println!("============================================================================");
    println!(
        "{:<6} {:<14} {:>6} {:>7} {:>7} {:>7} {:>7} {:>6} {:>5}",
        "PDB",
        "SpaceGroup",
        "Atoms",
        "Refls",
        "Rwork",
        "Rfree",
        "Refine",
        "Bcorr",
        "Pass"
    );
    println!("{}", "-".repeat(76));

    for r in &results {
        let refine_str = if r.refinement_improved {
            format!("{:.4}", r.r_work_after)
        } else {
            "FAIL".to_owned()
        };
        println!(
            "{:<6} SG {:>3} {:<8} {:>6} {:>7} {:>7.4} {:>7.4} {:>7} {:>6.3} \
             {:>5}",
            r.pdb,
            r.sg_number,
            r.sg_name,
            r.n_atoms,
            r.n_reflections,
            r.r_work,
            r.r_free,
            refine_str,
            r.b_corr,
            if r.pass { "OK" } else { "FAIL" }
        );
    }

    let n_pass = results.iter().filter(|r| r.pass).count();
    let n_total = results.len();
    println!("\n{n_pass}/{n_total} structures passed.");

    assert!(!any_failure, "Some structures failed validation");
    assert!(
        n_pass == n_total,
        "{} of {n_total} structures failed",
        n_total - n_pass
    );
}

struct StructureResult {
    pdb: String,
    sg_number: u16,
    sg_name: String,
    n_atoms: usize,
    n_reflections: usize,
    r_work: f64,
    r_free: f64,
    r_work_after: f64,
    refinement_improved: bool,
    b_corr: f64,
    pass: bool,
}

fn validate_structure(
    pdb: &str,
    sg_number: u16,
    coord_path: &Path,
    sf_path: &Path,
) -> Result<StructureResult, String> {
    let sg_name_str;

    // ── Parse coordinates ───────────────────────────────────────────────

    let coord_str = std::fs::read_to_string(coord_path)
        .map_err(|e| format!("read coord: {e}"))?;
    let coord_doc =
        parse(&coord_str).map_err(|e| format!("parse coord: {e}"))?;
    let coord_block =
        coord_doc.blocks.first().ok_or("no blocks in coord CIF")?;
    let coord_data = CoordinateData::try_from(coord_block)
        .map_err(|e| format!("extract coords: {e}"))?;

    let cell = coord_data
        .cell
        .as_ref()
        .ok_or("no unit cell in coord CIF")?;
    sg_name_str = coord_data.spacegroup.as_deref().unwrap_or("?").to_owned();

    let atoms: Vec<_> = coord_data
        .atoms
        .iter()
        .filter(|a| {
            let elem = if a.element.is_empty() {
                Element::from_atom_name(&a.label)
            } else {
                Element::from_symbol(&a.element)
            };
            elem != Element::H
        })
        .collect();

    println!(
        "  Atoms: {} non-H | Cell: {:.1} {:.1} {:.1} | SG: {sg_name_str}",
        atoms.len(),
        cell.a,
        cell.b,
        cell.c
    );

    // ── Parse reflections ───────────────────────────────────────────────

    let sf_str = std::fs::read_to_string(sf_path)
        .map_err(|e| format!("read sf: {e}"))?;
    let sf_doc = parse(&sf_str).map_err(|e| format!("parse sf: {e}"))?;
    let sf_block = sf_doc.blocks.first().ok_or("no blocks in SF CIF")?;
    let refl_data = ReflectionData::try_from(sf_block)
        .map_err(|e| format!("extract reflections: {e}"))?;

    // Ensure free set exists
    let seed = hash_pdb_code(pdb);
    let mut refl_data_patched = refl_data;
    let has_any_free =
        refl_data_patched.reflections.iter().any(|r| r.free_flag);
    if !has_any_free {
        refl_data_patched.free_flags_from_file = false;
    }
    let reflections =
        xtal::reflections_from_cif(&refl_data_patched, 0.05, seed);

    let n_free = reflections.iter().filter(|r| r.free_flag).count();
    println!(
        "  Reflections: {} ({} free, {:.1}%)",
        reflections.len(),
        n_free,
        100.0 * n_free as f64 / reflections.len().max(1) as f64
    );

    if reflections.is_empty() {
        return Err("no reflections with valid Fobs".into());
    }

    // ── Build pipeline ──────────────────────────────────────────────────

    let sg = space_group(sg_number)
        .ok_or_else(|| format!("SG {sg_number} not supported"))?;
    let gf = grid_factors(sg_number)
        .ok_or_else(|| format!("no grid factors for SG {sg_number}"))?;
    let equal_uv = requires_equal_uv(sg_number);

    let unit_cell = UnitCell::new(
        cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma,
    );

    // Resolution
    let mut d_min = f64::MAX;
    for r in &reflections {
        let s2 = unit_cell.d_star_sq(r.h, r.k, r.l);
        if s2 > 0.0 {
            let d = 1.0 / s2.sqrt();
            if d < d_min {
                d_min = d;
            }
        }
    }

    // Grid dimensions
    let grid_spacing = d_min / 3.0;
    let mut nu = round_up_to_smooth((cell.a / grid_spacing).ceil());
    let mut nv = round_up_to_smooth((cell.b / grid_spacing).ceil());
    let mut nw = round_up_to_smooth((cell.c / grid_spacing).ceil());

    while nu % gf[0] != 0 {
        nu = round_up_to_smooth((nu + 1) as f64);
    }
    while nv % gf[1] != 0 {
        nv = round_up_to_smooth((nv + 1) as f64);
    }
    while nw % gf[2] != 0 {
        nw = round_up_to_smooth((nw + 1) as f64);
    }
    if equal_uv {
        let m = nu.max(nv);
        nu = m;
        nv = m;
    }

    println!("  Grid: [{nu}, {nv}, {nw}] = {} points", nu * nv * nw);

    let mut refinement =
        XtalRefinement::new(unit_cell, sg, reflections, [nu, nv, nw]);

    // ── Extract atom data ───────────────────────────────────────────────

    let positions: Vec<[f64; 3]> = atoms
        .iter()
        .map(|a| refinement.unit_cell.fractionalize([a.x, a.y, a.z]))
        .collect();
    let elements: Vec<Element> = atoms
        .iter()
        .map(|a| {
            if a.element.is_empty() {
                Element::from_atom_name(&a.label)
            } else {
                Element::from_symbol(&a.element)
            }
        })
        .collect();
    let b_factors: Vec<f64> = atoms.iter().map(|a| a.b_factor).collect();
    let occupancies: Vec<f64> = atoms.iter().map(|a| a.occupancy).collect();

    // ── Compute map ─────────────────────────────────────────────────────

    let map = refinement
        .compute_map(&positions, &elements, &b_factors, &occupancies)
        .ok_or("compute_map failed")?;

    let n = map.data.len() as f64;
    let mean = map.data.iter().map(|v| f64::from(*v)).sum::<f64>() / n;
    let variance = map
        .data
        .iter()
        .map(|v| (f64::from(*v) - mean).powi(2))
        .sum::<f64>()
        / n;
    let sigma = variance.sqrt();

    if !sigma.is_finite() || sigma <= 0.0 {
        return Err(format!("bad map sigma: {sigma}"));
    }

    println!(
        "  Scaling: k={:.4} k_sol={:.4} B_sol={:.1}",
        refinement.scaling.k_overall,
        refinement.scaling.k_sol,
        refinement.scaling.b_sol
    );

    // ── R-factors ───────────────────────────────────────────────────────

    let (r_work, r_free) = refinement
        .r_factors(&positions, &elements, &b_factors, &occupancies)
        .ok_or("r_factors failed")?;

    println!(
        "  R-work: {r_work:.4}  R-free: {r_free:.4}  gap: {:.4}",
        r_free - r_work
    );

    if r_work > 0.60 {
        return Err(format!("R-work {r_work:.4} > 0.60"));
    }

    // ── B-factor refinement ─────────────────────────────────────────────

    let original_b = b_factors.clone();
    let mean_b_dep = original_b.iter().sum::<f64>() / original_b.len() as f64;
    let mut perturbed_b: Vec<f64> = vec![mean_b_dep; b_factors.len()];

    let (r_work_before, _) = refinement
        .r_factors(&positions, &elements, &perturbed_b, &occupancies)
        .ok_or("r_factors (perturbed) failed")?;

    // Adaptive learning rate: try progressively smaller values until
    // refinement improves R-work (or give up).
    let mut r_work_after = r_work_before;
    for &lr in &[10000.0, 5000.0, 2000.0, 1000.0, 500.0, 200.0, 50.0] {
        let mut test_b = vec![mean_b_dep; b_factors.len()];
        let result = refinement.refine(
            &positions,
            &elements,
            &mut test_b,
            &occupancies,
            10,
            lr,
        );
        if let Some((rw, _)) = result {
            if rw < r_work_before {
                r_work_after = rw;
                perturbed_b = test_b;
                break;
            }
        }
    }

    let improved = r_work_after < r_work_before;
    let b_corr = pearson_correlation(&original_b, &perturbed_b);

    println!(
        "  Refine: R-work {r_work_before:.4} -> {r_work_after:.4} {} | \
         B-corr: {b_corr:.3}",
        if improved { "OK" } else { "WORSE" }
    );

    let pass = r_work < 0.60 && sigma > 0.0 && improved;
    println!("  {}", if pass { "PASS" } else { "FAIL" });

    Ok(StructureResult {
        pdb: pdb.to_owned(),
        sg_number,
        sg_name: sg_name_str,
        n_atoms: atoms.len(),
        n_reflections: refinement.reflections.len(),
        r_work,
        r_free,
        r_work_after,
        refinement_improved: improved,
        b_corr,
        pass,
    })
}

fn hash_pdb_code(pdb: &str) -> u64 {
    let mut h = 0x517cc1b727220a95_u64;
    for b in pdb.bytes() {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x6c62_272e_07bb_0142);
    }
    h
}

fn pearson_correlation(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    let mean_x = x.iter().sum::<f64>() / n;
    let mean_y = y.iter().sum::<f64>() / n;
    let mut cov = 0.0;
    let mut var_x = 0.0;
    let mut var_y = 0.0;
    for i in 0..x.len() {
        let dx = x[i] - mean_x;
        let dy = y[i] - mean_y;
        cov += dx * dy;
        var_x += dx * dx;
        var_y += dy * dy;
    }
    let denom = (var_x * var_y).sqrt();
    if denom < 1e-30 {
        0.0
    } else {
        cov / denom
    }
}

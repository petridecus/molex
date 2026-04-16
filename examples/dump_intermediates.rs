//! Dump molex xtal pipeline intermediate values for cross-validation.
//!
//! Outputs CSV-formatted Fc amplitudes and pipeline parameters that can be
//! compared with gemmi reference values.
//!
//! Usage:
//!   cargo run --example dump_intermediates --features xtal --release --
//! <coord.cif> <sf.cif> <sg_number>
#![allow(
    clippy::print_stdout,
    clippy::expect_used,
    clippy::unwrap_used,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::missing_docs_in_private_items,
    clippy::too_many_lines,
    missing_docs,
    unused_results
)]

use std::env;
use std::path::Path;

use molex::adapters::cif::{parse, CoordinateData, ReflectionData};
use molex::element::Element;
use molex::xtal::{
    self, grid_factors, requires_equal_uv, round_up_to_smooth, space_group,
    UnitCell, XtalRefinement,
};

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 4 {
        eprintln!("Usage: {} <coord.cif> <sf.cif> <sg_number>", args[0]);
        std::process::exit(1);
    }

    let coord_path = Path::new(&args[1]);
    let sf_path = Path::new(&args[2]);
    let sg_number: u16 = args[3].parse().expect("invalid SG number");

    // Parse coordinates
    let coord_str =
        std::fs::read_to_string(coord_path).expect("failed to read coord file");
    let coord_doc = parse(&coord_str).expect("failed to parse coord CIF");
    let coord_block = coord_doc.blocks.first().expect("no blocks");
    let coord_data = CoordinateData::try_from(coord_block)
        .expect("failed to extract coords");

    let cell = coord_data.cell.as_ref().expect("no unit cell");

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

    // Parse reflections
    let sf_str =
        std::fs::read_to_string(sf_path).expect("failed to read SF file");
    let sf_doc = parse(&sf_str).expect("failed to parse SF CIF");
    let sf_block = sf_doc.blocks.first().expect("no blocks");
    let refl_data = ReflectionData::try_from(sf_block)
        .expect("failed to extract reflections");

    let mut refl_patched = refl_data;
    let has_free = refl_patched.reflections.iter().any(|r| r.free_flag);
    if !has_free {
        refl_patched.free_flags_from_file = false;
    }
    let reflections = xtal::reflections_from_cif(&refl_patched, 0.05, 42);

    // Build pipeline
    let sg = space_group(sg_number).expect("unsupported SG");
    let gf = grid_factors(sg_number).expect("no grid factors");
    let equal_uv = requires_equal_uv(sg_number);

    let unit_cell = UnitCell::new(
        cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma,
    );

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

    let mut refinement =
        XtalRefinement::new(unit_cell, sg, reflections, [nu, nv, nw]);

    // Extract atom data
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

    // Compute map
    let _map = refinement
        .compute_map(&positions, &elements, &b_factors, &occupancies)
        .expect("compute_map failed");

    // R-factors
    let (r_work, r_free) = refinement
        .r_factors(&positions, &elements, &b_factors, &occupancies)
        .expect("r_factors failed");

    // Output summary as key=value
    println!("MOLEX_GRID={},{},{}", nu, nv, nw);
    println!("MOLEX_N_ATOMS={}", atoms.len());
    println!("MOLEX_N_REFLECTIONS={}", refinement.reflections.len());
    println!("MOLEX_K_OVERALL={:.6}", refinement.scaling.k_overall);
    println!("MOLEX_K_SOL={:.6}", refinement.scaling.k_sol);
    println!("MOLEX_B_SOL={:.4}", refinement.scaling.b_sol);
    println!("MOLEX_R_WORK={:.6}", r_work);
    println!("MOLEX_R_FREE={:.6}", r_free);

    if let Some(ref sa) = refinement.sigma_a {
        println!(
            "MOLEX_D_BINS={}",
            sa.d_bins
                .iter()
                .map(|d| format!("{d:.6}"))
                .collect::<Vec<_>>()
                .join(",")
        );
        println!(
            "MOLEX_SIGMA_SQ_BINS={}",
            sa.sigma_sq_bins
                .iter()
                .map(|s| format!("{s:.6}"))
                .collect::<Vec<_>>()
                .join(",")
        );
    }
}

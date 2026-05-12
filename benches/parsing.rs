//! Benchmarks for file format parsing.
#![allow(
    missing_docs,
    unused_results,
    clippy::unwrap_used,
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_possible_wrap,
    clippy::suboptimal_flops,
    clippy::format_push_string,
    clippy::uninlined_format_args
)]

use criterion::{
    black_box, criterion_group, criterion_main, BenchmarkId, Criterion,
};
use molex::adapters::cif::mmcif_str_to_entities;
use molex::adapters::pdb::pdb_str_to_entities;

fn generate_pdb(n_residues: usize) -> String {
    let mut pdb = String::new();
    let mut serial = 1usize;
    for res_num in 1..=n_residues {
        for &(name, dx, dy, dz, elem) in &[
            ("N  ", 0.0f64, 0.0, 0.0, "N"),
            ("CA ", 1.47, 0.0, 0.0, "C"),
            ("C  ", 2.45, 1.0, 0.0, "C"),
            ("O  ", 2.45, 2.2, 0.0, "O"),
            ("CB ", 1.47, -1.5, 0.0, "C"),
        ] {
            let x = dx + (res_num as f64) * 3.8;
            pdb.push_str(&format!(
                "ATOM  {:>5} {:<4}ALA A{:>4}    {:>8.3}{:>8.3}{:>8.3}  1.00  \
                 0.00          {:>2}  \n",
                serial, name, res_num, x, dy, dz, elem
            ));
            serial += 1;
        }
    }
    pdb.push_str("END\n");
    pdb
}

fn generate_mmcif(n_residues: usize) -> String {
    let mut cif = String::from(concat!(
        "data_test\n",
        "loop_\n",
        "_atom_site.group_PDB\n",
        "_atom_site.label_atom_id\n",
        "_atom_site.label_comp_id\n",
        "_atom_site.label_asym_id\n",
        "_atom_site.label_seq_id\n",
        "_atom_site.Cartn_x\n",
        "_atom_site.Cartn_y\n",
        "_atom_site.Cartn_z\n",
        "_atom_site.type_symbol\n",
        "_atom_site.occupancy\n",
        "_atom_site.B_iso_or_equiv\n",
    ));
    for res_num in 1..=n_residues {
        for &(name, elem, dx, dy, dz) in &[
            ("N", "N", 0.0f64, 0.0, 0.0),
            ("CA", "C", 1.47, 0.0, 0.0),
            ("C", "C", 2.45, 1.0, 0.0),
            ("O", "O", 2.45, 2.2, 0.0),
            ("CB", "C", 1.47, -1.5, 0.0),
        ] {
            let x = dx + (res_num as f64) * 3.8;
            cif.push_str(&format!(
                "ATOM {name} ALA A {res_num} {x:.3} {dy:.3} {dz:.3} {elem} \
                 1.00 0.00\n"
            ));
        }
    }
    cif
}

fn bench_pdb_parsing(c: &mut Criterion) {
    let mut group = c.benchmark_group("pdb_parse");
    for n_residues in [10, 50, 200, 1000] {
        let pdb = generate_pdb(n_residues);
        let n_atoms = n_residues * 5;
        group.bench_with_input(
            BenchmarkId::new("to_entities", format!("{n_atoms}_atoms")),
            &pdb,
            |b, pdb| b.iter(|| pdb_str_to_entities(black_box(pdb)).unwrap()),
        );
    }
    group.finish();
}

fn bench_mmcif_parsing(c: &mut Criterion) {
    let mut group = c.benchmark_group("mmcif_parse");
    for n_residues in [10, 50, 200, 1000] {
        let cif = generate_mmcif(n_residues);
        let n_atoms = n_residues * 5;
        group.bench_with_input(
            BenchmarkId::new("to_entities", format!("{n_atoms}_atoms")),
            &cif,
            |b, cif| b.iter(|| mmcif_str_to_entities(black_box(cif)).unwrap()),
        );
    }
    group.finish();
}

criterion_group!(benches, bench_pdb_parsing, bench_mmcif_parsing);
criterion_main!(benches);

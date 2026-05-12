//! Benchmarks for structural analysis operations.
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
use glam::Vec3;
use molex::ops::transform::{extract_ca_positions, kabsch_alignment};

fn bench_kabsch_alignment(c: &mut Criterion) {
    let mut group = c.benchmark_group("kabsch_alignment");

    for n_points in [10, 50, 200, 1000] {
        let reference: Vec<Vec3> = (0..n_points)
            .map(|i| {
                let t = i as f32 * 0.1;
                Vec3::new(t.cos() * 10.0, t.sin() * 10.0, t * 0.5)
            })
            .collect();

        // Rotate + translate target
        let rotation =
            glam::Mat3::from_rotation_z(0.3) * glam::Mat3::from_rotation_x(0.2);
        let translation = Vec3::new(5.0, -3.0, 2.0);
        let target: Vec<Vec3> = reference
            .iter()
            .map(|p| rotation * *p + translation)
            .collect();

        group.bench_with_input(
            BenchmarkId::new("kabsch", format!("{n_points}_points")),
            &(reference.clone(), target.clone()),
            |b, (reference, target)| {
                b.iter(|| {
                    kabsch_alignment(black_box(reference), black_box(target))
                });
            },
        );
    }

    group.finish();
}

fn bench_ca_extraction(c: &mut Criterion) {
    use molex::adapters::pdb::pdb_str_to_entities;

    let mut group = c.benchmark_group("ca_extraction");

    for n_residues in [50usize, 200, 1000] {
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
                    "ATOM  {:>5} {:<4}ALA A{:>4}    {:>8.3}{:>8.3}{:>8.3}  \
                     1.00  0.00          {:>2}  \n",
                    serial, name, res_num, x, dy, dz, elem
                ));
                serial += 1;
            }
        }
        pdb.push_str("END\n");
        let entities = pdb_str_to_entities(&pdb).unwrap();

        group.bench_with_input(
            BenchmarkId::new("extract_ca", format!("{n_residues}_residues")),
            &entities,
            |b, entities| {
                b.iter(|| extract_ca_positions(black_box(entities)));
            },
        );
    }

    group.finish();
}

criterion_group!(benches, bench_kabsch_alignment, bench_ca_extraction);
criterion_main!(benches);

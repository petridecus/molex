//! Benchmarks for structural analysis operations.
#![allow(
    missing_docs,
    unused_results,
    clippy::unwrap_used,
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_possible_wrap
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
    use molex::element::Element;
    use molex::ops::codec::{split_into_entities, Coords, CoordsAtom};

    let mut group = c.benchmark_group("ca_extraction");

    for n_residues in [50, 200, 1000] {
        let n_atoms = n_residues * 5;
        let atom_names_cycle: &[[u8; 4]] =
            &[*b"N   ", *b"CA  ", *b"C   ", *b"O   ", *b"CB  "];

        let coords = Coords {
            num_atoms: n_atoms,
            atoms: (0..n_atoms)
                .map(|i| CoordsAtom {
                    x: (i as f32) * 1.5,
                    y: 0.0,
                    z: 0.0,
                    occupancy: 1.0,
                    b_factor: 0.0,
                })
                .collect(),
            chain_ids: vec![b'A'; n_atoms],
            res_names: vec![*b"ALA"; n_atoms],
            res_nums: (0..n_atoms).map(|i| (i / 5 + 1) as i32).collect(),
            atom_names: (0..n_atoms).map(|i| atom_names_cycle[i % 5]).collect(),
            elements: vec![Element::N; n_atoms],
        };
        let entities = split_into_entities(&coords);

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

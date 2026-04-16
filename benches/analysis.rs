//! Benchmarks for structural analysis operations.
#![allow(
    missing_docs,
    unused_results,
    unused_imports,
    unused_variables,
    deprecated,
    clippy::unwrap_used,
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_possible_wrap
)]

use criterion::{
    black_box, criterion_group, criterion_main, BenchmarkId, Criterion,
};
use glam::Vec3;
use molex::analysis::bonds::hydrogen::{detect_hbonds, HBond};
use molex::analysis::ss::dssp::classify;
use molex::entity::molecule::protein::ResidueBackbone;
use molex::ops::transform::{extract_ca_positions, kabsch_alignment};

/// Generate a synthetic alpha-helix backbone with `n` residues.
/// Positions approximate ideal helix geometry (3.6 residues/turn, 1.5 Å rise).
fn make_helix_backbone(n: usize) -> Vec<ResidueBackbone> {
    (0..n)
        .map(|i| {
            let t = i as f32;
            let angle = t * std::f32::consts::TAU / 3.6;
            let rise = t * 1.5;
            let r = 2.3; // helix radius

            ResidueBackbone {
                n: Vec3::new(r * angle.cos(), r * angle.sin(), rise),
                ca: Vec3::new(
                    r * (angle + 0.3).cos(),
                    r * (angle + 0.3).sin(),
                    rise + 0.5,
                ),
                c: Vec3::new(
                    r * (angle + 0.6).cos(),
                    r * (angle + 0.6).sin(),
                    rise + 1.0,
                ),
                o: Vec3::new(
                    (r + 1.2) * (angle + 0.6).cos(),
                    (r + 1.2) * (angle + 0.6).sin(),
                    rise + 1.0,
                ),
            }
        })
        .collect()
}

fn bench_hbond_detection(c: &mut Criterion) {
    let mut group = c.benchmark_group("hbond_detection");

    for n_residues in [20, 50, 100, 300] {
        let backbone = make_helix_backbone(n_residues);

        group.bench_with_input(
            BenchmarkId::new("detect_hbonds", format!("{n_residues}_residues")),
            &backbone,
            |b, backbone| {
                b.iter(|| detect_hbonds(black_box(backbone)));
            },
        );
    }

    group.finish();
}

fn bench_dssp_classification(c: &mut Criterion) {
    let mut group = c.benchmark_group("dssp_classification");

    for n_residues in [20, 50, 100, 300] {
        let backbone = make_helix_backbone(n_residues);
        let hbonds = detect_hbonds(&backbone);

        group.bench_with_input(
            BenchmarkId::new("classify", format!("{n_residues}_residues")),
            &(hbonds.clone(), n_residues),
            |b, (hbonds, n)| {
                b.iter(|| classify(black_box(hbonds), black_box(*n)));
            },
        );
    }

    group.finish();
}

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

criterion_group!(
    benches,
    bench_hbond_detection,
    bench_dssp_classification,
    bench_kabsch_alignment,
    bench_ca_extraction
);
criterion_main!(benches);

//! Benchmarks for serialization/deserialization (COORDS01, ASSEM01).
#![allow(
    missing_docs,
    unused_results,
    clippy::unwrap_used,
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_possible_wrap,
    unused_variables
)]

use criterion::{
    black_box, criterion_group, criterion_main, BenchmarkId, Criterion,
};
use molex::element::Element;
use molex::ops::codec::{
    deserialize, deserialize_assembly, merge_entities, serialize,
    serialize_assembly, split_into_entities, Coords, CoordsAtom,
};

/// Build a synthetic Coords with `n` atoms across `n/5` ALA residues.
fn make_coords(n_atoms: usize) -> Coords {
    let n_residues = (n_atoms / 5).max(1);
    let atom_names_cycle: &[[u8; 4]] =
        &[*b"N   ", *b"CA  ", *b"C   ", *b"O   ", *b"CB  "];

    Coords {
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
        res_names: (0..n_atoms).map(|_| *b"ALA").collect(),
        res_nums: (0..n_atoms).map(|i| (i / 5 + 1) as i32).collect(),
        atom_names: (0..n_atoms).map(|i| atom_names_cycle[i % 5]).collect(),
        elements: (0..n_atoms)
            .map(|i| match i % 5 {
                0 => Element::N,
                1 | 2 | 4 => Element::C,
                3 => Element::O,
                _ => Element::Unknown,
            })
            .collect(),
    }
}

fn bench_coords_roundtrip(c: &mut Criterion) {
    let mut group = c.benchmark_group("coords_roundtrip");

    for n_atoms in [50, 500, 5000, 50_000] {
        let coords = make_coords(n_atoms);
        let bytes = serialize(&coords).unwrap();

        group.bench_with_input(
            BenchmarkId::new("serialize", format!("{n_atoms}_atoms")),
            &coords,
            |b, coords| {
                b.iter(|| serialize(black_box(coords)).unwrap());
            },
        );

        group.bench_with_input(
            BenchmarkId::new("deserialize", format!("{n_atoms}_atoms")),
            &bytes,
            |b, bytes| {
                b.iter(|| deserialize(black_box(bytes)).unwrap());
            },
        );
    }

    group.finish();
}

fn bench_assembly_roundtrip(c: &mut Criterion) {
    let mut group = c.benchmark_group("assembly_roundtrip");

    for n_atoms in [50, 500, 5000] {
        let coords = make_coords(n_atoms);
        let entities = split_into_entities(&coords);
        let bytes = serialize_assembly(&entities).unwrap();

        group.bench_with_input(
            BenchmarkId::new("serialize", format!("{n_atoms}_atoms")),
            &entities,
            |b, entities| {
                b.iter(|| serialize_assembly(black_box(entities)).unwrap());
            },
        );

        group.bench_with_input(
            BenchmarkId::new("deserialize", format!("{n_atoms}_atoms")),
            &bytes,
            |b, bytes| {
                b.iter(|| deserialize_assembly(black_box(bytes)).unwrap());
            },
        );
    }

    group.finish();
}

fn bench_merge_entities(c: &mut Criterion) {
    let mut group = c.benchmark_group("merge_entities");

    for n_atoms in [50, 500, 5000] {
        let coords = make_coords(n_atoms);
        let entities = split_into_entities(&coords);

        group.bench_with_input(
            BenchmarkId::new("merge", format!("{n_atoms}_atoms")),
            &entities,
            |b, entities| {
                b.iter(|| merge_entities(black_box(entities)));
            },
        );
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_coords_roundtrip,
    bench_assembly_roundtrip,
    bench_merge_entities
);
criterion_main!(benches);

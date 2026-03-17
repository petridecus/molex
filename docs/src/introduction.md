# molex

**molex** is a Rust library for parsing, transforming, and serializing
molecular structure data. It provides a unified type system for working
with proteins, nucleic acids, ligands, and other biomolecules across
multiple file formats.

## What it does

- **Parse** PDB, mmCIF, BinaryCIF, DCD trajectory, and MRC density files
- **Convert** between formats with a canonical intermediate representation
- **Transform** coordinates: alignment, superposition, filtering, interpolation
- **Analyze** secondary structure (DSSP), bond inference, validation
- **Extract** render-ready data (backbone chains, sidechain atoms, bonds)
- **Serialize** to a compact binary format for IPC and storage

## Key types

| Type | Description |
|------|-------------|
| `Coords` | Flat atom arrays: positions, names, chains, residues, elements |
| `MoleculeEntity` | A classified molecule (protein, ligand, etc.) with its `Coords` |
| `RenderCoords` | Backbone chains + sidechain atoms, ready for GPU consumption |
| `SSType` | Secondary structure classification (helix, sheet, coil, turn) |
| `DensityMap` | 3D electron density grid from MRC/CCP4 files |

## Design principles

1. **Zero-copy where possible.** Parsing produces owned data, but
   transforms operate on slices and iterators.
2. **Format-agnostic core.** `Coords` and `MoleculeEntity` carry no
   format-specific metadata — adapters handle the translation.
3. **Embeddable.** No filesystem, network, or GPU dependencies in the
   core. Optional `python` feature adds PyO3 bindings.

## Crate structure

```
molex/
├── types/              Core data structures (Coords, Entity, Density)
├── adapters/           Format I/O (PDB, mmCIF, BinaryCIF, DCD, MRC, AtomWorks)
├── cif/                CIF/STAR parser and typed extractors
├── ops/                Coordinate transforms, validation, bond inference
├── render/             Render-ready data extraction
├── secondary_structure/ DSSP and SS type classification
├── ffi/                C-compatible FFI layer
└── python/             PyO3 bindings (feature-gated)
```

## API documentation

For the full Rust API reference, run:

```bash
cargo doc --open --document-private-items
```

# Architecture Overview

## Module layout

```
molex/src/
├── adapters/         File format parsers (PDB, mmCIF, BinaryCIF, MRC, DCD, AtomWorks)
├── analysis/         Structural analysis (bonds, secondary structure, AABB)
├── element.rs        Element enum (symbols, covalent radii, colors)
├── entity/           Entity system
│   ├── molecule/     MoleculeEntity enum + subtypes (protein, nucleic acid, small molecule, bulk)
│   └── surface/      Surface types (VoxelGrid, Density)
├── ops/              Operations
│   ├── codec/        Wire formats (COORDS01, ASSEM01), serialize/deserialize, split/merge
│   └── transform/    Kabsch alignment, CA extraction, backbone segments
├── ffi.rs            C FFI bindings
├── python.rs         PyO3 bindings
└── lib.rs            Crate root, re-exports
```

## Entity-first design

Entities (`Vec<MoleculeEntity>`) are the primary data model.

**Adapters** parse files into entities. The `*_to_entities` functions are the primary API. The `*_to_coords` functions parse to entities, then flatten to `Coords` via `merge_entities`.

**Analysis** operates on `&[Atom]`, `&[ResidueBackbone]`, or `&[MoleculeEntity]`.

**`Coords`** is a flat, column-oriented wire format for FFI and IPC (parallel arrays of x/y/z, chain IDs, residue names, etc.).

## Entity classification

When a file is parsed, atoms are grouped into entities by chain ID, residue name, and molecule type. The `classify_residue` function maps 3-letter residue codes to `MoleculeType` values using lookup tables:

| MoleculeType | Examples |
|---|---|
| `Protein` | Standard amino acids (ALA, GLY, ...) |
| `DNA` | DA, DT, DC, DG |
| `RNA` | A, U, C, G |
| `Ligand` | ATP, HEM, NAG, ... |
| `Ion` | ZN, MG, CA, FE, ... |
| `Water` | HOH, WAT, DOD |
| `Lipid` | OLC, PLM, ... |
| `Cofactor` | HEM, NAD, FAD, ... |
| `Solvent` | GOL, EDO, PEG, ... |

## Type hierarchy

```
MoleculeEntity (enum)
├── Protein(ProteinEntity)      -- chain with residues, segment breaks
├── NucleicAcid(NAEntity)       -- DNA/RNA chain with residues
├── SmallMolecule(SmallMoleculeEntity) -- single ligand, ion, cofactor, lipid
└── Bulk(BulkEntity)            -- grouped water or solvent molecules

Entity trait   -- id(), molecule_type(), atoms(), positions(), atom_count()
Polymer trait  -- residues(), segment_breaks(), segment_count(), segment_range()
```

`ProteinEntity` and `NAEntity` implement both `Entity` and `Polymer`. `SmallMoleculeEntity` and `BulkEntity` implement `Entity` only.

## Surface types

The `entity::surface` module provides volumetric data types:

- **`VoxelGrid`** -- a generic 3D grid (`ndarray::Array3<f32>`) with crystallographic cell metadata. Handles fractional-to-Cartesian coordinate conversion.
- **`Density`** -- wraps `VoxelGrid` with density-specific metadata. Constructed by the MRC adapter.

## Binary formats

molex defines two compact binary formats for IPC:

- **COORDS01** -- flat atom array with element data (magic: `COORDS01`)
- **ASSEM01** -- entity-aware format preserving molecule type metadata per entity (magic: `ASSEM01\0`)

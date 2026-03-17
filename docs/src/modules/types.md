# Types — Core Data Structures

The `types` module contains molex's canonical data representations.

## `types::coords` — Atom-Level Data

### `Coords`

Flat parallel arrays representing atoms in a molecular structure.
Every atom has a position (`CoordsAtom`), chain ID, residue name,
residue number, atom name, and element type.

This is the lowest-level representation — no hierarchy, no
classification. Use `split_into_entities()` to get grouped molecules.

### `CoordsAtom`

```rust,ignore
pub struct CoordsAtom {
    pub x: f32, pub y: f32, pub z: f32,
    pub occupancy: f32,
    pub b_factor: f32,
}
```

### `Element`

Chemical element enum with methods for covalent radius, VDW radius,
and CPK color. Supports lookup by symbol or by atom name heuristics.

### Binary serialization

- `serialize(&coords) → Vec<u8>` — COORDS01 format (26 bytes/atom)
- `deserialize(&bytes) → Coords` — inverse
- `serialize_assembly(&entities) → Vec<u8>` — ASSEM01 multi-entity format
- `deserialize_assembly(&bytes) → Vec<MoleculeEntity>` — inverse

## `types::entity` — Molecule Classification

### `MoleculeEntity`

A classified molecule with its atom data:

- `entity_id: u32` — unique identifier
- `molecule_type: MoleculeType` — Protein, DNA, RNA, Ligand, Ion, etc.
- `kind: EntityKind` — `Polymer(PolymerData)` or `AtomSet(AtomSet)`

### `MoleculeType`

```rust,ignore
pub enum MoleculeType {
    Protein, DNA, RNA, Ligand, Ion, Water, Lipid, Cofactor, Solvent,
}
```

### Key functions

- `split_into_entities(&coords) → Vec<MoleculeEntity>` — classify and group
- `merge_entities(&[MoleculeEntity]) → Coords` — flatten back to atoms
- `classify_residue(name) → MoleculeType` — single residue lookup

## `types::density` — Electron Density

### `DensityMap`

3D grid of electron density values from MRC/CCP4 files:

- Grid dimensions, cell parameters, origin
- `f32` voxel data
- Methods for coordinate-to-grid mapping

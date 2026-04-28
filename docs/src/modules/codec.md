# Codec

Wire formats, serialization, and entity splitting/merging live in `molex::ops::codec` (source: `src/ops/codec/`).

## Coords

The flat, column-oriented wire format:

```rust,ignore
pub struct Coords {
    pub num_atoms: usize,
    pub atoms: Vec<CoordsAtom>,       // x, y, z, occupancy, b_factor
    pub chain_ids: Vec<u8>,           // per-atom chain ID byte
    pub res_names: Vec<[u8; 3]>,      // per-atom residue name
    pub res_nums: Vec<i32>,           // per-atom residue number
    pub atom_names: Vec<[u8; 4]>,     // per-atom PDB atom name
    pub elements: Vec<Element>,       // per-atom element
}
```

`Coords` is a serialization/interop format for FFI and IPC.

## Wire formats

### COORDS01

Flat atom array binary format. Header magic: `COORDS01` (backward-compatible reader also handles `COORDS00` which omits elements).

```rust,ignore
use molex::ops::codec::{serialize, deserialize};

let bytes: Vec<u8> = serialize(&coords)?;
let coords: Coords = deserialize(&bytes)?;
```

### ASSEM01

Entity-aware binary format that preserves molecule type metadata per entity. Header magic: `ASSEM01\0`. Use this when you need to round-trip entities without re-running residue classification.

```rust,ignore
use molex::ops::codec::{serialize_assembly, deserialize_assembly};

let bytes: Vec<u8> = serialize_assembly(&entities)?;
let entities: Vec<MoleculeEntity> = deserialize_assembly(&bytes)?;
```

## Entity splitting and merging

```rust,ignore
use molex::ops::codec::{split_into_entities, merge_entities};

// Coords -> entities (classifies residues, groups by chain)
let entities: Vec<MoleculeEntity> = split_into_entities(&coords);

// Entities -> Coords (flattens back)
let coords: Coords = merge_entities(&entities);
```

`split_into_entities` groups atoms by:
- Chain ID + molecule type for polymers (one entity per chain)
- Chain ID + residue number for small molecules (one entity per molecule)
- All waters into a single `Bulk` entity
- All solvents into a single `Bulk` entity

## Transforms (`ops::transform`)

Structural alignment and extraction utilities:

```rust,ignore
use molex::ops::transform::*;

// Kabsch alignment (minimize RMSD between point sets) — returns Option
let (rotation, translation) = kabsch_alignment(&source, &target)?;
let (rotation, translation, scale) = kabsch_alignment_with_scale(&source, &target)?;

// Apply transform to entities (rotation/translation are Copy)
transform_entities(&mut entities, rotation, translation);
transform_entities_with_scale(&mut entities, rotation, translation, scale);

// Align entities to a reference structure
align_to_reference(&mut mobile, &reference);

// Extract alpha-carbon positions
let ca_positions: Vec<Vec3> = extract_ca_positions(&entities);
let ca_by_chain: Vec<Vec<Vec3>> = extract_ca_from_chains(&entities);

// Get continuous backbone segments
let segments: Vec<Vec<Vec3>> = extract_backbone_segments(&entities);

// Compute centroid
let center: Vec3 = centroid(&positions);
```

## ChainIdMapper

Maps multi-character chain ID strings (e.g. "AA", "AB" for structures with >26 chains) to unique `u8` values for the `Coords` format:

```rust,ignore
let mut mapper = ChainIdMapper::new();
let id = mapper.get_or_assign("AA"); // assigns a unique byte
```

## Error type

All codec operations return `CoordsError`:

```rust,ignore
pub enum CoordsError {
    InvalidFormat(String),
    PdbParseError(String),
    SerializationError(String),
}
```

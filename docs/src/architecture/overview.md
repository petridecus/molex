# Architecture Overview

molex is organized as a layered conversion pipeline. Raw file bytes
enter through adapters, pass through a canonical intermediate
representation, and exit as either transformed coordinates, render-ready
geometry, or binary-serialized bytes.

## Layer diagram

```
                          ┌─────────────┐
  PDB ──►┐               │   Coords    │               ┌──► Binary COORDS
  CIF ──►├─► adapters ──►│ (canonical) ├──► ops ───────►├──► PDB string
 BCIF ──►├─►             │  + Entity   │   transform    ├──► RenderCoords
  DCD ──►┘               │  classify   │   validate     └──► AtomArray (Py)
                          └──────┬──────┘   align
                                 │
                                 ▼
                          secondary_structure
                          (DSSP → SSType[])
```

## Core types

### `Coords`

The canonical atom-level representation. Flat parallel arrays:

```rust,ignore
pub struct Coords {
    pub num_atoms: usize,
    pub atoms: Vec<CoordsAtom>,      // x, y, z, occupancy, b_factor
    pub chain_ids: Vec<u8>,
    pub res_names: Vec<[u8; 3]>,
    pub res_nums: Vec<i32>,
    pub atom_names: Vec<[u8; 4]>,
    pub elements: Vec<Element>,
}
```

Flat arrays make iteration, slicing, and binary serialization cheap.
The tradeoff is no hierarchical chain→residue→atom tree — use
`MoleculeEntity` when you need entity-level grouping.

### `MoleculeEntity`

A classified molecule with its coordinates:

```rust,ignore
pub struct MoleculeEntity {
    pub entity_id: u32,
    pub molecule_type: MoleculeType,  // Protein, DNA, RNA, Ligand, ...
    pub kind: EntityKind,             // Polymer or AtomSet
}
```

`split_into_entities()` classifies residues by name and groups them
into entities. This is the primary input to viso's rendering engine.

### `RenderCoords`

Extracted backbone chains (N-CA-C triples) and sidechain atoms with
bond connectivity. Bridge between `Coords` and GPU renderers.

## Module responsibilities

| Module | Responsibility |
|--------|---------------|
| `types` | `Coords`, `MoleculeEntity`, `DensityMap`, binary serialization |
| `adapters` | Format I/O: PDB, mmCIF, BinaryCIF, DCD, MRC, AtomWorks |
| `cif` | Low-level CIF/STAR parser with typed extractors |
| `ops` | Coordinate transforms, Kabsch alignment, bond inference, validation |
| `render` | Backbone/sidechain extraction, sequence extraction |
| `secondary_structure` | DSSP algorithm, SS type classification |
| `ffi` | `extern "C"` functions for C/C++ integration |
| `python` | PyO3 bindings (feature-gated) |

## Binary formats

molex defines two compact binary formats for IPC:

- **COORDS01** — single molecule: magic + atom data (26 bytes/atom)
- **ASSEM01** — multi-entity assembly: magic + entity headers + atom data

Both use big-endian encoding and are designed for zero-overhead
round-tripping between Rust and C++ backends.

# Data Flow

## Overview

```
                       ┌──────────────┐
 PDB / mmCIF / BCIF ──>│              ├──> Vec<MoleculeEntity>
 MRC / CCP4         ──>│   Adapters   ├──> Density (wraps VoxelGrid)
 DCD                ──>│              ├──> Vec<DcdFrame>
                       └──────┬───────┘
                              │
                              v
                  ┌───────────────────────┐
                  │       Assembly        │
                  │  (entities + derived: │
                  │   ss, hbonds, S-S)    │
                  └──┬────────┬────────┬──┘
                     │        │        │
                     v        v        v
              ┌──────────┐ ┌────────┐ ┌─────────┐
              │ Analysis │ │  Ops   │ │  Codec  │
              │          │ │        │ │         │
              │ dssp     │ │ kabsch │ │serialize│
              │ bonds    │ │ align  │ │serialize│
              │ disulfide│ │extract │ │_assembly│
              │ aabb     │ │        │ │         │
              └──────────┘ └────────┘ └────┬────┘
                                           │
                                           v
                                  FFI / IPC / Python
```

Analysis, Transform, and Codec are independent — use any combination depending on what you need.

## 1. Parsing

Every structure adapter returns `Vec<MoleculeEntity>`:

```rust,ignore
let entities = pdb_file_to_entities(Path::new("1ubq.pdb"))?;
let entities = mmcif_file_to_entities(Path::new("3nez.cif"))?;
let entities = bcif_file_to_entities(Path::new("1ubq.bcif"))?;
```

Density and trajectory adapters return their own types:

```rust,ignore
let density = mrc_file_to_density(Path::new("emd_1234.map"))?;
let frames = dcd_file_to_frames(Path::new("trajectory.dcd"))?;
```

## 2. Entity splitting

`split_into_entities` groups atoms by:

1. Chain ID + molecule type for polymers (one entity per chain)
2. Chain ID + residue number for small molecules (one entity each)
3. All waters into a single `Bulk` entity
4. All solvents into a single `Bulk` entity

Each entity gets a unique `EntityId`.

## 3. Analysis

```rust,ignore
use molex::{detect_disulfides, Assembly};
use molex::analysis::{infer_bonds, DEFAULT_TOLERANCE};

// Top-level pipeline: build an Assembly to get DSSP, H-bonds,
// and cross-entity disulfides eagerly computed and kept in sync.
let assembly = Assembly::new(entities);
let hbonds = assembly.hbonds();
let ss     = assembly.ss_types(entity_id);

// Or call the building blocks directly:
let bonds      = infer_bonds(atoms, DEFAULT_TOLERANCE);   // distance-based
let disulfides = detect_disulfides(&entities);             // CYS SG-SG
let aabb       = entity.aabb();
```

## 4. Transforms

```rust,ignore
let (rotation, translation) = kabsch_alignment(&reference_ca, &target_ca)?;
transform_entities(&mut entities, rotation, translation);
let ca_positions = extract_ca_positions(&entities);
```

## 5. Serialization

For sending to C/C++/Python consumers:

```rust,ignore
// COORDS01 (flat atom array)
let bytes = serialize(&merge_entities(&entities))?;

// ASSEM01 (preserves entity types)
let bytes = serialize_assembly(&entities)?;
```

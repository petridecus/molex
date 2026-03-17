# Data Flow

This page traces the typical lifecycle of molecular data through
molex, from file input to render-ready output.

## 1. Parsing

Every supported format has a dedicated adapter that produces `Coords`:

```
PDB file  → adapters::pdb::pdb_to_coords()       → Vec<u8> (COORDS01)
CIF file  → adapters::pdb::mmcif_to_coords()      → Vec<u8> (COORDS01)
BCIF file → adapters::bcif::bcif_to_coords()       → Vec<u8> (COORDS01)
DCD file  → adapters::dcd::dcd_file_to_frames()    → Vec<DcdFrame>
MRC file  → adapters::mrc::mrc_to_density()        → DensityMap
```

The binary `Vec<u8>` is the COORDS01 format. Deserialize it to get
the `Coords` struct:

```rust,ignore
let coords = molex::types::coords::deserialize(&bytes)?;
```

## 2. Entity classification

```rust,ignore
let entities = molex::types::entity::split_into_entities(&coords);
// entities: Vec<MoleculeEntity>
// Each entity has: entity_id, molecule_type, kind (Polymer or AtomSet)
```

Classification is based on residue name lookup — standard amino acids
become `Protein`, nucleotides become `DNA`/`RNA`, `HOH` becomes
`Water`, and everything else is classified as `Ligand`, `Ion`,
`Cofactor`, etc.

## 3. Transforms

The `ops::transform` module provides coordinate manipulation:

```rust,ignore
// Filter to protein-only atoms
let protein = ops::transform::protein_only(&coords);

// Kabsch superposition onto a reference
let (aligned, rmsd) = ops::transform::kabsch_alignment(&mobile, &target);

// Extract backbone chains as Vec<Vec3> (N-CA-C triples)
let chains = ops::transform::extract_backbone_chains(&coords);
```

## 4. Secondary structure

DSSP classification from backbone geometry:

```rust,ignore
let ss_types: Vec<SSType> = secondary_structure::dssp::from_entity(&entity);
// SSType::Helix, SSType::Sheet, SSType::Coil, SSType::Turn
```

## 5. Render extraction

`RenderCoords` splits atoms into backbone and sidechain data suitable
for GPU rendering:

```rust,ignore
let render = RenderCoords::from_entity(&entity, is_hydrophobic, get_bonds);
// render.backbone_chains: Vec<Vec<Vec3>>  (N-CA-C per chain)
// render.sidechain_atoms: Vec<RenderSidechainAtom>
// render.sidechain_bonds: Vec<(u32, u32)>
```

## 6. Serialization

For IPC with C++ backends or storage:

```rust,ignore
// Single molecule
let bytes = molex::types::coords::serialize(&coords)?;

// Multi-entity assembly
let bytes = molex::types::coords::serialize_assembly(&entities)?;
```

## Pipeline summary

```
File → Adapter → Coords → split_into_entities → MoleculeEntity[]
                    │                                   │
                    ├──► ops::transform (align, filter)  │
                    ├──► ops::validation (completeness)  │
                    └──► serialize (binary IPC)          │
                                                        ▼
                                            secondary_structure::dssp
                                                        │
                                                        ▼
                                              RenderCoords::from_entity
                                              (backbone + sidechain data)
```

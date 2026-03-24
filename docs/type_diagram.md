# Molex Type System

## High-Level Data Flow

```
┌─────────────────────────────────┐
│       Raw File Bytes            │
│   (PDB / mmCIF / BCIF / MRC)   │
└───────────────┬─────────────────┘
                │
                │  adapters/{pdb,bcif,mrc}.rs
                │
    ┌───────────┼───────────────┐
    ▼           ▼               ▼
 Coords    Vec<Molecule-    DensityMap
 (SoA,      Entity>         (volumetric
  flat)    (structured)      grid)
    │           │
    │  split_   │
    │  into_    │
    │  entities │
    │     │     │
    └─────┘     │
          │     │
          ▼     ▼
    ┌───────────────────────────────────────────┐
    │           MoleculeEntity                  │
    │       (canonical source of truth)         │
    │                                           │
    │  entity_id: u32                           │
    │  molecule_type: MoleculeType              │
    │  kind: EntityKind ─┐                      │
    │                    │                      │
    │    ┌───────────────┼───────────────┐      │
    │    ▼               ▼               ▼      │
    │  Polymer       SmallMolecule     Bulk     │
    │  (PolymerData) (AtomSet + name)  (AtomSet)│
    └───────────────────────────────────────────┘
                │
                │  Derivation methods
                │
    ┌───────────┼──────────────┬──────────────┐
    ▼           ▼              ▼              ▼
to_residues() to_backbone() to_coords()  extract_*()
    │           │              │          (legacy)
    ▼           ▼              ▼
Vec<Protein-  Vec<Residue-   Coords
 Residue>      Backbone>     (SoA for
                             serialization)
```

## Entity Type Hierarchy

```
MoleculeType (enum)
├── Protein  ──→ EntityKind::Polymer(PolymerData)
├── DNA      ──→ EntityKind::Polymer(PolymerData)
├── RNA      ──→ EntityKind::Polymer(PolymerData)
├── Ligand   ──→ EntityKind::SmallMolecule
├── Ion      ──→ EntityKind::SmallMolecule
├── Cofactor ──→ EntityKind::SmallMolecule
├── Lipid    ──→ EntityKind::SmallMolecule
├── Water    ──→ EntityKind::Bulk
└── Solvent  ──→ EntityKind::Bulk
```

## Polymer Data Hierarchy

```
PolymerData
├── atoms: AtomSet
│   ├── atoms: Vec<CoordsAtom>   (x, y, z, occupancy, b_factor)
│   ├── atom_names: Vec<[u8;4]>  ("N   ", "CA  ", "CB  ", ...)
│   └── elements: Vec<Element>   (C, N, O, S, ...)
│
└── chains: Vec<PolymerChain>
    └── PolymerChain
        ├── chain_id: u8
        └── residues: Vec<Residue>
            └── Residue
                ├── name: [u8;3]          ("ALA", "CYS", ...)
                ├── number: i32           (PDB sequence number)
                └── atom_range: Range     (into PolymerData.atoms)
```

## Protein-Specific Types (types/protein.rs)

These are the primary extraction types from protein entities.

```
                       MoleculeEntity
                            │
                    to_residues(is_hydrophobic, get_bonds)
                            │
                            ▼
                   Vec<ProteinResidue>
                            │
              ┌─────────────┼──────────────┐
              ▼             │              ▼
      ResidueBackbone       │         Sidechain
      ├── n: Vec3           │         ├── atoms: AtomSet
      ├── ca: Vec3          │         ├── bonds: Vec<(usize, usize)>
      ├── c: Vec3           │         │   (local indices)
      └── o: Vec3           │         └── is_hydrophobic: bool
                            │
                     ProteinResidue
                     ├── name: [u8;3]
                     ├── number: i32
                     ├── chain_id: u8
                     ├── backbone: ResidueBackbone
                     └── sidechain: Sidechain
```

### Convenience derivations

```
to_backbone()    → Vec<ResidueBackbone>     (just backbone from each residue)
to_sidechains()  → (not yet implemented)
to_coords()      → Coords                  (SoA flat format for serialization)
```

## Analysis Module (analysis/)

Structural analysis: hydrogen bond detection and secondary structure
classification. Replaces the old `secondary_structure` module.

```
              Vec<ResidueBackbone>
              (from one or more entities,
               concatenated for cross-entity analysis)
                        │
           ┌────────────┼──────────────┐
           ▼            │              ▼
   detect_hbonds()      │       detect_dssp()
   (H-bonds only)       │       (both at once)
           │            │         │         │
           ▼            │         ▼         ▼
      Vec<HBond>        │    Vec<SSType>  Vec<HBond>
      ├── donor         │
      ├── acceptor      │
      └── energy        │
                        │
                   resolve_ss(override, residues)
                        │
                        ▼
                   Vec<SSType>
                   (Helix / Sheet / Coil)
```

### Module layout

```
analysis/
├── mod.rs        SSType, HBond, detect_dssp(), resolve_ss(),
│                 merge_short_segments()
├── hbond.rs      detect_hbonds(&[ResidueBackbone]) → Vec<HBond>
│                 Kabsch-Sander energy: E = 27.888 * (1/rON + 1/rCH
│                 - 1/rOH - 1/rCN), threshold < -0.5 kcal/mol
├── ss.rs         classify(&[HBond], n) → Vec<SSType>
│                 Helix: i→i+4, i→i+3, i→i+5 turn patterns
│                 Sheet: parallel/antiparallel bridge patterns
└── ss_string.rs  from_string("HHHEEECCC") → Vec<SSType>
```

## Coords (SoA Serialization Format)

Flat struct-of-arrays layout for binary serialization and GPU upload.
Derivable from entities via `to_coords()`, not the working data model.

```
Coords
├── num_atoms: usize
├── atoms: Vec<CoordsAtom>     (x, y, z, occupancy, b_factor)
├── chain_ids: Vec<u8>
├── res_names: Vec<[u8;3]>
├── res_nums: Vec<i32>
├── atom_names: Vec<[u8;4]>
└── elements: Vec<Element>

Binary formats:
  COORDS01 — per-atom: 3×f32 pos + f32 occ + f32 bfac + ... (26 bytes)
  ASSEM01  — entity headers + COORDS01 per entity
```

## Legacy Types (deprecated, backward compat)

These types exist for backward compatibility and will be removed
in a future release.

```
secondary_structure/     → re-exports from analysis/
  BackboneResidue        → type alias for ResidueBackbone
  DetectionInput         → wrapper calling analysis functions
  auto::detect           → CA-distance heuristic (to be removed)

render/
  RenderBackboneResidue  → parallel to ResidueBackbone (to be consolidated)
  RenderCoords           → bootstrap-only snapshot (to be simplified)
  ProteinBackbone        → replaced by Vec<ResidueBackbone>
  BackboneChain          → replaced by direct Vec<Vec3> usage
  SidechainAtoms         → replaced by Vec<Sidechain>
  SidechainAtomData      → replaced by AtomSet inside Sidechain
```

## Operations (ops/)

```
ops/
├── bond_inference.rs    infer_bonds(&Coords) → Vec<InferredBond>
│                        Distance-based, for small molecules
├── transform/
│   ├── mod.rs           protein_only(), backbone_only(),
│   │                    heavy_atoms_only(), filter_atoms()
│   ├── extract.rs       extract_backbone_chains(), extract_ca_positions()
│   ├── alignment.rs     kabsch_alignment(), transform_coords()
│   └── interpolate.rs   interpolate_coords() (animation lerp)
└── validation.rs        validate_completeness(), has_complete_backbone()
```

## Adapters (file format converters)

```
adapters/
├── pdb.rs       PDB + mmCIF → Coords or Vec<MoleculeEntity>
├── bcif/        BinaryCIF → Coords or Vec<MoleculeEntity>
├── dcd.rs       DCD trajectory → Vec<DcdFrame>
├── mrc/         MRC/CCP4 → DensityMap
└── atomworks/   Python NumPy ↔ entities (feature-gated)
```

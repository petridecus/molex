# Analysis

Structural analysis lives in `molex::analysis` (source: `src/analysis/`). All analysis functions operate on entity-level types (`&[Atom]`, `&[ResidueBackbone]`).

## Secondary structure (`analysis::ss`)

DSSP-based secondary structure classification.

```rust,ignore
use molex::analysis::{detect_dssp, resolve_ss, SSType};

// Full DSSP: detect H-bonds, then classify
let (ss_types, hbonds) = detect_dssp(&backbone_residues);
// ss_types: Vec<SSType> -- one per residue (Helix, Sheet, or Coil)
// hbonds: Vec<HBond> -- backbone H-bond pairs that produced the assignment

// With optional override (e.g. from mmCIF annotation)
let ss = resolve_ss(Some(&override_ss), &backbone_residues);
// Falls back to DSSP if override is None
```

`SSType` is a Q3 classification:

```rust,ignore
pub enum SSType { Helix, Sheet, Coil }
```

Each variant has a `.color()` method returning an RGB `[f32; 3]` for rendering.

Short isolated segments (1-residue helix/sheet runs) are automatically merged to `Coil` by `merge_short_segments`.

### SS from string

`analysis::ss::from_string` parses secondary structure strings (e.g. `"HHHCCCEEE"`) into `Vec<SSType>`.

## Bond detection (`analysis::bonds`)

### Covalent bonds

```rust,ignore
use molex::analysis::{infer_bonds, InferredBond, BondOrder, DEFAULT_TOLERANCE};

let bonds: Vec<InferredBond> = infer_bonds(atoms, tolerance);
// InferredBond { atom_a: usize, atom_b: usize, order: BondOrder }
// BondOrder: Single, Double, Triple, Aromatic
```

Distance-based inference using element covalent radii with a configurable tolerance (default: `DEFAULT_TOLERANCE`).

### Hydrogen bonds

```rust,ignore
use molex::analysis::detect_hbonds;

let hbonds: Vec<HBond> = detect_hbonds(&backbone_residues);
// HBond { donor: usize, acceptor: usize } -- residue indices
```

DSSP-style backbone N-H...O=C hydrogen bond detection using electrostatic energy criteria.

### Disulfide bonds

```rust,ignore
use molex::analysis::{detect_disulfide_bonds, DisulfideBond};

let disulfides: Vec<DisulfideBond> = detect_disulfide_bonds(atoms);
```

Detects CYS SG-SG bonds by distance.

## Bounding box (`analysis::aabb`)

```rust,ignore
use molex::analysis::Aabb;

let aabb = Aabb::from_positions(&positions)?;
aabb.center();   // Vec3 -- geometric center
aabb.extents();  // Vec3 -- size along each axis
aabb.radius();   // f32 -- bounding sphere radius

let merged = aabb.union(&other_aabb);
let combined = Aabb::from_aabbs(&[aabb1, aabb2, aabb3])?;
```

Also available directly on entities: `entity.aabb()`.

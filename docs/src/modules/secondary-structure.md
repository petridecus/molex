# Secondary Structure

The `secondary_structure` module classifies protein residues into
helix, sheet, coil, and turn based on backbone geometry.

## `SSType`

```rust,ignore
pub enum SSType {
    Helix,   // alpha-helix (H)
    Sheet,   // beta-strand (E)
    Coil,    // unstructured (C)
    Turn,    // hydrogen-bonded turn (T)
}
```

Implements `From<char>` for Q8-style single-letter codes and
`Display` for the reverse mapping.

## DSSP algorithm (`dssp` submodule)

The primary classification method. Computes hydrogen bond energies
from backbone N-H...O=C geometry and assigns secondary structure
based on the standard DSSP criteria.

```rust,ignore
use molex::secondary_structure::dssp;

// From a MoleculeEntity
let ss: Vec<SSType> = dssp::from_entity(&entity);

// From a Q8 string (e.g., from PDB HELIX/SHEET records)
let ss: Vec<SSType> = dssp::from_string("HHHHCCCEEEECCC");
```

## Auto-detection (`auto` submodule)

Fallback method using dihedral angle (phi/psi) ranges when full
DSSP is not needed:

```rust,ignore
use molex::secondary_structure::auto::detect;

let ss: Vec<SSType> = detect(&backbone_residues);
```

## Resolution and override

The `resolve()` function merges DSSP output with optional per-entity
Q8 overrides:

```rust,ignore
use molex::secondary_structure::{resolve, DetectionInput};

let input = DetectionInput { entity, override_q8: Some("HHHCCCEEE") };
let ss = resolve(&input);
```

This is the entry point used by viso's scene metadata pipeline.

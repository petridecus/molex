# Quick Start

## Parse a PDB file into entities

```rust,ignore
use molex::adapters::pdb::pdb_file_to_entities;
use std::path::Path;

let entities = pdb_file_to_entities(Path::new("1ubq.pdb"))?;
for e in &entities {
    println!("{}: {} atoms", e.label(), e.atom_count());
}
// Output:
//   Protein A: 660 atoms
//   Water (58 molecules): 58 atoms
```

## Auto-detect format by extension

`structure_file_to_entities` dispatches on `.pdb`/`.ent` vs mmCIF:

```rust,ignore
use molex::adapters::pdb::structure_file_to_entities;
use std::path::Path;

let entities = structure_file_to_entities(Path::new("3nez.cif"))?;
```

## Work with entities

```rust,ignore
use molex::{MoleculeEntity, MoleculeType};

// Filter to protein chains
let proteins: Vec<_> = entities.iter()
    .filter(|e| e.molecule_type() == MoleculeType::Protein)
    .collect();

// Access protein-specific data
for entity in &proteins {
    let protein = entity.as_protein().unwrap();
    let backbone = protein.to_backbone();
    println!(
        "Chain {}: {} residues, {} segments",
        protein.pdb_chain_id as char,
        protein.residues.len(),
        protein.segment_count(),
    );
}
```

## Run DSSP secondary structure assignment

```rust,ignore
use molex::analysis::{detect_dssp, SSType};

let protein = entities[0].as_protein().unwrap();
let backbone = protein.to_backbone();
let (ss_types, hbonds) = detect_dssp(&backbone);
for (i, ss) in ss_types.iter().enumerate() {
    println!("Residue {}: {:?}", i, ss); // Helix, Sheet, or Coil
}
```

## Serialize to COORDS binary (for FFI/IPC)

```rust,ignore
use molex::ops::codec::{serialize, merge_entities};

let coords = merge_entities(&entities);
let bytes = serialize(&coords)?;
// Send `bytes` over FFI, IPC, or network
```

## Python usage

```python
import molex

# PDB round-trip
coords_bytes = molex.pdb_to_coords(pdb_string)
pdb_back = molex.coords_to_pdb(coords_bytes)

# Entity-aware AtomArray conversion (for ML pipelines)
atom_array = molex.entities_to_atom_array(assembly_bytes)
assembly_bytes = molex.atom_array_to_entities(atom_array)
```

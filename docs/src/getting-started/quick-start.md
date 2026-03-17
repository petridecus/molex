# Quick Start

## Parse a PDB file

```rust,ignore
use molex::adapters::pdb::structure_file_to_coords;
use molex::types::entity::split_into_entities;

// Parse any PDB or mmCIF file
let coords = structure_file_to_coords("1ubq.pdb")?;

// Split into classified entities (protein, ligand, water, etc.)
let entities = split_into_entities(&coords);
for entity in &entities {
    println!(
        "Entity {}: {:?} ({} atoms)",
        entity.entity_id,
        entity.molecule_type,
        entity.atom_count(),
    );
}
```

## Extract backbone chains

```rust,ignore
use molex::ops::transform::{protein_only, extract_backbone_chains};

let protein_coords = protein_only(&coords);
let chains = extract_backbone_chains(&protein_coords);
for (i, chain) in chains.iter().enumerate() {
    println!("Chain {}: {} CA positions", i, chain.len());
}
```

## Compute secondary structure

```rust,ignore
use molex::secondary_structure::dssp::from_entity;

let ss_types = from_entity(&entities[0]);
for (i, ss) in ss_types.iter().enumerate() {
    println!("Residue {}: {:?}", i, ss);
}
```

## Convert between formats

```rust,ignore
use molex::adapters::pdb::{pdb_to_coords, coords_to_pdb};
use molex::types::coords::{serialize, deserialize};

// PDB string → binary COORDS
let coords_bytes = pdb_to_coords(pdb_string)?;

// Binary COORDS → Coords struct
let coords = deserialize(&coords_bytes)?;

// Coords struct → PDB string
let pdb_out = coords_to_pdb(&coords_bytes)?;
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

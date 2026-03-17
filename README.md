# molex

**Mol**ecular **ex**change — a Rust library for parsing, transforming,
and serializing molecular structure data. Provides a unified type system
for proteins, nucleic acids, ligands, and other biomolecules across
multiple file formats.

## Features

- **Parse** PDB, mmCIF, BinaryCIF, DCD trajectory, and MRC density files
- **Convert** between formats via a canonical intermediate representation
- **Transform** coordinates: Kabsch alignment, filtering, interpolation
- **Analyze** secondary structure (DSSP), bond inference, validation
- **Extract** render-ready geometry (backbone chains, sidechains, bonds)
- **Serialize** to compact binary formats for IPC and storage

## Quick start

```rust
use molex::{MoleculeEntity, MoleculeType};
use molex::adapters::pdb::structure_file_to_entities;

let entities = structure_file_to_entities("1ubq.pdb".as_ref())?;
for e in &entities {
    println!("{:?}: {} atoms", e.molecule_type, e.atom_count());
}
```

## Optional features

| Feature  | Description |
|----------|-------------|
| `python` | PyO3 bindings for use from Python (requires `pyo3`, `numpy`) |

## Documentation

- **API docs:** `cargo doc --open`
- **Book:** `cd docs && mdbook serve`

## License

MIT

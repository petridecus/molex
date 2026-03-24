# molex

**molex** (molecular exchange) is a Rust library for parsing, analyzing, transforming, and serializing molecular structure data. It supports PDB, mmCIF, BinaryCIF, MRC/CCP4 density maps, and DCD trajectories.

## Key concepts

- **`MoleculeEntity`** represents a single molecule: a protein chain, a DNA/RNA strand, a ligand, an ion, or a group of waters. Parsing a structure file produces a `Vec<MoleculeEntity>`.

- **`Atom`** holds a position, element, atom name, occupancy, and B-factor. Residue and chain context live on the entity that contains the atom.

- **`Coords`** is a binary serialization format used for FFI and IPC (e.g. iceoryx zero-copy between processes).

- **Analysis** includes covalent bond inference, DSSP hydrogen bond detection, disulfide bridges, and secondary structure classification.

- **`VoxelGrid`** and **`Density`** represent 3D volumetric data (electron density, cryo-EM maps).

## Crate features

| Feature  | Description |
|----------|-------------|
| `default` | Core Rust library (no Python) |
| `python` | PyO3 bindings + AtomWorks interop |

## API documentation

For the full Rust API reference, run:

```bash
cargo doc --open --document-private-items
```

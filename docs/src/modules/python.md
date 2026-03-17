# Python Bindings

molex provides Python bindings via PyO3, enabled with the `python`
feature flag. The module is built with maturin.

## Installation

```bash
cd crates/molex
maturin develop --release --features python
```

## Core functions

### `pdb_to_coords(pdb_string) → bytes`

Parse a PDB string and return COORDS01 binary bytes.

### `mmcif_to_coords(cif_string) → bytes`

Parse an mmCIF string and return COORDS01 binary bytes.

### `coords_to_pdb(coords_bytes) → str`

Convert COORDS01 bytes back to a PDB-format string.

### `deserialize_coords(coords_bytes) → dict`

Deserialize COORDS01 bytes into a Python dictionary with NumPy arrays:

```python
result = molex.deserialize_coords(coords_bytes)
# result["num_atoms"]: int
# result["x"], result["y"], result["z"]: np.ndarray[f32]
```

## AtomWorks adapters

For ML model pipelines that use Biotite `AtomArray` objects:

### `entities_to_atom_array(assembly_bytes) → AtomArray`

Convert ASSEM01 bytes to a Biotite AtomArray with standard and
AtomWorks-specific annotations (entity_id, mol_type, chain_type).

### `entities_to_atom_array_plus(assembly_bytes) → AtomArray`

Like `entities_to_atom_array` but also populates bonds via distance
inference.

### `atom_array_to_entities(atom_array) → bytes`

Convert a Biotite AtomArray back to ASSEM01 bytes.

### `entities_to_atom_array_parsed(assembly_bytes, filename) → AtomArray`

Convert via the full AtomWorks cleaning pipeline (leaving group
removal, charge correction, missing atom imputation).

### `parse_file_to_entities(path) → bytes`

Parse a structure file (PDB/mmCIF) directly to ASSEM01 bytes.

### `parse_file_full(path) → AtomArray`

Parse a structure file through the full AtomWorks pipeline.

### `coords_to_atom_array(coords_bytes) → AtomArray`

Convert single-molecule COORDS01 bytes to AtomArray.

### `coords_to_atom_array_plus(coords_bytes) → AtomArray`

Like `coords_to_atom_array` with bond inference.

### `atom_array_to_coords(atom_array) → bytes`

Convert AtomArray to single-molecule COORDS01 bytes.

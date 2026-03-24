# Python Bindings

Python bindings are available when molex is built with the `python` feature. The module is exposed as `import molex` via PyO3.

## Installation

```bash
cd crates/molex
maturin develop --release --features python
```

## Core COORDS functions (`python.rs`)

These operate on serialized COORDS bytes:

```python
import molex

# Parse PDB string to COORDS bytes
coords_bytes = molex.pdb_to_coords(pdb_string)

# Parse mmCIF string to COORDS bytes
coords_bytes = molex.mmcif_to_coords(cif_string)

# Convert COORDS bytes to PDB string
pdb_string = molex.coords_to_pdb(coords_bytes)

# Round-trip validation
validated = molex.deserialize_coords_py(coords_bytes)
```

## AtomWorks interop (`adapters::atomworks`)

Entity-aware conversions between molex's ASSEM01 binary format and Biotite `AtomArray` objects with AtomWorks annotations (`entity_id`, `mol_type`, `pn_unit_iid`, `chain_type`).

### molex -> AtomWorks

```python
import molex

# Convert ASSEM01 bytes to AtomArray (basic annotations)
atom_array = molex.entities_to_atom_array(assembly_bytes)

# Convert with full bond list and chain type annotations
atom_array = molex.entities_to_atom_array_plus(assembly_bytes)

# Convert with AtomWorks cleaning pipeline (leaving group removal,
# charge correction, missing atom imputation)
atom_array = molex.entities_to_atom_array_parsed(assembly_bytes, "3nez.cif.gz")
```

### AtomWorks -> molex

```python
# Convert AtomArray back to ASSEM01 bytes (preserves entity metadata)
assembly_bytes = molex.atom_array_to_entities(atom_array)
```

### File-based shortcuts

```python
# Parse structure file directly to ASSEM01 bytes via AtomWorks
assembly_bytes = molex.parse_file_to_entities("3nez.cif.gz")

# Parse to AtomArray with full AtomWorks pipeline
atom_array = molex.parse_file_full("3nez.cif.gz")
```

### Flat Coords functions

Flat Coords-based conversions are available for working with the COORDS01 format directly:

```python
atom_array = molex.coords_to_atom_array(coords_bytes)
atom_array = molex.coords_to_atom_array_plus(coords_bytes)
assembly_bytes = molex.atom_array_to_coords(atom_array)
```

These operate on the flat COORDS01 format and do not include entity metadata.

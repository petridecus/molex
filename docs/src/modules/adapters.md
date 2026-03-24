# Adapters

Format adapters live in `molex::adapters` (source: `src/adapters/`). Each adapter's primary API returns `Vec<MoleculeEntity>`. Coords-based variants are available for FFI/IPC consumers.

## PDB (`adapters::pdb`)

Parses PDB files using the `pdbtbx` crate. Handles non-standard lines (GROMACS/MemProtMD output) via sanitization.

```rust,ignore
// Primary API
pdb_file_to_entities(path: &Path) -> Result<Vec<MoleculeEntity>, CoordsError>
pdb_str_to_entities(pdb_str: &str) -> Result<Vec<MoleculeEntity>, CoordsError>
structure_file_to_entities(path: &Path) -> Result<Vec<MoleculeEntity>, CoordsError>

// Derived Coords API
pdb_file_to_coords(path: &Path) -> Result<Coords, CoordsError>
pdb_str_to_coords(pdb_str: &str) -> Result<Coords, CoordsError>
pdb_to_coords(pdb_str: &str) -> Result<Vec<u8>, CoordsError>  // serialized bytes

// Coords -> PDB
coords_to_pdb(coords_bytes: &[u8]) -> Result<String, CoordsError>
```

`structure_file_to_entities` auto-detects PDB vs mmCIF by file extension (`.pdb`/`.ent` -> PDB, everything else -> mmCIF).

## mmCIF (`adapters::cif`)

Custom DOM-based parser (no external crate). Parses CIF text into a DOM, then extracts atom sites via typed extractors.

```rust,ignore
mmcif_file_to_entities(path: &Path) -> Result<Vec<MoleculeEntity>, CoordsError>
mmcif_str_to_entities(cif_str: &str) -> Result<Vec<MoleculeEntity>, CoordsError>

mmcif_file_to_coords(path: &Path) -> Result<Coords, CoordsError>
mmcif_str_to_coords(cif_str: &str) -> Result<Coords, CoordsError>
```

## BinaryCIF (`adapters::bcif`)

Decodes BinaryCIF (MessagePack-encoded CIF) with column-level codecs.

```rust,ignore
bcif_file_to_entities(path: &Path) -> Result<Vec<MoleculeEntity>, CoordsError>
bcif_to_entities(data: &[u8]) -> Result<Vec<MoleculeEntity>, CoordsError>

bcif_file_to_coords(path: &Path) -> Result<Coords, CoordsError>
bcif_to_coords(data: &[u8]) -> Result<Coords, CoordsError>
```

## MRC/CCP4 density maps (`adapters::mrc`)

Parses MRC/CCP4 format electron density maps into `Density` (which wraps `VoxelGrid`).

```rust,ignore
mrc_file_to_density(path: &Path) -> Result<Density, DensityError>
mrc_to_density(data: &[u8]) -> Result<Density, DensityError>
```

## DCD trajectories (`adapters::dcd`)

Reads DCD binary trajectory files (CHARMM/NAMD format).

```rust,ignore
dcd_file_to_frames(path: &Path) -> Result<Vec<DcdFrame>, CoordsError>

pub struct DcdHeader { /* timestep, n_atoms, etc. */ }
pub struct DcdFrame { pub x: Vec<f32>, pub y: Vec<f32>, pub z: Vec<f32> }
pub struct DcdReader<R> { /* streaming reader */ }
```

## AtomWorks (`adapters::atomworks`, feature = `python`)

Bidirectional conversion between molex entities and Biotite `AtomArray` objects with AtomWorks annotations. Requires the `python` feature.

Entity-aware functions (preserve molecule type, entity ID, chain grouping):

```rust,ignore
entities_to_atom_array(assembly_bytes: Vec<u8>) -> PyResult<PyObject>
entities_to_atom_array_plus(assembly_bytes: Vec<u8>) -> PyResult<PyObject>
atom_array_to_entities(atom_array: PyObject) -> PyResult<Vec<u8>>
entities_to_atom_array_parsed(assembly_bytes: Vec<u8>, filename: &str) -> PyResult<PyObject>
parse_file_to_entities(path: &str) -> PyResult<Vec<u8>>
parse_file_full(path: &str) -> PyResult<PyObject>
```

Flat Coords functions:

```rust,ignore
coords_to_atom_array(py: Python, coords_bytes: Vec<u8>) -> PyResult<PyObject>
coords_to_atom_array_plus(py: Python, coords_bytes: Vec<u8>) -> PyResult<PyObject>
atom_array_to_coords(atom_array: PyObject) -> PyResult<Vec<u8>>
```

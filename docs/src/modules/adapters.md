# Adapters

Format adapters live in `molex::adapters` (source: `src/adapters/`). Each adapter's primary API returns `Vec<MoleculeEntity>`. The legacy `Coords`-returning variants were retired with the COORDS01 wire format; see `docs/COORDS_RETIREMENT_PLAN.md`.

## PDB (`adapters::pdb`)

Parses PDB files using the `pdbtbx` crate. Handles non-standard lines (GROMACS/MemProtMD output) via sanitization.

```rust,ignore
pdb_file_to_entities(path: &Path) -> Result<Vec<MoleculeEntity>, AdapterError>
pdb_str_to_entities(pdb_str: &str) -> Result<Vec<MoleculeEntity>, AdapterError>
structure_file_to_entities(path: &Path) -> Result<Vec<MoleculeEntity>, AdapterError>

// Writers (entity-direct; no Coords intermediate on the public path)
assembly_to_pdb(assembly: &Assembly) -> String
entities_to_pdb<E: Borrow<MoleculeEntity>>(entities: &[E]) -> String
```

`structure_file_to_entities` auto-detects PDB vs mmCIF by file extension (`.pdb`/`.ent` -> PDB, everything else -> mmCIF).

## mmCIF (`adapters::cif`)

Custom DOM-based parser (no external crate). Parses CIF text into a DOM, then extracts atom sites via typed extractors.

```rust,ignore
mmcif_file_to_entities(path: &Path) -> Result<Vec<MoleculeEntity>, AdapterError>
mmcif_str_to_entities(cif_str: &str) -> Result<Vec<MoleculeEntity>, AdapterError>
```

## BinaryCIF (`adapters::bcif`)

Decodes BinaryCIF (MessagePack-encoded CIF) with column-level codecs.

```rust,ignore
bcif_file_to_entities(path: &Path) -> Result<Vec<MoleculeEntity>, AdapterError>
bcif_to_entities(data: &[u8]) -> Result<Vec<MoleculeEntity>, AdapterError>
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
dcd_file_to_frames(path: &Path) -> Result<Vec<DcdFrame>, AdapterError>

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

assembly_bytes_to_atom_array(bytes: Vec<u8>) -> PyResult<PyObject>
assembly_bytes_to_atom_array_plus(bytes: Vec<u8>) -> PyResult<PyObject>
```

The legacy `coords_to_atom_array` / `coords_to_atom_array_plus` / `atom_array_to_coords` functions were renamed to the `assembly_bytes_*` / `atom_array_to_entities` forms when the COORDS01 wire format was retired.

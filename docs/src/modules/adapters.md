# Adapters — Format I/O

The `adapters` module provides parsers and serializers for common
molecular structure file formats.

## PDB / mmCIF (`adapters::pdb`)

The primary entry point for structure files:

```rust,ignore
// Auto-detect format from extension
let coords = structure_file_to_coords("path/to/file.pdb")?;

// Explicit format
let bytes = pdb_to_coords(pdb_string)?;      // PDB → COORDS01 bytes
let bytes = mmcif_to_coords(cif_string)?;     // mmCIF → COORDS01 bytes
let pdb = coords_to_pdb(&coords_bytes)?;      // COORDS01 → PDB string
```

Backed by the `pdbtbx` crate for robust PDB/mmCIF parsing.

## BinaryCIF (`adapters::bcif`)

RCSB's compressed binary CIF format. Handles gzip-compressed files
automatically.

```rust,ignore
let bytes = bcif_to_coords(&bcif_bytes)?;
let bytes = bcif_file_to_coords("path/to/file.bcif")?;
```

## DCD Trajectories (`adapters::dcd`)

CHARMM/NAMD trajectory format — a sequence of coordinate frames
sharing the same topology.

```rust,ignore
let frames: Vec<DcdFrame> = dcd_file_to_frames("trajectory.dcd")?;
// Each frame: Vec<f32> positions (x,y,z interleaved)
```

Also provides streaming via `DcdReader` for large trajectories.

## MRC / CCP4 Density (`adapters::mrc`)

Electron density maps:

```rust,ignore
let density = mrc_to_density(&bytes)?;
let density = mrc_file_to_density("map.mrc")?;
// density: DensityMap with 3D grid + cell parameters
```

## AtomWorks (`adapters::atomworks`)

Bidirectional conversion between molex entities and Biotite
`AtomArray` objects (via PyO3). Used by ML model pipelines
(RF3, RFdiffusion3, LigandMPNN).

```python
import molex

# molex → AtomArray (for model inference)
atom_array = molex.entities_to_atom_array(assembly_bytes)

# AtomArray → molex (after prediction)
assembly_bytes = molex.atom_array_to_entities(atom_array)
```

Feature-gated behind `python`.

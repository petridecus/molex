//! Format adapters for parsing molecular structure files.

pub mod bcif;
pub mod cif;
pub mod dcd;
pub mod mrc;
pub mod pdb;

#[cfg(feature = "python")]
pub mod atomworks;

// Entity-first re-exports (primary API)
// Coords/serialization re-exports (for FFI/IPC consumers)
pub use bcif::{
    bcif_file_to_coords, bcif_file_to_entities, bcif_to_coords,
    bcif_to_entities,
};
pub use cif::{
    mmcif_file_to_coords, mmcif_file_to_entities, mmcif_str_to_coords,
    mmcif_str_to_entities, mmcif_to_coords as mmcif_to_coords_internal,
};
pub use dcd::{dcd_file_to_frames, DcdFrame, DcdHeader, DcdReader};
pub use mrc::{mrc_file_to_density, mrc_to_density};
pub use pdb::{
    coords_to_pdb as coords_bytes_to_pdb, pdb_file_to_coords,
    pdb_file_to_entities, pdb_str_to_coords, pdb_str_to_entities,
    pdb_to_coords as pdb_to_coords_internal, structure_file_to_coords,
    structure_file_to_entities,
};

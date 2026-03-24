//! Binary wire formats, their type definitions, and entity bridge.
//!
//! Two wire formats:
//! - **COORDS01**: flat atom array with element data (backward-compatible with
//!   COORDS00 which omits elements).
//! - **ASSEM01**: entity-aware format that preserves molecule type metadata
//!   alongside atom data.

mod assembly;
mod bridge;
mod deserialize;
mod serialize;
mod types;

use crate::entity::molecule::MoleculeType;

// Wire format constants
/// Magic header bytes identifying the COORDS01 binary format.
pub(crate) const COORDS_MAGIC: &[u8; 8] = b"COORDS01";
pub(crate) const COORDS_MAGIC_V0: &[u8; 8] = b"COORDS00";
/// Magic header bytes identifying the ASSEM01 assembly binary format.
pub const ASSEMBLY_MAGIC: &[u8; 8] = b"ASSEM01\0";

// ASSEM01 molecule type wire encoding
/// Encode a `MoleculeType` to its ASSEM01 wire byte.
pub(crate) fn molecule_type_to_wire(mol_type: MoleculeType) -> u8 {
    match mol_type {
        MoleculeType::Protein => 0,
        MoleculeType::DNA => 1,
        MoleculeType::RNA => 2,
        MoleculeType::Ligand => 3,
        MoleculeType::Ion => 4,
        MoleculeType::Water => 5,
        MoleculeType::Lipid => 6,
        MoleculeType::Cofactor => 7,
        MoleculeType::Solvent => 8,
    }
}

/// Decode an ASSEM01 wire byte to a `MoleculeType`.
pub(crate) fn molecule_type_from_wire(b: u8) -> Option<MoleculeType> {
    match b {
        0 => Some(MoleculeType::Protein),
        1 => Some(MoleculeType::DNA),
        2 => Some(MoleculeType::RNA),
        3 => Some(MoleculeType::Ligand),
        4 => Some(MoleculeType::Ion),
        5 => Some(MoleculeType::Water),
        6 => Some(MoleculeType::Lipid),
        7 => Some(MoleculeType::Cofactor),
        8 => Some(MoleculeType::Solvent),
        _ => None,
    }
}

// Public API re-exports
pub use assembly::{
    assembly_bytes, ca_positions, protein_coords, residue_count,
    update_protein_entities,
};
pub use bridge::{
    coords_to_molecule_entity, extract_atom_set_and_residues, extract_by_type,
    merge_entities, split_into_entities,
};
pub use deserialize::{deserialize, deserialize_assembly};
pub use serialize::{serialize, serialize_assembly};
pub use types::{ChainIdMapper, Coords, CoordsAtom, CoordsError};

#[cfg(test)]
#[path = "tests.rs"]
mod tests;

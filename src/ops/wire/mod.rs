//! ASSEM01 binary wire format: encoder, decoder, header constants.
//!
//! The wire format carries entity headers (molecule type + atom count) plus
//! 26-byte atom rows. See [`crate::ops::wire::serialize_assembly`] for the
//! byte layout.
//!
//! Wire encoding/decoding internally materializes a `Coords` intermediate
//! for the per-atom byte block; that intermediate is private to the crate.

pub(crate) mod deserialize;
pub(crate) mod serialize;

use crate::entity::molecule::{MoleculeEntity, MoleculeType};
use crate::ops::codec::AdapterError;

/// Magic header bytes identifying the ASSEM01 assembly binary format.
pub const ASSEMBLY_MAGIC: &[u8; 8] = b"ASSEM01\0";

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

/// ASSEM01-encode a raw entity slice (includes molecule type metadata).
///
/// # Errors
///
/// Returns `AdapterError` if serialization fails.
pub fn assembly_bytes(
    entities: &[MoleculeEntity],
) -> Result<Vec<u8>, AdapterError> {
    serialize::serialize_entities(entities)
}

pub use deserialize::deserialize_assembly;
pub use serialize::serialize_assembly;

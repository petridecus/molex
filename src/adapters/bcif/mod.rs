//! BinaryCIF (.bcif) format decoder.
//!
//! BinaryCIF is a column-oriented binary encoding of mmCIF used by RCSB PDB
//! as the standard binary format. Files are MessagePack-encoded with optional
//! gzip compression. Spec:
//! <https://github.com/molstar/BinaryCIF/blob/master/encoding.md>.

mod codec;
mod decode;
mod hint;
mod refuse;

use std::path::Path;

use crate::entity::molecule::MoleculeEntity;
use crate::ops::codec::AdapterError;

/// Decode BinaryCIF bytes to entity list.
///
/// Returns the model whose `pdbx_PDB_model_num` matches the smallest value
/// present; files without a model column collapse to a single implicit model.
///
/// # Errors
///
/// Returns [`AdapterError`] if the bytes cannot be parsed as valid
/// BinaryCIF, contain multiple data blocks, or describe more chains than
/// the printable-byte mapper can hold.
pub fn bcif_to_entities(
    bytes: &[u8],
) -> Result<Vec<MoleculeEntity>, AdapterError> {
    decode::decode_to_entities(bytes)
}

/// Load a BinaryCIF file and convert to entity list.
///
/// # Errors
///
/// Returns [`AdapterError`] if the file cannot be read or parsing fails.
pub fn bcif_file_to_entities(
    path: &Path,
) -> Result<Vec<MoleculeEntity>, AdapterError> {
    let bytes = std::fs::read(path).map_err(|e| {
        AdapterError::InvalidFormat(format!("Failed to read file: {e}"))
    })?;
    bcif_to_entities(&bytes)
}

/// Decode BinaryCIF bytes into one entity list per `pdbx_PDB_model_num`.
///
/// Files without a model column collapse to a single bucket.
///
/// # Errors
///
/// Returns [`AdapterError`] if parsing fails.
pub fn bcif_to_all_models(
    bytes: &[u8],
) -> Result<Vec<Vec<MoleculeEntity>>, AdapterError> {
    decode::decode_to_all_models(bytes)
}

/// Load a BinaryCIF file and return one entity list per `pdbx_PDB_model_num`.
///
/// # Errors
///
/// Returns [`AdapterError`] if the file cannot be read or parsing fails.
pub fn bcif_file_to_all_models(
    path: &Path,
) -> Result<Vec<Vec<MoleculeEntity>>, AdapterError> {
    let bytes = std::fs::read(path).map_err(|e| {
        AdapterError::InvalidFormat(format!("Failed to read file: {e}"))
    })?;
    bcif_to_all_models(&bytes)
}

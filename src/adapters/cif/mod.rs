//! mmCIF text-format adapter.
//!
//! Two parser paths feed `EntityBuilder`: a hand-rolled streaming fast
//! scanner and a DOM-backed fallback. The DOM types remain exposed for
//! typed extractors ([`extract::CoordinateData`],
//! [`extract::ReflectionData`]).

pub mod dom;
mod dom_build;
pub mod extract;
mod fast;
mod fast_row;
mod hint;
pub mod parse;
mod refuse;

use std::path::Path;

pub use dom::{Block, ColumnIter, Columns, Document, Loop, RowIter, Value};
pub use extract::{
    AtomSite, CifContent, CoordinateData, ExtractionError, ObsDataType,
    Reflection, ReflectionData, UnitCell,
};
pub use parse::{parse, CifParseError};

use crate::entity::molecule::MoleculeEntity;
use crate::ops::codec::AdapterError;

/// Parse mmCIF format string to entity list.
///
/// Returns the model whose `pdbx_PDB_model_num` matches the smallest
/// value present; files without a model column collapse to a single
/// implicit model.
///
/// # Errors
///
/// Returns [`AdapterError`] if parsing fails.
pub fn mmcif_str_to_entities(
    cif_str: &str,
) -> Result<Vec<MoleculeEntity>, AdapterError> {
    if let Some(result) = fast::parse_mmcif_fast(cif_str) {
        return result;
    }
    dom_build::parse_mmcif_dom_to_entities(cif_str)
}

/// Load mmCIF file to entity list.
///
/// # Errors
///
/// Returns [`AdapterError`] if the file cannot be read or parsing fails.
pub fn mmcif_file_to_entities(
    path: &Path,
) -> Result<Vec<MoleculeEntity>, AdapterError> {
    let content = std::fs::read_to_string(path).map_err(|e| {
        AdapterError::InvalidFormat(format!("Failed to read file: {e}"))
    })?;
    mmcif_str_to_entities(&content)
}

/// Parse an mmCIF string and return one entity list per
/// `pdbx_PDB_model_num`. Files without a model column collapse to a
/// single bucket.
///
/// # Errors
///
/// Returns [`AdapterError`] if parsing fails.
pub fn mmcif_str_to_all_models(
    cif_str: &str,
) -> Result<Vec<Vec<MoleculeEntity>>, AdapterError> {
    if let Some(result) = fast::parse_mmcif_fast_to_all_models(cif_str) {
        return result;
    }
    dom_build::parse_mmcif_dom_to_all_models(cif_str)
}

/// Load an mmCIF file and return one entity list per `pdbx_PDB_model_num`.
///
/// # Errors
///
/// Returns [`AdapterError`] if the file cannot be read or parsing fails.
pub fn mmcif_file_to_all_models(
    path: &Path,
) -> Result<Vec<Vec<MoleculeEntity>>, AdapterError> {
    let content = std::fs::read_to_string(path).map_err(|e| {
        AdapterError::InvalidFormat(format!("Failed to read file: {e}"))
    })?;
    mmcif_str_to_all_models(&content)
}

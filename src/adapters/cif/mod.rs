//! CIF/STAR parser, typed extractors, and mmCIF adapter.

pub mod dom;
pub mod extract;
mod fast;
pub mod parse;

use std::path::Path;

// DOM types
pub use dom::{Block, ColumnIter, Columns, Document, Loop, RowIter, Value};
// Typed extractors
pub use extract::{
    AtomSite, CifContent, CoordinateData, ExtractionError, ObsDataType,
    Reflection, ReflectionData, UnitCell,
};
// Parser
pub use parse::{parse, CifParseError};

use crate::entity::molecule::MoleculeEntity;
use crate::ops::codec::{split_into_entities, AdapterError, Coords};

// ---------------------------------------------------------------------------
// Entity-first API (primary)
// ---------------------------------------------------------------------------

/// Parse mmCIF format string to entity list.
///
/// # Errors
///
/// Returns [`AdapterError`] if parsing fails.
pub fn mmcif_str_to_entities(
    cif_str: &str,
) -> Result<Vec<MoleculeEntity>, AdapterError> {
    parse_mmcif_to_entities(cif_str)
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
    parse_mmcif_to_entities(&content)
}

// ---------------------------------------------------------------------------
// Internal parsing
// ---------------------------------------------------------------------------

/// Parse mmCIF string into entities via temporary Coords + split.
fn parse_mmcif_to_entities(
    input: &str,
) -> Result<Vec<MoleculeEntity>, AdapterError> {
    let coords = parse_mmcif_to_coords(input)?;
    Ok(split_into_entities(&coords))
}

/// Parse mmCIF string into flat Coords (internal).
/// Tries the fast path (direct to Coords, no DOM) first, falls back to DOM.
fn parse_mmcif_to_coords(input: &str) -> Result<Coords, AdapterError> {
    if let Some(result) = fast::parse_mmcif_fast(input) {
        return result;
    }
    parse_mmcif_to_coords_dom(input)
}

/// DOM-based mmCIF parsing fallback.
fn parse_mmcif_to_coords_dom(input: &str) -> Result<Coords, AdapterError> {
    use crate::element::Element;
    use crate::ops::codec::{ChainIdMapper, CoordsAtom};

    let doc = parse(input).map_err(|e| {
        AdapterError::InvalidFormat(format!("CIF parse error: {e}"))
    })?;
    let block = doc.blocks.first().ok_or_else(|| {
        AdapterError::InvalidFormat("No data blocks in CIF".into())
    })?;
    let data = CoordinateData::try_from(block).map_err(|e| {
        AdapterError::InvalidFormat(format!("CIF extraction error: {e}"))
    })?;

    if data.atoms.is_empty() {
        return Err(AdapterError::InvalidFormat(
            "No atoms found in CIF".into(),
        ));
    }

    let n = data.atoms.len();
    let mut atoms = Vec::with_capacity(n);
    let mut chain_ids = Vec::with_capacity(n);
    let mut res_names = Vec::with_capacity(n);
    let mut res_nums = Vec::with_capacity(n);
    let mut atom_names = Vec::with_capacity(n);
    let mut elements = Vec::with_capacity(n);
    let mut chain_mapper = ChainIdMapper::new();

    #[allow(clippy::cast_possible_truncation)]
    for a in &data.atoms {
        atoms.push(CoordsAtom {
            x: a.x as f32,
            y: a.y as f32,
            z: a.z as f32,
            occupancy: a.occupancy as f32,
            b_factor: a.b_factor as f32,
        });
        chain_ids.push(chain_mapper.get_or_assign(&a.chain));
        res_names.push(name_to_bytes::<3>(&a.residue));
        res_nums.push(a.seq_id.unwrap_or(0));
        atom_names.push(name_to_bytes::<4>(&a.label));
        let elem = if a.element.is_empty() {
            Element::from_atom_name(&a.label)
        } else {
            Element::from_symbol(&a.element)
        };
        elements.push(elem);
    }

    Ok(Coords {
        num_atoms: atoms.len(),
        atoms,
        chain_ids,
        res_names,
        res_nums,
        atom_names,
        elements,
    })
}

/// Convert a string to a space-padded byte array of length N.
fn name_to_bytes<const N: usize>(name: &str) -> [u8; N] {
    let mut buf = [b' '; N];
    for (i, b) in name.bytes().take(N).enumerate() {
        buf[i] = b;
    }
    buf
}

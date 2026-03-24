//! BinaryCIF (.bcif) format decoder.
//!
//! BinaryCIF is a column-oriented binary encoding of mmCIF, used by RCSB PDB
//! as the standard binary format (replacing MMTF). Files are
//! MessagePack-encoded with optional gzip compression.
//!
//! Reference: <https://github.com/molstar/BinaryCIF/blob/master/encoding.md>

mod codec;

use std::collections::HashMap;
use std::io::Read;
use std::path::Path;

use codec::{decode_column, decode_msgpack, ColData, MsgVal};

use crate::element::Element;
use crate::entity::molecule::MoleculeEntity;
use crate::ops::codec::{
    merge_entities, split_into_entities, ChainIdMapper, Coords, CoordsAtom,
    CoordsError,
};

// ---------------------------------------------------------------------------
// Entity-first API (primary)
// ---------------------------------------------------------------------------

/// Decode BinaryCIF bytes to entity list.
///
/// # Errors
///
/// Returns `CoordsError::InvalidFormat` if the bytes cannot be parsed as
/// valid BinaryCIF.
pub fn bcif_to_entities(
    bytes: &[u8],
) -> Result<Vec<MoleculeEntity>, CoordsError> {
    let coords = parse_bcif_to_coords(bytes)?;
    Ok(split_into_entities(&coords))
}

/// Load a BinaryCIF file and convert to entity list.
///
/// # Errors
///
/// Returns `CoordsError::InvalidFormat` if the file cannot be read or parsed
/// as valid BinaryCIF.
pub fn bcif_file_to_entities(
    path: &Path,
) -> Result<Vec<MoleculeEntity>, CoordsError> {
    let bytes = std::fs::read(path).map_err(|e| {
        CoordsError::InvalidFormat(format!("Failed to read file: {e}"))
    })?;
    bcif_to_entities(&bytes)
}

// ---------------------------------------------------------------------------
// Coords/serialization API (derived from entities)
// ---------------------------------------------------------------------------

/// Load a BinaryCIF file and convert to Coords.
///
/// # Errors
///
/// Returns `CoordsError::InvalidFormat` if the file cannot be read,
/// decompressed, or parsed as valid BinaryCIF.
pub fn bcif_file_to_coords(path: &Path) -> Result<Coords, CoordsError> {
    let entities = bcif_file_to_entities(path)?;
    Ok(merge_entities(&entities))
}

/// Decode BinaryCIF bytes (possibly gzipped) into a Coords struct.
///
/// # Errors
///
/// Returns `CoordsError::InvalidFormat` if the bytes cannot be decompressed
/// or parsed as valid BinaryCIF.
pub fn bcif_to_coords(bytes: &[u8]) -> Result<Coords, CoordsError> {
    let entities = bcif_to_entities(bytes)?;
    Ok(merge_entities(&entities))
}

// ---------------------------------------------------------------------------
// Internal parsing
// ---------------------------------------------------------------------------

fn decompress_if_gzip(bytes: &[u8]) -> Result<Vec<u8>, CoordsError> {
    if bytes.len() >= 2 && bytes[0] == 0x1f && bytes[1] == 0x8b {
        let mut decoder = flate2::read::GzDecoder::new(bytes);
        let mut out = Vec::new();
        let _bytes_read = decoder.read_to_end(&mut out).map_err(|e| {
            CoordsError::InvalidFormat(format!(
                "Gzip decompression failed: {e}"
            ))
        })?;
        Ok(out)
    } else {
        Ok(bytes.to_vec())
    }
}

/// Parse BinaryCIF bytes into flat Coords (internal).
#[allow(
    clippy::too_many_lines,
    reason = "coordinate extraction from BinaryCIF requires processing many \
              columns"
)]
fn parse_bcif_to_coords(raw_bytes: &[u8]) -> Result<Coords, CoordsError> {
    let data = decompress_if_gzip(raw_bytes)?;
    let root = decode_msgpack(&data)?;

    let data_blocks = root
        .get("dataBlocks")
        .and_then(MsgVal::as_array)
        .ok_or_else(|| {
            CoordsError::InvalidFormat("Missing 'dataBlocks'".into())
        })?;

    if data_blocks.is_empty() {
        return Err(CoordsError::InvalidFormat("No data blocks found".into()));
    }

    let block = &data_blocks[0];
    let categories = block
        .get("categories")
        .and_then(MsgVal::as_array)
        .ok_or_else(|| {
            CoordsError::InvalidFormat(
                "Missing 'categories' in data block".into(),
            )
        })?;

    let atom_site = categories
        .iter()
        .find(|cat| {
            cat.get("name").and_then(MsgVal::as_str) == Some("_atom_site")
        })
        .ok_or_else(|| {
            CoordsError::InvalidFormat("No '_atom_site' category found".into())
        })?;

    #[allow(
        clippy::cast_possible_truncation,
        reason = "row count fits in usize"
    )]
    #[allow(clippy::cast_sign_loss, reason = "row count is non-negative")]
    let row_count = atom_site
        .get("rowCount")
        .and_then(MsgVal::as_i64)
        .ok_or_else(|| {
            CoordsError::InvalidFormat("Missing 'rowCount'".into())
        })? as usize;

    let columns = atom_site
        .get("columns")
        .and_then(MsgVal::as_array)
        .ok_or_else(|| {
            CoordsError::InvalidFormat("Missing 'columns'".into())
        })?;

    let col_map: HashMap<&str, &MsgVal> = columns
        .iter()
        .filter_map(|col| {
            let name = col.get("name")?.as_str()?;
            Some((name, col))
        })
        .collect();

    let cartn_x = decode_float_col(&col_map, "Cartn_x", row_count)?;
    let cartn_y = decode_float_col(&col_map, "Cartn_y", row_count)?;
    let cartn_z = decode_float_col(&col_map, "Cartn_z", row_count)?;
    let label_atom_id =
        decode_string_col(&col_map, "label_atom_id", row_count)?;
    let label_comp_id =
        decode_string_col(&col_map, "label_comp_id", row_count)?;
    let label_asym_id =
        decode_string_col(&col_map, "label_asym_id", row_count)?;
    let label_seq_id = decode_int_col(&col_map, "label_seq_id", row_count)?;

    let occupancy = decode_float_col_opt(&col_map, "occupancy", row_count, 1.0);
    let b_factor =
        decode_float_col_opt(&col_map, "B_iso_or_equiv", row_count, 0.0);

    let type_symbol =
        decode_string_col(&col_map, "type_symbol", row_count).ok();

    let mut atoms = Vec::with_capacity(row_count);
    let mut chain_ids = Vec::with_capacity(row_count);
    let mut res_names = Vec::with_capacity(row_count);
    let mut res_nums = Vec::with_capacity(row_count);
    let mut atom_names = Vec::with_capacity(row_count);
    let mut elements = Vec::with_capacity(row_count);
    let mut chain_mapper = ChainIdMapper::new();

    #[allow(
        clippy::cast_possible_truncation,
        reason = "f64→f32 truncation is intentional for coordinate storage"
    )]
    for i in 0..row_count {
        atoms.push(CoordsAtom {
            x: cartn_x[i] as f32,
            y: cartn_y[i] as f32,
            z: cartn_z[i] as f32,
            occupancy: occupancy[i] as f32,
            b_factor: b_factor[i] as f32,
        });

        chain_ids.push(chain_mapper.get_or_assign(&label_asym_id[i]));

        let mut rn = [b' '; 3];
        for (j, b) in label_comp_id[i].bytes().take(3).enumerate() {
            rn[j] = b;
        }
        res_names.push(rn);

        res_nums.push(label_seq_id[i]);

        let mut an = [b' '; 4];
        for (j, b) in label_atom_id[i].bytes().take(4).enumerate() {
            an[j] = b;
        }
        atom_names.push(an);

        let elem = type_symbol.as_ref().map_or_else(
            || Element::from_atom_name(&label_atom_id[i]),
            |ts| Element::from_symbol(&ts[i]),
        );
        elements.push(elem);
    }

    if atoms.is_empty() {
        return Err(CoordsError::InvalidFormat(
            "No ATOM records found in BinaryCIF".into(),
        ));
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

// ---------------------------------------------------------------------------
// Column extraction helpers
// ---------------------------------------------------------------------------

fn get_col_data<'a>(
    col_map: &HashMap<&str, &'a MsgVal>,
    name: &str,
) -> Result<&'a MsgVal, CoordsError> {
    let col = col_map.get(name).ok_or_else(|| {
        CoordsError::InvalidFormat(format!("Missing column '{name}'"))
    })?;
    col.get("data").ok_or_else(|| {
        CoordsError::InvalidFormat(format!("Column '{name}' has no data"))
    })
}

fn decode_float_col(
    col_map: &HashMap<&str, &MsgVal>,
    name: &str,
    expected: usize,
) -> Result<Vec<f64>, CoordsError> {
    let data = get_col_data(col_map, name)?;
    match decode_column(data)? {
        ColData::FloatArray(v) if v.len() == expected => Ok(v),
        ColData::FloatArray(v) => Err(CoordsError::InvalidFormat(format!(
            "Column '{name}': expected {expected} rows, got {}",
            v.len()
        ))),
        ColData::IntArray(v) if v.len() == expected => {
            Ok(v.iter().map(|&x| f64::from(x)).collect())
        }
        _ => Err(CoordsError::InvalidFormat(format!(
            "Column '{name}': expected float array"
        ))),
    }
}

fn decode_float_col_opt(
    col_map: &HashMap<&str, &MsgVal>,
    name: &str,
    count: usize,
    default: f64,
) -> Vec<f64> {
    decode_float_col(col_map, name, count)
        .unwrap_or_else(|_| vec![default; count])
}

fn decode_int_col(
    col_map: &HashMap<&str, &MsgVal>,
    name: &str,
    expected: usize,
) -> Result<Vec<i32>, CoordsError> {
    let data = get_col_data(col_map, name)?;
    match decode_column(data)? {
        ColData::IntArray(v) if v.len() == expected => Ok(v),
        ColData::IntArray(v) => Err(CoordsError::InvalidFormat(format!(
            "Column '{name}': expected {expected} rows, got {}",
            v.len()
        ))),
        _ => Err(CoordsError::InvalidFormat(format!(
            "Column '{name}': expected int array"
        ))),
    }
}

fn decode_string_col(
    col_map: &HashMap<&str, &MsgVal>,
    name: &str,
    expected: usize,
) -> Result<Vec<String>, CoordsError> {
    let data = get_col_data(col_map, name)?;
    match decode_column(data)? {
        ColData::StringArray(v) if v.len() == expected => Ok(v),
        ColData::StringArray(v) => Err(CoordsError::InvalidFormat(format!(
            "Column '{name}': expected {expected} rows, got {}",
            v.len()
        ))),
        _ => Err(CoordsError::InvalidFormat(format!(
            "Column '{name}': expected string array"
        ))),
    }
}

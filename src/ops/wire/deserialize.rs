//! Deserialization for the ASSEM01 binary wire format.

use super::{molecule_type_from_wire, ASSEMBLY_MAGIC};
use crate::element::Element;
use crate::entity::molecule::id::EntityIdAllocator;
use crate::entity::molecule::{MoleculeEntity, MoleculeType};
use crate::ops::codec::{
    coords_to_molecule_entity, AdapterError, Coords, CoordsAtom,
};

/// One atom's worth of parsed ASSEM01 binary data.
struct ParsedAtom<'a> {
    atom: CoordsAtom,
    chain_id: u8,
    res_name: [u8; 3],
    res_num: i32,
    atom_name: [u8; 4],
    element: Element,
    rest: &'a [u8],
}

/// Read one 26-byte ASSEM01 atom row from a cursor slice.
fn read_atom_from_cursor(
    cursor: &[u8],
) -> Result<ParsedAtom<'_>, AdapterError> {
    let x = f32::from_be_bytes(cursor[0..4].try_into().map_err(|_| {
        AdapterError::SerializationError("Invalid x coordinate".to_owned())
    })?);
    let y = f32::from_be_bytes(cursor[4..8].try_into().map_err(|_| {
        AdapterError::SerializationError("Invalid y coordinate".to_owned())
    })?);
    let z = f32::from_be_bytes(cursor[8..12].try_into().map_err(|_| {
        AdapterError::SerializationError("Invalid z coordinate".to_owned())
    })?);
    let atom = CoordsAtom {
        x,
        y,
        z,
        occupancy: 1.0,
        b_factor: 0.0,
    };
    let rest = &cursor[12..];

    let chain_id = rest[0];
    let rest = &rest[1..];

    let mut res_name = [0u8; 3];
    res_name.copy_from_slice(&rest[0..3]);
    let rest = &rest[3..];

    let res_num = i32::from_be_bytes(rest[0..4].try_into().map_err(|_| {
        AdapterError::SerializationError("Invalid residue number".to_owned())
    })?);
    let rest = &rest[4..];

    let mut atom_name = [0u8; 4];
    atom_name.copy_from_slice(&rest[0..4]);
    let rest = &rest[4..];

    let sym_str = std::str::from_utf8(&rest[0..2])
        .unwrap_or("")
        .trim_matches('\0')
        .trim();
    let element = Element::from_symbol(sym_str);
    let rest = &rest[2..];

    Ok(ParsedAtom {
        atom,
        chain_id,
        res_name,
        res_num,
        atom_name,
        element,
        rest,
    })
}

/// Read `num_atoms` 26-byte rows from a cursor, returning a `Coords` and
/// the remaining bytes.
fn read_atoms_to_coords(
    mut cursor: &[u8],
    num_atoms: usize,
) -> Result<(Coords, &[u8]), AdapterError> {
    let mut atoms = Vec::with_capacity(num_atoms);
    let mut chain_ids = Vec::with_capacity(num_atoms);
    let mut res_names = Vec::with_capacity(num_atoms);
    let mut res_nums = Vec::with_capacity(num_atoms);
    let mut atom_names = Vec::with_capacity(num_atoms);
    let mut elements = Vec::with_capacity(num_atoms);

    for _ in 0..num_atoms {
        let parsed = read_atom_from_cursor(cursor)?;
        atoms.push(parsed.atom);
        chain_ids.push(parsed.chain_id);
        res_names.push(parsed.res_name);
        res_nums.push(parsed.res_num);
        atom_names.push(parsed.atom_name);
        elements.push(parsed.element);
        cursor = parsed.rest;
    }

    let coords = Coords {
        num_atoms,
        atoms,
        chain_ids,
        res_names,
        res_nums,
        atom_names,
        elements,
    };
    Ok((coords, cursor))
}

/// Parse entity headers from ASSEM01 binary, returning
/// (molecule_type, atom_count) pairs and the offset past all headers.
fn parse_entity_headers(
    bytes: &[u8],
    entity_count: usize,
) -> Result<(Vec<(MoleculeType, usize)>, usize), AdapterError> {
    let mut headers = Vec::with_capacity(entity_count);
    let mut offset = 12;
    for _ in 0..entity_count {
        let mol_type =
            molecule_type_from_wire(bytes[offset]).ok_or_else(|| {
                AdapterError::InvalidFormat(format!(
                    "Unknown molecule type byte: {}",
                    bytes[offset]
                ))
            })?;
        offset += 1;
        let atom_count = u32::from_be_bytes(
            bytes[offset..offset + 4].try_into().map_err(|_| {
                AdapterError::InvalidFormat(
                    "Invalid atom count in entity header".to_owned(),
                )
            })?,
        ) as usize;
        offset += 4;
        headers.push((mol_type, atom_count));
    }
    Ok((headers, offset))
}

/// Deserialize ASSEM01 binary format into an [`Assembly`] with
/// derived data populated.
///
/// Runs [`Assembly::new`] over the decoded entities so callers see
/// `ss_types`, `hbonds`, and `cross_entity_bonds` without a follow-up
/// step.
///
/// # Errors
///
/// Returns `AdapterError::InvalidFormat` if the magic header, entity headers,
/// or atom data are malformed or truncated. Returns
/// `AdapterError::SerializationError` if individual atom fields cannot be
/// parsed.
///
/// [`Assembly`]: crate::Assembly
/// [`Assembly::new`]: crate::Assembly::new
pub fn deserialize_assembly(
    bytes: &[u8],
) -> Result<crate::Assembly, AdapterError> {
    let entities = deserialize_assembly_entities(bytes)?;
    Ok(crate::Assembly::new(entities))
}

/// Deserialize ASSEM01 bytes into a raw entity vec without building
/// an [`Assembly`].
///
/// Internal helper for in-crate paths that mutate entities in place
/// before re-serializing and don't need the derived-data recompute
/// [`Assembly::new`] performs.
///
/// # Errors
///
/// Same as [`deserialize_assembly`].
///
/// [`Assembly`]: crate::Assembly
/// [`Assembly::new`]: crate::Assembly::new
pub(crate) fn deserialize_assembly_entities(
    bytes: &[u8],
) -> Result<Vec<MoleculeEntity>, AdapterError> {
    if bytes.len() < 12 {
        return Err(AdapterError::InvalidFormat(
            "Data too short for ASSEM01 header".to_owned(),
        ));
    }

    let magic = &bytes[0..8];
    if magic != ASSEMBLY_MAGIC {
        return Err(AdapterError::InvalidFormat(
            "Invalid magic number for ASSEM01".to_owned(),
        ));
    }

    let entity_count =
        u32::from_be_bytes(bytes[8..12].try_into().map_err(|_| {
            AdapterError::InvalidFormat("Invalid entity count".to_owned())
        })?) as usize;

    let (entity_headers, headers_end) =
        parse_entity_headers(bytes, entity_count)?;

    let total_atoms: usize = entity_headers.iter().map(|(_, c)| c).sum();
    if bytes.len() < headers_end + total_atoms * 26 {
        return Err(AdapterError::InvalidFormat(
            "Data too short for atom data".to_owned(),
        ));
    }

    let mut cursor = &bytes[headers_end..];
    let mut entities = Vec::with_capacity(entity_count);
    let mut allocator = EntityIdAllocator::new();

    for (mol_type, atom_count) in entity_headers {
        let (coords, rest) = read_atoms_to_coords(cursor, atom_count)?;
        cursor = rest;
        let id = allocator.allocate();
        entities.push(coords_to_molecule_entity(id, mol_type, &coords));
    }

    Ok(entities)
}

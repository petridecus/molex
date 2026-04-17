//! Deserialization for COORDS01/COORDS00 and ASSEM01 binary formats.

use super::{
    coords_to_molecule_entity, molecule_type_from_wire, Coords, CoordsAtom,
    CoordsError, ASSEMBLY_MAGIC, COORDS_MAGIC, COORDS_MAGIC_V0,
};
use crate::element::Element;
use crate::entity::molecule::id::EntityIdAllocator;
use crate::entity::molecule::{MoleculeEntity, MoleculeType};

/// One atom's worth of parsed COORDS binary data.
struct ParsedAtom<'a> {
    atom: CoordsAtom,
    chain_id: u8,
    res_name: [u8; 3],
    res_num: i32,
    atom_name: [u8; 4],
    element: Element,
    rest: &'a [u8],
}

/// Read one atom's worth of COORDS binary data from a cursor slice.
/// `has_elements` controls whether a 2-byte element symbol follows
/// (COORDS01) or elements are inferred from atom name (COORDS00).
fn read_atom_from_cursor(
    cursor: &[u8],
    has_elements: bool,
) -> Result<ParsedAtom<'_>, CoordsError> {
    let x = f32::from_be_bytes(cursor[0..4].try_into().map_err(|_| {
        CoordsError::SerializationError("Invalid x coordinate".to_owned())
    })?);
    let y = f32::from_be_bytes(cursor[4..8].try_into().map_err(|_| {
        CoordsError::SerializationError("Invalid y coordinate".to_owned())
    })?);
    let z = f32::from_be_bytes(cursor[8..12].try_into().map_err(|_| {
        CoordsError::SerializationError("Invalid z coordinate".to_owned())
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
        CoordsError::SerializationError("Invalid residue number".to_owned())
    })?);
    let rest = &rest[4..];

    let mut atom_name = [0u8; 4];
    atom_name.copy_from_slice(&rest[0..4]);
    let rest = &rest[4..];

    let (elem, rest) = if has_elements {
        let sym_str = std::str::from_utf8(&rest[0..2])
            .unwrap_or("")
            .trim_matches('\0')
            .trim();
        (Element::from_symbol(sym_str), &rest[2..])
    } else {
        let aname = std::str::from_utf8(&atom_name).unwrap_or("");
        (Element::from_atom_name(aname), rest)
    };

    Ok(ParsedAtom {
        atom,
        chain_id,
        res_name,
        res_num,
        atom_name,
        element: elem,
        rest,
    })
}

/// Read `num_atoms` atoms from a cursor, returning a `Coords` struct and the
/// remaining bytes.
fn read_atoms_to_coords(
    mut cursor: &[u8],
    num_atoms: usize,
    has_elements: bool,
) -> Result<(Coords, &[u8]), CoordsError> {
    let mut atoms = Vec::with_capacity(num_atoms);
    let mut chain_ids = Vec::with_capacity(num_atoms);
    let mut res_names = Vec::with_capacity(num_atoms);
    let mut res_nums = Vec::with_capacity(num_atoms);
    let mut atom_names = Vec::with_capacity(num_atoms);
    let mut elements = Vec::with_capacity(num_atoms);

    for _ in 0..num_atoms {
        let parsed = read_atom_from_cursor(cursor, has_elements)?;
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

/// Deserialize COORDS binary format to `Coords` struct.
/// Supports both COORDS00 (no element data) and COORDS01 (with element data).
///
/// # Errors
///
/// Returns `CoordsError::InvalidFormat` if the magic header is wrong or data is
/// truncated. Returns `CoordsError::SerializationError` if individual atom
/// fields cannot be parsed.
pub fn deserialize(coords_bytes: &[u8]) -> Result<Coords, CoordsError> {
    if coords_bytes.len() < 8 {
        return Err(CoordsError::InvalidFormat(
            "Data too short to be valid COORDS".to_owned(),
        ));
    }

    let magic = &coords_bytes[0..8];
    let has_elements = magic == COORDS_MAGIC;
    if magic != COORDS_MAGIC && magic != COORDS_MAGIC_V0 {
        return Err(CoordsError::InvalidFormat(
            "Invalid magic number in COORDS header".to_owned(),
        ));
    }

    let cursor = &coords_bytes[8..];
    let num_atoms = u32::from_be_bytes(
        cursor
            .get(0..4)
            .ok_or_else(|| {
                CoordsError::InvalidFormat("Missing num_atoms field".to_owned())
            })?
            .try_into()
            .map_err(|_| {
                CoordsError::SerializationError(
                    "Invalid num_atoms size".to_owned(),
                )
            })?,
    ) as usize;

    let per_atom = if has_elements { 26 } else { 24 };
    if cursor.len() - 4 < num_atoms * per_atom {
        return Err(CoordsError::InvalidFormat(
            "Data too short for declared number of atoms".to_owned(),
        ));
    }

    let (coords, _) =
        read_atoms_to_coords(&cursor[4..], num_atoms, has_elements)?;
    Ok(coords)
}

/// Parse entity headers from ASSEM01 binary, returning
/// (molecule_type, atom_count) pairs and the offset past all headers.
fn parse_entity_headers(
    bytes: &[u8],
    entity_count: usize,
) -> Result<(Vec<(MoleculeType, usize)>, usize), CoordsError> {
    let mut headers = Vec::with_capacity(entity_count);
    let mut offset = 12;
    for _ in 0..entity_count {
        let mol_type =
            molecule_type_from_wire(bytes[offset]).ok_or_else(|| {
                CoordsError::InvalidFormat(format!(
                    "Unknown molecule type byte: {}",
                    bytes[offset]
                ))
            })?;
        offset += 1;
        let atom_count = u32::from_be_bytes(
            bytes[offset..offset + 4].try_into().map_err(|_| {
                CoordsError::InvalidFormat(
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
/// Returns `CoordsError::InvalidFormat` if the magic header, entity headers,
/// or atom data are malformed or truncated. Returns
/// `CoordsError::SerializationError` if individual atom fields cannot be
/// parsed.
///
/// [`Assembly`]: crate::Assembly
/// [`Assembly::new`]: crate::Assembly::new
pub fn deserialize_assembly(
    bytes: &[u8],
) -> Result<crate::Assembly, CoordsError> {
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
) -> Result<Vec<MoleculeEntity>, CoordsError> {
    if bytes.len() < 12 {
        return Err(CoordsError::InvalidFormat(
            "Data too short for ASSEM01 header".to_owned(),
        ));
    }

    let magic = &bytes[0..8];
    if magic != ASSEMBLY_MAGIC {
        return Err(CoordsError::InvalidFormat(
            "Invalid magic number for ASSEM01".to_owned(),
        ));
    }

    let entity_count =
        u32::from_be_bytes(bytes[8..12].try_into().map_err(|_| {
            CoordsError::InvalidFormat("Invalid entity count".to_owned())
        })?) as usize;

    let (entity_headers, headers_end) =
        parse_entity_headers(bytes, entity_count)?;

    let total_atoms: usize = entity_headers.iter().map(|(_, c)| c).sum();
    if bytes.len() < headers_end + total_atoms * 26 {
        return Err(CoordsError::InvalidFormat(
            "Data too short for atom data".to_owned(),
        ));
    }

    let mut cursor = &bytes[headers_end..];
    let mut entities = Vec::with_capacity(entity_count);
    let mut allocator = EntityIdAllocator::new();

    for (mol_type, atom_count) in entity_headers {
        let (coords, rest) = read_atoms_to_coords(cursor, atom_count, true)?;
        cursor = rest;
        let id = allocator.allocate();
        entities.push(coords_to_molecule_entity(id, mol_type, &coords));
    }

    Ok(entities)
}

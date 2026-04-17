//! Serialization for COORDS01 and ASSEM01 binary formats.

use super::{
    molecule_type_to_wire, Coords, CoordsError, ASSEMBLY_MAGIC, COORDS_MAGIC,
};
use crate::assembly::Assembly;
use crate::entity::molecule::MoleculeEntity;

/// Serialize `Coords` struct to COORDS01 binary format.
///
/// # Errors
///
/// Currently infallible but returns `Result` for API consistency.
pub fn serialize(coords: &Coords) -> Result<Vec<u8>, CoordsError> {
    let mut buffer =
        Vec::with_capacity(8 + 4 + coords.num_atoms * (12 + 1 + 3 + 4 + 4 + 2));

    buffer.extend_from_slice(COORDS_MAGIC);

    #[allow(clippy::cast_possible_truncation)] // atom counts fit in u32
    let num_atoms_u32 = coords.num_atoms as u32;
    buffer.extend_from_slice(&num_atoms_u32.to_be_bytes());

    for i in 0..coords.num_atoms {
        write_atom(&mut buffer, coords, i);
    }

    Ok(buffer)
}

/// Serialize an [`Assembly`] to ASSEM01 binary format.
///
/// Writes the entity list only; derived fields (`ss_types`, `hbonds`,
/// `cross_entity_bonds`, `generation`) are rebuilt by
/// [`deserialize_assembly`] via [`Assembly::new`].
///
/// Format:
/// - 8 bytes: magic "ASSEM01\0"
/// - 4 bytes: entity_count (u32 BE)
/// - Per entity header (5 bytes each):
///   - 1 byte: `molecule_type` wire byte
///   - 4 bytes: `atom_count` (`u32` BE)
/// - Per atom (26 bytes, same as COORDS01):
///   - 12 bytes: x,y,z (`f32` BE)
///   - 1 byte: `chain_id`
///   - 3 bytes: `res_name`
///   - 4 bytes: `res_num` (`i32` BE)
///   - 4 bytes: `atom_name`
///   - 2 bytes: element symbol
///
/// # Errors
///
/// Currently infallible but returns `Result` for API consistency.
pub fn serialize_assembly(
    assembly: &Assembly,
) -> Result<Vec<u8>, CoordsError> {
    serialize_entities(assembly.entities())
}

/// Serialize a raw entity slice to ASSEM01 binary format.
///
/// Internal helper for in-crate paths that operate on owned entity
/// slices and don't need an [`Assembly`] wrapper.
#[allow(
    clippy::unnecessary_wraps,
    reason = "mirrors `serialize_assembly` which keeps the `Result` for \
              API consistency"
)]
pub(crate) fn serialize_entities(
    entities: &[MoleculeEntity],
) -> Result<Vec<u8>, CoordsError> {
    let total_atoms: usize =
        entities.iter().map(MoleculeEntity::atom_count).sum();
    let header_size = 8 + 4 + entities.len() * 5;
    let atom_size = total_atoms * 26;
    let mut buffer = Vec::with_capacity(header_size + atom_size);

    // Magic
    buffer.extend_from_slice(ASSEMBLY_MAGIC);

    // Entity count
    #[allow(clippy::cast_possible_truncation)] // entity count fits in u32
    buffer.extend_from_slice(&(entities.len() as u32).to_be_bytes());

    // Per-entity headers
    for entity in entities {
        buffer.push(molecule_type_to_wire(entity.molecule_type()));
        #[allow(clippy::cast_possible_truncation)] // atom count fits in u32
        buffer.extend_from_slice(&(entity.atom_count() as u32).to_be_bytes());
    }

    // Atom data (same layout as COORDS01)
    for entity in entities {
        let c = entity.to_coords();
        for i in 0..c.num_atoms {
            write_atom(&mut buffer, &c, i);
        }
    }

    Ok(buffer)
}

/// Write a single atom's data in the COORDS01 wire layout (26 bytes).
fn write_atom(buffer: &mut Vec<u8>, coords: &Coords, i: usize) {
    let atom = &coords.atoms[i];
    buffer.extend_from_slice(&atom.x.to_be_bytes());
    buffer.extend_from_slice(&atom.y.to_be_bytes());
    buffer.extend_from_slice(&atom.z.to_be_bytes());
    buffer.push(coords.chain_ids[i]);
    buffer.extend_from_slice(&coords.res_names[i]);
    buffer.extend_from_slice(&coords.res_nums[i].to_be_bytes());
    buffer.extend_from_slice(&coords.atom_names[i]);
    let sym = coords.elements.get(i).map_or("X", |e| e.symbol());
    let sym_bytes = sym.as_bytes();
    buffer.push(sym_bytes.first().copied().unwrap_or(b'X'));
    buffer.push(sym_bytes.get(1).copied().unwrap_or(0));
}

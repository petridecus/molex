//! Serialization for the ASSEM01 binary wire format.

use super::{molecule_type_to_wire, ASSEMBLY_MAGIC};
use crate::assembly::Assembly;
use crate::entity::molecule::atom::Atom;
use crate::entity::molecule::polymer::Residue;
use crate::entity::molecule::MoleculeEntity;
use crate::ops::codec::AdapterError;

/// Serialize an [`Assembly`] to ASSEM01 binary format.
///
/// Writes the entity list only; derived fields (`ss_types`, `hbonds`,
/// `cross_entity_bonds`, `generation`) are rebuilt by
/// [`deserialize_assembly`](super::deserialize::deserialize_assembly) via
/// [`Assembly::new`].
///
/// Format:
/// - 8 bytes: magic "ASSEM01\0"
/// - 4 bytes: entity_count (u32 BE)
/// - Per entity header (5 bytes each):
///   - 1 byte: `molecule_type` wire byte
///   - 4 bytes: `atom_count` (`u32` BE)
/// - Per atom (26 bytes):
///   - 12 bytes: x, y, z (`f32` BE × 3)
///   - 1 byte:   `chain_id`
///   - 3 bytes:  `res_name`
///   - 4 bytes:  `res_num` (`i32` BE)
///   - 4 bytes:  `atom_name`
///   - 2 bytes:  element symbol (byte 0, byte 1 or 0)
///
/// Occupancy and b_factor are not preserved on the wire; deserialize
/// resets them to 1.0 and 0.0 respectively.
///
/// # Errors
///
/// Currently infallible but returns `Result` for API consistency.
pub fn serialize_assembly(
    assembly: &Assembly,
) -> Result<Vec<u8>, AdapterError> {
    serialize_entities(assembly.entities())
}

/// Serialize a raw entity slice to ASSEM01 binary format.
///
/// Internal helper for in-crate paths that operate on owned entity
/// slices and don't need an [`Assembly`] wrapper.
#[allow(
    clippy::unnecessary_wraps,
    reason = "mirrors `serialize_assembly` which keeps the `Result` for API \
              consistency"
)]
pub(crate) fn serialize_entities<E: std::borrow::Borrow<MoleculeEntity>>(
    entities: &[E],
) -> Result<Vec<u8>, AdapterError> {
    let total_atoms: usize =
        entities.iter().map(|e| e.borrow().atom_count()).sum();
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
        let entity = entity.borrow();
        buffer.push(molecule_type_to_wire(entity.molecule_type()));
        #[allow(clippy::cast_possible_truncation)] // atom count fits in u32
        buffer.extend_from_slice(&(entity.atom_count() as u32).to_be_bytes());
    }

    // Atom data per entity, walking residues directly.
    for entity in entities {
        write_entity_atoms(entity.borrow(), &mut buffer);
    }

    Ok(buffer)
}

fn write_entity_atoms(entity: &MoleculeEntity, buffer: &mut Vec<u8>) {
    match entity {
        MoleculeEntity::Protein(e) => {
            write_polymer_atoms(&e.atoms, &e.residues, e.pdb_chain_id, buffer);
        }
        MoleculeEntity::NucleicAcid(e) => {
            write_polymer_atoms(&e.atoms, &e.residues, e.pdb_chain_id, buffer);
        }
        MoleculeEntity::SmallMolecule(e) => {
            for atom in &e.atoms {
                write_atom_row(atom, b' ', e.residue_name, 1, buffer);
            }
        }
        MoleculeEntity::Bulk(e) =>
        {
            #[allow(
                clippy::cast_possible_truncation,
                clippy::cast_possible_wrap,
                reason = "atom count fits in i32 for valid structures"
            )]
            for (i, atom) in e.atoms.iter().enumerate() {
                write_atom_row(
                    atom,
                    b' ',
                    e.residue_name,
                    (i as i32) + 1,
                    buffer,
                );
            }
        }
    }
}

fn write_polymer_atoms(
    atoms: &[Atom],
    residues: &[Residue],
    chain_id: u8,
    buffer: &mut Vec<u8>,
) {
    for residue in residues {
        for idx in residue.atom_range.clone() {
            write_atom_row(
                &atoms[idx],
                chain_id,
                residue.name,
                residue.number,
                buffer,
            );
        }
    }
}

fn write_atom_row(
    atom: &Atom,
    chain_id: u8,
    res_name: [u8; 3],
    res_num: i32,
    buffer: &mut Vec<u8>,
) {
    buffer.extend_from_slice(&atom.position.x.to_be_bytes());
    buffer.extend_from_slice(&atom.position.y.to_be_bytes());
    buffer.extend_from_slice(&atom.position.z.to_be_bytes());
    buffer.push(chain_id);
    buffer.extend_from_slice(&res_name);
    buffer.extend_from_slice(&res_num.to_be_bytes());
    buffer.extend_from_slice(&atom.name);
    let sym = atom.element.symbol();
    let sym_bytes = sym.as_bytes();
    buffer.push(sym_bytes.first().copied().unwrap_or(b'X'));
    buffer.push(sym_bytes.get(1).copied().unwrap_or(0));
}

//! Deserialization for the ASSEM01 binary wire format.

use std::collections::HashSet;

use glam::Vec3;

use super::{molecule_type_from_wire, ASSEMBLY_MAGIC};
use crate::element::Element;
use crate::entity::molecule::atom::Atom;
use crate::entity::molecule::bulk::BulkEntity;
use crate::entity::molecule::id::{EntityId, EntityIdAllocator};
use crate::entity::molecule::nucleic_acid::NAEntity;
use crate::entity::molecule::polymer::Residue;
use crate::entity::molecule::protein::ProteinEntity;
use crate::entity::molecule::small_molecule::SmallMoleculeEntity;
use crate::entity::molecule::{MoleculeEntity, MoleculeType};
use crate::ops::codec::AdapterError;

/// Width of one ASSEM01 atom row in bytes.
const ATOM_ROW_BYTES: usize = 26;

/// One atom row decoded from the wire, paired with the per-atom
/// residue/chain context used to group atoms into residues.
struct AtomRow {
    atom: Atom,
    chain_id: u8,
    res_name: [u8; 3],
    res_num: i32,
}

fn read_atom_row(cursor: &[u8]) -> Result<AtomRow, AdapterError> {
    let x = f32::from_be_bytes(cursor[0..4].try_into().map_err(|_| {
        AdapterError::SerializationError("Invalid x coordinate".to_owned())
    })?);
    let y = f32::from_be_bytes(cursor[4..8].try_into().map_err(|_| {
        AdapterError::SerializationError("Invalid y coordinate".to_owned())
    })?);
    let z = f32::from_be_bytes(cursor[8..12].try_into().map_err(|_| {
        AdapterError::SerializationError("Invalid z coordinate".to_owned())
    })?);

    let chain_id = cursor[12];

    let mut res_name = [0u8; 3];
    res_name.copy_from_slice(&cursor[13..16]);

    let res_num =
        i32::from_be_bytes(cursor[16..20].try_into().map_err(|_| {
            AdapterError::SerializationError(
                "Invalid residue number".to_owned(),
            )
        })?);

    let mut atom_name = [0u8; 4];
    atom_name.copy_from_slice(&cursor[20..24]);

    let sym_str = std::str::from_utf8(&cursor[24..26])
        .unwrap_or("")
        .trim_matches('\0')
        .trim();
    let element = Element::from_symbol(sym_str);

    Ok(AtomRow {
        atom: Atom {
            position: Vec3::new(x, y, z),
            occupancy: 1.0,
            b_factor: 0.0,
            element,
            name: atom_name,
            formal_charge: 0,
        },
        chain_id,
        res_name,
        res_num,
    })
}

/// Read `atom_count` atom rows from a cursor, returning the rows and
/// the remaining bytes.
fn read_atom_rows(
    mut cursor: &[u8],
    atom_count: usize,
) -> Result<(Vec<AtomRow>, &[u8]), AdapterError> {
    let mut rows = Vec::with_capacity(atom_count);
    for _ in 0..atom_count {
        rows.push(read_atom_row(cursor)?);
        cursor = &cursor[ATOM_ROW_BYTES..];
    }
    Ok((rows, cursor))
}

/// Split atom rows into `(Vec<Atom>, Vec<Residue>)` by grouping
/// consecutive rows sharing the same residue number.
fn into_atoms_and_residues(rows: Vec<AtomRow>) -> (Vec<Atom>, Vec<Residue>) {
    let mut atoms = Vec::with_capacity(rows.len());
    let mut residues = Vec::new();
    let mut current_res_num: Option<i32> = None;
    let mut current_res_name: [u8; 3] = [b' '; 3];
    let mut current_start = 0usize;

    for row in rows {
        let new_residue = current_res_num.is_none_or(|n| n != row.res_num);
        if new_residue {
            if let Some(rn) = current_res_num {
                residues.push(Residue {
                    name: current_res_name,
                    label_seq_id: rn,
                    auth_seq_id: None,
                    auth_comp_id: None,
                    ins_code: None,
                    atom_range: current_start..atoms.len(),
                });
            }
            current_start = atoms.len();
            current_res_num = Some(row.res_num);
            current_res_name = row.res_name;
        }
        atoms.push(row.atom);
    }
    if let Some(rn) = current_res_num {
        residues.push(Residue {
            name: current_res_name,
            label_seq_id: rn,
            auth_seq_id: None,
            auth_comp_id: None,
            ins_code: None,
            atom_range: current_start..atoms.len(),
        });
    }

    (atoms, residues)
}

fn build_entity(
    id: EntityId,
    mol_type: MoleculeType,
    rows: Vec<AtomRow>,
) -> MoleculeEntity {
    let pdb_chain_id = rows.first().map_or(b' ', |r| r.chain_id);

    match mol_type {
        MoleculeType::Protein => {
            let (atoms, residues) = into_atoms_and_residues(rows);
            MoleculeEntity::Protein(ProteinEntity::new(
                id,
                atoms,
                residues,
                pdb_chain_id,
                None,
            ))
        }
        MoleculeType::DNA | MoleculeType::RNA => {
            let (atoms, residues) = into_atoms_and_residues(rows);
            MoleculeEntity::NucleicAcid(NAEntity::new(
                id,
                mol_type,
                atoms,
                residues,
                pdb_chain_id,
                None,
            ))
        }
        MoleculeType::Water | MoleculeType::Solvent => {
            let residue_name = rows.first().map_or([b' '; 3], |r| r.res_name);
            let mut seen = HashSet::new();
            for row in &rows {
                let _ = seen.insert((row.chain_id, row.res_num));
            }
            let molecule_count = seen.len();
            let atoms = rows.into_iter().map(|r| r.atom).collect();
            MoleculeEntity::Bulk(BulkEntity::new(
                id,
                mol_type,
                atoms,
                residue_name,
                molecule_count,
            ))
        }
        MoleculeType::Ligand
        | MoleculeType::Ion
        | MoleculeType::Cofactor
        | MoleculeType::Lipid => {
            let residue_name = rows.first().map_or([b' '; 3], |r| r.res_name);
            let atoms = rows.into_iter().map(|r| r.atom).collect();
            MoleculeEntity::SmallMolecule(SmallMoleculeEntity::new(
                id,
                mol_type,
                atoms,
                residue_name,
            ))
        }
    }
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
    if bytes.len() < headers_end + total_atoms * ATOM_ROW_BYTES {
        return Err(AdapterError::InvalidFormat(
            "Data too short for atom data".to_owned(),
        ));
    }

    let mut cursor = &bytes[headers_end..];
    let mut entities = Vec::with_capacity(entity_count);
    let mut allocator = EntityIdAllocator::new();

    for (mol_type, atom_count) in entity_headers {
        let (rows, rest) = read_atom_rows(cursor, atom_count)?;
        cursor = rest;
        let id = allocator.allocate();
        entities.push(build_entity(id, mol_type, rows));
    }

    Ok(entities)
}

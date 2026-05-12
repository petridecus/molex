//! Coords ↔ Entity bridge: splitting, merging, and conversion.

use std::collections::BTreeMap;

use glam::Vec3;

use super::Coords;
use crate::element::Element;
use crate::entity::molecule::atom::Atom;
use crate::entity::molecule::bulk::BulkEntity;
use crate::entity::molecule::classify::classify_residue;
use crate::entity::molecule::id::{EntityId, EntityIdAllocator};
use crate::entity::molecule::nucleic_acid::NAEntity;
use crate::entity::molecule::protein::ProteinEntity;
use crate::entity::molecule::small_molecule::SmallMoleculeEntity;
use crate::entity::molecule::{MoleculeEntity, MoleculeType, Residue};

/// Extract a `Vec<Atom>` and `Vec<Residue>` from flat `Coords` indices.
///
/// Groups atoms by residue number, preserving order. Used by concrete
/// entity `from_coords_indices` constructors.
#[must_use]
pub fn extract_atom_set_and_residues(
    indices: &[usize],
    coords: &Coords,
) -> (Vec<Atom>, Vec<Residue>) {
    let mut atoms = Vec::with_capacity(indices.len());
    let mut residues = Vec::new();

    let mut current_res_num: Option<i32> = None;
    let mut current_start = 0usize;

    for &idx in indices {
        let res_num = coords.res_nums[idx];
        if current_res_num.is_some_and(|r| r != res_num) {
            let end = atoms.len();
            if let Some(rn) = current_res_num {
                residues.push(Residue {
                    name: coords.res_names[indices[current_start]],
                    label_seq_id: rn,
                    auth_seq_id: None,
                    auth_comp_id: None,
                    ins_code: None,
                    atom_range: current_start..end,
                });
            }
            current_start = end;
        }
        current_res_num = Some(res_num);
        let ca = &coords.atoms[idx];
        atoms.push(Atom {
            position: Vec3::new(ca.x, ca.y, ca.z),
            occupancy: ca.occupancy,
            b_factor: ca.b_factor,
            element: coords
                .elements
                .get(idx)
                .copied()
                .unwrap_or(Element::Unknown),
            name: coords.atom_names[idx],
            formal_charge: 0,
        });
    }
    if let Some(rn) = current_res_num {
        residues.push(Residue {
            name: coords.res_names[indices[current_start]],
            label_seq_id: rn,
            auth_seq_id: None,
            auth_comp_id: None,
            ins_code: None,
            atom_range: current_start..atoms.len(),
        });
    }

    (atoms, residues)
}

/// Build a `MoleculeEntity` from a `Coords` + known `MoleculeType`.
///
/// Used by ASSEM01 deserialization where each entity's type is stored in the
/// header.
#[must_use]
pub fn coords_to_molecule_entity(
    id: EntityId,
    mol_type: MoleculeType,
    coords: &Coords,
) -> MoleculeEntity {
    let indices: Vec<usize> = (0..coords.num_atoms).collect();
    let pdb_id = coords.chain_ids.first().copied().unwrap_or(b' ');
    match mol_type {
        MoleculeType::Protein => MoleculeEntity::Protein(
            ProteinEntity::from_coords_indices(id, &indices, coords, pdb_id),
        ),
        MoleculeType::DNA | MoleculeType::RNA => {
            MoleculeEntity::NucleicAcid(NAEntity::from_coords_indices(
                id, mol_type, &indices, coords, pdb_id,
            ))
        }
        MoleculeType::Water | MoleculeType::Solvent => MoleculeEntity::Bulk(
            BulkEntity::from_coords_indices(id, mol_type, &indices, coords),
        ),
        _ => MoleculeEntity::SmallMolecule(
            SmallMoleculeEntity::from_coords_indices(
                id, mol_type, &indices, coords,
            ),
        ),
    }
}

// ---------------------------------------------------------------------------
// Entity splitting / merging
// ---------------------------------------------------------------------------

/// Key for grouping atoms into entities.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum EntityKey {
    Chain(u8, MoleculeTypeOrd),
    Water,
    Solvent,
    SmallMolecule(u8, i32, MoleculeTypeOrd),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct MoleculeTypeOrd(MoleculeType);

impl PartialOrd for MoleculeTypeOrd {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for MoleculeTypeOrd {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        (self.0 as u8).cmp(&(other.0 as u8))
    }
}

/// Check whether a set of atom indices contains backbone atoms N, CA, C.
fn residue_has_backbone(indices: &[usize], coords: &Coords) -> bool {
    let mut has_n = false;
    let mut has_ca = false;
    let mut has_c = false;
    for &idx in indices {
        let name = &coords.atom_names[idx];
        match name {
            [b' ', b'N', b' ', b' '] | [b'N', b' ', b' ', b' '] => has_n = true,
            [b' ', b'C', b'A', b' '] | [b'C', b'A', b' ', b' '] => {
                has_ca = true;
            }
            [b' ', b'C', b' ', b' '] | [b'C', b' ', b' ', b' '] => has_c = true,
            _ => {}
        }
    }
    has_n && has_ca && has_c
}

/// Split a flat `Coords` into per-entity `MoleculeEntity` groups.
#[must_use]
#[allow(
    clippy::excessive_nesting,
    reason = "grouping logic with nested match arms is natural"
)]
#[allow(
    clippy::too_many_lines,
    reason = "entity splitting with modified-residue merging is a single \
              logical operation"
)]
pub fn split_into_entities(coords: &Coords) -> Vec<MoleculeEntity> {
    let mut groups: BTreeMap<EntityKey, Vec<usize>> = BTreeMap::new();

    for i in 0..coords.num_atoms {
        let res_name = std::str::from_utf8(&coords.res_names[i])
            .unwrap_or("")
            .trim();
        let mol_type = classify_residue(res_name);
        let chain_id = coords.chain_ids[i];

        let key = match mol_type {
            MoleculeType::Water => EntityKey::Water,
            MoleculeType::Solvent => EntityKey::Solvent,
            MoleculeType::Protein | MoleculeType::DNA | MoleculeType::RNA => {
                EntityKey::Chain(chain_id, MoleculeTypeOrd(mol_type))
            }
            MoleculeType::Ligand
            | MoleculeType::Ion
            | MoleculeType::Cofactor
            | MoleculeType::Lipid => {
                let res_num = coords.res_nums[i];
                EntityKey::SmallMolecule(
                    chain_id,
                    res_num,
                    MoleculeTypeOrd(mol_type),
                )
            }
        };

        groups.entry(key).or_default().push(i);
    }

    // Merge modified amino acids back into their protein chain.
    let protein_chains: Vec<u8> = groups
        .keys()
        .filter_map(|k| match k {
            EntityKey::Chain(cid, mt) if mt.0 == MoleculeType::Protein => {
                Some(*cid)
            }
            _ => None,
        })
        .collect();

    let merge_keys: Vec<EntityKey> = groups
        .iter()
        .filter_map(|(key, indices)| {
            if let EntityKey::SmallMolecule(chain_id, _, _) = key {
                if protein_chains.contains(chain_id)
                    && residue_has_backbone(indices, coords)
                {
                    return Some(key.clone());
                }
            }
            None
        })
        .collect();

    for key in merge_keys {
        if let EntityKey::SmallMolecule(chain_id, _, _) = &key {
            let chain_key = EntityKey::Chain(
                *chain_id,
                MoleculeTypeOrd(MoleculeType::Protein),
            );
            if let Some(indices) = groups.remove(&key) {
                groups.entry(chain_key).or_default().extend(indices);
            }
        }
    }

    let mut allocator = EntityIdAllocator::new();

    groups
        .into_iter()
        .map(|(key, indices)| {
            let (mol_type, pdb_chain_id) = match &key {
                EntityKey::Chain(cid, mt)
                | EntityKey::SmallMolecule(cid, _, mt) => (mt.0, *cid),
                EntityKey::Water => (MoleculeType::Water, b' '),
                EntityKey::Solvent => (MoleculeType::Solvent, b' '),
            };

            let id = allocator.allocate();

            match mol_type {
                MoleculeType::Protein => {
                    MoleculeEntity::Protein(ProteinEntity::from_coords_indices(
                        id,
                        &indices,
                        coords,
                        pdb_chain_id,
                    ))
                }
                MoleculeType::DNA | MoleculeType::RNA => {
                    MoleculeEntity::NucleicAcid(NAEntity::from_coords_indices(
                        id,
                        mol_type,
                        &indices,
                        coords,
                        pdb_chain_id,
                    ))
                }
                MoleculeType::Water | MoleculeType::Solvent => {
                    MoleculeEntity::Bulk(BulkEntity::from_coords_indices(
                        id, mol_type, &indices, coords,
                    ))
                }
                _ => MoleculeEntity::SmallMolecule(
                    SmallMoleculeEntity::from_coords_indices(
                        id, mol_type, &indices, coords,
                    ),
                ),
            }
        })
        .collect()
}

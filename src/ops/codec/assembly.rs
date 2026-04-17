//! Assembly helpers: convenience functions for working with
//! `Vec<MoleculeEntity>`.

use glam::Vec3;

use super::serialize::serialize_entities;
use super::{merge_entities, split_into_entities, Coords, CoordsError};
use crate::entity::molecule::id::EntityIdAllocator;
use crate::entity::molecule::{MoleculeEntity, MoleculeType};

/// Protein-only Coords derived from entities (for serialization).
#[must_use]
pub fn protein_coords(entities: &[MoleculeEntity]) -> Coords {
    let proteins: Vec<MoleculeEntity> = entities
        .iter()
        .filter(|e| e.molecule_type() == MoleculeType::Protein)
        .cloned()
        .collect();
    merge_entities(&proteins)
}

/// ASSEM01-encode a raw entity slice (includes molecule type metadata).
///
/// # Errors
///
/// Returns `CoordsError` if serialization fails.
pub fn assembly_bytes(
    entities: &[MoleculeEntity],
) -> Result<Vec<u8>, CoordsError> {
    serialize_entities(entities)
}

/// CA positions from protein entities.
#[must_use]
pub fn ca_positions(entities: &[MoleculeEntity]) -> Vec<Vec3> {
    entities
        .iter()
        .filter_map(MoleculeEntity::as_protein)
        .flat_map(|p| p.to_backbone().into_iter().map(|bb| bb.ca))
        .collect()
}

/// Number of protein residues.
#[must_use]
pub fn residue_count(entities: &[MoleculeEntity]) -> usize {
    entities
        .iter()
        .filter_map(MoleculeEntity::as_protein)
        .map(|p| p.residues.len())
        .sum()
}

/// Replace protein entity coords (keeps non-protein entities).
/// Splits the incoming combined protein coords by chain ID so each
/// entity only receives its own chain's atoms, avoiding duplication.
pub fn update_protein_entities(
    entities: &mut Vec<MoleculeEntity>,
    protein: &Coords,
) {
    let new_protein = split_into_entities(protein);

    entities.retain(|e| e.molecule_type() != MoleculeType::Protein);

    let mut updated: Vec<MoleculeEntity> = new_protein
        .into_iter()
        .filter(|e| e.molecule_type() == MoleculeType::Protein)
        .collect();
    updated.append(entities);

    let mut allocator = EntityIdAllocator::new();
    for entity in &mut updated {
        entity.set_id(allocator.allocate());
    }

    *entities = updated;
}

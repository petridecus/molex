//! Entity-utility helpers (CA positions, residue count) and the
//! `update_protein_entities` live-position-update wire (Coords-shaped).

use glam::Vec3;

use super::{split_into_entities, Coords};
use crate::entity::molecule::id::EntityIdAllocator;
use crate::entity::molecule::{MoleculeEntity, MoleculeType};

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
#[allow(
    dead_code,
    reason = "test-only caller in molex; production callers live in viso"
)]
pub(crate) fn update_protein_entities(
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

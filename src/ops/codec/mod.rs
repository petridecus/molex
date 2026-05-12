//! Adapter / wire error type and protein-CA helper.

mod types;

use glam::Vec3;
pub use types::AdapterError;

use crate::entity::molecule::MoleculeEntity;

/// CA positions across every protein entity, residue-ordered per chain.
#[must_use]
pub fn ca_positions(entities: &[MoleculeEntity]) -> Vec<Vec3> {
    entities
        .iter()
        .filter_map(MoleculeEntity::as_protein)
        .flat_map(|p| p.to_backbone().into_iter().map(|bb| bb.ca))
        .collect()
}

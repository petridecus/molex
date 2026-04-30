//! Bulk entity — many identical small molecules (water, solvent).

use std::collections::HashSet;

use glam::Vec3;

use super::atom::Atom;
use super::id::EntityId;
use super::traits::Entity;
use super::MoleculeType;
use crate::element::Element;
use crate::ops::codec::Coords;

/// A group of identical small molecules (water, solvent).
#[derive(Debug, Clone)]
pub struct BulkEntity {
    /// Unique entity identifier.
    pub id: EntityId,
    /// Molecule type (Water or Solvent).
    pub mol_type: MoleculeType,
    /// Atom data for all molecules in this group.
    pub atoms: Vec<Atom>,
    /// 3-character residue code (e.g. b"HOH").
    pub residue_name: [u8; 3],
    /// Number of individual molecules in this group.
    pub molecule_count: usize,
}

impl BulkEntity {
    /// Construct from a list of atoms. `molecule_count` is supplied by
    /// the caller (typically the number of residues in the source).
    #[must_use]
    pub fn new(
        id: EntityId,
        mol_type: MoleculeType,
        atoms: Vec<Atom>,
        residue_name: [u8; 3],
        molecule_count: usize,
    ) -> Self {
        Self {
            id,
            mol_type,
            atoms,
            residue_name,
            molecule_count,
        }
    }

    /// Construct from flat `Coords` atom indices during entity splitting.
    #[must_use]
    pub fn from_coords_indices(
        id: EntityId,
        mol_type: MoleculeType,
        indices: &[usize],
        coords: &Coords,
    ) -> Self {
        let mut atoms = Vec::with_capacity(indices.len());
        let mut seen = HashSet::new();
        for &idx in indices {
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
            });
            let _ = seen.insert((coords.chain_ids[idx], coords.res_nums[idx]));
        }
        let residue_name = coords.res_names[indices[0]];
        Self {
            id,
            mol_type,
            atoms,
            residue_name,
            molecule_count: seen.len(),
        }
    }
}

impl Entity for BulkEntity {
    fn id(&self) -> EntityId {
        self.id
    }
    fn molecule_type(&self) -> MoleculeType {
        self.mol_type
    }
    fn atoms(&self) -> &[Atom] {
        &self.atoms
    }
}

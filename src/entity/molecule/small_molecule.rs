//! Small molecule entity — a single non-polymer molecule (ligand, ion,
//! cofactor, lipid).

use glam::Vec3;

use super::atom::Atom;
use super::id::EntityId;
use super::traits::Entity;
use super::MoleculeType;
use crate::element::Element;
use crate::ops::codec::Coords;

/// A single non-polymer molecule.
#[derive(Debug, Clone)]
pub struct SmallMoleculeEntity {
    /// Unique entity identifier.
    pub id: EntityId,
    /// Molecule type (Ligand, Ion, Cofactor, or Lipid).
    pub mol_type: MoleculeType,
    /// Atom data.
    pub atoms: Vec<Atom>,
    /// 3-character residue code (e.g. b"ATP").
    pub residue_name: [u8; 3],
    /// Human-readable name for display (e.g. "Chlorophyll A").
    pub display_name: String,
}

impl SmallMoleculeEntity {
    /// Construct from flat `Coords` atom indices during entity splitting.
    #[must_use]
    pub fn from_coords_indices(
        id: EntityId,
        mol_type: MoleculeType,
        indices: &[usize],
        coords: &Coords,
    ) -> Self {
        let atoms: Vec<Atom> = indices
            .iter()
            .map(|&idx| {
                let ca = &coords.atoms[idx];
                Atom {
                    position: Vec3::new(ca.x, ca.y, ca.z),
                    occupancy: ca.occupancy,
                    b_factor: ca.b_factor,
                    element: coords
                        .elements
                        .get(idx)
                        .copied()
                        .unwrap_or(Element::Unknown),
                    name: coords.atom_names[idx],
                }
            })
            .collect();
        let residue_name = coords.res_names[indices[0]];
        let rn_str = std::str::from_utf8(&residue_name).unwrap_or("???").trim();
        let display_name =
            super::classify::small_molecule_display_name(mol_type, rn_str);
        Self {
            id,
            mol_type,
            atoms,
            residue_name,
            display_name,
        }
    }
}

impl Entity for SmallMoleculeEntity {
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

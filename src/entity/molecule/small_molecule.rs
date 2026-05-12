//! Small molecule entity — a single non-polymer molecule (ligand, ion,
//! cofactor, lipid).

use glam::Vec3;

use super::atom::Atom;
use super::id::EntityId;
use super::traits::Entity;
use super::MoleculeType;
use crate::analysis::{infer_bonds, DEFAULT_TOLERANCE};
use crate::atom_id::AtomId;
use crate::bond::CovalentBond;
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
    /// Intra-entity covalent bonds with `AtomId` endpoints. Populated
    /// at construction by distance-based [`infer_bonds`] over the
    /// entity's atoms.
    pub bonds: Vec<CovalentBond>,
}

impl SmallMoleculeEntity {
    /// Construct from a list of atoms. Bonds are inferred via
    /// distance-based [`infer_bonds`]. The display name is derived from
    /// `residue_name` and `mol_type`.
    #[must_use]
    #[allow(
        clippy::cast_possible_truncation,
        reason = "atom indices bounded by small-molecule atom count"
    )]
    pub fn new(
        id: EntityId,
        mol_type: MoleculeType,
        atoms: Vec<Atom>,
        residue_name: [u8; 3],
    ) -> Self {
        let rn_str = std::str::from_utf8(&residue_name).unwrap_or("???").trim();
        let display_name =
            super::classify::small_molecule_display_name(mol_type, rn_str);
        let bonds = infer_bonds(&atoms, DEFAULT_TOLERANCE)
            .into_iter()
            .map(|b| CovalentBond {
                a: AtomId {
                    entity: id,
                    index: b.atom_a as u32,
                },
                b: AtomId {
                    entity: id,
                    index: b.atom_b as u32,
                },
                order: b.order,
            })
            .collect();
        Self {
            id,
            mol_type,
            atoms,
            residue_name,
            display_name,
            bonds,
        }
    }

    /// Construct from flat `Coords` atom indices during entity splitting.
    #[must_use]
    #[allow(
        clippy::cast_possible_truncation,
        reason = "atom indices bounded by small-molecule atom count"
    )]
    pub(crate) fn from_coords_indices(
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
        let bonds = infer_bonds(&atoms, DEFAULT_TOLERANCE)
            .into_iter()
            .map(|b| CovalentBond {
                a: AtomId {
                    entity: id,
                    index: b.atom_a as u32,
                },
                b: AtomId {
                    entity: id,
                    index: b.atom_b as u32,
                },
                order: b.order,
            })
            .collect();
        Self {
            id,
            mol_type,
            atoms,
            residue_name,
            display_name,
            bonds,
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

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::float_cmp)]
mod tests {
    use super::*;
    use crate::ops::codec::{split_into_entities, Coords, CoordsAtom};

    fn make_atom(x: f32, y: f32, z: f32) -> CoordsAtom {
        CoordsAtom {
            x,
            y,
            z,
            occupancy: 1.0,
            b_factor: 0.0,
        }
    }

    fn res_name(s: &str) -> [u8; 3] {
        let mut n = [b' '; 3];
        for (i, b) in s.bytes().take(3).enumerate() {
            n[i] = b;
        }
        n
    }

    fn atom_name(s: &str) -> [u8; 4] {
        let mut n = [b' '; 4];
        for (i, b) in s.bytes().take(4).enumerate() {
            n[i] = b;
        }
        n
    }

    /// A three-atom water-like ligand: O with two Hs at bonding distance.
    #[test]
    fn small_molecule_populates_bonds_from_positions() {
        let coords = Coords {
            num_atoms: 3,
            atoms: vec![
                make_atom(0.0, 0.0, 0.0),
                make_atom(0.95, 0.0, 0.0),
                make_atom(-0.24, 0.92, 0.0),
            ],
            chain_ids: vec![b'L', b'L', b'L'],
            res_names: vec![res_name("HOL"); 3],
            res_nums: vec![1, 1, 1],
            atom_names: vec![atom_name("O1"), atom_name("H1"), atom_name("H2")],
            elements: vec![Element::O, Element::H, Element::H],
        };
        let entities = split_into_entities(&coords);
        assert_eq!(entities.len(), 1);
        let sm = entities[0].as_small_molecule().unwrap();
        assert_eq!(sm.bonds.len(), 2);
        for b in &sm.bonds {
            assert_eq!(b.a.entity, sm.id);
            assert_eq!(b.b.entity, sm.id);
        }
    }
}

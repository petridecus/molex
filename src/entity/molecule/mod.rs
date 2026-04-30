//! Entity types and molecule classification.

/// Core atom type.
pub mod atom;
/// Bulk entity (water, solvent).
pub mod bulk;
pub(crate) mod classify;
/// Opaque entity ID with controlled allocation.
pub mod id;
/// Nucleic acid entity (DNA, RNA).
pub mod nucleic_acid;
/// Polymer residue (shared by protein and nucleic acid entities).
pub mod polymer;
/// Protein entity with residues and segment breaks.
pub mod protein;
/// Small molecule entity (ligand, ion, cofactor).
pub mod small_molecule;
/// Entity and Polymer traits.
pub mod traits;

pub use atom::Atom;
pub use classify::classify_residue;
use glam::Vec3;
pub use id::{EntityId, EntityIdAllocator};
pub use nucleic_acid::NucleotideRing;
pub use polymer::Residue;

use self::bulk::BulkEntity;
use self::nucleic_acid::NAEntity;
use self::protein::ProteinEntity;
use self::small_molecule::SmallMoleculeEntity;
use self::traits::Entity;
use crate::analysis::aabb::Aabb;
// Coords↔Entity bridge functions live in `ops::codec`.
pub use crate::ops::codec::{
    coords_to_molecule_entity, extract_atom_set_and_residues, extract_by_type,
    merge_entities, split_into_entities,
};
use crate::ops::codec::{Coords, CoordsAtom};
/// Classification of molecule types found in structural biology files.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MoleculeType {
    /// Amino acid polymer.
    Protein,
    /// Deoxyribonucleic acid polymer.
    DNA,
    /// Ribonucleic acid polymer.
    RNA,
    /// Non-polymer small molecule (drug, substrate, etc.).
    Ligand,
    /// Single-atom metal or halide ion.
    Ion,
    /// Water molecule.
    Water,
    /// Lipid or detergent molecule.
    Lipid,
    /// Enzyme cofactor (heme, NAD, FAD, Fe-S cluster, etc.).
    Cofactor,
    /// Crystallization solvent or buffer artifact.
    Solvent,
}

// ---------------------------------------------------------------------------
// MoleculeEntity enum
// ---------------------------------------------------------------------------

/// A single entity: one logical molecule (a protein chain, a ligand, waters,
/// etc.).
///
/// This is an enum wrapping concrete entity types. Each variant owns its
/// entity data directly. Use the accessor methods to work with entities
/// polymorphically.
#[derive(Debug, Clone)]
pub enum MoleculeEntity {
    /// A single protein chain.
    Protein(ProteinEntity),
    /// A single DNA or RNA chain.
    NucleicAcid(NAEntity),
    /// A single non-polymer molecule (ligand, ion, cofactor, lipid).
    SmallMolecule(SmallMoleculeEntity),
    /// A group of identical small molecules (water, solvent).
    Bulk(BulkEntity),
}

impl MoleculeEntity {
    // -- Entity trait delegation --

    /// Unique entity identifier.
    #[must_use]
    pub fn id(&self) -> EntityId {
        match self {
            MoleculeEntity::Protein(e) => e.id(),
            MoleculeEntity::NucleicAcid(e) => e.id(),
            MoleculeEntity::SmallMolecule(e) => e.id(),
            MoleculeEntity::Bulk(e) => e.id(),
        }
    }

    /// Classification of this entity's molecule type.
    #[must_use]
    pub fn molecule_type(&self) -> MoleculeType {
        match self {
            MoleculeEntity::Protein(e) => e.molecule_type(),
            MoleculeEntity::NucleicAcid(e) => e.molecule_type(),
            MoleculeEntity::SmallMolecule(e) => e.molecule_type(),
            MoleculeEntity::Bulk(e) => e.molecule_type(),
        }
    }

    /// Reference to the underlying `Vec<Atom>`.
    #[must_use]
    pub fn atom_set(&self) -> &[Atom] {
        match self {
            MoleculeEntity::Protein(e) => e.atoms(),
            MoleculeEntity::NucleicAcid(e) => e.atoms(),
            MoleculeEntity::SmallMolecule(e) => e.atoms(),
            MoleculeEntity::Bulk(e) => e.atoms(),
        }
    }

    /// Mutable reference to the underlying `Vec<Atom>`.
    pub fn atom_set_mut(&mut self) -> &mut [Atom] {
        match self {
            MoleculeEntity::Protein(e) => &mut e.atoms,
            MoleculeEntity::NucleicAcid(e) => &mut e.atoms,
            MoleculeEntity::SmallMolecule(e) => &mut e.atoms,
            MoleculeEntity::Bulk(e) => &mut e.atoms,
        }
    }

    /// All atom positions as Vec3.
    #[must_use]
    pub fn positions(&self) -> Vec<Vec3> {
        self.atom_set().iter().map(|a| a.position).collect()
    }

    /// Number of atoms in this entity.
    #[must_use]
    pub fn atom_count(&self) -> usize {
        self.atom_set().len()
    }

    // -- Variant-specific accessors --

    /// If this entity is a protein, return it.
    #[must_use]
    pub fn as_protein(&self) -> Option<&ProteinEntity> {
        match self {
            MoleculeEntity::Protein(e) => Some(e),
            _ => None,
        }
    }

    /// If this entity is a nucleic acid, return it.
    #[must_use]
    pub fn as_nucleic_acid(&self) -> Option<&NAEntity> {
        match self {
            MoleculeEntity::NucleicAcid(e) => Some(e),
            _ => None,
        }
    }

    /// PDB chain identifier byte for polymer entities, `None` for others.
    #[must_use]
    pub fn pdb_chain_id(&self) -> Option<u8> {
        match self {
            MoleculeEntity::Protein(e) => Some(e.pdb_chain_id),
            MoleculeEntity::NucleicAcid(e) => Some(e.pdb_chain_id),
            _ => None,
        }
    }

    /// If this entity is a small molecule, return it.
    #[must_use]
    pub fn as_small_molecule(&self) -> Option<&SmallMoleculeEntity> {
        match self {
            MoleculeEntity::SmallMolecule(e) => Some(e),
            _ => None,
        }
    }

    /// If this entity is a bulk entity, return it.
    #[must_use]
    pub fn as_bulk(&self) -> Option<&BulkEntity> {
        match self {
            MoleculeEntity::Bulk(e) => Some(e),
            _ => None,
        }
    }

    /// Convert to a flat `Coords` for serialization or interop.
    #[must_use]
    pub fn to_coords(&self) -> Coords {
        match self {
            MoleculeEntity::Protein(e) => {
                polymer_entity_to_coords(&e.atoms, &e.residues, e.pdb_chain_id)
            }
            MoleculeEntity::NucleicAcid(e) => {
                polymer_entity_to_coords(&e.atoms, &e.residues, e.pdb_chain_id)
            }
            MoleculeEntity::SmallMolecule(e) => atoms_to_coords(
                &e.atoms,
                e.residue_name,
                vec![1; e.atoms.len()],
            ),
            MoleculeEntity::Bulk(e) => {
                let n = e.atoms.len();
                #[allow(
                    clippy::cast_possible_truncation,
                    clippy::cast_possible_wrap,
                    reason = "atom count fits in i32 for valid structures"
                )]
                let res_nums = (1..=n as i32).collect();
                atoms_to_coords(&e.atoms, e.residue_name, res_nums)
            }
        }
    }

    /// Compute the axis-aligned bounding box for this entity's atoms.
    #[must_use]
    pub fn aabb(&self) -> Option<Aabb> {
        Aabb::from_positions(&self.positions())
    }

    /// Human-readable label (e.g. "Protein Chain A", "Ligand (ATP)", "Zn2+
    /// Ion").
    #[must_use]
    #[allow(
        clippy::too_many_lines,
        reason = "match arms per molecule type are straightforward"
    )]
    pub fn label(&self) -> String {
        let mol_type = self.molecule_type();
        match mol_type {
            MoleculeType::Protein => polymer_label(self, "Protein"),
            MoleculeType::DNA => polymer_label(self, "DNA"),
            MoleculeType::RNA => polymer_label(self, "RNA"),
            MoleculeType::Ligand => {
                if let MoleculeEntity::SmallMolecule(e) = self {
                    format!("Ligand ({})", e.display_name)
                } else {
                    "Ligand".to_owned()
                }
            }
            MoleculeType::Ion => {
                if let MoleculeEntity::SmallMolecule(e) = self {
                    format!("{} Ion", e.display_name)
                } else {
                    "Ion".to_owned()
                }
            }
            MoleculeType::Water => {
                if let MoleculeEntity::Bulk(e) = self {
                    format!("Water ({} molecules)", e.molecule_count)
                } else {
                    "Water".to_owned()
                }
            }
            MoleculeType::Lipid => {
                if let MoleculeEntity::SmallMolecule(e) = self {
                    format!("Lipid ({})", e.display_name)
                } else {
                    format!("Lipid ({} molecules)", self.residue_count())
                }
            }
            MoleculeType::Cofactor => {
                if let MoleculeEntity::SmallMolecule(e) = self {
                    e.display_name.clone()
                } else {
                    "Cofactor".to_owned()
                }
            }
            MoleculeType::Solvent => {
                if let MoleculeEntity::Bulk(e) = self {
                    format!("Solvent ({} molecules)", e.molecule_count)
                } else {
                    "Solvent".to_owned()
                }
            }
        }
    }

    /// Whether this entity type participates in tab-cycling focus.
    /// Protein: no (focused at group level). Water, Ion: no (ambient).
    /// Ligand, DNA, RNA: yes.
    #[must_use]
    pub fn is_focusable(&self) -> bool {
        !matches!(
            self.molecule_type(),
            MoleculeType::Water | MoleculeType::Ion | MoleculeType::Solvent
        )
    }

    /// Number of residues (for polymer/nucleic) or molecules (for small
    /// mol/ion/water).
    #[must_use]
    pub fn residue_count(&self) -> usize {
        match self {
            MoleculeEntity::Protein(e) => e.residues.len(),
            MoleculeEntity::NucleicAcid(e) => e.residues.len(),
            MoleculeEntity::SmallMolecule(_) => 1,
            MoleculeEntity::Bulk(e) => e.molecule_count,
        }
    }

    /// Set the entity ID. Used during re-assignment (e.g. after
    /// `update_protein_entities`).
    pub fn set_id(&mut self, new_id: EntityId) {
        match self {
            MoleculeEntity::Protein(e) => e.id = new_id,
            MoleculeEntity::NucleicAcid(e) => e.id = new_id,
            MoleculeEntity::SmallMolecule(e) => e.id = new_id,
            MoleculeEntity::Bulk(e) => e.id = new_id,
        }
    }
}

/// Format a polymer entity label from its PDB chain ID.
fn polymer_label(entity: &MoleculeEntity, type_name: &str) -> String {
    entity.pdb_chain_id().map_or_else(
        || type_name.to_owned(),
        |id| format!("{type_name} {}", id as char),
    )
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Convert polymer entity fields to flat `Coords`.
fn polymer_entity_to_coords(
    atoms: &[Atom],
    residues: &[Residue],
    pdb_chain_id: u8,
) -> Coords {
    let n = atoms.len();
    let mut res_names = Vec::with_capacity(n);
    let mut res_nums = Vec::with_capacity(n);
    for residue in residues {
        for _ in residue.atom_range.clone() {
            res_names.push(residue.name);
            res_nums.push(residue.number);
        }
    }
    Coords {
        num_atoms: n,
        atoms: atoms
            .iter()
            .map(|a| CoordsAtom {
                x: a.position.x,
                y: a.position.y,
                z: a.position.z,
                occupancy: a.occupancy,
                b_factor: a.b_factor,
            })
            .collect(),
        chain_ids: vec![pdb_chain_id; n],
        res_names,
        res_nums,
        atom_names: atoms.iter().map(|a| a.name).collect(),
        elements: atoms.iter().map(|a| a.element).collect(),
    }
}

/// Convert a flat `Vec<Atom>` to `Coords` with uniform residue metadata.
fn atoms_to_coords(
    atoms: &[Atom],
    residue_name: [u8; 3],
    res_nums: Vec<i32>,
) -> Coords {
    let n = atoms.len();
    Coords {
        num_atoms: n,
        atoms: atoms
            .iter()
            .map(|a| CoordsAtom {
                x: a.position.x,
                y: a.position.y,
                z: a.position.z,
                occupancy: a.occupancy,
                b_factor: a.b_factor,
            })
            .collect(),
        chain_ids: vec![b' '; n],
        res_names: vec![residue_name; n],
        res_nums,
        atom_names: atoms.iter().map(|a| a.name).collect(),
        elements: atoms.iter().map(|a| a.element).collect(),
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::float_cmp)]
mod tests {
    use super::*;
    use crate::element::Element;
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

    /// Build a 2-residue protein (ALA-GLY) on chain A with backbone atoms
    /// at known positions.
    fn make_two_residue_protein() -> Coords {
        Coords {
            num_atoms: 8,
            atoms: vec![
                make_atom(1.0, 2.0, 3.0),    // res1 N
                make_atom(4.0, 5.0, 6.0),    // res1 CA
                make_atom(7.0, 8.0, 9.0),    // res1 C
                make_atom(10.0, 11.0, 12.0), // res1 O
                make_atom(13.0, 14.0, 15.0), // res2 N
                make_atom(16.0, 17.0, 18.0), // res2 CA
                make_atom(19.0, 20.0, 21.0), // res2 C
                make_atom(22.0, 23.0, 24.0), // res2 O
            ],
            chain_ids: vec![b'A'; 8],
            res_names: vec![
                res_name("ALA"),
                res_name("ALA"),
                res_name("ALA"),
                res_name("ALA"),
                res_name("GLY"),
                res_name("GLY"),
                res_name("GLY"),
                res_name("GLY"),
            ],
            res_nums: vec![1, 1, 1, 1, 2, 2, 2, 2],
            atom_names: vec![
                atom_name("N"),
                atom_name("CA"),
                atom_name("C"),
                atom_name("O"),
                atom_name("N"),
                atom_name("CA"),
                atom_name("C"),
                atom_name("O"),
            ],
            elements: vec![
                Element::N,
                Element::C,
                Element::C,
                Element::O,
                Element::N,
                Element::C,
                Element::C,
                Element::O,
            ],
        }
    }

    #[test]
    fn split_produces_protein_entity() {
        let coords = make_two_residue_protein();
        let entities = split_into_entities(&coords);
        assert_eq!(entities.len(), 1);
        assert_eq!(entities[0].molecule_type(), MoleculeType::Protein);
    }

    #[test]
    fn entity_id_is_set() {
        let coords = make_two_residue_protein();
        let entities = split_into_entities(&coords);
        // ID should be valid (non-zero value from allocator)
        let _id = entities[0].id();
    }

    #[test]
    fn atom_set_and_atom_count() {
        let coords = make_two_residue_protein();
        let entities = split_into_entities(&coords);
        assert_eq!(entities[0].atom_count(), 8);
        assert_eq!(entities[0].atom_set().len(), 8);
    }

    #[test]
    fn as_protein_returns_some_for_protein() {
        let coords = make_two_residue_protein();
        let entities = split_into_entities(&coords);
        assert!(entities[0].as_protein().is_some());
        assert!(entities[0].as_nucleic_acid().is_none());
        assert!(entities[0].as_small_molecule().is_none());
        assert!(entities[0].as_bulk().is_none());
    }

    #[test]
    fn label_for_protein() {
        let coords = make_two_residue_protein();
        let entities = split_into_entities(&coords);
        let label = entities[0].label();
        assert!(label.contains("Protein"), "label={label}");
        assert!(label.contains('A'), "label should contain chain A: {label}");
    }

    #[test]
    fn residue_count_for_protein() {
        let coords = make_two_residue_protein();
        let entities = split_into_entities(&coords);
        assert_eq!(entities[0].residue_count(), 2);
    }

    #[test]
    fn to_coords_roundtrip_preserves_atom_count() {
        let coords = make_two_residue_protein();
        let entities = split_into_entities(&coords);
        let recovered = entities[0].to_coords();
        assert_eq!(recovered.num_atoms, 8);
    }

    #[test]
    fn to_coords_roundtrip_preserves_positions() {
        let coords = make_two_residue_protein();
        let entities = split_into_entities(&coords);
        let recovered = entities[0].to_coords();
        assert!((recovered.atoms[0].x - 1.0).abs() < 1e-6);
        assert!((recovered.atoms[0].y - 2.0).abs() < 1e-6);
        assert!((recovered.atoms[0].z - 3.0).abs() < 1e-6);
        assert!((recovered.atoms[7].x - 22.0).abs() < 1e-6);
    }

    #[test]
    fn to_coords_roundtrip_preserves_chain_ids() {
        let coords = make_two_residue_protein();
        let entities = split_into_entities(&coords);
        let recovered = entities[0].to_coords();
        for &cid in &recovered.chain_ids {
            assert_eq!(cid, b'A');
        }
    }

    #[test]
    fn to_coords_roundtrip_preserves_res_names() {
        let coords = make_two_residue_protein();
        let entities = split_into_entities(&coords);
        let recovered = entities[0].to_coords();
        assert_eq!(recovered.res_names[0], res_name("ALA"));
        assert_eq!(recovered.res_names[4], res_name("GLY"));
    }

    #[test]
    fn aabb_is_some_for_nonempty_entity() {
        let coords = make_two_residue_protein();
        let entities = split_into_entities(&coords);
        let aabb = entities[0].aabb();
        assert!(aabb.is_some());
        let bb = aabb.unwrap();
        assert!(bb.min.x <= 1.0);
        assert!(bb.max.x >= 22.0);
    }

    #[test]
    fn positions_returns_all_atom_positions() {
        let coords = make_two_residue_protein();
        let entities = split_into_entities(&coords);
        let positions = entities[0].positions();
        assert_eq!(positions.len(), 8);
        assert!((positions[0].x - 1.0).abs() < 1e-6);
    }

    #[test]
    fn water_entity_accessors() {
        let coords = Coords {
            num_atoms: 2,
            atoms: vec![make_atom(1.0, 2.0, 3.0), make_atom(4.0, 5.0, 6.0)],
            chain_ids: vec![b' ', b' '],
            res_names: vec![res_name("HOH"), res_name("HOH")],
            res_nums: vec![100, 101],
            atom_names: vec![atom_name("O"), atom_name("O")],
            elements: vec![Element::O; 2],
        };
        let entities = split_into_entities(&coords);
        assert_eq!(entities.len(), 1);
        let water = &entities[0];
        assert_eq!(water.molecule_type(), MoleculeType::Water);
        assert!(water.as_bulk().is_some());
        assert_eq!(water.atom_count(), 2);
        assert_eq!(water.residue_count(), 2);
        assert!(water.label().contains("Water"));
    }

    #[test]
    fn ion_entity_accessors() {
        let coords = Coords {
            num_atoms: 1,
            atoms: vec![make_atom(5.0, 6.0, 7.0)],
            chain_ids: vec![b'Z'],
            res_names: vec![res_name("ZN")],
            res_nums: vec![200],
            atom_names: vec![atom_name("ZN")],
            elements: vec![Element::Zn],
        };
        let entities = split_into_entities(&coords);
        assert_eq!(entities.len(), 1);
        let ion = &entities[0];
        assert_eq!(ion.molecule_type(), MoleculeType::Ion);
        assert!(ion.as_small_molecule().is_some());
        assert_eq!(ion.residue_count(), 1);
    }

    #[test]
    fn is_focusable_for_protein_and_water() {
        let coords = make_two_residue_protein();
        let entities = split_into_entities(&coords);
        // Protein is focusable
        assert!(entities[0].is_focusable());

        // Water is not focusable
        let water_coords = Coords {
            num_atoms: 1,
            atoms: vec![make_atom(0.0, 0.0, 0.0)],
            chain_ids: vec![b' '],
            res_names: vec![res_name("HOH")],
            res_nums: vec![1],
            atom_names: vec![atom_name("O")],
            elements: vec![Element::O],
        };
        let water_entities = split_into_entities(&water_coords);
        assert!(!water_entities[0].is_focusable());
    }

    #[test]
    fn pdb_chain_id_for_protein() {
        let coords = make_two_residue_protein();
        let entities = split_into_entities(&coords);
        assert_eq!(entities[0].pdb_chain_id(), Some(b'A'));
    }
}

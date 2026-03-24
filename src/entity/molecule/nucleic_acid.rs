//! Nucleic acid entity — a single DNA or RNA chain instance.

use std::collections::HashMap;

use glam::Vec3;

use super::atom::Atom;
use super::id::EntityId;
use super::traits::{Entity, Polymer};
use super::{MoleculeType, Residue};
use crate::ops::codec::Coords;

/// Pre-extracted ring geometry for a single nucleotide base.
#[derive(Debug, Clone)]
pub struct NucleotideRing {
    /// Hexagonal ring atom positions in order: N1, C2, N3, C4, C5, C6.
    pub hex_ring: Vec<Vec3>,
    /// Pentagonal ring for purines: C4, C5, N7, C8, N9 (empty for
    /// pyrimidines).
    pub pent_ring: Vec<Vec3>,
    /// NDB color for this base.
    pub color: [f32; 3],
    /// C1' sugar carbon position (for anchoring stem to backbone spline).
    pub c1_prime: Option<Vec3>,
}

/// A single DNA or RNA entity (one PDB chain instance).
#[derive(Debug, Clone)]
pub struct NAEntity {
    /// Unique entity identifier.
    pub id: EntityId,
    /// Molecule type (DNA or RNA).
    pub na_type: MoleculeType,
    /// All atoms in this entity.
    pub atoms: Vec<Atom>,
    /// Ordered residues (nucleotides).
    pub residues: Vec<Residue>,
    /// Backbone segment break indices.
    pub segment_breaks: Vec<usize>,
    /// PDB chain identifier byte.
    pub pdb_chain_id: u8,
}

impl Entity for NAEntity {
    fn id(&self) -> EntityId {
        self.id
    }
    fn molecule_type(&self) -> MoleculeType {
        self.na_type
    }
    fn atoms(&self) -> &[Atom] {
        &self.atoms
    }
}

const HEX_RING_ATOMS: &[&str] = &["N1", "C2", "N3", "C4", "C5", "C6"];
const PENT_RING_ATOMS: &[&str] = &["C4", "C5", "N7", "C8", "N9"];

fn ndb_base_color(res_name: &str) -> Option<[f32; 3]> {
    match res_name {
        "DA" | "A" | "ADE" | "RAD" => Some([0.85, 0.20, 0.20]),
        "DG" | "G" | "GUA" | "RGU" => Some([0.20, 0.80, 0.20]),
        "DC" | "C" | "CYT" | "RCY" => Some([0.90, 0.90, 0.20]),
        "DT" | "THY" => Some([0.20, 0.20, 0.85]),
        "DU" | "U" | "URA" => Some([0.20, 0.85, 0.85]),
        _ => None,
    }
}

fn is_purine(res_name: &str) -> bool {
    matches!(
        res_name,
        "DA" | "DG" | "DI" | "A" | "G" | "ADE" | "GUA" | "I" | "RAD" | "RGU"
    )
}

impl NAEntity {
    /// Construct from flat `Coords` atom indices during entity splitting.
    #[must_use]
    pub fn from_coords_indices(
        id: EntityId,
        na_type: MoleculeType,
        indices: &[usize],
        coords: &Coords,
        pdb_chain_id: u8,
    ) -> Self {
        let (atoms, residues) =
            super::extract_atom_set_and_residues(indices, coords);
        Self {
            id,
            na_type,
            atoms,
            residues,
            segment_breaks: Vec::new(), // TODO: compute from P-P distances
            pdb_chain_id,
        }
    }

    /// Extract phosphorus (P) atom positions, split into segments at
    /// gaps where consecutive P-P distance exceeds ~8 Å.
    #[must_use]
    #[allow(
        clippy::excessive_nesting,
        reason = "gap-splitting logic is inherently nested"
    )]
    pub fn extract_p_atom_segments(&self) -> Vec<Vec<Vec3>> {
        const MAX_PP_DIST_SQ: f32 = 8.0 * 8.0;

        let mut p_positions = Vec::new();
        for residue in &self.residues {
            for idx in residue.atom_range.clone() {
                let name = std::str::from_utf8(&self.atoms[idx].name)
                    .unwrap_or("")
                    .trim();
                if name == "P" {
                    let a = &self.atoms[idx];
                    p_positions.push(a.position);
                }
            }
        }

        // Split at large gaps (missing residues)
        let mut result = Vec::new();
        let mut segment = Vec::new();
        for pos in p_positions {
            if let Some(&prev) = segment.last() {
                if pos.distance_squared(prev) > MAX_PP_DIST_SQ {
                    if segment.len() >= 2 {
                        result.push(std::mem::take(&mut segment));
                    } else {
                        segment.clear();
                    }
                }
            }
            segment.push(pos);
        }
        if segment.len() >= 2 {
            result.push(segment);
        }

        result
    }

    /// Extract base ring geometry for each nucleotide residue.
    #[must_use]
    #[allow(clippy::too_many_lines)]
    pub fn extract_base_rings(&self) -> Vec<NucleotideRing> {
        let mut rings = Vec::new();

        for residue in &self.residues {
            let res_name =
                std::str::from_utf8(&residue.name).unwrap_or("").trim();

            let Some(color) = ndb_base_color(res_name) else {
                continue;
            };

            let mut atom_map: HashMap<String, Vec3> = HashMap::new();
            for idx in residue.atom_range.clone() {
                let name = std::str::from_utf8(&self.atoms[idx].name)
                    .unwrap_or("")
                    .trim()
                    .trim_matches('\0')
                    .to_owned();
                let a = &self.atoms[idx];
                let _ = atom_map.insert(name, a.position);
            }

            let hex_ring: Vec<Vec3> = HEX_RING_ATOMS
                .iter()
                .filter_map(|name| atom_map.get(*name).copied())
                .collect();
            if hex_ring.len() != 6 {
                continue;
            }

            let pent_ring = if is_purine(res_name) {
                let pent: Vec<Vec3> = PENT_RING_ATOMS
                    .iter()
                    .filter_map(|name| atom_map.get(*name).copied())
                    .collect();
                if pent.len() == 5 {
                    pent
                } else {
                    Vec::new()
                }
            } else {
                Vec::new()
            };

            let c1_prime =
                atom_map.get("C1'").or_else(|| atom_map.get("C1*")).copied();

            rings.push(NucleotideRing {
                hex_ring,
                pent_ring,
                color,
                c1_prime,
            });
        }

        rings
    }
}

impl Polymer for NAEntity {
    fn residues(&self) -> &[Residue] {
        &self.residues
    }
    fn segment_breaks(&self) -> &[usize] {
        &self.segment_breaks
    }
}

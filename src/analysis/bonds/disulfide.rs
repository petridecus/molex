//! Disulfide bond detection between cysteine residues.
//!
//! Detects Cys-Cys disulfide bridges by finding SG-SG atom pairs within
//! the expected distance range (~2.05 Å).

use crate::entity::molecule::atom::Atom;

/// Maximum SG-SG distance (angstroms) for a disulfide bond.
const MAX_SS_DISTANCE: f32 = 2.5;
/// Minimum SG-SG distance (angstroms) to avoid clashes.
const MIN_SS_DISTANCE: f32 = 1.5;

/// A disulfide bond between two cysteine residues.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DisulfideBond {
    /// Atom index of the first SG sulfur.
    pub sg_a: usize,
    /// Atom index of the second SG sulfur.
    pub sg_b: usize,
    /// SG-SG distance in angstroms.
    pub distance: f32,
}

/// Detect disulfide bonds from SG atom positions.
///
/// Finds all pairs of atoms named " SG " within the expected disulfide
/// distance range (1.5–2.5 Å).
#[must_use]
pub fn detect_disulfide_bonds(atoms: &[Atom]) -> Vec<DisulfideBond> {
    let sg_indices: Vec<usize> = atoms
        .iter()
        .enumerate()
        .filter(|(_, a)| a.name == *b" SG ")
        .map(|(i, _)| i)
        .collect();

    let mut bonds = Vec::new();
    for (ai, &i) in sg_indices.iter().enumerate() {
        for &j in &sg_indices[ai + 1..] {
            let dist = atoms[i].position.distance(atoms[j].position);
            if (MIN_SS_DISTANCE..=MAX_SS_DISTANCE).contains(&dist) {
                bonds.push(DisulfideBond {
                    sg_a: i,
                    sg_b: j,
                    distance: dist,
                });
            }
        }
    }
    bonds
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use glam::Vec3;

    use super::*;
    use crate::element::Element;

    fn make_sg(x: f32, y: f32, z: f32) -> Atom {
        Atom {
            position: Vec3::new(x, y, z),
            occupancy: 1.0,
            b_factor: 0.0,
            element: Element::S,
            name: *b" SG ",
        }
    }

    #[test]
    fn detect_single_disulfide() {
        let atoms = vec![make_sg(0.0, 0.0, 0.0), make_sg(2.03, 0.0, 0.0)];
        let bonds = detect_disulfide_bonds(&atoms);
        assert_eq!(bonds.len(), 1);
        assert!((bonds[0].distance - 2.03).abs() < 1e-4);
    }

    #[test]
    fn no_disulfide_too_far() {
        let atoms = vec![make_sg(0.0, 0.0, 0.0), make_sg(5.0, 0.0, 0.0)];
        assert!(detect_disulfide_bonds(&atoms).is_empty());
    }

    #[test]
    fn ignores_non_sg() {
        let atoms = vec![
            Atom {
                position: Vec3::new(0.0, 0.0, 0.0),
                occupancy: 1.0,
                b_factor: 0.0,
                element: Element::S,
                name: *b" SD ",
            },
            make_sg(2.03, 0.0, 0.0),
        ];
        assert!(detect_disulfide_bonds(&atoms).is_empty());
    }
}

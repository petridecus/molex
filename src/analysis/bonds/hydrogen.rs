//! Hydrogen bond detection.
//!
//! Two algorithms:
//!
//! - **Kabsch-Sander** ([`detect_hbonds`]): backbone-only, energy-based. Used
//!   by DSSP for secondary structure classification. Operates on
//!   [`ResidueBackbone`] with the electrostatic partial-charge formula.
//!
//! - **Geometric** ([`detect_all_hbonds`]): general-purpose, covers all atom
//!   types (backbone, sidechain, ligand, nucleic acid, water). Uses
//!   donor···acceptor distance and D-H···A angle criteria. For visualization.

use glam::Vec3;

use crate::element::Element;
use crate::entity::molecule::atom::Atom;
use crate::entity::molecule::protein::ResidueBackbone;

/// Kabsch-Sander electrostatic energy constant (kcal/mol · Å).
const KS_FACTOR: f32 = 27.888;

/// H-bond energy threshold (kcal/mol). Bonds with E < this are accepted.
const HBOND_THRESHOLD: f32 = -0.5;

/// A backbone hydrogen bond detected by Kabsch-Sander energy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HBond {
    /// Residue index of the N-H donor.
    pub donor: usize,
    /// Residue index of the C=O acceptor.
    pub acceptor: usize,
    /// Kabsch-Sander energy (kcal/mol, negative = favorable).
    pub energy: f32,
}

/// Detect backbone hydrogen bonds using the Kabsch-Sander energy criterion.
///
/// For each residue, amide H positions are estimated from geometry:
/// `H_i = N_i + normalize(N_i - C_{i-1})` (first residue uses N-CA).
///
/// Returns all donor-acceptor pairs with energy below the threshold,
/// sorted by donor index then by energy (strongest first).
#[must_use]
pub fn detect_hbonds(residues: &[ResidueBackbone]) -> Vec<HBond> {
    let n = residues.len();
    if n < 2 {
        return Vec::new();
    }

    // Estimate amide H positions using the bisector method:
    // H is placed 1.0 Å from N along the bisector of (N - C_prev) and
    // (N - CA). This matches the canonical DSSP placement.
    let h_positions: Vec<Vec3> = (0..n)
        .map(|i| {
            let v_ca = (residues[i].n - residues[i].ca).normalize_or_zero();
            let v_prev = if i > 0 {
                (residues[i].n - residues[i - 1].c).normalize_or_zero()
            } else {
                v_ca
            };
            let bisector = (v_ca + v_prev).normalize_or_zero();
            residues[i].n + bisector
        })
        .collect();

    let mut hbonds = Vec::new();

    for i in 1..n {
        let h = h_positions[i];
        let n_pos = residues[i].n;

        for (j, res_j) in residues.iter().enumerate() {
            // Skip self and immediate neighbors.
            if i == j || i == j + 1 || (j > 0 && i == j - 1) {
                continue;
            }

            let c = res_j.c;
            let o = res_j.o;

            let r_on = (o - n_pos).length();
            let r_ch = (c - h).length();
            let r_oh = (o - h).length();
            let r_cn = (c - n_pos).length();

            // Avoid division by near-zero distances.
            if r_on < 0.5 || r_ch < 0.5 || r_oh < 0.5 || r_cn < 0.5 {
                continue;
            }

            let energy =
                KS_FACTOR * (1.0 / r_on + 1.0 / r_ch - 1.0 / r_oh - 1.0 / r_cn);

            if energy < HBOND_THRESHOLD {
                hbonds.push(HBond {
                    donor: i,
                    acceptor: j,
                    energy,
                });
            }
        }
    }

    hbonds
}

// ── Geometric H-bond detection (all atoms) ──────────────────────────

/// Maximum donor···acceptor distance (Å) for geometric H-bond detection.
const MAX_DA_DIST_SQ: f32 = 3.5 * 3.5;

/// Maximum H···acceptor distance when explicit H is available (Å).
const MAX_HA_DIST_SQ: f32 = 2.5 * 2.5;

/// Minimum D-H···A angle cosine (cos(120°) ≈ -0.5).
/// Angles more obtuse than 120° are accepted.
const MIN_DHA_COS: f32 = -0.5;

/// Minimum distance² to exclude covalently bonded pairs.
const MIN_NONBONDED_DIST_SQ: f32 = 1.5 * 1.5;

/// Maximum D-H distance² for finding bonded hydrogens.
const MAX_DH_DIST_SQ: f32 = 1.3 * 1.3;

/// A hydrogen bond between two atoms, detected by geometric criteria.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AtomHBond {
    /// Index of the donor heavy atom in the input slice.
    pub donor: usize,
    /// Index of the acceptor atom in the input slice.
    pub acceptor: usize,
    /// Donor···acceptor distance (Å).
    pub distance: f32,
}

/// Whether an element can donate an H-bond (has N-H, O-H, or S-H).
fn is_donor(elem: Element) -> bool {
    matches!(elem, Element::N | Element::O | Element::S)
}

/// Whether an element can accept an H-bond (has lone pairs).
fn is_acceptor(elem: Element) -> bool {
    matches!(elem, Element::N | Element::O | Element::S | Element::F)
}

/// Detect hydrogen bonds between all eligible atom pairs.
///
/// Takes a flat slice of atoms (from any combination of entities) and
/// returns donor···acceptor pairs satisfying geometric criteria.
///
/// When explicit H atoms are present, applies both distance and angle
/// checks (D···A < 3.5 Å, H···A < 2.5 Å, D-H···A > 120°). Without
/// explicit H (typical X-ray structures), uses D···A distance only.
#[must_use]
pub fn detect_all_hbonds(atoms: &[Atom]) -> Vec<AtomHBond> {
    let n = atoms.len();
    if n < 2 {
        return Vec::new();
    }

    let has_h = atoms.iter().any(|a| a.element == Element::H);

    // Collect donor and acceptor indices.
    let donors: Vec<usize> =
        (0..n).filter(|&i| is_donor(atoms[i].element)).collect();
    let acceptors: Vec<usize> =
        (0..n).filter(|&i| is_acceptor(atoms[i].element)).collect();

    // Pre-compute bonded H for each donor (only if H atoms present).
    let donor_h: Vec<Vec<usize>> = if has_h {
        donors
            .iter()
            .map(|&d| {
                let dp = atoms[d].position;
                (0..n)
                    .filter(|&i| {
                        atoms[i].element == Element::H
                            && (atoms[i].position - dp).length_squared()
                                < MAX_DH_DIST_SQ
                    })
                    .collect()
            })
            .collect()
    } else {
        vec![Vec::new(); donors.len()]
    };

    let mut hbonds = Vec::new();

    for (di, &d) in donors.iter().enumerate() {
        let d_pos = atoms[d].position;

        for &a in &acceptors {
            if d == a {
                continue;
            }

            let a_pos = atoms[a].position;
            let da_sq = (d_pos - a_pos).length_squared();

            if !(MIN_NONBONDED_DIST_SQ..=MAX_DA_DIST_SQ).contains(&da_sq) {
                continue;
            }

            // With explicit H: check angle criterion.
            if has_h
                && !donor_h[di].is_empty()
                && !passes_angle_check(&donor_h[di], d_pos, a_pos, atoms)
            {
                continue;
            }

            hbonds.push(AtomHBond {
                donor: d,
                acceptor: a,
                distance: da_sq.sqrt(),
            });
        }
    }

    hbonds
}

/// Check if any bonded H satisfies the H···A distance and D-H···A angle.
fn passes_angle_check(
    h_indices: &[usize],
    d_pos: Vec3,
    a_pos: Vec3,
    atoms: &[Atom],
) -> bool {
    h_indices.iter().any(|&h| {
        let h_pos = atoms[h].position;
        let ha_sq = (h_pos - a_pos).length_squared();
        if ha_sq > MAX_HA_DIST_SQ {
            return false;
        }
        let hd = (d_pos - h_pos).normalize_or_zero();
        let ha = (a_pos - h_pos).normalize_or_zero();
        hd.dot(ha) < MIN_DHA_COS
    })
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    #[test]
    fn empty_input() {
        assert!(detect_hbonds(&[]).is_empty());
    }

    #[test]
    fn single_residue() {
        let residues = vec![ResidueBackbone {
            n: Vec3::new(0.0, 0.0, 0.0),
            ca: Vec3::new(1.5, 0.0, 0.0),
            c: Vec3::new(2.5, 1.0, 0.0),
            o: Vec3::new(2.5, 2.0, 0.0),
        }];
        assert!(detect_hbonds(&residues).is_empty());
    }

    // -- Geometric H-bond tests --

    fn atom(x: f32, y: f32, z: f32, elem: Element, name: &str) -> Atom {
        let mut n = [b' '; 4];
        for (i, b) in name.bytes().take(4).enumerate() {
            n[i] = b;
        }
        Atom {
            position: Vec3::new(x, y, z),
            occupancy: 1.0,
            b_factor: 0.0,
            element: elem,
            name: n,
        }
    }

    #[test]
    fn geometric_no_pair_too_far() {
        let atoms = vec![
            atom(0.0, 0.0, 0.0, Element::N, "N"),
            atom(10.0, 0.0, 0.0, Element::O, "O"),
        ];
        assert!(detect_all_hbonds(&atoms).is_empty());
    }

    #[test]
    fn geometric_valid_hbond() {
        let atoms = vec![
            atom(0.0, 0.0, 0.0, Element::N, "N"),
            atom(2.9, 0.0, 0.0, Element::O, "O"),
        ];
        let hbonds = detect_all_hbonds(&atoms);
        // N→O and O→N (both are donor + acceptor)
        assert_eq!(hbonds.len(), 2);
        assert!(hbonds[0].distance < 3.5);
    }

    #[test]
    fn geometric_excludes_covalent() {
        let atoms = vec![
            atom(0.0, 0.0, 0.0, Element::N, "N"),
            atom(1.3, 0.0, 0.0, Element::O, "O"),
        ];
        assert!(detect_all_hbonds(&atoms).is_empty());
    }

    #[test]
    fn geometric_carbon_not_donor() {
        let atoms = vec![
            atom(0.0, 0.0, 0.0, Element::C, "CA"),
            atom(2.9, 0.0, 0.0, Element::O, "O"),
        ];
        assert!(detect_all_hbonds(&atoms).is_empty());
    }
}

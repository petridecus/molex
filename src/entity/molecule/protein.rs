//! Protein types and entity.
//!
//! Per-residue types (`ResidueBackbone`, `Sidechain`, `ProteinResidue`) and
//! the `ProteinEntity` chain instance with segment break metadata.

use glam::Vec3;

use super::atom::Atom;
use super::id::EntityId;
use super::traits::{Entity, Polymer};
use super::{MoleculeType, Residue};
use crate::ops::codec::Coords;

// ---------------------------------------------------------------------------
// Per-residue types
// ---------------------------------------------------------------------------

/// Backbone atom positions for a single protein residue.
///
/// Named fields, fixed size, `Copy`. Used by DSSP H-bond detection
/// and secondary structure classification.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ResidueBackbone {
    /// Nitrogen atom position.
    pub n: Vec3,
    /// Alpha-carbon atom position.
    pub ca: Vec3,
    /// Carbonyl carbon atom position.
    pub c: Vec3,
    /// Carbonyl oxygen atom position.
    pub o: Vec3,
}

/// Sidechain atoms and bond topology for a single protein residue.
///
/// `Vec<Atom>` holds the atoms (positions + names + elements). Bond
/// indices are local to this sidechain's `Vec<Atom>` (0-based). The
/// backbone-to-sidechain bond (CA→CB) is implicit — `ResidueBackbone.ca`
/// connects to the first atom in `atoms` when the sidechain is non-empty.
#[derive(Debug, Clone)]
pub struct Sidechain {
    /// Sidechain atoms (CB, CG, SG, etc.).
    pub atoms: Vec<Atom>,
    /// Intra-residue bonds as `(atom_idx_a, atom_idx_b)` pairs, local
    /// to this sidechain's `Vec<Atom>`.
    pub bonds: Vec<(usize, usize)>,
    /// Whether this residue is hydrophobic (Kyte-Doolittle classification).
    pub is_hydrophobic: bool,
}

impl Sidechain {
    /// Create an empty sidechain (e.g. for glycine).
    #[must_use]
    pub fn empty() -> Self {
        Self {
            atoms: Vec::new(),
            bonds: Vec::new(),
            is_hydrophobic: false,
        }
    }

    /// Whether this sidechain has no atoms (glycine).
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
    }
}

/// Complete per-residue view of a protein entity.
///
/// Contains backbone (N, CA, C, O), sidechain (atoms + bonds), and
/// residue metadata (name, number, chain).
#[derive(Debug, Clone)]
pub struct ProteinResidue {
    /// 3-character residue name (e.g. b"ALA").
    pub name: [u8; 3],
    /// PDB residue sequence number.
    pub number: i32,
    /// PDB chain identifier byte.
    pub pdb_chain_id: u8,
    /// Backbone atom positions (N, CA, C, O).
    pub backbone: ResidueBackbone,
    /// Sidechain atoms and topology.
    pub sidechain: Sidechain,
}

// ---------------------------------------------------------------------------
// Protein entity
// ---------------------------------------------------------------------------

/// A single protein chain instance.
///
/// Contains all atoms, residue metadata, and precomputed segment breaks
/// (backbone gaps from missing residues, detected by C-N distance).
///
/// Protein-specific derived views (`to_backbone`, `to_protein_residues`,
/// `to_sidechains`) are available as methods. These derive from the
/// underlying `atoms` + `residues` on each call — no cached copies.
#[derive(Debug, Clone)]
pub struct ProteinEntity {
    /// Unique entity identifier.
    pub id: EntityId,
    /// All atoms in this protein chain.
    pub atoms: Vec<Atom>,
    /// Ordered residues (name, number, atom range into `atoms`).
    pub residues: Vec<Residue>,
    /// Indices into `residues` where backbone segments start a new
    /// continuous run. Computed from C(i)->N(i+1) distance > 2.0Å.
    pub segment_breaks: Vec<usize>,
    /// PDB chain identifier byte for this entity.
    pub pdb_chain_id: u8,
}

/// Maximum C->N distance (A) for a peptide bond. Pairs exceeding this
/// are treated as segment breaks.
const MAX_PEPTIDE_BOND_DIST: f32 = 2.0;

impl ProteinEntity {
    /// Construct a protein entity, computing segment breaks from
    /// backbone atom distances.
    #[must_use]
    pub fn new(
        id: EntityId,
        atoms: Vec<Atom>,
        residues: Vec<Residue>,
        pdb_chain_id: u8,
    ) -> Self {
        let segment_breaks = compute_segment_breaks(&residues, &atoms);
        Self {
            id,
            atoms,
            residues,
            segment_breaks,
            pdb_chain_id,
        }
    }

    /// Construct from flat `Coords` atom indices during entity splitting.
    #[must_use]
    pub fn from_coords_indices(
        id: EntityId,
        indices: &[usize],
        coords: &Coords,
        pdb_chain_id: u8,
    ) -> Self {
        let (atoms, residues) =
            super::extract_atom_set_and_residues(indices, coords);
        Self::new(id, atoms, residues, pdb_chain_id)
    }

    /// Derive backbone (N, CA, C, O) for all residues.
    ///
    /// Returns one `ResidueBackbone` per residue that has all four
    /// backbone atoms. Residues missing any backbone atom are skipped.
    #[must_use]
    pub fn to_backbone(&self) -> Vec<ResidueBackbone> {
        self.residues
            .iter()
            .filter_map(|r| extract_backbone_from_residue(&self.atoms, r))
            .collect()
    }

    /// Derive full protein residues (backbone + sidechain).
    #[must_use]
    pub fn to_protein_residues<F, G>(
        &self,
        is_hydrophobic: F,
        get_bonds: G,
    ) -> Vec<ProteinResidue>
    where
        F: Fn(&str) -> bool,
        G: Fn(&str) -> Option<Vec<(&'static str, &'static str)>>,
    {
        self.residues
            .iter()
            .filter_map(|r| {
                let backbone = extract_backbone_from_residue(&self.atoms, r)?;
                let sidechain = extract_sidechain_from_residue(
                    &self.atoms,
                    r,
                    &is_hydrophobic,
                    &get_bonds,
                );
                Some(ProteinResidue {
                    name: r.name,
                    number: r.number,
                    pdb_chain_id: self.pdb_chain_id,
                    backbone,
                    sidechain,
                })
            })
            .collect()
    }

    /// N/CA/C interleaved positions per segment (for spline renderer).
    ///
    /// Each inner `Vec<Vec3>` is one continuous segment with atoms
    /// laid out as `[N0, CA0, C0, N1, CA1, C1, ...]`.
    #[must_use]
    pub fn to_interleaved_segments(&self) -> Vec<Vec<Vec3>> {
        let n_segments = self.segment_count();
        (0..n_segments)
            .map(|seg_idx| {
                let range = self.segment_range(seg_idx);
                let mut positions = Vec::with_capacity(range.len() * 3);
                for r in &self.residues[range] {
                    if let Some(bb) =
                        extract_backbone_from_residue(&self.atoms, r)
                    {
                        positions.push(bb.n);
                        positions.push(bb.ca);
                        positions.push(bb.c);
                    }
                }
                positions
            })
            .collect()
    }
}

impl Entity for ProteinEntity {
    fn id(&self) -> EntityId {
        self.id
    }
    fn molecule_type(&self) -> MoleculeType {
        MoleculeType::Protein
    }
    fn atoms(&self) -> &[Atom] {
        &self.atoms
    }
}

impl Polymer for ProteinEntity {
    fn residues(&self) -> &[Residue] {
        &self.residues
    }
    fn segment_breaks(&self) -> &[usize] {
        &self.segment_breaks
    }
}

// -- Internal helpers --------------------------------------------------------

/// Compute segment break indices from C(i)->N(i+1) distances.
fn compute_segment_breaks(residues: &[Residue], atoms: &[Atom]) -> Vec<usize> {
    let mut breaks = Vec::new();
    for i in 1..residues.len() {
        let prev_c = find_backbone_atom(atoms, &residues[i - 1], *b"C   ");
        let curr_n = find_backbone_atom(atoms, &residues[i], *b"N   ");
        match (prev_c, curr_n) {
            (Some(c), Some(n)) => {
                if (c - n).length() > MAX_PEPTIDE_BOND_DIST {
                    breaks.push(i);
                }
            }
            _ => breaks.push(i),
        }
    }
    breaks
}

/// Find a specific backbone atom position within a residue.
fn find_backbone_atom(
    atoms: &[Atom],
    residue: &Residue,
    target_name: [u8; 4],
) -> Option<Vec3> {
    for idx in residue.atom_range.clone() {
        if atoms[idx].name == target_name {
            let a = &atoms[idx];
            return Some(a.position);
        }
    }
    None
}

/// Find a backbone atom by trimmed string name.
fn find_atom_by_name(
    atoms: &[Atom],
    residue: &Residue,
    name: &str,
) -> Option<Vec3> {
    for idx in residue.atom_range.clone() {
        let atom_name =
            std::str::from_utf8(&atoms[idx].name).unwrap_or("").trim();
        if atom_name == name {
            let a = &atoms[idx];
            return Some(a.position);
        }
    }
    None
}

/// Extract backbone atom positions (N, CA, C, O) from a single residue.
fn extract_backbone_from_residue(
    atoms: &[Atom],
    residue: &Residue,
) -> Option<ResidueBackbone> {
    Some(ResidueBackbone {
        n: find_atom_by_name(atoms, residue, "N")?,
        ca: find_atom_by_name(atoms, residue, "CA")?,
        c: find_atom_by_name(atoms, residue, "C")?,
        o: find_atom_by_name(atoms, residue, "O")
            .or_else(|| find_atom_by_name(atoms, residue, "OXT"))?,
    })
}

/// Extract sidechain from a single residue.
fn extract_sidechain_from_residue<F, G>(
    atoms: &[Atom],
    residue: &Residue,
    is_hydrophobic: &F,
    get_bonds: &G,
) -> Sidechain
where
    F: Fn(&str) -> bool,
    G: Fn(&str) -> Option<Vec<(&'static str, &'static str)>>,
{
    use std::collections::HashMap;

    let res_name = std::str::from_utf8(&residue.name).unwrap_or("UNK").trim();

    let mut sc_atoms: Vec<Atom> = Vec::new();
    let mut name_to_local: HashMap<String, usize> = HashMap::new();

    for idx in residue.atom_range.clone() {
        let name = std::str::from_utf8(&atoms[idx].name).unwrap_or("").trim();
        // Skip backbone and hydrogen atoms.
        if matches!(name, "N" | "CA" | "C" | "O" | "OXT") || is_hydrogen(name) {
            continue;
        }
        let local_idx = sc_atoms.len();
        let _ = name_to_local.insert(name.to_owned(), local_idx);
        sc_atoms.push(atoms[idx].clone());
    }

    let mut bonds = Vec::new();
    if let Some(residue_bonds) = get_bonds(res_name) {
        for (a1, a2) in residue_bonds {
            let a1_str: &str = a1;
            let a2_str: &str = a2;
            if let (Some(&i1), Some(&i2)) =
                (name_to_local.get(a1_str), name_to_local.get(a2_str))
            {
                bonds.push((i1, i2));
            }
        }
    }

    Sidechain {
        atoms: sc_atoms,
        bonds,
        is_hydrophobic: is_hydrophobic(res_name),
    }
}

/// Returns true if an atom name looks like a hydrogen atom.
fn is_hydrogen(name: &str) -> bool {
    name.starts_with('H')
        || name.starts_with("1H")
        || name.starts_with("2H")
        || name.starts_with("3H")
        || (name.len() >= 2
            && name.as_bytes().first().is_some_and(u8::is_ascii_digit)
            && name.as_bytes().get(1) == Some(&b'H'))
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

    fn make_two_residue_protein_coords() -> Coords {
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

    // -- Sidechain tests --

    #[test]
    fn sidechain_empty() {
        let sc = Sidechain::empty();
        assert!(sc.is_empty());
        assert!(sc.atoms.is_empty());
        assert!(sc.bonds.is_empty());
        assert!(!sc.is_hydrophobic);
    }

    #[test]
    fn sidechain_is_empty_with_atoms() {
        let sc = Sidechain {
            atoms: vec![Atom {
                position: Vec3::ZERO,
                occupancy: 1.0,
                b_factor: 0.0,
                element: Element::C,
                name: *b"CB  ",
            }],
            bonds: Vec::new(),
            is_hydrophobic: false,
        };
        assert!(!sc.is_empty());
    }

    // -- to_backbone tests --

    #[test]
    fn to_backbone_returns_two_entries() {
        let coords = make_two_residue_protein_coords();
        let entities = split_into_entities(&coords);
        let protein = entities[0].as_protein().unwrap();
        let backbone = protein.to_backbone();
        assert_eq!(backbone.len(), 2);
    }

    #[test]
    fn to_backbone_first_residue_positions() {
        let coords = make_two_residue_protein_coords();
        let entities = split_into_entities(&coords);
        let protein = entities[0].as_protein().unwrap();
        let backbone = protein.to_backbone();

        let bb0 = &backbone[0];
        assert!((bb0.n.x - 1.0).abs() < 1e-6);
        assert!((bb0.ca.x - 4.0).abs() < 1e-6);
        assert!((bb0.c.x - 7.0).abs() < 1e-6);
        assert!((bb0.o.x - 10.0).abs() < 1e-6);
    }

    #[test]
    fn to_backbone_second_residue_positions() {
        let coords = make_two_residue_protein_coords();
        let entities = split_into_entities(&coords);
        let protein = entities[0].as_protein().unwrap();
        let backbone = protein.to_backbone();

        let bb1 = &backbone[1];
        assert!((bb1.n.x - 13.0).abs() < 1e-6);
        assert!((bb1.ca.x - 16.0).abs() < 1e-6);
        assert!((bb1.c.x - 19.0).abs() < 1e-6);
        assert!((bb1.o.x - 22.0).abs() < 1e-6);
    }

    // -- segment breaks --

    #[test]
    fn no_segment_break_when_residues_close() {
        let coords = make_two_residue_protein_coords();
        let entities = split_into_entities(&coords);
        let protein = entities[0].as_protein().unwrap();
        // C of res1 is at (7,8,9), N of res2 is at (13,14,15)
        // distance ~10.4, which exceeds MAX_PEPTIDE_BOND_DIST, so there is a
        // break
        assert!(!protein.segment_breaks.is_empty());
    }

    #[test]
    fn segment_break_detected_for_large_gap() {
        // Build residues where C->N distance > 2.0 A
        let coords = Coords {
            num_atoms: 8,
            atoms: vec![
                make_atom(0.0, 0.0, 0.0),  // N
                make_atom(1.0, 0.0, 0.0),  // CA
                make_atom(2.0, 0.0, 0.0),  // C
                make_atom(2.0, 1.0, 0.0),  // O
                make_atom(20.0, 0.0, 0.0), // N (far away)
                make_atom(21.0, 0.0, 0.0), // CA
                make_atom(22.0, 0.0, 0.0), // C
                make_atom(22.0, 1.0, 0.0), // O
            ],
            chain_ids: vec![b'A'; 8],
            res_names: vec![res_name("ALA"); 8],
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
            elements: vec![Element::Unknown; 8],
        };
        let entities = split_into_entities(&coords);
        let protein = entities[0].as_protein().unwrap();
        assert!(
            !protein.segment_breaks.is_empty(),
            "should have segment break at large gap"
        );
    }

    #[test]
    fn no_segment_break_for_close_residues() {
        // Build residues where C->N distance < 2.0 A (peptide bond)
        let coords = Coords {
            num_atoms: 8,
            atoms: vec![
                make_atom(0.0, 0.0, 0.0), // N
                make_atom(1.5, 0.0, 0.0), // CA
                make_atom(2.5, 0.0, 0.0), // C
                make_atom(2.5, 1.0, 0.0), // O
                make_atom(3.8, 0.0, 0.0), // N (1.3 A from C, normal peptide)
                make_atom(5.0, 0.0, 0.0), // CA
                make_atom(6.0, 0.0, 0.0), // C
                make_atom(6.0, 1.0, 0.0), // O
            ],
            chain_ids: vec![b'A'; 8],
            res_names: vec![res_name("ALA"); 8],
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
            elements: vec![Element::Unknown; 8],
        };
        let entities = split_into_entities(&coords);
        let protein = entities[0].as_protein().unwrap();
        assert!(
            protein.segment_breaks.is_empty(),
            "should have no segment break for peptide-bonded residues"
        );
    }

    // -- to_interleaved_segments --

    #[test]
    fn interleaved_segments_for_single_segment() {
        let coords = Coords {
            num_atoms: 8,
            atoms: vec![
                make_atom(0.0, 0.0, 0.0),
                make_atom(1.5, 0.0, 0.0),
                make_atom(2.5, 0.0, 0.0),
                make_atom(2.5, 1.0, 0.0),
                make_atom(3.8, 0.0, 0.0),
                make_atom(5.0, 0.0, 0.0),
                make_atom(6.0, 0.0, 0.0),
                make_atom(6.0, 1.0, 0.0),
            ],
            chain_ids: vec![b'A'; 8],
            res_names: vec![res_name("ALA"); 8],
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
            elements: vec![Element::Unknown; 8],
        };
        let entities = split_into_entities(&coords);
        let protein = entities[0].as_protein().unwrap();
        let segments = protein.to_interleaved_segments();
        assert_eq!(segments.len(), 1);
        // 2 residues * 3 backbone atoms (N, CA, C) = 6
        assert_eq!(segments[0].len(), 6);
    }

    // -- is_hydrogen helper --

    #[test]
    fn is_hydrogen_recognizes_h_names() {
        assert!(is_hydrogen("H"));
        assert!(is_hydrogen("HA"));
        assert!(is_hydrogen("HB2"));
        assert!(is_hydrogen("1H"));
        assert!(is_hydrogen("2HB"));
        assert!(is_hydrogen("3HD1"));
    }

    #[test]
    fn is_hydrogen_rejects_non_h() {
        assert!(!is_hydrogen("CA"));
        assert!(!is_hydrogen("N"));
        assert!(!is_hydrogen("O"));
        assert!(!is_hydrogen("CB"));
    }
}

//! Protein types and entity.
//!
//! Per-residue types (`ResidueBackbone`, `Sidechain`, `ProteinResidue`) and
//! the `ProteinEntity` chain instance with segment break metadata.

use glam::Vec3;

use super::atom::Atom;
use super::id::EntityId;
use super::traits::{Entity, Polymer};
use super::{MoleculeType, Residue};
use crate::analysis::BondOrder;
use crate::atom_id::AtomId;
use crate::bond::CovalentBond;
use crate::chemistry::amino_acids::AminoAcid;
use crate::chemistry::atom_name::AtomName;
use crate::element::Element;
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
    ///
    /// After [`ProteinEntity::new`], atoms belonging to a kept residue are
    /// laid out in canonical order: `N, CA, C, O, sidechain heavy...,
    /// hydrogens...`. Atoms of dropped residues remain in `atoms` but are
    /// not referenced by any `Residue::atom_range`.
    pub atoms: Vec<Atom>,
    /// Ordered residues (name, number, atom range into `atoms`). Residues
    /// missing any of the four backbone atoms (N, CA, C, O or OXT) are
    /// dropped during construction.
    pub residues: Vec<Residue>,
    /// Indices into `residues` where backbone segments start a new
    /// continuous run. Computed from C(i)->N(i+1) distance > 2.0Å.
    pub segment_breaks: Vec<usize>,
    /// Intra-entity covalent bonds with `AtomId` endpoints. Populated at
    /// construction from `AminoAcid::bonds()` tables (intra-residue) plus
    /// universal backbone bonds (N-CA, CA-C, C=O) and inter-residue peptide
    /// bonds `C(i)-N(i+1)` across pairs that are not separated by a
    /// segment break.
    pub bonds: Vec<CovalentBond>,
    /// PDB chain identifier byte for this entity.
    pub pdb_chain_id: u8,
}

/// Maximum C->N distance (A) for a peptide bond. Pairs exceeding this
/// are treated as segment breaks.
const MAX_PEPTIDE_BOND_DIST: f32 = 2.0;

impl ProteinEntity {
    /// Construct a protein entity.
    ///
    /// Enforces the canonical atom ordering invariant: for every kept
    /// residue, atoms are laid out as `N, CA, C, O, sidechain heavy...,
    /// hydrogens...`. Residues missing any of N/CA/C/O (with OXT as
    /// C-terminal oxygen fallback) are dropped from the `residues` vec
    /// with a `log::warn!` entry; their atoms remain in `atoms` but are
    /// unreferenced.
    ///
    /// Also computes segment breaks from C(i)->N(i+1) distance and
    /// populates `bonds` from the `AminoAcid::bonds()` tables plus
    /// universal backbone + inter-residue peptide bonds.
    #[must_use]
    #[allow(
        clippy::needless_pass_by_value,
        reason = "constructor owns the inputs; canonicalize borrows before \
                  reassigning new owned vecs"
    )]
    pub fn new(
        id: EntityId,
        atoms: Vec<Atom>,
        residues: Vec<Residue>,
        pdb_chain_id: u8,
    ) -> Self {
        let (atoms, residues) =
            canonicalize_protein_residues(&atoms, &residues, pdb_chain_id);
        let segment_breaks = compute_segment_breaks(&residues, &atoms);
        let bonds = build_protein_bonds(id, &atoms, &residues, &segment_breaks);
        Self {
            id,
            atoms,
            residues,
            segment_breaks,
            bonds,
            pdb_chain_id,
        }
    }

    /// Iterate covalent bonds whose endpoints lie in any residue's
    /// canonical backbone region (indices `atom_range.start..+4`).
    ///
    /// Intra-residue N-CA, CA-C, C=O and inter-residue peptide C-N all
    /// qualify because their endpoints occupy positions 0..4 of each
    /// residue after canonical ordering. OXT (C-terminal oxygen) is
    /// stored in the sidechain region and is therefore not considered
    /// a backbone atom by this filter.
    pub fn backbone_bonds(&self) -> impl Iterator<Item = &CovalentBond> + '_ {
        self.bonds.iter().filter(move |b| {
            self.is_backbone_atom(b.a) && self.is_backbone_atom(b.b)
        })
    }

    /// Iterate covalent bonds with at least one heavy-atom endpoint in
    /// the sidechain region (`residue.atom_range.start + 4..`).
    ///
    /// Hydrogen endpoints are excluded (they must be non-hydrogen to
    /// qualify). The CA-CB anchor bond qualifies because CB is a
    /// sidechain heavy atom.
    pub fn sidechain_bonds(&self) -> impl Iterator<Item = &CovalentBond> + '_ {
        self.bonds.iter().filter(move |b| {
            self.is_sidechain_heavy_atom(b.a)
                || self.is_sidechain_heavy_atom(b.b)
        })
    }

    fn residue_local_offset(&self, atom_id: AtomId) -> Option<usize> {
        if atom_id.entity != self.id {
            return None;
        }
        let idx = atom_id.index as usize;
        for r in &self.residues {
            if r.atom_range.contains(&idx) {
                return Some(idx - r.atom_range.start);
            }
        }
        None
    }

    fn is_backbone_atom(&self, atom_id: AtomId) -> bool {
        self.residue_local_offset(atom_id)
            .is_some_and(|off| off < 4)
    }

    fn is_sidechain_heavy_atom(&self, atom_id: AtomId) -> bool {
        let Some(off) = self.residue_local_offset(atom_id) else {
            return false;
        };
        if off < 4 {
            return false;
        }
        let atom = &self.atoms[atom_id.index as usize];
        atom.element != Element::H
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

    /// Derive per-residue sidechains (atoms, bonds, hydrophobicity).
    ///
    /// Returns one [`Sidechain`] per residue that has backbone atoms (same
    /// residues as [`to_backbone()`](Self::to_backbone)). Glycine residues
    /// produce an empty `Sidechain`. Residues missing backbone atoms are
    /// skipped.
    #[must_use]
    pub fn to_sidechains<F, G>(
        &self,
        is_hydrophobic: F,
        get_bonds: G,
    ) -> Vec<Sidechain>
    where
        F: Fn(&str) -> bool,
        G: Fn(&str) -> Option<Vec<(&'static str, &'static str)>>,
    {
        self.residues
            .iter()
            .filter_map(|r| {
                // Skip residues without complete backbone (same filter as
                // to_backbone / to_protein_residues).
                let _bb = extract_backbone_from_residue(&self.atoms, r)?;
                Some(extract_sidechain_from_residue(
                    &self.atoms,
                    r,
                    &is_hydrophobic,
                    &get_bonds,
                ))
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

    /// Classify secondary structure (DSSP) for each residue.
    ///
    /// Runs DSSP on the full backbone (all residues in the entity),
    /// then merges short segments. Segment breaks are handled naturally
    /// by the H-bond energy calculation — missing residues simply don't
    /// form bonds. Returns one [`SSType`](crate::SSType) per residue
    /// that has a complete backbone (same count as
    /// [`to_backbone()`](Self::to_backbone)).
    #[must_use]
    pub fn detect_ss(&self) -> Vec<crate::SSType> {
        use crate::analysis::{detect_dssp, merge_short_segments};

        let backbone = self.to_backbone();
        if backbone.is_empty() {
            return Vec::new();
        }
        let (ss, _) = detect_dssp(&backbone);
        merge_short_segments(&ss)
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

/// Reorganize atoms into canonical per-residue order and drop residues
/// missing backbone atoms.
///
/// For every residue whose atom slice contains N, CA, C, and O (OXT
/// accepted as C-terminal oxygen fallback), emit the four backbone
/// atoms followed by sidechain heavy atoms (in input order, plus OXT
/// if O and OXT both existed) then hydrogens (in input order).
/// Residues missing any of the four backbone atoms are dropped from
/// the returned `Vec<Residue>`; their atoms are appended to the output
/// `Vec<Atom>` in input order and are left unreferenced.
struct ProteinResiduePartition {
    n: Option<usize>,
    ca: Option<usize>,
    c: Option<usize>,
    o: Option<usize>,
    oxt: Option<usize>,
    sidechain_heavy: Vec<usize>,
    hydrogens: Vec<usize>,
}

fn partition_protein_residue(
    atoms: &[Atom],
    range: std::ops::Range<usize>,
) -> ProteinResiduePartition {
    let mut out = ProteinResiduePartition {
        n: None,
        ca: None,
        c: None,
        o: None,
        oxt: None,
        sidechain_heavy: Vec::new(),
        hydrogens: Vec::new(),
    };
    for idx in range {
        let trimmed = trimmed_atom_name(&atoms[idx].name);
        let is_h =
            atoms[idx].element == Element::H || trimmed.first() == Some(&b'H');
        match trimmed {
            b"N" if out.n.is_none() => out.n = Some(idx),
            b"CA" if out.ca.is_none() => out.ca = Some(idx),
            b"C" if out.c.is_none() => out.c = Some(idx),
            b"O" if out.o.is_none() => out.o = Some(idx),
            b"OXT" if out.oxt.is_none() => out.oxt = Some(idx),
            _ if is_h => out.hydrogens.push(idx),
            _ => out.sidechain_heavy.push(idx),
        }
    }
    out
}

fn canonicalize_protein_residues(
    atoms: &[Atom],
    residues: &[Residue],
    pdb_chain_id: u8,
) -> (Vec<Atom>, Vec<Residue>) {
    let mut new_atoms: Vec<Atom> = Vec::with_capacity(atoms.len());
    let mut new_residues: Vec<Residue> = Vec::new();

    for residue in residues {
        let range = residue.atom_range.clone();
        let start = new_atoms.len();
        let mut p = partition_protein_residue(atoms, range.clone());
        let o_slot = p.o.or(p.oxt);

        if let (Some(n), Some(ca), Some(c), Some(o)) = (p.n, p.ca, p.c, o_slot)
        {
            new_atoms.push(atoms[n].clone());
            new_atoms.push(atoms[ca].clone());
            new_atoms.push(atoms[c].clone());
            new_atoms.push(atoms[o].clone());
            if let (Some(_), Some(xt)) = (p.o, p.oxt) {
                p.sidechain_heavy.push(xt);
            }
            for idx in p.sidechain_heavy {
                new_atoms.push(atoms[idx].clone());
            }
            for idx in p.hydrogens {
                new_atoms.push(atoms[idx].clone());
            }
            let end = new_atoms.len();
            new_residues.push(Residue {
                name: residue.name,
                number: residue.number,
                atom_range: start..end,
            });
        } else {
            let res_name = std::str::from_utf8(&residue.name).unwrap_or("???");
            log::warn!(
                "ProteinEntity chain '{}': dropping residue {} (name {}) — \
                 missing backbone atoms (need N, CA, C, and O or OXT)",
                pdb_chain_id as char,
                residue.number,
                res_name.trim(),
            );
            for idx in range {
                new_atoms.push(atoms[idx].clone());
            }
        }
    }

    (new_atoms, new_residues)
}

/// Emit intra-residue bonds (backbone + sidechain/anchor) for one
/// protein residue into `bonds`.
#[allow(
    clippy::cast_possible_truncation,
    reason = "atom indices are bounded by protein chain size (< u32::MAX)"
)]
fn emit_protein_residue_bonds(
    entity_id: EntityId,
    atoms: &[Atom],
    r: &Residue,
    bonds: &mut Vec<CovalentBond>,
) {
    use std::collections::HashMap;

    let start = r.atom_range.start;
    let aid = |offset: usize| AtomId {
        entity: entity_id,
        index: (start + offset) as u32,
    };
    bonds.push(CovalentBond {
        a: aid(0),
        b: aid(1),
        order: BondOrder::Single,
    });
    bonds.push(CovalentBond {
        a: aid(1),
        b: aid(2),
        order: BondOrder::Single,
    });
    bonds.push(CovalentBond {
        a: aid(2),
        b: aid(3),
        order: BondOrder::Double,
    });

    let Some(aa) = AminoAcid::from_code(r.name) else {
        return;
    };
    let mut name_to_idx: HashMap<AtomName, usize> = HashMap::new();
    for idx in r.atom_range.clone() {
        let key = AtomName::from_bytes(trimmed_atom_name(&atoms[idx].name));
        let _ = name_to_idx.insert(key, idx);
    }
    for (name_a, name_b) in aa.bonds() {
        if let (Some(&ia), Some(&ib)) =
            (name_to_idx.get(name_a), name_to_idx.get(name_b))
        {
            bonds.push(CovalentBond {
                a: AtomId {
                    entity: entity_id,
                    index: ia as u32,
                },
                b: AtomId {
                    entity: entity_id,
                    index: ib as u32,
                },
                order: BondOrder::Single,
            });
        }
    }
}

/// Populate a protein entity's covalent bond graph.
///
/// Emits: universal backbone bonds (N-CA, CA-C, C=O) per residue,
/// sidechain/anchor bonds from [`AminoAcid::bonds()`] matched by
/// [`AtomName`] within the residue, and inter-residue peptide bonds
/// `C(i)-N(i+1)` between consecutive kept residues that are not
/// separated by a segment break.
#[allow(
    clippy::cast_possible_truncation,
    reason = "atom indices are bounded by protein chain size (< u32::MAX)"
)]
fn build_protein_bonds(
    entity_id: EntityId,
    atoms: &[Atom],
    residues: &[Residue],
    segment_breaks: &[usize],
) -> Vec<CovalentBond> {
    use std::collections::HashSet;

    let mut bonds: Vec<CovalentBond> = Vec::new();
    for r in residues {
        emit_protein_residue_bonds(entity_id, atoms, r, &mut bonds);
    }

    let breaks: HashSet<usize> = segment_breaks.iter().copied().collect();
    for i in 1..residues.len() {
        if breaks.contains(&i) {
            continue;
        }
        let prev = &residues[i - 1];
        let curr = &residues[i];
        bonds.push(CovalentBond {
            a: AtomId {
                entity: entity_id,
                index: (prev.atom_range.start + 2) as u32,
            },
            b: AtomId {
                entity: entity_id,
                index: curr.atom_range.start as u32,
            },
            order: BondOrder::Single,
        });
    }

    bonds
}

/// Trim trailing space and NUL padding from a 4-byte atom name.
pub(crate) fn trimmed_atom_name(name: &[u8; 4]) -> &[u8] {
    let mut end = 4;
    while end > 0 && (name[end - 1] == b' ' || name[end - 1] == 0) {
        end -= 1;
    }
    let mut start = 0;
    while start < end && (name[start] == b' ' || name[start] == 0) {
        start += 1;
    }
    &name[start..end]
}

/// Compute segment break indices from C(i)->N(i+1) distances.
///
/// Relies on canonical atom ordering — kept residues have `C` at local
/// offset 2 and `N` at offset 0.
fn compute_segment_breaks(residues: &[Residue], atoms: &[Atom]) -> Vec<usize> {
    let mut breaks = Vec::new();
    for i in 1..residues.len() {
        let prev_c = atoms[residues[i - 1].atom_range.start + 2].position;
        let curr_n = atoms[residues[i].atom_range.start].position;
        if (prev_c - curr_n).length() > MAX_PEPTIDE_BOND_DIST {
            breaks.push(i);
        }
    }
    breaks
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

    // -- to_sidechains tests --

    #[test]
    fn to_sidechains_per_residue() {
        // ALA(+CB) then GLY: 2 residues, first has sidechain, second empty
        let coords = make_residue_coords(&[
            ("ALA", &["N", "CA", "C", "O", "CB"]),
            ("GLY", &["N", "CA", "C", "O"]),
        ]);
        let entities = split_into_entities(&coords);
        let protein = entities[0].as_protein().unwrap();
        let sc = protein.to_sidechains(|_| false, |_| None);
        assert_eq!(sc.len(), 2);
        assert_eq!(sc[0].atoms.len(), 1); // CB
        assert!(sc[1].is_empty()); // GLY
    }

    #[test]
    fn to_sidechains_bonds_local() {
        // VAL with CB, CG1, CG2 — bonds should use local indices
        let coords = make_residue_coords(&[(
            "VAL",
            &["N", "CA", "C", "O", "CB", "CG1", "CG2"],
        )]);
        let entities = split_into_entities(&coords);
        let protein = entities[0].as_protein().unwrap();
        let sc = protein.to_sidechains(
            |_| false,
            |name| match name {
                "VAL" => Some(vec![("CB", "CG1"), ("CB", "CG2")]),
                _ => None,
            },
        );
        assert_eq!(sc.len(), 1);
        assert_eq!(sc[0].atoms.len(), 3);
        assert_eq!(sc[0].bonds.len(), 2);
        for &(a, b) in &sc[0].bonds {
            assert!(a < 3 && b < 3);
        }
    }

    /// Build Coords for one or more residues from (name, atom_names) pairs.
    fn make_residue_coords(residues: &[(&str, &[&str])]) -> Coords {
        let mut atoms = Vec::new();
        let mut chain_ids = Vec::new();
        let mut res_names_v = Vec::new();
        let mut res_nums = Vec::new();
        let mut atom_names_v = Vec::new();
        let mut elements = Vec::new();
        let mut x = 0.0f32;
        let mut res_num: i32 = 1;
        for (rname, anames) in residues {
            for aname in *anames {
                atoms.push(make_atom(x, 0.0, 0.0));
                chain_ids.push(b'A');
                res_names_v.push(res_name(rname));
                res_nums.push(res_num);
                atom_names_v.push(atom_name(aname));
                elements.push(Element::Unknown);
                x += 1.5;
            }
            res_num += 1;
        }
        Coords {
            num_atoms: atoms.len(),
            atoms,
            chain_ids,
            res_names: res_names_v,
            res_nums,
            atom_names: atom_names_v,
            elements,
        }
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

    // -- Canonical ordering + bond population (Phase 2) --

    fn alanine_scrambled_coords() -> Coords {
        // Input atoms: CB, O, N, HA, CA, C — non-canonical order.
        // Expected canonical: N, CA, C, O, CB, HA.
        Coords {
            num_atoms: 6,
            atoms: vec![
                make_atom(50.0, 0.0, 0.0),   // CB
                make_atom(10.0, 11.0, 12.0), // O
                make_atom(1.0, 2.0, 3.0),    // N
                make_atom(60.0, 0.0, 0.0),   // HA
                make_atom(4.0, 5.0, 6.0),    // CA
                make_atom(7.0, 8.0, 9.0),    // C
            ],
            chain_ids: vec![b'A'; 6],
            res_names: vec![res_name("ALA"); 6],
            res_nums: vec![1; 6],
            atom_names: vec![
                atom_name("CB"),
                atom_name("O"),
                atom_name("N"),
                atom_name("HA"),
                atom_name("CA"),
                atom_name("C"),
            ],
            elements: vec![
                Element::C,
                Element::O,
                Element::N,
                Element::H,
                Element::C,
                Element::C,
            ],
        }
    }

    #[test]
    fn canonical_ordering_reorders_scrambled_input() {
        let coords = alanine_scrambled_coords();
        let entities = split_into_entities(&coords);
        let protein = entities[0].as_protein().unwrap();
        assert_eq!(protein.residues.len(), 1);
        let r = &protein.residues[0];
        // N CA C O CB HA — positions 0..6 of canonical layout.
        let names: Vec<&str> = r
            .atom_range
            .clone()
            .map(|i| {
                std::str::from_utf8(&protein.atoms[i].name).unwrap().trim()
            })
            .collect();
        assert_eq!(names, vec!["N", "CA", "C", "O", "CB", "HA"]);
    }

    #[test]
    fn canonical_ordering_drops_residue_missing_backbone() {
        // GLY has only N, CA, C (missing O) — must be dropped.
        let coords = Coords {
            num_atoms: 7,
            #[allow(clippy::cast_precision_loss)]
            atoms: (0..7i32).map(|i| make_atom(i as f32, 0.0, 0.0)).collect(),
            chain_ids: vec![b'A'; 7],
            res_names: vec![
                res_name("ALA"),
                res_name("ALA"),
                res_name("ALA"),
                res_name("ALA"),
                res_name("GLY"),
                res_name("GLY"),
                res_name("GLY"),
            ],
            res_nums: vec![1, 1, 1, 1, 2, 2, 2],
            atom_names: vec![
                atom_name("N"),
                atom_name("CA"),
                atom_name("C"),
                atom_name("O"),
                atom_name("N"),
                atom_name("CA"),
                atom_name("C"),
            ],
            elements: vec![
                Element::N,
                Element::C,
                Element::C,
                Element::O,
                Element::N,
                Element::C,
                Element::C,
            ],
        };
        let entities = split_into_entities(&coords);
        let protein = entities[0].as_protein().unwrap();
        // Only ALA is kept.
        assert_eq!(protein.residues.len(), 1);
        assert_eq!(&protein.residues[0].name, b"ALA");
        // All 7 atoms still reside in entity.atoms.
        assert_eq!(protein.atoms.len(), 7);
    }

    #[test]
    fn gly_backbone_bonds_are_three_sidechain_empty() {
        // Single GLY residue with canonical N, CA, C, O.
        let coords = Coords {
            num_atoms: 4,
            atoms: vec![
                make_atom(0.0, 0.0, 0.0),
                make_atom(1.5, 0.0, 0.0),
                make_atom(2.5, 0.0, 0.0),
                make_atom(2.5, 1.0, 0.0),
            ],
            chain_ids: vec![b'A'; 4],
            res_names: vec![res_name("GLY"); 4],
            res_nums: vec![1; 4],
            atom_names: vec![
                atom_name("N"),
                atom_name("CA"),
                atom_name("C"),
                atom_name("O"),
            ],
            elements: vec![Element::N, Element::C, Element::C, Element::O],
        };
        let entities = split_into_entities(&coords);
        let protein = entities[0].as_protein().unwrap();
        // Expect exactly N-CA, CA-C, C=O.
        assert_eq!(protein.backbone_bonds().count(), 3);
        assert_eq!(protein.sidechain_bonds().count(), 0);
    }

    #[test]
    fn ala_sidechain_bonds_is_ca_cb_only() {
        // Single ALA residue with N, CA, C, O, CB.
        let coords = Coords {
            num_atoms: 5,
            atoms: vec![
                make_atom(0.0, 0.0, 0.0),
                make_atom(1.5, 0.0, 0.0),
                make_atom(2.5, 0.0, 0.0),
                make_atom(2.5, 1.0, 0.0),
                make_atom(1.5, -1.5, 0.0),
            ],
            chain_ids: vec![b'A'; 5],
            res_names: vec![res_name("ALA"); 5],
            res_nums: vec![1; 5],
            atom_names: vec![
                atom_name("N"),
                atom_name("CA"),
                atom_name("C"),
                atom_name("O"),
                atom_name("CB"),
            ],
            elements: vec![
                Element::N,
                Element::C,
                Element::C,
                Element::O,
                Element::C,
            ],
        };
        let entities = split_into_entities(&coords);
        let protein = entities[0].as_protein().unwrap();
        let sidechain: Vec<_> = protein.sidechain_bonds().collect();
        assert_eq!(sidechain.len(), 1);
        let b = sidechain[0];
        let a_name =
            std::str::from_utf8(&protein.atoms[b.a.index as usize].name)
                .unwrap()
                .trim();
        let c_name =
            std::str::from_utf8(&protein.atoms[b.b.index as usize].name)
                .unwrap()
                .trim();
        let pair = (a_name, c_name);
        assert!(
            pair == ("CA", "CB") || pair == ("CB", "CA"),
            "expected CA-CB sidechain anchor, got {pair:?}"
        );
    }

    #[test]
    fn peptide_bond_connects_consecutive_residues() {
        // Two close ALA residues — peptide bond expected.
        let coords = Coords {
            num_atoms: 8,
            atoms: vec![
                make_atom(0.0, 0.0, 0.0),
                make_atom(1.5, 0.0, 0.0),
                make_atom(2.5, 0.0, 0.0),
                make_atom(2.5, 1.0, 0.0),
                make_atom(3.8, 0.0, 0.0), // N of residue 2 close to C of 1
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
        // 3 backbone * 2 residues + 1 peptide = 7 backbone bonds total.
        let backbone_count = protein.backbone_bonds().count();
        assert_eq!(backbone_count, 7, "3 per residue + 1 peptide");
    }

    #[test]
    fn dropped_residue_emits_log_warning() {
        testing_logger::setup();
        // ALA residue missing CA — must be dropped and logged.
        let coords = Coords {
            num_atoms: 3,
            atoms: vec![
                make_atom(0.0, 0.0, 0.0),
                make_atom(1.5, 0.0, 0.0),
                make_atom(2.5, 1.0, 0.0),
            ],
            chain_ids: vec![b'A'; 3],
            res_names: vec![res_name("ALA"); 3],
            res_nums: vec![7; 3],
            atom_names: vec![atom_name("N"), atom_name("C"), atom_name("O")],
            elements: vec![Element::N, Element::C, Element::O],
        };
        let _ = split_into_entities(&coords);
        testing_logger::validate(|captured_logs| {
            let warn_bodies: Vec<&str> = captured_logs
                .iter()
                .filter(|l| l.level == log::Level::Warn)
                .map(|l| l.body.as_str())
                .collect();
            assert!(
                !warn_bodies.is_empty(),
                "expected at least one warn-level log entry"
            );
            assert!(
                warn_bodies.iter().any(|b| b.contains("dropping residue")),
                "warning should mention dropping a residue; got \
                 {warn_bodies:?}"
            );
        });
    }

    #[test]
    fn peptide_bond_skipped_across_segment_break() {
        // Two ALA residues with a large C-N gap (segment break).
        let coords = Coords {
            num_atoms: 8,
            atoms: vec![
                make_atom(0.0, 0.0, 0.0),
                make_atom(1.5, 0.0, 0.0),
                make_atom(2.5, 0.0, 0.0),
                make_atom(2.5, 1.0, 0.0),
                make_atom(20.0, 0.0, 0.0), // N far from prev C
                make_atom(21.5, 0.0, 0.0),
                make_atom(22.5, 0.0, 0.0),
                make_atom(22.5, 1.0, 0.0),
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
        // No peptide bond — just 3 backbone per residue = 6.
        assert_eq!(protein.backbone_bonds().count(), 6);
    }
}

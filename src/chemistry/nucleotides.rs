//! Canonical nucleotides: enum + 3-letter code parsing.
//!
//! The bond table is a skeleton in Phase 1 — `bonds()` returns an
//! empty slice. The real per-base table lands in Phase 2 if a concrete
//! `NAEntity::new` consumer needs it.

use super::atom_name::AtomName;

/// Canonical nucleotide bases. The DNA/RNA distinction lives on the
/// owning [`crate::entity::molecule::MoleculeType`]; this enum names
/// only the base.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Nucleotide {
    /// Adenine.
    A,
    /// Cytosine.
    C,
    /// Guanine.
    G,
    /// Thymine (DNA).
    T,
    /// Uracil (RNA, occasionally DNA).
    U,
}

impl Nucleotide {
    /// Parse a 1-, 2-, or 3-letter PDB nucleotide code.
    ///
    /// Accepts the canonical mmCIF and PDB aliases that
    /// [`crate::entity::molecule::classify_residue`] already
    /// recognizes:
    ///
    /// - 1-letter RNA: `A`, `C`, `G`, `U`, and `T`
    /// - 2-letter DNA: `DA`, `DC`, `DG`, `DT`, `DU`
    /// - 3-letter aliases: `ADE`, `CYT`, `GUA`, `URA`, `THY`, `RAD`, `RCY`,
    ///   `RGU`
    ///
    /// Comparison is ASCII-case-insensitive. Trailing space or NUL
    /// padding is stripped before matching. Returns `None` for unknown
    /// codes.
    #[must_use]
    pub fn from_code(code: [u8; 3]) -> Option<Self> {
        let upper = [
            code[0].to_ascii_uppercase(),
            code[1].to_ascii_uppercase(),
            code[2].to_ascii_uppercase(),
        ];
        let mut end = 3usize;
        while end > 0 && (upper[end - 1] == 0 || upper[end - 1] == b' ') {
            end -= 1;
        }
        match &upper[..end] {
            b"A" | b"DA" | b"ADE" | b"RAD" => Some(Self::A),
            b"C" | b"DC" | b"CYT" | b"RCY" => Some(Self::C),
            b"G" | b"DG" | b"GUA" | b"RGU" => Some(Self::G),
            b"T" | b"DT" | b"THY" => Some(Self::T),
            b"U" | b"DU" | b"URA" => Some(Self::U),
            _ => None,
        }
    }

    /// Per-base bond table. Empty in Phase 1 — populated in Phase 2 if
    /// a concrete `NAEntity::new` consumer requires it.
    #[must_use]
    pub const fn bonds(self) -> &'static [(AtomName, AtomName)] {
        &[]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_code_canonical_one_letter() {
        assert_eq!(Nucleotide::from_code(*b"A  "), Some(Nucleotide::A));
        assert_eq!(Nucleotide::from_code(*b"C  "), Some(Nucleotide::C));
        assert_eq!(Nucleotide::from_code(*b"G  "), Some(Nucleotide::G));
        assert_eq!(Nucleotide::from_code(*b"T  "), Some(Nucleotide::T));
        assert_eq!(Nucleotide::from_code(*b"U  "), Some(Nucleotide::U));
    }

    #[test]
    fn from_code_dna_two_letter() {
        assert_eq!(Nucleotide::from_code(*b"DA\0"), Some(Nucleotide::A));
        assert_eq!(Nucleotide::from_code(*b"DC\0"), Some(Nucleotide::C));
        assert_eq!(Nucleotide::from_code(*b"DG\0"), Some(Nucleotide::G));
        assert_eq!(Nucleotide::from_code(*b"DT\0"), Some(Nucleotide::T));
    }

    #[test]
    fn from_code_three_letter_aliases() {
        assert_eq!(Nucleotide::from_code(*b"ADE"), Some(Nucleotide::A));
        assert_eq!(Nucleotide::from_code(*b"GUA"), Some(Nucleotide::G));
        assert_eq!(Nucleotide::from_code(*b"THY"), Some(Nucleotide::T));
    }

    #[test]
    fn from_code_unknown_returns_none() {
        assert_eq!(Nucleotide::from_code(*b"XXX"), None);
        assert_eq!(Nucleotide::from_code(*b"DI\0"), None);
    }

    #[test]
    fn from_code_handles_padding() {
        assert_eq!(Nucleotide::from_code(*b"A  "), Some(Nucleotide::A));
    }

    #[test]
    fn bonds_skeleton_is_empty() {
        for nt in [
            Nucleotide::A,
            Nucleotide::C,
            Nucleotide::G,
            Nucleotide::T,
            Nucleotide::U,
        ] {
            assert!(nt.bonds().is_empty());
        }
    }
}

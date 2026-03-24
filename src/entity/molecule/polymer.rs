//! Polymer residue type.

use std::ops::Range;

/// A single residue within a polymer entity.
#[derive(Debug, Clone)]
pub struct Residue {
    /// 3-character residue name (e.g. b"ALA").
    pub name: [u8; 3],
    /// Residue sequence number.
    pub number: i32,
    /// Index range into the parent entity's atom list.
    pub atom_range: Range<usize>,
}

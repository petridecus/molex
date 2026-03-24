//! Wire format types: Coords, CoordsAtom, CoordsError, ChainIdMapper.

use thiserror::Error;

use crate::element::Element;

/// Errors that can occur during COORDS operations.
#[derive(Error, Debug)]
pub enum CoordsError {
    /// The binary data does not conform to the expected COORDS/ASSEM format.
    #[error("Invalid COORDS format: {0}")]
    InvalidFormat(String),
    /// A PDB file could not be parsed.
    #[error("Failed to parse PDB: {0}")]
    PdbParseError(String),
    /// An error occurred during binary serialization or deserialization.
    #[error("Serialization error: {0}")]
    SerializationError(String),
}

/// Maps multi-character chain ID strings to unique `u8` values.
///
/// Structures with >26 chains (ribosomes, virus capsids) use multi-character
/// chain IDs in mmCIF format (e.g., "AA", "AB"). Since `Coords.chain_ids`
/// stores a single `u8` per atom, this mapper assigns a unique byte to each
/// distinct chain string, preventing collisions that cause cross-chain
/// rendering artifacts.
///
/// Assigns printable ASCII characters (A-Z, a-z, 0-9, then other printable
/// chars) so that PDB export produces valid chain ID columns.
pub struct ChainIdMapper {
    map: std::collections::HashMap<String, u8>,
    next_idx: usize,
}

/// Printable chain ID characters in conventional order.
const CHAIN_CHARS: &[u8] =
    b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789!#$%&()*+,-./:;<=>?@[]^_`{|}~";

impl Default for ChainIdMapper {
    fn default() -> Self {
        Self::new()
    }
}

impl ChainIdMapper {
    /// Create a new empty mapper with no chain assignments.
    #[must_use]
    pub fn new() -> Self {
        Self {
            map: std::collections::HashMap::new(),
            next_idx: 0,
        }
    }

    /// Get or assign a unique `u8` for the given chain ID string.
    pub fn get_or_assign(&mut self, chain_id: &str) -> u8 {
        if let Some(&id) = self.map.get(chain_id) {
            return id;
        }
        let byte = if self.next_idx < CHAIN_CHARS.len() {
            CHAIN_CHARS[self.next_idx]
        } else {
            #[allow(clippy::cast_possible_truncation)]
            {
                (self.next_idx - CHAIN_CHARS.len()) as u8
            }
        };
        self.next_idx += 1;
        let _ = self.map.insert(chain_id.to_owned(), byte);
        byte
    }
}

/// Single atom with coordinates and crystallographic factors.
#[derive(Debug, Clone)]
pub struct CoordsAtom {
    /// X coordinate in angstroms.
    pub x: f32,
    /// Y coordinate in angstroms.
    pub y: f32,
    /// Z coordinate in angstroms.
    pub z: f32,
    /// Crystallographic occupancy (0.0 to 1.0).
    pub occupancy: f32,
    /// Temperature factor (B-factor) in square angstroms.
    pub b_factor: f32,
}

/// Complete coordinate structure with atom metadata.
#[derive(Debug, Clone)]
pub struct Coords {
    /// Total number of atoms in this structure.
    pub num_atoms: usize,
    /// Per-atom position and crystallographic data.
    pub atoms: Vec<CoordsAtom>,
    /// Per-atom chain identifier byte (e.g., b'A').
    pub chain_ids: Vec<u8>,
    /// Per-atom 3-character residue name (e.g., b"ALA").
    pub res_names: Vec<[u8; 3]>,
    /// Per-atom residue sequence number.
    pub res_nums: Vec<i32>,
    /// Per-atom 4-character PDB atom name (e.g., b" CA ").
    pub atom_names: Vec<[u8; 4]>,
    /// Chemical element per atom.
    pub elements: Vec<Element>,
}

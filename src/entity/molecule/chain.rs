//! `ChainIdMapper`: shared helper used by [`EntityBuilder`] and the
//! AtomWorks bridge to map multi-character chain ID strings down to a
//! single printable byte that the entity types and PDB writer expect.
//!
//! Structures with >26 chains (ribosomes, virus capsids) use multi-character
//! chain IDs in mmCIF format (e.g., "AA", "AB"). Assigning a printable
//! ASCII byte (A-Z, a-z, 0-9, then other printable characters) per
//! distinct chain string prevents collisions that would otherwise show
//! up as cross-chain rendering artifacts.
//!
//! [`EntityBuilder`]: super::builder::EntityBuilder

use std::collections::HashMap;

/// Maps multi-character chain ID strings to unique `u8` values.
pub(crate) struct ChainIdMapper {
    map: HashMap<String, u8>,
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
            map: HashMap::new(),
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

//! Adapter-internal `Coords` parser intermediate and entity-utility
//! helpers.
//!
//! The `Coords` struct is a flat parallel-array atom record used by the
//! PDB/CIF/BCIF parsers and the ASSEM01 wire decoder as an intermediate
//! step on the way to / from `Vec<MoleculeEntity>`. It is not part of the
//! public API (see `docs/COORDS_RETIREMENT_PLAN.md`).
//!
//! The ASSEM01 wire format itself lives in `crate::ops::wire`.

mod assembly;
mod bridge;
mod types;

#[allow(
    unused_imports,
    reason = "test-only caller post-Phase 4a; function retained pending Phase \
              3/6 cleanup"
)]
pub(crate) use assembly::update_protein_entities;
pub use assembly::{ca_positions, residue_count};
pub(crate) use bridge::{
    coords_to_molecule_entity, extract_atom_set_and_residues,
    split_into_entities,
};
pub use types::AdapterError;
pub(crate) use types::{ChainIdMapper, Coords, CoordsAtom};

#[cfg(test)]
#[path = "tests.rs"]
mod tests;

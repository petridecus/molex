//! Entity system: molecules and surfaces.

pub mod molecule;
pub mod surface;

// Re-export commonly used molecule types at the entity level.
pub use molecule::atom::Atom;
pub use molecule::id::{EntityId, EntityIdAllocator};
#[allow(
    deprecated,
    reason = "legacy re-export; removed alongside ProteinResidue in Phase 5"
)]
pub use molecule::protein::{ProteinResidue, ResidueBackbone, Sidechain};
pub use molecule::{
    classify_residue, MoleculeEntity, MoleculeType, NucleotideRing,
};
// Re-export surface types.
pub use surface::{Density, DensityError, VoxelGrid};

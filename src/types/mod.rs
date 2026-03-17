//! Core data types for the molex crate.

pub mod assembly;
pub mod coords;
pub mod density;
pub mod element;
pub mod entity;

// Re-export commonly used items
pub use assembly::{
    assembly_bytes,
    ca_positions,
    // Free functions (replace Assembly struct methods)
    protein_coords,
    residue_count,
    update_protein_entities,
};
pub use coords::{
    atom_count,
    // Binary serialization
    deserialize,
    deserialize_assembly,
    serialize,
    serialize_assembly,
    AtomMetadata,
    ChainIdMapper,
    Coords,
    CoordsAtom,
    CoordsError,
    Element,
    ResidueAtoms,
    ValidationResult,
    ASSEMBLY_MAGIC,
};
pub use density::{DensityError, DensityMap};
pub use entity::{
    classify_residue, coords_to_entity_kind, extract_by_type, merge_entities,
    split_into_entities, AtomSet, EntityKind, MoleculeEntity, MoleculeType,
    PolymerChain, PolymerData, Residue,
};

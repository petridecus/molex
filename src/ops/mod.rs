//! Operations on molecular data.

/// Wire formats (COORDS01, ASSEM01), their types, and assembly helpers.
pub mod codec;
pub mod transform;

pub use transform::{
    align_coords_bytes, align_to_reference, centroid,
    extract_backbone_segments, extract_ca_from_chains, extract_ca_positions,
    kabsch_alignment, kabsch_alignment_with_scale, transform_entities,
    transform_entities_with_scale,
};

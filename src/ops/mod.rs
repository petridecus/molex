//! Operations on molecular data.

/// Adapter / wire error type and protein-CA helper.
pub mod codec;
pub mod transform;
/// ASSEM01 binary wire format encoder/decoder.
pub mod wire;

pub use transform::{
    align_to_reference, centroid, extract_backbone_segments,
    extract_ca_from_chains, extract_ca_positions, kabsch_alignment,
    kabsch_alignment_with_scale, transform_entities,
    transform_entities_with_scale,
};

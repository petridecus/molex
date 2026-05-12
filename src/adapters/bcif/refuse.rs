//! Refusal helpers for the BinaryCIF adapter.

use crate::ops::codec::AdapterError;

/// Refusal for files containing multiple `dataBlocks` entries.
///
/// Production coordinate files are single-block; multi-block inputs come from
/// chemical-component dictionaries which the coordinate adapter intentionally
/// does not handle.
pub(super) fn multi_block_error(count: usize) -> AdapterError {
    AdapterError::InvalidFormat(format!(
        "BinaryCIF input has {count} data blocks; the coordinate adapter \
         consumes single-block files only. Use a chemical-component \
         dictionary adapter for components.bcif-style inputs."
    ))
}

/// Refusal for inputs whose chain count exhausts the printable-byte mapper.
pub(super) fn too_many_chains_error(limit: usize) -> AdapterError {
    AdapterError::InvalidFormat(format!(
        "BinaryCIF input declares more than {limit} distinct chains; molex's \
         chain-byte mapper is exhausted. Re-ingest via the mmCIF route on a \
         smaller subset."
    ))
}

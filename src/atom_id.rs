//! Stable cross-entity atom identifier.
//!
//! [`AtomId`] is the cross-cutting identity used by bond endpoints,
//! cross-entity references, and any other API that names an atom from
//! outside its owning entity.

use crate::entity::EntityId;

/// Identifier for an atom owned by a specific entity.
///
/// Layout: `EntityId` (4-byte newtype around `u32`) + `u32` index, for
/// 8 bytes total.
///
/// # Stability
///
/// An [`AtomId`] remains valid across `Assembly` mutations that do not
/// remove the atom's owning entity. Mutations that reorder or remove
/// atoms within an entity, or that remove the entity itself, invalidate
/// existing ids; callers must drop them. There is no automatic
/// invalidation mechanism.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct AtomId {
    /// Owning entity.
    pub entity: EntityId,
    /// Index into the owning entity's `atoms` slice.
    pub index: u32,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::entity::EntityIdAllocator;

    #[test]
    fn layout_is_eight_bytes() {
        assert_eq!(size_of::<AtomId>(), 8);
        assert_eq!(align_of::<AtomId>(), 4);
    }

    #[test]
    fn round_trip_through_components() {
        let mut alloc = EntityIdAllocator::new();
        let entity = alloc.allocate();
        let id = AtomId { entity, index: 42 };
        assert_eq!(id.entity, entity);
        assert_eq!(id.index, 42);
    }

    #[test]
    fn copy_and_eq() {
        let mut alloc = EntityIdAllocator::new();
        let e = alloc.allocate();
        let a = AtomId {
            entity: e,
            index: 7,
        };
        let b = a;
        assert_eq!(a, b);
    }
}

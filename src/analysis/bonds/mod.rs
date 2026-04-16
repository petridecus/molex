//! Bond detection: covalent, hydrogen, and disulfide.

pub mod covalent;
pub mod disulfide;
pub mod hydrogen;

pub use covalent::{infer_bonds, BondOrder, InferredBond, DEFAULT_TOLERANCE};
pub use disulfide::detect_disulfides;
#[allow(
    deprecated,
    reason = "legacy re-export; removed in Phase 5 of assembly migration"
)]
pub use disulfide::{detect_disulfide_bonds, DisulfideBond};
#[allow(
    deprecated,
    reason = "legacy re-export; detect_hbonds removed in Phase 5, HBond stays"
)]
pub use hydrogen::detect_hbonds;
pub use hydrogen::HBond;

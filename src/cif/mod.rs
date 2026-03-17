//! CIF/STAR parser and typed extractors.
//!
//! Two-layer design:
//! - **Layer 1 (DOM)**: [`mod@parse`] any CIF/STAR file into an untyped
//!   [`Document`] tree.
//! - **Layer 2 (Extractors)**: Pull typed data out via `TryFrom<&Block>` —
//!   [`CoordinateData`], [`ReflectionData`], or auto-detect with
//!   [`CifContent`].
//!
//! ```ignore
//! let doc = molex::cif::parse(input)?;
//! let block = &doc.blocks[0];
//!
//! // Caller knows what they have:
//! let coords = CoordinateData::try_from(block)?;
//!
//! // Or auto-detect:
//! let content = CifContent::from(block.clone());
//! ```

pub mod dom;
pub mod extract;
pub mod parse;

// DOM types
pub use dom::{Block, ColumnIter, Columns, Document, Loop, RowIter, Value};
// Typed extractors
pub use extract::{
    AtomSite, CifContent, CoordinateData, ExtractionError, Reflection,
    ReflectionData, UnitCell,
};
// Parser
pub use parse::{parse, CifParseError};

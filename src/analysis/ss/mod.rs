//! Secondary structure detection and parsing.

pub mod dssp;
pub mod string;

pub use dssp::classify;
pub use string::from_string;

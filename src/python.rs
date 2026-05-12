//! Python bindings for parsing structure files into ASSEM01 bytes and
//! emitting PDB from ASSEM01 bytes. AtomArray interop lives in
//! `adapters::atomworks`.

use pyo3::prelude::*;

use crate::adapters::{cif, pdb};
use crate::ops::wire::{
    assembly_bytes, deserialize_assembly, serialize_assembly,
};

/// Parse a PDB string and emit ASSEM01 binary bytes.
///
/// # Errors
///
/// Returns a `PyValueError` if the PDB string cannot be parsed.
#[pyfunction]
#[allow(clippy::needless_pass_by_value)]
pub fn pdb_to_assembly_bytes(pdb_str: String) -> PyResult<Vec<u8>> {
    let entities = pdb::pdb_str_to_entities(&pdb_str).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string())
    })?;
    assembly_bytes(&entities).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string())
    })
}

/// Parse an mmCIF string and emit ASSEM01 binary bytes.
///
/// # Errors
///
/// Returns a `PyValueError` if the mmCIF string cannot be parsed.
#[pyfunction]
#[allow(clippy::needless_pass_by_value)]
pub fn mmcif_to_assembly_bytes(cif_str: String) -> PyResult<Vec<u8>> {
    let entities = cif::mmcif_str_to_entities(&cif_str).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string())
    })?;
    assembly_bytes(&entities).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string())
    })
}

/// Decode ASSEM01 bytes and emit a PDB-format string.
///
/// # Errors
///
/// Returns a `PyValueError` if the bytes cannot be deserialized.
#[pyfunction]
#[allow(clippy::needless_pass_by_value)]
pub fn assembly_bytes_to_pdb(bytes: Vec<u8>) -> PyResult<String> {
    let assembly = deserialize_assembly(&bytes).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string())
    })?;
    pdb::assembly_to_pdb(&assembly).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string())
    })
}

/// Round-trip ASSEM01 bytes through `Assembly` and back (validation).
///
/// # Errors
///
/// Returns a `PyValueError` if the bytes cannot be deserialized or
/// re-serialized.
#[pyfunction]
#[allow(clippy::needless_pass_by_value)]
pub fn deserialize_assembly_bytes(bytes: Vec<u8>) -> PyResult<Vec<u8>> {
    deserialize_assembly(&bytes)
        .and_then(|assembly| serialize_assembly(&assembly))
        .map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string())
        })
}

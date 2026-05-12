//! Python bindings for parsing structure files into ASSEM01 bytes and
//! emitting PDB from ASSEM01 bytes. AtomArray interop lives in
//! `adapters::atomworks`.
//!
//! Handle-based surface (`Assembly`, `EditList`, `Variant`, `AtomRow`,
//! `ProtonationState`) parallels the C API in `c_api::edit` per
//! `feedback_molex_api_parity.md`.

use glam::Vec3;
use pyo3::prelude::*;

use crate::adapters::{cif, pdb};
use crate::assembly::Assembly;
use crate::chemistry::variant::{ProtonationState, VariantTag};
use crate::element::Element;
use crate::entity::molecule::atom::Atom as MolexAtom;
use crate::entity::molecule::id::EntityIdAllocator;
use crate::ops::edit::AssemblyEdit;
use crate::ops::wire::delta::{
    deserialize_edits, serialize_edits, DeltaSerializeError,
};
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

// ---------------------------------------------------------------------------
// Pyclass: Assembly
// ---------------------------------------------------------------------------

fn value_err<E: std::fmt::Display>(e: E) -> PyErr {
    PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string())
}

/// Python handle wrapping a `molex::Assembly`. Parallels
/// `molex_Assembly` in the C API and `molex::Assembly` in C++.
#[pyclass(name = "Assembly", module = "molex")]
pub struct PyAssembly {
    inner: Assembly,
}

#[pymethods]
impl PyAssembly {
    /// Decode ASSEM01 / ASSEM02 bytes into an `Assembly`.
    ///
    /// # Errors
    ///
    /// `PyValueError` if the magic header is wrong or the payload is
    /// truncated / malformed.
    #[staticmethod]
    #[pyo3(name = "from_assem01")]
    #[allow(clippy::needless_pass_by_value)]
    pub fn from_assem01(bytes: Vec<u8>) -> PyResult<Self> {
        let inner = deserialize_assembly(&bytes).map_err(value_err)?;
        Ok(Self { inner })
    }

    /// Emit the assembly as ASSEM02 bytes (latest wire format).
    ///
    /// # Errors
    ///
    /// `PyValueError` propagated from the serializer (currently
    /// infallible in practice; the `Result` is kept for API stability).
    pub fn to_assem01(&self) -> PyResult<Vec<u8>> {
        serialize_assembly(&self.inner).map_err(value_err)
    }

    /// Monotonic generation counter.
    #[must_use]
    #[getter]
    pub fn generation(&self) -> u64 {
        self.inner.generation()
    }

    /// Apply every edit in `edits` to the assembly in order.
    ///
    /// # Errors
    ///
    /// `PyValueError` when an edit references an unknown entity /
    /// residue, when a coord-count doesn't match, or when an edit
    /// targets a non-polymer entity. Edits applied before the failing
    /// one stay applied; the generation counter advances per
    /// successful edit.
    pub fn apply_edits(&mut self, edits: &PyEditList) -> PyResult<()> {
        self.inner.apply_edits(&edits.inner).map_err(value_err)
    }

    /// Decode DELTA01 bytes and apply them in one call.
    ///
    /// # Errors
    ///
    /// `PyValueError` for invalid DELTA01 bytes (bad magic / truncated
    /// payload / unknown tag) or any of the apply failures documented
    /// on `apply_edits`.
    #[allow(clippy::needless_pass_by_value)]
    pub fn apply_delta01(&mut self, bytes: Vec<u8>) -> PyResult<()> {
        let edits = deserialize_edits(&bytes).map_err(value_err)?;
        self.inner.apply_edits(&edits).map_err(value_err)
    }
}

// ---------------------------------------------------------------------------
// Pyclass: EditList
// ---------------------------------------------------------------------------

/// Python handle for an ordered list of typed Assembly edits.
/// Parallels `molex_EditList` in C and `molex::EditList` in C++.
#[pyclass(name = "EditList", module = "molex")]
#[derive(Default)]
pub struct PyEditList {
    inner: Vec<AssemblyEdit>,
}

fn make_entity_id(raw: u32) -> crate::entity::molecule::id::EntityId {
    EntityIdAllocator::new().from_raw(raw)
}

fn coords_to_vec3(coords: Vec<(f32, f32, f32)>) -> Vec<Vec3> {
    coords
        .into_iter()
        .map(|(x, y, z)| Vec3::new(x, y, z))
        .collect()
}

#[pymethods]
impl PyEditList {
    /// Construct an empty edit list.
    #[must_use]
    #[new]
    pub fn new() -> Self {
        Self::default()
    }

    /// Number of edits currently in the list.
    #[must_use]
    pub fn __len__(&self) -> usize {
        self.inner.len()
    }

    /// Append a `SetEntityCoords` edit. `coords` is a sequence of
    /// `(x, y, z)` tuples.
    #[allow(clippy::needless_pass_by_value)]
    pub fn push_set_entity_coords(
        &mut self,
        entity_id: u32,
        coords: Vec<(f32, f32, f32)>,
    ) {
        self.inner.push(AssemblyEdit::SetEntityCoords {
            entity: make_entity_id(entity_id),
            coords: coords_to_vec3(coords),
        });
    }

    /// Append a `SetResidueCoords` edit.
    #[allow(clippy::needless_pass_by_value)]
    pub fn push_set_residue_coords(
        &mut self,
        entity_id: u32,
        residue_idx: usize,
        coords: Vec<(f32, f32, f32)>,
    ) {
        self.inner.push(AssemblyEdit::SetResidueCoords {
            entity: make_entity_id(entity_id),
            residue_idx,
            coords: coords_to_vec3(coords),
        });
    }

    /// Append a `MutateResidue` edit. `new_name` is a 3-byte residue
    /// name (e.g. `b"HIS"`); `atoms` is a list of `AtomRow`;
    /// `variants` is a list of `Variant`.
    #[allow(
        clippy::needless_pass_by_value,
        clippy::too_many_arguments,
        reason = "pyo3 method signature mirrors AssemblyEdit::MutateResidue \
                  variant's fields; Python callers use keyword arguments so \
                  the flat shape is the most ergonomic surface."
    )]
    pub fn push_mutate_residue(
        &mut self,
        entity_id: u32,
        residue_idx: usize,
        new_name: [u8; 3],
        atoms: Vec<PyAtomRow>,
        variants: Vec<PyVariant>,
    ) {
        let new_atoms = atoms.into_iter().map(PyAtomRow::into_atom).collect();
        let new_variants =
            variants.into_iter().map(|v| v.inner).collect::<Vec<_>>();
        self.inner.push(AssemblyEdit::MutateResidue {
            entity: make_entity_id(entity_id),
            residue_idx,
            new_name,
            new_atoms,
            new_variants,
        });
    }

    /// Append a `SetVariants` edit.
    #[allow(clippy::needless_pass_by_value)]
    pub fn push_set_variants(
        &mut self,
        entity_id: u32,
        residue_idx: usize,
        variants: Vec<PyVariant>,
    ) {
        let variants =
            variants.into_iter().map(|v| v.inner).collect::<Vec<_>>();
        self.inner.push(AssemblyEdit::SetVariants {
            entity: make_entity_id(entity_id),
            residue_idx,
            variants,
        });
    }

    /// Serialize this list as DELTA01 bytes.
    ///
    /// # Errors
    ///
    /// `PyValueError` when the list contains a topology edit
    /// (`AddEntity` / `RemoveEntity`) which is not representable in
    /// DELTA01; callers should broadcast a full ASSEM02 snapshot in
    /// that case.
    pub fn to_delta01(&self) -> PyResult<Vec<u8>> {
        serialize_edits(&self.inner)
            .map_err(|e: DeltaSerializeError| value_err(format!("{e}")))
    }

    /// Decode DELTA01 bytes into a new edit list.
    ///
    /// # Errors
    ///
    /// `PyValueError` if the magic header is wrong, the payload is
    /// truncated, or an unknown edit tag is encountered.
    #[staticmethod]
    #[allow(clippy::needless_pass_by_value)]
    pub fn from_delta01(bytes: Vec<u8>) -> PyResult<Self> {
        let inner = deserialize_edits(&bytes).map_err(value_err)?;
        Ok(Self { inner })
    }
}

// ---------------------------------------------------------------------------
// Pyclass: Variant
// ---------------------------------------------------------------------------

/// Python handle wrapping a single `VariantTag`.
#[pyclass(name = "Variant", module = "molex")]
#[derive(Clone)]
pub struct PyVariant {
    inner: VariantTag,
}

#[pymethods]
impl PyVariant {
    /// Chain N-terminus patch.
    #[must_use]
    #[staticmethod]
    pub fn n_terminus() -> Self {
        Self {
            inner: VariantTag::NTerminus,
        }
    }

    /// Chain C-terminus patch.
    #[must_use]
    #[staticmethod]
    pub fn c_terminus() -> Self {
        Self {
            inner: VariantTag::CTerminus,
        }
    }

    /// Participates in a disulfide bond.
    #[must_use]
    #[staticmethod]
    pub fn disulfide() -> Self {
        Self {
            inner: VariantTag::Disulfide,
        }
    }

    /// `HID` — delta-protonated histidine.
    #[must_use]
    #[staticmethod]
    pub fn protonation_his_delta() -> Self {
        Self {
            inner: VariantTag::Protonation(ProtonationState::HisDelta),
        }
    }

    /// `HIE` — epsilon-protonated histidine.
    #[must_use]
    #[staticmethod]
    pub fn protonation_his_epsilon() -> Self {
        Self {
            inner: VariantTag::Protonation(ProtonationState::HisEpsilon),
        }
    }

    /// `HIP` — doubly-protonated histidine.
    #[must_use]
    #[staticmethod]
    pub fn protonation_his_doubly() -> Self {
        Self {
            inner: VariantTag::Protonation(ProtonationState::HisDoubly),
        }
    }

    /// Custom protonation state, e.g. `"ASP_PROTONATED"`.
    #[must_use]
    #[staticmethod]
    pub fn protonation_custom(name: String) -> Self {
        Self {
            inner: VariantTag::Protonation(ProtonationState::Custom(name)),
        }
    }

    /// Open-ended variant tag carrying an arbitrary string.
    #[must_use]
    #[staticmethod]
    pub fn other(name: String) -> Self {
        Self {
            inner: VariantTag::Other(name),
        }
    }
}

// ---------------------------------------------------------------------------
// Pyclass: AtomRow
// ---------------------------------------------------------------------------

/// Python handle wrapping a single atom row (position + name +
/// element). Parallels `molex_AtomRow` in C.
#[pyclass(name = "AtomRow", module = "molex")]
#[derive(Clone)]
pub struct PyAtomRow {
    /// X coordinate.
    #[pyo3(get, set)]
    pub position_x: f32,
    /// Y coordinate.
    #[pyo3(get, set)]
    pub position_y: f32,
    /// Z coordinate.
    #[pyo3(get, set)]
    pub position_z: f32,
    /// PDB atom name (4 bytes, space-padded; e.g. `b"CA  "`).
    #[pyo3(get, set)]
    pub name: [u8; 4],
    /// Element symbol (e.g. `"C"`, `"Fe"`).
    #[pyo3(get, set)]
    pub element: String,
}

#[pymethods]
impl PyAtomRow {
    /// Construct an `AtomRow` from positional fields.
    #[must_use]
    #[new]
    pub fn new(
        position_x: f32,
        position_y: f32,
        position_z: f32,
        name: [u8; 4],
        element: String,
    ) -> Self {
        Self {
            position_x,
            position_y,
            position_z,
            name,
            element,
        }
    }
}

impl PyAtomRow {
    fn into_atom(self) -> MolexAtom {
        let element = Element::from_symbol(self.element.trim());
        MolexAtom {
            position: Vec3::new(
                self.position_x,
                self.position_y,
                self.position_z,
            ),
            occupancy: 1.0,
            b_factor: 0.0,
            element,
            name: self.name,
            formal_charge: 0,
        }
    }
}

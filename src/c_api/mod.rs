//! C ABI surface for molex.
//!
//! Gated behind the `c-api` cargo feature. Bridges `Assembly` and its
//! associated read accessors plus PDB/CIF/BCIF read + PDB write across the
//! Rust/C FFI boundary via opaque-pointer handles and heap-allocated
//! output buffers.
//!
//! # Naming
//!
//! Types exposed here use bare names (`Assembly`, `Entity`, `Residue`,
//! `Atom`, `molex_MoleculeType`). C++ consumers wrap them in `namespace molex`
//! for disambiguation; pure C consumers accept the global-namespace
//! placement. Free functions keep a `molex_` prefix because C has no
//! namespace mechanism and the prefix is part of the function name
//! itself rather than a per-layer transformation.
//!
//! # Ownership model
//!
//! - `*mut molex_Assembly` returned by parser entry points is owned by the
//!   caller and must be released with [`molex_assembly_free`].
//! - `*const molex_Entity`, `*const molex_Residue`, `*const molex_Atom` handles
//!   are non-owning views into the parent `Assembly`. They are valid for the
//!   lifetime of the parent and must not outlive it.
//! - Byte buffers returned via `out_buf` / `out_len` (e.g. from
//!   [`molex_assembly_to_pdb`]) are heap-allocated by Rust and must be freed
//!   with [`molex_free_bytes`].
//! - `const char*` + length pairs passed *into* parsers are borrowed for the
//!   duration of the call only.
//!
//! # Error reporting
//!
//! Fallible entry points either return a null handle (parsers) or a
//! nonzero status code (writers). On failure, a human-readable message
//! is recorded in thread-local storage and can be retrieved with
//! [`molex_last_error_message`]. Successful calls do not clear the
//! error state, so consumers should only read it immediately after a
//! failed call.

#![allow(
    clippy::missing_safety_doc,
    reason = "FFI surface: every fn that takes raw pointers documents its \
              preconditions in prose; safety contract is on the caller side \
              rather than each unsafe block."
)]
#![allow(
    clippy::not_unsafe_ptr_arg_deref,
    reason = "FFI entry points receive raw pointers from C/C++ callers; \
              dereferencing them after a null check is the whole point."
)]
#![allow(
    non_camel_case_types,
    reason = "FFI types intentionally match their C-side spelling \
              (`molex_Assembly` etc.) so the same identifier appears in Rust \
              source, cbindgen output, and C consumer code."
)]

use std::cell::RefCell;
use std::ffi::c_char;
use std::sync::Arc;

use crate::adapters::{bcif, cif, pdb};
use crate::assembly::Assembly;
use crate::element::Element;
use crate::entity::molecule::{Atom, MoleculeEntity, MoleculeType, Residue};

thread_local! {
    static LAST_ERROR: RefCell<Option<Vec<u8>>> = const { RefCell::new(None) };
}

fn set_last_error<E: std::fmt::Display>(err: &E) {
    let msg = format!("{err}").into_bytes();
    LAST_ERROR.with(|slot| {
        *slot.borrow_mut() = Some(msg);
    });
}

fn clear_last_error() {
    LAST_ERROR.with(|slot| {
        *slot.borrow_mut() = None;
    });
}

// ---------------------------------------------------------------------------
// Opaque handle types
// ---------------------------------------------------------------------------
//
// These are zero-sized tag types. Pointers in the FFI are tagged with
// these but always cast back to the real inner type before dereference.
// Defining them as 0-sized with no fields stops cbindgen from drilling
// into the inner Rust types and leaking molex internals into `molex.h`.

/// Owned, top-level assembly handle. Free with [`molex_assembly_free`].
pub struct molex_Assembly;

/// Non-owning view of a single entity within an assembly.
pub struct molex_Entity;

/// Non-owning view of a single polymer residue within an entity.
pub struct molex_Residue;

/// Non-owning view of a single atom within an entity.
pub struct molex_Atom;

/// Discriminant for an entity's molecule classification across the FFI
/// boundary. Stable integer codes so C consumers can pattern-match
/// without depending on Rust's enum layout.
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum molex_MoleculeType {
    /// Amino acid polymer.
    Protein = 0,
    /// Deoxyribonucleic acid polymer.
    Dna = 1,
    /// Ribonucleic acid polymer.
    Rna = 2,
    /// Non-polymer small molecule (drug, substrate).
    Ligand = 3,
    /// Single-atom metal or halide ion.
    Ion = 4,
    /// Water molecule.
    Water = 5,
    /// Lipid or detergent molecule.
    Lipid = 6,
    /// Enzyme cofactor (heme, NAD, FAD, Fe-S cluster).
    Cofactor = 7,
    /// Crystallization solvent or buffer artifact.
    Solvent = 8,
}

impl From<MoleculeType> for molex_MoleculeType {
    fn from(value: MoleculeType) -> Self {
        match value {
            MoleculeType::Protein => molex_MoleculeType::Protein,
            MoleculeType::DNA => molex_MoleculeType::Dna,
            MoleculeType::RNA => molex_MoleculeType::Rna,
            MoleculeType::Ligand => molex_MoleculeType::Ligand,
            MoleculeType::Ion => molex_MoleculeType::Ion,
            MoleculeType::Water => molex_MoleculeType::Water,
            MoleculeType::Lipid => molex_MoleculeType::Lipid,
            MoleculeType::Cofactor => molex_MoleculeType::Cofactor,
            MoleculeType::Solvent => molex_MoleculeType::Solvent,
        }
    }
}

// ---------------------------------------------------------------------------
// Status codes
// ---------------------------------------------------------------------------

/// Success status returned by writer-style entry points.
pub const MOLEX_OK: i32 = 0;
/// Generic failure status returned by writer-style entry points.
pub const MOLEX_ERR: i32 = -1;
/// One of the input pointers was null.
pub const MOLEX_ERR_NULL: i32 = -2;

// ---------------------------------------------------------------------------
// Error reporting
// ---------------------------------------------------------------------------

/// Pointer to the most recent error message recorded on this thread.
///
/// `out_len` receives the length in bytes (excluding any null terminator;
/// the returned pointer is *not* guaranteed to be null-terminated). When
/// no error has been recorded, returns null and writes 0 to `out_len`.
///
/// The returned pointer is valid until the next FFI call that records or
/// clears the thread-local error state on the same thread. Copy it
/// immediately if the caller needs to retain the message.
#[no_mangle]
pub extern "C" fn molex_last_error_message(
    out_len: *mut usize,
) -> *const c_char {
    LAST_ERROR.with(|slot| {
        let borrow = slot.borrow();
        match borrow.as_ref() {
            None => {
                if !out_len.is_null() {
                    unsafe {
                        *out_len = 0;
                    }
                }
                std::ptr::null()
            }
            Some(msg) => {
                if !out_len.is_null() {
                    unsafe {
                        *out_len = msg.len();
                    }
                }
                msg.as_ptr().cast::<c_char>()
            }
        }
    })
}

// ---------------------------------------------------------------------------
// Buffer / handle lifecycle
// ---------------------------------------------------------------------------

/// Free a buffer previously returned via `out_buf` / `out_len` from a
/// writer-style entry point.
///
/// Safe to call with a null pointer (no-op). `len` must match the value
/// originally written to `out_len`.
#[no_mangle]
pub extern "C" fn molex_free_bytes(bytes: *mut u8, len: usize) {
    if bytes.is_null() {
        return;
    }
    drop(unsafe {
        Box::from_raw(std::ptr::slice_from_raw_parts_mut(bytes, len))
    });
}

/// Free an assembly handle returned by a parser entry point.
///
/// Safe to call with a null pointer (no-op). Borrowed sub-handles
/// (`*const molex_Entity`, `*const molex_Residue`, `*const molex_Atom`)
/// returned from the assembly are invalidated and must not be used after this
/// call.
#[no_mangle]
pub extern "C" fn molex_assembly_free(assembly: *mut molex_Assembly) {
    if assembly.is_null() {
        return;
    }
    drop(unsafe { Box::from_raw(assembly.cast::<Assembly>()) });
}

// ---------------------------------------------------------------------------
// Parser entry points (PDB / CIF / BCIF -> Assembly)
// ---------------------------------------------------------------------------

fn assembly_from_entities(
    entities: Vec<MoleculeEntity>,
) -> *mut molex_Assembly {
    let inner = Assembly::new(entities);
    Box::into_raw(Box::new(inner)).cast::<molex_Assembly>()
}

fn assembly_inner<'a>(assembly: *const molex_Assembly) -> Option<&'a Assembly> {
    if assembly.is_null() {
        return None;
    }
    Some(unsafe { &*assembly.cast::<Assembly>() })
}

fn slice_from_raw<'a, T>(ptr: *const T, len: usize) -> Option<&'a [T]> {
    if ptr.is_null() {
        return None;
    }
    Some(unsafe { std::slice::from_raw_parts(ptr, len) })
}

/// Parse a PDB-format string into an `Assembly`.
///
/// Returns null on failure with the error message available via
/// [`molex_last_error_message`]. The caller owns the returned handle
/// and must free it with [`molex_assembly_free`].
#[no_mangle]
pub extern "C" fn molex_pdb_str_to_assembly(
    str_ptr: *const c_char,
    len: usize,
) -> *mut molex_Assembly {
    let Some(bytes) = slice_from_raw(str_ptr.cast::<u8>(), len) else {
        set_last_error(&"molex_pdb_str_to_assembly: null input pointer");
        return std::ptr::null_mut();
    };
    let s = match std::str::from_utf8(bytes) {
        Ok(s) => s,
        Err(e) => {
            set_last_error(&format!("invalid UTF-8 in PDB input: {e}"));
            return std::ptr::null_mut();
        }
    };
    match pdb::pdb_str_to_entities(s) {
        Ok(entities) => {
            clear_last_error();
            assembly_from_entities(entities)
        }
        Err(e) => {
            set_last_error(&e);
            std::ptr::null_mut()
        }
    }
}

/// Parse an mmCIF-format string into an `Assembly`.
///
/// Returns null on failure with the error message available via
/// [`molex_last_error_message`]. The caller owns the returned handle
/// and must free it with [`molex_assembly_free`].
#[no_mangle]
pub extern "C" fn molex_cif_str_to_assembly(
    str_ptr: *const c_char,
    len: usize,
) -> *mut molex_Assembly {
    let Some(bytes) = slice_from_raw(str_ptr.cast::<u8>(), len) else {
        set_last_error(&"molex_cif_str_to_assembly: null input pointer");
        return std::ptr::null_mut();
    };
    let s = match std::str::from_utf8(bytes) {
        Ok(s) => s,
        Err(e) => {
            set_last_error(&format!("invalid UTF-8 in mmCIF input: {e}"));
            return std::ptr::null_mut();
        }
    };
    match cif::mmcif_str_to_entities(s) {
        Ok(entities) => {
            clear_last_error();
            assembly_from_entities(entities)
        }
        Err(e) => {
            set_last_error(&e);
            std::ptr::null_mut()
        }
    }
}

/// Decode BinaryCIF bytes into an `Assembly`.
///
/// Returns null on failure with the error message available via
/// [`molex_last_error_message`]. The caller owns the returned handle
/// and must free it with [`molex_assembly_free`].
#[no_mangle]
pub extern "C" fn molex_bcif_to_assembly(
    bytes_ptr: *const u8,
    len: usize,
) -> *mut molex_Assembly {
    let Some(bytes) = slice_from_raw(bytes_ptr, len) else {
        set_last_error(&"molex_bcif_to_assembly: null input pointer");
        return std::ptr::null_mut();
    };
    match bcif::bcif_to_entities(bytes) {
        Ok(entities) => {
            clear_last_error();
            assembly_from_entities(entities)
        }
        Err(e) => {
            set_last_error(&e);
            std::ptr::null_mut()
        }
    }
}

// ---------------------------------------------------------------------------
// Writer entry points (Assembly -> PDB)
// ---------------------------------------------------------------------------

fn vec_to_out_buffer(
    bytes: Vec<u8>,
    out_buf: *mut *mut u8,
    out_len: *mut usize,
) -> i32 {
    let len = bytes.len();
    let mut boxed = bytes.into_boxed_slice();
    let ptr = boxed.as_mut_ptr();
    std::mem::forget(boxed);
    unsafe {
        *out_buf = ptr;
        *out_len = len;
    }
    MOLEX_OK
}

/// Emit an `Assembly` as a PDB-format byte buffer.
///
/// On success returns [`MOLEX_OK`] and writes the heap-allocated buffer
/// pointer + length to `out_buf` / `out_len`; the caller frees with
/// [`molex_free_bytes`]. On failure returns a nonzero status and the
/// error message is available via [`molex_last_error_message`].
#[no_mangle]
pub extern "C" fn molex_assembly_to_pdb(
    assembly: *const molex_Assembly,
    out_buf: *mut *mut u8,
    out_len: *mut usize,
) -> i32 {
    if out_buf.is_null() || out_len.is_null() {
        set_last_error(&"molex_assembly_to_pdb: null pointer argument");
        return MOLEX_ERR_NULL;
    }
    let Some(assembly) = assembly_inner(assembly) else {
        set_last_error(&"molex_assembly_to_pdb: null pointer argument");
        return MOLEX_ERR_NULL;
    };
    match pdb::assembly_to_pdb(assembly) {
        Ok(s) => {
            clear_last_error();
            vec_to_out_buffer(s.into_bytes(), out_buf, out_len)
        }
        Err(e) => {
            set_last_error(&e);
            MOLEX_ERR
        }
    }
}

// ---------------------------------------------------------------------------
// Assembly walk accessors
// ---------------------------------------------------------------------------

/// Monotonic generation counter; increments on every mutation.
///
/// Use this for cheap change detection on consumer side without
/// snapshotting the full atom set.
#[no_mangle]
pub extern "C" fn molex_assembly_generation(
    assembly: *const molex_Assembly,
) -> u64 {
    assembly_inner(assembly).map_or(0, Assembly::generation)
}

/// Number of entities in the assembly. Returns 0 if `assembly` is null.
#[no_mangle]
pub extern "C" fn molex_assembly_num_entities(
    assembly: *const molex_Assembly,
) -> usize {
    assembly_inner(assembly).map_or(0, |a| a.entities().len())
}

/// Borrow a non-owning view of the i-th entity. Returns null when
/// `assembly` is null or `i` is out of bounds.
#[no_mangle]
pub extern "C" fn molex_assembly_entity(
    assembly: *const molex_Assembly,
    i: usize,
) -> *const molex_Entity {
    let Some(a) = assembly_inner(assembly) else {
        return std::ptr::null();
    };
    a.entities().get(i).map_or(std::ptr::null(), |entity_arc| {
        let entity_ref: &MoleculeEntity = Arc::as_ref(entity_arc);
        std::ptr::from_ref(entity_ref).cast::<molex_Entity>()
    })
}

// ---------------------------------------------------------------------------
// Entity walk accessors
// ---------------------------------------------------------------------------

fn entity_inner<'a>(entity: *const molex_Entity) -> Option<&'a MoleculeEntity> {
    if entity.is_null() {
        return None;
    }
    Some(unsafe { &*entity.cast::<MoleculeEntity>() })
}

/// Raw entity id (`u32`). Returns 0 if `entity` is null.
#[no_mangle]
pub extern "C" fn molex_entity_id(entity: *const molex_Entity) -> u32 {
    entity_inner(entity).map_or(0, |e| e.id().raw())
}

/// Molecule type discriminant. Returns [`molex_MoleculeType::Solvent`]'s
/// integer value (8) as a placeholder if `entity` is null - callers
/// should null-check the handle first.
#[no_mangle]
pub extern "C" fn molex_entity_molecule_type(
    entity: *const molex_Entity,
) -> molex_MoleculeType {
    entity_inner(entity)
        .map_or(molex_MoleculeType::Solvent, |e| e.molecule_type().into())
}

/// PDB chain identifier byte for polymer entities. Returns -1 when the
/// entity has no chain id (small molecule / bulk) or when `entity` is null.
#[no_mangle]
pub extern "C" fn molex_entity_pdb_chain_id(
    entity: *const molex_Entity,
) -> i32 {
    entity_inner(entity)
        .and_then(MoleculeEntity::pdb_chain_id)
        .map_or(-1, i32::from)
}

/// Total atom count in this entity. Returns 0 if `entity` is null.
#[no_mangle]
pub extern "C" fn molex_entity_num_atoms(entity: *const molex_Entity) -> usize {
    entity_inner(entity).map_or(0, MoleculeEntity::atom_count)
}

/// Borrow a non-owning view of the i-th atom in this entity's flat atom
/// list. Returns null when `entity` is null or `i` is out of bounds.
#[no_mangle]
pub extern "C" fn molex_entity_atom(
    entity: *const molex_Entity,
    i: usize,
) -> *const molex_Atom {
    let Some(e) = entity_inner(entity) else {
        return std::ptr::null();
    };
    e.atom_set().get(i).map_or(std::ptr::null(), |atom: &Atom| {
        std::ptr::from_ref(atom).cast::<molex_Atom>()
    })
}

fn polymer_residues(entity: &MoleculeEntity) -> Option<&[Residue]> {
    match entity {
        MoleculeEntity::Protein(p) => Some(&p.residues),
        MoleculeEntity::NucleicAcid(n) => Some(&n.residues),
        MoleculeEntity::SmallMolecule(_) | MoleculeEntity::Bulk(_) => None,
    }
}

/// Number of indexable residues in this entity.
///
/// Returns the residue count for protein and nucleic acid entities; 0 for
/// small-molecule and bulk entities (which do not expose individual
/// residue records). Returns 0 if `entity` is null.
#[no_mangle]
pub extern "C" fn molex_entity_num_residues(
    entity: *const molex_Entity,
) -> usize {
    entity_inner(entity)
        .and_then(polymer_residues)
        .map_or(0, <[Residue]>::len)
}

/// Borrow a non-owning view of the i-th residue in this entity.
///
/// Returns null when `entity` is null, the entity is not a polymer, or
/// `i` is out of bounds.
#[no_mangle]
pub extern "C" fn molex_entity_residue(
    entity: *const molex_Entity,
    i: usize,
) -> *const molex_Residue {
    let Some(e) = entity_inner(entity) else {
        return std::ptr::null();
    };
    polymer_residues(e)
        .and_then(|residues| residues.get(i))
        .map_or(std::ptr::null(), |residue: &Residue| {
            std::ptr::from_ref(residue).cast::<molex_Residue>()
        })
}

/// Number of atoms in the i-th residue of this entity. Returns 0 if
/// `entity` is null, the entity has no residue records, or `i` is out
/// of bounds.
#[no_mangle]
pub extern "C" fn molex_entity_residue_num_atoms(
    entity: *const molex_Entity,
    i: usize,
) -> usize {
    let Some(e) = entity_inner(entity) else {
        return 0;
    };
    polymer_residues(e)
        .and_then(|residues| residues.get(i))
        .map_or(0, |residue| residue.atom_range.len())
}

/// Borrow a non-owning view of the j-th atom in the i-th residue of
/// this entity. Returns null on any out-of-bounds access or null
/// `entity` argument.
#[no_mangle]
pub extern "C" fn molex_entity_residue_atom(
    entity: *const molex_Entity,
    residue_idx: usize,
    atom_idx: usize,
) -> *const molex_Atom {
    let Some(e) = entity_inner(entity) else {
        return std::ptr::null();
    };
    let Some(residues) = polymer_residues(e) else {
        return std::ptr::null();
    };
    let Some(residue) = residues.get(residue_idx) else {
        return std::ptr::null();
    };
    let atoms = e.atom_set();
    let flat_idx = residue.atom_range.start + atom_idx;
    if flat_idx >= residue.atom_range.end {
        return std::ptr::null();
    }
    atoms.get(flat_idx).map_or(std::ptr::null(), |atom: &Atom| {
        std::ptr::from_ref(atom).cast::<molex_Atom>()
    })
}

// ---------------------------------------------------------------------------
// Residue accessors
// ---------------------------------------------------------------------------

fn residue_inner<'a>(residue: *const molex_Residue) -> Option<&'a Residue> {
    if residue.is_null() {
        return None;
    }
    Some(unsafe { &*residue.cast::<Residue>() })
}

/// Pointer to this residue's 3-byte name (e.g. b"ALA"). Writes 3 to
/// `out_len` on success. Returns null and writes 0 if `residue` is null.
///
/// The buffer is space-padded to 3 bytes; callers that want a trimmed
/// name should strip trailing ASCII spaces.
#[no_mangle]
pub extern "C" fn molex_residue_name(
    residue: *const molex_Residue,
    out_len: *mut usize,
) -> *const u8 {
    let Some(r) = residue_inner(residue) else {
        if !out_len.is_null() {
            unsafe {
                *out_len = 0;
            }
        }
        return std::ptr::null();
    };
    if !out_len.is_null() {
        unsafe {
            *out_len = r.name.len();
        }
    }
    r.name.as_ptr()
}

/// Author-side sequence id, falling back to the structural-side
/// (`label_seq_id`) when the author id is absent.
#[no_mangle]
pub extern "C" fn molex_residue_seq_id(residue: *const molex_Residue) -> i32 {
    residue_inner(residue)
        .map_or(0, |r| r.auth_seq_id.unwrap_or(r.label_seq_id))
}

/// Structural-side sequence id (`label_seq_id`). Use this for stable
/// internal ordering rather than display.
#[no_mangle]
pub extern "C" fn molex_residue_label_seq_id(
    residue: *const molex_Residue,
) -> i32 {
    residue_inner(residue).map_or(0, |r| r.label_seq_id)
}

/// PDB insertion code byte (`iCode` / `pdbx_PDB_ins_code`). Returns 0
/// when the residue has no insertion code or `residue` is null.
#[no_mangle]
pub extern "C" fn molex_residue_ins_code(residue: *const molex_Residue) -> u8 {
    residue_inner(residue).and_then(|r| r.ins_code).unwrap_or(0)
}

// ---------------------------------------------------------------------------
// Atom accessors
// ---------------------------------------------------------------------------

fn atom_inner<'a>(atom: *const molex_Atom) -> Option<&'a Atom> {
    if atom.is_null() {
        return None;
    }
    Some(unsafe { &*atom.cast::<Atom>() })
}

/// Pointer to this atom's 4-byte PDB-style name (e.g. b"CA  "). Writes
/// 4 to `out_len` on success. Returns null and writes 0 if `atom` is null.
///
/// The buffer is space-padded; callers that want a trimmed atom name
/// should strip ASCII spaces.
#[no_mangle]
pub extern "C" fn molex_atom_name(
    atom: *const molex_Atom,
    out_len: *mut usize,
) -> *const u8 {
    let Some(a) = atom_inner(atom) else {
        if !out_len.is_null() {
            unsafe {
                *out_len = 0;
            }
        }
        return std::ptr::null();
    };
    if !out_len.is_null() {
        unsafe {
            *out_len = a.name.len();
        }
    }
    a.name.as_ptr()
}

/// Atomic number for this atom's element, or 0 for
/// [`Element::Unknown`] / a null `atom`.
#[no_mangle]
pub extern "C" fn molex_atom_atomic_number(atom: *const molex_Atom) -> u8 {
    atom_inner(atom).map_or(0, |a| Element::atomic_number(a.element))
}

/// Write this atom's `(x, y, z)` position into the 3-float output array.
/// No-op if either pointer is null.
#[no_mangle]
pub extern "C" fn molex_atom_position(
    atom: *const molex_Atom,
    out_xyz: *mut f32,
) {
    if out_xyz.is_null() {
        return;
    }
    let Some(a) = atom_inner(atom) else { return };
    unsafe {
        *out_xyz.add(0) = a.position.x;
        *out_xyz.add(1) = a.position.y;
        *out_xyz.add(2) = a.position.z;
    }
}

/// Crystallographic occupancy (0.0 to 1.0). Returns 0 if `atom` is null.
#[no_mangle]
pub extern "C" fn molex_atom_occupancy(atom: *const molex_Atom) -> f32 {
    atom_inner(atom).map_or(0.0, |a| a.occupancy)
}

/// Temperature factor (B-factor) in square angstroms. Returns 0 if
/// `atom` is null.
#[no_mangle]
pub extern "C" fn molex_atom_b_factor(atom: *const molex_Atom) -> f32 {
    atom_inner(atom).map_or(0.0, |a| a.b_factor)
}

/// Signed formal charge (0 means neutral). Returns 0 if `atom` is null.
#[no_mangle]
pub extern "C" fn molex_atom_formal_charge(atom: *const molex_Atom) -> i8 {
    atom_inner(atom).map_or(0, |a| a.formal_charge)
}

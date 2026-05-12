//! C ABI walk accessors: assembly -> entity -> residue / atom.
//!
//! Split out of `c_api/mod.rs` for file-length reasons. All entries
//! here are `#[no_mangle] pub extern "C"` and so emit C symbols
//! identically to functions declared at the parent module scope.
//! cbindgen discovers them via the parent module's `pub mod walk;`.

#![allow(
    clippy::missing_safety_doc,
    clippy::not_unsafe_ptr_arg_deref,
    reason = "Inherits the FFI-surface conventions documented at the top of \
              c_api/mod.rs."
)]

use std::sync::Arc;

use super::{
    assembly_inner, molex_Assembly, molex_Atom, molex_Entity,
    molex_MoleculeType, molex_Residue,
};
use crate::assembly::Assembly;
use crate::element::Element;
use crate::entity::molecule::{Atom, MoleculeEntity, Residue};

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

/// Pointer to the single 3-byte residue name carried by a non-polymer
/// entity (`SmallMolecule` / `Bulk`). Writes 3 to `out_len` on success;
/// returns null and writes 0 for polymers or a null `entity`.
///
/// The buffer is space-padded to 3 bytes; callers should strip trailing
/// ASCII spaces if needed.
#[no_mangle]
pub extern "C" fn molex_entity_residue_name_single(
    entity: *const molex_Entity,
    out_len: *mut usize,
) -> *const u8 {
    let write_len = |len: usize| {
        if !out_len.is_null() {
            unsafe {
                *out_len = len;
            }
        }
    };
    let Some(e) = entity_inner(entity) else {
        write_len(0);
        return std::ptr::null();
    };
    let name: &[u8; 3] = match e {
        MoleculeEntity::SmallMolecule(s) => &s.residue_name,
        MoleculeEntity::Bulk(b) => &b.residue_name,
        MoleculeEntity::Protein(_) | MoleculeEntity::NucleicAcid(_) => {
            write_len(0);
            return std::ptr::null();
        }
    };
    write_len(name.len());
    name.as_ptr()
}

/// Number of equal-sized molecule chunks the atom set should be split into.
///
/// For non-polymer entities only: returns 1 for `SmallMolecule`,
/// `BulkEntity::molecule_count` for `Bulk`, and 0 for polymers or a null
/// `entity`.
#[no_mangle]
pub extern "C" fn molex_entity_molecule_count(
    entity: *const molex_Entity,
) -> usize {
    let Some(e) = entity_inner(entity) else {
        return 0;
    };
    match e {
        MoleculeEntity::SmallMolecule(_) => 1,
        MoleculeEntity::Bulk(b) => b.molecule_count,
        MoleculeEntity::Protein(_) | MoleculeEntity::NucleicAcid(_) => 0,
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

//! C-compatible FFI layer for COORDS conversion functions.

use std::ffi::{c_char, CString};
use std::ptr;

/// FFI-safe result type for COORDS conversion functions.
#[repr(C)]
pub struct CoordsResult {
    /// Pointer to the output byte buffer, or null on error.
    pub data: *const u8,
    /// Length of the data buffer in bytes.
    pub len: usize,
    /// Allocated capacity of the data buffer.
    pub data_len: usize,
    /// Pointer to a null-terminated error string, or null on success.
    pub error: *const c_char,
}

impl CoordsResult {
    /// Create a successful result from owned bytes.
    #[must_use]
    pub fn success(data: Vec<u8>) -> Self {
        let len = data.len();
        let data_ptr = if data.is_empty() {
            ptr::null()
        } else {
            let boxed = data.into_boxed_slice();
            let ptr = boxed.as_ptr();
            std::mem::forget(boxed);
            ptr
        };
        Self {
            data: data_ptr,
            len,
            data_len: len,
            error: ptr::null(),
        }
    }

    /// Create an error result with a message string.
    ///
    /// # Panics
    ///
    /// Panics if `msg` contains an interior null byte.
    #[must_use]
    #[allow(clippy::expect_used)]
    pub fn error(msg: &str) -> Self {
        Self {
            data: ptr::null(),
            len: 0,
            data_len: 0,
            error: CString::new(msg)
                .expect("error message must not contain null bytes")
                .into_raw(),
        }
    }
}

/// # Safety
/// `result` must be a valid pointer returned by a `coords_*` FFI function, or
/// null.
#[no_mangle]
pub unsafe extern "C" fn coords_free_result(result: *const CoordsResult) {
    if result.is_null() {
        return;
    }
    let r = &*result;
    if !r.data.is_null() && r.len > 0 {
        let _: Box<[u8]> =
            Box::from(core::slice::from_raw_parts(r.data, r.len));
    }
    if !r.error.is_null() {
        let _ = CString::from_raw(r.error.cast_mut().cast());
    }
}

/// Free a C string returned by a `coords_*` FFI function.
#[no_mangle]
pub extern "C" fn coords_free_string(s: *const c_char) {
    if !s.is_null() {
        unsafe {
            let _ = CString::from_raw(s.cast_mut());
        }
    }
}

/// Parse a PDB string into COORDS binary format.
#[no_mangle]
pub extern "C" fn pdb_to_coords_bytes(
    pdb_ptr: *const c_char,
    pdb_len: usize,
) -> CoordsResult {
    if pdb_ptr.is_null() {
        return CoordsResult::error("PDB string is null");
    }
    unsafe {
        let pdb_slice =
            std::slice::from_raw_parts(pdb_ptr.cast::<u8>(), pdb_len);
        let Ok(pdb_str) = std::str::from_utf8(pdb_slice) else {
            return CoordsResult::error("Invalid UTF-8 in PDB string");
        };
        match crate::adapters::pdb::pdb_to_coords(pdb_str) {
            Ok(coords_bytes) => CoordsResult::success(coords_bytes),
            Err(e) => CoordsResult::error(&e.to_string()),
        }
    }
}

/// # Safety
/// `coords_ptr` must point to `coords_len` valid bytes. `out_len` must be a
/// valid, non-null pointer.
///
/// # Panics
///
/// Panics if the resulting PDB string or error message contains interior null
/// bytes.
#[no_mangle]
#[allow(clippy::expect_used)]
pub unsafe extern "C" fn coords_to_pdb(
    coords_ptr: *const u8,
    coords_len: usize,
    out_len: *mut usize,
) -> *const c_char {
    if coords_ptr.is_null() {
        return CString::new("COORDS data is null")
            .expect("static string has no null bytes")
            .into_raw();
    }
    if out_len.is_null() {
        return CString::new("out_len is null")
            .expect("static string has no null bytes")
            .into_raw();
    }
    let coords_slice = std::slice::from_raw_parts(coords_ptr, coords_len);
    match crate::adapters::pdb::coords_to_pdb(coords_slice) {
        Ok(pdb_string) => {
            *out_len = pdb_string.len();
            CString::new(pdb_string)
                .expect("PDB string must not contain null bytes")
                .into_raw()
        }
        Err(e) => CString::new(e.to_string())
            .expect("error string must not contain null bytes")
            .into_raw(),
    }
}

/// # Safety
/// `coords_ptr` must point to `coords_len` valid bytes.
#[no_mangle]
pub unsafe extern "C" fn coords_from_coords(
    coords_ptr: *const u8,
    coords_len: usize,
) -> CoordsResult {
    if coords_ptr.is_null() {
        return CoordsResult::error("Coords data is null");
    }
    let coords_slice = std::slice::from_raw_parts(coords_ptr, coords_len);
    match crate::ops::codec::deserialize(coords_slice)
        .and_then(|coords| crate::ops::codec::serialize(&coords))
    {
        Ok(coords_bytes) => CoordsResult::success(coords_bytes),
        Err(e) => CoordsResult::error(&e.to_string()),
    }
}

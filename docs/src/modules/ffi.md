# C FFI

C-compatible bindings live in `molex::ffi` (source: `src/ffi.rs`). These expose COORDS conversion functions for consumption from C, C++, Swift, or any language with C FFI support.

## Result type

```c
typedef struct {
    const uint8_t *data;    // output bytes, or NULL on error
    size_t len;             // data length
    size_t data_len;        // allocated capacity
    const char *error;      // error string, or NULL on success
} CoordsResult;
```

## Functions

### `pdb_to_coords_bytes`

Parse a PDB string into COORDS binary format.

```c
CoordsResult pdb_to_coords_bytes(const char *pdb_ptr, size_t pdb_len);
```

### `coords_to_pdb`

Convert COORDS binary to a PDB-format string. Returns a null-terminated C string. The caller must free the string with `coords_free_string`.

```c
const char *coords_to_pdb(const uint8_t *coords_ptr, size_t coords_len, size_t *out_len);
```

### `coords_from_coords`

Deserialize and re-serialize COORDS bytes (round-trip validation).

```c
CoordsResult coords_from_coords(const uint8_t *coords_ptr, size_t coords_len);
```

## Memory management

```c
void coords_free_result(const CoordsResult *result);
void coords_free_string(const char *s);
```

All pointers returned by FFI functions must be freed using the corresponding free function. Do not use `free()` directly.

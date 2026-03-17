# FFI — C Bindings

The `ffi` module provides `extern "C"` functions for integration with
C and C++ applications. These are re-exported by `foldit-runner` for
inclusion in its generated C header.

## Functions

### `coords_from_backbone`

```c
CoordsResult coords_from_backbone(
    const float* positions,  // [x,y,z] × n_atoms
    uint32_t n_atoms,
    const uint8_t* chain_ids,
    const int32_t* res_nums,
    const char (*res_names)[4],  // 3-letter codes, null-padded
    const char (*atom_names)[5]  // 4-letter codes, null-padded
);
```

Constructs COORDS01 bytes from raw arrays. Returns a `CoordsResult`
containing either serialized bytes or an error string.

### `coords_free_result`

Free memory allocated by `coords_from_backbone`.

### `coords_free_string`

Free an error string from a `CoordsResult`.

## `CoordsResult`

```c
typedef struct {
    const uint8_t* data;   // COORDS01 bytes (null on error)
    uint32_t len;          // byte length
    const char* error;     // error message (null on success)
} CoordsResult;
```

## Usage from C++

```cpp
#include "molex.h"

auto result = coords_from_backbone(
    positions.data(), n_atoms,
    chain_ids.data(), res_nums.data(),
    res_names.data(), atom_names.data()
);

if (result.error) {
    fprintf(stderr, "Error: %s\n", result.error);
    coords_free_string(result.error);
} else {
    // Use result.data[0..result.len]
    coords_free_result(result);
}
```

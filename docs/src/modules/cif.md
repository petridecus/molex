# CIF Parser

The `cif` module is a standalone CIF/STAR parser with typed data
extraction. It operates in two layers.

## Layer 1 — DOM parsing (`cif::parse`)

```rust,ignore
let doc = molex::cif::parse(input)?;
```

Parses any CIF or STAR file into an untyped `Document` tree:

- `Document` contains `Vec<Block>`
- Each `Block` has `name`, `categories` (loop tables), and `pairs`
  (key-value entries)
- `Category` holds `Column` data (string values, one per row)

This layer makes no assumptions about the content — it works for
mmCIF, CCD, reflection data, or any STAR-format file.

## Layer 2 — Typed extractors (`cif::extract`)

Pull structured data from a parsed `Block`:

```rust,ignore
use molex::cif::extract::{CoordinateData, CifContent};

// Caller knows the content type:
let coords = CoordinateData::try_from(&block)?;

// Or auto-detect:
match CifContent::try_from(&block)? {
    CifContent::Coordinates(data) => { /* atom_site data */ },
    CifContent::Reflections(data) => { /* refln data */ },
    CifContent::Dictionary(data)  => { /* CCD entry */ },
}
```

### `CoordinateData`

Extracted from `_atom_site` loops:

- Atom positions, names, elements, B-factors, occupancy
- Chain IDs, residue names, residue numbers
- Entity ID and molecule type annotations
- Cell parameters and space group (if present)

### `ReflectionData`

Extracted from `_refln` loops (X-ray diffraction data).

### `DictionaryEntry`

Extracted from CCD (Chemical Component Dictionary) entries — ideal
coordinates, bond tables, and chemical metadata.

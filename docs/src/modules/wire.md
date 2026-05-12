# Wire Format

`molex::ops::wire` is the home of the ASSEM01 binary wire format —
the on-disk / IPC format used to round-trip an `Assembly` through
bytes.

ASSEM01 is entity-aware: it preserves per-entity molecule-type
metadata so the decoder can reconstruct `MoleculeEntity` variants
without re-running residue classification. It replaced the legacy
COORDS01 format (which was retired entirely in the Phase 1/2 work
documented in `docs/COORDS_RETIREMENT_PLAN.md`).

## Public surface

```rust,ignore
use molex::ops::wire::{
    assembly_bytes, deserialize_assembly, serialize_assembly, ASSEMBLY_MAGIC,
};

// Serialize an Assembly
let bytes: Vec<u8> = serialize_assembly(&assembly)?;

// Serialize a raw entity slice (no Assembly wrapper)
let bytes: Vec<u8> = assembly_bytes(&entities)?;

// Deserialize back to an Assembly (recomputes derived data via Assembly::new)
let assembly: Assembly = deserialize_assembly(&bytes)?;

// Magic header for ASSEM01 detection
assert_eq!(&bytes[0..8], ASSEMBLY_MAGIC);
```

## Byte layout

```
8 bytes:   magic "ASSEM01\0"
4 bytes:   entity_count (u32 BE)
Per entity header (5 bytes each):
  1 byte:  molecule_type wire byte
  4 bytes: atom_count (u32 BE)
Per atom (26 bytes):
  12 bytes: x, y, z (f32 BE × 3)
  1 byte:   chain_id
  3 bytes:  res_name
  4 bytes:  res_num (i32 BE)
  4 bytes:  atom_name
  2 bytes:  element symbol (byte 0, byte 1 or 0)
```

Occupancy and b_factor are **not** preserved on the wire; deserialize
resets them to 1.0 and 0.0 respectively. Derived data (`ss_types`,
`hbonds`, `cross_entity_bonds`) is recomputed by `Assembly::new` at
deserialize time.

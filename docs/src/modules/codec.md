# Codec

`molex::ops::codec` is the home of the crate-internal `Coords` parser
intermediate and a small set of entity-utility helpers. It is **not**
the ASSEM01 wire format — that lives in `molex::ops::wire` (see
[Wire](wire.md)).

`Coords` is a flat parallel-array atom record used inside molex by
the PDB/CIF/BCIF parsers and the ASSEM01 decoder to bridge raw bytes
into `Vec<MoleculeEntity>`. It is `pub` only because
`update_protein_entities` still takes `&Coords`; that signature is
scheduled for retirement in Phase 4a of
`docs/COORDS_RETIREMENT_PLAN.md`. Do not grow new external usage of
`Coords` or `CoordsAtom`.

## Public surface

```rust,ignore
use molex::ops::codec::{ca_positions, residue_count, update_protein_entities};

// CA positions across all protein entities in an entity slice.
let ca: Vec<Vec3> = ca_positions(&entities);

// Total protein residue count.
let n: usize = residue_count(&entities);

// Replace the protein-typed entities in `entities` with the chains
// reconstructed from `protein_coords`. Currently the only consumer is
// viso's live-update wire; this function (and its `&Coords` argument)
// is slated to be replaced by an entity-shaped signature.
update_protein_entities(&mut entities, &protein_coords);
```

## Error type

All adapter and wire operations return `AdapterError`:

```rust,ignore
pub enum AdapterError {
    InvalidFormat(String),
    PdbParseError(String),
    SerializationError(String),
}
```

## Retired surface

The following functions used to live in `ops::codec` and have been
removed. Migration targets are listed for any out-of-tree callers:

| Removed                 | Replacement                                   |
| ----------------------- | --------------------------------------------- |
| `serialize`             | `ops::wire::serialize_assembly` (ASSEM01 only) |
| `deserialize`           | `ops::wire::deserialize_assembly`              |
| `serialize_assembly`    | `ops::wire::serialize_assembly`                |
| `deserialize_assembly`  | `ops::wire::deserialize_assembly`              |
| `assembly_bytes`        | `ops::wire::assembly_bytes`                    |
| `ASSEMBLY_MAGIC`        | `ops::wire::ASSEMBLY_MAGIC`                    |
| `merge_entities`        | (none; entity walks replace this)              |
| `extract_by_type`       | (none; filter `assembly.entities()` directly)  |
| `protein_coords`        | (none)                                         |
| `align_coords_bytes`    | (none; deleted with COORDS01)                  |

See `docs/COORDS_RETIREMENT_PLAN.md` for full context.

## Transforms (`ops::transform`)

Structural alignment and extraction utilities:

```rust,ignore
use molex::ops::transform::*;

// Kabsch alignment (minimize RMSD between point sets)
let (rotation, translation) = kabsch_alignment(&source, &target)?;
let (rotation, translation, scale) =
    kabsch_alignment_with_scale(&source, &target)?;

// Apply transform to entities
transform_entities(&mut entities, rotation, translation);
transform_entities_with_scale(&mut entities, rotation, translation, scale);

// Align entities to a reference structure
align_to_reference(&mut mobile, &reference);

// Extract alpha-carbon positions
let ca_positions: Vec<Vec3> = extract_ca_positions(&entities);
let ca_by_chain: Vec<Vec<Vec3>> = extract_ca_from_chains(&entities);

// Get continuous backbone segments
let segments: Vec<Vec<Vec3>> = extract_backbone_segments(&entities);

// Compute centroid
let center: Vec3 = centroid(&positions);
```

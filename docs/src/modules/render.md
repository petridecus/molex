# Render — Visualization Data

The `render` module extracts GPU-ready data from `Coords` and
`MoleculeEntity` values. It serves as the bridge between molex's
canonical types and rendering engines like viso.

## `RenderCoords`

The primary output type, containing separated backbone and sidechain
data with bond connectivity:

```rust,ignore
let render = RenderCoords::from_entity(&entity, is_hydrophobic, get_bonds);
```

### Fields

| Field | Type | Description |
|-------|------|-------------|
| `backbone_chains` | `Vec<Vec<Vec3>>` | N-CA-C position triples per chain |
| `backbone_chain_ids` | `Vec<u8>` | Chain ID for each backbone chain |
| `backbone_residue_chains` | `Vec<Vec<RenderBackboneResidue>>` | Full N/CA/C/O per residue |
| `sidechain_atoms` | `Vec<RenderSidechainAtom>` | Non-backbone heavy atoms |
| `sidechain_bonds` | `Vec<(u32, u32)>` | Sidechain bond pairs (indices into `sidechain_atoms`) |
| `backbone_sidechain_bonds` | `Vec<(Vec3, u32)>` | CA→CB connections |
| `all_positions` | `Vec<Vec3>` | Every atom position (for bounding box, picking) |

### Construction methods

- `from_entity()` — from a `MoleculeEntity` with topology callbacks
- `from_coords_with_topology()` — from raw `Coords` with callbacks
- `from_coords()` — minimal extraction (no bonds or hydrophobicity)

### Queries

```rust,ignore
render.get_atom_position(residue_idx, "CB")  // Option<Vec3>
render.find_closest_atom(residue_idx, point) // Option<(Vec3, String)>
render.ca_positions()                        // Vec<Vec3>
render.residue_count()                       // usize
```

## `render::backbone`

Backbone-specific extraction from `MoleculeEntity`:

- `BackboneData` — chains as flat CA arrays + chain IDs
- `ca_positions_from_chains()` — extract CA positions from N-CA-C chains

## `render::sidechain`

Sidechain extraction with bond topology:

- `SidechainAtoms` — atom positions, names, bonds, hydrophobicity
- `SidechainAtomData` — per-atom data struct

## `extract_sequences()`

Extract amino acid sequences from `Coords`:

```rust,ignore
let (full_sequence, chain_sequences) = extract_sequences(&coords);
// full_sequence: "MQIFVKTL..."
// chain_sequences: vec![(b'A', "MQIFVKTL..."), (b'B', "...")]
```

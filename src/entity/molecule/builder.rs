//! Builder that turns parser-output atom rows into typed `MoleculeEntity`
//! values. Used by the PDB, mmCIF, and BinaryCIF adapters.

#![allow(
    dead_code,
    reason = "internals are wired in by adapters during the Phase 3a/3b/3c \
              cutover"
)]

use std::collections::HashMap;

use glam::Vec3;

use super::atom::Atom;
use super::bulk::BulkEntity;
use super::classify::classify_residue;
use super::id::EntityIdAllocator;
use super::nucleic_acid::NAEntity;
use super::polymer::Residue;
use super::protein::ProteinEntity;
use super::small_molecule::SmallMoleculeEntity;
use super::{MoleculeEntity, MoleculeType};
use crate::element::Element;
use crate::ops::codec::ChainIdMapper;

/// Maximum number of distinct chains the builder accepts. Matches the
/// printable-byte capacity of `ChainIdMapper`.
const MAX_CHAINS: usize = 90;
const DEFAULT_WATER_RESNAME: [u8; 3] = *b"HOH";
const DEFAULT_SOLVENT_RESNAME: [u8; 3] = *b"GOL";

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Hint about an entity's classification, sourced from mmCIF
/// `_entity.type` joined with `_entity_poly.type`. The PDB ingest path
/// supplies `Unknown` for every atom; mmCIF/BCIF supply explicit hints.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(
    clippy::upper_case_acronyms,
    reason = "DNA / RNA mirror MoleculeType's variant names"
)]
pub(crate) enum ExpectedEntityType {
    /// `_entity_poly.type` = `polypeptide(L)` or `polypeptide(D)`.
    Protein,
    /// `_entity_poly.type` = `polydeoxyribonucleotide`.
    DNA,
    /// `_entity_poly.type` = `polyribonucleotide`.
    RNA,
    /// `_entity.type` = `non-polymer` or `branched`. Builder
    /// sub-classifies each residue via the 3-letter heuristic.
    NonPolymer,
    /// `_entity.type` = `water`.
    Water,
    /// Entity tables absent (PDB path) or unrecognized polymer subtype.
    /// Builder runs the full heuristic plus modified-residue merge.
    Unknown,
}

/// One parser-emitted atom row.
pub(crate) struct AtomRow {
    /// `label_asym_id` (mmCIF) or chain letter (PDB).
    pub label_asym_id: String,
    /// `label_seq_id` (mmCIF) or `resSeq` (PDB). For non-polymer
    /// entities with `.` in the source, parsers pass `1`.
    pub label_seq_id: i32,
    /// `label_comp_id` (mmCIF) or `resName` (PDB), 3 chars, space-padded.
    pub label_comp_id: [u8; 3],
    /// `label_atom_id` (mmCIF) or atom name (PDB), 4 chars.
    pub label_atom_id: [u8; 4],
    /// `label_entity_id` (mmCIF). Joined against registered hints to
    /// pick the entity type. `None` on the PDB path.
    pub label_entity_id: Option<String>,

    /// `auth_asym_id`. `None` defaults to `label_asym_id`.
    pub auth_asym_id: Option<String>,
    /// `auth_seq_id`. `None` defaults to `label_seq_id`.
    pub auth_seq_id: Option<i32>,
    /// `auth_comp_id`. `None` defaults to `label_comp_id`.
    pub auth_comp_id: Option<[u8; 3]>,
    /// `auth_atom_id`. `None` defaults to `label_atom_id`.
    pub auth_atom_id: Option<[u8; 4]>,

    /// Alternate location indicator. `None` = blank (present in every
    /// conformer).
    pub alt_loc: Option<u8>,
    /// Insertion code. `None` = blank.
    pub ins_code: Option<u8>,
    /// Chemical element.
    pub element: Element,
    /// X coordinate in angstroms.
    pub x: f32,
    /// Y coordinate in angstroms.
    pub y: f32,
    /// Z coordinate in angstroms.
    pub z: f32,
    /// Crystallographic occupancy (0.0..=1.0).
    pub occupancy: f32,
    /// Temperature factor (B-factor) in square angstroms.
    pub b_factor: f32,
    /// Formal charge. 0 default.
    pub formal_charge: i8,
}

/// Builder-side errors. Adapters wrap these into `AdapterError`.
#[derive(Debug, thiserror::Error)]
pub(crate) enum BuildError {
    /// An atom row contains NaN or infinite coordinates.
    #[error("invalid coordinate: {axis} = {value} (atom {label_atom_id})")]
    InvalidCoordinate {
        /// Offending axis ('x', 'y', or 'z').
        axis: char,
        /// Offending value.
        value: f32,
        /// Trimmed atom name from the offending row.
        label_atom_id: String,
    },
    /// Chain count exceeded printable-byte capacity. Adapters convert
    /// this to the mmCIF redirect error.
    #[error("structure exceeds {limit} chains; use mmCIF instead")]
    TooManyChains {
        /// Capacity that was exceeded.
        limit: usize,
    },
}

/// Builds typed `MoleculeEntity` values from a stream of atom rows.
///
/// Usage:
///
/// 1. [`EntityBuilder::new`].
/// 2. (mmCIF/BCIF only) [`EntityBuilder::register_entity`] once per entity from
///    the `_entity` / `_entity_poly` pre-pass.
/// 3. [`EntityBuilder::push_atom`] for each ATOM/HETATM row in parser order.
/// 4. [`EntityBuilder::finish`] to consume the builder.
pub(crate) struct EntityBuilder {
    allocator: EntityIdAllocator,
    hints: HashMap<String, ExpectedEntityType>,
    chains: HashMap<String, ChainState>,
    chain_mapper: ChainIdMapper,
    chain_order: Vec<String>,
}

impl EntityBuilder {
    /// Create an empty builder.
    pub(crate) fn new() -> Self {
        Self {
            allocator: EntityIdAllocator::new(),
            hints: HashMap::new(),
            chains: HashMap::new(),
            chain_mapper: ChainIdMapper::new(),
            chain_order: Vec::new(),
        }
    }

    /// Register an entity-type hint sourced from mmCIF `_entity` +
    /// `_entity_poly`. Idempotent on repeat with the same hint;
    /// first-wins on conflicting hints (with a debug log).
    #[allow(
        clippy::needless_pass_by_value,
        reason = "ExpectedEntityType is Copy; by-value matches the surface"
    )]
    pub(crate) fn register_entity(
        &mut self,
        label_entity_id: &str,
        hint: ExpectedEntityType,
    ) {
        if let Some(&existing) = self.hints.get(label_entity_id) {
            if existing != hint {
                log::debug!(
                    "EntityBuilder: ignoring conflicting hint {hint:?} for \
                     entity {label_entity_id} (first registration: \
                     {existing:?})",
                );
            }
            return;
        }
        let _ = self.hints.insert(label_entity_id.to_owned(), hint);
    }

    /// Push an atom row. Handles altLoc dedup, residue grouping, and
    /// chain bucketing.
    ///
    /// Errors only on invalid coordinates or chain-mapper overflow;
    /// rows with unknown elements or unrecognized residue names are
    /// accepted.
    #[allow(
        clippy::needless_pass_by_value,
        reason = "AtomRow is moved into builder state in 3b/3c"
    )]
    pub(crate) fn push_atom(&mut self, row: AtomRow) -> Result<(), BuildError> {
        validate_coords(&row)?;
        self.ensure_chain(&row)?;
        let Some(chain) = self.chains.get_mut(&row.label_asym_id) else {
            unreachable!("chain inserted by ensure_chain");
        };
        let res_idx = chain.locate_or_create_residue(&row);
        chain.residues[res_idx].apply_altloc_dedup(&row);
        Ok(())
    }

    /// Consume the builder and produce `MoleculeEntity` values.
    /// Chain insertion order; global Water and Solvent bulks
    /// (accumulated under Unknown / NonPolymer hints) appended last.
    #[allow(
        clippy::unnecessary_wraps,
        reason = "Result is part of the API contract; future error paths land \
                  here"
    )]
    pub(crate) fn finish(self) -> Result<Vec<MoleculeEntity>, BuildError> {
        let Self {
            mut allocator,
            hints,
            mut chains,
            chain_mapper: _,
            chain_order,
        } = self;

        let mut out: Vec<MoleculeEntity> = Vec::new();
        let mut water = GlobalBulk::default();
        let mut solvent = GlobalBulk::default();

        {
            let mut ctx = ChainCtx {
                allocator: &mut allocator,
                out: &mut out,
                water: &mut water,
                solvent: &mut solvent,
            };
            for chain_key in &chain_order {
                let Some(mut chain) = chains.remove(chain_key) else {
                    unreachable!("chain in chain_order");
                };
                chain.residues.sort_by_key(|r| (r.label_seq_id, r.ins_code));
                let hint = chain
                    .entity_hint_key
                    .as_deref()
                    .and_then(|k| hints.get(k).copied())
                    .unwrap_or(ExpectedEntityType::Unknown);
                classify_chain(
                    hint,
                    chain.pdb_chain_id,
                    &chain.residues,
                    &mut ctx,
                );
            }
        }

        water.flush(
            MoleculeType::Water,
            DEFAULT_WATER_RESNAME,
            &mut allocator,
            &mut out,
        );
        solvent.flush(
            MoleculeType::Solvent,
            DEFAULT_SOLVENT_RESNAME,
            &mut allocator,
            &mut out,
        );

        Ok(out)
    }

    fn ensure_chain(&mut self, row: &AtomRow) -> Result<(), BuildError> {
        if self.chains.contains_key(&row.label_asym_id) {
            return Ok(());
        }
        if self.chain_order.len() >= MAX_CHAINS {
            return Err(BuildError::TooManyChains { limit: MAX_CHAINS });
        }
        let pdb_chain_id = self.chain_mapper.get_or_assign(&row.label_asym_id);
        let state = ChainState {
            pdb_chain_id,
            entity_hint_key: row.label_entity_id.clone(),
            residues: Vec::new(),
            residue_index: HashMap::new(),
        };
        // Cloning the chain key: needed both as HashMap key and as the
        // chain_order entry that drives finish() ordering.
        let key = row.label_asym_id.clone();
        self.chain_order.push(key.clone());
        let _ = self.chains.insert(key, state);
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Internal state
// ---------------------------------------------------------------------------

struct ChainState {
    pdb_chain_id: u8,
    entity_hint_key: Option<String>,
    residues: Vec<ResidueAccum>,
    residue_index: HashMap<(i32, Option<u8>), usize>,
}

impl ChainState {
    fn locate_or_create_residue(&mut self, row: &AtomRow) -> usize {
        let key = (row.label_seq_id, row.ins_code);
        if let Some(&idx) = self.residue_index.get(&key) {
            return idx;
        }
        let idx = self.residues.len();
        self.residues.push(ResidueAccum {
            label_seq_id: row.label_seq_id,
            label_comp_id: row.label_comp_id,
            ins_code: row.ins_code,
            atoms: HashMap::new(),
            atom_order: Vec::new(),
        });
        let _ = self.residue_index.insert(key, idx);
        idx
    }
}

struct ResidueAccum {
    label_seq_id: i32,
    label_comp_id: [u8; 3],
    ins_code: Option<u8>,
    atoms: HashMap<[u8; 4], AtomChoice>,
    atom_order: Vec<[u8; 4]>,
}

impl ResidueAccum {
    fn apply_altloc_dedup(&mut self, row: &AtomRow) {
        let candidate = AtomChoice::from_row(row);
        match self.atoms.get_mut(&row.label_atom_id) {
            None => {
                let _ = self.atoms.insert(row.label_atom_id, candidate);
                self.atom_order.push(row.label_atom_id);
            }
            Some(existing) => {
                if candidate_should_replace(existing, &candidate) {
                    *existing = candidate;
                }
            }
        }
    }
}

struct AtomChoice {
    alt_loc: Option<u8>,
    occupancy: f32,
    label_atom_id: [u8; 4],
    element: Element,
    position: Vec3,
    b_factor: f32,
}

impl AtomChoice {
    fn from_row(row: &AtomRow) -> Self {
        Self {
            alt_loc: row.alt_loc,
            occupancy: row.occupancy,
            label_atom_id: row.label_atom_id,
            element: row.element,
            position: Vec3::new(row.x, row.y, row.z),
            b_factor: row.b_factor,
        }
    }

    fn to_atom(&self) -> Atom {
        Atom {
            position: self.position,
            occupancy: self.occupancy,
            b_factor: self.b_factor,
            element: self.element,
            name: self.label_atom_id,
        }
    }
}

/// AltLoc dedup precedence: blank trumps any non-blank; among non-blank
/// higher occupancy wins; ties go to the alphabetically-earlier byte.
fn candidate_should_replace(
    existing: &AtomChoice,
    candidate: &AtomChoice,
) -> bool {
    match (existing.alt_loc, candidate.alt_loc) {
        (None, _) => false,
        (Some(_), None) => true,
        (Some(_), Some(_)) => {
            if candidate.occupancy > existing.occupancy {
                true
            } else if candidate.occupancy < existing.occupancy {
                false
            } else {
                candidate.alt_loc < existing.alt_loc
            }
        }
    }
}

// ---------------------------------------------------------------------------
// finish() context
// ---------------------------------------------------------------------------

struct ChainCtx<'a> {
    allocator: &'a mut EntityIdAllocator,
    out: &'a mut Vec<MoleculeEntity>,
    water: &'a mut GlobalBulk,
    solvent: &'a mut GlobalBulk,
}

#[derive(Default)]
struct GlobalBulk {
    atoms: Vec<Atom>,
    residue_count: usize,
    residue_name: Option<[u8; 3]>,
}

impl GlobalBulk {
    fn ingest(&mut self, r: &ResidueAccum) {
        for name in &r.atom_order {
            self.atoms.push(r.atoms[name].to_atom());
        }
        self.residue_count += 1;
        if self.residue_name.is_none() {
            self.residue_name = Some(r.label_comp_id);
        }
    }

    fn flush(
        self,
        mol_type: MoleculeType,
        default_resname: [u8; 3],
        allocator: &mut EntityIdAllocator,
        out: &mut Vec<MoleculeEntity>,
    ) {
        if self.atoms.is_empty() {
            return;
        }
        let id = allocator.allocate();
        out.push(MoleculeEntity::Bulk(BulkEntity::new(
            id,
            mol_type,
            self.atoms,
            self.residue_name.unwrap_or(default_resname),
            self.residue_count,
        )));
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn validate_coords(row: &AtomRow) -> Result<(), BuildError> {
    for (axis, value) in [('x', row.x), ('y', row.y), ('z', row.z)] {
        if !value.is_finite() {
            return Err(BuildError::InvalidCoordinate {
                axis,
                value,
                label_atom_id: trim_atom_name(row.label_atom_id),
            });
        }
    }
    Ok(())
}

fn trim_atom_name(name: [u8; 4]) -> String {
    std::str::from_utf8(&name)
        .unwrap_or("")
        .trim_matches(|c: char| c == ' ' || c == '\0')
        .to_owned()
}

fn trim_res_name(name: &[u8; 3]) -> &str {
    std::str::from_utf8(name).unwrap_or("").trim()
}

/// True if the residue has PDB-style N, CA, and C atoms — the trigger
/// for modified-residue merging into a Protein chain. Matches the
/// existing `bridge::residue_has_backbone` check byte-for-byte.
fn residue_has_protein_backbone(residue: &ResidueAccum) -> bool {
    let mut has_n = false;
    let mut has_ca = false;
    let mut has_c = false;
    for name in residue.atoms.keys() {
        match name {
            [b' ', b'N', b' ', b' '] | [b'N', b' ', b' ', b' '] => has_n = true,
            [b' ', b'C', b'A', b' '] | [b'C', b'A', b' ', b' '] => {
                has_ca = true;
            }
            [b' ', b'C', b' ', b' '] | [b'C', b' ', b' ', b' '] => has_c = true,
            _ => {}
        }
    }
    has_n && has_ca && has_c
}

fn residue_to_atoms(r: &ResidueAccum) -> Vec<Atom> {
    r.atom_order
        .iter()
        .map(|name| r.atoms[name].to_atom())
        .collect()
}

fn flatten_residues<'a>(
    residues: impl IntoIterator<Item = &'a ResidueAccum>,
) -> (Vec<Atom>, Vec<Residue>) {
    let mut atoms: Vec<Atom> = Vec::new();
    let mut out_residues: Vec<Residue> = Vec::new();
    for r in residues {
        let start = atoms.len();
        for name in &r.atom_order {
            atoms.push(r.atoms[name].to_atom());
        }
        let end = atoms.len();
        out_residues.push(Residue {
            name: r.label_comp_id,
            number: r.label_seq_id,
            atom_range: start..end,
        });
    }
    (atoms, out_residues)
}

// ---------------------------------------------------------------------------
// Classification dispatch
// ---------------------------------------------------------------------------

fn classify_chain(
    hint: ExpectedEntityType,
    pdb_chain_id: u8,
    residues: &[ResidueAccum],
    ctx: &mut ChainCtx,
) {
    match hint {
        ExpectedEntityType::Protein => {
            emit_polymer_chain(
                residues,
                MoleculeType::Protein,
                pdb_chain_id,
                ctx,
            );
        }
        ExpectedEntityType::DNA => {
            emit_polymer_chain(residues, MoleculeType::DNA, pdb_chain_id, ctx);
        }
        ExpectedEntityType::RNA => {
            emit_polymer_chain(residues, MoleculeType::RNA, pdb_chain_id, ctx);
        }
        ExpectedEntityType::Water => emit_chain_bulk(
            residues,
            MoleculeType::Water,
            DEFAULT_WATER_RESNAME,
            ctx,
        ),
        ExpectedEntityType::NonPolymer => {
            emit_non_polymer_chain(residues, ctx);
        }
        ExpectedEntityType::Unknown => {
            emit_unknown_chain(pdb_chain_id, residues, ctx);
        }
    }
}

fn emit_polymer_chain(
    residues: &[ResidueAccum],
    mol_type: MoleculeType,
    pdb_chain_id: u8,
    ctx: &mut ChainCtx,
) {
    if residues.is_empty() {
        return;
    }
    let (atoms, res_vec) = flatten_residues(residues);
    let id = ctx.allocator.allocate();
    match mol_type {
        MoleculeType::Protein => ctx.out.push(MoleculeEntity::Protein(
            ProteinEntity::new(id, atoms, res_vec, pdb_chain_id),
        )),
        MoleculeType::DNA | MoleculeType::RNA => {
            ctx.out.push(MoleculeEntity::NucleicAcid(NAEntity::new(
                id,
                mol_type,
                atoms,
                res_vec,
                pdb_chain_id,
            )));
        }
        _ => unreachable!("emit_polymer_chain called with non-polymer type"),
    }
}

fn emit_chain_bulk(
    residues: &[ResidueAccum],
    mol_type: MoleculeType,
    default_resname: [u8; 3],
    ctx: &mut ChainCtx,
) {
    if residues.is_empty() {
        return;
    }
    let mut atoms: Vec<Atom> = Vec::new();
    for r in residues {
        for name in &r.atom_order {
            atoms.push(r.atoms[name].to_atom());
        }
    }
    let residue_name = residues
        .first()
        .map_or(default_resname, |r| r.label_comp_id);
    let id = ctx.allocator.allocate();
    ctx.out.push(MoleculeEntity::Bulk(BulkEntity::new(
        id,
        mol_type,
        atoms,
        residue_name,
        residues.len(),
    )));
}

fn emit_single_residue_small_molecule(
    r: &ResidueAccum,
    mol_type: MoleculeType,
    ctx: &mut ChainCtx,
) {
    let atoms = residue_to_atoms(r);
    let id = ctx.allocator.allocate();
    ctx.out
        .push(MoleculeEntity::SmallMolecule(SmallMoleculeEntity::new(
            id,
            mol_type,
            atoms,
            r.label_comp_id,
        )));
}

fn emit_non_polymer_chain(residues: &[ResidueAccum], ctx: &mut ChainCtx) {
    for r in residues {
        let mol_type = classify_residue(trim_res_name(&r.label_comp_id));
        match mol_type {
            MoleculeType::Water => ctx.water.ingest(r),
            MoleculeType::Solvent => ctx.solvent.ingest(r),
            MoleculeType::Protein | MoleculeType::DNA | MoleculeType::RNA => {
                // Hint says non-polymer; backbone-bearing residues fall
                // back to Ligand rather than building a polymer.
                emit_single_residue_small_molecule(
                    r,
                    MoleculeType::Ligand,
                    ctx,
                );
            }
            MoleculeType::Ligand
            | MoleculeType::Ion
            | MoleculeType::Cofactor
            | MoleculeType::Lipid => {
                emit_single_residue_small_molecule(r, mol_type, ctx);
            }
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
#[allow(
    clippy::upper_case_acronyms,
    reason = "DNA / RNA mirror MoleculeType's variant names"
)]
enum UnknownBucket {
    Protein,
    DNA,
    RNA,
    Water,
    Solvent,
    Small(MoleculeType),
}

fn assign_unknown_bucket(r: &ResidueAccum, has_protein: bool) -> UnknownBucket {
    let mol_type = classify_residue(trim_res_name(&r.label_comp_id));
    match mol_type {
        MoleculeType::Protein => UnknownBucket::Protein,
        MoleculeType::DNA => UnknownBucket::DNA,
        MoleculeType::RNA => UnknownBucket::RNA,
        MoleculeType::Water => UnknownBucket::Water,
        MoleculeType::Solvent => UnknownBucket::Solvent,
        MoleculeType::Ligand
        | MoleculeType::Ion
        | MoleculeType::Cofactor
        | MoleculeType::Lipid => {
            if has_protein && residue_has_protein_backbone(r) {
                UnknownBucket::Protein
            } else {
                UnknownBucket::Small(mol_type)
            }
        }
    }
}

/// Full residue-name heuristic with the protein-merge logic preserved
/// from `bridge::split_into_entities`: a residue that classifies as
/// Ligand/Ion/Cofactor/Lipid but carries protein backbone atoms gets
/// folded into the chain's protein bucket when one exists.
fn emit_unknown_chain(
    pdb_chain_id: u8,
    residues: &[ResidueAccum],
    ctx: &mut ChainCtx,
) {
    let has_protein = residues.iter().any(|r| {
        classify_residue(trim_res_name(&r.label_comp_id))
            == MoleculeType::Protein
    });
    let buckets: Vec<UnknownBucket> = residues
        .iter()
        .map(|r| assign_unknown_bucket(r, has_protein))
        .collect();

    emit_unknown_polymer(
        UnknownBucket::Protein,
        residues,
        &buckets,
        pdb_chain_id,
        ctx,
    );
    emit_unknown_polymer(
        UnknownBucket::DNA,
        residues,
        &buckets,
        pdb_chain_id,
        ctx,
    );
    emit_unknown_polymer(
        UnknownBucket::RNA,
        residues,
        &buckets,
        pdb_chain_id,
        ctx,
    );

    for (r, bucket) in residues.iter().zip(buckets.iter()) {
        match bucket {
            UnknownBucket::Water => ctx.water.ingest(r),
            UnknownBucket::Solvent => ctx.solvent.ingest(r),
            UnknownBucket::Small(mt) => {
                emit_single_residue_small_molecule(r, *mt, ctx);
            }
            UnknownBucket::Protein
            | UnknownBucket::DNA
            | UnknownBucket::RNA => {}
        }
    }
}

fn emit_unknown_polymer(
    target: UnknownBucket,
    residues: &[ResidueAccum],
    buckets: &[UnknownBucket],
    pdb_chain_id: u8,
    ctx: &mut ChainCtx,
) {
    let selected: Vec<&ResidueAccum> = residues
        .iter()
        .zip(buckets.iter())
        .filter(|(_, b)| **b == target)
        .map(|(r, _)| r)
        .collect();
    if selected.is_empty() {
        return;
    }
    let (atoms, res_vec) = flatten_residues(selected.iter().copied());
    let id = ctx.allocator.allocate();
    match target {
        UnknownBucket::Protein => ctx.out.push(MoleculeEntity::Protein(
            ProteinEntity::new(id, atoms, res_vec, pdb_chain_id),
        )),
        UnknownBucket::DNA => ctx.out.push(MoleculeEntity::NucleicAcid(
            NAEntity::new(id, MoleculeType::DNA, atoms, res_vec, pdb_chain_id),
        )),
        UnknownBucket::RNA => ctx.out.push(MoleculeEntity::NucleicAcid(
            NAEntity::new(id, MoleculeType::RNA, atoms, res_vec, pdb_chain_id),
        )),
        UnknownBucket::Water
        | UnknownBucket::Solvent
        | UnknownBucket::Small(_) => {
            unreachable!("emit_unknown_polymer called with non-polymer bucket")
        }
    }
}

#[cfg(test)]
#[path = "builder_tests.rs"]
mod tests;

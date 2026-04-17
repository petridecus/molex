//! Top-level `Assembly` container: entities, cross-entity bonds, and
//! eagerly-recomputed derived data (secondary structure + backbone H-bonds).
//!
//! `Assembly` is the host-owned structural source of truth. Every
//! `&mut Assembly` mutation bumps the generation counter and recomputes
//! all derived data before returning, so readers of a given snapshot see
//! consistent `ss_types`, `hbonds`, and `cross_entity_bonds`.

use std::collections::HashMap;

use glam::Vec3;

use crate::analysis::bonds::disulfide::detect_disulfides;
#[allow(
    deprecated,
    reason = "Assembly is the replacement path; the free function stays as an \
              implementation detail through Phase 4 and gets deleted in Phase \
              5."
)]
use crate::analysis::bonds::hydrogen::detect_hbonds;
use crate::analysis::bonds::hydrogen::HBond;
use crate::analysis::ss::dssp::classify;
use crate::analysis::SSType;
use crate::atom_id::AtomId;
use crate::bond::CovalentBond;
use crate::entity::molecule::id::EntityId;
use crate::entity::molecule::MoleculeEntity;

/// Top-level container of entities plus eagerly-computed derived data.
///
/// Layout matches `Assembly migration` decisions #2, #5, #7, #11: entities
/// own their intra-entity bonds; `Assembly` owns cross-entity bonds
/// (disulfides only in this migration) and the full derived-data set
/// (`ss_types`, `hbonds`). The generation counter increments on every
/// mutation so consumers can detect snapshots cheaply.
#[derive(Debug, Clone)]
pub struct Assembly {
    entities: Vec<MoleculeEntity>,
    cross_entity_bonds: Vec<CovalentBond>,
    ss_types: HashMap<EntityId, Vec<SSType>>,
    hbonds: Vec<HBond>,
    generation: u64,
}

/// Snapshot of per-entity atom positions for whole-assembly replacement.
///
/// Entities present in the snapshot but missing from the target
/// `Assembly` are ignored; entities in the target but missing from the
/// snapshot retain their current positions.
#[derive(Debug, Clone, Default)]
pub struct CoordinateSnapshot {
    per_entity: HashMap<EntityId, Vec<Vec3>>,
}

impl CoordinateSnapshot {
    /// Create a snapshot from a pre-built per-entity map.
    #[must_use]
    pub fn new(per_entity: HashMap<EntityId, Vec<Vec3>>) -> Self {
        Self { per_entity }
    }

    /// Capture the current positions of every entity in an `Assembly`.
    ///
    /// Intended for "reset to original" flows: take a snapshot after
    /// load, run mutations, then pass the snapshot back to
    /// [`Assembly::set_coordinate_snapshot`] to restore.
    #[must_use]
    pub fn from_assembly(assembly: &Assembly) -> Self {
        let per_entity = assembly
            .entities
            .iter()
            .map(|e| (e.id(), e.positions()))
            .collect();
        Self { per_entity }
    }

    /// Positions for a specific entity, if present.
    #[must_use]
    pub fn positions(&self, id: EntityId) -> Option<&[Vec3]> {
        self.per_entity.get(&id).map(Vec::as_slice)
    }

    /// Number of entities covered by this snapshot.
    #[must_use]
    pub fn len(&self) -> usize {
        self.per_entity.len()
    }

    /// Whether the snapshot covers no entities.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.per_entity.is_empty()
    }
}

impl Assembly {
    /// Build an `Assembly` from a collection of entities.
    ///
    /// Runs disulfide detection, per-entity DSSP classification, and
    /// flat-backbone H-bond detection to populate all derived data.
    /// `generation` is initialized to 0.
    #[must_use]
    pub fn new(entities: Vec<MoleculeEntity>) -> Self {
        let mut this = Self {
            entities,
            cross_entity_bonds: Vec::new(),
            ss_types: HashMap::new(),
            hbonds: Vec::new(),
            generation: 0,
        };
        this.recompute_derived();
        this
    }

    // ── Read accessors ──────────────────────────────────────────────

    /// All entities in declaration order.
    #[must_use]
    pub fn entities(&self) -> &[MoleculeEntity] {
        &self.entities
    }

    /// Look up an entity by id.
    #[must_use]
    pub fn entity(&self, id: EntityId) -> Option<&MoleculeEntity> {
        self.entities.iter().find(|e| e.id() == id)
    }

    /// Monotonic counter incremented on every mutation.
    #[must_use]
    pub fn generation(&self) -> u64 {
        self.generation
    }

    /// Backbone hydrogen bonds (Kabsch-Sander) across all protein
    /// entities. Donor and acceptor indices refer to the flattened
    /// per-protein-entity backbone sequence produced by concatenating
    /// `ProteinEntity::to_backbone()` in entity order.
    #[must_use]
    pub fn hbonds(&self) -> &[HBond] {
        &self.hbonds
    }

    /// Secondary structure classification for an entity.
    ///
    /// Returns an empty slice for entities that don't have an SS
    /// assignment (non-protein, or proteins with fewer than two complete
    /// backbone residues).
    #[must_use]
    pub fn ss_types(&self, id: EntityId) -> &[SSType] {
        self.ss_types.get(&id).map_or(&[], Vec::as_slice)
    }

    /// Cross-entity covalent bonds. In this migration this contains
    /// disulfides only (decision #12); same-chain SG-SG bridges are
    /// also collected here rather than on the owning `ProteinEntity`.
    #[must_use]
    pub fn cross_entity_bonds(&self) -> &[CovalentBond] {
        &self.cross_entity_bonds
    }

    /// All atoms bonded to `atom`, yielding the far endpoint of each
    /// matching bond. Walks the owning entity's intra-entity bond list
    /// and the assembly's `cross_entity_bonds`.
    pub fn bonds_touching(
        &self,
        atom: AtomId,
    ) -> impl Iterator<Item = AtomId> + '_ {
        let intra = self
            .entity(atom.entity)
            .and_then(entity_bonds)
            .unwrap_or(&[])
            .iter()
            .filter_map(move |b| other_endpoint(b, atom));
        let cross = self
            .cross_entity_bonds
            .iter()
            .filter_map(move |b| other_endpoint(b, atom));
        intra.chain(cross)
    }

    /// Cross-entity disulfide bonds — pairs where both endpoints
    /// resolve to an SG atom inside a CYS residue.
    pub fn disulfides(&self) -> impl Iterator<Item = &CovalentBond> + '_ {
        self.cross_entity_bonds
            .iter()
            .filter(|b| self.is_cys_sg(b.a) && self.is_cys_sg(b.b))
    }

    // ── Mutation methods ────────────────────────────────────────────

    /// Append an entity. Bumps the generation counter and recomputes
    /// all derived data.
    pub fn add_entity(&mut self, entity: MoleculeEntity) {
        self.entities.push(entity);
        self.after_mutation();
    }

    /// Remove an entity by id. Any `cross_entity_bonds` touching the
    /// removed entity are purged as part of the derived-data
    /// recomputation.
    pub fn remove_entity(&mut self, id: EntityId) {
        self.entities.retain(|e| e.id() != id);
        let _ = self.ss_types.remove(&id);
        self.after_mutation();
    }

    /// Replace the positions of a single entity. If `coords.len()` does
    /// not match the entity's atom count, the mutation is abandoned
    /// (logged at error level) and the generation counter is not
    /// advanced — readers keep seeing the previous snapshot.
    pub fn update_positions(&mut self, entity: EntityId, coords: &[Vec3]) {
        let Some(idx) = self.entities.iter().position(|e| e.id() == entity)
        else {
            log::error!(
                "Assembly::update_positions: unknown entity id {entity}"
            );
            return;
        };
        let atoms = self.entities[idx].atom_set_mut();
        if atoms.len() != coords.len() {
            log::error!(
                "Assembly::update_positions: entity {entity} has {} atoms but \
                 {} coords were supplied; skipping update",
                atoms.len(),
                coords.len(),
            );
            return;
        }
        for (atom, &pos) in atoms.iter_mut().zip(coords.iter()) {
            atom.position = pos;
        }
        self.after_mutation();
    }

    /// Apply positions for every entity covered by `snapshot`.
    ///
    /// Entities in the snapshot whose atom counts don't match their
    /// target entity are skipped (logged at error level) while other
    /// entities still get their positions applied. Entities absent from
    /// the snapshot retain their current positions. Bumps the
    /// generation counter and recomputes derived data exactly once,
    /// regardless of how many entities were touched.
    #[allow(
        clippy::needless_pass_by_value,
        reason = "phase_3.md locks this signature — `snapshot` is the \
                  one-shot \"apply this whole-assembly snapshot\" input; \
                  callers hand over ownership rather than keeping the \
                  snapshot alive."
    )]
    pub fn set_coordinate_snapshot(&mut self, snapshot: CoordinateSnapshot) {
        for entity in &mut self.entities {
            let entity_id = entity.id();
            let Some(coords) = snapshot.per_entity.get(&entity_id) else {
                continue;
            };
            let atoms = entity.atom_set_mut();
            if atoms.len() != coords.len() {
                log::error!(
                    "Assembly::set_coordinate_snapshot: entity {entity_id} \
                     has {} atoms but snapshot carries {}; skipping this \
                     entity",
                    atoms.len(),
                    coords.len(),
                );
                continue;
            }
            for (atom, &pos) in atoms.iter_mut().zip(coords.iter()) {
                atom.position = pos;
            }
        }
        self.after_mutation();
    }

    // ── Internal helpers ────────────────────────────────────────────

    fn after_mutation(&mut self) {
        self.generation = self.generation.saturating_add(1);
        self.recompute_derived();
    }

    fn recompute_derived(&mut self) {
        self.cross_entity_bonds = detect_disulfides(&self.entities);
        self.ss_types = compute_per_entity_ss(&self.entities);
        self.hbonds = compute_flat_hbonds(&self.entities);
    }

    fn is_cys_sg(&self, atom: AtomId) -> bool {
        let Some(entity) = self.entity(atom.entity) else {
            return false;
        };
        let Some(protein) = entity.as_protein() else {
            return false;
        };
        let idx = atom.index as usize;
        let Some(residue) = protein
            .residues
            .iter()
            .find(|r| r.atom_range.contains(&idx))
        else {
            return false;
        };
        if trimmed(&residue.name) != b"CYS" {
            return false;
        }
        let Some(atom) = protein.atoms.get(idx) else {
            return false;
        };
        trimmed_atom_name(&atom.name) == b"SG"
    }
}

fn entity_bonds(entity: &MoleculeEntity) -> Option<&[CovalentBond]> {
    match entity {
        MoleculeEntity::Protein(e) => Some(&e.bonds),
        MoleculeEntity::NucleicAcid(e) => Some(&e.bonds),
        MoleculeEntity::SmallMolecule(e) => Some(&e.bonds),
        MoleculeEntity::Bulk(_) => None,
    }
}

fn other_endpoint(bond: &CovalentBond, atom: AtomId) -> Option<AtomId> {
    if bond.a == atom {
        Some(bond.b)
    } else if bond.b == atom {
        Some(bond.a)
    } else {
        None
    }
}

#[allow(
    deprecated,
    reason = "Assembly owns the per-entity DSSP recompute — it wraps the \
              legacy detect_hbonds, which stays as the implementation detail \
              until the final Phase 5 deletion."
)]
fn compute_per_entity_ss(
    entities: &[MoleculeEntity],
) -> HashMap<EntityId, Vec<SSType>> {
    let mut out = HashMap::new();
    for entity in entities {
        let Some(protein) = entity.as_protein() else {
            continue;
        };
        let backbone = protein.to_backbone();
        if backbone.is_empty() {
            continue;
        }
        let hbonds = detect_hbonds(&backbone);
        let ss = classify(&hbonds, backbone.len());
        let _ = out.insert(protein.id, ss);
    }
    out
}

#[allow(
    deprecated,
    reason = "Assembly owns the flat-backbone H-bond recompute — it wraps the \
              legacy detect_hbonds, which stays as the implementation detail \
              until the final Phase 5 deletion."
)]
fn compute_flat_hbonds(entities: &[MoleculeEntity]) -> Vec<HBond> {
    let mut flat = Vec::new();
    for entity in entities {
        let Some(protein) = entity.as_protein() else {
            continue;
        };
        flat.extend(protein.to_backbone());
    }
    detect_hbonds(&flat)
}

fn trimmed(name: &[u8; 3]) -> &[u8] {
    let mut end = name.len();
    while end > 0 && (name[end - 1] == b' ' || name[end - 1] == 0) {
        end -= 1;
    }
    &name[..end]
}

fn trimmed_atom_name(name: &[u8; 4]) -> &[u8] {
    let mut end = name.len();
    while end > 0 && (name[end - 1] == b' ' || name[end - 1] == 0) {
        end -= 1;
    }
    let mut start = 0;
    while start < end && (name[start] == b' ' || name[start] == 0) {
        start += 1;
    }
    &name[start..end]
}

#[cfg(test)]
#[allow(
    clippy::unwrap_used,
    clippy::float_cmp,
    clippy::cast_precision_loss,
    deprecated
)]
mod tests {
    use std::path::Path;

    use glam::Vec3;

    use super::*;
    use crate::adapters::pdb::structure_file_to_entities;
    use crate::analysis::bonds::disulfide::detect_disulfides;
    use crate::analysis::bonds::hydrogen::detect_hbonds;
    use crate::analysis::ss::dssp::classify;
    use crate::element::Element;
    use crate::entity::molecule::atom::Atom;
    use crate::entity::molecule::id::EntityIdAllocator;
    use crate::entity::molecule::protein::ProteinEntity;
    use crate::entity::molecule::Residue;

    fn mk_atom(name: [u8; 4], el: Element, pos: Vec3) -> Atom {
        Atom {
            position: pos,
            occupancy: 1.0,
            b_factor: 0.0,
            element: el,
            name,
        }
    }

    /// Build a minimal two-residue protein (ALA-GLY) with backbone +
    /// one sidechain heavy atom, laid out in canonical order after
    /// `ProteinEntity::new`.
    fn make_dipeptide(
        alloc: &mut EntityIdAllocator,
        chain: u8,
        origin: Vec3,
    ) -> MoleculeEntity {
        let id = alloc.allocate();
        let atoms = vec![
            mk_atom(*b"N   ", Element::N, origin),
            mk_atom(*b"CA  ", Element::C, origin + Vec3::new(1.0, 0.0, 0.0)),
            mk_atom(*b"C   ", Element::C, origin + Vec3::new(2.0, 0.0, 0.0)),
            mk_atom(*b"O   ", Element::O, origin + Vec3::new(2.0, 1.0, 0.0)),
            mk_atom(*b"CB  ", Element::C, origin + Vec3::new(1.0, -1.0, 0.0)),
            mk_atom(*b"N   ", Element::N, origin + Vec3::new(3.2, 0.0, 0.0)),
            mk_atom(*b"CA  ", Element::C, origin + Vec3::new(4.2, 0.0, 0.0)),
            mk_atom(*b"C   ", Element::C, origin + Vec3::new(5.2, 0.0, 0.0)),
            mk_atom(*b"O   ", Element::O, origin + Vec3::new(5.2, 1.0, 0.0)),
        ];
        let residues = vec![
            Residue {
                name: *b"ALA",
                number: 1,
                atom_range: 0..5,
            },
            Residue {
                name: *b"GLY",
                number: 2,
                atom_range: 5..9,
            },
        ];
        MoleculeEntity::Protein(ProteinEntity::new(id, atoms, residues, chain))
    }

    /// A single cysteine residue with a bondable SG at a given position.
    fn cys_residue_with_sg(
        alloc: &mut EntityIdAllocator,
        chain: u8,
        sg_pos: Vec3,
    ) -> MoleculeEntity {
        let id = alloc.allocate();
        let atoms = vec![
            mk_atom(*b"N   ", Element::N, Vec3::new(0.0, 0.0, 0.0)),
            mk_atom(*b"CA  ", Element::C, Vec3::new(1.0, 0.0, 0.0)),
            mk_atom(*b"C   ", Element::C, Vec3::new(2.0, 0.0, 0.0)),
            mk_atom(*b"O   ", Element::O, Vec3::new(2.0, 1.0, 0.0)),
            mk_atom(*b"CB  ", Element::C, Vec3::new(1.0, -1.0, 0.0)),
            mk_atom(*b"SG  ", Element::S, sg_pos),
        ];
        let residues = vec![Residue {
            name: *b"CYS",
            number: 1,
            atom_range: 0..atoms.len(),
        }];
        MoleculeEntity::Protein(ProteinEntity::new(id, atoms, residues, chain))
    }

    // -- Construction + generation --

    #[test]
    fn new_starts_at_generation_zero() {
        let mut alloc = EntityIdAllocator::new();
        let dipep = make_dipeptide(&mut alloc, b'A', Vec3::ZERO);
        let assembly = Assembly::new(vec![dipep]);
        assert_eq!(assembly.generation(), 0);
    }

    #[test]
    fn new_exposes_all_entities() {
        let mut alloc = EntityIdAllocator::new();
        let a = make_dipeptide(&mut alloc, b'A', Vec3::ZERO);
        let b = make_dipeptide(&mut alloc, b'B', Vec3::new(20.0, 0.0, 0.0));
        let assembly = Assembly::new(vec![a, b]);
        assert_eq!(assembly.entities().len(), 2);
    }

    // -- Mutation generation + recompute --

    #[test]
    fn add_entity_bumps_generation_exactly_once() {
        let mut alloc = EntityIdAllocator::new();
        let mut assembly =
            Assembly::new(vec![make_dipeptide(&mut alloc, b'A', Vec3::ZERO)]);
        let before = assembly.generation();
        assembly.add_entity(make_dipeptide(
            &mut alloc,
            b'B',
            Vec3::new(20.0, 0.0, 0.0),
        ));
        assert_eq!(assembly.generation(), before + 1);
    }

    #[test]
    fn remove_entity_bumps_generation_exactly_once() {
        let mut alloc = EntityIdAllocator::new();
        let a = make_dipeptide(&mut alloc, b'A', Vec3::ZERO);
        let b_id;
        let b = {
            let e = make_dipeptide(&mut alloc, b'B', Vec3::new(20.0, 0.0, 0.0));
            b_id = e.id();
            e
        };
        let mut assembly = Assembly::new(vec![a, b]);
        let before = assembly.generation();
        assembly.remove_entity(b_id);
        assert_eq!(assembly.generation(), before + 1);
        assert_eq!(assembly.entities().len(), 1);
    }

    #[test]
    fn update_positions_bumps_generation_and_moves_atoms() {
        let mut alloc = EntityIdAllocator::new();
        let entity = make_dipeptide(&mut alloc, b'A', Vec3::ZERO);
        let id = entity.id();
        let n_atoms = entity.atom_count();
        let mut assembly = Assembly::new(vec![entity]);
        let before = assembly.generation();

        let shifted: Vec<Vec3> = (0..n_atoms)
            .map(|i| Vec3::new(i as f32, 100.0, 0.0))
            .collect();
        assembly.update_positions(id, &shifted);

        assert_eq!(assembly.generation(), before + 1);
        let updated = assembly.entity(id).unwrap();
        assert!(
            (updated.atom_set()[0].position - Vec3::new(0.0, 100.0, 0.0))
                .length()
                < 1e-6
        );
    }

    #[test]
    fn update_positions_length_mismatch_is_a_no_op() {
        let mut alloc = EntityIdAllocator::new();
        let entity = make_dipeptide(&mut alloc, b'A', Vec3::ZERO);
        let id = entity.id();
        let mut assembly = Assembly::new(vec![entity]);
        let before = assembly.generation();

        // Too few coords — must be rejected without bumping generation.
        assembly.update_positions(id, &[Vec3::ZERO]);
        assert_eq!(assembly.generation(), before);
    }

    #[test]
    fn mutation_result_matches_fresh_build() {
        // A fresh Assembly over the same final entity set must produce
        // the same derived outputs as a mutated Assembly: that is the
        // contract "derived data is recomputed on every mutation".
        let mut alloc = EntityIdAllocator::new();
        let a = make_dipeptide(&mut alloc, b'A', Vec3::ZERO);
        let b = make_dipeptide(&mut alloc, b'B', Vec3::new(20.0, 0.0, 0.0));
        let a_id = a.id();

        // Mutation path: start with just `a`, add `b`.
        let mut mutated = Assembly::new(vec![a.clone()]);
        mutated.add_entity(b.clone());

        // Fresh path: build from both at once.
        let fresh = Assembly::new(vec![a, b]);

        assert_eq!(mutated.hbonds().len(), fresh.hbonds().len());
        assert_eq!(
            mutated.cross_entity_bonds().len(),
            fresh.cross_entity_bonds().len(),
        );
        assert_eq!(mutated.ss_types(a_id).len(), fresh.ss_types(a_id).len());
    }

    // -- Disulfide handling --

    #[test]
    fn disulfides_filtered_and_cross_bonds_populated() {
        let mut alloc = EntityIdAllocator::new();
        // Two cysteine entities with SGs at ~2.03 Å separation.
        let sg_a = Vec3::new(0.0, 0.0, 0.0);
        let sg_b = Vec3::new(2.03, 0.0, 0.0);
        let ca = cys_residue_with_sg(&mut alloc, b'A', sg_a);
        let cb = cys_residue_with_sg(&mut alloc, b'B', sg_b);

        let assembly = Assembly::new(vec![ca, cb]);
        assert_eq!(assembly.cross_entity_bonds().len(), 1);
        assert_eq!(assembly.disulfides().count(), 1);
    }

    #[test]
    fn remove_entity_purges_touching_cross_bonds() {
        let mut alloc = EntityIdAllocator::new();
        let sg_a = Vec3::new(0.0, 0.0, 0.0);
        let sg_b = Vec3::new(2.03, 0.0, 0.0);
        let ca = cys_residue_with_sg(&mut alloc, b'A', sg_a);
        let cb = cys_residue_with_sg(&mut alloc, b'B', sg_b);
        let cb_id = cb.id();

        let mut assembly = Assembly::new(vec![ca, cb]);
        assert_eq!(assembly.cross_entity_bonds().len(), 1);

        assembly.remove_entity(cb_id);
        assert!(assembly.cross_entity_bonds().is_empty());
        assert_eq!(assembly.disulfides().count(), 0);
    }

    // -- bonds_touching --

    #[test]
    fn bonds_touching_walks_both_intra_and_cross() {
        let mut alloc = EntityIdAllocator::new();
        let sg_a = Vec3::new(0.0, 0.0, 0.0);
        let sg_b = Vec3::new(2.03, 0.0, 0.0);
        let ca = cys_residue_with_sg(&mut alloc, b'A', sg_a);
        let ca_id = ca.id();
        let cb = cys_residue_with_sg(&mut alloc, b'B', sg_b);
        let assembly = Assembly::new(vec![ca, cb]);

        // SG in chain A sits at index 5 after canonical ordering
        // (N, CA, C, O, CB, SG).
        let sg_atom = AtomId {
            entity: ca_id,
            index: 5,
        };
        let neighbors: Vec<AtomId> = assembly.bonds_touching(sg_atom).collect();

        // Must include the CB neighbor (intra) and the SG partner on
        // chain B (cross). That is at least two.
        assert!(
            neighbors.len() >= 2,
            "expected >=2 neighbors for SG, got {}: {neighbors:?}",
            neighbors.len()
        );
    }

    // -- CoordinateSnapshot roundtrip --

    #[test]
    fn coordinate_snapshot_roundtrip_restores_positions() {
        let mut alloc = EntityIdAllocator::new();
        let entity = make_dipeptide(&mut alloc, b'A', Vec3::ZERO);
        let id = entity.id();
        let mut assembly = Assembly::new(vec![entity]);
        let original = CoordinateSnapshot::from_assembly(&assembly);

        // Shift every atom by +10 in z via set_coordinate_snapshot.
        let shifted_positions: Vec<Vec3> = assembly
            .entity(id)
            .unwrap()
            .atom_set()
            .iter()
            .map(|a| a.position + Vec3::new(0.0, 0.0, 10.0))
            .collect();
        let mut per_entity = HashMap::new();
        let _ = per_entity.insert(id, shifted_positions);
        assembly.set_coordinate_snapshot(CoordinateSnapshot::new(per_entity));
        assert!(
            (assembly.entity(id).unwrap().atom_set()[0].position.z - 10.0)
                .abs()
                < 1e-6
        );

        // Restore and confirm we're back.
        assembly.set_coordinate_snapshot(original);
        assert!(
            assembly.entity(id).unwrap().atom_set()[0].position.z.abs() < 1e-6
        );
    }

    #[test]
    fn coordinate_snapshot_skips_length_mismatch_but_applies_matching() {
        let mut alloc = EntityIdAllocator::new();
        let a = make_dipeptide(&mut alloc, b'A', Vec3::ZERO);
        let b = make_dipeptide(&mut alloc, b'B', Vec3::new(20.0, 0.0, 0.0));
        let a_id = a.id();
        let b_id = b.id();
        let b_original_pos = b.atom_set()[0].position;
        let mut assembly = Assembly::new(vec![a, b]);
        let gen_before = assembly.generation();

        let a_count = assembly.entity(a_id).unwrap().atom_count();
        let shifted_a: Vec<Vec3> =
            (0..a_count).map(|i| Vec3::new(i as f32, 50.0, 0.0)).collect();
        let too_short_b = vec![Vec3::ZERO];

        let mut per_entity = HashMap::new();
        let _ = per_entity.insert(a_id, shifted_a);
        let _ = per_entity.insert(b_id, too_short_b);
        assembly.set_coordinate_snapshot(CoordinateSnapshot::new(per_entity));

        assert_eq!(assembly.generation(), gen_before + 1);
        assert!(
            (assembly.entity(a_id).unwrap().atom_set()[0].position
                - Vec3::new(0.0, 50.0, 0.0))
            .length()
                < 1e-6
        );
        assert!(
            (assembly.entity(b_id).unwrap().atom_set()[0].position
                - b_original_pos)
                .length()
                < 1e-6
        );
    }

    #[test]
    fn coordinate_snapshot_with_unknown_ids_still_bumps_generation() {
        let mut alloc = EntityIdAllocator::new();
        let a = make_dipeptide(&mut alloc, b'A', Vec3::ZERO);
        let a_id = a.id();
        let original_pos = a.atom_set()[0].position;
        let mut assembly = Assembly::new(vec![a]);
        let gen_before = assembly.generation();

        let mut other_alloc = EntityIdAllocator::new();
        let unknown_id = other_alloc.from_raw(10_000);

        let mut per_entity = HashMap::new();
        let _ = per_entity.insert(unknown_id, vec![Vec3::ZERO; 9]);
        assembly.set_coordinate_snapshot(CoordinateSnapshot::new(per_entity));

        assert_eq!(assembly.generation(), gen_before + 1);
        assert!(
            (assembly.entity(a_id).unwrap().atom_set()[0].position
                - original_pos)
                .length()
                < 1e-6
        );
    }

    // -- Byte-identity vs the legacy per-slice path (viso's Phase-3 gate) --

    /// On a real protein, `Assembly::ss_types` / `Assembly::hbonds`
    /// must match what the legacy free-function path produces on the
    /// same inputs. This is the "byte-identity" gate that lets Phase 4
    /// swap the data source without rendering regressions.
    #[test]
    fn byte_identity_against_legacy_path_1ubq() {
        let path = Path::new("../viso/assets/models/1ubq.cif");
        assert!(
            path.exists(),
            "byte-identity fixture missing at {}",
            path.display(),
        );

        let entities = structure_file_to_entities(path).unwrap();
        let protein = entities.iter().find_map(|e| e.as_protein()).unwrap();
        let protein_id = protein.id;

        // Legacy path (still used by viso pre-Phase-4):
        let legacy_backbone: Vec<_> = entities
            .iter()
            .filter_map(MoleculeEntity::as_protein)
            .flat_map(ProteinEntity::to_backbone)
            .collect();
        let legacy_hbonds = detect_hbonds(&legacy_backbone);
        let per_entity_backbone = protein.to_backbone();
        let legacy_ss = classify(
            &detect_hbonds(&per_entity_backbone),
            per_entity_backbone.len(),
        );

        // Assembly path:
        let assembly = Assembly::new(entities.clone());
        assert_eq!(
            assembly.hbonds(),
            legacy_hbonds.as_slice(),
            "flat-backbone H-bonds must match the legacy viso path"
        );
        assert_eq!(
            assembly.ss_types(protein_id),
            legacy_ss.as_slice(),
            "per-entity SS must match the legacy per-entity DSSP path"
        );
        // Disulfide coverage lives in `disulfides_filtered_and_cross_bonds_populated`;
        // asserting it here against `detect_disulfides` would be tautological
        // because `Assembly::new` calls the same function internally.
    }
}

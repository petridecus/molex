//! Entity-based extraction utilities.

use glam::Vec3;

use crate::entity::molecule::protein::ProteinEntity;
use crate::entity::molecule::MoleculeEntity;

/// Extract CA positions from protein entities.
#[must_use]
pub fn extract_ca_positions(entities: &[MoleculeEntity]) -> Vec<Vec3> {
    entities
        .iter()
        .filter_map(MoleculeEntity::as_protein)
        .flat_map(|p| p.to_backbone().into_iter().map(|bb| bb.ca))
        .collect()
}

/// Extract backbone segments (N-CA-C interleaved) from protein entities.
///
/// Returns one `Vec<Vec3>` per continuous backbone segment across all
/// protein entities.
#[must_use]
pub fn extract_backbone_segments(
    entities: &[MoleculeEntity],
) -> Vec<Vec<Vec3>> {
    entities
        .iter()
        .filter_map(MoleculeEntity::as_protein)
        .flat_map(ProteinEntity::to_interleaved_segments)
        .collect()
}

/// Extract CA positions from backbone chains (every 2nd element in N-CA-C
/// pattern).
#[must_use]
pub fn extract_ca_from_chains(chains: &[Vec<Vec3>]) -> Vec<Vec3> {
    let mut ca_positions = Vec::new();
    for chain in chains {
        for (i, pos) in chain.iter().enumerate() {
            if i % 3 == 1 {
                ca_positions.push(*pos);
            }
        }
    }
    ca_positions
}

/// Compute centroid of a point set.
#[must_use]
#[allow(clippy::cast_precision_loss, reason = "point count fits in f32")]
pub fn centroid(points: &[Vec3]) -> Vec3 {
    if points.is_empty() {
        return Vec3::ZERO;
    }
    let sum: Vec3 = points.iter().copied().sum();
    sum / points.len() as f32
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::float_cmp)]
mod tests {
    use super::*;
    use crate::element::Element;
    use crate::ops::codec::{split_into_entities, Coords, CoordsAtom};

    fn make_atom(x: f32, y: f32, z: f32) -> CoordsAtom {
        CoordsAtom {
            x,
            y,
            z,
            occupancy: 1.0,
            b_factor: 0.0,
        }
    }
    fn res_name(s: &str) -> [u8; 3] {
        let mut n = [b' '; 3];
        for (i, b) in s.bytes().take(3).enumerate() {
            n[i] = b;
        }
        n
    }
    fn atom_name(s: &str) -> [u8; 4] {
        let mut n = [b' '; 4];
        for (i, b) in s.bytes().take(4).enumerate() {
            n[i] = b;
        }
        n
    }

    fn make_two_residue_protein_entities() -> Vec<MoleculeEntity> {
        let coords = Coords {
            num_atoms: 8,
            atoms: vec![
                make_atom(1.0, 2.0, 3.0), // N
                make_atom(4.0, 5.0, 6.0), // CA
                make_atom(5.3, 5.0, 6.0), /* C (close to next N for no
                                           * break) */
                make_atom(5.3, 6.0, 6.0),    // O
                make_atom(6.5, 5.0, 6.0),    // N (1.2A from prev C)
                make_atom(16.0, 17.0, 18.0), // CA
                make_atom(19.0, 20.0, 21.0), // C
                make_atom(22.0, 23.0, 24.0), // O
            ],
            chain_ids: vec![b'A'; 8],
            res_names: vec![
                res_name("ALA"),
                res_name("ALA"),
                res_name("ALA"),
                res_name("ALA"),
                res_name("GLY"),
                res_name("GLY"),
                res_name("GLY"),
                res_name("GLY"),
            ],
            res_nums: vec![1, 1, 1, 1, 2, 2, 2, 2],
            atom_names: vec![
                atom_name("N"),
                atom_name("CA"),
                atom_name("C"),
                atom_name("O"),
                atom_name("N"),
                atom_name("CA"),
                atom_name("C"),
                atom_name("O"),
            ],
            elements: vec![
                Element::N,
                Element::C,
                Element::C,
                Element::O,
                Element::N,
                Element::C,
                Element::C,
                Element::O,
            ],
        };
        split_into_entities(&coords)
    }

    // -- extract_ca_positions --

    #[test]
    fn extract_ca_positions_from_protein() {
        let entities = make_two_residue_protein_entities();
        let ca_pos = extract_ca_positions(&entities);
        assert_eq!(ca_pos.len(), 2);
        assert!((ca_pos[0].x - 4.0).abs() < 1e-6);
        assert!((ca_pos[1].x - 16.0).abs() < 1e-6);
    }

    #[test]
    fn extract_ca_positions_empty_for_non_protein() {
        let coords = Coords {
            num_atoms: 1,
            atoms: vec![make_atom(1.0, 2.0, 3.0)],
            chain_ids: vec![b' '],
            res_names: vec![res_name("HOH")],
            res_nums: vec![100],
            atom_names: vec![atom_name("O")],
            elements: vec![Element::O],
        };
        let entities = split_into_entities(&coords);
        let ca_pos = extract_ca_positions(&entities);
        assert!(ca_pos.is_empty());
    }

    #[test]
    fn extract_ca_positions_empty_for_empty_entities() {
        let ca_pos = extract_ca_positions(&[]);
        assert!(ca_pos.is_empty());
    }

    // -- extract_backbone_segments --

    #[test]
    fn extract_backbone_segments_from_protein() {
        let entities = make_two_residue_protein_entities();
        let segments = extract_backbone_segments(&entities);
        assert!(!segments.is_empty());
        // Each segment should have positions that are multiples of 3 (N,CA,C
        // per residue)
        for seg in &segments {
            assert_eq!(seg.len() % 3, 0);
        }
    }

    #[test]
    fn extract_backbone_segments_empty_for_non_protein() {
        let coords = Coords {
            num_atoms: 1,
            atoms: vec![make_atom(1.0, 2.0, 3.0)],
            chain_ids: vec![b' '],
            res_names: vec![res_name("HOH")],
            res_nums: vec![100],
            atom_names: vec![atom_name("O")],
            elements: vec![Element::O],
        };
        let entities = split_into_entities(&coords);
        let segments = extract_backbone_segments(&entities);
        assert!(segments.is_empty());
    }

    // -- extract_ca_from_chains --

    #[test]
    fn extract_ca_from_chains_picks_every_third() {
        let chains = vec![vec![
            Vec3::new(0.0, 0.0, 0.0), // N
            Vec3::new(1.0, 1.0, 1.0), // CA
            Vec3::new(2.0, 2.0, 2.0), // C
            Vec3::new(3.0, 3.0, 3.0), // N
            Vec3::new(4.0, 4.0, 4.0), // CA
            Vec3::new(5.0, 5.0, 5.0), // C
        ]];
        let cas = extract_ca_from_chains(&chains);
        assert_eq!(cas.len(), 2);
        assert!((cas[0].x - 1.0).abs() < 1e-6);
        assert!((cas[1].x - 4.0).abs() < 1e-6);
    }

    #[test]
    fn extract_ca_from_chains_empty() {
        let cas = extract_ca_from_chains(&[]);
        assert!(cas.is_empty());
    }

    // -- centroid --

    #[test]
    fn centroid_single_point() {
        let c = centroid(&[Vec3::new(3.0, 4.0, 5.0)]);
        assert!((c.x - 3.0).abs() < 1e-6);
        assert!((c.y - 4.0).abs() < 1e-6);
        assert!((c.z - 5.0).abs() < 1e-6);
    }

    #[test]
    fn centroid_two_points() {
        let c = centroid(&[Vec3::new(0.0, 0.0, 0.0), Vec3::new(4.0, 6.0, 8.0)]);
        assert!((c.x - 2.0).abs() < 1e-6);
        assert!((c.y - 3.0).abs() < 1e-6);
        assert!((c.z - 4.0).abs() < 1e-6);
    }

    #[test]
    fn centroid_empty() {
        let c = centroid(&[]);
        assert_eq!(c, Vec3::ZERO);
    }
}

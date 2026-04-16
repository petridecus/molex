//! Tests for `protein.rs`. Split into a sibling file to keep the primary
//! module under the 800-line source cap enforced by `just file-lengths`.

#![allow(clippy::unwrap_used, clippy::float_cmp)]

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

fn make_two_residue_protein_coords() -> Coords {
    Coords {
        num_atoms: 8,
        atoms: vec![
            make_atom(1.0, 2.0, 3.0),    // res1 N
            make_atom(4.0, 5.0, 6.0),    // res1 CA
            make_atom(7.0, 8.0, 9.0),    // res1 C
            make_atom(10.0, 11.0, 12.0), // res1 O
            make_atom(13.0, 14.0, 15.0), // res2 N
            make_atom(16.0, 17.0, 18.0), // res2 CA
            make_atom(19.0, 20.0, 21.0), // res2 C
            make_atom(22.0, 23.0, 24.0), // res2 O
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
    }
}

// -- Sidechain tests --

#[test]
fn sidechain_empty() {
    let sc = Sidechain::empty();
    assert!(sc.is_empty());
    assert!(sc.atoms.is_empty());
    assert!(sc.bonds.is_empty());
    assert!(!sc.is_hydrophobic);
}

#[test]
fn sidechain_is_empty_with_atoms() {
    let sc = Sidechain {
        atoms: vec![Atom {
            position: Vec3::ZERO,
            occupancy: 1.0,
            b_factor: 0.0,
            element: Element::C,
            name: *b"CB  ",
        }],
        bonds: Vec::new(),
        is_hydrophobic: false,
    };
    assert!(!sc.is_empty());
}

// -- to_backbone tests --

#[test]
fn to_backbone_returns_two_entries() {
    let coords = make_two_residue_protein_coords();
    let entities = split_into_entities(&coords);
    let protein = entities[0].as_protein().unwrap();
    let backbone = protein.to_backbone();
    assert_eq!(backbone.len(), 2);
}

#[test]
fn to_backbone_first_residue_positions() {
    let coords = make_two_residue_protein_coords();
    let entities = split_into_entities(&coords);
    let protein = entities[0].as_protein().unwrap();
    let backbone = protein.to_backbone();

    let bb0 = &backbone[0];
    assert!((bb0.n.x - 1.0).abs() < 1e-6);
    assert!((bb0.ca.x - 4.0).abs() < 1e-6);
    assert!((bb0.c.x - 7.0).abs() < 1e-6);
    assert!((bb0.o.x - 10.0).abs() < 1e-6);
}

#[test]
fn to_backbone_second_residue_positions() {
    let coords = make_two_residue_protein_coords();
    let entities = split_into_entities(&coords);
    let protein = entities[0].as_protein().unwrap();
    let backbone = protein.to_backbone();

    let bb1 = &backbone[1];
    assert!((bb1.n.x - 13.0).abs() < 1e-6);
    assert!((bb1.ca.x - 16.0).abs() < 1e-6);
    assert!((bb1.c.x - 19.0).abs() < 1e-6);
    assert!((bb1.o.x - 22.0).abs() < 1e-6);
}

// -- segment breaks --

#[test]
fn no_segment_break_when_residues_close() {
    let coords = make_two_residue_protein_coords();
    let entities = split_into_entities(&coords);
    let protein = entities[0].as_protein().unwrap();
    // C of res1 is at (7,8,9), N of res2 is at (13,14,15)
    // distance ~10.4, which exceeds MAX_PEPTIDE_BOND_DIST, so there is a
    // break
    assert!(!protein.segment_breaks.is_empty());
}

#[test]
fn segment_break_detected_for_large_gap() {
    // Build residues where C->N distance > 2.0 A
    let coords = Coords {
        num_atoms: 8,
        atoms: vec![
            make_atom(0.0, 0.0, 0.0),  // N
            make_atom(1.0, 0.0, 0.0),  // CA
            make_atom(2.0, 0.0, 0.0),  // C
            make_atom(2.0, 1.0, 0.0),  // O
            make_atom(20.0, 0.0, 0.0), // N (far away)
            make_atom(21.0, 0.0, 0.0), // CA
            make_atom(22.0, 0.0, 0.0), // C
            make_atom(22.0, 1.0, 0.0), // O
        ],
        chain_ids: vec![b'A'; 8],
        res_names: vec![res_name("ALA"); 8],
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
        elements: vec![Element::Unknown; 8],
    };
    let entities = split_into_entities(&coords);
    let protein = entities[0].as_protein().unwrap();
    assert!(
        !protein.segment_breaks.is_empty(),
        "should have segment break at large gap"
    );
}

#[test]
fn no_segment_break_for_close_residues() {
    // Build residues where C->N distance < 2.0 A (peptide bond)
    let coords = Coords {
        num_atoms: 8,
        atoms: vec![
            make_atom(0.0, 0.0, 0.0), // N
            make_atom(1.5, 0.0, 0.0), // CA
            make_atom(2.5, 0.0, 0.0), // C
            make_atom(2.5, 1.0, 0.0), // O
            make_atom(3.8, 0.0, 0.0), // N (1.3 A from C, normal peptide)
            make_atom(5.0, 0.0, 0.0), // CA
            make_atom(6.0, 0.0, 0.0), // C
            make_atom(6.0, 1.0, 0.0), // O
        ],
        chain_ids: vec![b'A'; 8],
        res_names: vec![res_name("ALA"); 8],
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
        elements: vec![Element::Unknown; 8],
    };
    let entities = split_into_entities(&coords);
    let protein = entities[0].as_protein().unwrap();
    assert!(
        protein.segment_breaks.is_empty(),
        "should have no segment break for peptide-bonded residues"
    );
}

// -- to_interleaved_segments --

#[test]
fn interleaved_segments_for_single_segment() {
    let coords = Coords {
        num_atoms: 8,
        atoms: vec![
            make_atom(0.0, 0.0, 0.0),
            make_atom(1.5, 0.0, 0.0),
            make_atom(2.5, 0.0, 0.0),
            make_atom(2.5, 1.0, 0.0),
            make_atom(3.8, 0.0, 0.0),
            make_atom(5.0, 0.0, 0.0),
            make_atom(6.0, 0.0, 0.0),
            make_atom(6.0, 1.0, 0.0),
        ],
        chain_ids: vec![b'A'; 8],
        res_names: vec![res_name("ALA"); 8],
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
        elements: vec![Element::Unknown; 8],
    };
    let entities = split_into_entities(&coords);
    let protein = entities[0].as_protein().unwrap();
    let segments = protein.to_interleaved_segments();
    assert_eq!(segments.len(), 1);
    // 2 residues * 3 backbone atoms (N, CA, C) = 6
    assert_eq!(segments[0].len(), 6);
}

// -- to_sidechains tests --

#[test]
fn to_sidechains_per_residue() {
    // ALA(+CB) then GLY: 2 residues, first has sidechain, second empty
    let coords = make_residue_coords(&[
        ("ALA", &["N", "CA", "C", "O", "CB"]),
        ("GLY", &["N", "CA", "C", "O"]),
    ]);
    let entities = split_into_entities(&coords);
    let protein = entities[0].as_protein().unwrap();
    let sc = protein.to_sidechains(|_| false, |_| None);
    assert_eq!(sc.len(), 2);
    assert_eq!(sc[0].atoms.len(), 1); // CB
    assert!(sc[1].is_empty()); // GLY
}

#[test]
fn to_sidechains_bonds_local() {
    // VAL with CB, CG1, CG2 — bonds should use local indices
    let coords = make_residue_coords(&[(
        "VAL",
        &["N", "CA", "C", "O", "CB", "CG1", "CG2"],
    )]);
    let entities = split_into_entities(&coords);
    let protein = entities[0].as_protein().unwrap();
    let sc = protein.to_sidechains(
        |_| false,
        |name| match name {
            "VAL" => Some(vec![("CB", "CG1"), ("CB", "CG2")]),
            _ => None,
        },
    );
    assert_eq!(sc.len(), 1);
    assert_eq!(sc[0].atoms.len(), 3);
    assert_eq!(sc[0].bonds.len(), 2);
    for &(a, b) in &sc[0].bonds {
        assert!(a < 3 && b < 3);
    }
}

/// Build Coords for one or more residues from (name, atom_names) pairs.
fn make_residue_coords(residues: &[(&str, &[&str])]) -> Coords {
    let mut atoms = Vec::new();
    let mut chain_ids = Vec::new();
    let mut res_names_v = Vec::new();
    let mut res_nums = Vec::new();
    let mut atom_names_v = Vec::new();
    let mut elements = Vec::new();
    let mut x = 0.0f32;
    let mut res_num: i32 = 1;
    for (rname, anames) in residues {
        for aname in *anames {
            atoms.push(make_atom(x, 0.0, 0.0));
            chain_ids.push(b'A');
            res_names_v.push(res_name(rname));
            res_nums.push(res_num);
            atom_names_v.push(atom_name(aname));
            elements.push(Element::Unknown);
            x += 1.5;
        }
        res_num += 1;
    }
    Coords {
        num_atoms: atoms.len(),
        atoms,
        chain_ids,
        res_names: res_names_v,
        res_nums,
        atom_names: atom_names_v,
        elements,
    }
}

// -- is_hydrogen helper --

#[test]
fn is_hydrogen_recognizes_h_names() {
    assert!(is_hydrogen("H"));
    assert!(is_hydrogen("HA"));
    assert!(is_hydrogen("HB2"));
    assert!(is_hydrogen("1H"));
    assert!(is_hydrogen("2HB"));
    assert!(is_hydrogen("3HD1"));
}

#[test]
fn is_hydrogen_rejects_non_h() {
    assert!(!is_hydrogen("CA"));
    assert!(!is_hydrogen("N"));
    assert!(!is_hydrogen("O"));
    assert!(!is_hydrogen("CB"));
}

// -- Canonical ordering + bond population (Phase 2) --

fn alanine_scrambled_coords() -> Coords {
    // Input atoms: CB, O, N, HA, CA, C — non-canonical order.
    // Expected canonical: N, CA, C, O, CB, HA.
    Coords {
        num_atoms: 6,
        atoms: vec![
            make_atom(50.0, 0.0, 0.0),   // CB
            make_atom(10.0, 11.0, 12.0), // O
            make_atom(1.0, 2.0, 3.0),    // N
            make_atom(60.0, 0.0, 0.0),   // HA
            make_atom(4.0, 5.0, 6.0),    // CA
            make_atom(7.0, 8.0, 9.0),    // C
        ],
        chain_ids: vec![b'A'; 6],
        res_names: vec![res_name("ALA"); 6],
        res_nums: vec![1; 6],
        atom_names: vec![
            atom_name("CB"),
            atom_name("O"),
            atom_name("N"),
            atom_name("HA"),
            atom_name("CA"),
            atom_name("C"),
        ],
        elements: vec![
            Element::C,
            Element::O,
            Element::N,
            Element::H,
            Element::C,
            Element::C,
        ],
    }
}

#[test]
fn canonical_ordering_reorders_scrambled_input() {
    let coords = alanine_scrambled_coords();
    let entities = split_into_entities(&coords);
    let protein = entities[0].as_protein().unwrap();
    assert_eq!(protein.residues.len(), 1);
    let r = &protein.residues[0];
    // N CA C O CB HA — positions 0..6 of canonical layout.
    let names: Vec<&str> = r
        .atom_range
        .clone()
        .map(|i| std::str::from_utf8(&protein.atoms[i].name).unwrap().trim())
        .collect();
    assert_eq!(names, vec!["N", "CA", "C", "O", "CB", "HA"]);
}

#[test]
fn canonical_ordering_drops_residue_missing_backbone() {
    // GLY has only N, CA, C (missing O) — must be dropped.
    let coords = Coords {
        num_atoms: 7,
        #[allow(clippy::cast_precision_loss)]
        atoms: (0..7i32).map(|i| make_atom(i as f32, 0.0, 0.0)).collect(),
        chain_ids: vec![b'A'; 7],
        res_names: vec![
            res_name("ALA"),
            res_name("ALA"),
            res_name("ALA"),
            res_name("ALA"),
            res_name("GLY"),
            res_name("GLY"),
            res_name("GLY"),
        ],
        res_nums: vec![1, 1, 1, 1, 2, 2, 2],
        atom_names: vec![
            atom_name("N"),
            atom_name("CA"),
            atom_name("C"),
            atom_name("O"),
            atom_name("N"),
            atom_name("CA"),
            atom_name("C"),
        ],
        elements: vec![
            Element::N,
            Element::C,
            Element::C,
            Element::O,
            Element::N,
            Element::C,
            Element::C,
        ],
    };
    let entities = split_into_entities(&coords);
    let protein = entities[0].as_protein().unwrap();
    // Only ALA is kept.
    assert_eq!(protein.residues.len(), 1);
    assert_eq!(&protein.residues[0].name, b"ALA");
    // All 7 atoms still reside in entity.atoms.
    assert_eq!(protein.atoms.len(), 7);
}

#[test]
fn gly_backbone_bonds_are_three_sidechain_empty() {
    // Single GLY residue with canonical N, CA, C, O.
    let coords = Coords {
        num_atoms: 4,
        atoms: vec![
            make_atom(0.0, 0.0, 0.0),
            make_atom(1.5, 0.0, 0.0),
            make_atom(2.5, 0.0, 0.0),
            make_atom(2.5, 1.0, 0.0),
        ],
        chain_ids: vec![b'A'; 4],
        res_names: vec![res_name("GLY"); 4],
        res_nums: vec![1; 4],
        atom_names: vec![
            atom_name("N"),
            atom_name("CA"),
            atom_name("C"),
            atom_name("O"),
        ],
        elements: vec![Element::N, Element::C, Element::C, Element::O],
    };
    let entities = split_into_entities(&coords);
    let protein = entities[0].as_protein().unwrap();
    // Expect exactly N-CA, CA-C, C=O.
    assert_eq!(protein.backbone_bonds().count(), 3);
    assert_eq!(protein.sidechain_bonds().count(), 0);
}

#[test]
fn ala_sidechain_bonds_is_ca_cb_only() {
    // Single ALA residue with N, CA, C, O, CB.
    let coords = Coords {
        num_atoms: 5,
        atoms: vec![
            make_atom(0.0, 0.0, 0.0),
            make_atom(1.5, 0.0, 0.0),
            make_atom(2.5, 0.0, 0.0),
            make_atom(2.5, 1.0, 0.0),
            make_atom(1.5, -1.5, 0.0),
        ],
        chain_ids: vec![b'A'; 5],
        res_names: vec![res_name("ALA"); 5],
        res_nums: vec![1; 5],
        atom_names: vec![
            atom_name("N"),
            atom_name("CA"),
            atom_name("C"),
            atom_name("O"),
            atom_name("CB"),
        ],
        elements: vec![
            Element::N,
            Element::C,
            Element::C,
            Element::O,
            Element::C,
        ],
    };
    let entities = split_into_entities(&coords);
    let protein = entities[0].as_protein().unwrap();
    let sidechain: Vec<_> = protein.sidechain_bonds().collect();
    assert_eq!(sidechain.len(), 1);
    let b = sidechain[0];
    let a_name = std::str::from_utf8(&protein.atoms[b.a.index as usize].name)
        .unwrap()
        .trim();
    let c_name = std::str::from_utf8(&protein.atoms[b.b.index as usize].name)
        .unwrap()
        .trim();
    let pair = (a_name, c_name);
    assert!(
        pair == ("CA", "CB") || pair == ("CB", "CA"),
        "expected CA-CB sidechain anchor, got {pair:?}"
    );
}

#[test]
fn peptide_bond_connects_consecutive_residues() {
    // Two close ALA residues — peptide bond expected.
    let coords = Coords {
        num_atoms: 8,
        atoms: vec![
            make_atom(0.0, 0.0, 0.0),
            make_atom(1.5, 0.0, 0.0),
            make_atom(2.5, 0.0, 0.0),
            make_atom(2.5, 1.0, 0.0),
            make_atom(3.8, 0.0, 0.0), // N of residue 2 close to C of 1
            make_atom(5.0, 0.0, 0.0),
            make_atom(6.0, 0.0, 0.0),
            make_atom(6.0, 1.0, 0.0),
        ],
        chain_ids: vec![b'A'; 8],
        res_names: vec![res_name("ALA"); 8],
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
        elements: vec![Element::Unknown; 8],
    };
    let entities = split_into_entities(&coords);
    let protein = entities[0].as_protein().unwrap();
    // 3 backbone * 2 residues + 1 peptide = 7 backbone bonds total.
    let backbone_count = protein.backbone_bonds().count();
    assert_eq!(backbone_count, 7, "3 per residue + 1 peptide");
}

#[test]
fn dropped_residue_emits_log_warning() {
    testing_logger::setup();
    // ALA residue missing CA — must be dropped and logged.
    let coords = Coords {
        num_atoms: 3,
        atoms: vec![
            make_atom(0.0, 0.0, 0.0),
            make_atom(1.5, 0.0, 0.0),
            make_atom(2.5, 1.0, 0.0),
        ],
        chain_ids: vec![b'A'; 3],
        res_names: vec![res_name("ALA"); 3],
        res_nums: vec![7; 3],
        atom_names: vec![atom_name("N"), atom_name("C"), atom_name("O")],
        elements: vec![Element::N, Element::C, Element::O],
    };
    let _ = split_into_entities(&coords);
    testing_logger::validate(|captured_logs| {
        let warn_bodies: Vec<&str> = captured_logs
            .iter()
            .filter(|l| l.level == log::Level::Warn)
            .map(|l| l.body.as_str())
            .collect();
        assert!(
            !warn_bodies.is_empty(),
            "expected at least one warn-level log entry"
        );
        assert!(
            warn_bodies.iter().any(|b| b.contains("dropping residue")),
            "warning should mention dropping a residue; got {warn_bodies:?}"
        );
    });
}

#[test]
fn peptide_bond_skipped_across_segment_break() {
    // Two ALA residues with a large C-N gap (segment break).
    let coords = Coords {
        num_atoms: 8,
        atoms: vec![
            make_atom(0.0, 0.0, 0.0),
            make_atom(1.5, 0.0, 0.0),
            make_atom(2.5, 0.0, 0.0),
            make_atom(2.5, 1.0, 0.0),
            make_atom(20.0, 0.0, 0.0), // N far from prev C
            make_atom(21.5, 0.0, 0.0),
            make_atom(22.5, 0.0, 0.0),
            make_atom(22.5, 1.0, 0.0),
        ],
        chain_ids: vec![b'A'; 8],
        res_names: vec![res_name("ALA"); 8],
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
        elements: vec![Element::Unknown; 8],
    };
    let entities = split_into_entities(&coords);
    let protein = entities[0].as_protein().unwrap();
    // No peptide bond — just 3 backbone per residue = 6.
    assert_eq!(protein.backbone_bonds().count(), 6);
}

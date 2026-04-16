#![allow(
    clippy::unwrap_used,
    clippy::cast_precision_loss,
    clippy::too_many_lines,
    clippy::float_cmp
)]

use super::*;
use crate::entity::molecule::{MoleculeEntity, MoleculeType};

fn make_atom(x: f32) -> CoordsAtom {
    CoordsAtom {
        x,
        y: 0.0,
        z: 0.0,
        occupancy: 1.0,
        b_factor: 0.0,
    }
}

fn res_name(s: &str) -> [u8; 3] {
    let mut name = [b' '; 3];
    for (i, b) in s.bytes().take(3).enumerate() {
        name[i] = b;
    }
    name
}

fn atom_name(s: &str) -> [u8; 4] {
    let mut name = [b' '; 4];
    for (i, b) in s.bytes().take(4).enumerate() {
        name[i] = b;
    }
    name
}

fn make_full_atom(x: f32, y: f32, z: f32) -> CoordsAtom {
    CoordsAtom {
        x,
        y,
        z,
        occupancy: 1.0,
        b_factor: 0.0,
    }
}

use crate::element::Element;

// -- Assembly helpers --

#[test]
fn test_update_protein_entities_no_duplication() {
    let coords = Coords {
        num_atoms: 7,
        atoms: (0..7).map(|i| make_atom(i as f32)).collect(),
        chain_ids: vec![b'A', b'A', b'A', b'B', b'B', b'B', b'C'],
        res_names: vec![
            res_name("ALA"),
            res_name("ALA"),
            res_name("ALA"),
            res_name("GLY"),
            res_name("GLY"),
            res_name("GLY"),
            res_name("HOH"),
        ],
        res_nums: vec![1, 1, 1, 1, 1, 1, 100],
        atom_names: vec![
            atom_name("N"),
            atom_name("CA"),
            atom_name("C"),
            atom_name("N"),
            atom_name("CA"),
            atom_name("C"),
            atom_name("O"),
        ],
        elements: vec![Element::Unknown; 7],
    };

    let mut entities = split_into_entities(&coords);
    assert_eq!(entities.len(), 3);
    let total: usize = entities.iter().map(MoleculeEntity::atom_count).sum();
    assert_eq!(total, 7);

    let updated_protein = Coords {
        num_atoms: 6,
        atoms: (10..16).map(|i| make_atom(i as f32)).collect(),
        chain_ids: vec![b'A', b'A', b'A', b'B', b'B', b'B'],
        res_names: vec![
            res_name("ALA"),
            res_name("ALA"),
            res_name("ALA"),
            res_name("GLY"),
            res_name("GLY"),
            res_name("GLY"),
        ],
        res_nums: vec![1, 1, 1, 1, 1, 1],
        atom_names: vec![
            atom_name("N"),
            atom_name("CA"),
            atom_name("C"),
            atom_name("N"),
            atom_name("CA"),
            atom_name("C"),
        ],
        elements: vec![Element::Unknown; 6],
    };

    update_protein_entities(&mut entities, &updated_protein);

    let total: usize = entities.iter().map(MoleculeEntity::atom_count).sum();
    assert_eq!(total, 7);
    assert_eq!(entities.len(), 3);

    let prot = protein_coords(&entities);
    assert_eq!(prot.num_atoms, 6);
    assert!(prot.atoms[0].x >= 10.0);
}

// -- Assembly byte layout --

#[test]
fn test_assembly_bytes_roundtrip_mixed() {
    let coords = Coords {
        num_atoms: 6,
        atoms: vec![
            make_atom(1.0),
            make_atom(2.0),
            make_atom(3.0),
            make_atom(4.0),
            make_atom(10.0),
            make_atom(20.0),
        ],
        chain_ids: vec![b'A', b'A', b'A', b'A', b'B', b'C'],
        res_names: vec![
            res_name("ALA"),
            res_name("ALA"),
            res_name("ALA"),
            res_name("ALA"),
            res_name("ATP"),
            res_name("ZN"),
        ],
        res_nums: vec![1, 1, 1, 1, 1, 1],
        atom_names: vec![
            atom_name("N"),
            atom_name("CA"),
            atom_name("C"),
            atom_name("O"),
            atom_name("C1"),
            atom_name("ZN"),
        ],
        elements: vec![
            Element::N,
            Element::C,
            Element::C,
            Element::O,
            Element::C,
            Element::Zn,
        ],
    };

    let entities = split_into_entities(&coords);
    assert!(entities.len() >= 3);

    let bytes = assembly_bytes(&entities).unwrap();
    assert_eq!(&bytes[0..8], b"ASSEM01\0");

    let roundtripped = deserialize_assembly(&bytes).unwrap();
    assert_eq!(roundtripped.len(), entities.len());

    for (orig, rt) in entities.iter().zip(roundtripped.iter()) {
        assert_eq!(orig.molecule_type(), rt.molecule_type());
        assert_eq!(orig.atom_count(), rt.atom_count());
    }
}

#[test]
fn test_assembly_bytes_protein_only() {
    let coords = Coords {
        num_atoms: 4,
        atoms: vec![
            make_atom(1.0),
            make_atom(2.0),
            make_atom(3.0),
            make_atom(4.0),
        ],
        chain_ids: vec![b'A'; 4],
        res_names: vec![res_name("ALA"); 4],
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
    let bytes = assembly_bytes(&entities).unwrap();
    let roundtripped = deserialize_assembly(&bytes).unwrap();
    assert_eq!(roundtripped.len(), 1);
    assert_eq!(roundtripped[0].molecule_type(), MoleculeType::Protein);
    assert_eq!(roundtripped[0].atom_count(), 4);
}

#[test]
fn test_assembly_bytes_single_atom_ion() {
    let coords = Coords {
        num_atoms: 1,
        atoms: vec![make_atom(5.5)],
        chain_ids: vec![b'X'],
        res_names: vec![res_name("ZN")],
        res_nums: vec![99],
        atom_names: vec![atom_name("ZN")],
        elements: vec![Element::Zn],
    };

    let entities = split_into_entities(&coords);
    let bytes = assembly_bytes(&entities).unwrap();
    let roundtripped = deserialize_assembly(&bytes).unwrap();
    assert_eq!(roundtripped.len(), 1);
    assert_eq!(roundtripped[0].molecule_type(), MoleculeType::Ion);
    assert!((roundtripped[0].atom_set()[0].position.x - 5.5).abs() < 1e-6);
}

#[test]
fn test_assembly_bytes_empty_entities() {
    let entities: Vec<MoleculeEntity> = Vec::new();
    let bytes = serialize_assembly(&entities).unwrap();
    let roundtripped = deserialize_assembly(&bytes).unwrap();
    assert!(roundtripped.is_empty());
}

#[test]
fn test_assembly_byte_layout() {
    let coords = Coords {
        num_atoms: 4,
        atoms: vec![
            CoordsAtom {
                x: 1.0,
                y: 2.0,
                z: 3.0,
                occupancy: 1.0,
                b_factor: 0.0,
            },
            CoordsAtom {
                x: 1.5,
                y: 2.5,
                z: 3.5,
                occupancy: 1.0,
                b_factor: 0.0,
            },
            CoordsAtom {
                x: 2.0,
                y: 3.0,
                z: 4.0,
                occupancy: 1.0,
                b_factor: 0.0,
            },
            CoordsAtom {
                x: 2.5,
                y: 3.5,
                z: 4.5,
                occupancy: 1.0,
                b_factor: 0.0,
            },
        ],
        chain_ids: vec![b'A'; 4],
        res_names: vec![res_name("ALA"); 4],
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
    let bytes = serialize_assembly(&entities).unwrap();

    // 8 magic + 4 count + 5 per-entity header + 4 atoms * 26 = 121.
    assert_eq!(&bytes[0..8], b"ASSEM01\0");
    assert_eq!(u32::from_be_bytes(bytes[8..12].try_into().unwrap()), 1);
    assert_eq!(bytes[12], 0); // Protein
    assert_eq!(u32::from_be_bytes(bytes[13..17].try_into().unwrap()), 4);
    // First atom starts at offset 17; canonical ordering may reorder
    // within residue, so we only check the header + total length.
    assert_eq!(bytes.len(), 8 + 4 + 5 + 4 * 26);
}

// -- Split/merge --

#[test]
fn test_split_protein_only() {
    let coords = Coords {
        num_atoms: 6,
        atoms: (0..6).map(|i| make_atom(i as f32)).collect(),
        chain_ids: vec![b'A'; 6],
        res_names: vec![
            res_name("ALA"),
            res_name("ALA"),
            res_name("ALA"),
            res_name("GLY"),
            res_name("GLY"),
            res_name("GLY"),
        ],
        res_nums: vec![1, 1, 1, 2, 2, 2],
        atom_names: vec![
            atom_name("N"),
            atom_name("CA"),
            atom_name("C"),
            atom_name("N"),
            atom_name("CA"),
            atom_name("C"),
        ],
        elements: vec![Element::Unknown; 6],
    };
    let entities = split_into_entities(&coords);
    assert_eq!(entities.len(), 1);
    assert_eq!(entities[0].molecule_type(), MoleculeType::Protein);
    assert_eq!(entities[0].atom_count(), 6);
}

#[test]
fn test_split_mixed() {
    let coords = Coords {
        num_atoms: 5,
        atoms: (0..5).map(|i| make_atom(i as f32)).collect(),
        chain_ids: vec![b'A', b'A', b'A', b'B', b'C'],
        res_names: vec![
            res_name("ALA"),
            res_name("ALA"),
            res_name("HOH"),
            res_name("ATP"),
            res_name("HOH"),
        ],
        res_nums: vec![1, 1, 100, 1, 200],
        atom_names: vec![
            atom_name("N"),
            atom_name("CA"),
            atom_name("O"),
            atom_name("C1"),
            atom_name("O"),
        ],
        elements: vec![Element::Unknown; 5],
    };
    let entities = split_into_entities(&coords);
    assert_eq!(entities.len(), 3);
}

#[test]
fn test_merge_roundtrip() {
    let coords = Coords {
        num_atoms: 4,
        atoms: (0..4).map(|i| make_atom(i as f32)).collect(),
        chain_ids: vec![b'A', b'A', b'B', b'B'],
        res_names: vec![
            res_name("ALA"),
            res_name("ALA"),
            res_name("HOH"),
            res_name("HOH"),
        ],
        res_nums: vec![1, 1, 100, 101],
        atom_names: vec![
            atom_name("N"),
            atom_name("CA"),
            atom_name("O"),
            atom_name("O"),
        ],
        elements: vec![Element::Unknown; 4],
    };
    let entities = split_into_entities(&coords);
    let merged = merge_entities(&entities);
    assert_eq!(merged.num_atoms, coords.num_atoms);
}

#[test]
fn test_extract_by_type() {
    let coords = Coords {
        num_atoms: 3,
        atoms: (0..3).map(|i| make_atom(i as f32)).collect(),
        chain_ids: vec![b'A', b'A', b'A'],
        res_names: vec![res_name("ALA"), res_name("HOH"), res_name("ZN")],
        res_nums: vec![1, 100, 200],
        atom_names: vec![atom_name("CA"), atom_name("O"), atom_name("ZN")],
        elements: vec![Element::Unknown; 3],
    };
    let entities = split_into_entities(&coords);
    let protein = extract_by_type(&entities, MoleculeType::Protein);
    assert!(protein.is_some());
    assert_eq!(protein.unwrap().num_atoms, 1);
    assert!(extract_by_type(&entities, MoleculeType::DNA).is_none());
}

#[test]
fn test_split_cofactor_grouping() {
    let coords = Coords {
        num_atoms: 4,
        atoms: (0..4).map(|i| make_atom(i as f32)).collect(),
        chain_ids: vec![b'A', b'A', b'D', b'D'],
        res_names: vec![res_name("CLA"); 4],
        res_nums: vec![1, 1, 2, 2],
        atom_names: vec![
            atom_name("MG"),
            atom_name("NA"),
            atom_name("MG"),
            atom_name("NA"),
        ],
        elements: vec![Element::Unknown; 4],
    };
    let entities = split_into_entities(&coords);
    assert_eq!(entities.len(), 2);
    assert_eq!(entities[0].molecule_type(), MoleculeType::Cofactor);
}

#[test]
fn test_split_solvent_consolidated() {
    let coords = Coords {
        num_atoms: 3,
        atoms: (0..3).map(|i| make_atom(i as f32)).collect(),
        chain_ids: vec![b'A', b'B', b'C'],
        res_names: vec![res_name("GOL"), res_name("SO4"), res_name("GOL")],
        res_nums: vec![1, 2, 3],
        atom_names: vec![atom_name("O1"), atom_name("S"), atom_name("O1")],
        elements: vec![Element::Unknown; 3],
    };
    let entities = split_into_entities(&coords);
    assert_eq!(entities.len(), 1);
    assert_eq!(entities[0].molecule_type(), MoleculeType::Solvent);
}

#[test]
fn test_polymer_structure() {
    let coords = Coords {
        num_atoms: 8,
        atoms: (0..8).map(|i| make_atom(i as f32)).collect(),
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
        elements: vec![Element::Unknown; 8],
    };
    let entities = split_into_entities(&coords);
    assert_eq!(entities.len(), 1);
    let protein = entities[0].as_protein().unwrap();
    assert_eq!(protein.pdb_chain_id, b'A');
    assert_eq!(protein.residues.len(), 2);
}

#[test]
fn test_small_molecule_no_chain_residue() {
    let coords = Coords {
        num_atoms: 1,
        atoms: vec![make_atom(1.0)],
        chain_ids: vec![b'A'],
        res_names: vec![res_name("ZN")],
        res_nums: vec![300],
        atom_names: vec![atom_name("ZN")],
        elements: vec![Element::Zn],
    };
    let entities = split_into_entities(&coords);
    assert_eq!(entities.len(), 1);
    assert_eq!(entities[0].molecule_type(), MoleculeType::Ion);
}

#[test]
fn test_to_coords_roundtrip() {
    let coords = Coords {
        num_atoms: 6,
        atoms: (0..6).map(|i| make_atom(i as f32)).collect(),
        chain_ids: vec![b'A'; 6],
        res_names: vec![
            res_name("ALA"),
            res_name("ALA"),
            res_name("ALA"),
            res_name("GLY"),
            res_name("GLY"),
            res_name("GLY"),
        ],
        res_nums: vec![1, 1, 1, 2, 2, 2],
        atom_names: vec![
            atom_name("N"),
            atom_name("CA"),
            atom_name("C"),
            atom_name("N"),
            atom_name("CA"),
            atom_name("C"),
        ],
        elements: vec![Element::Unknown; 6],
    };
    let entities = split_into_entities(&coords);
    let recovered = entities[0].to_coords();
    assert_eq!(recovered.num_atoms, 6);
    assert_eq!(recovered.chain_ids, vec![b'A'; 6]);
}

#[test]
fn test_split_modified_amino_acid_merges_into_protein() {
    let coords = Coords {
        num_atoms: 9,
        atoms: (0..9).map(|i| make_atom(i as f32)).collect(),
        chain_ids: vec![b'A'; 9],
        res_names: vec![
            res_name("ALA"),
            res_name("ALA"),
            res_name("ALA"),
            res_name("SEP"),
            res_name("SEP"),
            res_name("SEP"),
            res_name("GLY"),
            res_name("GLY"),
            res_name("GLY"),
        ],
        res_nums: vec![1, 1, 1, 2, 2, 2, 3, 3, 3],
        atom_names: vec![
            atom_name("N"),
            atom_name("CA"),
            atom_name("C"),
            atom_name("N"),
            atom_name("CA"),
            atom_name("C"),
            atom_name("N"),
            atom_name("CA"),
            atom_name("C"),
        ],
        elements: vec![Element::Unknown; 9],
    };
    let entities = split_into_entities(&coords);
    assert_eq!(entities.len(), 1);
    assert_eq!(entities[0].molecule_type(), MoleculeType::Protein);
    assert_eq!(entities[0].atom_count(), 9);
}

// -- COORDS01 serialize/deserialize round-trip --

#[test]
fn test_coords_serialize_deserialize_roundtrip() {
    let coords = Coords {
        num_atoms: 3,
        atoms: vec![
            make_full_atom(1.5, 2.5, 3.5),
            make_full_atom(4.0, 5.0, 6.0),
            make_full_atom(7.0, 8.0, 9.0),
        ],
        chain_ids: vec![b'A', b'A', b'B'],
        res_names: vec![res_name("ALA"), res_name("ALA"), res_name("GLY")],
        res_nums: vec![1, 1, 2],
        atom_names: vec![atom_name("N"), atom_name("CA"), atom_name("C")],
        elements: vec![Element::N, Element::C, Element::C],
    };
    let bytes = serialize(&coords).unwrap();
    let recovered = deserialize(&bytes).unwrap();

    assert_eq!(recovered.num_atoms, 3);
    assert!((recovered.atoms[0].x - 1.5).abs() < 1e-6);
    assert_eq!(recovered.chain_ids, vec![b'A', b'A', b'B']);
    assert_eq!(recovered.elements[0], Element::N);
}

#[test]
fn test_coords_serialize_preserves_elements() {
    let coords = Coords {
        num_atoms: 2,
        atoms: vec![
            make_full_atom(0.0, 0.0, 0.0),
            make_full_atom(1.0, 1.0, 1.0),
        ],
        chain_ids: vec![b'A'; 2],
        res_names: vec![res_name("ALA"); 2],
        res_nums: vec![1; 2],
        atom_names: vec![atom_name("FE"), atom_name("ZN")],
        elements: vec![Element::Fe, Element::Zn],
    };
    let bytes = serialize(&coords).unwrap();
    let recovered = deserialize(&bytes).unwrap();
    assert_eq!(recovered.elements[0], Element::Fe);
    assert_eq!(recovered.elements[1], Element::Zn);
}

#[test]
fn test_coords_single_atom_roundtrip() {
    let coords = Coords {
        num_atoms: 1,
        atoms: vec![make_full_atom(42.0, -17.5, 0.001)],
        chain_ids: vec![b'Z'],
        res_names: vec![res_name("ZN")],
        res_nums: vec![999],
        atom_names: vec![atom_name("ZN")],
        elements: vec![Element::Zn],
    };
    let bytes = serialize(&coords).unwrap();
    let recovered = deserialize(&bytes).unwrap();
    assert_eq!(recovered.num_atoms, 1);
    assert!((recovered.atoms[0].x - 42.0).abs() < 1e-6);
    assert_eq!(recovered.elements[0], Element::Zn);
}

// -- Error cases --

#[test]
fn test_deserialize_empty_bytes() {
    assert!(deserialize(&[]).is_err());
}

#[test]
fn test_deserialize_wrong_magic() {
    assert!(deserialize(b"BADMAGIC\x00\x00\x00\x01").is_err());
}

#[test]
fn test_deserialize_truncated_data() {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(b"COORDS01");
    bytes.extend_from_slice(&1u32.to_be_bytes());
    assert!(deserialize(&bytes).is_err());
}

#[test]
fn test_deserialize_assembly_empty_bytes() {
    assert!(deserialize_assembly(&[]).is_err());
}

#[test]
fn test_deserialize_assembly_wrong_magic() {
    assert!(deserialize_assembly(b"BADMAGIC\x00\x00\x00\x00").is_err());
}

#[test]
fn test_deserialize_assembly_truncated_atom_data() {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(b"ASSEM01\0");
    bytes.extend_from_slice(&1u32.to_be_bytes()); // 1 entity
    bytes.push(0); // Protein type
    bytes.extend_from_slice(&1u32.to_be_bytes()); // 1 atom, no atom data
    assert!(deserialize_assembly(&bytes).is_err());
}

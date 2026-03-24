//! PDB format parsing and writing.

use std::fmt::Write as _;
use std::io::BufReader;

use pdbtbx::{
    ContainsAtomConformer, ContainsAtomConformerResidue,
    ContainsAtomConformerResidueChain, ReadOptions, StrictnessLevel,
};

use crate::entity::molecule::MoleculeEntity;
use crate::ops::codec::{
    deserialize, merge_entities, serialize, split_into_entities, Coords,
    CoordsError,
};

// ---------------------------------------------------------------------------
// Entity-first API (primary)
// ---------------------------------------------------------------------------

/// Parse PDB format string to entity list.
///
/// # Errors
///
/// Returns [`CoordsError`] if parsing fails.
pub fn pdb_str_to_entities(
    pdb_str: &str,
) -> Result<Vec<MoleculeEntity>, CoordsError> {
    parse_pdb_to_entities(pdb_str)
}

/// Load PDB file to entity list.
///
/// Sanitizes non-standard lines (e.g. GROMACS/MemProtMD output) before
/// parsing.
///
/// # Errors
///
/// Returns [`CoordsError`] if the file cannot be read or parsing fails.
pub fn pdb_file_to_entities(
    path: &std::path::Path,
) -> Result<Vec<MoleculeEntity>, CoordsError> {
    let content = std::fs::read_to_string(path).map_err(|e| {
        CoordsError::PdbParseError(format!("Failed to read file: {e}"))
    })?;
    let sanitized = sanitize_pdb(&content);
    parse_pdb_to_entities(&sanitized)
}

/// Load structure file (PDB or mmCIF, detected by extension) to entity list.
///
/// # Errors
///
/// Returns [`CoordsError`] if the file cannot be read or parsing fails.
pub fn structure_file_to_entities(
    path: &std::path::Path,
) -> Result<Vec<MoleculeEntity>, CoordsError> {
    let ext = path
        .extension()
        .and_then(|e| e.to_str())
        .unwrap_or("")
        .to_lowercase();
    match ext.as_str() {
        "pdb" | "ent" => pdb_file_to_entities(path),
        _ => super::cif::mmcif_file_to_entities(path),
    }
}

// ---------------------------------------------------------------------------
// Coords/serialization API (derived from entities)
// ---------------------------------------------------------------------------

/// Parse PDB format string to COORDS binary format.
///
/// # Errors
///
/// Returns [`CoordsError`] if parsing or serialization fails.
pub fn pdb_to_coords(pdb_str: &str) -> Result<Vec<u8>, CoordsError> {
    let entities = parse_pdb_to_entities(pdb_str)?;
    serialize(&merge_entities(&entities))
}

/// Parse PDB format string directly to Coords struct.
///
/// # Errors
///
/// Returns [`CoordsError`] if parsing fails.
pub fn pdb_str_to_coords(pdb_str: &str) -> Result<Coords, CoordsError> {
    let entities = parse_pdb_to_entities(pdb_str)?;
    Ok(merge_entities(&entities))
}

/// Load PDB file directly to Coords struct.
///
/// # Errors
///
/// Returns [`CoordsError`] if the file cannot be read or parsing fails.
pub fn pdb_file_to_coords(
    path: &std::path::Path,
) -> Result<Coords, CoordsError> {
    let entities = pdb_file_to_entities(path)?;
    Ok(merge_entities(&entities))
}

/// Load a structure file (PDB or mmCIF) by detecting format from extension.
///
/// # Errors
///
/// Returns [`CoordsError`] if the file cannot be read or parsing fails.
pub fn structure_file_to_coords(
    path: &std::path::Path,
) -> Result<Coords, CoordsError> {
    let entities = structure_file_to_entities(path)?;
    Ok(merge_entities(&entities))
}

/// Convert COORDS binary to PDB format string.
///
/// # Errors
///
/// Returns [`CoordsError`] if deserialization fails.
pub fn coords_to_pdb(coords_bytes: &[u8]) -> Result<String, CoordsError> {
    let coords = deserialize(coords_bytes)?;

    let mut pdb_string = String::new();

    for i in 0..coords.num_atoms {
        let atom = &coords.atoms[i];
        let chain_id = coords.chain_ids[i] as char;
        let res_num = coords.res_nums[i];

        let atom_name =
            std::str::from_utf8(&coords.atom_names[i]).unwrap_or("X   ");
        let res_name =
            std::str::from_utf8(&coords.res_names[i]).unwrap_or("UNK");

        let _ = writeln!(
            pdb_string,
            "ATOM  {:>5} {:<4} {:>3} {}{:>4}    \
             {:>8.3}{:>8.3}{:>8.3}{:>6.2}{:>6.2}",
            i + 1,
            atom_name,
            res_name,
            chain_id,
            res_num,
            atom.x,
            atom.y,
            atom.z,
            atom.occupancy,
            atom.b_factor
        );
    }

    pdb_string.push_str("END\n");

    Ok(pdb_string)
}

// ---------------------------------------------------------------------------
// Internal parsing
// ---------------------------------------------------------------------------

/// Strip or fix PDB lines that cause pdbtbx to error even in Loose mode.
fn sanitize_pdb(content: &str) -> String {
    content
        .lines()
        .map(|line| {
            line.strip_prefix("REMARK").map_or_else(
                || line.to_owned(),
                |after| {
                    let trimmed = after.trim_start();
                    if trimmed.is_empty()
                        || !trimmed.as_bytes()[0].is_ascii_digit()
                    {
                        format!("REMARK   0 {trimmed}")
                    } else {
                        line.to_owned()
                    }
                },
            )
        })
        .collect::<Vec<_>>()
        .join("\n")
}

/// Parse PDB string into entities via temporary Coords + split.
fn parse_pdb_to_entities(
    input: &str,
) -> Result<Vec<MoleculeEntity>, CoordsError> {
    let coords = parse_pdb_to_coords(input)?;
    Ok(split_into_entities(&coords))
}

/// Parse PDB string into flat Coords (internal, used by entity pipeline).
fn parse_pdb_to_coords(input: &str) -> Result<Coords, CoordsError> {
    use crate::element::Element;
    use crate::ops::codec::{ChainIdMapper, CoordsAtom};

    let reader = BufReader::new(input.as_bytes());

    let (pdb, _errors) = ReadOptions::new()
        .set_format(pdbtbx::Format::Pdb)
        .set_level(StrictnessLevel::Loose)
        .read_raw(reader)
        .map_err(|errs| {
            CoordsError::PdbParseError(
                errs.iter()
                    .map(ToString::to_string)
                    .collect::<Vec<_>>()
                    .join("; "),
            )
        })?;

    let mut atoms = Vec::new();
    let mut chain_ids = Vec::new();
    let mut res_names = Vec::new();
    let mut res_nums = Vec::new();
    let mut atom_names = Vec::new();
    let mut elements = Vec::new();
    let mut chain_mapper = ChainIdMapper::new();

    for hier in pdb.atoms_with_hierarchy() {
        let atom = hier.atom();
        let chain = hier.chain();
        let residue = hier.residue();
        let conformer = hier.conformer();

        #[allow(clippy::cast_possible_truncation)]
        atoms.push(CoordsAtom {
            x: atom.x() as f32,
            y: atom.y() as f32,
            z: atom.z() as f32,
            occupancy: atom.occupancy() as f32,
            b_factor: atom.b_factor() as f32,
        });

        chain_ids.push(chain_mapper.get_or_assign(chain.id()));
        res_names.push(name_to_bytes::<3>(conformer.name()));
        #[allow(clippy::cast_possible_truncation)]
        res_nums.push(residue.serial_number() as i32);
        let aname = atom.name();
        atom_names.push(name_to_bytes::<4>(aname));
        elements.push(atom.element().map_or_else(
            || Element::from_atom_name(aname),
            |e| Element::from_symbol(e.symbol()),
        ));
    }

    if atoms.is_empty() {
        return Err(CoordsError::PdbParseError(
            "No atoms found in structure".to_owned(),
        ));
    }

    Ok(Coords {
        num_atoms: atoms.len(),
        atoms,
        chain_ids,
        res_names,
        res_nums,
        atom_names,
        elements,
    })
}

/// Convert a string to a space-padded byte array of length N.
fn name_to_bytes<const N: usize>(name: &str) -> [u8; N] {
    let mut buf = [b' '; N];
    for (i, b) in name.bytes().take(N).enumerate() {
        buf[i] = b;
    }
    buf
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::float_cmp)]
mod tests {
    use super::*;
    use crate::entity::molecule::MoleculeType;

    /// Minimal PDB with one residue (N, CA, C, O) for alanine.
    const MINIMAL_PDB: &str = "\
ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       4.000   5.000   6.000  1.00  0.00           C
ATOM      3  C   ALA A   1       7.000   8.000   9.000  1.00  0.00           C
ATOM      4  O   ALA A   1      10.000  11.000  12.000  1.00  0.00           O
END
";

    #[test]
    fn pdb_str_to_entities_minimal() {
        let entities = pdb_str_to_entities(MINIMAL_PDB).unwrap();
        assert!(!entities.is_empty());
        let protein = entities
            .iter()
            .find(|e| e.molecule_type() == MoleculeType::Protein);
        assert!(protein.is_some());
        assert_eq!(protein.unwrap().atom_count(), 4);
    }

    #[test]
    fn pdb_str_to_entities_preserves_positions() {
        let entities = pdb_str_to_entities(MINIMAL_PDB).unwrap();
        let protein = entities
            .iter()
            .find(|e| e.molecule_type() == MoleculeType::Protein)
            .unwrap();
        let pos = &protein.atom_set()[0].position;
        assert!((pos.x - 1.0).abs() < 0.01);
        assert!((pos.y - 2.0).abs() < 0.01);
        assert!((pos.z - 3.0).abs() < 0.01);
    }

    #[test]
    fn pdb_str_to_entities_two_chains() {
        let pdb = "\
ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       4.000   5.000   6.000  1.00  0.00           C
ATOM      3  C   ALA A   1       7.000   8.000   9.000  1.00  0.00           C
ATOM      4  N   GLY B   1       1.000   2.000   3.000  1.00  0.00           N
ATOM      5  CA  GLY B   1       4.000   5.000   6.000  1.00  0.00           C
ATOM      6  C   GLY B   1       7.000   8.000   9.000  1.00  0.00           C
END
";
        let entities = pdb_str_to_entities(pdb).unwrap();
        let protein_count = entities
            .iter()
            .filter(|e| e.molecule_type() == MoleculeType::Protein)
            .count();
        assert_eq!(protein_count, 2);
    }

    #[test]
    fn pdb_str_to_entities_water() {
        let pdb = "\
ATOM      1  O   HOH A 100       1.000   2.000   3.000  1.00  0.00           O
ATOM      2  O   HOH A 101       4.000   5.000   6.000  1.00  0.00           O
END
";
        let entities = pdb_str_to_entities(pdb).unwrap();
        let water = entities
            .iter()
            .find(|e| e.molecule_type() == MoleculeType::Water);
        assert!(water.is_some());
        assert_eq!(water.unwrap().atom_count(), 2);
    }

    #[test]
    fn coords_to_pdb_produces_valid_output() {
        let entities = pdb_str_to_entities(MINIMAL_PDB).unwrap();
        let merged = merge_entities(&entities);
        let bytes = serialize(&merged).unwrap();
        let pdb_output = coords_to_pdb(&bytes).unwrap();

        assert!(pdb_output.contains("ATOM"));
        assert!(pdb_output.contains("ALA"));
        assert!(pdb_output.ends_with("END\n"));
        // Should have 4 ATOM lines + END
        let atom_line_count =
            pdb_output.lines().filter(|l| l.starts_with("ATOM")).count();
        assert_eq!(atom_line_count, 4);
    }

    #[test]
    fn coords_to_pdb_preserves_coordinates_in_output() {
        let entities = pdb_str_to_entities(MINIMAL_PDB).unwrap();
        let merged = merge_entities(&entities);
        let bytes = serialize(&merged).unwrap();
        let pdb_output = coords_to_pdb(&bytes).unwrap();

        // Verify coordinate values appear in the output
        assert!(pdb_output.contains("1.000"));
        assert!(pdb_output.contains("2.000"));
        assert!(pdb_output.contains("3.000"));
        // Verify chain ID present
        assert!(pdb_output.contains('A'));
        // Verify residue numbering
        assert!(pdb_output.contains('1'));
    }

    #[test]
    fn pdb_str_empty_produces_error() {
        let result = pdb_str_to_entities("");
        assert!(result.is_err());
    }

    #[test]
    fn sanitize_pdb_fixes_remark_lines() {
        let input = "REMARK some info\nATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N\nEND\n";
        let sanitized = sanitize_pdb(input);
        assert!(sanitized.contains("REMARK   0"));
    }

    #[test]
    fn name_to_bytes_pads_with_spaces() {
        let result: [u8; 4] = name_to_bytes::<4>("CA");
        assert_eq!(result, [b'C', b'A', b' ', b' ']);
    }

    #[test]
    fn name_to_bytes_truncates_long_names() {
        let result: [u8; 3] = name_to_bytes::<3>("ALAAA");
        assert_eq!(result, [b'A', b'L', b'A']);
    }
}

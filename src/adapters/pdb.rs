//! PDB format parsing and writing.

use std::fmt::Write as _;
use std::io::BufReader;

use pdbtbx::{
    ContainsAtomConformer, ContainsAtomConformerResidue,
    ContainsAtomConformerResidueChain, ReadOptions, StrictnessLevel,
};

use crate::assembly::Assembly;
use crate::entity::molecule::atom::Atom;
use crate::entity::molecule::polymer::Residue;
use crate::entity::molecule::MoleculeEntity;
use crate::ops::codec::{split_into_entities, AdapterError, Coords};

// ---------------------------------------------------------------------------
// Entity-first API (primary)
// ---------------------------------------------------------------------------

/// Parse PDB format string to entity list.
///
/// # Errors
///
/// Returns [`AdapterError`] if parsing fails.
pub fn pdb_str_to_entities(
    pdb_str: &str,
) -> Result<Vec<MoleculeEntity>, AdapterError> {
    parse_pdb_to_entities(pdb_str)
}

/// Load PDB file to entity list.
///
/// Sanitizes non-standard lines (e.g. GROMACS/MemProtMD output) before
/// parsing.
///
/// # Errors
///
/// Returns [`AdapterError`] if the file cannot be read or parsing fails.
pub fn pdb_file_to_entities(
    path: &std::path::Path,
) -> Result<Vec<MoleculeEntity>, AdapterError> {
    let content = std::fs::read_to_string(path).map_err(|e| {
        AdapterError::PdbParseError(format!("Failed to read file: {e}"))
    })?;
    let sanitized = sanitize_pdb(&content);
    parse_pdb_to_entities(&sanitized)
}

/// Load structure file (PDB or mmCIF, detected by extension) to entity list.
///
/// # Errors
///
/// Returns [`AdapterError`] if the file cannot be read or parsing fails.
pub fn structure_file_to_entities(
    path: &std::path::Path,
) -> Result<Vec<MoleculeEntity>, AdapterError> {
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
// PDB writers
// ---------------------------------------------------------------------------

/// Emit an `Assembly` as a PDB-format string. All entities are flattened
/// into a single atom stream in declaration order.
#[must_use]
pub fn assembly_to_pdb(assembly: &Assembly) -> String {
    entities_to_pdb(assembly.entities())
}

/// Emit an entity slice as a PDB-format string.
#[must_use]
pub fn entities_to_pdb<E: std::borrow::Borrow<MoleculeEntity>>(
    entities: &[E],
) -> String {
    let mut pdb = String::new();
    let mut serial: usize = 0;
    for entity in entities {
        write_entity_atoms(entity.borrow(), &mut serial, &mut pdb);
    }
    pdb.push_str("END\n");
    pdb
}

fn write_entity_atoms(
    entity: &MoleculeEntity,
    serial: &mut usize,
    out: &mut String,
) {
    match entity {
        MoleculeEntity::Protein(e) => {
            write_polymer_atoms(
                &e.atoms,
                &e.residues,
                e.pdb_chain_id,
                serial,
                out,
            );
        }
        MoleculeEntity::NucleicAcid(e) => {
            write_polymer_atoms(
                &e.atoms,
                &e.residues,
                e.pdb_chain_id,
                serial,
                out,
            );
        }
        MoleculeEntity::SmallMolecule(e) => {
            for atom in &e.atoms {
                *serial += 1;
                write_atom_line(
                    *serial,
                    atom,
                    ResidueCtx {
                        chain_id: b' ',
                        res_name: e.residue_name,
                        res_num: 1,
                    },
                    out,
                );
            }
        }
        MoleculeEntity::Bulk(e) =>
        {
            #[allow(
                clippy::cast_possible_truncation,
                clippy::cast_possible_wrap,
                reason = "atom count fits in i32 for valid structures"
            )]
            for (i, atom) in e.atoms.iter().enumerate() {
                *serial += 1;
                write_atom_line(
                    *serial,
                    atom,
                    ResidueCtx {
                        chain_id: b' ',
                        res_name: e.residue_name,
                        res_num: (i as i32) + 1,
                    },
                    out,
                );
            }
        }
    }
}

fn write_polymer_atoms(
    atoms: &[Atom],
    residues: &[Residue],
    chain_id: u8,
    serial: &mut usize,
    out: &mut String,
) {
    for residue in residues {
        for idx in residue.atom_range.clone() {
            *serial += 1;
            write_atom_line(
                *serial,
                &atoms[idx],
                ResidueCtx {
                    chain_id,
                    res_name: residue.name,
                    res_num: residue.number,
                },
                out,
            );
        }
    }
}

#[derive(Copy, Clone)]
struct ResidueCtx {
    chain_id: u8,
    res_name: [u8; 3],
    res_num: i32,
}

fn write_atom_line(
    serial: usize,
    atom: &Atom,
    ctx: ResidueCtx,
    out: &mut String,
) {
    let atom_name = std::str::from_utf8(&atom.name).unwrap_or("X   ");
    let res_name_str = std::str::from_utf8(&ctx.res_name).unwrap_or("UNK");
    let _ = writeln!(
        out,
        "ATOM  {:>5} {:<4} {:>3} {}{:>4}    \
         {:>8.3}{:>8.3}{:>8.3}{:>6.2}{:>6.2}",
        serial,
        atom_name,
        res_name_str,
        ctx.chain_id as char,
        ctx.res_num,
        atom.position.x,
        atom.position.y,
        atom.position.z,
        atom.occupancy,
        atom.b_factor,
    );
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
) -> Result<Vec<MoleculeEntity>, AdapterError> {
    let coords = parse_pdb_to_coords(input)?;
    Ok(split_into_entities(&coords))
}

/// Parse PDB string into flat Coords (internal, used by entity pipeline).
fn parse_pdb_to_coords(input: &str) -> Result<Coords, AdapterError> {
    use crate::element::Element;
    use crate::ops::codec::{ChainIdMapper, CoordsAtom};

    let reader = BufReader::new(input.as_bytes());

    let (pdb, _errors) = ReadOptions::new()
        .set_format(pdbtbx::Format::Pdb)
        .set_level(StrictnessLevel::Loose)
        .read_raw(reader)
        .map_err(|errs| {
            AdapterError::PdbParseError(
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
        return Err(AdapterError::PdbParseError(
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
    fn entities_to_pdb_produces_valid_output() {
        let entities = pdb_str_to_entities(MINIMAL_PDB).unwrap();
        let pdb_output = entities_to_pdb(&entities);

        assert!(pdb_output.contains("ATOM"));
        assert!(pdb_output.contains("ALA"));
        assert!(pdb_output.ends_with("END\n"));
        // Should have 4 ATOM lines + END
        let atom_line_count =
            pdb_output.lines().filter(|l| l.starts_with("ATOM")).count();
        assert_eq!(atom_line_count, 4);
    }

    #[test]
    fn entities_to_pdb_preserves_coordinates_in_output() {
        let entities = pdb_str_to_entities(MINIMAL_PDB).unwrap();
        let pdb_output = entities_to_pdb(&entities);

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

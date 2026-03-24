//! Residue classification: mapping residue names to molecule types.

use super::MoleculeType;

/// Standard amino acid residue names.
pub const PROTEIN_RESIDUES: &[&str] = &[
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    // Non-standard but protein-like
    "MSE", "SEC", "PYL",
];

/// Standard DNA residue names (mmCIF convention + common alias THY).
const DNA_RESIDUES: &[&str] = &["DA", "DC", "DG", "DT", "DU", "DI", "THY"];

/// Standard RNA residue names.
const RNA_RESIDUES: &[&str] = &[
    "A", "C", "G", "U", "ADE", "CYT", "GUA", "URA", "I", "RAD", "RCY", "RGU",
];

/// Water residue names.
pub(crate) const WATER_RESIDUES: &[&str] = &[
    "HOH", "WAT", "H2O", "DOD",
    // MD simulation water models (GROMACS, AMBER, CHARMM, etc.)
    "SOL", "TIP", "TP3", "TIP3", "T3P", "SPC", "TP4", "TIP4", "T4P", "TP5",
    "TIP5",
];

/// Known ion residue names.
const ION_RESIDUES: &[&str] = &[
    "ZN", "MG", "NA", "CL", "FE", "MN", "CO", "NI", "CU", "K", "CA", "BR", "I",
    "F", "LI", "CD", "SR", "BA", "CS", "RB", "PB", "HG", "PT", "AU", "AG",
];

/// Known lipid residue 3-char truncated names.
const LIPID_RESIDUES: &[&str] = &[
    "DPP", "POP", "DOP", "DMP", "DSP", "DLP", "PPE", "DPE", "PPG", "DPG",
    "PPS", "DPS", "CHO", "CHL", "SPH", "CER", "PAL", "OLE", "STE", "MYR",
    "LAU", "LHG", "LMG", "DGD", "SQD", "LMT", "HTG",
];

/// Known cofactor residue names.
pub(crate) const COFACTOR_RESIDUES: &[&str] = &[
    "HEM", "HEC", "HEA", "HEB", "CLA", "CHL", "PHO", "BCR", "BCB", "PL9",
    "PLQ", "UQ1", "UQ2", "MQ7", "NAD", "NAP", "NAI", "NDP", "FAD", "FMN",
    "ATP", "ADP", "AMP", "ANP", "GTP", "GDP", "GMP", "GNP", "SAM", "SAH",
    "COA", "ACO", "PLP", "PMP", "TPP", "TDP", "BTN", "BIO", "H4B", "BH4",
    "SF4", "FES", "F3S",
];

/// Known solvent / crystallization artifact residue names.
const SOLVENT_RESIDUES: &[&str] = &[
    "GOL", "EDO", "PEG", "1PE", "P6G", "PG4", "PGE", "SO4", "SUL", "PO4",
    "ACT", "ACE", "CIT", "FMT", "TRS", "MES", "EPE", "IMD", "MPD", "DMS",
    "BME", "IPA", "EOH",
];

/// Human-readable display name for a cofactor residue code.
fn cofactor_display_name(res_name: &str) -> &str {
    match res_name {
        "CLA" => "Chlorophyll A",
        "CHL" => "Chlorophyll B",
        "BCR" => "Beta-Carotene",
        "BCB" => "Beta-Carotene B",
        "HEM" => "Heme",
        "HEC" => "Heme C",
        "HEA" => "Heme A",
        "HEB" => "Heme B",
        "PHO" => "Pheophytin",
        "PL9" | "PLQ" => "Plastoquinone",
        "UQ1" | "UQ2" => "Ubiquinone",
        "MQ7" => "Menaquinone",
        "NAD" | "NAP" | "NAI" | "NDP" => "NAD",
        "SAM" | "SAH" => "SAM/SAH",
        "COA" | "ACO" => "Coenzyme A",
        "PLP" | "PMP" => "PLP",
        "TPP" | "TDP" => "Thiamine PP",
        "BTN" | "BIO" => "Biotin",
        "H4B" | "BH4" => "Tetrahydrobiopterin",
        "SF4" => "[4Fe-4S] Cluster",
        "FES" => "[2Fe-2S] Cluster",
        "F3S" => "[3Fe-4S] Cluster",
        _ => res_name,
    }
}

/// Display name for a small molecule, dispatching by molecule type.
pub(super) fn small_molecule_display_name(
    mol_type: MoleculeType,
    res_name: &str,
) -> String {
    match mol_type {
        MoleculeType::Cofactor => cofactor_display_name(res_name).to_owned(),
        _ => res_name.to_owned(),
    }
}

/// Classify a residue name into a `MoleculeType`.
///
/// The name should be trimmed of whitespace before calling.
#[must_use]
pub fn classify_residue(name: &str) -> MoleculeType {
    if PROTEIN_RESIDUES.contains(&name) {
        return MoleculeType::Protein;
    }
    if WATER_RESIDUES.contains(&name) {
        return MoleculeType::Water;
    }
    if DNA_RESIDUES.contains(&name) {
        return MoleculeType::DNA;
    }
    if RNA_RESIDUES.contains(&name) {
        return MoleculeType::RNA;
    }
    if ION_RESIDUES.contains(&name) {
        return MoleculeType::Ion;
    }
    if COFACTOR_RESIDUES.contains(&name) {
        return MoleculeType::Cofactor;
    }
    if SOLVENT_RESIDUES.contains(&name) {
        return MoleculeType::Solvent;
    }
    let truncated = if name.len() > 3 { &name[..3] } else { name };
    if LIPID_RESIDUES.contains(&truncated) {
        return MoleculeType::Lipid;
    }
    MoleculeType::Ligand
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_classify_protein() {
        assert_eq!(classify_residue("ALA"), MoleculeType::Protein);
        assert_eq!(classify_residue("GLY"), MoleculeType::Protein);
        assert_eq!(classify_residue("MSE"), MoleculeType::Protein);
    }

    #[test]
    fn test_classify_nucleic() {
        assert_eq!(classify_residue("DA"), MoleculeType::DNA);
        assert_eq!(classify_residue("DT"), MoleculeType::DNA);
        assert_eq!(classify_residue("A"), MoleculeType::RNA);
        assert_eq!(classify_residue("U"), MoleculeType::RNA);
    }

    #[test]
    fn test_classify_water_ion_ligand() {
        assert_eq!(classify_residue("HOH"), MoleculeType::Water);
        assert_eq!(classify_residue("WAT"), MoleculeType::Water);
        assert_eq!(classify_residue("ZN"), MoleculeType::Ion);
        assert_eq!(classify_residue("MG"), MoleculeType::Ion);
        assert_eq!(classify_residue("ATP"), MoleculeType::Cofactor);
        assert_eq!(classify_residue("HEM"), MoleculeType::Cofactor);
        assert_eq!(classify_residue("UNL"), MoleculeType::Ligand);
    }

    #[test]
    fn test_classify_cofactor() {
        assert_eq!(classify_residue("CLA"), MoleculeType::Cofactor);
        assert_eq!(classify_residue("HEM"), MoleculeType::Cofactor);
        assert_eq!(classify_residue("FAD"), MoleculeType::Cofactor);
        assert_eq!(classify_residue("NAD"), MoleculeType::Cofactor);
        assert_eq!(classify_residue("SF4"), MoleculeType::Cofactor);
        assert_eq!(classify_residue("BCR"), MoleculeType::Cofactor);
        assert_eq!(classify_residue("PL9"), MoleculeType::Cofactor);
    }

    #[test]
    fn test_classify_solvent() {
        assert_eq!(classify_residue("GOL"), MoleculeType::Solvent);
        assert_eq!(classify_residue("EDO"), MoleculeType::Solvent);
        assert_eq!(classify_residue("SO4"), MoleculeType::Solvent);
        assert_eq!(classify_residue("PEG"), MoleculeType::Solvent);
        assert_eq!(classify_residue("MPD"), MoleculeType::Solvent);
        assert_eq!(classify_residue("DMS"), MoleculeType::Solvent);
    }
}

//! Fast-path mmCIF parser that goes directly to Coords, skipping the DOM.
//!
//! Reads the `_atom_site` loop inline and builds flat column vectors without
//! intermediate String allocations. Falls back to the DOM path for anything
//! the fast path can't handle.

use crate::element::Element;
use crate::ops::codec::{AdapterError, ChainIdMapper, Coords, CoordsAtom};

/// Column indices within the `_atom_site` loop.
struct AtomSiteColumns {
    label_atom_id: Option<usize>,
    label_comp_id: Option<usize>,
    label_asym_id: Option<usize>,
    label_seq_id: Option<usize>,
    cartn_x: Option<usize>,
    cartn_y: Option<usize>,
    cartn_z: Option<usize>,
    type_symbol: Option<usize>,
    occupancy: Option<usize>,
    b_iso: Option<usize>,
    ncols: usize,
}

impl AtomSiteColumns {
    fn new() -> Self {
        Self {
            label_atom_id: None,
            label_comp_id: None,
            label_asym_id: None,
            label_seq_id: None,
            cartn_x: None,
            cartn_y: None,
            cartn_z: None,
            type_symbol: None,
            occupancy: None,
            b_iso: None,
            ncols: 0,
        }
    }

    fn map_tag(&mut self, tag: &str, index: usize) {
        match tag.to_ascii_lowercase().as_str() {
            "_atom_site.label_atom_id" => self.label_atom_id = Some(index),
            "_atom_site.label_comp_id" => self.label_comp_id = Some(index),
            "_atom_site.label_asym_id" => self.label_asym_id = Some(index),
            "_atom_site.label_seq_id" => self.label_seq_id = Some(index),
            "_atom_site.cartn_x" => self.cartn_x = Some(index),
            "_atom_site.cartn_y" => self.cartn_y = Some(index),
            "_atom_site.cartn_z" => self.cartn_z = Some(index),
            "_atom_site.type_symbol" => self.type_symbol = Some(index),
            "_atom_site.occupancy" => self.occupancy = Some(index),
            "_atom_site.b_iso_or_equiv" => self.b_iso = Some(index),
            _ => {}
        }
        self.ncols = self.ncols.max(index + 1);
    }

    fn is_valid(&self) -> bool {
        self.label_atom_id.is_some()
            && self.label_comp_id.is_some()
            && self.label_asym_id.is_some()
            && self.cartn_x.is_some()
            && self.cartn_y.is_some()
            && self.cartn_z.is_some()
    }
}

/// Parse mmCIF text directly to Coords, bypassing the DOM.
///
/// Returns `None` if the fast path can't handle the input (no `_atom_site`
/// loop found). The caller should fall back to the DOM path.
#[allow(clippy::too_many_lines)]
pub(crate) fn parse_mmcif_fast(
    input: &str,
) -> Option<Result<Coords, AdapterError>> {
    let bytes = input.as_bytes();
    let len = bytes.len();
    let mut pos = 0;

    // Scan for the _atom_site loop
    let cols = find_atom_site_loop(bytes, &mut pos)?;

    // Now read the atom data rows directly into vectors
    let mut atoms = Vec::new();
    let mut chain_ids = Vec::new();
    let mut res_names = Vec::new();
    let mut res_nums = Vec::new();
    let mut atom_names = Vec::new();
    let mut elements = Vec::new();
    let mut chain_mapper = ChainIdMapper::new();

    loop {
        skip_whitespace_and_comments(bytes, &mut pos);
        if pos >= len {
            break;
        }

        // Check if we hit a new keyword (loop_, data_, save_, or a tag _)
        if bytes[pos] == b'_'
            || (pos + 4 < len
                && bytes[pos..pos + 5].eq_ignore_ascii_case(b"loop_"))
            || (pos + 4 < len
                && bytes[pos..pos + 5].eq_ignore_ascii_case(b"data_"))
            || (pos + 4 < len
                && bytes[pos..pos + 5].eq_ignore_ascii_case(b"save_"))
        {
            break;
        }

        // Read ncols values for this row
        let mut row_vals: Vec<&str> = Vec::with_capacity(cols.ncols);
        for _ in 0..cols.ncols {
            skip_whitespace_and_comments(bytes, &mut pos);
            if pos >= len {
                return Some(Err(AdapterError::InvalidFormat(
                    "Unexpected end of CIF data in _atom_site loop".into(),
                )));
            }
            row_vals.push(scan_value(input, bytes, &mut pos));
        }

        // Extract fields from row
        let x = parse_cif_float(row_vals[cols.cartn_x?]);
        let y = parse_cif_float(row_vals[cols.cartn_y?]);
        let z = parse_cif_float(row_vals[cols.cartn_z?]);

        let (Some(x), Some(y), Some(z)) = (x, y, z) else {
            continue; // skip atoms with inapplicable coords
        };

        #[allow(clippy::cast_possible_truncation)]
        atoms.push(CoordsAtom {
            x: x as f32,
            y: y as f32,
            z: z as f32,
            occupancy: cols
                .occupancy
                .and_then(|i| parse_cif_float(row_vals[i]))
                .unwrap_or(1.0) as f32,
            b_factor: cols
                .b_iso
                .and_then(|i| parse_cif_float(row_vals[i]))
                .unwrap_or(0.0) as f32,
        });

        let asym_id = row_vals[cols.label_asym_id?];
        chain_ids.push(chain_mapper.get_or_assign(asym_id));

        let comp_id = row_vals[cols.label_comp_id?];
        res_names.push(name_to_bytes::<3>(comp_id));

        let seq_id = cols
            .label_seq_id
            .and_then(|i| row_vals[i].parse::<i32>().ok())
            .unwrap_or(0);
        res_nums.push(seq_id);

        let atom_id = row_vals[cols.label_atom_id?];
        atom_names.push(name_to_bytes::<4>(atom_id));

        let elem = cols.type_symbol.map_or_else(
            || Element::from_atom_name(atom_id),
            |i| {
                let s = row_vals[i];
                if s == "." || s == "?" || s.is_empty() {
                    Element::from_atom_name(atom_id)
                } else {
                    Element::from_symbol(s)
                }
            },
        );
        elements.push(elem);
    }

    if atoms.is_empty() {
        return Some(Err(AdapterError::InvalidFormat(
            "No atoms found in CIF".into(),
        )));
    }

    Some(Ok(Coords {
        num_atoms: atoms.len(),
        atoms,
        chain_ids,
        res_names,
        res_nums,
        atom_names,
        elements,
    }))
}

// ---------------------------------------------------------------------------
// Inline helpers (no allocations)
// ---------------------------------------------------------------------------

fn skip_whitespace_and_comments(bytes: &[u8], pos: &mut usize) {
    while *pos < bytes.len() {
        match bytes[*pos] {
            b' ' | b'\t' | b'\r' | b'\n' => *pos += 1,
            b'#' => {
                while *pos < bytes.len() && bytes[*pos] != b'\n' {
                    *pos += 1;
                }
            }
            _ => break,
        }
    }
}

/// Scan an unquoted token, returning a &str slice of the input.
fn scan_unquoted_token<'a>(bytes: &'a [u8], pos: &mut usize) -> &'a str {
    let start = *pos;
    while *pos < bytes.len() {
        let c = bytes[*pos];
        if c.is_ascii_whitespace() || c == b'#' {
            break;
        }
        *pos += 1;
    }
    // SAFETY: CIF is text, input was &str, so this slice is valid UTF-8
    unsafe { std::str::from_utf8_unchecked(&bytes[start..*pos]) }
}

/// Scan a value (handles quoting), returning a &str slice from `input`.
fn scan_value<'a>(input: &'a str, bytes: &'a [u8], pos: &mut usize) -> &'a str {
    let b = bytes[*pos];

    // Quoted string
    if b == b'\'' || b == b'"' {
        let quote = b;
        *pos += 1;
        let start = *pos;
        loop {
            if *pos >= bytes.len() {
                break;
            }
            if bytes[*pos] == quote {
                let end = *pos;
                *pos += 1;
                return &input[start..end];
            }
            *pos += 1;
        }
        return &input[start..*pos];
    }

    // Unquoted
    scan_unquoted_token(bytes, pos)
}

/// Parse a CIF float, handling uncertainty notation like `50.123(4)`.
fn parse_cif_float(s: &str) -> Option<f64> {
    if s == "." || s == "?" {
        return None;
    }
    let s = s.find('(').map_or(s, |idx| &s[..idx]);
    s.parse().ok()
}

fn name_to_bytes<const N: usize>(name: &str) -> [u8; N] {
    let mut buf = [b' '; N];
    for (i, b) in name.bytes().take(N).enumerate() {
        buf[i] = b;
    }
    buf
}

/// Scan forward to find the `_atom_site` loop and collect column mappings.
fn find_atom_site_loop(
    bytes: &[u8],
    pos: &mut usize,
) -> Option<AtomSiteColumns> {
    let len = bytes.len();
    loop {
        if *pos >= len {
            return None;
        }
        skip_whitespace_and_comments(bytes, pos);
        if *pos >= len {
            return None;
        }

        let token = scan_unquoted_token(bytes, pos);
        if !token.eq_ignore_ascii_case("loop_") {
            continue;
        }

        // Check if the next tag starts with _atom_site
        skip_whitespace_and_comments(bytes, pos);
        if *pos >= len || bytes[*pos] != b'_' {
            continue;
        }

        let peek_start = *pos;
        let peek = scan_unquoted_token(bytes, pos);
        if !peek.to_ascii_lowercase().starts_with("_atom_site.") {
            *pos = peek_start;
            continue;
        }

        // Found the atom_site loop — collect all tags
        let mut cols = AtomSiteColumns::new();
        cols.map_tag(peek, 0);
        let mut tag_idx = 1;

        loop {
            skip_whitespace_and_comments(bytes, pos);
            if *pos >= len || bytes[*pos] != b'_' {
                break;
            }
            let tag = scan_unquoted_token(bytes, pos);
            cols.map_tag(tag, tag_idx);
            tag_idx += 1;
        }

        return if cols.is_valid() { Some(cols) } else { None };
    }
}

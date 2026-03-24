//! DSSP secondary structure classification from hydrogen bond pairs.
//!
//! Takes H-bond pairs (from [`crate::analysis::bonds::hydrogen`]) and
//! classifies residues into Helix, Sheet, or Coil using DSSP pattern
//! matching rules.

use crate::analysis::bonds::hydrogen::HBond;
use crate::analysis::SSType;

/// Classify secondary structure from H-bond pairs.
///
/// Detects helices (α i→i+4, 3₁₀ i→i+3, π i→i+5 turn patterns) and
/// sheets (parallel/antiparallel bridge patterns) using standard DSSP
/// rules.
#[must_use]
#[allow(
    clippy::too_many_lines,
    reason = "DSSP helix + sheet detection is inherently sequential"
)]
pub fn classify(hbonds: &[HBond], n_residues: usize) -> Vec<SSType> {
    if n_residues < 2 {
        return vec![SSType::Coil; n_residues];
    }

    let has_hbond = |donor: usize, acceptor: usize| -> bool {
        hbonds
            .iter()
            .any(|h| h.donor == donor && h.acceptor == acceptor)
    };

    let mut raw_ss = vec![SSType::Coil; n_residues];

    // Detect helices: n-turn at residue i means H-bond from i+n to i.
    for turn_size in [4usize, 3, 5] {
        let mut consecutive_turns = 0usize;
        let min_consecutive = if turn_size == 4 { 4 } else { 3 };

        for i in 0..n_residues {
            if !(i + turn_size < n_residues && has_hbond(i + turn_size, i)) {
                consecutive_turns = 0;
                continue;
            }

            consecutive_turns += 1;
            if consecutive_turns < min_consecutive {
                continue;
            }

            let start = if consecutive_turns == min_consecutive {
                i + 1 - (min_consecutive - 1)
            } else {
                i
            };
            for ss_entry in raw_ss
                .iter_mut()
                .take(n_residues)
                .skip(start)
                .take(i + turn_size + 1 - start)
            {
                *ss_entry = SSType::Helix;
            }
        }
    }

    // Detect sheets: bridge patterns between non-adjacent residues.
    for i in 1..n_residues.saturating_sub(1) {
        if raw_ss[i] == SSType::Helix {
            continue;
        }
        for j in (i + 2)..n_residues {
            if raw_ss[j] == SSType::Helix {
                continue;
            }

            let parallel = (i > 0
                && j + 1 < n_residues
                && has_hbond(i, j)
                && has_hbond(j + 1, i))
                || (j > 0
                    && i + 1 < n_residues
                    && has_hbond(j, i)
                    && has_hbond(i + 1, j));

            let antiparallel = (has_hbond(i, j) && has_hbond(j, i))
                || (i > 0
                    && j + 1 < n_residues
                    && j > 0
                    && i + 1 < n_residues
                    && has_hbond(i, j)
                    && has_hbond(j, i));

            if parallel || antiparallel {
                raw_ss[i] = SSType::Sheet;
                raw_ss[j] = SSType::Sheet;
            }
        }
    }

    raw_ss
}

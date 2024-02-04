// Copyright Ryangguk Kim @ Oak Bioinformatics, LLC
// 
// This software is available under a dual licensing model, offering users the choice between the Affero General Public License version 3 (AGPL-3) for open-source use and a commercial license for proprietary or commercial use. 
// 
// To obtain a commercial license, please contact info@oakbioinformatics.com.

//! Example:
//!
//! ```
//! use anyhow::Result;
//! 
//! fn example() -> Result<String> {
//!     let s = spdi::SPDI::new("path/to/2bit/file")?;
//!     let spdi_string = s.get_spdi_string("chr1", 99092, "C", "CT")?;
//!     println!("{}", spdi_string);
//! }
//! ```

use anyhow::Result;
use noodles::vcf::record::reference_bases::Base;
mod error;
mod util;
mod trim;
mod grow;
mod tests;
use trim::{trim_left, trim_right};
use grow::Grower;
use util::{get_bases_of_string, get_string_of_bases};

pub struct SPDI {
    grower: Grower,
}

impl SPDI {
    pub fn new(twobit_path: &str) -> Result<SPDI> {
        let grower = Grower::new(twobit_path)?;
        Ok(SPDI { grower })
    }

    pub fn get_spdi_conversion(
        &mut self,
        chrom: &str,
        pos: usize,
        ref_bases: &Vec<Base>,
        alt_bases: &Vec<Base>,
    ) -> Result<(usize, Vec<Base>, Vec<Base>)> {
        let ref_bases_len = ref_bases.len();
        let alt_bases_len = alt_bases.len();
        let ref_start = 0;
        let ref_end = match ref_bases_len {
            0 => 0,
            _ => ref_bases_len,
        };
        let alt_start = 0;
        let alt_end = match alt_bases_len {
            0 => 0,
            _ => alt_bases_len,
        };
        let new_ref_start: usize;
        let new_ref_end: usize;
        let new_alt_start: usize;
        let new_alt_end: usize;
        // shrink
        (new_ref_end, new_alt_end) =
            trim_right(ref_bases, alt_bases, ref_start, ref_end, alt_start, alt_end);
        (new_ref_start, new_alt_start) = trim_left(
            ref_bases,
            alt_bases,
            ref_start,
            new_ref_end,
            alt_start,
            new_alt_end,
        );
        let shrunk_ref_bases: Vec<Base>;
        let shrunk_alt_bases: Vec<Base>;
        shrunk_ref_bases = ref_bases[new_ref_start..new_ref_end].to_vec();
        shrunk_alt_bases = alt_bases[new_alt_start..new_alt_end].to_vec();
        let shrunk_pos = pos + new_ref_start;
        let shrunk_ref_bases_len = shrunk_ref_bases.len();
        let shrunk_alt_bases_len = shrunk_alt_bases.len();
        match shrunk_ref_bases_len {
            0 => {
                match shrunk_alt_bases_len {
                    // same
                    0 => {
                        let base = ref_bases[0..1].to_vec();
                        Ok((pos, base.clone(), base))
                    }
                    // insertion
                    _ => self.grower.grow_pos_ref_alt(
                        chrom,
                        shrunk_pos,
                        &shrunk_ref_bases,
                        &shrunk_alt_bases,
                    ),
                }
            }
            _ => {
                match shrunk_alt_bases_len {
                    // deletion
                    0 => self.grower.grow_pos_ref_alt(
                        chrom,
                        shrunk_pos,
                        &shrunk_ref_bases,
                        &shrunk_alt_bases,
                    ),
                    // ambiguous
                    _ => Ok((pos, shrunk_ref_bases, shrunk_alt_bases)),
                }
            }
        }
    }

    pub fn get_spdi_string_components(
        &mut self,
        chrom: &str,
        pos: usize,
        ref_bases: &Vec<Base>,
        alt_bases: &Vec<Base>,
    ) -> Result<(usize, String, String)> {
        let new_pos: usize;
        let new_ref_bases: Vec<Base>;
        let new_alt_bases: Vec<Base>;
        (new_pos, new_ref_bases, new_alt_bases) =
            self.get_spdi_conversion(chrom, pos, ref_bases, alt_bases)?;
        let new_ref_bases_string = get_string_of_bases(&new_ref_bases);
        let new_alt_bases_string = get_string_of_bases(&new_alt_bases);
        Ok((new_pos, new_ref_bases_string, new_alt_bases_string))
    }

    pub fn get_spdi_string(
        &mut self,
        chrom: &str,
        pos: usize,
        ref_bases_s: &str,
        alt_bases_s: &str,
    ) -> Result<String> {
        let ref_bases = get_bases_of_string(ref_bases_s)?;
        let alt_bases = get_bases_of_string(alt_bases_s)?;
        let new_pos: usize;
        let new_ref_bases_s: String;
        let new_alt_bases_s: String;
        (new_pos, new_ref_bases_s, new_alt_bases_s) =
            self.get_spdi_string_components(chrom, pos, &ref_bases, &alt_bases)?;
        Ok(format!(
            "{}:{}:{}:{}",
            chrom, new_pos, new_ref_bases_s, new_alt_bases_s
        ))
    }
}

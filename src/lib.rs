// Copyright Ryangguk Kim @ Oak Bioinformatics, LLC
//
// This software is available under a dual licensing model, offering users the choice between the Affero General Public License version 3 (AGPL-3) for open-source use and a commercial license for proprietary or commercial use.
//
// To obtain a commercial license, please contact info@oakbioinformatics.com.

//! Example:
//!
//! ```
//! use spdi::error::Error;
//!
//! fn example() -> Result<(), Error> {
//!     let mut s = spdi::SPDI::new("path/to/2bit/file")?;
//!     let spdi_string = s.get_spdi_string("chr1".as_bytes(), 99092, "C".as_bytes(),
//!     "CT".as_bytes())?;
//!     assert_eq!(spdi_string, "chr1:99092:C:CT");
//!     Ok(())
//! }
//! ```

pub mod error;
mod grow;
mod tests;
mod trim;
pub mod util;
use grow::Grower;
pub use noodles::vcf;
pub type Base = vcf::record::reference_bases::base::Base;
use trim::{trim_left, trim_right};
use util::{get_bases_of_vu8, get_string_of_bases};
use error::Error;
use std::path::PathBuf;

pub struct SPDI {
    grower: Grower,
}

impl SPDI {
    pub fn new(twobit_path: &PathBuf) -> std::result::Result<SPDI, Error> {
        let grower = Grower::new(twobit_path)?;
        Ok(SPDI { grower })
    }

    pub fn get_spdi_conversion_str(
        &mut self,
        chrom: &[u8],
        pos: usize,
        ref_bases: &[u8],
        alt_bases: &[u8],
    ) -> std::result::Result<(usize, Box<[Base]>, Box<[Base]>), Error> {
        let ref_bases_v: Vec<Base>;
        let alt_bases_v: Vec<Base>;
        match util::get_bases_of_vu8(ref_bases) {
            Ok(v) => ref_bases_v = v,
            Err(e) => {
                return Err(e);
            }
        }
        let alt_bases_q: &[u8];
        if alt_bases[0] == '.' as u8 {
            alt_bases_q = ref_bases;
        } else {
            alt_bases_q = alt_bases;
        }
        match util::get_bases_of_vu8(alt_bases_q) {
            Ok(v) => alt_bases_v = v,
            Err(e) => {
                return Err(e);
            }
        }
        self.get_spdi_conversion(chrom, pos, &ref_bases_v, &alt_bases_v)
    }

    pub fn get_spdi_conversion(
        &mut self,
        chrom: &[u8],
        pos: usize,
        ref_bases: &[Base],
        alt_bases: &[Base],
    ) -> std::result::Result<(usize, Box<[Base]>, Box<[Base]>), Error> {
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
        let shrunk_ref_bases = &ref_bases[new_ref_start..new_ref_end];
        let shrunk_alt_bases = &alt_bases[new_alt_start..new_alt_end];
        let shrunk_pos = pos + new_ref_start;
        let shrunk_ref_bases_len = shrunk_ref_bases.len();
        let shrunk_alt_bases_len = shrunk_alt_bases.len();
        match shrunk_ref_bases_len {
            0 => {
                match shrunk_alt_bases_len {
                    // same
                    0 => {
                        let base = &ref_bases[0..1];
                        Ok((pos, base.to_vec().into_boxed_slice(), base.to_vec().into_boxed_slice()))
                    }
                    // insertion
                    _ => {
                        match self.grower.grow(
                            chrom,
                            shrunk_pos,
                            shrunk_ref_bases,
                            shrunk_alt_bases,
                        ) {
                            Ok(v) => Ok(v),
                            Err(e) => {
                                eprintln!(
                                    "Error: {}. {}:{}:{:?}:{:?}",
                                    e, std::str::from_utf8(chrom).unwrap(), pos, ref_bases, alt_bases
                                );
                                Err(e)
                            }
                        }
                    }
                }
            }
            _ => {
                match shrunk_alt_bases_len {
                    // deletion
                    0 => match self.grower.grow(
                        chrom,
                        shrunk_pos,
                        shrunk_ref_bases,
                        shrunk_alt_bases,
                    ) {
                        Ok(v) => Ok(v),
                        Err(e) => {
                            eprintln!(
                                "Error: {}. {}:{}:{:?}:{:?}",
                                e, std::str::from_utf8(chrom).unwrap(), pos, ref_bases, alt_bases
                            );
                            Err(e)
                        }
                    },
                    // ambiguous
                    _ => Ok((pos, shrunk_ref_bases.to_vec().into_boxed_slice(), shrunk_alt_bases.to_vec().into_boxed_slice())),
                }
            }
        }
    }

    pub fn get_spdi_string_components(
        &mut self,
        chrom: &[u8],
        pos: usize,
        ref_bases: &Vec<Base>,
        alt_bases: &Vec<Base>,
    ) -> Result<(usize, String, String), Error> {
        let new_pos: usize;
        let new_ref_bases: Box<[Base]>;
        let new_alt_bases: Box<[Base]>;
        (new_pos, new_ref_bases, new_alt_bases) =
            self.get_spdi_conversion(chrom, pos, ref_bases, alt_bases)?;
        let new_ref_bases_string = get_string_of_bases(&new_ref_bases);
        let new_alt_bases_string = get_string_of_bases(&new_alt_bases);
        Ok((new_pos, new_ref_bases_string, new_alt_bases_string))
    }

    pub fn get_spdi_string(
        &mut self,
        chrom: &[u8],
        pos: usize,
        ref_bases_s: &[u8],
        alt_bases_s: &[u8],
    ) -> Result<String, Error> {
        let ref_bases = get_bases_of_vu8(ref_bases_s)?;
        let alt_bases = get_bases_of_vu8(alt_bases_s)?;
        let new_pos: usize;
        let new_ref_bases_s: String;
        let new_alt_bases_s: String;
        (new_pos, new_ref_bases_s, new_alt_bases_s) =
            self.get_spdi_string_components(chrom, pos, &ref_bases, &alt_bases)?;
        Ok(format!(
            "{}:{}:{}:{}",
            std::str::from_utf8(chrom).unwrap(), new_pos, new_ref_bases_s, new_alt_bases_s
        ))
    }
}


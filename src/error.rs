// Copyright Ryangguk Kim @ Oak Bioinformatics, LLC
//
// This software is available under a dual licensing model, offering users the choice between the Affero General Public License version 3 (AGPL-3) for open-source use and a commercial license for proprietary or commercial use.
//
// To obtain a commercial license, please contact info@oakbioinformatics.com.

#[allow(dead_code)]
#[derive(Debug)]
pub enum Error {
    InvalidBase { base: String },
    NoNonRepeat {
        chrom: String,
        search_len: usize,
        search_start: usize,
    },
    NotIndel {
        chrom: String,
        pos: usize,
        ref_base: String,
        alt_base: String,
    },
    InvalidPosition {
        chrom: String,
        pos: usize,
    },
    EmptyVariant {
        chrom: String,
        pos: usize,
        ref_base: String,
        alt_base: String,
    },
    TwoBitError(twobit::Error),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Error::InvalidBase { base } => write!(f, "Wrong base: {}", base),
            Error::NoNonRepeat { chrom, search_len, search_start } => 
                write!(
                    f,
                    "Could not find a non-repeat sequence after {} bases from {}:{}",
                    search_len, chrom, search_start
                ),
            Error::NotIndel { chrom, pos, ref_base, alt_base } => write!(f, "Not an indel: {}:{}:{}:{}", chrom, pos, ref_base, alt_base),
            Error::InvalidPosition { chrom, pos } => write!(f, "Invalid genomic positions: {}:{}", chrom, pos),
            Error::EmptyVariant { chrom, pos, ref_base, alt_base } => write!(f, "Empty variant: {}:{}:{}:{}", chrom, pos, ref_base, alt_base),
            Error::TwoBitError(e) => write!(f, "TwoBitError: {}", e),
        }
    }
}

impl std::error::Error for Error {}


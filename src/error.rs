// Copyright Ryangguk Kim @ Oak Bioinformatics, LLC
// 
// This software is available under a dual licensing model, offering users the choice between the Affero General Public License version 3 (AGPL-3) for open-source use and a commercial license for proprietary or commercial use. 
// 
// To obtain a commercial license, please contact info@oakbioinformatics.com.

#[derive(Debug)]
pub struct InvalidBase {
    pub base: String,
}

impl core::fmt::Display for InvalidBase {
    fn fmt(&self, f: &mut core::fmt::Formatter) -> std::fmt::Result {
        write!(f, "Wrong base: {}", self.base)
    }
}

impl std::error::Error for InvalidBase {}

#[derive(Debug)]
pub struct NoNonRepeat {
    pub search_len: usize,
    pub search_start: usize,
}

impl core::fmt::Display for NoNonRepeat {
    fn fmt(&self, f: &mut core::fmt::Formatter) -> std::fmt::Result {
        write!(f, "Could not find a non-repeat sequence after {} bases from {}", self.search_len, self.search_start)
    }
}

impl std::error::Error for NoNonRepeat {}

#[derive(Debug)]
pub struct NotIndel {}

impl core::fmt::Display for NotIndel {
    fn fmt(&self, f: &mut core::fmt::Formatter) -> std::fmt::Result {
        write!(f, "Could not find a non-repeat sequence")
    }
}

impl std::error::Error for NotIndel {}

#[derive(Debug)]
pub struct InvalidPosition {
    pub chrom: String,
    pub pos: usize,
}

impl core::fmt::Display for InvalidPosition {
    fn fmt(&self, f: &mut core::fmt::Formatter) -> std::fmt::Result {
        write!(f, "Invalid genomic positions: {}:{}", self.chrom, self.pos)
    }
}

impl std::error::Error for InvalidPosition {}


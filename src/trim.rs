// Copyright Ryangguk Kim @ Oak Bioinformatics, LLC
// 
// This software is available under a dual licensing model, offering users the choice between the Affero General Public License version 3 (AGPL-3) for open-source use and a commercial license for proprietary or commercial use. 
// 
// To obtain a commercial license, please contact info@oakbioinformatics.com.

use noodles::vcf::record::reference_bases::Base;

pub fn trim_left(
    ref_bases: &Vec<Base>,
    alt_bases: &Vec<Base>,
    ref_start: usize,
    ref_end: usize,
    alt_start: usize,
    alt_end: usize,
) -> (usize, usize) {
    let ref_span = ref_end - ref_start;
    let alt_span = alt_end - alt_start;
    let span = std::cmp::min(ref_span, alt_span);
    let mut new_ref_start = ref_start;
    let mut new_alt_start = alt_start;
    for i in 0..span {
        let ref_pos = ref_start + i;
        let ref_base = ref_bases.get(ref_pos).unwrap();
        let alt_pos = alt_start + i;
        let alt_base = alt_bases.get(alt_pos).unwrap();
        if ref_base == alt_base {
            new_ref_start += 1;
            new_alt_start += 1;
        } else {
            break;
        }
    }
    (new_ref_start, new_alt_start)
}

pub fn trim_right(
    ref_bases: &Vec<Base>,
    alt_bases: &Vec<Base>,
    ref_start: usize,
    ref_end: usize,
    alt_start: usize,
    alt_end: usize,
) -> (usize, usize) {
    let ref_span = ref_end - ref_start;
    let alt_span = alt_end - alt_start;
    let span = std::cmp::min(ref_span, alt_span);
    let mut new_ref_end = ref_end;
    let mut new_alt_end = alt_end;
    for i in 1..(span + 1) {
        let ref_pos = ref_end - i;
        let ref_base = ref_bases.get(ref_pos).unwrap();
        let alt_pos = alt_end - i;
        let alt_base = alt_bases.get(alt_pos).unwrap();
        if ref_base == alt_base {
            new_ref_end -= 1;
            new_alt_end -= 1;
        } else {
            break;
        }
    }
    (new_ref_end, new_alt_end)
}


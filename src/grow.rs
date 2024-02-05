// Copyright Ryangguk Kim @ Oak Bioinformatics, LLC
// 
// This software is available under a dual licensing model, offering users the choice between the Affero General Public License version 3 (AGPL-3) for open-source use and a commercial license for proprietary or commercial use. 
// 
// To obtain a commercial license, please contact info@oakbioinformatics.com.

use anyhow::Result;
use crate::error::{InvalidPosition, NoNonRepeat, NotIndel};
use crate::util::is_base_same_as_char;

static GROW_RIGHT_MAX: usize = 30;

pub struct Grower {
    pub tb: twobit::TwoBitFile<std::io::BufReader<std::fs::File>>,
}

impl Grower {
    pub fn new(twobit_path: &str) -> Result<Grower> {
        let tb = twobit::TwoBitFile::open(twobit_path)?;
        Ok(Grower { tb })
    }

    pub fn grow_right(
        &mut self,
        chrom: &str,
        pos: usize,
        bases: &Vec<Base>,
    ) -> Result<(usize, Vec<Base>)> {
        let bases_len = bases.len();
        let mut expansion: Vec<Base> = Vec::with_capacity(bases_len);
        match bases_len {
            0 => return Ok((pos, expansion)),
            _ => {
                // 0-based
                let mut probe_start = pos;
                let mut probe_end = pos + bases_len;
                let mut growth_end = probe_start;
                let mut grow_c: usize = 0;
                loop {
                    let frag = self
                        .tb
                        .read_sequence(chrom, probe_start - 1..probe_end - 1)?;
                    if frag.len() == 0 {
                        return Err(anyhow::Error::new(InvalidPosition {
                            chrom: chrom.to_string(),
                            pos
                        }));
                    }
                    let mut chars = frag.chars();
                    let mut bases_iter = bases.iter();
                    let mut diff_found = false;
                    for offset in 0..bases_len {
                        let base = bases_iter.next().unwrap();
                        let c = chars.next().unwrap();
                        if !is_base_same_as_char(base, c) {
                            growth_end = probe_start + offset;
                            diff_found = true;
                            break;
                        } else {
                            expansion.push(base.clone());
                        }
                    }
                    if diff_found {
                        break;
                    }
                    probe_start = probe_start + bases_len;
                    probe_end = probe_start + bases_len;
                    grow_c += 1;
                    if grow_c > GROW_RIGHT_MAX {
                        return Err(anyhow::Error::new(NoNonRepeat { search_len: GROW_RIGHT_MAX, search_start: probe_start }));
                    }
                }
                Ok((growth_end, expansion))
            }
        }
    }

    pub fn grow_left(
        &mut self,
        chrom: &str,
        pos: usize,
        bases: &Vec<Base>,
    ) -> Result<(usize, Vec<Base>)> {
        let bases_len = bases.len();
        let mut probe_end = pos;
        let mut probe_start = probe_end - bases_len;
        let mut growth_start = probe_end;
        let mut expansion: Vec<Base> = Vec::with_capacity(bases_len);
        loop {
            let frag = self
                .tb
                .read_sequence(chrom, probe_start - 1..probe_end - 1)?;
            let mut chars = frag.chars().rev();
            let mut bases_iter = bases.iter().rev();
            let mut diff_found = false;
            for offset in 0..bases_len {
                let base = bases_iter.next().unwrap();
                let c = chars.next().unwrap();
                if !is_base_same_as_char(base, c) {
                    growth_start = probe_end - offset;
                    diff_found = true;
                    break;
                } else {
                    expansion.push(base.clone());
                }
            }
            if diff_found {
                break;
            }
            probe_end = probe_start;
            probe_start = probe_start - bases_len;
        }
        expansion.reverse();
        Ok((growth_start, expansion))
    }

    pub fn grow_pos_ref_alt(
        &mut self,
        chrom: &str,
        pos: usize,
        ref_bases: &Vec<Base>,
        alt_bases: &Vec<Base>,
    ) -> Result<(usize, Vec<Base>, Vec<Base>)> {
        let growth_left_start: usize;
        let growth_right_end: usize;
        let growth_left_bases: Vec<Base>;
        let growth_right_bases: Vec<Base>;
        let growth_right_size: usize;
        let ref_bases_len = ref_bases.len();
        let alt_bases_len = alt_bases.len();
        match ref_bases_len {
            0 => match alt_bases_len {
                0 => {
                    return Err(anyhow::Error::new(NotIndel {}));
                }
                _ => {
                    (growth_right_end, growth_right_bases) =
                        self.grow_right(chrom, pos, alt_bases)?;
                    (growth_left_start, growth_left_bases) =
                        self.grow_left(chrom, pos, alt_bases)?;
                    growth_right_size = growth_right_end - pos;
                }
            },
            _ => match alt_bases_len {
                0 => {
                    let growth_right_start: usize = pos + ref_bases_len;
                    (growth_right_end, growth_right_bases) =
                        self.grow_right(chrom, growth_right_start, ref_bases)?;
                    (growth_left_start, growth_left_bases) =
                        self.grow_left(chrom, pos, ref_bases)?;
                    growth_right_size = growth_right_end - pos - ref_bases_len + 1;
                }
                _ => {
                    return Err(anyhow::Error::new(NotIndel {}));
                }
            },
        }
        let growth_left_size = pos - growth_left_start;
        let new_ref_bases_len = growth_left_size + ref_bases_len + growth_right_size;
        let new_alt_bases_len = growth_left_size + alt_bases_len + growth_right_size;
        let mut new_ref_bases: Vec<Base> = Vec::with_capacity(new_ref_bases_len);
        let mut new_alt_bases: Vec<Base> = Vec::with_capacity(new_alt_bases_len);
        new_ref_bases.extend(growth_left_bases.iter());
        new_ref_bases.extend(ref_bases.iter());
        new_ref_bases.extend(growth_right_bases.iter());
        new_alt_bases.extend(growth_left_bases.iter());
        new_alt_bases.extend(alt_bases.iter());
        new_alt_bases.extend(growth_right_bases.iter());
        Ok((growth_left_start, new_ref_bases, new_alt_bases))
    }
}

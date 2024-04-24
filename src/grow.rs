// Copyright Ryangguk Kim @ Oak Bioinformatics, LLC
//
// This software is available under a dual licensing model, offering users the choice between the Affero General Public License version 3 (AGPL-3) for open-source use and a commercial license for proprietary or commercial use.
//
// To obtain a commercial license, please contact info@oakbioinformatics.com.

use crate::error::Error;
use crate::util::is_base_same_as_char;
use crate::Base;
use std::path::PathBuf;

static GROW_RIGHT_MAX: usize = 10000;

pub struct Grower {
    pub tb: twobit::TwoBitMemoryFile
}

impl Grower {
    pub fn new(twobit_path: &PathBuf) -> Result<Grower, Error> {
        match twobit::TwoBitFile::open_and_read(twobit_path) {
            Ok(tb) => Ok(Grower { tb }),
            Err(e) => Err(Error::TwoBitError(e)),
        }
    }

    pub fn grow_right(
        &mut self,
        chrom: &[u8],
        pos: usize,
        bases: &[Base],
    ) -> Result<(usize, Box<[Base]>), Error> {
        let bases_len = bases.len();
        let mut expansion: Vec<Base> = Vec::with_capacity(bases_len);
        match bases_len {
            0 => Ok((pos, expansion.into_boxed_slice())),
            _ => {
                // 0-based
                let mut probe_start = pos;
                let mut probe_end = pos + bases_len;
                let mut growth_end = probe_start;
                let mut grow_c: usize = 0;
                let chrom_str = std::str::from_utf8(chrom).unwrap();
                loop {
                    let frag = self.tb.read_sequence(chrom_str, probe_start - 1..probe_end - 1).unwrap_or("".to_string());
                    if frag.len() == 0 {
                        return Ok((growth_end, expansion.into_boxed_slice()));
                    }
                    let check_len: usize = std::cmp::min(frag.len(), bases_len);
                    let mut chars = frag.chars();
                    let mut bases_iter = bases.iter();
                    let mut diff_found = false;
                    for offset in 0..check_len {
                        let base = bases_iter.next().unwrap();
                        match chars.next() {
                            Some(c) => {
                                if !is_base_same_as_char(base, c) {
                                    growth_end = probe_start + offset;
                                    diff_found = true;
                                    break;
                                }
                                expansion.push(base.clone());
                            }
                            None => {
                                return Err(Error::InvalidPosition {
                                    chrom: chrom_str.to_string(),
                                    pos: probe_start + offset,
                                });
                            }
                        };
                    }
                    if diff_found {
                        break;
                    }
                    probe_start = probe_start + bases_len;
                    probe_end = probe_start + bases_len;
                    grow_c += 1;
                    if grow_c > GROW_RIGHT_MAX {
                        return Err(Error::NoNonRepeat {
                            chrom: chrom_str.to_string(),
                            search_len: GROW_RIGHT_MAX,
                            search_start: probe_start,
                        });
                    }
                }
                Ok((growth_end, expansion.into_boxed_slice()))
            }
        }
    }

    pub fn grow_left(
        &mut self,
        chrom: &[u8],
        pos: usize,
        bases: &[Base],
    ) -> Result<(usize, Box<[Base]>), Error> {
        let bases_len = bases.len();
        let mut growth_start = pos;
        let mut expansion: Vec<Base> = Vec::with_capacity(bases_len);
        let mut probe_end = pos;
        if bases_len == 0 {
            return Ok((growth_start, expansion.into_boxed_slice()));
        }
        if probe_end <= 1 {
            return Ok((growth_start, expansion.into_boxed_slice()));
        }
        let mut probe_start: usize;
        if probe_end <= bases_len {
            probe_start = 1;
        } else {
            probe_start = probe_end - bases_len;
        }
        let chrom_str = std::str::from_utf8(chrom).unwrap();
        loop {
            match self.tb.read_sequence(chrom_str, probe_start - 1..probe_end - 1) {
                Ok(frag) => {
                    if frag.len() == 0 {
                        return Ok((growth_start, expansion.into_boxed_slice()));
                    }
                    let mut chars = frag.chars().rev();
                    let mut bases_iter = bases.iter().rev();
                    let mut diff_found = false;
                    let check_len: usize = std::cmp::min(frag.len(), bases_len);
                    for offset in 0..check_len {
                        let base = bases_iter.next().unwrap();
                        match chars.next() {
                            Some(c) => {
                                if !is_base_same_as_char(base, c) {
                                    growth_start = probe_end - offset;
                                    diff_found = true;
                                    break;
                                }
                                expansion.push(base.clone());
                            }
                            None => {
                                return Err(Error::InvalidPosition {
                                    chrom: chrom_str.to_string(),
                                    pos: probe_end - offset,
                                });
                            }
                        }
                    }
                    if diff_found {
                        break;
                    }
                    if probe_start <= 1 || probe_start <= bases_len {
                        return Ok((growth_start, expansion.into_boxed_slice()));
                    }
                    probe_end = probe_start;
                    if probe_end <= 1 {
                        break;
                    }
                    if probe_start <= bases_len {
                        probe_start = 1;
                    } else {
                        probe_start = probe_start - bases_len;
                    }
                },
                Err(_) => {
                    return Err(Error::InvalidPosition {
                        chrom: chrom_str.to_string(),
                        pos: probe_start,
                    });
                }
            }
        }
        expansion.reverse();
        Ok((growth_start, expansion.into_boxed_slice()))
    }

    pub fn grow(
        &mut self,
        chrom: &[u8],
        pos: usize,
        ref_bases: &[Base],
        alt_bases: &[Base],
    ) -> Result<(usize, Box<[Base]>, Box<[Base]>), Error> {
        let growth_left_start: usize;
        let growth_right_end: usize;
        let growth_left_bases: Box<[Base]>;
        let growth_right_bases: Box<[Base]>;
        let growth_right_size: usize;
        let ref_bases_len = ref_bases.len();
        let alt_bases_len = alt_bases.len();
        match ref_bases_len {
            0 => match alt_bases_len {
                0 => {
                    return Err(Error::EmptyVariant { chrom: String::from_utf8(chrom.to_vec()).unwrap(), pos, ref_base: format!("{:?}", ref_bases), alt_base: format!("{:?}", alt_bases) });
                }
                _ => {
                    match self.grow_right(chrom, pos, alt_bases) {
                        Ok((end, bases)) => {
                            growth_right_end = end;
                            growth_right_bases = bases;
                        }
                        Err(e) => {
                            eprintln!(
                                "Error: {}. {}:{}:{:?}:{:?}",
                                e, String::from_utf8(chrom.to_vec()).unwrap(), pos, ref_bases, alt_bases
                            );
                            return Err(e);
                        }
                    }
                    match self.grow_left(chrom, pos, alt_bases) {
                        Ok((start, bases)) => {
                            growth_left_start = start;
                            growth_left_bases = bases;
                        }
                        Err(e) => {
                            eprintln!(
                                "Error: {}. {}:{}:{:?}:{:?}",
                                e, String::from_utf8(chrom.to_vec()).unwrap(), pos, ref_bases, alt_bases
                            );
                            return Err(e);
                        }
                    }
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
                    return Err(Error::NotIndel { chrom: String::from_utf8(chrom.to_vec()).unwrap(), pos, ref_base: format!("{:?}", ref_bases), alt_base: format!("{:?}", alt_bases) });
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
        Ok((growth_left_start, new_ref_bases.into_boxed_slice(), new_alt_bases.into_boxed_slice()))
    }
}

#[cfg(test)]
mod tests_grow {
    use super::*;

    #[test]
    fn test_grow_left() {
        let twobit_fname = std::env::var("TWOBIT_FNAME").unwrap();
        let mut grower = Grower::new(&std::path::PathBuf::from(twobit_fname)).unwrap();
        let result = grower.grow_left("chr19_GL383575v2_alt".as_bytes(), 1, &vec![Base::C, Base::A, Base::C, Base::A]).unwrap();
        assert_eq!(result, (1, vec![].into_boxed_slice()));
        let result = grower.grow_left("chr19_GL383575v2_alt".as_bytes(), 3, &vec![Base::A, Base::G, Base::C, Base::C]).unwrap();
        assert_eq!(result, (3, vec![Base::C, Base::C].into_boxed_slice()));
    }

    #[test]
    fn test_grow_right() {
        let twobit_fname = std::env::var("TWOBIT_FNAME").unwrap();
        let mut grower = Grower::new(&std::path::PathBuf::from(twobit_fname)).unwrap();
        let result = grower.grow_right("chr19_GL383576v1_alt".as_bytes(), 188023, &vec![Base::C, Base::A, Base::C, Base::A]).unwrap();
        assert_eq!(result, (188023, vec![].into_boxed_slice()));
        let result = grower.grow_right("chr19_GL383576v1_alt".as_bytes(), 188023, &vec![Base::T]).unwrap();
        assert_eq!(result, (188023, vec![Base::T, Base::T].into_boxed_slice()));
        let result = grower.grow_right("chr19_GL383576v1_alt".as_bytes(), 188023, &vec![Base::T, Base::T, Base::T]).unwrap();
        assert_eq!(result, (188023, vec![Base::T, Base::T].into_boxed_slice()));
    }
}

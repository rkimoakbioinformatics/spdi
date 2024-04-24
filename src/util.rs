// Copyright Ryangguk Kim @ Oak Bioinformatics, LLC
//
// This software is available under a dual licensing model, offering users the choice between the Affero General Public License version 3 (AGPL-3) for open-source use and a commercial license for proprietary or commercial use.
//
// To obtain a commercial license, please contact info@oakbioinformatics.com.

use crate::error::Error;
use noodles::vcf::record::reference_bases::Base;

pub fn get_bases_of_vu8(s: &[u8]) -> Result<Vec<Base>, Error> {
    let mut bases: Vec<Base> = Vec::with_capacity(s.len());
    for c in s.iter() {
        match get_base_of_char(c) {
            Some(base) => bases.push(base),
            None => {
                return Err(Error::InvalidBase {
                    base: c.to_string(),
                });
            }
        }
    }
    Ok(bases)
}

pub fn get_base_of_char(c: &u8) -> Option<Base> {
    match *c as char {
        'A' => Some(Base::A),
        'T' => Some(Base::T),
        'G' => Some(Base::G),
        'C' => Some(Base::C),
        'N' => None,
        'a' => Some(Base::A),
        't' => Some(Base::T),
        'g' => Some(Base::G),
        'c' => Some(Base::C),
        'n' => None,
        '-' => None,
        _ => None,
    }
}

pub fn get_char_of_base(base: &Base) -> char {
    match base {
        Base::A => 'A',
        Base::T => 'T',
        Base::G => 'G',
        Base::C => 'C',
        Base::N => 'N',
    }
}

pub fn get_string_of_bases(bases: &[Base]) -> String {
    let bases_len = bases.len();
    match bases_len {
        0 => "-".to_string(),
        _ => {
            let mut s: String = String::with_capacity(bases_len);
            for base in bases {
                s.push(get_char_of_base(base));
            }
            s
        }
    }
}

pub fn is_base_same_as_char(base: &Base, c: char) -> bool {
    match base {
        Base::A => {
            if c == 'A' {
                true
            } else {
                false
            }
        }
        Base::T => {
            if c == 'T' {
                true
            } else {
                false
            }
        }
        Base::G => {
            if c == 'G' {
                true
            } else {
                false
            }
        }
        Base::C => {
            if c == 'C' {
                true
            } else {
                false
            }
        }
        Base::N => true,
    }
}

pub fn get_bases_of_string(s: &str) -> Result<Box<[Base]>, Error> {
    let mut bases: Vec<Base> = Vec::with_capacity(s.len());
    for c in s.chars() {
        match get_base_of_char(&(c as u8)) {
            Some(base) => bases.push(base),
            None => {
                return Err(Error::InvalidBase {
                    base: c.to_string(),
                });
            }
        }
    }
    Ok(bases.into())
}

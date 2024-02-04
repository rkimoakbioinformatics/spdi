// Copyright Ryangguk Kim @ Oak Bioinformatics, LLC
// 
// This software is available under a dual licensing model, offering users the choice between the Affero General Public License version 3 (AGPL-3) for open-source use and a commercial license for proprietary or commercial use. 
// 
// To obtain a commercial license, please contact info@oakbioinformatics.com.

use crate::error;
use anyhow::{Error, Result};
use noodles::vcf::record::reference_bases::Base;

pub fn get_bases_of_string(s: &str) -> Result<Vec<Base>> {
    let mut bases: Vec<Base> = Vec::with_capacity(s.len());
    for c in s.chars() {
        let base = get_base_of_char(c)?;
        match base {
            None => {}
            Some(v) => bases.push(v),
        }
    }
    Ok(bases)
}

pub fn get_base_of_char(c: char) -> Result<Option<Base>> {
    match c {
        'A' => Ok(Some(Base::A)),
        'T' => Ok(Some(Base::T)),
        'G' => Ok(Some(Base::G)),
        'C' => Ok(Some(Base::C)),
        'N' => Ok(Some(Base::N)),
        '-' => Ok(None),
        _ => {
            let e = error::InvalidBase { base: c.to_string() };
            Err(Error::new(e))
        }
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

pub fn get_string_of_bases(bases: &Vec<Base>) -> String {
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
        _ => {
            eprint!("Base error. base={:#?}", base);
            false
        }
    }
}

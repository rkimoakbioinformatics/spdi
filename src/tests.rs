// Copyright Ryangguk Kim @ Oak Bioinformatics, LLC
//
// This software is available under a dual licensing model, offering users the choice between the Affero General Public License version 3 (AGPL-3) for open-source use and a commercial license for proprietary or commercial use.
//
// To obtain a commercial license, please contact info@oakbioinformatics.com.

#[test]
fn get_spdi_string_del() {
    use crate::SPDI;
    use std::path::PathBuf;

    let twobit_fname = std::env::var("TWOBIT_FNAME").unwrap().to_string();
    let mut spdi = SPDI::new(&PathBuf::from(twobit_fname)).unwrap();
    assert_eq!(
        "chr1:141454:TTATTATTATTATT:TTATTATTATT".to_string(),
        spdi.get_spdi_string("chr1".as_bytes(), 141457, "TTA".as_bytes(), "".as_bytes()).unwrap()
    );
    assert_eq!(
        "chr1:141454:TTATTATTATTATT:TTATTATTATT".to_string(),
        spdi.get_spdi_string("chr1".as_bytes(), 141455, "TAT".as_bytes(), "".as_bytes()).unwrap()
    );
}

#[test]
fn get_spdi_string_ins() {
    use crate::SPDI;
    use std::path::PathBuf;

    let twobit_fname = std::env::var("TWOBIT_FNAME").unwrap().to_string();
    let mut spdi = SPDI::new(&PathBuf::from(twobit_fname)).unwrap();
    assert_eq!(
        "chr1:141454:TTATTATTATTATT:TTATTATTATTATTATT".to_string(),
        spdi.get_spdi_string("chr1".as_bytes(), 141458, "".as_bytes(), "TAT".as_bytes()).unwrap()
    );
    assert_eq!(
        "chr1:141454:TTATTATTATTATT:TTATTATTATTATTATT".to_string(),
        spdi.get_spdi_string("chr1".as_bytes(), 141457, "".as_bytes(), "TTA".as_bytes()).unwrap()
    );
    assert_eq!(
        "chr1:4950532:AAAATAAAATAAAATAAAATAAAATAA:AAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAA".to_string(),
        spdi.get_spdi_string("chr1".as_bytes(), 4950531, "AAAAATAAAAT".as_bytes(), "AAAAATAAAATAAAATAAAATAAAATAAAATAAAAT".as_bytes()).unwrap()
    );
}

#[test]
fn get_spdi_string_snv() {
    use crate::SPDI;
    use std::path::PathBuf;

    let twobit_fname = std::env::var("TWOBIT_FNAME").unwrap().to_string();
    let mut spdi = SPDI::new(&PathBuf::from(twobit_fname)).unwrap();
    assert_eq!(
        "chr1:141453:G:C".to_string(),
        spdi.get_spdi_string("chr1".as_bytes(), 141453, "G".as_bytes(), "C".as_bytes()).unwrap()
    );
    assert_eq!(
        "chr1:30025797:A:T".to_string(),
        spdi.get_spdi_string("chr1".as_bytes(), 30025797, "A".as_bytes(), "T".as_bytes()).unwrap()
    );
}

#[test]
fn trim_right() {
    use crate::trim::trim_right;
    use crate::util::get_bases_of_string;
    use noodles::vcf::record::reference_bases::Base;

    // ACTGTC
    // AC  TC
    let ref_bases: Box<[Base]> = get_bases_of_string("ACTGTC").unwrap();
    let alt_bases: Box<[Base]> = get_bases_of_string("ACTC").unwrap();
    assert_eq!((4, 2), trim_right(&ref_bases, &alt_bases, 0, 6, 0, 4));
    // A
    // G
    let ref_bases: Box<[Base]> = get_bases_of_string("A").unwrap();
    let alt_bases: Box<[Base]> = get_bases_of_string("G").unwrap();
    assert_eq!((1, 1), trim_right(&ref_bases, &alt_bases, 0, 1, 0, 1));
}

#[test]
fn trim_left() {
    use crate::trim::trim_left;
    use crate::util::get_bases_of_string;
    use noodles::vcf::record::reference_bases::Base;

    // A
    // G
    let ref_bases: Box<[Base]> = get_bases_of_string("A").unwrap();
    let alt_bases: Box<[Base]> = get_bases_of_string("G").unwrap();
    assert_eq!((0, 0), trim_left(&ref_bases, &alt_bases, 0, 1, 0, 1));
    // ACTG
    // AC
    let ref_bases: Box<[Base]> = get_bases_of_string("ACTG").unwrap();
    let alt_bases: Box<[Base]> = get_bases_of_string("AC").unwrap();
    assert_eq!((2, 2), trim_left(&ref_bases, &alt_bases, 0, 4, 0, 2));
}

#[test]
fn grow_left() {
    use crate::util::get_bases_of_string;
    use crate::SPDI;
    use noodles::vcf::record::reference_bases::Base;
    use std::path::PathBuf;

    // TTA(TTA)TTA
    //     ---
    let twobit_fname = std::env::var("TWOBIT_FNAME").unwrap().to_string();
    let mut spdi = SPDI::new(&PathBuf::from(twobit_fname)).unwrap();
    let ref_bases: Box<[Base]> = get_bases_of_string("TTA").unwrap();
    assert_eq!(
        (141454, get_bases_of_string("TTA").unwrap()),
        spdi.grower.grow_left("chr1".as_bytes(), 141457, &ref_bases).unwrap()
    );
}

#[test]
fn grow_right() {
    use crate::util::get_bases_of_string;
    use crate::SPDI;
    use noodles::vcf::record::reference_bases::Base;
    use std::path::PathBuf;

    // TTA(TTA)TTA
    //     ---
    let twobit_fname = std::env::var("TWOBIT_FNAME").unwrap().to_string();
    let mut spdi = SPDI::new(&PathBuf::from(twobit_fname)).unwrap();
    let bases: Box<[Base]> = get_bases_of_string("TTA").unwrap();
    assert_eq!(
        (141468, get_bases_of_string("TTATTATTATT").unwrap()),
        spdi.grower.grow_right("chr1".as_bytes(), 141457, &bases).unwrap()
    );
}

#[test]
fn grow() {
    use crate::util::get_bases_of_string;
    use crate::SPDI;
    use noodles::vcf::record::reference_bases::Base;
    use std::path::PathBuf;

    let twobit_fname = std::env::var("TWOBIT_FNAME").unwrap().to_string();
    let mut spdi = SPDI::new(&PathBuf::from(twobit_fname)).unwrap();
    // TTA(TTA)TTA
    //     ---
    let ref_bases: Box<[Base]> = get_bases_of_string("TTA").unwrap();
    let alt_bases: Box<[Base]> = get_bases_of_string("").unwrap();
    assert_eq!(
        (
            141454,
            get_bases_of_string("TTATTATTATTATT").unwrap(),
            get_bases_of_string("TTATTATTATT").unwrap()
        ),
        spdi.grower
            .grow("chr1".as_bytes(), 141457, &ref_bases, &alt_bases)
            .unwrap()
    );
    let ref_bases: Box<[Base]> = get_bases_of_string("").unwrap();
    let alt_bases: Box<[Base]> = get_bases_of_string("TTA").unwrap();
    assert_eq!(
        (
            141454,
            get_bases_of_string("TTATTATTATTATT").unwrap(),
            get_bases_of_string("TTATTATTATTATTATT").unwrap()
        ),
        spdi.grower
            .grow("chr1".as_bytes(), 141457, &ref_bases, &alt_bases)
            .unwrap()
    );
}

#[test]
fn long_conversion() {
    use crate::util::get_bases_of_string;
    use crate::SPDI;
    use crate::util::get_char_of_base;
    use noodles::vcf::record::reference_bases::Base;
    use std::path::PathBuf;


    let twobit_fname = std::env::var("TWOBIT_FNAME").unwrap().to_string();
    let mut spdi = SPDI::new(&PathBuf::from(&twobit_fname)).unwrap();
    let ref_bases: Box<[Base]> = get_bases_of_string("GATTC").unwrap();
    let alt_bases: Box<[Base]> = get_bases_of_string("GATTCTATTC").unwrap();
    let result = spdi.get_spdi_conversion("chr22".as_bytes(), 45795354, &ref_bases, &alt_bases)
            .unwrap();
    println!("converted ref_base_len={} alt_base_len={} total={}", result.1.len(), result.2.len(), result.1.len() + result.2.len());
    println!("orig ref_base_len={} alt_base_len={} total={}", ref_bases.len(), alt_bases.len(), ref_bases.len() + alt_bases.len());
    println!("converted ref_bases={}", result.1.iter().map(|b| get_char_of_base(b)).collect::<String>());
    println!("converted alt_bases={}", result.2.iter().map(|b| get_char_of_base(b)).collect::<String>());
    assert_eq!(
        (
            45795355,
            get_bases_of_string("ATTCTATTCTATTCTATTCTATTCTATTCTATTCTATTCTATTCTATTCTATTCTATTCTATTCTATTCT").unwrap(),
            get_bases_of_string("ATTCTATTCTATTCTATTCTATTCTATTCTATTCTATTCTATTCTATTCTATTCTATTCTATTCTATTCTATTCT").unwrap(),
        ),
        result
    );
}

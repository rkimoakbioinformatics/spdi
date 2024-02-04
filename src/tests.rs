// Copyright Ryangguk Kim @ Oak Bioinformatics, LLC
// 
// This software is available under a dual licensing model, offering users the choice between the Affero General Public License version 3 (AGPL-3) for open-source use and a commercial license for proprietary or commercial use. 
// 
// To obtain a commercial license, please contact info@oakbioinformatics.com.

#[test]
fn get_spdi_string_del() {
    use crate::SPDI;

    let mut spdi = SPDI::new("/Users/rick/Downloads/hg38.2bit").unwrap();
    assert_eq!(
        "chr1:141454:TTATTATTATTATT:TTATTATTATT".to_string(),
        spdi.get_spdi_string("chr1", 141457, "TTA", "-").unwrap()
    );
    assert_eq!(
        "chr1:141454:TTATTATTATTATT:TTATTATTATT".to_string(),
        spdi.get_spdi_string("chr1", 141455, "TAT", "-").unwrap()
    );
}

#[test]
fn get_spdi_string_ins() {
    use crate::SPDI;

    let mut spdi = SPDI::new("/Users/rick/Downloads/hg38.2bit").unwrap();
    assert_eq!(
        "chr1:141454:TTATTATTATTATT:TTATTATTATTATTATT".to_string(),
        spdi.get_spdi_string("chr1", 141458, "-", "TAT").unwrap()
    );
    assert_eq!(
        "chr1:141454:TTATTATTATTATT:TTATTATTATTATTATT".to_string(),
        spdi.get_spdi_string("chr1", 141457, "-", "TTA").unwrap()
    );
}

#[test]
fn get_spdi_string_snv() {
    use crate::SPDI;

    let mut spdi = SPDI::new("/Users/rick/Downloads/hg38.2bit").unwrap();
    assert_eq!(
        "chr1:141453:G:C".to_string(),
        spdi.get_spdi_string("chr1", 141453, "G", "C").unwrap()
    );
    assert_eq!(
        "chr1:30025797:A:T".to_string(),
        spdi.get_spdi_string("chr1", 30025797, "A", "T").unwrap()
    );
}

#[test]
fn trim_right() {
    use noodles::vcf::record::reference_bases::Base;
    use crate::util::get_bases_of_string;
    use crate::trim::trim_right;

    // ACTGTC
    // AC  TC
    let ref_bases: Vec<Base> = get_bases_of_string("ACTGTC").unwrap();
    let alt_bases: Vec<Base> = get_bases_of_string("ACTC").unwrap();
    assert_eq!((4, 2), trim_right(&ref_bases, &alt_bases, 0, 6, 0, 4));
    // A
    // G
    let ref_bases: Vec<Base> = get_bases_of_string("A").unwrap();
    let alt_bases: Vec<Base> = get_bases_of_string("G").unwrap();
    assert_eq!((1, 1), trim_right(&ref_bases, &alt_bases, 0, 1, 0, 1));
}

#[test]
fn trim_left() {
    use noodles::vcf::record::reference_bases::Base;
    use crate::util::get_bases_of_string;
    use crate::trim::trim_left;

    // A
    // G
    let ref_bases: Vec<Base> = get_bases_of_string("A").unwrap();
    let alt_bases: Vec<Base> = get_bases_of_string("G").unwrap();
    assert_eq!((0, 0), trim_left(&ref_bases, &alt_bases, 0, 1, 0, 1));
    // ACTG
    // AC
    let ref_bases: Vec<Base> = get_bases_of_string("ACTG").unwrap();
    let alt_bases: Vec<Base> = get_bases_of_string("AC").unwrap();
    assert_eq!((2, 2), trim_left(&ref_bases, &alt_bases, 0, 4, 0, 2));
}

#[test]
fn grow_left() {
    use noodles::vcf::record::reference_bases::Base;
    use crate::SPDI;
    use crate::util::get_bases_of_string;

    // TTA(TTA)TTA
    //     ---
    let mut spdi = SPDI::new("/Users/rick/Downloads/hg38.2bit").unwrap();
    let ref_bases: Vec<Base> = get_bases_of_string("TTA").unwrap();
    assert_eq!(
        (141454, get_bases_of_string("TTA").unwrap()),
        spdi.grower.grow_left("chr1", 141457, &ref_bases).unwrap()
    );
}

#[test]
fn grow_right() {
    use noodles::vcf::record::reference_bases::Base;
    use crate::SPDI;
    use crate::util::get_bases_of_string;

    // TTA(TTA)TTA
    //     ---
    let mut spdi = SPDI::new("/Users/rick/Downloads/hg38.2bit").unwrap();
    let bases: Vec<Base> = get_bases_of_string("TTA").unwrap();
    assert_eq!(
        (141468, get_bases_of_string("TTATTATTATT").unwrap()),
        spdi.grower.grow_right("chr1", 141457, &bases).unwrap()
    );
}

#[test]
fn grow_pos_ref_alt() {
    use noodles::vcf::record::reference_bases::Base;
    use crate::SPDI;
    use crate::util::get_bases_of_string;

    let mut spdi = SPDI::new("/Users/rick/Downloads/hg38.2bit").unwrap();
    // TTA(TTA)TTA
    //     ---
    let ref_bases: Vec<Base> = get_bases_of_string("TTA").unwrap();
    let alt_bases: Vec<Base> = get_bases_of_string("-").unwrap();
    assert_eq!(
        (
            141454,
            get_bases_of_string("TTATTATTATTATT").unwrap(),
            get_bases_of_string("TTATTATTATT").unwrap()
        ),
        spdi.grower.grow_pos_ref_alt("chr1", 141457, &ref_bases, &alt_bases)
            .unwrap()
    );
    let ref_bases: Vec<Base> = get_bases_of_string("-").unwrap();
    let alt_bases: Vec<Base> = get_bases_of_string("TTA").unwrap();
    assert_eq!(
        (
            141454,
            get_bases_of_string("TTATTATTATTATT").unwrap(),
            get_bases_of_string("TTATTATTATTATTATT").unwrap()
        ),
        spdi.grower.grow_pos_ref_alt("chr1", 141457, &ref_bases, &alt_bases)
            .unwrap()
    );
}


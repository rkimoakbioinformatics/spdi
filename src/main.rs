// Copyright Ryangguk Kim @ Oak Bioinformatics, LLC
// 
// This software is available under a dual licensing model, offering users the choice between the Affero General Public License version 3 (AGPL-3) for open-source use and a commercial license for proprietary or commercial use. 
// 
// To obtain a commercial license, please contact info@oakbioinformatics.com.

use clap::Parser;
use std::io::BufRead;
use lazy_static::lazy_static;
use ahash::HashMap;
mod error;
use spdi::SPDI;

lazy_static! {
    static ref CHROMS: HashMap<&'static str, &'static str> = {
        let mut m = HashMap::default();
        m.insert("KI270721.1", "chr11_KI270721v1_random");
        m.insert("GL000009.2", "chr14_GL000009v2_random");
        m.insert("GL000225.1", "chr14_GL000225v1_random");
        m.insert("KI270722.1", "chr14_KI270722v1_random");
        m.insert("GL000194.1", "chr14_GL000194v1_random");
        m.insert("KI270723.1", "chr14_KI270723v1_random");
        m.insert("KI270724.1", "chr14_KI270724v1_random");
        m.insert("KI270725.1", "chr14_KI270725v1_random");
        m.insert("KI270726.1", "chr14_KI270726v1_random");
        m.insert("KI270727.1", "chr15_KI270727v1_random");
        m.insert("KI270728.1", "chr16_KI270728v1_random");
        m.insert("GL000205.2", "chr17_GL000205v2_random");
        m.insert("KI270729.1", "chr17_KI270729v1_random");
        m.insert("KI270730.1", "chr17_KI270730v1_random");
        m.insert("KI270706.1", "chr1_KI270706v1_random");
        m.insert("KI270707.1", "chr1_KI270707v1_random");
        m.insert("KI270708.1", "chr1_KI270708v1_random");
        m.insert("KI270709.1", "chr1_KI270709v1_random");
        m.insert("KI270710.1", "chr1_KI270710v1_random");
        m.insert("KI270711.1", "chr1_KI270711v1_random");
        m.insert("KI270712.1", "chr1_KI270712v1_random");
        m.insert("KI270713.1", "chr1_KI270713v1_random");
        m.insert("KI270714.1", "chr1_KI270714v1_random");
        m.insert("KI270731.1", "chr22_KI270731v1_random");
        m.insert("KI270732.1", "chr22_KI270732v1_random");
        m.insert("KI270733.1", "chr22_KI270733v1_random");
        m.insert("KI270734.1", "chr22_KI270734v1_random");
        m.insert("KI270735.1", "chr22_KI270735v1_random");
        m.insert("KI270736.1", "chr22_KI270736v1_random");
        m.insert("KI270737.1", "chr22_KI270737v1_random");
        m.insert("KI270738.1", "chr22_KI270738v1_random");
        m.insert("KI270739.1", "chr22_KI270739v1_random");
        m.insert("KI270715.1", "chr2_KI270715v1_random");
        m.insert("KI270716.1", "chr2_KI270716v1_random");
        m.insert("GL000221.1", "chr3_GL000221v1_random");
        m.insert("GL000008.2", "chr4_GL000008v2_random");
        m.insert("GL000208.1", "chr5_GL000208v1_random");
        m.insert("KI270717.1", "chr9_KI270717v1_random");
        m.insert("KI270718.1", "chr9_KI270718v1_random");
        m.insert("KI270719.1", "chr9_KI270719v1_random");
        m.insert("KI270720.1", "chr9_KI270720v1_random");
        m.insert("KI270762.1", "chr1_KI270762v1_alt");
        m.insert("KI270766.1", "chr1_KI270766v1_alt");
        m.insert("KI270760.1", "chr1_KI270760v1_alt");
        m.insert("KI270765.1", "chr1_KI270765v1_alt");
        m.insert("GL383518.1", "chr1_GL383518v1_alt");
        m.insert("GL383519.1", "chr1_GL383519v1_alt");
        m.insert("GL383520.2", "chr1_GL383520v2_alt");
        m.insert("KI270764.1", "chr1_KI270764v1_alt");
        m.insert("KI270763.1", "chr1_KI270763v1_alt");
        m.insert("KI270759.1", "chr1_KI270759v1_alt");
        m.insert("KI270761.1", "chr1_KI270761v1_alt");
        m.insert("KI270770.1", "chr2_KI270770v1_alt");
        m.insert("KI270773.1", "chr2_KI270773v1_alt");
        m.insert("KI270774.1", "chr2_KI270774v1_alt");
        m.insert("KI270769.1", "chr2_KI270769v1_alt");
        m.insert("GL383521.1", "chr2_GL383521v1_alt");
        m.insert("KI270772.1", "chr2_KI270772v1_alt");
        m.insert("KI270775.1", "chr2_KI270775v1_alt");
        m.insert("KI270771.1", "chr2_KI270771v1_alt");
        m.insert("KI270768.1", "chr2_KI270768v1_alt");
        m.insert("GL582966.2", "chr2_GL582966v2_alt");
        m.insert("GL383522.1", "chr2_GL383522v1_alt");
        m.insert("KI270776.1", "chr2_KI270776v1_alt");
        m.insert("KI270767.1", "chr2_KI270767v1_alt");
        m.insert("JH636055.2", "chr3_JH636055v2_alt");
        m.insert("KI270783.1", "chr3_KI270783v1_alt");
        m.insert("KI270780.1", "chr3_KI270780v1_alt");
        m.insert("GL383526.1", "chr3_GL383526v1_alt");
        m.insert("KI270777.1", "chr3_KI270777v1_alt");
        m.insert("KI270778.1", "chr3_KI270778v1_alt");
        m.insert("KI270781.1", "chr3_KI270781v1_alt");
        m.insert("KI270779.1", "chr3_KI270779v1_alt");
        m.insert("KI270782.1", "chr3_KI270782v1_alt");
        m.insert("KI270784.1", "chr3_KI270784v1_alt");
        m.insert("KI270790.1", "chr4_KI270790v1_alt");
        m.insert("GL383528.1", "chr4_GL383528v1_alt");
        m.insert("KI270787.1", "chr4_KI270787v1_alt");
        m.insert("GL000257.2", "chr4_GL000257v2_alt");
        m.insert("KI270788.1", "chr4_KI270788v1_alt");
        m.insert("GL383527.1", "chr4_GL383527v1_alt");
        m.insert("KI270785.1", "chr4_KI270785v1_alt");
        m.insert("KI270789.1", "chr4_KI270789v1_alt");
        m.insert("KI270786.1", "chr4_KI270786v1_alt");
        m.insert("KI270793.1", "chr5_KI270793v1_alt");
        m.insert("KI270792.1", "chr5_KI270792v1_alt");
        m.insert("KI270791.1", "chr5_KI270791v1_alt");
        m.insert("GL383532.1", "chr5_GL383532v1_alt");
        m.insert("GL949742.1", "chr5_GL949742v1_alt");
        m.insert("KI270794.1", "chr5_KI270794v1_alt");
        m.insert("GL339449.2", "chr5_GL339449v2_alt");
        m.insert("GL383530.1", "chr5_GL383530v1_alt");
        m.insert("KI270796.1", "chr5_KI270796v1_alt");
        m.insert("GL383531.1", "chr5_GL383531v1_alt");
        m.insert("KI270795.1", "chr5_KI270795v1_alt");
        m.insert("GL000250.2", "chr6_GL000250v2_alt");
        m.insert("KI270800.1", "chr6_KI270800v1_alt");
        m.insert("KI270799.1", "chr6_KI270799v1_alt");
        m.insert("GL383533.1", "chr6_GL383533v1_alt");
        m.insert("KI270801.1", "chr6_KI270801v1_alt");
        m.insert("KI270802.1", "chr6_KI270802v1_alt");
        m.insert("KB021644.2", "chr6_KB021644v2_alt");
        m.insert("KI270797.1", "chr6_KI270797v1_alt");
        m.insert("KI270798.1", "chr6_KI270798v1_alt");
        m.insert("KI270804.1", "chr7_KI270804v1_alt");
        m.insert("KI270809.1", "chr7_KI270809v1_alt");
        m.insert("KI270806.1", "chr7_KI270806v1_alt");
        m.insert("GL383534.2", "chr7_GL383534v2_alt");
        m.insert("KI270803.1", "chr7_KI270803v1_alt");
        m.insert("KI270808.1", "chr7_KI270808v1_alt");
        m.insert("KI270807.1", "chr7_KI270807v1_alt");
        m.insert("KI270805.1", "chr7_KI270805v1_alt");
        m.insert("KI270818.1", "chr8_KI270818v1_alt");
        m.insert("KI270812.1", "chr8_KI270812v1_alt");
        m.insert("KI270811.1", "chr8_KI270811v1_alt");
        m.insert("KI270821.1", "chr8_KI270821v1_alt");
        m.insert("KI270813.1", "chr8_KI270813v1_alt");
        m.insert("KI270822.1", "chr8_KI270822v1_alt");
        m.insert("KI270814.1", "chr8_KI270814v1_alt");
        m.insert("KI270810.1", "chr8_KI270810v1_alt");
        m.insert("KI270819.1", "chr8_KI270819v1_alt");
        m.insert("KI270820.1", "chr8_KI270820v1_alt");
        m.insert("KI270817.1", "chr8_KI270817v1_alt");
        m.insert("KI270816.1", "chr8_KI270816v1_alt");
        m.insert("KI270815.1", "chr8_KI270815v1_alt");
        m.insert("GL383539.1", "chr9_GL383539v1_alt");
        m.insert("GL383540.1", "chr9_GL383540v1_alt");
        m.insert("GL383541.1", "chr9_GL383541v1_alt");
        m.insert("GL383542.1", "chr9_GL383542v1_alt");
        m.insert("KI270823.1", "chr9_KI270823v1_alt");
        m.insert("GL383545.1", "chr10_GL383545v1_alt");
        m.insert("KI270824.1", "chr10_KI270824v1_alt");
        m.insert("GL383546.1", "chr10_GL383546v1_alt");
        m.insert("KI270825.1", "chr10_KI270825v1_alt");
        m.insert("KI270832.1", "chr11_KI270832v1_alt");
        m.insert("KI270830.1", "chr11_KI270830v1_alt");
        m.insert("KI270831.1", "chr11_KI270831v1_alt");
        m.insert("KI270829.1", "chr11_KI270829v1_alt");
        m.insert("GL383547.1", "chr11_GL383547v1_alt");
        m.insert("JH159136.1", "chr11_JH159136v1_alt");
        m.insert("JH159137.1", "chr11_JH159137v1_alt");
        m.insert("KI270827.1", "chr11_KI270827v1_alt");
        m.insert("KI270826.1", "chr11_KI270826v1_alt");
        m.insert("GL877875.1", "chr12_GL877875v1_alt");
        m.insert("GL877876.1", "chr12_GL877876v1_alt");
        m.insert("KI270837.1", "chr12_KI270837v1_alt");
        m.insert("GL383549.1", "chr12_GL383549v1_alt");
        m.insert("KI270835.1", "chr12_KI270835v1_alt");
        m.insert("GL383550.2", "chr12_GL383550v2_alt");
        m.insert("GL383552.1", "chr12_GL383552v1_alt");
        m.insert("GL383553.2", "chr12_GL383553v2_alt");
        m.insert("KI270834.1", "chr12_KI270834v1_alt");
        m.insert("GL383551.1", "chr12_GL383551v1_alt");
        m.insert("KI270833.1", "chr12_KI270833v1_alt");
        m.insert("KI270836.1", "chr12_KI270836v1_alt");
        m.insert("KI270840.1", "chr13_KI270840v1_alt");
        m.insert("KI270839.1", "chr13_KI270839v1_alt");
        m.insert("KI270843.1", "chr13_KI270843v1_alt");
        m.insert("KI270841.1", "chr13_KI270841v1_alt");
        m.insert("KI270838.1", "chr13_KI270838v1_alt");
        m.insert("KI270842.1", "chr13_KI270842v1_alt");
        m.insert("KI270844.1", "chr14_KI270844v1_alt");
        m.insert("KI270847.1", "chr14_KI270847v1_alt");
        m.insert("KI270845.1", "chr14_KI270845v1_alt");
        m.insert("KI270846.1", "chr14_KI270846v1_alt");
        m.insert("KI270852.1", "chr15_KI270852v1_alt");
        m.insert("KI270851.1", "chr15_KI270851v1_alt");
        m.insert("KI270848.1", "chr15_KI270848v1_alt");
        m.insert("GL383554.1", "chr15_GL383554v1_alt");
        m.insert("KI270849.1", "chr15_KI270849v1_alt");
        m.insert("GL383555.2", "chr15_GL383555v2_alt");
        m.insert("KI270850.1", "chr15_KI270850v1_alt");
        m.insert("KI270854.1", "chr16_KI270854v1_alt");
        m.insert("KI270856.1", "chr16_KI270856v1_alt");
        m.insert("KI270855.1", "chr16_KI270855v1_alt");
        m.insert("KI270853.1", "chr16_KI270853v1_alt");
        m.insert("GL383556.1", "chr16_GL383556v1_alt");
        m.insert("GL383557.1", "chr16_GL383557v1_alt");
        m.insert("GL383563.3", "chr17_GL383563v3_alt");
        m.insert("KI270862.1", "chr17_KI270862v1_alt");
        m.insert("KI270861.1", "chr17_KI270861v1_alt");
        m.insert("KI270857.1", "chr17_KI270857v1_alt");
        m.insert("JH159146.1", "chr17_JH159146v1_alt");
        m.insert("JH159147.1", "chr17_JH159147v1_alt");
        m.insert("GL383564.2", "chr17_GL383564v2_alt");
        m.insert("GL000258.2", "chr17_GL000258v2_alt");
        m.insert("GL383565.1", "chr17_GL383565v1_alt");
        m.insert("KI270858.1", "chr17_KI270858v1_alt");
        m.insert("KI270859.1", "chr17_KI270859v1_alt");
        m.insert("GL383566.1", "chr17_GL383566v1_alt");
        m.insert("KI270860.1", "chr17_KI270860v1_alt");
        m.insert("KI270864.1", "chr18_KI270864v1_alt");
        m.insert("GL383567.1", "chr18_GL383567v1_alt");
        m.insert("GL383570.1", "chr18_GL383570v1_alt");
        m.insert("GL383571.1", "chr18_GL383571v1_alt");
        m.insert("GL383568.1", "chr18_GL383568v1_alt");
        m.insert("GL383569.1", "chr18_GL383569v1_alt");
        m.insert("GL383572.1", "chr18_GL383572v1_alt");
        m.insert("KI270863.1", "chr18_KI270863v1_alt");
        m.insert("KI270868.1", "chr19_KI270868v1_alt");
        m.insert("KI270865.1", "chr19_KI270865v1_alt");
        m.insert("GL383573.1", "chr19_GL383573v1_alt");
        m.insert("GL383575.2", "chr19_GL383575v2_alt");
        m.insert("GL383576.1", "chr19_GL383576v1_alt");
        m.insert("GL383574.1", "chr19_GL383574v1_alt");
        m.insert("KI270866.1", "chr19_KI270866v1_alt");
        m.insert("KI270867.1", "chr19_KI270867v1_alt");
        m.insert("GL949746.1", "chr19_GL949746v1_alt");
        m.insert("GL383577.2", "chr20_GL383577v2_alt");
        m.insert("KI270869.1", "chr20_KI270869v1_alt");
        m.insert("KI270871.1", "chr20_KI270871v1_alt");
        m.insert("KI270870.1", "chr20_KI270870v1_alt");
        m.insert("GL383578.2", "chr21_GL383578v2_alt");
        m.insert("KI270874.1", "chr21_KI270874v1_alt");
        m.insert("KI270873.1", "chr21_KI270873v1_alt");
        m.insert("GL383579.2", "chr21_GL383579v2_alt");
        m.insert("GL383580.2", "chr21_GL383580v2_alt");
        m.insert("GL383581.2", "chr21_GL383581v2_alt");
        m.insert("KI270872.1", "chr21_KI270872v1_alt");
        m.insert("KI270875.1", "chr22_KI270875v1_alt");
        m.insert("KI270878.1", "chr22_KI270878v1_alt");
        m.insert("KI270879.1", "chr22_KI270879v1_alt");
        m.insert("KI270876.1", "chr22_KI270876v1_alt");
        m.insert("KI270877.1", "chr22_KI270877v1_alt");
        m.insert("GL383583.2", "chr22_GL383583v2_alt");
        m.insert("GL383582.2", "chr22_GL383582v2_alt");
        m.insert("KI270880.1", "chrX_KI270880v1_alt");
        m.insert("KI270881.1", "chrX_KI270881v1_alt");
        m.insert("KI270882.1", "chr19_KI270882v1_alt");
        m.insert("KI270883.1", "chr19_KI270883v1_alt");
        m.insert("KI270884.1", "chr19_KI270884v1_alt");
        m.insert("KI270885.1", "chr19_KI270885v1_alt");
        m.insert("KI270886.1", "chr19_KI270886v1_alt");
        m.insert("KI270887.1", "chr19_KI270887v1_alt");
        m.insert("KI270888.1", "chr19_KI270888v1_alt");
        m.insert("KI270889.1", "chr19_KI270889v1_alt");
        m.insert("KI270890.1", "chr19_KI270890v1_alt");
        m.insert("KI270891.1", "chr19_KI270891v1_alt");
        m.insert("KI270892.1", "chr1_KI270892v1_alt");
        m.insert("KI270894.1", "chr2_KI270894v1_alt");
        m.insert("KI270893.1", "chr2_KI270893v1_alt");
        m.insert("KI270895.1", "chr3_KI270895v1_alt");
        m.insert("KI270896.1", "chr4_KI270896v1_alt");
        m.insert("KI270897.1", "chr5_KI270897v1_alt");
        m.insert("KI270898.1", "chr5_KI270898v1_alt");
        m.insert("GL000251.2", "chr6_GL000251v2_alt");
        m.insert("KI270899.1", "chr7_KI270899v1_alt");
        m.insert("KI270901.1", "chr8_KI270901v1_alt");
        m.insert("KI270900.1", "chr8_KI270900v1_alt");
        m.insert("KI270902.1", "chr11_KI270902v1_alt");
        m.insert("KI270903.1", "chr11_KI270903v1_alt");
        m.insert("KI270904.1", "chr12_KI270904v1_alt");
        m.insert("KI270906.1", "chr15_KI270906v1_alt");
        m.insert("KI270905.1", "chr15_KI270905v1_alt");
        m.insert("KI270907.1", "chr17_KI270907v1_alt");
        m.insert("KI270910.1", "chr17_KI270910v1_alt");
        m.insert("KI270909.1", "chr17_KI270909v1_alt");
        m.insert("JH159148.1", "chr17_JH159148v1_alt");
        m.insert("KI270908.1", "chr17_KI270908v1_alt");
        m.insert("KI270912.1", "chr18_KI270912v1_alt");
        m.insert("KI270911.1", "chr18_KI270911v1_alt");
        m.insert("GL949747.2", "chr19_GL949747v2_alt");
        m.insert("KB663609.1", "chr22_KB663609v1_alt");
        m.insert("KI270913.1", "chrX_KI270913v1_alt");
        m.insert("KI270914.1", "chr19_KI270914v1_alt");
        m.insert("KI270915.1", "chr19_KI270915v1_alt");
        m.insert("KI270916.1", "chr19_KI270916v1_alt");
        m.insert("KI270917.1", "chr19_KI270917v1_alt");
        m.insert("KI270918.1", "chr19_KI270918v1_alt");
        m.insert("KI270919.1", "chr19_KI270919v1_alt");
        m.insert("KI270920.1", "chr19_KI270920v1_alt");
        m.insert("KI270921.1", "chr19_KI270921v1_alt");
        m.insert("KI270922.1", "chr19_KI270922v1_alt");
        m.insert("KI270923.1", "chr19_KI270923v1_alt");
        m.insert("KI270924.1", "chr3_KI270924v1_alt");
        m.insert("KI270925.1", "chr4_KI270925v1_alt");
        m.insert("GL000252.2", "chr6_GL000252v2_alt");
        m.insert("KI270926.1", "chr8_KI270926v1_alt");
        m.insert("KI270927.1", "chr11_KI270927v1_alt");
        m.insert("GL949748.2", "chr19_GL949748v2_alt");
        m.insert("KI270928.1", "chr22_KI270928v1_alt");
        m.insert("KI270929.1", "chr19_KI270929v1_alt");
        m.insert("KI270930.1", "chr19_KI270930v1_alt");
        m.insert("KI270931.1", "chr19_KI270931v1_alt");
        m.insert("KI270932.1", "chr19_KI270932v1_alt");
        m.insert("KI270933.1", "chr19_KI270933v1_alt");
        m.insert("GL000209.2", "chr19_GL000209v2_alt");
        m.insert("KI270934.1", "chr3_KI270934v1_alt");
        m.insert("GL000253.2", "chr6_GL000253v2_alt");
        m.insert("GL949749.2", "chr19_GL949749v2_alt");
        m.insert("KI270935.1", "chr3_KI270935v1_alt");
        m.insert("GL000254.2", "chr6_GL000254v2_alt");
        m.insert("GL949750.2", "chr19_GL949750v2_alt");
        m.insert("KI270936.1", "chr3_KI270936v1_alt");
        m.insert("GL000255.2", "chr6_GL000255v2_alt");
        m.insert("GL949751.2", "chr19_GL949751v2_alt");
        m.insert("KI270937.1", "chr3_KI270937v1_alt");
        m.insert("GL000256.2", "chr6_GL000256v2_alt");
        m.insert("GL949752.1", "chr19_GL949752v1_alt");
        m.insert("KI270758.1", "chr6_KI270758v1_alt");
        m.insert("GL949753.2", "chr19_GL949753v2_alt");
        m.insert("KI270938.1", "chr19_KI270938v1_alt");
        m.insert("KI270302.1", "chrUn_KI270302v1");
        m.insert("KI270304.1", "chrUn_KI270304v1");
        m.insert("KI270303.1", "chrUn_KI270303v1");
        m.insert("KI270305.1", "chrUn_KI270305v1");
        m.insert("KI270322.1", "chrUn_KI270322v1");
        m.insert("KI270320.1", "chrUn_KI270320v1");
        m.insert("KI270310.1", "chrUn_KI270310v1");
        m.insert("KI270316.1", "chrUn_KI270316v1");
        m.insert("KI270315.1", "chrUn_KI270315v1");
        m.insert("KI270312.1", "chrUn_KI270312v1");
        m.insert("KI270311.1", "chrUn_KI270311v1");
        m.insert("KI270317.1", "chrUn_KI270317v1");
        m.insert("KI270412.1", "chrUn_KI270412v1");
        m.insert("KI270411.1", "chrUn_KI270411v1");
        m.insert("KI270414.1", "chrUn_KI270414v1");
        m.insert("KI270419.1", "chrUn_KI270419v1");
        m.insert("KI270418.1", "chrUn_KI270418v1");
        m.insert("KI270420.1", "chrUn_KI270420v1");
        m.insert("KI270424.1", "chrUn_KI270424v1");
        m.insert("KI270417.1", "chrUn_KI270417v1");
        m.insert("KI270422.1", "chrUn_KI270422v1");
        m.insert("KI270423.1", "chrUn_KI270423v1");
        m.insert("KI270425.1", "chrUn_KI270425v1");
        m.insert("KI270429.1", "chrUn_KI270429v1");
        m.insert("KI270442.1", "chrUn_KI270442v1");
        m.insert("KI270466.1", "chrUn_KI270466v1");
        m.insert("KI270465.1", "chrUn_KI270465v1");
        m.insert("KI270467.1", "chrUn_KI270467v1");
        m.insert("KI270435.1", "chrUn_KI270435v1");
        m.insert("KI270438.1", "chrUn_KI270438v1");
        m.insert("KI270468.1", "chrUn_KI270468v1");
        m.insert("KI270510.1", "chrUn_KI270510v1");
        m.insert("KI270509.1", "chrUn_KI270509v1");
        m.insert("KI270518.1", "chrUn_KI270518v1");
        m.insert("KI270508.1", "chrUn_KI270508v1");
        m.insert("KI270516.1", "chrUn_KI270516v1");
        m.insert("KI270512.1", "chrUn_KI270512v1");
        m.insert("KI270519.1", "chrUn_KI270519v1");
        m.insert("KI270522.1", "chrUn_KI270522v1");
        m.insert("KI270511.1", "chrUn_KI270511v1");
        m.insert("KI270515.1", "chrUn_KI270515v1");
        m.insert("KI270507.1", "chrUn_KI270507v1");
        m.insert("KI270517.1", "chrUn_KI270517v1");
        m.insert("KI270529.1", "chrUn_KI270529v1");
        m.insert("KI270528.1", "chrUn_KI270528v1");
        m.insert("KI270530.1", "chrUn_KI270530v1");
        m.insert("KI270539.1", "chrUn_KI270539v1");
        m.insert("KI270538.1", "chrUn_KI270538v1");
        m.insert("KI270544.1", "chrUn_KI270544v1");
        m.insert("KI270548.1", "chrUn_KI270548v1");
        m.insert("KI270583.1", "chrUn_KI270583v1");
        m.insert("KI270587.1", "chrUn_KI270587v1");
        m.insert("KI270580.1", "chrUn_KI270580v1");
        m.insert("KI270581.1", "chrUn_KI270581v1");
        m.insert("KI270579.1", "chrUn_KI270579v1");
        m.insert("KI270589.1", "chrUn_KI270589v1");
        m.insert("KI270590.1", "chrUn_KI270590v1");
        m.insert("KI270584.1", "chrUn_KI270584v1");
        m.insert("KI270582.1", "chrUn_KI270582v1");
        m.insert("KI270588.1", "chrUn_KI270588v1");
        m.insert("KI270593.1", "chrUn_KI270593v1");
        m.insert("KI270591.1", "chrUn_KI270591v1");
        m.insert("KI270330.1", "chrUn_KI270330v1");
        m.insert("KI270329.1", "chrUn_KI270329v1");
        m.insert("KI270334.1", "chrUn_KI270334v1");
        m.insert("KI270333.1", "chrUn_KI270333v1");
        m.insert("KI270335.1", "chrUn_KI270335v1");
        m.insert("KI270338.1", "chrUn_KI270338v1");
        m.insert("KI270340.1", "chrUn_KI270340v1");
        m.insert("KI270336.1", "chrUn_KI270336v1");
        m.insert("KI270337.1", "chrUn_KI270337v1");
        m.insert("KI270363.1", "chrUn_KI270363v1");
        m.insert("KI270364.1", "chrUn_KI270364v1");
        m.insert("KI270362.1", "chrUn_KI270362v1");
        m.insert("KI270366.1", "chrUn_KI270366v1");
        m.insert("KI270378.1", "chrUn_KI270378v1");
        m.insert("KI270379.1", "chrUn_KI270379v1");
        m.insert("KI270389.1", "chrUn_KI270389v1");
        m.insert("KI270390.1", "chrUn_KI270390v1");
        m.insert("KI270387.1", "chrUn_KI270387v1");
        m.insert("KI270395.1", "chrUn_KI270395v1");
        m.insert("KI270396.1", "chrUn_KI270396v1");
        m.insert("KI270388.1", "chrUn_KI270388v1");
        m.insert("KI270394.1", "chrUn_KI270394v1");
        m.insert("KI270386.1", "chrUn_KI270386v1");
        m.insert("KI270391.1", "chrUn_KI270391v1");
        m.insert("KI270383.1", "chrUn_KI270383v1");
        m.insert("KI270393.1", "chrUn_KI270393v1");
        m.insert("KI270384.1", "chrUn_KI270384v1");
        m.insert("KI270392.1", "chrUn_KI270392v1");
        m.insert("KI270381.1", "chrUn_KI270381v1");
        m.insert("KI270385.1", "chrUn_KI270385v1");
        m.insert("KI270382.1", "chrUn_KI270382v1");
        m.insert("KI270376.1", "chrUn_KI270376v1");
        m.insert("KI270374.1", "chrUn_KI270374v1");
        m.insert("KI270372.1", "chrUn_KI270372v1");
        m.insert("KI270373.1", "chrUn_KI270373v1");
        m.insert("KI270375.1", "chrUn_KI270375v1");
        m.insert("KI270371.1", "chrUn_KI270371v1");
        m.insert("KI270448.1", "chrUn_KI270448v1");
        m.insert("KI270521.1", "chrUn_KI270521v1");
        m.insert("GL000195.1", "chrUn_GL000195v1");
        m.insert("GL000219.1", "chrUn_GL000219v1");
        m.insert("GL000220.1", "chrUn_GL000220v1");
        m.insert("GL000224.1", "chrUn_GL000224v1");
        m.insert("KI270741.1", "chrUn_KI270741v1");
        m.insert("GL000226.1", "chrUn_GL000226v1");
        m.insert("GL000213.1", "chrUn_GL000213v1");
        m.insert("KI270743.1", "chrUn_KI270743v1");
        m.insert("KI270744.1", "chrUn_KI270744v1");
        m.insert("KI270745.1", "chrUn_KI270745v1");
        m.insert("KI270746.1", "chrUn_KI270746v1");
        m.insert("KI270747.1", "chrUn_KI270747v1");
        m.insert("KI270748.1", "chrUn_KI270748v1");
        m.insert("KI270749.1", "chrUn_KI270749v1");
        m.insert("KI270750.1", "chrUn_KI270750v1");
        m.insert("KI270751.1", "chrUn_KI270751v1");
        m.insert("KI270752.1", "chrUn_KI270752v1");
        m.insert("KI270753.1", "chrUn_KI270753v1");
        m.insert("KI270754.1", "chrUn_KI270754v1");
        m.insert("KI270755.1", "chrUn_KI270755v1");
        m.insert("KI270756.1", "chrUn_KI270756v1");
        m.insert("KI270757.1", "chrUn_KI270757v1");
        m.insert("GL000214.1", "chrUn_GL000214v1");
        m.insert("KI270742.1", "chrUn_KI270742v1");
        m.insert("GL000216.2", "chrUn_GL000216v2");
        m.insert("GL000218.1", "chrUn_GL000218v1");
        m.insert("KI270740.1", "chrY_KI270740v1_random");
        m
    };
}

#[derive(Parser)]
#[command(name = "SPDI")]
#[command(author = "Ryangguk Kim <rkim@oakbioinformatics.com>")]
#[command(version = "0.1.0")]
#[command(about = "SPDI: SPDI format converter")]
#[command(
    after_help = "<variant> is in the form \"chrom:position:reference base:alternate base\" without quotation marks.\nFor example,\n\"chr1:398239:A:C\" for SNV\n\"chr1:26748347:GAC:TA\" for MNV\n\"chr1:2378233:-:A\" for insertion\n\"chr1:72378854:T:-\" for deletion\n\nOutput is an SPDI format string. For example, \"chr1:8734834:GTGT:GT\"\n\nSPDI paper: https://doi.org/10.1093%2Fbioinformatics%2Fbtz856\n\nCopyright 2024 Ryangguk Kim @ Oak Bioinformatics, LLC. Licensed under AGPL-3 and commercial license terms"
)]
struct Cli {
    #[arg(id = "twobit_path")]
    #[arg(
        help = "Path to a 2bit file. What is a 2bit file? See https://genome.ucsc.edu/goldenPath/help/twoBit.html. 2bit files can be downloaded at for example https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/"
    )]
    #[arg(short = 't')]
    twobit_path: String,
    #[arg(help = "Input variant")]
    #[arg(id = "variant")]
    #[arg(short = 'v')]
    #[arg(default_value = "")]
    variant: String,
    #[arg(help = "Variant VCF file")]
    #[arg(id = "input_file")]
    #[arg(short = 'f')]
    #[arg(default_value = "")]
    input_file: String,
}

fn process_variant(variant: &String, spdi: &mut SPDI) {
    let words: Vec<&str> = variant.split(':').collect();
    if words.len() != 4 {
        eprintln!("\nWrong input format: [{}]\n", variant);
        std::process::exit(1);
    }
    let chrom: String = words[0].into();
    let pos_r = words[1].parse::<usize>();
    let pos: usize;
    match pos_r {
        Ok(v) => pos = v,
        Err(_) => {
            eprintln!("\n[{}] is not a valid position.\n", words[1]);
            std::process::exit(1);
        }
    }
    let ref_bases_s: String = words[2].into();
    let alt_bases_s: String = words[3].into();
    let ret = spdi.get_spdi_string(&chrom, pos, &ref_bases_s, &alt_bases_s);
    match ret {
        Err(e) => {
            eprintln!("Error: {:#?}", e);
        }
        Ok(v) => {
            println!("{}", v);
        }
    }
}

fn header_has_sample(line: &String) -> bool {
    let words: Vec<&str> = line.split("\t").collect();
    words.len() > 8
}

fn process_input_file(input_file: &String, spdi: &mut SPDI) {
    let f_r = std::fs::File::open(input_file);
    let f: std::fs::File;
    match f_r {
        Err(_) => {
            eprintln!("Cannot open input file: [{}]", input_file);
            std::process::exit(1);
        }
        Ok(v) => {
            f = v;
        }
    }
    let reader = std::io::BufReader::new(f);
    let mut has_sample: bool = false;
    let mut line: String;
    for line_r in reader.lines() {
        match line_r {
            Err(_) => {
                eprintln!("Error while reading a line from input file");
                std::process::exit(1);
            }
            Ok(v) => {
                line = v;
            }
        }
        if line.starts_with("#CHROM") {
            has_sample = header_has_sample(&line);
            println!("##INFO=<ID=OV_SPDI_IDS,Number=A,Type=String,Description=\"SPDI notation of each alternate allele\">");
            println!("{}", line);
            continue;
        }
        if line.starts_with("#") {
            println!("{}", line);
            continue;
        }
        let words: Vec<&str> = line.split("\t").collect();
        let words_len = words.len();
        if words_len < 8 {
            println!("{}", line);
        }
        let chrom_s = words[0];
        let chrom_1st_c: char = chrom_s.chars().next().unwrap();
        let chrom: String;
        if (chrom_1st_c >= '1' && chrom_1st_c <= '9') || chrom_1st_c == 'X' || chrom_1st_c == 'Y' {
            chrom = format!("chr{}", chrom_s);
        } else if chrom_1st_c == 'M' {
            chrom = "chrM".to_string();
        } else {
            let chrom_o = CHROMS.get(&chrom_s);
            match chrom_o {
                None => {
                    eprintln!("Chromosome [{}] not supported: {}", chrom_s, line);
                    continue;
                },
                Some(v) => {
                    chrom = v.to_string();
                }
            }
        }
        let pos_r = words[1].parse::<usize>();
        let pos: usize;
        match pos_r {
            Err(_) => {
                eprintln!("Invalid POS: {}", line);
                println!("{}", line);
                continue;
            }
            Ok(v) => {
                pos = v;
            }
        }
        let ref_base = words[3];
        let alt_bases = words[4].split(",");
        let mut spdi_strings: Vec<String> = Vec::with_capacity(4);
        for alt_base in alt_bases {
            match spdi.get_spdi_string(&chrom, pos, ref_base, alt_base) {
                Err(e) => {
                    eprintln!("{}: {}", e, line);
                    spdi_strings.push(".".to_string());
                }
                Ok(v) => {
                    spdi_strings.push(v);
                }
            }
        }
        let spdi_string = spdi_strings.join(",");
        let first = words[0..8].join("\t");
        if !has_sample {
            println!("{};OV_SPDI_IDS={}", first, spdi_string);
        } else {
            let last = words[8..].join("\t");
            println!("{};OV_SPDI_IDS={}\t{}", first, spdi_string, last);
        }
    }
}

fn main() {
    let cli = Cli::parse();
    let spdi_r = SPDI::new(&cli.twobit_path);
    let mut spdi: SPDI;
    match spdi_r {
        Err(_) => {
            eprintln!("Cannot open a 2bit file at [{}].", cli.twobit_path);
            std::process::exit(1);
        }
        Ok(v) => spdi = v,
    }
    let variant_len = cli.variant.len();
    let input_file_len = cli.input_file.len();
    match variant_len {
        0 => match input_file_len {
            0 => {
                eprintln!("-v <variant> or -f <input_file> should be given.");
                std::process::exit(1);
            }
            _ => {
                process_input_file(&cli.input_file, &mut spdi);
            }
        },
        _ => process_variant(&cli.variant, &mut spdi),
    }
}

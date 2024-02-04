SPDI is a Rust library and a command-line utility for getting the SPDI-format representation of genomic variants.

## As a command-line utility

```
git clone git@github.com:rkimoakbioinformatics/spdi.git
cd spdi
cargo build --release
# Add SPDI representation as OV_SPDI_IDS field in INFO of a VCF file.
./target/release/spdi -t <2bit file path> -f <VCF file path> 1>out.vcf 2>err.txt
# Get SPDI representation of a single variant.
./target/release/spdi -t <2bit file path> -v chr1:99092:C:CT
```

## As a library
To add:
```
cargo add spdi
cargo add anyhow
```
To use:
```
use anyhow::Result;

fn example() -> Result<String> {
    let s = spdi::SPDI::new("path/to/2bit/file")?;
    let spdi_string = s.get_spdi_string("chr1", 99092, "C", "CT")?;
    println!("{}", spdi_string);
}
```

# Count FASTQ Read barcodes

The `fqgrep` processed a `fastq.gz` file by matching each read to a specified regular expression. If the expression matches, the read is assigned a barcode based on the
matched barcode. Once all reads have been processed, the frequency of all barcodes is returned.

## Usage

~~~plain
Usage: fqgrep [OPTIONS] <REGEX> <EXPR> <FILE>

Arguments:
  <REGEX>  Search expresion
  <EXPR>   Replacement expression
  <FILE>   Input fastq.gz file

Options:
  -v, --verbose...        Show log messages. Multiple -v options increase the verbosity
  -c, --sort              Sort returned barcodes by count
  -n, --unmatched <FILE>  Write non-barcoded sequences to file
  -h, --help              Print help
~~~

Each read is matched against the search regular expression `REGEX`. If no match is found, the read is classed as `no_barcode`. If a match is found, the read barcode is calculated by using the replacement expression `EXPR` on the match.

For example, if using the `REGEX` match expression of `ATAATACGACTCACTATAAAACTGGAAG(.{20})` to look for 20-mer barcodes immediately following a truseq adapter, we would use the `EXPR` replacement expression `${1}` to define the returned barcode as the captured 20-mer sequence.  This would label a read as follows:

~~~plain
INPUT:   ATAATACGACTCACTATAAAACTGGAAGAACGCTGACCACAAGTTCGACCTCCCTTTAGTGAGGGTTAATTTGGAAACCGACGCCCCAGCACTCGTCCGAGGGCAAAGGAATAGTACTGTCTCTTATACACATCTCCGAGCCCACGAGAC
CAPTURE: ATAATACGACTCACTATAAAACTGGAAGAACGCTGACCACAAGTTCGA------------------------------------------------------------------------------------------------------
BARCODE: ----------------------------AACGCTGACCACAAGTTCGA------------------------------------------------------------------------------------------------------
~~~

If `-n` is specified, the sequences (*not* the `.fastq` reads) of non-matching reads are written to the specified file.

## Installation from Source

Before installation, you'll need to install [Rust](https://www.rust-lang.org/).

Once you have rust installed, installation is as follows:

~~~bash
git clone https://github.com/uoy-research/fqgrep.git
cd fqgrep
cargo install --path .
~~~

## Licence

This tool is released under the [MIT License](https://opensource.org/licenses/MIT).

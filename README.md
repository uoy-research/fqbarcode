# Count FASTQ Read barcodes

The `fqbarcode` processes a `fastq.gz` file by matching each read to a specified regular expression. If the expression matches, the read is assigned a barcode based on the
matched regular expression capture groups. Once all reads have been processed, the frequency of all barcodes is returned.

## Usage

~~~plain
Usage: fqbarcode [OPTIONS] <REGEX> <FILE>

Arguments:
  <REGEX>  Search expresion
  <FILE>   Input fastq.gz file

Options:
  -v, --verbose...              Show log messages. Multiple -v options increase the verbosity
  -c, --sort                    Sort returned barcodes by count
  -n, --unmatched <FILE>        Write non-barcoded sequences to file
  -r, --replacement <EXPR>      Replacement expression [default: ${1}]
  -m, --merge-count <N>         Threshold count for merging [default: 0]
  -t, --threshold-distance <D>  Threshold edit distance for merging [default: 1]
  -h, --help                    Print help
~~~

Each read is matched against the search regular expression `REGEX`. If no match is found, the read is classed as `no_barcode`. If a match is found, the read barcode is calculated by using the replacement expression `EXPR` on the match.

For example, if using the `REGEX` match expression of `ATAATACGACTCACTATAAAACTGGAAG(.{20})` to look for 20-mer barcodes immediately following a truseq adapter, we would use the `EXPR` replacement expression `${1}` to define the returned barcode as the captured 20-mer sequence.  This would label a read as follows:

~~~plain
INPUT:   ATAATACGACTCACTATAAAACTGGAAGAACGCTGACCACAAGTTCGACCTCCCTTTAGTGAGGGTTAATTTGGAAACCGACGCCCCAGCACTCGTCCGAGGGCAAAGGAATAGTACTGTCTCTTATACACATCTCCGAGCCCACGAGAC
CAPTURE: ATAATACGACTCACTATAAAACTGGAAGAACGCTGACCACAAGTTCGA------------------------------------------------------------------------------------------------------
BARCODE: ----------------------------AACGCTGACCACAAGTTCGA------------------------------------------------------------------------------------------------------
~~~

## Barcode Merging

As sequencing is an imperfect system, it is common to get a long tail of low-count barcodes that are simple transversions of other more common barcodes.  To address this, the optional `--merge-count` (`-m`) argument can be used to preform barcode merging.

If `-m` is assigned a value greater than zero, then barcode merging will be attempted.  The algorithm is as follows:

1. All (barcode, count) pairs are ordered based on descending count refequency;
2. All barcodes with more than `merge-count` reads are marked as endpoints, and will not be subject to merging;
3. All remaining barcodes are processed in ascending order of frequency:
   1. The [Levenshtein](https://en.wikipedia.org/wiki/Levenshtein_distance) edit distance between the barcode and all endpoints is calculated;
   2. The minimum edit distance is calculated (if the minimum edit distance is > `threshold-distance`, then do not merge);
   3. A single endpoint is selected at random from the set of endpoints having the minimum edit distance;
   4. The barcode counts are added to the selected endpoint counts, and the original barcode is removed.

## Outputs

* The tab-delimited read counts and barcodes are returned in decreasing order of frequency;
* The number of reads that did not match the input `REGEX` is returned as `no_barcode`
* If `-n` is specified, the sequences (*not* the `.fastq` reads) of non-matching reads are written to the specified file.

## Installation from Source

Before installation, you'll need to install [Rust](https://www.rust-lang.org/).

Once you have rust installed, installation is as follows:

~~~bash
git clone https://github.com/uoy-research/fqbarcode.git
cd fqbarcode
cargo install --path .
~~~

## Licence

This tool is released under the [MIT License](https://opensource.org/licenses/MIT).

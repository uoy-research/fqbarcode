use clap::Parser;
use flate2::read::MultiGzDecoder;
use log::*;
use regex::Regex;
use simple_eyre::eyre::Report;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

#[derive(Parser)]
struct Args {
    /// Show log messages. Multiple -v options increase the verbosity
    #[clap(short='v', long="verbose", action=clap::ArgAction::Count)]
    verbose: u8,
    /// Sort returned barcodes by count
    #[clap(short = 'c', long = "sort")]
    sort_barcodes: bool,
    /// Write non-barcoded sequences to file
    #[clap(short = 'n', long = "unmatched", value_name = "FILE")]
    unmatched_path: Option<PathBuf>,
    /// Search expresion
    #[clap(value_name = "REGEX")]
    barcode_expression: String,
    /// Replacement expression
    #[clap(value_name = "EXPR")]
    barcode_replacement: String,
    /// Input fastq.gz file
    #[clap(value_name = "FILE")]
    file_path: PathBuf,
}

fn main() -> Result<(), Report> {
    // Register the Eyre handler:
    simple_eyre::install()?;
    // Parse the CLI arguments:
    let args = Args::parse();
    // Build the logger:
    stderrlog::new()
        .module(module_path!())
        .verbosity(args.verbose as usize)
        .timestamp(stderrlog::Timestamp::Millisecond)
        .init()?;
    // Load the input file:
    info!("parsing reads from {}", args.file_path.to_string_lossy());
    let input_file = BufReader::new(File::open(&args.file_path)?);
    let input_buffer = BufReader::new(MultiGzDecoder::new(input_file));
    // Build the regular expression:
    debug!("building barcode regular expression");
    trace!("barcode regular expression is {}", args.barcode_expression);
    let barcode_re = Regex::new(&args.barcode_expression)?;
    // Define the barcode counts:
    let mut barcode_label = String::new();
    let mut barcode_counts: BTreeMap<String, u64> = BTreeMap::new();
    let mut total_count: u64 = 0_u64;
    let mut no_barcode: u64 = 0_u64;
    // Loop over all reads:
    debug!("processing reads");
    for read in input_buffer.lines().skip(1).step_by(4) {
        let read = read?;
        match barcode_re.captures(&read) {
            Some(c) => {
                c.expand(&args.barcode_replacement, &mut barcode_label);
                // let barcode = c.get(1).ok_or(eyre!("invalid capture group"))?.as_str().to_string();
                trace!("read {} barcode label is {}", read, barcode_label);
                let i = barcode_counts
                    .entry(barcode_label.to_owned())
                    .or_insert(0_u64);
                *i += 1_u64;
                barcode_label.clear();
            }
            None => {
                trace!("no barcode detected in read {}", read);
                no_barcode += 1_u64;
            }
        }
        total_count += 1_u64;
    }
    // Print out the results:
    info!("{} reads processed", total_count);
    if args.sort_barcodes {
        let mut barcodes_sorted: Vec<(&String, &u64)> = barcode_counts.iter().collect();
        barcodes_sorted.sort_by(|a, b| b.1.cmp(a.1));
        for (barcode, count) in barcodes_sorted {
            println!("{}\t{}", count, barcode);
        }
    } else {
        for (barcode, count) in barcode_counts.iter() {
            println!("{}\t{}", count, barcode);
        }
    }
    println!("{}\tno_barcode", no_barcode);
    Ok(())
}

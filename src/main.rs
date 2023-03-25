use clap::Parser;
use flate2::read::MultiGzDecoder;
use levenshtein::levenshtein;
use log::*;
use rand::seq::SliceRandom;
use rand::thread_rng;
use regex::Regex;
use simple_eyre::eyre::Report;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;

#[derive(Parser)]
#[command(version)]
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
    #[clap(
        short = 'r',
        long = "replacement",
        value_name = "EXPR",
        default_value = "${1}"
    )]
    barcode_replacement: String,
    /// Threshold count for merging
    #[clap(
        short = 'm',
        long = "merge-count",
        value_name = "N",
        default_value = "0"
    )]
    threshold_count: u64,
    /// Threshold edit distance for merging
    #[clap(
        short = 't',
        long = "threshold-distance",
        value_name = "D",
        default_value = "1"
    )]
    threshold_distance: usize,
    /// Input fastq.gz file
    #[clap(value_name = "FILE")]
    file_path: PathBuf,
}

fn sort_barcodes(s: &mut [(String, u64)]) {
    s.sort_by(|a, b| b.1.cmp(&a.1));
}

fn count_barcodes(m: &HashMap<String, u64>) -> u64 {
    m.iter().map(|(_, i)| *i).sum()
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
    // Set up the RNG:
    let mut rng = thread_rng();
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
    let mut barcodes: HashMap<String, u64> = HashMap::new();
    let mut total_reads: u64 = 0_u64;
    let mut no_barcode: u64 = 0_u64;
    // If requested, open the non-barcode sequence file:
    let mut unmatched_buffer = match args.unmatched_path {
        Some(unmatched_path) => {
            info!(
                "writing non-barcoded sequences to {}",
                unmatched_path.to_string_lossy()
            );
            Some(BufWriter::new(File::create(unmatched_path)?))
        }
        None => None,
    };
    // Loop over all reads:
    debug!("processing reads");
    for read in input_buffer.lines().skip(1).step_by(4) {
        let read = read?;
        match barcode_re.captures(&read) {
            Some(c) => {
                c.expand(&args.barcode_replacement, &mut barcode_label);
                trace!("read {} barcode label is {}", read, barcode_label);
                let i = barcodes.entry(barcode_label.to_owned()).or_insert(0_u64);
                *i += 1_u64;
                barcode_label.clear();
            }
            None => {
                trace!("no barcode detected in read {}", read);
                if let Some(ref mut buffer) = unmatched_buffer {
                    writeln!(buffer, "{read}")?;
                }
                no_barcode += 1_u64;
            }
        }
        total_reads += 1_u64;
    }
    info!("processed {total_reads} reads");
    info!(
        "{}/{} ({:0.2}%) reads did not match barcode",
        no_barcode,
        total_reads,
        (no_barcode as f32 / total_reads as f32) * 100_f32
    );
    info!("{} barcodes detected", barcodes.len());

    // Now we have all the barcodes we can extract a list of the "endpoint" barcodes, i.e.
    // those that can accept merged barcodes. A barcode is an endpoint if it currently has
    // at least args.threshold_count reads associated to it.
    let endpoint_barcodes: HashSet<String> = barcodes
        .iter()
        .filter_map(|(barcode, count)| match count > &args.threshold_count {
            false => None,
            true => Some(barcode.to_owned()),
        })
        .collect();
    debug!("{} barcodes pass threshold count", endpoint_barcodes.len());

    if !endpoint_barcodes.is_empty() {
        // Get a list of the non-endpoint barcodes sorted by their count (lowest first):
        let mut non_endpoint_barcodes: Vec<(String, u64)> = barcodes
            .iter()
            .filter_map(|(barcode, count)| {
                let barcode = barcode.to_owned();
                match endpoint_barcodes.contains(&barcode) {
                    false => Some((barcode, *count)),
                    true => None,
                }
            })
            .collect();
        non_endpoint_barcodes.sort_by(|a, b| a.1.cmp(&b.1));

        // Iterate through each of the non-endpoint barcodes, and attempt to mege it into a single on of the endpoints.
        for (barcode, count) in non_endpoint_barcodes.iter() {
            let barcode = barcode.to_owned();
            debug!(
                "barcode {barcode} count {count} <= {}; attempting to merge",
                &args.threshold_count
            );
            // Get the edit distances between this barcode and all endpoints:
            let mut end_point_distances: HashMap<String, usize> =
                HashMap::with_capacity(endpoint_barcodes.len());
            for endpoint_barcode in endpoint_barcodes.iter() {
                end_point_distances.insert(
                    endpoint_barcode.to_owned(),
                    levenshtein(&barcode, endpoint_barcode),
                );
            }
            // Find the minimum edit distance:
            let min_endpoint_distance: usize = end_point_distances
                .values()
                .copied()
                .min()
                .unwrap_or(0_usize);
            // Only continue if the minimum distance is suitably low:
            if min_endpoint_distance <= args.threshold_distance {
                // Get a set of all the endpoints with the minimum distance:
                let min_distance_endpoint_barcodes: Vec<String> = end_point_distances
                    .iter()
                    .filter_map(
                        |(barcode, distance)| match distance == &min_endpoint_distance {
                            false => None,
                            true => Some(barcode.to_owned()),
                        },
                    )
                    .collect();
                // Select a single endpoint from the available options:
                if let Some(selected_endpoint) = min_distance_endpoint_barcodes.choose(&mut rng) {
                    let selected_endpoint = selected_endpoint.to_owned();
                    debug!("merging barcode {barcode} (count={count}) into {selected_endpoint} (distance is {min_endpoint_distance})");
                    // Move across the merged counts to the endpoint:
                    *barcodes.entry(selected_endpoint).or_insert(0) += count;
                    // Delete this now-merged barcode:
                    barcodes.remove(&barcode);
                }
            } else {
                debug!("barcode {barcode} minumum edit distance ({min_endpoint_distance}) is too great; not merging");
            }
        }
    } else {
        info!(
            "no barcodes have counts > {}; merging not performed",
            args.threshold_count
        );
    }

    // Print out the results:
    info!("{} reads assigned a barcode", count_barcodes(&barcodes));
    info!("{} barcodes remain after merging", barcodes.len());
    let mut barcodes: Vec<(String, u64)> =
        barcodes.iter().map(|(s, i)| (s.to_owned(), *i)).collect();
    sort_barcodes(&mut barcodes);
    for (barcode, count) in barcodes {
        println!("{}\t{}", count, barcode);
    }
    println!("{no_barcode}\tno_barcode");
    Ok(())
}

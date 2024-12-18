use clap::{Parser, ValueEnum};
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::io::{BufReader, BufRead};
use std::fs::OpenOptions;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::io::BufWriter;
use std::io::{self, Read};
use std::path::Path;
use std::time::Instant;
use smitten::{Identifier, IDVersion};
use bio::io::fasta;
use regex::Regex;
use memmap2::Mmap;

/// Command-line arguments
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Path to the reference sequence file [twobit, or fasta]
    #[arg(short, long)]
    reference: String,

    /// Enable mapping (Boyer-Moore) for invalid identifiers
    #[arg(short, long, default_value = "false")]
    map_sequences: bool,

    // TODO: Change this to Option<String>
    /// Optional output file path
    #[arg(short, long)]
    output: Option<String>,

    /// Log level (Summary, PerRecord or Detailed)
    #[arg(short, long, default_value = "summary")]
    log_level: LogLevel,

    /// Input file path (Fasta, Stockholm or Tab/Comma Delimited file)
    #[arg()]
    input: String,

    /// Threads to use for parallel processing
    #[arg(short = 'x', long)]
    threads: Option<usize>,
}

// Define a structure to store parsed ID components and sequences
#[derive(Debug)]
struct SequenceRecord {
    original_id: Option<String>,          // optional
    assembly_id: Option<String>,          // optional
    sequence_id: String,                  // Mandatory
    start: Option<u64>,                   // optional
    end: Option<u64>,                     // optional
    orient: Option<char>,                 // optional
    inferred_version: Option<IDVersion>,  // optional
    sequence: Vec<u8>,                    // Mandatory
    aligned_seq: Option<Vec<u8>>,         // optional
    validated: Option<String>,            // optional
}

impl SequenceRecord {
    fn print_record(&self) {
        println!(
            "Smitten::Identifier: original_id: {}, assembly_id: {}, sequence_id: {}, start: {}, end: {}, orient: {}, inferred_version: {:?}, validated: {}",
            self.original_id.as_deref().unwrap_or("Unknown"),
            self.assembly_id.as_deref().unwrap_or("None"),
            self.sequence_id,
            self.start.unwrap_or(0),  // Default to 0 if None
            self.end.unwrap_or(0),    // Default to 0 if None
            self.orient.unwrap_or('?'),
            self.inferred_version,
            self.validated.as_deref().unwrap_or(""),
        );
    }
}

#[derive(Clone, Debug, ValueEnum, PartialEq)]
enum LogLevel {
    Summary,
    PerRecord,
    Detailed,
}


fn parse_and_validate_stockholm(file_path: &str, genome_map: &HashMap<String, Vec<u8>>, output_path: &Option<String>, is_gzip: bool, log_level: LogLevel, map_seqs: bool) -> io::Result<()> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    let mut current_record: HashMap<String, String> = HashMap::new();
    let mut metadata: Vec<String> = Vec::new();

    // For summary stats
    let mut combined_status_counts: HashMap<String, usize> = HashMap::new();
    let mut total_seq_records = 0;

    // Regex to match the ID/ACC lines
    let re_id = Regex::new(r"^#=GF ID\s+(.+)$").unwrap();
    let re_ac = Regex::new(r"^#=GF AC\s+(.+)$").unwrap();


    let mut ident = String::new();
    let mut accession = String::new();
    let mut record_number = 1;
    for line in reader.lines() {
        let line = line?;

        if line.is_empty() {
            continue;
        }

        if line.starts_with("#") {
            metadata.push(line.clone());
            if re_id.is_match(&line) {
                ident = re_id.captures(&line).unwrap().get(1).unwrap().as_str().to_string();
            } else if re_ac.is_match(&line) {
                accession = re_ac.captures(&line).unwrap().get(1).unwrap().as_str().to_string();
            }
            continue;
        }

        if line.starts_with("//") {
            let label = if !accession.is_empty() {
                accession.clone()
            } else if !ident.is_empty() { 
                ident.clone()
            }else {
                format!("record_{}", record_number)
            };
            //println!("Processing record: {}", label);
            let (rec_seq_count, rec_status_counts) = process_stockholm_record(&mut current_record, &metadata, genome_map, output_path, is_gzip, &log_level, map_seqs, label);
            for (key, value) in rec_status_counts {
                *combined_status_counts.entry(key.clone()).or_insert(0) += value; 
            }  
            total_seq_records += rec_seq_count;
            current_record.clear();
            record_number += 1;
        } else {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() == 2 {
                let (name, seq) = (parts[0].to_string(), parts[1].to_string());
                current_record.entry(name).or_insert(String::new()).push_str(&seq);
            }
        }
    }

    if !current_record.is_empty() {
        panic!("Incomplete record found in Stockholm file");
    }

    // Output summary stats
    println!("-------------------------------------------------");
    println!("Summary for {}:", file_path);
    println!("  Total Sequences: {}", total_seq_records);
    if combined_status_counts.get("valid").is_some() {
        println!("     Accurate Coordinates: {}", combined_status_counts.get("valid").unwrap());
    }else {
        println!("     Accurate Coordinates: 0");
    }
    let mut repaired_total = 0;
    for (fix_type, count) in combined_status_counts.iter() {
        if fix_type != "valid" && fix_type != "invalid" {
            repaired_total += count;
        }
    }
    println!("     Repaired Coordinates: {}", repaired_total);
    for (fix_type, count) in combined_status_counts.iter() {
        if fix_type == "valid" || fix_type == "invalid" {
            continue;
        }
        println!("        {}: {}", fix_type, count);
    }
    if combined_status_counts.get("invalid").is_some() {
        println!("     Invalid Coordinates: {}", combined_status_counts.get("invalid").unwrap());
    }else {
        println!("     Invalid Coordinates: 0");
    }


    Ok(())
}

    
fn process_stockholm_record(
    current_record: &mut HashMap<String, String>,
    metadata: &[String],
    genome_map: &HashMap<String, Vec<u8>>,
    output_path: &Option<String>, 
    is_gzip: bool,
    log_level: &LogLevel,
    map_seqs: bool,
    label: String,
) -> (usize, HashMap<String, usize>){
    let mut records = Vec::new();

    for (name, seq) in current_record.drain() {
        let (v2_id, inferred_version) = Identifier::from_unknown_format(&name, false, true)
            .expect("Failed to convert ID to V2");
        let normalized_parsed_id = v2_id.normalize().unwrap();

        let assembly_id = normalized_parsed_id.assembly_id.clone();
        let sequence_id = normalized_parsed_id.sequence_id.clone();
        let orient = normalized_parsed_id.ranges.get(0).map(|r| r.orientation.clone());
        let start = normalized_parsed_id.ranges.get(0).map(|r| r.start as u64);
        let end = normalized_parsed_id.ranges.get(0).map(|r| r.end as u64);

        let raw_seq: String = seq.chars().filter(|&c| c != '.' && c != '-').collect();

        records.push(SequenceRecord {
            original_id: Some(name),
            assembly_id,
            sequence_id,
            start,
            end,
            orient,
            inferred_version: Some(inferred_version),
            sequence: raw_seq.into_bytes(),
            aligned_seq: Some(seq.into_bytes()),
            validated: None,
        });
    }

    // Validate the records
    validate_sequences(&mut records, genome_map, false);
    
    if map_seqs {
        // Run Boyer-Moore search on unvalidated records
        boyer_moore_search_with_validation(&mut records, genome_map, false);
    }

    let total_records = records.len();
    let mut status_counts = HashMap::new();
    for record in records.iter_mut() {
        if record.validated.is_none() {
            record.validated = Some("invalid".to_string());
        }
        *status_counts.entry(record.validated.clone().unwrap()).or_insert(0) += 1;
    }

    // Log what we did
    if *log_level == LogLevel::PerRecord || *log_level == LogLevel::Detailed{
        output_results(&records, log_level.clone(), label.clone());
    }

    // Write the output in append mode

    if let Some(outpath) = output_path {
        write_stockholm_output(&records, metadata, outpath, is_gzip, true).expect("Failed to write Stockholm output");
    }

    (total_records, status_counts)
}

fn write_stockholm_output(
    records: &[SequenceRecord],
    metadata: &[String],
    output_path: &str,
    is_gzip: bool,
    append: bool
) -> io::Result<()> {
    let file = OpenOptions::new()
        .write(true)
        .create(true)
        .append(append)
        .open(output_path)?;

    let mut writer: Box<dyn Write> = if is_gzip {
        Box::new(BufWriter::new(GzEncoder::new(file, Compression::default())))
    } else {
        Box::new(BufWriter::new(file))
    };

    for line in metadata {
        writeln!(writer, "{}", line)?;
    }
    for record in records {
        let v2_id = if let Some(assembly) = &record.assembly_id {
            format!("{}:{}", assembly, record.sequence_id)
        } else {
            record.sequence_id.clone()
        };
        if record.start.is_some() && record.end.is_some() && record.orient.is_some() {
            writeln!(writer, "{}:{}-{}_{} {}", v2_id, record.start.unwrap(), record.end.unwrap(), 
                    record.orient.unwrap(), String::from_utf8_lossy(&record.aligned_seq.as_ref().unwrap()))?;
        } else {
            writeln!(writer, "{} {}", v2_id, String::from_utf8_lossy(&record.aligned_seq.as_ref().unwrap()))?;
        }
    }
    writeln!(writer, "//")?;

    writer.flush()?;
    Ok(())
}

fn write_fasta_output(
    records: &[SequenceRecord],
    output_path: &str,
    is_gzip: bool,
    append: bool
) -> io::Result<()> {
    let file = OpenOptions::new()
        .write(true)
        .create(true)
        .append(append)
        .open(output_path)?;

    let mut writer: Box<dyn Write> = if is_gzip {
        Box::new(BufWriter::new(GzEncoder::new(file, Compression::default())))
    } else {
        Box::new(BufWriter::new(file))
    };

    for record in records {
        // TODO: Use display functionality now
        let v2_id = if let Some(assembly) = &record.assembly_id {
            format!("{}:{}", assembly, record.sequence_id)
        } else {
            record.sequence_id.clone()
        };
        if record.start.is_some() && record.end.is_some() && record.orient.is_some() {
            writeln!(writer, ">{}:{}-{}_{}", v2_id, record.start.unwrap(), record.end.unwrap(), 
                    record.orient.unwrap())?;
        } else {
            writeln!(writer, ">{}", v2_id)?;
        }
        writeln!(writer, "{}", String::from_utf8_lossy(&record.sequence))?;
    }

    writer.flush()?;
    Ok(())
}

fn write_delimited_output(
    records: &[SequenceRecord],
    output_path: &str,
    is_gzip: bool,
    append: bool,
    format: &str,
) -> io::Result<()> {
    // Determine the delimiter based on the format
    let delimiter = match format {
        "TabDelimited" => '\t',
        "CommaDelimited" => ',',
        _ => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "Unsupported format. Use 'TabDelimited' or 'CommaDelimited'.",
            ))
        }
    };

    let file = OpenOptions::new()
        .write(true)
        .create(true)
        .append(append)
        .open(output_path)?;

    let mut writer: Box<dyn Write> = if is_gzip {
        Box::new(BufWriter::new(GzEncoder::new(file, Compression::default())))
    } else {
        Box::new(BufWriter::new(file))
    };

    for record in records {
        // Extract components for delimited output
        let assembly_id = record.assembly_id.clone().unwrap_or_default();
        let sequence_id = record.sequence_id.clone();
        let start = record.start.map(|v| v.to_string()).unwrap_or_default();
        let end = record.end.map(|v| v.to_string()).unwrap_or_default();
        let orient = record.orient.clone().unwrap_or_default();
        let sequence = String::from_utf8_lossy(&record.sequence);

        // Write the record as a single line with the selected delimiter
        writeln!(
            writer,
            "{}{}{}{}{}{}{}{}{}{}{}",
            assembly_id, delimiter,
            sequence_id, delimiter,
            start, delimiter,
            end, delimiter,
            orient, delimiter,
            sequence
        )?;
    }

    writer.flush()?;
    Ok(())
}


fn detect_format_and_compression(path: &str) -> io::Result<(bool, &str)> {
    let file = File::open(path)?;
    let mut buf_reader = BufReader::new(file);

    // Check if the file is gzip-compressed
    let is_gzip = match buf_reader.fill_buf() {
        Ok(data) => data.starts_with(b"\x1F\x8B"), // Gzip magic number
        Err(_) => false,
    };

    let mut reader: Box<dyn BufRead> = if is_gzip {
        let file = File::open(path)?; // Re-open the file for gzip decompression
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(buf_reader)
    };

    let mut magic_buffer = [0; 4];
    reader.read_exact(&mut magic_buffer)?;
    let magic_be = u32::from_be_bytes(magic_buffer);
    let magic_le = u32::from_le_bytes(magic_buffer);
    if magic_be == 0x1A412743 || magic_le == 0x1A412743 {
        return Ok((is_gzip, "TwoBit"));
    }

    // Convert the magic bytes into a string
    let mut first_line = String::from_utf8_lossy(&magic_buffer).to_string();
    let mut format = None;
    for _ in 0..15 {
        let mut line = String::new();

        // If this is the first iteration, start with the prepended magic bytes
        if !first_line.is_empty() {
            line.push_str(&first_line);
            first_line.clear(); // Clear after the first use
        } 
        if reader.read_line(&mut line)? == 0 {
            break; // EOF reached
        }

        // Skip blank or whitespace-only lines
        if line.trim().is_empty() {
            continue;
        }

        // Check for FASTA or Stockholm format based on line content
        if line.starts_with(">") {
            format = Some("Fasta");
            break;
        } else if line.starts_with("# STOCKHOLM 1.0") {
            format = Some("Stockholm");
            break;
        }

        // Detect TabDelimited or CommaDelimited files
        if line.contains('\t') && line.split('\t').count() > 1 {
            format = Some("TabDelimited");
            break;
        } else if line.contains(',') && line.split(',').count() > 1 {
            format = Some("CommaDelimited");
            break;
        }
    }

    match format {
        Some(fmt) => Ok((is_gzip, fmt)),
        None => Err(io::Error::new(io::ErrorKind::InvalidData, "Unknown file format")),
    }
}

fn validate_sequences(
    records: &mut [SequenceRecord],
    genome_map: &HashMap<String, Vec<u8>>,
    debug_mode: bool,
) {
    let mut fix_counts: HashMap<String, usize> = HashMap::new();

    for record in records.iter_mut() {
        //record.print_record();

        let genome_sequence = match genome_map.get(&record.sequence_id) {
            Some(seq) => seq,
            None => continue, // Skip if the sequence ID is not found in the genome
        };

        // Full-length sequences ( no range specified)
        if record.start.is_none() {
            // Is it the same as the database sequence?
            if genome_sequence == &record.sequence {
                record.validated = Some("valid".to_string());
                *fix_counts.entry(record.validated.clone().unwrap()).or_insert(0) += 1;
            }
            continue;
        }

        // Calculate the length of the range in the record
        let range_length = match (record.end, record.start) {
            (Some(end), Some(start)) => end - start,
             _ => 0, // or an appropriate default if either end or start is None
        };
        let fasta_sequence_length = record.sequence.len() as u64;

        // Check if the coordinates are half-open (range matches sequence length exactly)
        let mut validation_str = String::new();
        if range_length == fasta_sequence_length {
            if debug_mode {
                println!(
                    "Detected half-open coordinates for record {}. Converting to one-based fully closed.",
                    record.original_id.as_ref().unwrap()
                );
            }
            // Convert to one-based fully closed by incrementing the start position
            // Note: This is working under the assumption that it is probably zero-based-half-open
            //       Although, it won't matter as this will be followed up by looking at shifts.
            validation_str.push_str("_halfopen");
            record.start = record.start.map(|start| start + 1);
        }

        let start = record.start.unwrap() as usize - 1;
        let end = record.end.unwrap() as usize;
        let fasta_sequence = &record.sequence;
        let rev_complement = reverse_complement(fasta_sequence);
        let mut located = false;

        //println!("  Working on start {} and end {}", start, end);
        //println!("  Working on genome_sequence {:?}", genome_sequence[start..end].to_vec());
        //println!("  Working on fasta_sequence {:?}", fasta_sequence);
        //println!("  Working on REV fasta_sequence {:?}", rev_complement);

        // 1. Look for direct matches (either orientation)
        if start < genome_sequence.len() && end <= genome_sequence.len() {
            let mut direct_match_orient: Option<char> = None;
            // Direct match on fwd strand
            if &genome_sequence[start..end] == fasta_sequence {
                direct_match_orient = Some('+');
            // Direct match on rev strand
            }
            if &genome_sequence[start..end] == &rev_complement {
                if direct_match_orient.is_none() {
                    direct_match_orient = Some('-');
                }else {
                    // Both orientations match, so we can't determine the correct orientation
                    direct_match_orient = Some('B');
                }
            }
            if direct_match_orient.is_some() {
                located = true;
                if direct_match_orient == Some('B') || direct_match_orient == record.orient {
                    if debug_mode {
                        println!("Direct match validated for: {:?}", record);
                    }
                } else
                {
                    validation_str.push_str("_orient");
                    record.orient = direct_match_orient;
                }
            }
        }

        // 2. Shifted matches on both strands within +/-3bp range
        if !located {
            let shifts: [isize; 6] = [-3, -2, -1, 1, 2, 3];
            let orig_len = end.saturating_sub(start);
            for shift in shifts.iter() {
                let shifted_start = if *shift < 0 {
                    start.saturating_sub((-*shift) as usize)
                } else {
                    start.saturating_add(*shift as usize)
                };

                let shifted_end = if *shift < 0  {
                    end.saturating_sub((-*shift) as usize)
                } else {
                    end.saturating_add(*shift as usize)
                };
                let new_len = shifted_end.saturating_sub(shifted_start);

                if new_len == orig_len && shifted_end  <= genome_sequence.len() {
                    // Check forward orientation with shift
                    if &genome_sequence[shifted_start..shifted_end] == fasta_sequence {
                        validation_str.push_str(&format!("{}{}{}", 
                            if record.orient == Some('-') { "_orient" } else { "" },
                            if *shift >= 0 { "_plus" } else { "_minus" },
                            shift.abs()));
                        record.start = Some((shifted_start + 1) as u64);
                        record.end = Some(shifted_end as u64);
                        record.orient = if record.orient == Some('-') { Some('+') } else { Some('+') };
                        located = true;
                        break;
                    }

                    // Check reverse complement with shift
                    if &genome_sequence[shifted_start..shifted_end] == &rev_complement {
                        validation_str.push_str(&format!("{}{}{}", 
                            if record.orient == Some('+') { "_orient" } else { "" },
                            if *shift >= 0 { "_plus" } else { "_minus" },
                            shift.abs()));
                        record.start = Some((shifted_start + 1) as u64);
                        record.end = Some(shifted_end as u64);
                        record.orient = if record.orient == Some('+') { Some('-') } else { Some('-') };
                        located = true;
                        break;
                    }
                }
            }
        }

        if located {
            if validation_str.is_empty() {
                record.validated = Some("valid".to_string());
            } else {
                record.validated = Some(format!("fixed{}",validation_str));
            }
            //println!("Record validated: {:?}", record.validated);
            *fix_counts.entry(record.validated.clone().unwrap()).or_insert(0) += 1;
        }
    }
}



// Generates the reverse complement of a DNA sequence
fn reverse_complement(dna: &[u8]) -> Vec<u8> {
    dna.iter()
        .rev()
        .map(|&base| match base {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => base,
        })
        .collect()
}

fn load_and_parse_fasta(fasta_path: &str) -> Vec<SequenceRecord> {
    let reader = fasta::Reader::from_file(fasta_path).expect("Could not read FASTA file");
    let mut parsed_records = Vec::new();

    for rec in reader.records() {
        let rec = rec.expect("Could not parse FASTA record");
        let original_id = rec.id().to_string();

        // Convert the identifier to V2 format and capture inferred version
        let (v2_id, inferred_version) = Identifier::from_unknown_format(&original_id, false, true).expect("Failed to convert ID to V2");
        //println!("Original ID: {}, V2 ID: {}, Inferred Version: {}", original_id, v2_id, inferred_version);

        // Parse and normalize the V2 identifier
        let normalized_parsed_id = v2_id.normalize().unwrap();

        // Extract components
        let assembly_id = normalized_parsed_id.assembly_id.clone();
        let sequence_id = normalized_parsed_id.sequence_id.clone();
        let orient = normalized_parsed_id.ranges.get(0).map(|range| range.orientation.clone());
        let start = normalized_parsed_id.ranges.get(0).map(|range| range.start as u64);
        let end = normalized_parsed_id.ranges.get(0).map(|range| range.end as u64);

        // Store the parsed ID information along with the sequence and inferred version
        parsed_records.push(SequenceRecord {
            original_id: Some(original_id),                 // Original ID
            assembly_id,
            sequence_id,
            start: start,
            end: end,
            orient: orient,
            inferred_version: Some(inferred_version),  // Inferred version
            sequence: rec.seq().to_owned().to_ascii_uppercase(),
            aligned_seq: None,
            validated: None,
        });
    }
    parsed_records
}

fn load_and_parse_delimited_file(svs_path: &str) -> Vec<SequenceRecord> {
    // Open the file for reading
    let file = File::open(Path::new(svs_path)).expect("Could not read SVS file");
    let reader = io::BufReader::new(file);

    let mut parsed_records = Vec::new();

    // Process each line in the file
    for line in reader.lines() {
        let line = line.expect("Could not read line from SVS file");

        // Split the line by either tab or comma
        let fields: Vec<&str> = line.split(|c| c == '\t' || c == ',').collect();

        // Ensure the line has the expected number of fields
        if fields.len() < 6 {
            panic!("Invalid record format: expected 6 fields, found {}", fields.len());
        }

        // Extract additional components if needed (redundant here as they're directly loaded)
        let assembly_id = fields[0].to_string();
        let sequence_id = fields[1].to_string();
        let orient = fields.get(4).and_then(|field| field.chars().next());
        let start = fields[2].parse::<u64>().ok(); 
        let end = fields[3].parse::<u64>().ok();
        let sequence = fields[5].as_bytes().to_vec();
        // Generate an original ID (reconstruct from parts)
        let original_id = format!("{}:{}:{}-{}:{}", assembly_id, sequence_id, start.unwrap_or(0), end.unwrap_or(0), orient.unwrap_or('+'));

        // Store the parsed ID information along with the sequence and inferred version
        parsed_records.push(SequenceRecord {
            original_id: Some(original_id),
            assembly_id: Some(assembly_id),
            sequence_id,
            start,
            end,
            orient,
            inferred_version: None,
            sequence,
            aligned_seq: None,
            validated: None,
        });
    }

    parsed_records
}

fn boyer_moore_search_with_validation(
    records: &mut [SequenceRecord],
    genome_map: &HashMap<String, Vec<u8>>,
    debug_mode: bool,
) {
    records.par_iter_mut().filter(|r| r.validated.is_none()).for_each(|record| {
        let pattern = &record.sequence;
        let rev_complement_pattern = reverse_complement(pattern);
        let mut found_positions = vec![];
        let mut found_sequence_id = None;

        // Attempt search on the sequence specified by sequence_id
        if let Some(target_sequence) = genome_map.get(&record.sequence_id) {
            // Search forward strand
            found_positions.extend(boyer_moore_search(target_sequence, pattern).into_iter().map(|pos| (pos, '+')));

            // Search reverse strand
            found_positions.extend(boyer_moore_search(target_sequence, &rev_complement_pattern).into_iter().map(|pos| (pos, '-')));
        }

        // Check if there are matches on the specified sequence_id
        if found_positions.is_empty() {
            // Search other sequences in genome_map if no match is found on the specified sequence_id
            for (seq_name, genome_sequence) in genome_map {
                if seq_name == &record.sequence_id {
                    continue; // Skip the original sequence_id since it was already checked
                }
                // Search forward strand
                found_positions.extend(boyer_moore_search(genome_sequence, pattern).into_iter().map(|pos| (pos, '+')));

                // Search reverse strand
                found_positions.extend(boyer_moore_search(genome_sequence, &rev_complement_pattern).into_iter().map(|pos| (pos, '-')));

                if !found_positions.is_empty() {
                    found_sequence_id = Some(seq_name.clone());
                    break;
                }
            }
        } else {
            found_sequence_id = Some(record.sequence_id.clone());
        }

        // Process found positions
        if !found_positions.is_empty() {
            let closest = if found_positions.len() == 1 || record.start.is_none() {
                found_positions[0]
            } else  {
                // Find closest match to the start position in found_positions
                found_positions
                    .iter()
                    .min_by_key(|&&(pos, _)| (pos as isize - record.start.unwrap() as isize).abs())
                    .copied()
                    .unwrap_or(found_positions[0])
            };

            // Update start, end, orientation, sequence_id, and validated status
            record.start = Some(closest.0 as u64 + 1);
            record.end = Some((closest.0 + pattern.len()) as u64);
            record.orient = Some(closest.1);
            record.sequence_id = found_sequence_id.unwrap();
            record.validated = Some(if found_positions.len() == 1 {
                "fixed_remapped_unique".to_string()
            } else {
                "fixed_remapped_ambig".to_string()
            });
        }

        if debug_mode {
            if let Some(validation) = &record.validated {
                println!("Boyer-Moore {} fix for record: {:?}", validation, record);
            } else {
                println!("Boyer-Moore failed to fix for record: {:?}", record);
            }
        }
    });
}


fn bad_char_heuristic(pattern: &[u8], bad_char: &mut [isize; 256]) {
    for i in 0..256 {
        bad_char[i] = -1; // Initialize all occurrences as -1
    }

    for (i, &ch) in pattern.iter().enumerate() {
        bad_char[ch as usize] = i as isize; // Last occurrence of each character in the pattern
    }
}

fn boyer_moore_search(text: &[u8], pattern: &[u8]) -> Vec<usize> {
    let m = pattern.len();
    let n = text.len();
    if m == 0 || m > n {
        return vec![]; // Early return if pattern is empty or longer than text
    }

    let mut bad_char = [-1; 256];
    bad_char_heuristic(pattern, &mut bad_char);

    let mut positions = Vec::new();
    let mut s = 0; // s is the shift of the pattern with respect to text

    while s <= n - m {
        let mut j = (m - 1) as isize;

        // Decrease j while characters of pattern and text are matching at this shift s
        while j >= 0 && pattern[j as usize] == text[s + j as usize] {
            j -= 1;
        }

        if j < 0 {
            // Pattern found at current shift
            positions.push(s);

            // Shift pattern to align next character in text with last occurrence in pattern
            s += if s + m < n { m - bad_char[text[s + m] as usize] as usize } else { 1 };
        } else {
            // Shift pattern to align bad character in text with last occurrence in pattern
            s += (j - bad_char[text[s + j as usize] as usize]).max(1) as usize;
        }
    }

    positions
}



// LogLevel options
fn output_results(records: &[SequenceRecord], format: LogLevel, label: String) {
    match format {
        LogLevel::Summary | LogLevel::PerRecord => {
            let total_records = records.len();
            let mut fix_counts = HashMap::new();
            let mut fixed_count = 0;
            for record in records.iter() {
                if record.validated.is_some() && record.validated.as_deref() != Some("valid") &&
                   record.validated.as_deref() != Some("invalid") {
                    fixed_count = fixed_count + 1;
                    *fix_counts.entry(record.validated.clone().unwrap()).or_insert(0) += 1;
                }
            }
            let valid_count = records.iter().filter(|r| r.validated.as_deref() == Some("valid")).count();
            let invalid_count = records.iter().filter(|r| r.validated.as_deref() == Some("invalid")).count();

            println!("{}:", label);
            println!("  Total Sequences: {}", total_records);
            println!("     Accurate Coordinates: {}", valid_count);
            println!("     Repaired Coordinates: {}", fixed_count);
            for (fix_type, count) in fix_counts {
                println!("        {}: {}", fix_type, count);
            }
            println!("     Invalid Coordinates: {}", invalid_count);
        }
        LogLevel::Detailed => {
            println!("Detailed Report:");
            for record in records {
                record.print_record();
            }
        }
    }
}

fn main() {
    let args = Args::parse();
    let debug_mode = false;

    println!("##\n## DisCoord Version {}\n##", env!("CARGO_PKG_VERSION"));

    let output_file = &args.output;

    if let Some(n) = args.threads {
        ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .expect("Failed to build global thread pool");
        println!("## Threads: {}", n);
    } else {
        println!("## Threads: all available");
    }

    let mut start = Instant::now();
    let (is_gzip, format) = detect_format_and_compression(&args.input).expect("Failed to detect input format and compression");
    println!("## File: {}, Format: {}, Compression: {}", &args.input, format, if is_gzip { "Gzip" } else { "None" });

    let (ref_is_gzip, ref_format) = detect_format_and_compression(&args.reference).expect("Failed to detect reference format and compression");
    println!("## Reference: {}, format: {}, Compression: {}", &args.reference, ref_format, if ref_is_gzip { "Gzip" } else { "None" });
    let genome_map = if ref_format == "TwoBit" {
        load_genome_from_2bit_parallel(&args.reference).expect("Failed to load genome")
    } else {
        // TODO: Support compressed fasta files
        load_genome_from_fasta_parallel(&args.reference).expect("Failed to load genome")
    };

    let loading_duration = start.elapsed(); // Get the elapsed time

    // TODO: we should probably flag duplicate coordinates following fixes
   
    start = Instant::now();
    if  format == "Fasta" {
        let mut records = load_and_parse_fasta(&args.input);
        validate_sequences(&mut records, &genome_map, false);
        // Run Boyer-Moore search on unvalidated records
        if args.map_sequences {
            boyer_moore_search_with_validation(&mut records, &genome_map, debug_mode);
        }
        for record in records.iter_mut() {
            if record.validated.is_none() {
                record.validated = Some("invalid".to_string());
            }
        }
        output_results(&records, args.log_level, format!("Summary for {}",args.input));
        if let Some(output_file) = &args.output {
            write_fasta_output(&records, output_file, is_gzip, false).expect("Failed to write output");
        }
    } else if format == "Stockholm" {
        parse_and_validate_stockholm(&args.input, &genome_map, &args.output, is_gzip, args.log_level, args.map_sequences).expect("Failed to parse Stockholm file");
    } else if format == "TabDelimited" || format == "CommaDelimited" {
        let mut records = load_and_parse_delimited_file(&args.input);
        validate_sequences(&mut records, &genome_map, false);
        // Run Boyer-Moore search on unvalidated records
        if args.map_sequences {
            boyer_moore_search_with_validation(&mut records, &genome_map, debug_mode);
        }
        for record in records.iter_mut() {
            if record.validated.is_none() {
                record.validated = Some("invalid".to_string());
            }
        }
        output_results(&records, args.log_level, format!("Summary for {}",args.input));
        if let Some(output_file) = &args.output {
            write_delimited_output(&records, output_file, is_gzip, false, format).expect("Failed to write output");
        }
    } else {
        panic!("Unsupported format");
    }
    let validation_duration = start.elapsed(); // Get the elapsed time

    println!("\nRuntime stats:");
    println!("  Loading sequences: {:?} s", loading_duration);
    println!("         Validation: {:?} s", validation_duration);
}


pub fn load_genome_from_fasta_parallel(path: &str) -> io::Result<HashMap<String, Vec<u8>>> {
    // Open the FASTA file for reading
    let reader = fasta::Reader::from_file(path).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    
    // Collect sequences into a vector of tuples (name, sequence)
    let sequences: Vec<(String, Vec<u8>)> = reader
        .records()
        .par_bridge() // Convert to a parallel iterator
        .map(|result| {
            let record = result.map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            let name = record.id().to_string(); // Sequence ID
            let sequence = record.seq().to_ascii_uppercase().to_vec(); // Sequence data
            Ok((name, sequence))
        })
        .collect::<Result<_, io::Error>>()?;
    
    // Collect into a HashMap
    let genome_map: HashMap<String, Vec<u8>> = sequences.into_iter().collect();

    Ok(genome_map)
}



pub fn load_genome_from_2bit_parallel(path: &str) -> io::Result<HashMap<String, Vec<u8>>> {
    let file = File::open(path)?;
    let mmap = unsafe { Mmap::map(&file)? };
    let buffer = &mmap[..];

    // Step 1: Parse header
    let is_little_endian = match u32::from_be_bytes(buffer[0..4].try_into().unwrap()) {
        0x1A412743 => false,
        0x4327411A => true,
        _ => return Err(io::Error::new(io::ErrorKind::InvalidData, "Invalid 2bit signature")),
    };

    let read_u32 = |offset: usize| {
        let bytes: [u8; 4] = buffer[offset..offset + 4].try_into().unwrap();
        if is_little_endian {
            u32::from_le_bytes(bytes)
        } else {
            u32::from_be_bytes(bytes)
        }
    };

    let version = read_u32(4);
    if version > 1 {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "Unsupported 2bit version"));
    }

    let seq_count = read_u32(8) as usize;

    // Step 2: Parse sequence metadata
    let mut sequences = Vec::new();
    let mut offset = 16;

    for _ in 0..seq_count {
        let name_len = buffer[offset] as usize;
        offset += 1;

        let name = String::from_utf8(buffer[offset..offset + name_len].to_vec()).unwrap();
        offset += name_len;

        let seq_offset = if version == 0 {
            read_u32(offset) as u64
        } else {
            u64::from_be_bytes(buffer[offset..offset + 8].try_into().unwrap())
        };
        offset += if version == 0 { 4 } else { 8 };

        sequences.push((name, seq_offset));
    }

    // Step 3: Decode DNA in parallel
    let genome_map: HashMap<String, Vec<u8>> = sequences
        .into_par_iter()
        .map(|(name, seq_offset)| {
            let dna_size = read_u32(seq_offset as usize) as usize;

            // Read N block data
            let n_block_count = read_u32((seq_offset + 4) as usize) as usize;
            let mut n_block_starts = Vec::with_capacity(n_block_count);
            let mut n_block_sizes = Vec::with_capacity(n_block_count);

            let mut current_offset = (seq_offset + 8) as usize;

            for _ in 0..n_block_count {
                let start = read_u32(current_offset) as usize;
                n_block_starts.push(start);
                current_offset += 4;
            }

            for _ in 0..n_block_count {
                let size = read_u32(current_offset) as usize;
                n_block_sizes.push(size);
                current_offset += 4;
            }
            
            // Read mask block data
            let mask_block_count = read_u32((current_offset) as usize) as usize;
            // For now we do not need the mask data
            current_offset = current_offset + (mask_block_count * 8) + 4;

            // Skipped the reserved bytes
            current_offset += 4;

            // Read packed DNA
            let mut genome = vec![b'N'; dna_size];
            for i in 0..((dna_size + 3) / 4) {
                let byte = buffer[current_offset + i];
                for j in 0..4 {
                    let pos = i * 4 + j;
                    if pos >= dna_size {
                        break;
                    }
                    genome[pos] = match (byte >> ((3 - j) * 2)) & 0b11 {
                        0 => b'T',
                        1 => b'C',
                        2 => b'A',
                        3 => b'G',
                        _ => b'N',
                    };
                }
            }

            // Apply N blocks
            for (&start, &size) in n_block_starts.iter().zip(n_block_sizes.iter()) {
                for pos in start..(start + size) {
                    if pos < genome.len() {
                        genome[pos] = b'N';
                    }
                }
            }

            // Return the sequence name and genome data
            (name, genome)
        })
        .collect();

    Ok(genome_map)
}



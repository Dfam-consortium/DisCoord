use clap::{Parser, ValueEnum};
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::io::{BufReader, BufRead};
use std::fs::OpenOptions;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::io::BufWriter;
use std::io::{self, Read, Seek, SeekFrom};
use std::path::Path;
use std::time::Instant;
use smitten::{Identifier, IDVersion};
use bio::io::fasta;

/// Command-line arguments
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Path to the 2bit genome file
    #[arg(short, long)]
    twobit: String,

    /// Enable mapping (Boyer-Moore) for invalid identifiers
    #[arg(short, long, default_value = "false")]
    map_sequences: bool,

    /// Optional output file path
    #[arg(short, long, default_value = "")]
    output: String,

    /// Log level (Summary or Detailed)
    #[arg(short, long, default_value = "summary")]
    log_level: LogLevel,

    /// Input file path (Fasta or Stockholm format)
    #[arg()]
    input: String,
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
            self.validated.as_deref().unwrap_or("Not validated"),
        );
    }
}

#[derive(Clone, Debug, ValueEnum)]
enum LogLevel {
    Summary,
    Detailed,
}


fn parse_and_validate_stockholm(file_path: &str, genome_map: &HashMap<String, Vec<u8>>, output_path: &str, is_gzip: bool, log_level: LogLevel, map_seqs: bool) -> io::Result<()> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    let mut current_record: HashMap<String, String> = HashMap::new();
    let mut metadata: Vec<String> = Vec::new();

    for line in reader.lines() {
        let line = line?;

        if line.is_empty() {
            continue;
        }

        if line.starts_with("#") {
            metadata.push(line);
            continue;
        }

        if line.starts_with("//") {
            process_stockholm_record(&mut current_record, &metadata, genome_map, output_path, is_gzip, &log_level, map_seqs);
            current_record.clear();
        } else {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() == 2 {
                let (name, seq) = (parts[0].to_string(), parts[1].to_string());
                current_record.entry(name).or_insert(String::new()).push_str(&seq);
            }
        }
    }

    // Process any remaining record
    if !current_record.is_empty() {
        process_stockholm_record(
            &mut current_record,
            &metadata,
            genome_map,
            output_path,
            is_gzip,
            &log_level,
            map_seqs,
        );
    }

    Ok(())
}

    
fn process_stockholm_record(
    current_record: &mut HashMap<String, String>,
    metadata: &[String],
    genome_map: &HashMap<String, Vec<u8>>,
    output_path: &str,
    is_gzip: bool,
    log_level: &LogLevel,
    map_seqs: bool,
) {
    let mut records = Vec::new();

    for (name, seq) in current_record.drain() {
        let (v2_id, inferred_version) = Identifier::from_unknown_format(&name, false)
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

    // Log what we did
    output_results(&records, log_level.clone());

    // Write the output in append mode

    if ! output_path.is_empty() {
        write_stockholm_output(&records, metadata, output_path, is_gzip, true).expect("Failed to write Stockholm output");
    }
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


/// Detects format and compression type for input files (FASTA or Stockholm).
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

    let mut format = None;
    for _ in 0..15 {
        let mut line = String::new();
        if reader.read_line(&mut line)? == 0 {
            break; // EOF reached
        }

        // Skip blank or whitespace-only lines
        if line.trim().is_empty() {
            continue;
        }

        // Check for FASTA or Stockholm format based on line content
        if line.starts_with(">") {
            format = Some("FASTA");
            break;
        } else if line.starts_with("# STOCKHOLM 1.0") {
            format = Some("Stockholm");
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
        //println!("  Working on fasta_sequence {:?}", fasta_sequence);
        //println!("  Working on genome_sequence {:?}", genome_sequence[start..end].to_vec());

        // 1. Direct match on same strand
        if start < genome_sequence.len() && end <= genome_sequence.len() {
            if &genome_sequence[start..end] == fasta_sequence {
                located = true;
                if record.orient != Some('+') {
                    validation_str.push_str("_orient");
                    record.orient = Some('-');
                }
                if debug_mode {
                    println!("Direct match validated for: {:?}", record);
                }
            }
        }

        // 2. Reverse complement match on opposite strand
        if !located && start < genome_sequence.len() && end <= genome_sequence.len() {
            if &genome_sequence[start..end] == &rev_complement {
                located = true;
                if record.orient == Some('+') {
                    validation_str.push_str("_orient");
                    record.orient = Some('-');
                }
                if debug_mode {
                    println!("Reverse complement validated for: {:?}", record);
                }
            }
        }

        if !located {
            // 3. Shifted matches on both strands within +/-3bp range
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
        let (v2_id, inferred_version) = Identifier::from_unknown_format(&original_id, false).expect("Failed to convert ID to V2");
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
            sequence: rec.seq().to_owned(),
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
    //for record in records.iter_mut().filter(|r| r.validated.is_none()) {
        let pattern = &record.sequence;
        let rev_complement_pattern = reverse_complement(pattern);
        let mut found_positions = vec![];
        let mut found_sequence_id = None;
        //let mut found_orientation: Option<char> = None;

        //println!("Searching for pattern: {:?}", pattern);

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
            let closest = if found_positions.len() == 1 {
                found_positions[0]
            } else {
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
fn output_results(records: &[SequenceRecord], format: LogLevel) {
    match format {
        LogLevel::Summary => {
            let total_records = records.len();
            //let fixed_records: Vec<_> = records.iter().filter(|r| r.validated.is_some()).collect();
            let mut fix_counts = HashMap::new();
            let mut fixed_count = 0;
            for record in records.iter() {
                if record.validated.is_some() && record.validated.as_deref() != Some("valid") {
                    fixed_count = fixed_count + 1;
                    *fix_counts.entry(record.validated.clone().unwrap()).or_insert(0) += 1;
                }
            }
            let valid_count = records.iter().filter(|r| r.validated.as_deref() == Some("valid")).count();
            let unvalidated_count = records.iter().filter(|r| r.validated.is_none()).count();

            println!("Summary Report:");
            println!("  Total Sequences: {}", total_records);
            println!("     Accurate Coordinates: {}", valid_count);
            println!("     Repaired Coordinates: {}", fixed_count);
            for (fix_type, count) in fix_counts {
                println!("        {}: {}", fix_type, count);
            }
            println!("     Invalid Coordinates: {}", unvalidated_count);
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

    // Check if the output file already exists
    if !args.output.is_empty() && Path::new(&args.output).exists() {
        eprintln!("Error: The output file '{}' already exists. Please specify a different filename.", &args.output);
        std::process::exit(1); // Exit the program to avoid overwriting
    }

    let mut start = Instant::now();
    let (is_gzip, format) = detect_format_and_compression(&args.input).expect("Failed to detect format and compression");
    let genome_map = load_genome_from_2bit_parallel(&args.twobit).expect("Failed to load genome");
    let loading_duration = start.elapsed(); // Get the elapsed time

    // TODO: we should probably flag duplicate coordinates following fixes
   
    start = Instant::now();
    if  format == "FASTA" {
        let mut records = load_and_parse_fasta(&args.input);
        validate_sequences(&mut records, &genome_map, false);
        // Run Boyer-Moore search on unvalidated records
        if args.map_sequences {
            boyer_moore_search_with_validation(&mut records, &genome_map, debug_mode);
        }
        output_results(&records, args.log_level);
        if ! args.output.is_empty() {
            write_fasta_output(&records, &args.output, is_gzip, false).expect("Failed to write output");
        }
    } else if format == "Stockholm" {
        parse_and_validate_stockholm(&args.input, &genome_map, &args.output, is_gzip, args.log_level, args.map_sequences).expect("Failed to parse Stockholm file");
    } else {
        panic!("Unsupported format");
    }
    let validation_duration = start.elapsed(); // Get the elapsed time

    println!("\nRuntime stats:");
    println!("  Loading sequences: {:?} s", loading_duration);
    println!("         Validation: {:?} s", validation_duration);
}



fn load_genome_from_2bit_parallel(path: &str) -> io::Result<HashMap<String, Vec<u8>>> {
    let mut file = File::open(Path::new(path))?;
    let mut buffer = [0; 4];

    // Step 1: Read and check the file signature to determine endianness
    file.read_exact(&mut buffer)?;
    let signature_be = u32::from_be_bytes(buffer);
    let signature_le = u32::from_le_bytes(buffer);
    let is_little_endian = if signature_be == 0x1A412743 {
        false
    } else if signature_le == 0x1A412743 {
        true
    } else {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "Invalid 2bit file signature"));
    };

    // Step 2: Read the version number (4 bytes)
    file.read_exact(&mut buffer)?;
    let version = if is_little_endian {
        u32::from_le_bytes(buffer)
    } else {
        u32::from_be_bytes(buffer)
    };
    if version > 1 {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "Unsupported 2bit file version"));
    }

    // Step 3: Read the sequence count
    file.read_exact(&mut buffer)?;
    let seq_count = if is_little_endian {
        u32::from_le_bytes(buffer)
    } else {
        u32::from_be_bytes(buffer)
    } as usize;

    // Step 4: Skip the reserved 4 bytes
    file.seek(SeekFrom::Current(4))?;

    // Step 5: Read sequence names and offsets without seeking to sequence data yet
    let mut sequences = Vec::new();
    for _ in 0..seq_count {
        // Read the sequence name length (1 byte)
        let mut name_len_buf = [0; 1];
        file.read_exact(&mut name_len_buf)?;
        let name_len = name_len_buf[0] as usize;

        // Read the sequence name
        let mut name_buf = vec![0; name_len];
        file.read_exact(&mut name_buf)?;
        let name = String::from_utf8(name_buf).expect("Invalid UTF-8 in sequence name");

        // Read the sequence offset based on file version
        let offset = if version == 0 {
            file.read_exact(&mut buffer)?;
            if is_little_endian {
                u32::from_le_bytes(buffer) as u64
            } else {
                u32::from_be_bytes(buffer) as u64
            }
        } else {
            let mut offset_buf = [0; 8];
            file.read_exact(&mut offset_buf)?;
            if is_little_endian {
                u64::from_le_bytes(offset_buf)
            } else {
                u64::from_be_bytes(offset_buf)
            }
        };

        // Collect sequence metadata
        sequences.push((name, offset));
    }

    // Step 6: Load each sequence in parallel
    let genome_map: HashMap<String, Vec<u8>> = sequences
        .into_par_iter()
        .map(|(name, offset)| {
            let mut file = File::open(Path::new(path)).expect("Could not reopen 2bit file");
            let mut buffer = [0; 4];

            // Seek to the sequence offset
            file.seek(SeekFrom::Start(offset)).expect("Seek failed");

            // Step 6.1: Read dnaSize
            file.read_exact(&mut buffer).expect("Failed to read dnaSize");
            let dna_size = if is_little_endian {
                u32::from_le_bytes(buffer) as usize
            } else {
                u32::from_be_bytes(buffer) as usize
            };

            // Step 6.2: Read nBlockCount and nBlockStarts/nBlockSizes
            file.read_exact(&mut buffer).expect("Failed to read nBlockCount");
            let n_block_count = if is_little_endian {
                u32::from_le_bytes(buffer) as usize
            } else {
                u32::from_be_bytes(buffer) as usize
            };

            let mut n_block_starts = vec![0; n_block_count];
            for start in n_block_starts.iter_mut() {
                file.read_exact(&mut buffer).expect("Failed to read nBlockStart");
                *start = if is_little_endian {
                    u32::from_le_bytes(buffer) as usize
                } else {
                    u32::from_be_bytes(buffer) as usize
                };
            }

            let mut n_block_sizes = vec![0; n_block_count];
            for size in n_block_sizes.iter_mut() {
                file.read_exact(&mut buffer).expect("Failed to read nBlockSize");
                *size = if is_little_endian {
                    u32::from_le_bytes(buffer) as usize
                } else {
                    u32::from_be_bytes(buffer) as usize
                };
            }

            // Step 6.3: Read maskBlockCount and maskBlockStarts/maskBlockSizes
            file.read_exact(&mut buffer).expect("Failed to read maskBlockCount");
            let mask_block_count = if is_little_endian {
                u32::from_le_bytes(buffer) as usize
            } else {
                u32::from_be_bytes(buffer) as usize
            };

            let mut mask_block_starts = vec![0; mask_block_count];
            for start in mask_block_starts.iter_mut() {
                file.read_exact(&mut buffer).expect("Failed to read maskBlockStart");
                *start = if is_little_endian {
                    u32::from_le_bytes(buffer) as usize
                } else {
                    u32::from_be_bytes(buffer) as usize
                };
            }

            let mut mask_block_sizes = vec![0; mask_block_count];
            for size in mask_block_sizes.iter_mut() {
                file.read_exact(&mut buffer).expect("Failed to read maskBlockSize");
                *size = if is_little_endian {
                    u32::from_le_bytes(buffer) as usize
                } else {
                    u32::from_be_bytes(buffer) as usize
                };
            }

            // Step 6.4: Skip the reserved 4 bytes
            file.seek(SeekFrom::Current(4)).expect("Failed to skip reserved bytes");

            // Step 6.5: Read packed DNA
            let mut genome = vec![b'N'; dna_size];
            for i in 0..((dna_size + 3) / 4) {
                file.read_exact(&mut buffer[0..1]).expect("Failed to read packed DNA");
                let byte = buffer[0];
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

            // Apply mask blocks (convert to lowercase for masked regions)
            for (&start, &size) in mask_block_starts.iter().zip(mask_block_sizes.iter()) {
                for pos in start..(start + size) {
                    if pos < genome.len() {
                        genome[pos] = genome[pos].to_ascii_lowercase();
                    }
                }
            }

            // Return the sequence name and data
            (name, genome)
        })
        .collect();

    Ok(genome_map)
}


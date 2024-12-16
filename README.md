
# DisCoord : A tool for validating and fixing sequence coordinates in FASTA or Stockholm files

## Overview

DisCoord is a tool for validating and fixing sequence coordinates in FASTA or Stockholm files. 
In addition to detecting small errors in the coordinates, such as shifted ranges, mixed up
interval boundaries ( e.g half-open, fully closed), or strand errors, DisCoord can also map 
the sequences to the reference sequence to detect and fix larger errors.

DisCoord is designed to work with sequence identifiers that include coordinates and strand
information.  Many formats do not formally define these data elements (I am looking at you
FASTA), and as a consequence many ad-hoc encodings have flourished. We have our own custom
encoding (Smitten Format -- https://github.com/Dfam-consortium/Smitten), but any parsable 
format could be added to the code to support coordinate validations.

For example a FASTA record might look like this:

```
>hs1:chr1:100-200_+
taaccctaaccctaaccctaaccctaaccctaaccctaaccctaacccta
accctaaccctaaccctaaccctaacccctaaccctaaccctaaccctaa
c
```

Where this sequence identifier encodes the assembly ('hs1'), the sequence identifier ('chr1'),
the start and end coordinates (100-200 a one-based, fully-closed interval) and the strand ('+'). 

## Installation

* [Install Rust](https://www.rust-lang.org/tools/install)
* Run **`cargo build --release`** to build the release binary file.
* The binary 'discoord' will be in the **`target/release`** directory.


## Usage

```
% ./target/release/discoord -h

Command-line arguments

Usage: discoord [OPTIONS] --reference <REFERENCE> --output <OUTPUT> <INPUT>

Arguments:
  <INPUT>  Input file path (Fasta, Stockholm or Tab/Comma Delimited file)

Options:
  -r, --reference <REFERENCE>  Path to the reference sequence file [twobit, or fasta]
  -m, --map-sequences          Enable mapping (Boyer-Moore) for invalid identifiers
  -o, --output <OUTPUT>        Optional output file path
  -l, --log-level <LOG_LEVEL>  Log level (Summary, PerRecord or Detailed) [default: summary] [possible values: summary, per-record, detailed]
  -x, --threads <THREADS>      Threads to use for parallel processing
  -h, --help                   Print help
  -V, --version                Print version


# Validate a FASTA file, looking for minor errors in coordinates
% ./target/release/discoord -r hg38.2bit test.fasta

or using a FASTA reference rather than UCSC 2bit format:

% ./target/release/discoord -r hg38.fa test.fasta 

# Validate a Stockholm file, looking for minor errors in coordinates
% ./target/release/discoord -r hg38.2bit test.stk

# Identify more serious errors using the Boyer-Moore algorithm to map sequences to the reference
% ./target/release/discoord -r hg38.2bit -m test.fasta

# Identify and fix errors, generating a new sequence file
% ./target/release/discoord -r hg38.2bit -m -o fixed.fasta test.fasta

```

### Delimited Format

DisCoord can read tab or comma delimited files with the following format:

```
assembly_id: An assembly/reference identifier [or blank]
sequence_id: A unique sequence identifier in the reference
      start: The start coordinate (one-based) 
        end: The end coordinate (one-based-fully-closed)
   sequence: The sequence data
```

### Smitten Format

The Smitten format is a simple encoding for sequence identifiers that (as of V2) includes:

* An assembly identifier [optional]
* A sequence identifier
* A start coordinate (one-based) [optional
* An end coordinate  (one-based-fully-closed) [optional]
* A strand indicator [optional]

The field are seperated as follows: assembly_id:sequene_id:start-end_strand 

For example 'hs1:chr1:100-200_+' is a valid Smitten V2 identifier. Multiple subranges
are possible, for example: 'hs1:chr1:100-200_+:5-10_-' is also a valid identifier
and indicates that the sequence is base pairs 5-10 of the sequence parent sequence 
'hs1:chr1:100-200_+' and is on the negative strand.

The only reserved character is the colon (':') which is used to separate the assembly_id,
sequence_id and subranges. 

### TODO
- Add support for file globbing for a set of input files against a reference
- Setup a quiet mode to return overall status with a return code

# BioPack-FQ

**Base-4 Genomic Encoding Format**

BioPack-FQ is a high-performance, open-source binary encoding format designed to replace FASTA and FASTQ files for modern sequencing workloads. DNA is fundamentally a base-4 (quaternary) information system — each nucleotide carries exactly 2 bits of information. By encoding nucleotides at 2 bits per base and packing quality scores at 6 bits per value, BioPack-FQ maps genomic data to its natural information-theoretic density while preserving lossless roundtrip fidelity.

---

## Table of Contents

- [Abstract](#abstract)
- [1. Introduction](#1-introduction)
- [2. Methods](#2-methods)
  - [2.1 Nucleotide Encoding](#21-nucleotide-encoding)
  - [2.2 Ambiguous Base (N) Handling](#22-ambiguous-base-n-handling)
  - [2.3 Quality Score Encoding](#23-quality-score-encoding)
  - [2.4 File Architecture](#24-file-architecture)
  - [2.5 Metadata Header](#25-metadata-header)
  - [2.6 Payload Entropy Coding](#26-payload-entropy-coding)
  - [2.7 Multi-Threaded Encoding](#27-multi-threaded-encoding)
- [3. Results](#3-results)
  - [3.1 Encoding Efficiency](#31-encoding-efficiency)
  - [3.2 Lossless Roundtrip Validation](#32-lossless-roundtrip-validation)
  - [3.3 Pipeline Impact Assessment](#33-pipeline-impact-assessment)
- [4. Usage](#4-usage)
- [5. Build](#5-build)
- [6. Technology Stack](#6-technology-stack)
- [7. Future Work](#7-future-work)
- [8. License](#8-license)

---

## Abstract

Next-generation sequencing produces vast quantities of data, commonly stored in FASTA or FASTQ formats. These human-readable formats use 8-bit ASCII characters to represent a fundamentally 2-bit-per-symbol information source. BioPack-FQ addresses this by introducing a compact binary format (`.g` for sequences, `.gq` for sequences with quality scores) that performs **base-4 encoding** of DNA: each of the four nucleotides (A, C, G, T) is mapped to a 2-bit code, achieving the information-theoretic lower bound. Phred quality scores are packed at 6 bits per value (the range of biologically meaningful scores). Ambiguous bases (N) are stored via a position index rather than in-band. On a 16 MB Oxford Nanopore dataset comprising 29,170 reads and 5.99 million bases, BioPack-FQ produces a 6.1 MB `.gq` file from base-4 encoding alone (**62.9% reduction**), and 5.1 MB with optional ZSTD entropy coding on top (**69.2% reduction**).

---

## 1. Introduction

FASTA and FASTQ have been the cornerstone of genomic data storage for decades. Despite their ubiquity, their ASCII-based representation is information-theoretically wasteful:

- **Nucleotides** stored as 8-bit ASCII characters (A, C, G, T) when 2 bits suffice (log₂4 = 2)
- **Quality scores** stored as 8-bit ASCII characters (Phred+33 or Phred+64) when 6 bits cover the biologically relevant range (0–60)
- **Read headers** duplicated across millions of entries

The fundamental insight: DNA is a base-4 information channel. Each position carries exactly one of four symbols — 2 bits of information. ASCII encoding inflates this by 4×. BioPack-FQ restores the encoding to its natural density.

---

## 2. Methods

### 2.1 Nucleotide Encoding

Each of the four standard DNA bases is assigned a 2-bit code — the minimum possible for a 4-symbol alphabet:

| Base | Binary |
|------|--------|
| A    | 00     |
| C    | 01     |
| G    | 10     |
| T    | 11     |

Four bases are packed into a single byte:

```
Byte: [b0:2][b1:2][b2:2][b3:2]
       MSB                  LSB
```

This achieves a 4:1 encoding density improvement over ASCII FASTA.

### 2.2 Ambiguous Base (N) Handling

Ambiguous bases (N) fall outside the 2-bit code space. Rather than sacrificing a code point or using a bitmask, their global positions are recorded in a sorted **N-position list** (none_list) stored after the file header. During encoding, N bases are omitted from the 2-bit packed stream and their positions appended to the list. During decoding, the packed stream and N-position list are merged to reconstruct the original sequence.

This approach adds no overhead for datasets without Ns and only ~0.25 bits per N for datasets with Ns.

### 2.3 Quality Score Encoding

Quality scores (Phred) are auto-detected as Phred+33 or Phred+64 by scanning the minimum ASCII value across all reads.

Each raw quality score (0–63) is stored in **6 bits** — sufficient to cover the full Phred range used in modern sequencing. Four scores are packed into three bytes:

```
Byte 0: [s0:6][s1:2]
Byte 1: [s1:4][s2:4]
Byte 2: [s2:2][s3:6]
```

This achieves a 25% density improvement over the 8-bit ASCII representation.

Header metadata records `phred_offset`, `phred_min`, `phred_max`, `phred_avg`, `q20_count`, and `q30_count` for rapid quality assessment without scanning the packed data.

### 2.4 File Architecture

```
┌─────────────────────────────────┐
│      GQHeader (200 bytes)       │  Fixed metadata header
├─────────────────────────────────┤
│      Read Lengths (variable)     │  uint32_t per read
├─────────────────────────────────┤
│      N-Position List            │  uint32_t count + positions
├─────────────────────────────────┤
│     2-bit Packed Sequence Data   │  Base-4 encoded A/C/G/T
├─────────────────────────────────┤
│  6-bit Packed Quality (.gq)      │  Phred scores (optional)
├─────────────────────────────────┤
│  (Optional) ZSTD Entropy Coding │  On-disk payload compression
└─────────────────────────────────┘
```

The first five sections form the semantic payload (base-4 sequences + 6-bit qualities). When written to disk, the payload may optionally be entropy-coded with ZSTD. The header's `compression` field (0 or 1) records whether entropy coding was applied.

### 2.5 Metadata Header

The 200-byte fixed header stores comprehensive metadata:

| Field | Size | Description |
|-------|------|-------------|
| `version` | 4 B | Format version |
| `file_type` | 1 B | 0 = `.g`, 1 = `.gq` |
| `compression` | 1 B | 0 = none, 1 = ZSTD |
| `original_size` | 8 B | Uncompressed payload size |
| `compressed_size` | 8 B | Compressed payload size |
| `read_count` | 8 B | Number of reads |
| `nucleotide_count` | 8 B | Total bases (incl. Ns) |
| `min_read_length` | 4 B | Shortest read |
| `max_read_length` | 4 B | Longest read |
| `avg_read_length` | 4 B | Mean read length |
| `gc_content` | 8 B | GC percentage |
| `count_a/c/g/t/n` | 40 B | Nucleotide counts |
| `phred_offset/min/max/avg` | 4 B | Quality statistics |
| `q20_count/q30_count` | 16 B | Bases above quality thresholds |
| `checksum` | 16 B | MD5 binary checksum |
| `original_file_size` | 8 B | Original input file size |
| `created_timestamp` | 8 B | Unix timestamp |

### 2.6 Payload Entropy Coding

After the semantic base-4 encoding, the binary payload can be further entropy-coded with **ZSTD** using its advanced API (`ZSTD_CCtx` with `ZSTD_c_nbWorkers` set to the number of available CPU cores, compression level 3). This is strictly optional — the base-4 encoding is the format's core contribution; ZSTD is a standard entropy coder applied on top of the already-encoded representation.

### 2.7 Multi-Threaded Encoding

Encoding and decoding are parallelized at multiple levels:

- **Sequence encoding**: The concatenated sequence string is split into chunks; each thread independently 2-bit packs its chunk and records N positions, with results merged in order.
- **ZSTD entropy coding**: The ZSTD library's built-in multi-threading (`nbWorkers`) compresses the payload using all available cores.
- **Sequence decoding**: Reads are split into chunks; each thread independently decodes its chunk using binary search on the N-position list to determine its starting state.

---

## 3. Results

### 3.1 Encoding Efficiency

Test dataset: Oxford Nanopore *rbcLa* amplicon sequencing, 29,170 reads, 5,989,618 total bases, Phred+33 qualities.

| Format | Size | vs. FASTQ |
|--------|------|-----------|
| Original FASTQ | 16,450,568 B | — |
| BioPack-FQ base-4 encoding only (`.gq`) | 6,106,504 B | **62.9% smaller** |
| BioPack-FQ + ZSTD entropy coding (`.gq`) | 5,062,164 B | **69.2% smaller** |

The primary savings come from the base-4 encoding itself (75% reduction on sequence data vs. ASCII) and 6-bit quality packing (25% on quality data). ZSTD provides additional entropy reduction on the already-densified payload.

### 3.2 Lossless Roundtrip Validation

Full roundtrip testing (FASTQ → `.gq` → decoded FASTQ) produces byte-identical sequences and quality scores, verified by `diff` across all 29,170 reads and 5.99 million bases.

### 3.3 Pipeline Impact Assessment

BioPack-FQ is a **storage and transfer format**, not a compute accelerator. Core pipeline stages — alignment, sorting, variant calling — are CPU-bound and will not run faster because the input is stored in `.gq` instead of FASTQ. The format's value lies in the data logistics that surround computation.

#### What BioPack-FQ Does Not Do

- Does not accelerate alignment (BWA, minimap2, Bowtie)
- Does not accelerate variant calling (GATK, FreeBayes, DeepVariant)
- Does not provide random access to individual reads (future feature)
- Does not replace indexing (`.bai`, `.fai`)

#### Where BioPack-FQ Provides Real Benefit

| Domain | Benefit | Why |
|--------|---------|-----|
| **Cloud storage** | 69% less egress/ingress cost | Files are 69% smaller |
| **Cloud transfer** | 69% faster upload/download | Proportional to size reduction |
| **Multi-sample projects** | 100 samples: 15 TB → 4.6 TB | Storage cost savings at scale |
| **Archival (S3 Glacier, tape)** | Lower tier cost per genome | Pay for 48 GB instead of 150 GB |
| **HPC scratch space** | Reduced disk footprint per job | Important on shared filesystems |
| **Concurrent NAS access** | Less I/O contention per read | Multiple pipelines reading simultaneously |
| **Output of basecallers** | Encoded on-the-fly without gzip overhead | No re-encode cycle |

#### Projected I/O Scenarios

| Operation | Format | Data read | Relative I/O time |
|-----------|--------|-----------|-------------------|
| Load into memory | FASTQ | 150 GB | 1.0× |
| Load into memory | gzip FASTQ | 50 GB + decompress | 2–3× slower * |
| Load into memory | `.gq` | 46 GB | **0.31×** |
| FASTQ → aligner pipe | FASTQ | 150 GB read + write /dev/null | 1.0× |
| `.gq` → aligner pipe | `.gq` decode → stdout | 46 GB read + bit unpack | **~0.5×** |

* gzip decompression is CPU-bound and typically reads slower than raw I/O on modern SSDs.

#### Honest Assessment

For a single-sample variant-calling pipeline (GATK best practices, ~24 h total), switching from FASTQ to `.gq` would save **under 30 min** — most of the pipeline is CPU-bound alignment and calling. The real savings come when that sample is:

- Uploaded to the cloud for processing
- Stored for 6 months post-analysis
- One of 10,000 samples in a biobank
- Transferred between collaborators

#### Key Takeaway

BioPack-FQ optimizes **data logistics**, not computation. Use it to reduce storage costs, accelerate transfers, and ease I/O pressure in multi-sample environments. Do not expect faster alignment or variant calling.

---

## 4. Usage

```
biopack <input> [-d] [-v] [-o output]
```

| Input extension | Default action | Output |
|----------------|----------------|--------|
| `.fasta`, `.fa` | Encode | `.g` |
| `.fastq`, `.fq` | Encode | `.gq` |
| `.fasta.gz`, `.fa.gz` | Encode (auto-decompress) | `.g` |
| `.fastq.gz`, `.fq.gz` | Encode (auto-decompress) | `.gq` |
| `.g` | Show info | stdout |
| `.gq` | Show info | stdout |

| Flag | Action |
|------|--------|
| `-v` | View full content (FASTA/FASTQ to stdout) |
| `-d` | Decode to FASTA/FASTQ |
| `-o <path>` | Write output to file (all modes) |

**Examples:**

```bash
# Encode FASTA
./biopack sample.fasta

# Encode gzipped FASTQ
./biopack sample.fastq.gz

# View file info
./biopack sample.gq

# View full content
./biopack sample.gq -v

# Decode to FASTQ
./biopack sample.gq -d -o decoded.fastq

# Write info to file
./biopack sample.gq -o info.txt
```

---

## 5. Build

**Dependencies:** `g++` (C++20), `zlib` (development headers), `libzstd` (development headers)

```bash
git clone <repo>
cd fast_to_g
make
make test
```

The build produces a single `biopack` binary.

---

## 6. Technology Stack

| Component | Technology |
|-----------|------------|
| Language | C++20 |
| Build system | GNU Make |
| Base-4 encoding | 2-bit nucleotide packing |
| Quality encoding | 6-bit Phred score packing |
| Input decompression | zlib (gzread) |
| Payload entropy coding | ZSTD (multi-threaded) |
| Parallelism | std::thread / std::async |
| Hashing | Planned (OpenSSL / xxHash) |
| SIMD | Planned (AVX2 / AVX-512) |

---

## 7. Future Work

- **Per-read header storage** for full FASTQ roundtrip with read names
- **Random-access indexing** for querying specific reads without full decode
- **Streaming support** for pipelining with bioinformatics tools
- **k-mer frequency tables** and **Bloom filters** for rapid dataset comparison
- **SIMD-optimized** packing and unpacking kernels
- **Checksum validation** via MD5 / SHA-256
- **GZIP** as an additional entropy coder option

---

## 8. License

GNU General Public License v3.0 (GPLv3). Commercial licensing exceptions available upon request. See `LICENSE` for details.

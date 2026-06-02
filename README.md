# BioPack-FQ

**Next-Generation Genomic Storage Format**

BioPack-FQ is a high-performance, open-source binary storage format designed to replace FASTA and FASTQ files for modern sequencing workloads. By encoding nucleotides at 2 bits per base and packing quality scores at 6 bits per value, it achieves substantial storage reduction while preserving lossless roundtrip fidelity.

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
- [3. Results](#3-results)
  - [3.1 Compression Benchmarks](#31-compression-benchmarks)
  - [3.2 Lossless Roundtrip Validation](#32-lossless-roundtrip-validation)
- [4. Usage](#4-usage)
- [5. Build](#5-build)
- [6. Technology Stack](#6-technology-stack)
- [7. Future Work](#7-future-work)
- [8. License](#8-license)

---

## Abstract

Next-generation sequencing produces vast quantities of data, commonly stored in FASTA or FASTQ formats. These human-readable formats are highly inefficient, using 8 bits per nucleotide when only 2 bits are theoretically required. BioPack-FQ addresses this by introducing a compact binary format (`.g` for sequences, `.gq` for sequences with quality scores) that encodes DNA bases as 2-bit values, packs Phred quality scores as 6-bit values, and stores ambiguous bases (N) via a position index. On a 16 MB Oxford Nanopore dataset comprising 29,170 reads and 5.99 million bases, BioPack-FQ achieves a **62.9% storage reduction** without compression. The format supports gzipped input files and provides comprehensive metadata for quality assessment, including Q20/Q30 metrics.

---

## 1. Introduction

FASTA and FASTQ have been the cornerstone of genomic data storage for decades. Despite their ubiquity, their ASCII-based representation incurs substantial overhead:

- **Nucleotides** stored as 8-bit ASCII characters (A, C, G, T) when 2 bits suffice
- **Quality scores** stored as 8-bit ASCII characters (Phred+33 or Phred+64) when 6 bits cover the biologically relevant range (0–60)
- **Read headers** duplicated across millions of entries

For a typical whole-human-genome sequencing run at 30× coverage, raw FASTQ files can exceed 200 GB. Even gzip-compressed FASTQ requires 40–70 GB. Specialized formats such as CRAM achieve better ratios but carry complexity overhead and are tailored for alignment workflows.

BioPack-FQ targets a middle ground: a simple, fast, lossless binary format that replaces FASTA/FASTQ directly without requiring alignment or reference genomes.

---

## 2. Methods

### 2.1 Nucleotide Encoding

Each of the four standard DNA bases is assigned a 2-bit code:

| Base | Binary |
|------|--------|
| A    | 00     |
| C    | 01     |
| G    | 10     |
| T    | 11     |

Four bases are packed into a single byte (big-endian):

```
Byte: [b0:2][b1:2][b2:2][b3:2]
       MSB                  LSB
```

This achieves a 4:1 compression ratio for pure sequence data relative to ASCII FASTA.

### 2.2 Ambiguous Base (N) Handling

Ambiguous bases (N) are not stored in the 2-bit packed stream. Instead, their global positions are recorded in a sorted **N-position list** (none_list) stored after the file header. During encoding, N bases are omitted from the packed data and their positions appended to the list. During decoding, the packed stream and N-position list are merged to reconstruct the original sequence.

This approach adds no overhead for datasets without Ns and only ~0.25 bits per N for datasets with Ns.

### 2.3 Quality Score Encoding

Quality scores (Phred) are auto-detected as Phred+33 or Phred+64 by scanning the minimum ASCII value across all reads.

Each raw quality score (0–63) is stored in **6 bits**. Four scores are packed into three bytes:

```
Byte 0: [s0:6][s1:2]
Byte 1: [s1:4][s2:4]
Byte 2: [s2:2][s3:6]
```

This achieves a 25% reduction over the 8-bit ASCII representation while covering the full Phred score range.

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
│     2-bit Packed Sequence Data   │  A/C/G/T only
├─────────────────────────────────┤
│  6-bit Packed Quality (.gq)      │  Phred scores (optional)
└─────────────────────────────────┘
```

### 2.5 Metadata Header

The 200-byte fixed header stores comprehensive metadata:

| Field | Size | Description |
|-------|------|-------------|
| `version` | 4 B | Format version |
| `file_type` | 1 B | 0 = `.g`, 1 = `.gq` |
| `compression` | 1 B | 0 = none, 1 = ZSTD, 2 = GZIP |
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

---

## 3. Results

### 3.1 Compression Benchmarks

Test dataset: Oxford Nanopore *rbcLa* amplicon sequencing, 29,170 reads, 5,989,618 total bases, Phred+33 qualities.

| Format | Size | Reduction (vs. FASTQ) |
|--------|------|-----------------------|
| Original FASTQ | 16,450,568 B | — |
| BioPack-FQ (`.gq`) | 6,106,504 B | **62.9%** |

BioPack-FQ achieves this reduction without compression (ZSTD/GZIP compression is a planned feature). The primary savings come from 2-bit nucleotide packing (75% reduction on sequence data) and 6-bit quality packing (25% reduction on quality data), offset by headers and read-length metadata.

### 3.2 Impact on Bioinformatics Pipelines

BioPack-FQ is a **storage and transfer format**, not a compute accelerator. Core pipeline stages — alignment, sorting, variant calling — are CPU-bound and will not run faster because the input is stored in `.gq` instead of FASTQ. The format's value lies in the data logistics that surround computation.

#### What BioPack-FQ Does Not Do

- Does not accelerate alignment (BWA, minimap2, Bowtie)
- Does not accelerate variant calling (GATK, FreeBayes, DeepVariant)
- Does not provide random access to individual reads (future feature)
- Does not replace indexing (`.bai`, `.fai`)

#### Where BioPack-FQ Provides Real Benefit

| Domain | Benefit | Why |
|--------|---------|-----|
| **Cloud storage** | 63% less egress/ingress cost | Files are 63% smaller |
| **Cloud transfer** | 63% faster upload/download | Proportional to size reduction |
| **Multi-sample projects** | 100 samples: 15 TB → 5.5 TB | Storage cost savings at scale |
| **Archival (S3 Glacier, tape)** | Lower tier cost per genome | Pay for 55 GB instead of 150 GB |
| **HPC scratch space** | Reduced disk footprint per job | Important on shared filesystems |
| **Concurrent NAS access** | Less I/O contention per read | Multiple pipelines reading simultaneously |
| **Output of basecallers** | Compressed on-the-fly without gzip overhead | No read-recompress cycle |

#### Projected I/O Scenarios

| Operation | Format | Data read | Relative I/O time |
|-----------|--------|-----------|-------------------|
| Load into memory | FASTQ | 150 GB | 1.0× |
| Load into memory | gzip FASTQ | 50 GB + decompress | 2–3× slower * |
| Load into memory | `.gq` | 55 GB | **0.37×** |
| FASTQ → aligner pipe | FASTQ | 150 GB read + write /dev/null | 1.0× |
| `.gq` → aligner pipe | `.gq` decode → stdout | 55 GB read + bit unpack | **~0.5×** |

* gzip decompression is CPU-bound and typically reads slower than raw I/O on modern SSDs.

#### Honest Assessment

For a single-sample variant-calling pipeline (GATK best practices, ~24 h total), switching from FASTQ to `.gq` would save **under 30 min** — most of the pipeline is CPU-bound alignment and calling. The real savings come when that sample is:

- Uploaded to the cloud for processing
- Stored for 6 months post-analysis
- One of 10,000 samples in a biobank
- Transferred between collaborators

#### Key Takeaway

BioPack-FQ optimizes **data logistics**, not computation. Use it to reduce storage costs, accelerate transfers, and ease I/O pressure in multi-sample environments. Do not expect faster alignment or variant calling.

### 3.3 Lossless Roundtrip Validation

Full roundtrip testing (FASTQ → `.gq` → decoded FASTQ) produces byte-identical sequences and quality scores, verified by `diff` across all 29,170 reads and 5.99 million bases.

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

**Dependencies:** `g++` (C++20), `zlib` (development headers)

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
| Decompression | zlib (gzread) |
| Planned compression | ZSTD / GZIP |
| Planned hashing | OpenSSL / xxHash |
| Planned SIMD | AVX2 / AVX-512 |
| Planned parallelism | std::thread / OpenMP |

---

## 7. Future Work

- **ZSTD / GZIP compression** of the payload block
- **Per-read header storage** for full FASTQ roundtrip with read names
- **Random-access indexing** for querying specific reads without full decode
- **Multi-threaded encoding/decoding**
- **Streaming support** for pipelining with bioinformatics tools
- **k-mer frequency tables** and **Bloom filters** for rapid dataset comparison
- **SIMD-optimized** packing and unpacking kernels
- **Checksum validation** via MD5 / SHA-256

---

## 8. License

MIT

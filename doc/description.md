# BioPack-FQ: Next-Generation Genomic Storage Format

## Overview

BioPack-FQ is a high-performance genomic storage format designed to replace traditional FASTA and FASTQ files with a compact binary representation optimized for modern sequencing workloads.

The project introduces two new file formats:

* **.g** — Compact genomic sequence format (FASTA replacement)
* **.gq** — Compact genomic sequence + quality format (FASTQ replacement)

can also use fasta.gz or fastaq.gz as input

The primary objective is to significantly reduce storage requirements, improve I/O performance, accelerate data transfers, and provide rich metadata for large-scale genomic analysis.

## Core Innovation

Instead of storing nucleotides as ASCII characters, BioPack-FQ encodes DNA bases using 2 bits:

A → 00
C → 01
G → 10
T → 11

This allows:

* 4 nucleotides per byte
* ~75% storage reduction for sequence data
* Faster reading and writing operations
* Improved memory efficiency

Quality scores are stored using packed binary representations, reducing the storage footprint while preserving the original information.

## File Architecture

Each BioPack-FQ file begins with a structured metadata header containing information about the dataset before the compressed genomic data.

### Metadata Header

The header may contain:

* Format version
* File type (.g or .gq)
* Original filename
* Original file size
* MD5 checksum of original file
* SHA256 checksum (optional)
* Compression method
* Total number of reads
* Total nucleotide count
* Average read length
* Minimum read length
* Maximum read length
* Nucleotide frequencies (A, C, G, T, N)
* GC content percentage
* Quality score statistics
* Creation timestamp
* Software version
* User-defined metadata

### Example Structure

```cpp
struct GQHeader {
    uint32_t version;

    uint64_t original_size;
    uint64_t compressed_size;

    uint64_t read_count;
    uint64_t nucleotide_count;

    uint32_t min_read_length;
    uint32_t max_read_length;
    uint32_t avg_read_length;

    double gc_content;

    uint64_t count_a;
    uint64_t count_c;
    uint64_t count_g;
    uint64_t count_t;
    uint64_t count_n;

    char md5[33];
    char sha256[65];

    uint64_t created_timestamp;
};
```

## Advanced Features

Future versions may include:

* Random-access indexing
* Multi-threaded encoding/decoding
* Block-level compression
* Streaming support
* Embedded statistics
* Read filtering indexes
* k-mer frequency tables
* Bloom filters for rapid dataset comparison
* ZSTD and GZIP support
* Distributed storage compatibility

## Benchmark Targets

### Human Genome WGS 30x

| Format     | Estimated Size |
| ---------- | -------------- |
| FASTQ      | 120–200 GB     |
| gzip FASTQ | 40–70 GB       |
| CRAM       | 10–30 GB       |
| .gq        | 50–80 GB       |
| .gq.gz     | 30–60 GB       |
| .g         | 30–50 GB       |
| .g.gz      | 20–40 GB       |

## Technology Stack

* Language: C++20
* Build System: CMake
* Compression: ZSTD / GZIP
* Hashing: OpenSSL / xxHash
* SIMD Optimization: AVX2 / AVX-512
* Parallelism: std::thread / OpenMP

## Vision

BioPack-FQ aims to become an open, efficient, and extensible genomic storage standard that bridges the gap between raw FASTQ files and highly specialized formats such as CRAM, while remaining simple to implement, portable, and optimized for modern high-performance computing environments.

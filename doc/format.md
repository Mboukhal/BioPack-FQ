# BioPack-FQ Binary Format Specification

## File Types

| Extension | Content |
|-----------|---------|
| `.g`      | Genomic sequences only (FASTA replacement) |
| `.gq`     | Genomic sequences + quality scores (FASTQ replacement) |

## File Layout

```
┌──────────────────────────────────┐
│         GQHeader (200 bytes)      │  ← fixed metadata header
├──────────────────────────────────┤
│      Read Lengths (variable)      │  ← uint32_t per read
├──────────────────────────────────┤
│         N Position List           │  ← uint32_t count + positions
├──────────────────────────────────┤
│    Sequence Data (2-bit packed)   │  ← A/C/G/T only, no Ns
├──────────────────────────────────┤
│  Quality Data (6-bit packed, .gq) │  ← optional, Phred scores
└──────────────────────────────────┘
```

## 1. GQHeader (200 bytes fixed)

All multibyte fields are **little-endian**.

| Offset | Size | Field | Description |
|--------|------|-------|-------------|
| 0      | 4    | `version` | Format version (currently 1) |
| 4      | 1    | `file_type` | 0 = .g, 1 = .gq |
| 5      | 1    | `compression` | 0 = none, 1 = ZSTD, 2 = GZIP |
| 6      | 2    | _reserved_ | Padding |
| 8      | 8    | `original_size` | Uncompressed payload size |
| 16     | 8    | `compressed_size` | Compressed payload size (0 if uncompressed) |
| 24     | 8    | `read_count` | Number of reads |
| 32     | 8    | `nucleotide_count` | Total bases (including Ns) |
| 40     | 4    | `min_read_length` | Shortest read |
| 44     | 4    | `max_read_length` | Longest read |
| 48     | 4    | `avg_read_length` | Mean read length |
| 52     | 8    | `gc_content` | GC percentage (e.g. 44.2) |
| 60     | 8    | `count_a` | A count |
| 68     | 8    | `count_c` | C count |
| 76     | 8    | `count_g` | G count |
| 84     | 8    | `count_t` | T count |
| 92     | 8    | `count_n` | N count |
| 100    | 1    | `phred_offset` | Phred offset (33 or 64) |
| 101    | 1    | `phred_min` | Min quality score |
| 102    | 1    | `phred_max` | Max quality score |
| 103    | 1    | `phred_avg` | Mean quality score |
| 104    | 16   | `checksum` | MD5 (binary) |
| 120    | 8    | `created_timestamp` | Unix timestamp |
| 128    | 72   | _padding_ | Zero-filled to 200 |

**Total: 200 bytes**

## 2. Read Lengths (variable)

```
[read_len_1:4] [read_len_2:4] ... [read_len_N:4]
```

Each is a `uint32_t`. N = `read_count`.

## 3. N Position List (variable)

```
[n_count:4] [n_pos_1:4] [n_pos_2:4] ... [n_pos_M:4]
```

- `n_count` (`uint32_t`): number of N positions
- `n_pos_i` (`uint32_t`): global position of the i-th N (0-based across concatenated sequences)

The list is **sorted in ascending order**.

## 4. Sequence Data (2-bit packed)

Nucleotides A, C, G, T are encoded with 2 bits:

| Base | Code |
|------|------|
| A    | 00   |
| C    | 01   |
| G    | 10   |
| T    | 11   |

All reads are concatenated into one sequence. N bases are **omitted** from the packed data (their positions are tracked in the N Position List).

Four bases packed per byte:

```
Byte: [b0:2][b1:2][b2:2][b3:2]
       MSB                  LSB
```

Packed size = `ceil((nucleotide_count - n_count) * 2 / 8)` bytes.

## 5. Quality Data (.gq only)

### Phred Detection

- Scan all quality chars, find minimum ASCII value
- `min >= 64` → offset = 64 (Illumina 1.3–1.7)
- Otherwise → offset = 33 (Sanger / Illumina 1.8+)

### 6-bit Packing

Each raw score (`ascii - offset`) is stored in 6 bits.
Four scores packed into 3 bytes:

```
Byte 0: [s0:6    | s1:2]
Byte 1: [s1:4    | s2:4]
Byte 2: [s2:2    | s3:6]
```

Packed size = `ceil(nucleotide_count * 6 / 8)` bytes, but
computed as `((nucleotide_count + 3) / 4) * 3`.

## 6. Complete File Examples

### .gq
```
[200-byte GQHeader]
[read_count × 4 bytes read lengths]
[4 bytes N count  +  n_count × 4 bytes N positions]
[ceil((nuc_count − n_count) / 4) bytes packed A/C/G/T]
[((nuc_count + 3) / 4) × 3 bytes packed quality scores]
```

### .g
```
[200-byte GQHeader]
[read_count × 4 bytes read lengths]
[4 bytes N count  +  n_count × 4 bytes N positions]
[ceil((nuc_count − n_count) / 4) bytes packed A/C/G/T]
```

## Compression

When `compression != 0`, the payload (read lengths + N list + sequence data + quality data) is compressed as one block. The GQHeader stays uncompressed.

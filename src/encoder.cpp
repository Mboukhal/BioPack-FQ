#include "encoder.hpp"
#include <fstream>
#include <iostream>
#include <cstring>
#include <sstream>
#include <numeric>
#include <algorithm>
#include <thread>
#include <future>
#include <zlib.h>
#include "zstd.h"

GQHeader Encoder::build_header(const std::string& /*filename*/,
                                const std::vector<std::string>& reads,
                                const std::vector<std::string>& qualities) {
    GQHeader hdr;
    hdr.file_type = qualities.empty() ? 0 : 1;
    hdr.read_count = reads.size();

    uint64_t total_nucs = 0;
    uint32_t min_len = UINT32_MAX;
    uint32_t max_len = 0;

    for (const auto& read : reads) {
        auto len = static_cast<uint32_t>(read.size());
        total_nucs += len;
        if (len < min_len) min_len = len;
        if (len > max_len) max_len = len;
    }

    hdr.nucleotide_count = total_nucs;
    hdr.min_read_length = min_len;
    hdr.max_read_length = max_len;
    hdr.avg_read_length = reads.empty() ? 0 : static_cast<uint32_t>(total_nucs / reads.size());

    uint64_t freq[5] = {0};
    for (const auto& read : reads) {
        for (char c : read) {
            switch (c) {
                case 'A': case 'a': freq[0]++; break;
                case 'C': case 'c': freq[1]++; break;
                case 'G': case 'g': freq[2]++; break;
                case 'T': case 't': freq[3]++; break;
                case 'N': case 'n': freq[4]++; break;
            }
        }
    }

    hdr.count_a = freq[0];
    hdr.count_c = freq[1];
    hdr.count_g = freq[2];
    hdr.count_t = freq[3];
    hdr.count_n = freq[4];

    if (total_nucs > 0) {
        hdr.gc_content = 100.0 * static_cast<double>(freq[1] + freq[2]) / total_nucs;
    }

    if (!qualities.empty()) {
        hdr.phred_offset = static_cast<uint8_t>(detect_phred_offset(qualities));

        int min_q = 255, max_q = 0;
        uint64_t sum_q = 0;
        uint64_t total_q = 0;
        uint64_t q20 = 0, q30 = 0;
        for (const auto& q : qualities) {
            for (char c : q) {
                int score = static_cast<uint8_t>(c) - hdr.phred_offset;
                if (score < 0) score = 0;
                if (score > 255) score = 255;
                auto s = static_cast<uint8_t>(score);
                if (s < min_q) min_q = s;
                if (s > max_q) max_q = s;
                sum_q += s;
                total_q++;
                if (s >= 20) q20++;
                if (s >= 30) q30++;
            }
        }
        hdr.phred_min = static_cast<uint8_t>(min_q);
        hdr.phred_max = static_cast<uint8_t>(max_q);
        hdr.phred_avg = total_q > 0 ? static_cast<uint8_t>(sum_q / total_q) : 0;
        hdr.q20_count = q20;
        hdr.q30_count = q30;
    }

    return hdr;
}

EncodedPayload Encoder::build_payload(const std::vector<std::string>& sequences,
                                       const std::vector<std::string>& qualities) {
    EncodedPayload payload;
    payload.read_lengths.reserve(sequences.size());

    uint64_t total_bases = 0;
    for (const auto& seq : sequences) {
        payload.read_lengths.push_back(static_cast<uint32_t>(seq.size()));
        total_bases += seq.size();
    }
    std::string all_seqs;
    all_seqs.reserve(total_bases);
    for (const auto& seq : sequences) {
        all_seqs += seq;
    }

    // Parallel sequence encoding
    int num_threads = std::thread::hardware_concurrency();
    if (num_threads < 1) num_threads = 1;
    if (total_bases < 65536) num_threads = 1;

    size_t chunk = (total_bases + static_cast<uint64_t>(num_threads) - 1) / num_threads;
    std::vector<std::future<std::pair<std::vector<uint8_t>, std::vector<uint32_t>>>> futures;

    for (int t = 0; t < num_threads; ++t) {
        size_t start = static_cast<size_t>(t) * chunk;
        if (start >= static_cast<size_t>(total_bases)) break;
        size_t end = std::min(start + chunk, static_cast<size_t>(total_bases));
        futures.push_back(std::async(std::launch::async,
            [&all_seqs, start, end]() {
                std::vector<uint32_t> npos;
                auto data = encode_sequence(all_seqs.data() + start, end - start,
                                            static_cast<uint32_t>(start), npos);
                return std::make_pair(std::move(data), std::move(npos));
            }));
    }

    std::vector<std::vector<uint8_t>> seq_chunks;
    seq_chunks.reserve(futures.size());
    for (auto& f : futures) {
        auto result = f.get();
        seq_chunks.push_back(std::move(result.first));
        payload.n_positions.insert(payload.n_positions.end(),
                                   result.second.begin(), result.second.end());
    }

    size_t total_sb = 0;
    for (const auto& c : seq_chunks) total_sb += c.size();
    payload.seq_data.reserve(total_sb);
    for (const auto& c : seq_chunks) {
        payload.seq_data.insert(payload.seq_data.end(), c.begin(), c.end());
    }

    if (!qualities.empty()) {
        int offset = detect_phred_offset(qualities);
        payload.qual_data = pack_qualities(qualities, offset);
    }

    return payload;
}

std::vector<uint8_t> Encoder::encode(const GQHeader& hdr, const EncodedPayload& payload) {
    // Build payload bytes: read_lengths + n_count + n_positions + seq_data + qual_data
    std::vector<uint8_t> body;

    // read_lengths
    for (auto len : payload.read_lengths) {
        body.push_back(static_cast<uint8_t>(len >> 0));
        body.push_back(static_cast<uint8_t>(len >> 8));
        body.push_back(static_cast<uint8_t>(len >> 16));
        body.push_back(static_cast<uint8_t>(len >> 24));
    }

    // n_positions count
    uint32_t n_count = static_cast<uint32_t>(payload.n_positions.size());
    body.push_back(static_cast<uint8_t>(n_count >> 0));
    body.push_back(static_cast<uint8_t>(n_count >> 8));
    body.push_back(static_cast<uint8_t>(n_count >> 16));
    body.push_back(static_cast<uint8_t>(n_count >> 24));

    // n_positions
    for (auto pos : payload.n_positions) {
        body.push_back(static_cast<uint8_t>(pos >> 0));
        body.push_back(static_cast<uint8_t>(pos >> 8));
        body.push_back(static_cast<uint8_t>(pos >> 16));
        body.push_back(static_cast<uint8_t>(pos >> 24));
    }

    // seq_data
    body.insert(body.end(), payload.seq_data.begin(), payload.seq_data.end());

    // qual_data
    body.insert(body.end(), payload.qual_data.begin(), payload.qual_data.end());

    // Compress
    auto hdr_copy = hdr;
    hdr_copy.original_size = body.size();
    auto compressed = compress(body, hdr_copy);
    hdr_copy.compressed_size = compressed.size();
    auto hdr_bytes = serialize_header(hdr_copy);

    std::vector<uint8_t> output;
    output.reserve(hdr_bytes.size() + compressed.size());
    output.insert(output.end(), hdr_bytes.begin(), hdr_bytes.end());
    output.insert(output.end(), compressed.begin(), compressed.end());

    return output;
}

bool Encoder::write_file(const std::string& output_path, const std::vector<uint8_t>& data) {
    std::ofstream out(output_path, std::ios::binary);
    if (!out) return false;
    out.write(reinterpret_cast<const char*>(data.data()), static_cast<std::streamsize>(data.size()));
    return out.good();
}

std::vector<uint8_t> Encoder::serialize_header(const GQHeader& hdr) {
    std::vector<uint8_t> buf(sizeof(GQHeader));
    std::memcpy(buf.data(), &hdr, sizeof(GQHeader));
    return buf;
}

struct DecodeChunkResult {
    std::vector<FastQRead> reads;
    uint32_t non_n_offset;
    size_t n_idx;
};

static std::pair<uint32_t, size_t> decode_start_state(
    uint32_t global_pos, const std::vector<uint32_t>& n_positions) {
    auto it = std::lower_bound(n_positions.begin(), n_positions.end(), global_pos);
    size_t n_idx = static_cast<size_t>(it - n_positions.begin());
    uint32_t non_n_offset = global_pos - static_cast<uint32_t>(n_idx);
    return {non_n_offset, n_idx};
}

static DecodeChunkResult decode_chunk(
    const uint8_t* seq_data,
    const std::vector<uint32_t>& n_positions,
    const std::vector<uint32_t>& read_lengths,
    uint64_t start_read, uint64_t end_read,
    uint32_t start_global_pos,
    uint32_t start_non_n_offset,
    size_t start_n_idx) {
    DecodeChunkResult result;
    result.reads.reserve(static_cast<size_t>(end_read - start_read));
    uint32_t global_pos = start_global_pos;
    uint32_t non_n_offset = start_non_n_offset;
    size_t n_idx = start_n_idx;
    for (uint64_t i = start_read; i < end_read; ++i) {
        FastQRead read;
        read.sequence.reserve(read_lengths[static_cast<size_t>(i)]);
        uint32_t read_end = global_pos + read_lengths[static_cast<size_t>(i)];
        for (uint32_t gp = global_pos; gp < read_end; ++gp) {
            if (n_idx < n_positions.size() && n_positions[n_idx] == gp) {
                read.sequence += 'N';
                n_idx++;
            } else {
                uint8_t byte = seq_data[non_n_offset / 4];
                uint8_t shift = 6 - 2 * (non_n_offset % 4);
                uint8_t v = (byte >> shift) & 0x03;
                read.sequence += nuc_to_char(static_cast<Nucleotide>(v));
                non_n_offset++;
            }
        }
        global_pos = read_end;
        result.reads.push_back(std::move(read));
    }
    result.non_n_offset = non_n_offset;
    result.n_idx = n_idx;
    return result;
}

static std::vector<uint8_t> decompress_zstd(const uint8_t* data, size_t size, size_t original_size) {
    if (size == 0) return {};
    std::vector<uint8_t> result(original_size);
    size_t rc = ZSTD_decompress(result.data(), result.size(), data, size);
    if (ZSTD_isError(rc)) {
        std::cerr << "ZSTD decompression error: " << ZSTD_getErrorName(rc) << "\n";
        return {};
    }
    result.resize(rc);
    return result;
}

std::vector<uint8_t> Encoder::compress(const std::vector<uint8_t>& data, GQHeader& hdr) {
    size_t bound = ZSTD_compressBound(data.size());
    if (ZSTD_isError(bound)) {
        hdr.compression = 0;
        return data;
    }

    ZSTD_CCtx* cctx = ZSTD_createCCtx();
    if (!cctx) {
        hdr.compression = 0;
        return data;
    }

    int num_threads = std::thread::hardware_concurrency();
    if (num_threads < 1) num_threads = 1;
    ZSTD_CCtx_setParameter(cctx, ZSTD_c_nbWorkers, num_threads);

    std::vector<uint8_t> compressed(bound);
    size_t rc = ZSTD_compress2(cctx, compressed.data(), compressed.size(),
                                data.data(), data.size());
    ZSTD_freeCCtx(cctx);

    if (ZSTD_isError(rc) || rc >= data.size()) {
        hdr.compression = 0;
        return data;
    }
    hdr.compression = 1;
    compressed.resize(rc);
    return compressed;
}

bool Encoder::read_header(const std::string& path, GQHeader& hdr_out,
                           std::vector<uint8_t>& body_out) {
    std::ifstream in(path, std::ios::binary);
    if (!in) return false;

    in.read(reinterpret_cast<char*>(&hdr_out), sizeof(GQHeader));
    if (!in) return false;

    uint64_t total_file = static_cast<uint64_t>(in.seekg(0, std::ios::end).tellg());
    uint64_t body_size = total_file - sizeof(GQHeader);
    if (hdr_out.compressed_size > 0 && hdr_out.compressed_size < body_size)
        body_size = hdr_out.compressed_size;

    in.seekg(sizeof(GQHeader));
    std::vector<uint8_t> compressed(body_size);
    in.read(reinterpret_cast<char*>(compressed.data()), static_cast<std::streamsize>(body_size));
    if (!in.good()) return false;

    if (hdr_out.compression != 0) {
        body_out = decompress_zstd(compressed.data(), compressed.size(),
                                   static_cast<size_t>(hdr_out.original_size));
        return !body_out.empty();
    }
    body_out = std::move(compressed);
    return true;
}

static uint32_t read_u32_le(const uint8_t* p) {
    return static_cast<uint32_t>(p[0]) |
           (static_cast<uint32_t>(p[1]) << 8) |
           (static_cast<uint32_t>(p[2]) << 16) |
           (static_cast<uint32_t>(p[3]) << 24);
}

std::vector<FastQRead> Encoder::decode(const GQHeader& hdr, const std::vector<uint8_t>& body) {
    if (hdr.read_count == 0) return {};

    size_t pos = 0;

    // Read lengths
    std::vector<uint32_t> read_lengths(static_cast<size_t>(hdr.read_count));
    for (uint64_t i = 0; i < hdr.read_count; ++i) {
        read_lengths[static_cast<size_t>(i)] = read_u32_le(body.data() + pos);
        pos += 4;
    }

    // N count
    uint32_t n_count = read_u32_le(body.data() + pos);
    pos += 4;

    // N positions
    std::vector<uint32_t> n_positions(n_count);
    for (uint32_t i = 0; i < n_count; ++i) {
        n_positions[i] = read_u32_le(body.data() + pos);
        pos += 4;
    }

    // Seq data size = ceil(non-N bases / 4)
    uint64_t non_n = hdr.nucleotide_count - n_count;
    size_t seq_bytes = (non_n + 3) / 4;
    size_t qual_bytes = body.size() - pos - seq_bytes;

    const uint8_t* seq_data = body.data() + pos;
    pos += seq_bytes;

    const uint8_t* qual_data = hdr.file_type ? body.data() + pos : nullptr;

    // Pre-compute cumulative positions per read
    std::vector<uint32_t> cum_pos(static_cast<size_t>(hdr.read_count) + 1, 0);
    for (uint64_t i = 0; i < hdr.read_count; ++i) {
        cum_pos[static_cast<size_t>(i) + 1] =
            cum_pos[static_cast<size_t>(i)] + read_lengths[static_cast<size_t>(i)];
    }

    // Parallel sequence decode
    int num_threads = std::thread::hardware_concurrency();
    if (num_threads < 1) num_threads = 1;
    if (hdr.read_count < 1000) num_threads = 1;

    uint64_t chunk_reads = (hdr.read_count + static_cast<uint64_t>(num_threads) - 1)
                           / static_cast<uint64_t>(num_threads);
    std::vector<std::future<DecodeChunkResult>> futures;

    for (int t = 0; t < num_threads; ++t) {
        uint64_t chunk_start = static_cast<uint64_t>(t) * chunk_reads;
        if (chunk_start >= hdr.read_count) break;
        uint64_t chunk_end = std::min(chunk_start + chunk_reads, hdr.read_count);

        uint32_t chunk_global = cum_pos[static_cast<size_t>(chunk_start)];
        auto [non_n_off, n_idx_off] = decode_start_state(chunk_global, n_positions);

        futures.push_back(std::async(std::launch::async,
            [seq_data, &n_positions, &read_lengths, chunk_start, chunk_end,
             chunk_global, non_n_off, n_idx_off]() {
                return decode_chunk(seq_data, n_positions, read_lengths,
                                    chunk_start, chunk_end, chunk_global,
                                    non_n_off, n_idx_off);
            }));
    }

    std::vector<FastQRead> reads;
    reads.reserve(static_cast<size_t>(hdr.read_count));
    for (auto& f : futures) {
        auto cr = f.get();
        for (auto& r : cr.reads) {
            reads.push_back(std::move(r));
        }
    }

    // Decode qualities (single-threaded, fast)
    if (hdr.file_type && qual_data) {
        std::vector<size_t> lens;
        lens.reserve(reads.size());
        for (const auto& r : reads) lens.push_back(r.sequence.size());
        auto quals = unpack_qualities(qual_data, qual_bytes, lens, hdr.phred_offset);
        for (size_t i = 0; i < reads.size(); ++i) {
            reads[i].quality = std::move(quals[i]);
        }
    }

    return reads;
}

bool Encoder::write_fastq(const std::string& path, const std::vector<FastQRead>& reads) {
    std::ofstream out(path);
    if (!out) return false;
    for (const auto& r : reads) {
        out << r.header << "\n" << r.sequence << "\n+\n" << r.quality << "\n";
    }
    return out.good();
}

bool Encoder::write_fasta(const std::string& path, const std::vector<FastQRead>& reads) {
    std::ofstream out(path);
    if (!out) return false;
    for (const auto& r : reads) {
        out << ">" << r.header << "\n" << r.sequence << "\n";
    }
    return out.good();
}

std::string Encoder::read_file_content(const std::string& path) {
    if (path.size() >= 3 && path.substr(path.size() - 3) == ".gz") {
        gzFile gz = gzopen(path.c_str(), "rb");
        if (!gz) return {};
        std::string content;
        char buf[65536];
        int n;
        while ((n = gzread(gz, buf, sizeof(buf))) > 0) {
            content.append(buf, static_cast<size_t>(n));
        }
        gzclose(gz);
        return content;
    }
    std::ifstream in(path, std::ios::binary);
    if (!in) return {};
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

std::vector<FastQRead> Encoder::parse_fastq_from(std::istream& in) {
    std::vector<FastQRead> result;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        FastQRead read;
        read.header = std::move(line);
        if (!std::getline(in, read.sequence)) break;
        if (!std::getline(in, line)) break;
        if (!std::getline(in, read.quality)) break;
        result.push_back(std::move(read));
    }
    return result;
}

std::vector<FastQRead> Encoder::parse_fastq(const std::string& path) {
    auto content = read_file_content(path);
    if (content.empty()) return {};
    std::istringstream in(content);
    return parse_fastq_from(in);
}

std::vector<FastQRead> Encoder::parse_fasta_from(std::istream& in) {
    std::vector<FastQRead> result;
    std::string line;
    FastQRead current;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!current.sequence.empty()) {
                result.push_back(std::move(current));
                current = FastQRead();
            }
            current.header = line.substr(1);
        } else {
            current.sequence += line;
        }
    }
    if (!current.sequence.empty()) {
        result.push_back(std::move(current));
    }
    return result;
}

std::vector<FastQRead> Encoder::parse_fasta(const std::string& path) {
    auto content = read_file_content(path);
    if (content.empty()) return {};
    std::istringstream in(content);
    return parse_fasta_from(in);
}

std::vector<std::string> Encoder::extract_sequences(const std::vector<FastQRead>& reads) {
    std::vector<std::string> seqs;
    seqs.reserve(reads.size());
    for (const auto& r : reads) seqs.push_back(r.sequence);
    return seqs;
}

std::vector<std::string> Encoder::extract_qualities(const std::vector<FastQRead>& reads) {
    std::vector<std::string> quals;
    quals.reserve(reads.size());
    for (const auto& r : reads) quals.push_back(r.quality);
    return quals;
}

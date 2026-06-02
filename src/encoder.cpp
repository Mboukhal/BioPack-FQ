#include "encoder.hpp"
#include <fstream>
#include <cstring>
#include <sstream>
#include <numeric>

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

    // Concatenate all sequences into one string, tracking read lengths
    std::string all_seqs;
    uint64_t total_bases = 0;
    for (const auto& seq : sequences) {
        payload.read_lengths.push_back(static_cast<uint32_t>(seq.size()));
        total_bases += seq.size();
    }
    all_seqs.reserve(total_bases);
    for (const auto& seq : sequences) {
        all_seqs += seq;
    }

    // Encode all at once (no per-read padding waste)
    payload.seq_data = encode_sequence(all_seqs, 0, payload.n_positions);

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
    auto compressed = compress(body);

    // Header
    auto hdr_copy = hdr;
    hdr_copy.original_size = body.size();
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

std::vector<uint8_t> Encoder::compress(const std::vector<uint8_t>& data) {
    // TODO: implement ZSTD compression
    return data;
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
    body_out.resize(body_size);
    in.read(reinterpret_cast<char*>(body_out.data()), static_cast<std::streamsize>(body_size));
    return in.good();
}

static uint32_t read_u32_le(const uint8_t* p) {
    return static_cast<uint32_t>(p[0]) |
           (static_cast<uint32_t>(p[1]) << 8) |
           (static_cast<uint32_t>(p[2]) << 16) |
           (static_cast<uint32_t>(p[3]) << 24);
}

std::vector<FastQRead> Encoder::decode(const GQHeader& hdr, const std::vector<uint8_t>& body) {
    std::vector<FastQRead> reads;
    if (hdr.read_count == 0) return reads;
    reads.reserve(hdr.read_count);

    size_t pos = 0;

    // Read lengths
    std::vector<uint32_t> read_lengths(hdr.read_count);
    for (uint64_t i = 0; i < hdr.read_count; ++i) {
        read_lengths[i] = read_u32_le(body.data() + pos);
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

    // Decode sequences (single pass through data and n_positions)
    uint32_t global_pos = 0;
    uint32_t non_n_offset = 0;
    size_t n_idx = 0;
    for (uint64_t i = 0; i < hdr.read_count; ++i) {
        FastQRead read;
        read.sequence.reserve(read_lengths[i]);
        uint32_t read_end = global_pos + read_lengths[i];
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
        reads.push_back(std::move(read));
    }

    // Decode qualities
    if (hdr.file_type && qual_data) {
        std::vector<size_t> lens;
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

std::vector<FastQRead> Encoder::parse_fastq(const std::string& path) {
    std::vector<FastQRead> result;
    std::ifstream in(path);
    if (!in) return result;

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

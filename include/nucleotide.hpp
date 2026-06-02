#pragma once

#include <cstdint>
#include <cassert>
#include <string>
#include <vector>
#include <stdexcept>
#include <algorithm>

// 2-bit encoding: A=00, C=01, G=10, T=11
// N is recorded in a separate none_list and skipped in packed data
enum class Nucleotide : uint8_t {
    A = 0,
    C = 1,
    G = 2,
    T = 3,
    N = 4
};

inline Nucleotide char_to_nuc(char c) {
    switch (c) {
        case 'A': case 'a': return Nucleotide::A;
        case 'C': case 'c': return Nucleotide::C;
        case 'G': case 'g': return Nucleotide::G;
        case 'T': case 't': return Nucleotide::T;
        case 'N': case 'n': return Nucleotide::N;
        default: throw std::invalid_argument("Invalid nucleotide character");
    }
}

inline char nuc_to_char(Nucleotide n) {
    switch (n) {
        case Nucleotide::A: return 'A';
        case Nucleotide::C: return 'C';
        case Nucleotide::G: return 'G';
        case Nucleotide::T: return 'T';
        default: return 'N';
    }
}

// Pack non-N nucleotides (0-3) into a byte (4 per byte)
inline uint8_t pack_nucleotides(uint8_t v0, uint8_t v1, uint8_t v2, uint8_t v3) {
    return (v0 << 6) | (v1 << 4) | (v2 << 2) | v3;
}

// Encode a DNA string: skip Ns, record their positions, pack only A/C/G/T
// n_positions: global N positions appended here
inline std::vector<uint8_t> encode_sequence(const char* seq, size_t len,
                                             uint32_t base_offset,
                                             std::vector<uint32_t>& n_positions) {
    std::vector<uint8_t> result;
    std::vector<uint8_t> bases;
    bases.reserve(len);
    for (size_t i = 0; i < len; ++i) {
        char c = seq[i];
        if (c == 'N' || c == 'n') {
            n_positions.push_back(static_cast<uint32_t>(base_offset + i));
        } else {
            uint8_t v = static_cast<uint8_t>(char_to_nuc(c)) & 0x03;
            bases.push_back(v);
        }
    }
    result.reserve((bases.size() + 3) / 4);
    for (size_t i = 0; i < bases.size(); i += 4) {
        auto v0 = bases[i];
        auto v1 = (i + 1 < bases.size()) ? bases[i + 1] : 0;
        auto v2 = (i + 2 < bases.size()) ? bases[i + 2] : 0;
        auto v3 = (i + 3 < bases.size()) ? bases[i + 3] : 0;
        result.push_back(pack_nucleotides(v0, v1, v2, v3));
    }
    return result;
}

inline std::vector<uint8_t> encode_sequence(const std::string& seq,
                                             uint32_t base_offset,
                                             std::vector<uint32_t>& n_positions) {
    return encode_sequence(seq.data(), seq.size(), base_offset, n_positions);
}

// Decode: reconstruct from packed data + none_list (sorted)
// non_n_offset: how many non-N bases have been read from data before this call
// n_idx_offset: how many N positions have been consumed before this call
// Returns the total number of non-N bases consumed after this call
inline uint32_t decode_sequence(const uint8_t* data,
                                 const std::vector<uint32_t>& n_positions,
                                 uint32_t base_offset,
                                 uint32_t total_nucs,
                                 std::string& out,
                                 uint32_t non_n_offset = 0,
                                 size_t n_idx_offset = 0) {
    out.reserve(out.size() + total_nucs);

    size_t n_idx = n_idx_offset;
    uint32_t data_idx = non_n_offset;
    for (uint32_t i = 0; i < total_nucs; ++i) {
        uint32_t global_pos = base_offset + i;
        if (n_idx < n_positions.size() && n_positions[n_idx] == global_pos) {
            out += 'N';
            n_idx++;
        } else {
            uint8_t byte = data[data_idx / 4];
            uint8_t shift = 6 - 2 * (data_idx % 4);
            uint8_t v = (byte >> shift) & 0x03;
            out += nuc_to_char(static_cast<Nucleotide>(v));
            data_idx++;
        }
    }
    return data_idx;
}

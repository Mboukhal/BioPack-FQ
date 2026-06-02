#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include <algorithm>
#include <climits>

// Detect Phred offset from quality strings
// Returns 64 if all chars >= 64, otherwise 33
inline int detect_phred_offset(const std::vector<std::string>& qualities) {
    int min_char = INT_MAX;
    for (const auto& q : qualities) {
        for (char c : q) {
            if (static_cast<unsigned char>(c) < static_cast<unsigned>(min_char))
                min_char = static_cast<unsigned char>(c);
        }
    }
    return (min_char >= 64) ? 64 : 33;
}

// Pack 4 quality scores (0-63) into 3 bytes
// Layout: [s0:6][s1:2] [s1:4][s2:4] [s2:2][s3:6]
inline void pack_4_scores(uint8_t s0, uint8_t s1, uint8_t s2, uint8_t s3,
                          uint8_t* out) {
    out[0] = (s0 << 2) | (s1 >> 4);
    out[1] = (s1 << 4) | (s2 >> 2);
    out[2] = (s2 << 6) | s3;
}

// Unpack 3 bytes into 4 quality scores
inline void unpack_4_scores(const uint8_t* in, uint8_t* s0, uint8_t* s1,
                            uint8_t* s2, uint8_t* s3) {
    *s0 = in[0] >> 2;
    *s1 = ((in[0] & 0x03) << 4) | (in[1] >> 4);
    *s2 = ((in[1] & 0x0F) << 2) | (in[2] >> 6);
    *s3 = in[2] & 0x3F;
}

inline std::vector<uint8_t> pack_qualities(const std::vector<std::string>& qualities,
                                            int offset) {
    size_t total_scores = 0;
    for (const auto& q : qualities) total_scores += q.size();

    std::vector<uint8_t> result;
    if (total_scores == 0) return result;

    size_t groups = (total_scores + 3) / 4;
    result.resize(groups * 3, 0);

    size_t score_idx = 0;

    for (const auto& q : qualities) {
        for (char c : q) {
            uint8_t raw = static_cast<uint8_t>(c) - offset;
            if (raw > 63) raw = 63;

            size_t g = score_idx / 4;
            size_t o = score_idx % 4;
            size_t base = g * 3;

            switch (o) {
                case 0: result[base]     |= (raw << 2); break;
                case 1: result[base]     |= (raw >> 4);
                        result[base + 1]  = (raw << 4); break;
                case 2: result[base + 1] |= (raw >> 2);
                        result[base + 2]  = (raw << 6); break;
                case 3: result[base + 2] |= raw;        break;
            }

            score_idx++;
        }
    }

    return result;
}

inline std::vector<std::string> unpack_qualities(const uint8_t* data, size_t /*data_bytes*/,
                                                  const std::vector<size_t>& lengths,
                                                  int offset) {
    std::vector<std::string> result;
    result.reserve(lengths.size());

    size_t score_idx = 0;
    for (size_t len : lengths) {
        std::string q;
        q.reserve(len);
        for (size_t i = 0; i < len; ++i) {
            size_t g = score_idx / 4;
            size_t o = score_idx % 4;
            size_t base = g * 3;

            uint8_t raw = 0;
            switch (o) {
                case 0: raw = data[base] >> 2;       break;
                case 1: raw = ((data[base] & 0x03) << 4) | (data[base + 1] >> 4); break;
                case 2: raw = ((data[base + 1] & 0x0F) << 2) | (data[base + 2] >> 6); break;
                case 3: raw = data[base + 2] & 0x3F; break;
            }

            q += static_cast<char>(raw + offset);
            score_idx++;
        }
        result.push_back(std::move(q));
    }

    return result;
}

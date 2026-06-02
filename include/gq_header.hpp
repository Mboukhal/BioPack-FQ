#pragma once

#include <cstdint>
#include <ctime>

static constexpr uint32_t BIOPACK_VERSION = 1;

#pragma pack(push, 1)
struct GQHeader {
    uint32_t version = BIOPACK_VERSION;
    uint8_t file_type = 0;       // 0=.g, 1=.gq
    uint8_t compression = 0;     // 0=none, 1=ZSTD, 2=GZIP
    uint8_t reserved[2] = {0};

    uint64_t original_size = 0;
    uint64_t compressed_size = 0;

    uint64_t read_count = 0;
    uint64_t nucleotide_count = 0;

    uint32_t min_read_length = 0;
    uint32_t max_read_length = 0;
    uint32_t avg_read_length = 0;

    double gc_content = 0.0;

    uint64_t count_a = 0;
    uint64_t count_c = 0;
    uint64_t count_g = 0;
    uint64_t count_t = 0;
    uint64_t count_n = 0;

    uint8_t phred_offset = 33;
    uint8_t phred_min = 0;
    uint8_t phred_max = 0;
    uint8_t phred_avg = 0;

    uint8_t checksum[16] = {0};  // MD5 binary

    uint64_t created_timestamp = 0;

    uint64_t original_file_size = 0;  // original input file size
    uint64_t q20_count = 0;           // bases with qual score >= 20
    uint64_t q30_count = 0;           // bases with qual score >= 30

    uint8_t _pad[48] = {0};          // zero padding to 200 bytes

    GQHeader() {
        created_timestamp = static_cast<uint64_t>(std::time(nullptr));
    }
};
#pragma pack(pop)

static_assert(sizeof(GQHeader) == 200, "GQHeader must be 200 bytes");

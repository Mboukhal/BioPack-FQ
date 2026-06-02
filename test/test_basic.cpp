#include <iostream>
#include <cassert>
#include <cctype>
#include <cstring>
#include "nucleotide.hpp"
#include "quality.hpp"
#include "encoder.hpp"

static int tests_run = 0;
static int tests_passed = 0;

#define TEST(name) do { std::cout << "  " << name << "... "; tests_run++; } while(0)
#define PASS() do { std::cout << "PASS\n"; tests_passed++; } while(0)
#define FAIL(msg) do { std::cout << "FAIL: " << msg << "\n"; return 1; } while(0)

static int test_nucleotide_encoding() {
    TEST("char_to_nuc / nuc_to_char roundtrip");
    for (char c : {'A', 'C', 'G', 'T', 'a', 'c', 'g', 't'}) {
        if (nuc_to_char(char_to_nuc(c)) != static_cast<char>(std::toupper(static_cast<unsigned char>(c))))
            FAIL("roundtrip failed for " << c);
    }
    PASS();
    return 0;
}

static int test_pack_unpack() {
    TEST("encode_sequence / decode_sequence roundtrip");
    std::string seq = "ACGT";
    std::vector<uint32_t> npos;
    auto packed = encode_sequence(seq, 0, npos);
    if (packed.size() != 1) FAIL("expected 1 byte for 4 nucs");
    if (!npos.empty()) FAIL("no N expected");
    std::string decoded;
    decode_sequence(packed.data(), npos, 0, 4, decoded);
    if (decoded != "ACGT") FAIL("got " << decoded);
    PASS();
    return 0;
}

static int test_n_handling() {
    TEST("N positions recorded and reconstructed");
    std::string seq = "ACGTACGTN";
    std::vector<uint32_t> npos;
    auto packed = encode_sequence(seq, 0, npos);
    if (packed.size() != 2) FAIL("expected 2 bytes for 8 non-N bases");
    if (npos.size() != 1 || npos[0] != 8) FAIL("expected N at position 8");
    std::string decoded;
    decode_sequence(packed.data(), npos, 0, 9, decoded);
    if (decoded != "ACGTACGTN") FAIL("got " << decoded);
    PASS();
    return 0;
}

static int test_header_build() {
    TEST("build_header computes stats correctly");
    Encoder enc;
    std::vector<std::string> reads = {"ACGT", "CC", "G"};

    auto hdr = enc.build_header("test.fasta", reads);

    if (hdr.read_count != 3) FAIL("read_count");
    if (hdr.nucleotide_count != 7) FAIL("nucleotide_count");
    if (hdr.min_read_length != 1) FAIL("min_read_length");
    if (hdr.max_read_length != 4) FAIL("max_read_length");
    if (hdr.avg_read_length != 2) FAIL("avg_read_length");
    if (hdr.count_a != 1) FAIL("count_a");
    if (hdr.count_c != 3) FAIL("count_c");
    if (hdr.count_g != 2) FAIL("count_g");
    if (hdr.count_t != 1) FAIL("count_t");

    PASS();
    return 0;
}

static int test_phred_detection() {
    TEST("detect_phred_offset");
    std::vector<std::string> sanger_quals = {"ABC", "!@#"};
    if (detect_phred_offset(sanger_quals) != 33) FAIL("expected Phred+33");
    std::vector<std::string> illum_quals = {"@ABC", "xyz"};
    if (detect_phred_offset(illum_quals) != 64) FAIL("expected Phred+64");
    PASS();
    return 0;
}

static int test_qual_pack_roundtrip() {
    TEST("quality 6-bit pack/unpack roundtrip");
    std::vector<std::string> quals = {"!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ"};
    int offset = 33;
    auto packed = pack_qualities(quals, offset);
    std::vector<size_t> lens;
    for (const auto& q : quals) lens.push_back(q.size());
    auto unpacked = unpack_qualities(packed.data(), packed.size(), lens, offset);
    if (unpacked.size() != quals.size()) FAIL("wrong count");
    for (size_t i = 0; i < quals.size(); ++i) {
        if (unpacked[i] != quals[i]) FAIL("mismatch at " << i);
    }
    PASS();
    return 0;
}

static int test_full_encode() {
    TEST("full encode pipeline");
    Encoder enc;
    std::vector<std::string> seqs = {"ACGT", "NNN"};
    std::vector<std::string> quals = {"IIII", "!!!!"};

    auto hdr = enc.build_header("test.fastq", seqs, quals);
    if (hdr.file_type != 1) FAIL("should be .gq");
    if (hdr.nucleotide_count != 7) FAIL("nucleotide_count");

    auto payload = enc.build_payload(seqs, quals);
    auto data = enc.encode(hdr, payload);
    if (data.size() < 200) FAIL("data too small");

    // Header (200) + body
    // read_lengths: 2 * 4 = 8
    // n_count: 4
    // n_positions: 3 * 4 = 12
    // seq_data: 1 byte (only ACGT, NNN has no ACGT)
    // qual_data: 6 bytes
    // Total body: 8 + 4 + 12 + 1 + 6 = 31
    // Total: 231
    if (data.size() != 231) FAIL("expected 231 bytes, got " << data.size());

    // Verify header at start
    GQHeader readback;
    std::memcpy(&readback, data.data(), sizeof(GQHeader));
    if (readback.read_count != 2) FAIL("header readback read_count");

    PASS();
    return 0;
}

int main() {
    std::cout << "BioPack-FQ Tests\n"
              << "================\n";

    int failures = 0;
    failures += test_nucleotide_encoding();
    failures += test_pack_unpack();
    failures += test_n_handling();
    failures += test_header_build();
    failures += test_phred_detection();
    failures += test_qual_pack_roundtrip();
    failures += test_full_encode();

    std::cout << "================\n"
              << tests_passed << "/" << tests_run << " tests passed\n";

    return failures;
}

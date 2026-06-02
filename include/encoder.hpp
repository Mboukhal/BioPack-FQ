#pragma once

#include "gq_header.hpp"
#include "nucleotide.hpp"
#include "quality.hpp"
#include <string>
#include <vector>
#include <sstream>

struct FastQRead {
    std::string header;
    std::string sequence;
    std::string quality;
};

struct EncodedPayload {
    std::vector<uint32_t> read_lengths;
    std::vector<uint32_t> n_positions;  // global N positions
    std::vector<uint8_t> seq_data;      // 2-bit packed A/C/G/T
    std::vector<uint8_t> qual_data;     // 6-bit packed quality (empty for .g)
};

class Encoder {
public:
    Encoder() = default;

    GQHeader build_header(const std::string& filename,
                          const std::vector<std::string>& reads,
                          const std::vector<std::string>& qualities = {});

    EncodedPayload build_payload(const std::vector<std::string>& sequences,
                                  const std::vector<std::string>& qualities = {});

    std::vector<uint8_t> encode(const GQHeader& hdr, const EncodedPayload& payload);

    bool write_file(const std::string& output_path, const std::vector<uint8_t>& data);

    static bool read_header(const std::string& path, GQHeader& hdr_out,
                            std::vector<uint8_t>& body_out);

    static std::vector<FastQRead> decode(const GQHeader& hdr, const std::vector<uint8_t>& body);

    static std::string read_file_content(const std::string& path);

    static std::vector<FastQRead> parse_fastq(const std::string& path);
    static std::vector<FastQRead> parse_fastq_from(std::istream& in);

    static std::vector<FastQRead> parse_fasta(const std::string& path);
    static std::vector<FastQRead> parse_fasta_from(std::istream& in);
    static std::vector<std::string> extract_sequences(const std::vector<FastQRead>& reads);
    static std::vector<std::string> extract_qualities(const std::vector<FastQRead>& reads);
    static bool write_fastq(const std::string& path, const std::vector<FastQRead>& reads);
    static bool write_fasta(const std::string& path, const std::vector<FastQRead>& reads);

private:
    std::vector<uint8_t> serialize_header(const GQHeader& hdr);
    std::vector<uint8_t> compress(const std::vector<uint8_t>& data);
};

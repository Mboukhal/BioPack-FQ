#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <fstream>
#include <filesystem>
#include <cstring>
#include <ctime>
#include "encoder.hpp"
#include "nucleotide.hpp"

namespace fs = std::filesystem;

enum class FileClass { UNKNOWN, FASTA, FASTQ, G, GQ };

static FileClass classify(const fs::path& p) {
    auto ext = p.extension().string();
    if (ext == ".gz") {
        auto stem = p.stem();
        return classify(stem);
    }
    if (ext == ".fasta" || ext == ".fa") return FileClass::FASTA;
    if (ext == ".fastq" || ext == ".fq" || ext == ".fastaq") return FileClass::FASTQ;
    if (ext == ".g") return FileClass::G;
    if (ext == ".gq") return FileClass::GQ;
    return FileClass::UNKNOWN;
}

static fs::path strip_gz(const fs::path& p) {
    if (p.extension() == ".gz") return p.stem();
    return p;
}

static fs::path infer_output(const fs::path& input, FileClass fc, bool decode) {
    auto base = strip_gz(input);
    if (decode) {
        if (fc == FileClass::G || fc == FileClass::FASTA) {
            return base.replace_extension(".fasta");
        }
        return base.replace_extension(".fastq");
    }
    if (fc == FileClass::FASTA) {
        return base.replace_extension(".g");
    }
    return base.replace_extension(".gq");
}

static void print_usage(const char* prog) {
    std::cerr << "BioPack-FQ v" << BIOPACK_VERSION << " — Next-Gen Genomic Storage\n"
              << "Usage:\n"
              << "  " << prog << " <input> [-d] [-v] [-o output]\n"
              << "\n"
              << "  .fasta/.fa/.fasta.gz  → encode to .g\n"
              << "  .fastq/.fq/.fastq.gz  → encode to .gq\n"
              << "  .g/.g.gz              → show info\n"
              << "  .g/.g.gz -v           → view full content\n"
              << "  .g/.g.gz -d           → decode to .fasta\n"
              << "  .gq/.gq.gz -v         → view full content\n"
              << "  .gq/.gq.gz -d         → decode to .fastq\n";
}

int main(int argc, char** argv) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    fs::path input = argv[1];
    bool decode = false;
    bool verbose = false;
    fs::path output;

    for (int i = 2; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "-d") {
            decode = true;
        } else if (a == "-v") {
            verbose = true;
        } else if (a == "-o" && i + 1 < argc) {
            output = argv[++i];
        }
    }

    auto fc = classify(input);
    if (fc == FileClass::UNKNOWN) {
        std::cerr << "Error: unrecognized file type: " << input << "\n";
        return 1;
    }

    // .g/.gq with -v → view (decode to stdout)
    if ((fc == FileClass::G || fc == FileClass::GQ) && verbose) {
        GQHeader hdr;
        std::vector<uint8_t> body;
        if (!Encoder::read_header(input.string(), hdr, body)) {
            std::cerr << "Error: failed to read " << input << "\n";
            return 1;
        }
        auto reads = Encoder::decode(hdr, body);
        if (reads.empty()) {
            std::cerr << "Error: decode produced no reads\n";
            return 1;
        }
        for (const auto& r : reads) {
            if (!r.header.empty())
                std::cout << r.header << "\n";
            else
                std::cout << "@\n";
            std::cout << r.sequence << "\n+\n"
                      << r.quality << "\n";
        }
        return 0;
    }

    // .g/.gq without -d → info
    if ((fc == FileClass::G || fc == FileClass::GQ) && !decode) {
        std::ifstream in(input.string(), std::ios::binary);
        if (!in) {
            std::cerr << "Error: cannot open " << input << "\n";
            return 1;
        }
        GQHeader hdr;
        in.read(reinterpret_cast<char*>(&hdr), sizeof(GQHeader));
        if (!in) {
            std::cerr << "Error: failed to read header\n";
            return 1;
        }

        auto ts = static_cast<std::time_t>(hdr.created_timestamp);
        char time_buf[64];
        std::strftime(time_buf, sizeof(time_buf), "%Y-%m-%d %H:%M:%S", std::gmtime(&ts));

        double nuc_total = static_cast<double>(hdr.nucleotide_count);
        double n_pct = nuc_total > 0 ? 100.0 * hdr.count_n / nuc_total : 0;

        std::cout << "BioPack-FQ v" << hdr.version << "\n"
                  << "File:          " << input.filename() << "\n"
                  << "Type:          " << (hdr.file_type == 1 ? ".gq (genomic+qual)" : ".g (genomic)")
                  << "\n"
                  << "Created:       " << time_buf << "\n"
                  << "Compression:   " << (hdr.compression == 0 ? "none" : hdr.compression == 1 ? "ZSTD" : "GZIP")
                  << "\n"
                  << "── Reads ──\n"
                  << "  Count:       " << hdr.read_count << "\n"
                  << "  Length:      " << hdr.min_read_length << " – " << hdr.max_read_length
                  << " (avg " << hdr.avg_read_length << ")\n"
                  << "── Bases ──\n"
                  << "  Total:       " << hdr.nucleotide_count << "\n"
                  << "  GC content:  " << hdr.gc_content << "%\n"
                  << "  N content:   " << n_pct << "%"
                  << " (" << hdr.count_n << ")\n"
                  << "  A:           " << (100.0 * hdr.count_a / nuc_total) << "%\n"
                  << "  C:           " << (100.0 * hdr.count_c / nuc_total) << "%\n"
                  << "  G:           " << (100.0 * hdr.count_g / nuc_total) << "%\n"
                  << "  T:           " << (100.0 * hdr.count_t / nuc_total) << "%\n";

        if (hdr.file_type == 1) {
            double q20_pct = nuc_total > 0 ? 100.0 * hdr.q20_count / nuc_total : 0;
            double q30_pct = nuc_total > 0 ? 100.0 * hdr.q30_count / nuc_total : 0;
            std::cout << "── Quality ──\n"
                      << "  Phred:       +" << static_cast<int>(hdr.phred_offset) << "\n"
                      << "  Range:       " << static_cast<int>(hdr.phred_min) << " – "
                      << static_cast<int>(hdr.phred_max)
                      << " (avg " << static_cast<int>(hdr.phred_avg) << ")\n"
                      << "  Q20:         " << q20_pct << "%\n"
                      << "  Q30:         " << q30_pct << "%\n";
        }

        if (hdr.original_file_size > 0) {
            auto gq_size = fs::file_size(input);
            double ratio = 100.0 * gq_size / hdr.original_file_size;
            std::cout << "── Storage ──\n"
                      << "  Original:    " << hdr.original_file_size << " bytes\n"
                      << "  Compressed:  " << gq_size << " bytes\n"
                      << "  Savings:     " << (100.0 - ratio) << "%\n";
        }

        return 0;
    }

    // .g/.gq with -d → decode
    if ((fc == FileClass::G || fc == FileClass::GQ) && decode) {
        if (output.empty()) {
            output = infer_output(input, fc, true);
        }
        GQHeader hdr;
        std::vector<uint8_t> body;
        if (!Encoder::read_header(input.string(), hdr, body)) {
            std::cerr << "Error: failed to read " << input << "\n";
            return 1;
        }
        auto reads = Encoder::decode(hdr, body);
        if (reads.empty()) {
            std::cerr << "Error: decode produced no reads\n";
            return 1;
        }
        if (Encoder::write_fastq(output.string(), reads)) {
            std::cout << "Decoded " << input << " -> " << output << " ("
                      << reads.size() << " reads)\n";
        } else {
            std::cerr << "Error: failed to write " << output << "\n";
            return 1;
        }
        return 0;
    }

    // FASTA/FASTQ → encode
    if (output.empty()) {
        output = infer_output(input, fc, false);
    }

    std::cerr << "Parsing " << input << "...\n";
    auto reads = Encoder::parse_fastq(input.string());
    if (reads.empty()) {
        std::cerr << "Error: no reads found or file not found\n";
        return 1;
    }
    std::cerr << "Found " << reads.size() << " reads\n";

    auto sequences = Encoder::extract_sequences(reads);
    auto qualities = Encoder::extract_qualities(reads);

    Encoder enc;
    auto hdr = enc.build_header(input.string(), sequences, qualities);
    hdr.original_file_size = fs::file_size(input);
    auto payload = enc.build_payload(sequences, qualities);
    auto data = enc.encode(hdr, payload);
    if (enc.write_file(output.string(), data)) {
        std::cout << "Encoded " << input << " -> " << output << " ("
                  << data.size() << " bytes)\n";
    } else {
        std::cerr << "Error: failed to write " << output << "\n";
        return 1;
    }

    return 0;
}

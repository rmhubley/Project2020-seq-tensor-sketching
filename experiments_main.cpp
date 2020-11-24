#include "util/args.hpp"
#include "util/modules.hpp"
#include "util/multivec.hpp"
#include "util/seqgen.hpp"
#include "util/timer.hpp"
#include "util/utils.hpp"

#include <filesystem>
#include <fstream>
#include <memory>

namespace fs = std::filesystem;

using namespace ts;

struct KmerModule : public BasicModule {
    void override_module_params() override {
        sig_len = int_pow<size_t>(sig_len, kmer_size);
    }
};

template <class seq_type, class embed_type>
struct SeqGenModule {
    Vec2D<seq_type> seqs;
    Vec<std::string> seq_names;
    string test_id;
    Vec2D<seq_type> kmer_seqs;
    Vec2D<embed_type> mh_sketch;
    Vec2D<embed_type> wmh_sketch;
    Vec2D<embed_type> omh_sketch;
    Vec2D<embed_type> ten_sketch;
    Vec3D<embed_type> slide_sketch;
    Vec3D<embed_type> dists;

    BasicModule basicModules;
    KmerModule kmerModules;
    string output;

    void parse(int argc, char **argv) {
        basicModules.parse(argc, argv);
        basicModules.models_init();
        kmerModules.parse(argc, argv);
        kmerModules.models_init();
        output = basicModules.directory + basicModules.output;
    }

    void write_fasta(Vec2D<seq_type> &seq_vec, bool Abc = false) {
        std::ofstream fo;
        fo.open(output + "seqs.fa");
        test_id = "#" + std::to_string(random());
        fo << test_id << "\n";
        for (size_t si = 0; si < seq_vec.size(); si++) {
            fo << "> " << si << "\n";
            for (size_t i = 0; i < seq_vec[i].size(); i++) {
                if (Abc) {
                    fo << (char)(seq_vec[si][i] + (int)'A');
                } else {
                    fo << seq_vec[si][i] << ",";
                }
            }
            fo << "\n\n";
        }
        fo.close();
    }

    void read_fasta(Vec2D<seq_type> &seq_vec) {
        seq_vec.clear();
        string file = (output + "/seqs.fa");
        std::ifstream infile = std::ifstream(file);
        string line;

        std::getline(infile, line);
        if (line[0] == '#') {
            test_id = line;
            std::getline(infile, line);
        }
        while (line[0] != '>') {
            std::cout << line << "\n";
            std::getline(infile, line);
        }
        string name = line;
        Vec<seq_type> seq;
        while (std::getline(infile, line)) {
            if (line[0] == '>') {
                seq_vec.push_back(seq);
                seq_names.push_back(name);
                seq.clear();
                name = line;
            } else if (!line.empty()) {
                for (char c : line) {
                    int ic = c - (int)'A';
                    seq.push_back(ic);
                }
            }
        }
    }


    void generate_sequences() {
        if (basicModules.mutation_pattern == "pairs") {
            basicModules.seq_gen.genseqs_pairs(seqs);
        } else if (basicModules.mutation_pattern == "linear") {
            basicModules.seq_gen.genseqs_linear(seqs);
        } else if (basicModules.mutation_pattern == "tree") {
            basicModules.seq_gen.genseqs_tree(seqs, basicModules.sequence_seeds);
            //            basicModules.seq_gen.genseqs_tree2(seqs);
        } else {
            std::cerr << " mutation pattern `" << basicModules.mutation_pattern
                      << "` is not valid\n";
            exit(1);
        }
    }

    void compute_sketches() {
        size_t num_seqs = seqs.size();
        kmer_seqs.resize(num_seqs);
        wmh_sketch.resize(num_seqs);
        mh_sketch.resize(num_seqs);
        omh_sketch.resize(num_seqs);
        ten_sketch.resize(num_seqs);
        slide_sketch.resize(num_seqs);
        for (size_t si = 0; si < num_seqs; si++) {
            kmer_seqs[si] = seq2kmer<seq_type, seq_type>(seqs[si], basicModules.kmer_size, basicModules.sig_len);
            minhash(kmer_seqs[si], mh_sketch[si], kmerModules.mh_params);
            weighted_minhash(kmer_seqs[si], wmh_sketch[si], kmerModules.wmh_params);
            if (basicModules.tuple_on_kmer) {
                ordered_minhash_flat(kmer_seqs[si], omh_sketch[si], kmerModules.omh_params);
                //                tensor_sketch(kmer_seqs[si], ten_sketch[si],
                //                kmerModules.tensor_params); tensor_slide_sketch(kmer_seqs[si],
                //                slide_sketch[si], kmerModules.tensor_slide_params);
            } else {
                ordered_minhash_flat(seqs[si], omh_sketch[si], basicModules.omh_params);
            }
            tensor_sketch(seqs[si], ten_sketch[si], basicModules.tensor_params);
            tensor_slide_sketch(seqs[si], slide_sketch[si], basicModules.tensor_slide_params);
        }
    }

    void compute_pairwise_dists() {
        int num_seqs = seqs.size();
        if (basicModules.mutation_pattern == "pairs") {
            dists = new3D<double>(8, num_seqs, 1, -1);
            for (size_t i = 0; i < seqs.size(); i += 2) {
                int j = i + 1;
                dists[0][i][0] = edit_distance(seqs[i], seqs[j]);
                dists[1][i][0] = hamming_dist(mh_sketch[i], mh_sketch[j]);
                dists[2][i][0] = hamming_dist(wmh_sketch[i], wmh_sketch[j]);
                dists[3][i][0] = hamming_dist(omh_sketch[i], omh_sketch[j]);
                dists[4][i][0] = l1_dist(ten_sketch[i], ten_sketch[j]);
                dists[5][i][0] = l1_dist2D_minlen(slide_sketch[i], slide_sketch[j]);
            }
        } else {
            dists = new3D<double>(8, num_seqs, num_seqs, 0);
            for (size_t i = 0; i < seqs.size(); i++) {
                for (size_t j = i + 1; j < seqs.size(); j++) {
                    dists[0][i][j] = edit_distance(seqs[i], seqs[j]);
                    dists[1][i][j] = hamming_dist(mh_sketch[i], mh_sketch[j]);
                    dists[2][i][j] = hamming_dist(wmh_sketch[i], wmh_sketch[j]);
                    dists[3][i][j] = hamming_dist(omh_sketch[i], omh_sketch[j]);
                    dists[4][i][j] = l1_dist(ten_sketch[i], ten_sketch[j]);
                    dists[5][i][j] = l1_dist2D_minlen(slide_sketch[i], slide_sketch[j]);
                }
            }
        }
    }

    void save_output() {
        Vec<string> method_names
                = { "ED", "MH", "WMH", "OMH", "TenSketch", "TenSlide", "Ten2", "Ten2Slide" };
        std::ofstream fo;

        // std::filesystem::remove_all(std::filesystem::path(output));
        // std::filesystem::create_directories(std::filesystem::path(output + "/dists"));
        // std::filesystem::create_directories(std::filesystem::path(output + "/sketches"));
        fs::remove_all(fs::path(output));
        fs::create_directories(fs::path(output + "/dists"));
        fs::create_directories(fs::path(output + "/sketches"));

        fo.open(output + "conf.csv");
        assert(fo.is_open());
        fo << basicModules.config();
        fo.close();

        fo.open(output + "timing.csv");
        assert(fo.is_open());
        fo << Timer::summary();
        fo.close();

        write_fasta(seqs);

        size_t num_seqs = seqs.size();
        for (int m = 0; m < 6; m++) {
            fo.open(output + "dists/" + method_names[m] + ".txt");
            assert(fo.is_open());
            if (basicModules.mutation_pattern == "pairs") {
                for (size_t i = 0; i < num_seqs; i += 2) {
                    size_t j = i + 1;
                    fo << i << ", " << j << ", " << dists[m][i][0] << "\n";
                }
            } else {
                for (size_t i = 0; i < num_seqs; i++) {
                    for (size_t j = i + 1; j < seqs.size(); j++) {
                        fo << i << ", " << j << ", " << dists[m][i][j] << "\n";
                    }
                }
            }
            fo.close();
        }

        fo.open(output + "sketches/mh.txt");
        assert(fo.is_open());
        for (size_t si = 0; si < num_seqs; si++) {
            fo << ">> seq " << si << "\n";
            for (const auto &e : mh_sketch[si]) {
                fo << e << ", ";
            }
            fo << "\n";
        }
        fo.close();

        fo.open(output + "sketches/wmh.txt");
        assert(fo.is_open());
        for (size_t si = 0; si < num_seqs; si++) {
            fo << ">> seq " << si << "\n";
            for (const auto &e : wmh_sketch[si]) {
                fo << e << ", ";
            }
            fo << "\n";
        }
        fo.close();

        fo.open(output + "sketches/omh.txt");
        assert(fo.is_open());
        for (size_t si = 0; si < num_seqs; si++) {
            fo << ">> seq " << si << "\n";
            for (const auto &e : omh_sketch[si]) {
                fo << e << ", ";
            }
            fo << "\n";
        }
        fo.close();

        fo.open(output + "sketches/ten.txt");
        assert(fo.is_open());
        for (size_t si = 0; si < seqs.size(); si++) {
            fo << ">> seq " << si << "\n";
            for (const auto &e : ten_sketch[si]) {
                fo << e << ", ";
            }
            fo << "\n";
        }
        fo.close();

        fo.open(output + "sketches/ten_slide.txt");
        for (size_t si = 0; si < seqs.size(); si++) {
            auto &sk = slide_sketch[si];
            for (size_t dim = 0; dim < sk.size(); dim++) {
                fo << ">> seq: " << si << ", dim: " << dim << "\n";
                for (auto &item : sk[dim])
                    fo << item << ", ";
                fo << "\n";
            }
            fo << "\n";
        }
        fo.close();
    }
};

int main(int argc, char *argv[]) {
    SeqGenModule<int, double> experiment;
    experiment.parse(argc, argv);
    if (experiment.basicModules.show_help) {
        std::cout << experiment.basicModules.description();
    } else {
        experiment.generate_sequences();
        experiment.compute_sketches();
        experiment.compute_pairwise_dists();
        experiment.save_output();
    }
    return 0;
}

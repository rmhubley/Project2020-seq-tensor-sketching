#pragma once

#include "util/utils.hpp"

#include <algorithm>
#include <functional>
#include <iostream>
#include <memory>
#include <random>
#include <utility>

namespace ts { // ts = Tensor Sketch

class SeqGen {
  public:
    SeqGen(uint8_t alphabet_size,
           bool fix_len,
           uint32_t num_seqs,
           uint32_t seq_len,
           uint32_t group_size,
           double max_mutation_rate,
           double min_mutation_rate,
           std::string phylogeny_shape,
           uint64_t seed = 0,
           // RMH: Added to support more sophisticated sequence generation
           double gc_background = 0.5,
           double tandem_fraction = 0.0)
        : alphabet_size(alphabet_size),
          fix_len(fix_len),
          num_seqs(num_seqs),
          seq_len(seq_len),
          group_size(group_size),
          max_mutation_rate(max_mutation_rate),
          min_mutation_rate(min_mutation_rate),
          phylogeny_shape(std::move(phylogeny_shape)),
          seed(seed),
          gc_background(gc_background),
          tandem_fraction(tandem_fraction) {
        assert(group_size >= 2 && "group size<=1 leads to completely independent sequences");

        if ( seed == 0 ) {
          // RMH: The original code used a fixed seed of 341234, which supports only
          //      one set of random sequences (for reproducibility I assume).  I've
          //      replaced this with the ability to set the seed or let it randomize
          //      each time.  Seed value of 0 will be used as a flag to randomize.
          //std::mt19937 gen = std::mt19937(341234);
          std::random_device rd;
          seed = rd();
        }
        gen = std::mt19937(seed);
    }

    /**
     * Generate sequences divided into independent groups of size `group_size`, and store
     * ingroup_pairs within each group in `ingroup_pairs`
     * @tparam T character type
     * @tparam C index type
     * @param seqs generated sequences
     * @param pairs sequence ingroup_pairs within each group
     */
    template <class T>
    Vec2D<T> generate_seqs() {
        if (phylogeny_shape == "pair") { // shape=path: implemented as path & group_size=2
            phylogeny_shape = "path";
            group_size = 2;
        }
        Vec2D<T> seqs;
        seqs.reserve(num_seqs);
        while (seqs.size() < num_seqs) {
            Vec2D<T> group;

            // tree-like: g1->g2, add g2 to pool, g1->g3, g2->g4, add g3, g4 to pool
            if (phylogeny_shape == "tree") {
                group = Vec2D<T>(1);
                // RMH: Added more options for sequence generations
                //random_sequence(group[0], seq_len);
                random_sequence(group[0], seq_len, gc_background, tandem_fraction);
                Vec2D<T> children;
                while (group.size() < group_size) {
                    for (auto &seq : group) {
                        std::vector<T> ch;
                        mutate(seq, ch);
                        children.push_back(seq);
                        children.push_back(ch);
                    }
                    std::swap(group, children);
                    children.clear();
                }
            } else if (phylogeny_shape == "path") { // path-like: g0->g1->g2->g3->...
                group = Vec2D<T>(group_size);
                //random_sequence(group[0], seq_len);
                random_sequence(group[0], seq_len, gc_background, tandem_fraction);
                for (size_t i = 0; i < group_size - 1; i++) {
                    mutate(group[i], group[i + 1]);
                }
            } else if (phylogeny_shape == "star") { // star-like: g0->g1, g0->g2,g0->g3 ...
                group = Vec2D<T>(1);
                //random_sequence(group[0], seq_len);
                random_sequence(group[0], seq_len, gc_background, tandem_fraction);
                for (size_t i = 1; i < group_size; i++) {
                    mutate(group[0], group[i]);
                }
            }

            group.resize(group_size);
            seqs.insert(seqs.end(), group.begin(), group.end());
            if (seqs.size() > num_seqs) {
                seqs.resize(num_seqs);
            }
        }
        return seqs;
    }

    template <class T>
    void ingroup_pairs(std::vector<std::pair<T, T>> &pairs) {
        for (size_t go = 0; go < num_seqs; go += group_size) { // group-offset
            for (size_t i = 0; i < group_size && go + i < num_seqs; i++) { // group-member i
                for (size_t j = i + 1; j < group_size && go + j < num_seqs; j++) { // group-member j
                    pairs.push_back({ go + i, go + j });
                }
            }
        }
    }


  private:
    template <class T>
    void mutate(const std::vector<T> &ref, std::vector<T> &seq) {
        std::uniform_real_distribution<double> unif(min_mutation_rate, max_mutation_rate);
        // RMH: Added more sophisticated mutation with GC background
        //mutate(ref, seq, unif(gen));
        mutate(ref, seq, unif(gen), gc_background);
        if (fix_len)
            make_fix_len(seq);
    }

    /**
     * Mutate seq from ref, by mutating each position with the probability = `rate`
     * @tparam T element type in the sequence
     * @param ref
     * @param seq mutated sequence
     * @param rate probability of mutation at each index
     */
    template <class T>
    void mutate(const std::vector<T> &ref, std::vector<T> &seq, double rate) {
        assert((rate >= 0.0) && (rate <= 1.0) && " rate must be strictly in the range [0,1]");
        // probabilities for each index position: no mutation, insert, delete, substitute
        std::discrete_distribution<int> mut { 1 - rate, rate / 3, rate / 3, rate / 3 };
        // the range chosen such that (sub_char+ref % alphabet_size) will different from ref
        std::uniform_int_distribution<T> sub_char(1, alphabet_size - 1);
        // random character from the alphabet
        std::uniform_int_distribution<T> rand_char(0, alphabet_size - 1);
        for (size_t i = 0; i < ref.size(); i++) {
            switch (mut(gen)) {
                case 0: { // no mutation
                    seq.push_back(ref[i]);
                    break;
                }
                case 1: { // insert
                    seq.push_back(rand_char(gen));
                    i--; // init_tensor_slide_params negate the increment
                    break;
                }
                case 2: { // delete
                    break;
                }
                case 3: { // substitute
                    seq.push_back((sub_char(gen) + ref[i]) % alphabet_size);
                    break;
                }
            }
        }
    }

    /**
     * Mutate seq from ref, by mutating each position with the probability = `rate`
     * @tparam T element type in the sequence
     * @param ref
     * @param seq mutated sequence
     * @param rate probability of mutation at each index
     */
    template <class T>
    void mutate(const std::vector<T> &ref, std::vector<T> &seq, double rate, double gc_level) {
        assert((rate >= 0.0) && (rate <= 1.0) && " rate must be strictly in the range [0,1]");
        // probabilities for each index position: no mutation, insert, delete, substitute
        std::discrete_distribution<int> mut { 1 - rate, rate / 3, rate / 3, rate / 3 };

        // A=0,C=1,G=2,T=3, invalid=5
        // Sub probability for each base
        auto subDists = std::vector<std::discrete_distribution<T>>(4);
        // Probability for A or T
        double at_prob = (1.0 - gc_level) / 2.0;
        // Probability for C or G
        double gc_prob = gc_level / 2.0;
        double total_prob = gc_prob + gc_prob + at_prob;
        subDists[0] = std::discrete_distribution<T>({ 0, gc_prob/total_prob, gc_prob/total_prob, at_prob/total_prob });
        total_prob = at_prob + gc_prob + at_prob;
        subDists[1] = std::discrete_distribution<T>({ at_prob/total_prob, 0, gc_prob/total_prob, at_prob/total_prob });
        subDists[2] = std::discrete_distribution<T>({ at_prob/total_prob, gc_prob/total_prob, 0, at_prob/total_prob });
        total_prob = gc_prob + gc_prob + at_prob;
        subDists[3] = std::discrete_distribution<T>({ at_prob/total_prob, gc_prob/total_prob, gc_prob/total_prob, 0 });

        // random character from the alphabet
        std::discrete_distribution<T> insert_char({ at_prob, gc_prob, gc_prob, at_prob });

        for (size_t i = 0; i < ref.size(); i++) {
            switch (mut(gen)) {
                case 0: { // no mutation
                    seq.push_back(ref[i]);
                    break;
                }
                case 1: { // insert
                    seq.push_back(insert_char(gen));
                    i--; // init_tensor_slide_params negate the increment
                    break;
                }
                case 2: { // delete
                    break;
                }
                case 3: { // substitute
                    seq.push_back(subDists[ref[i]](gen));
                    break;
                }
            }
        }
    }

    template <class T>
    void make_fix_len(std::vector<T> &seq) {
        std::uniform_int_distribution<T> rand_char(0, alphabet_size - 1);
        if (seq.size() > seq_len) {
            // RMH: This appears to be a bug, as it doesn't actually truncate the sequence
            //seq = std::vector<T>(seq.begin(), seq.end());
            seq = std::vector<T>(seq.begin(), seq.end() - ( seq.size() - seq_len ));
        } else if (seq.size() < seq_len) {
            while (seq.size() < seq_len) {
                seq.push_back(rand_char(gen));
            }
        }
    }

    /**
     * Generate a random sequence of length `len`
     * @tparam T
     * @param seq : the result will be stored in `seq`
     * @param len : length of the random sequence
     */
    template <class T>
    void random_sequence(std::vector<T> &seq, size_t len) {
        seq.resize(len);
        std::uniform_int_distribution<T> rand_char(0, alphabet_size - 1);
        for (uint32_t i = 0; i < len; i++) {
            seq[i] = rand_char(gen);
        }
    }

    /**
     * Generate a random sequence of length `len` with a given GC fraction.
     * @tparam T
     * @param seq : the result will be stored in `seq`
     * @param len : length of the random sequence
     * @param gc_level : fraction of G+C in the sequence (0.0 to 1.0)
     *
     * RMH: Added this function to be able to generate sequences with a specific
     *      GC background.
     */
    template <class T>
    void random_sequence(std::vector<T> &seq, size_t len, double gc_level) {
        assert(gc_level >= 0.0 && gc_level <= 1.0 &&
               "gc_level should be between 0 and 1");

        seq.resize(len);

        // Probability for A or T
        double at_prob = (1.0 - gc_level) / 2.0;
        // Probability for C or G
        double gc_prob = gc_level / 2.0;
        
        // Set up the discrete distribution: A, C, G, T
        // Indices: 0 = A, 1 = C, 2 = G, 3 = T
        std::discrete_distribution<T> dist({ at_prob, gc_prob, gc_prob, at_prob });

        for (size_t i = 0; i < len; i++) {
            seq[i] = dist(gen);
        }
    }

    // Return {p(A), p(C), p(G), p(T)} consistent with gc_level
    // e.g. gc_level=0.6 => {0.2, 0.3, 0.3, 0.2}
    inline std::vector<double> getGcProb(double gc_level) {
        double at_prob = (1.0 - gc_level) / 2.0; // A,T share
        double gc_prob = gc_level / 2.0;        // C,G share
        // Order: A=0, C=1, G=2, T=3
        return {at_prob, gc_prob, gc_prob, at_prob};
    }

    /**
     * Generate a single random base (A, C, G, T) according to gc_level.
     */
    int random_base(double gc_level)
    {
        static thread_local std::discrete_distribution<int> dist; // We'll rebuild each call for demonstration
        static thread_local double last_gc = -1.0;

        // Build the distribution only if gc_level changed or first use
        if (gc_level != last_gc) {
            auto probs = getGcProb(gc_level);
            dist = std::discrete_distribution<int>(probs.begin(), probs.end());
            last_gc = gc_level;
        }

        return dist(gen);
    }

    /**
     * Fill `seq` of length `len` with random bases according to `gc_level`.
     * seq[i] will be in {0,1,2,3} => {A,C,G,T}.
     */
    template <class T>
    void fill_background_gc(std::vector<T> &seq, size_t len, double gc_level)
    {
        seq.resize(len);
        for (size_t i = 0; i < len; i++) {
            seq[i] = random_base(gc_level);
        }
    }

    /**
     * Generate a random pattern of length `L` (1..10) with GC bias.
     */
    template <class T>
    std::vector<T> make_pattern(size_t L, double gc_level)
    {
        std::vector<T> pattern(L);
        for (size_t i = 0; i < L; i++) {
            pattern[i] = random_base(gc_level);
        }
        return pattern;
    }

    /**
     * Randomly choose an integer in [minVal, maxVal], inclusive.
     */
    int rand_int(int minVal, int maxVal)
    {
        std::uniform_int_distribution<int> dist(minVal, maxVal);
        return dist(gen);
    }

    /**
     * Generate a random sequence of length `len` with a given GC level. 
     * Then embed up to `tandem_frac * len` bases of random tandem repeats 
     * (pattern length in [1..10], repeated [2..50] times), placed at random 
     * non-overlapping positions. 
     *
     * @param seq          output sequence in {0,1,2,3} = {A,C,G,T}
     * @param len          desired total length
     * @param gc_level     fraction of G+C in random picks, in [0..1]
     * @param tandem_frac  fraction of the total length that may be covered by tandem repeats
     */
    template <class T>
    void random_sequence(std::vector<T> &seq, size_t len, double gc_level, double tandem_frac)
    {
        // Basic checks
        assert(gc_level >= 0.0 && gc_level <= 1.0 && "gc_level out of [0,1]");
        assert(tandem_frac >= 0.0 && tandem_frac <= 1.0 && "tandem_frac out of [0,1]");

        // 1) Fill the entire sequence with random bases (GC distribution).
        fill_background_gc(seq, len, gc_level);

        // 2) We'll track which positions are "occupied" by a tandem repeat
        //    so we can avoid overlapping. Initially all false.
        std::vector<bool> occupied(len, false);

        // The max total of tandem bases we aim to place
        size_t target_tandem = static_cast<size_t>(std::floor(tandem_frac * len));
        if (target_tandem == 0) {
            return; // No tandem repeats needed
        }

        size_t placed_tandem = 0;

        // We'll attempt placing tandem repeats up to some maximum tries
        // to avoid infinite loops if there's no room left.
        const int MAX_ATTEMPTS = 10000;
        for (int attempt = 0; attempt < MAX_ATTEMPTS; attempt++) {
            // If we already placed enough tandem bases, stop
            if (placed_tandem >= target_tandem) {
                break;
            }

            // a) pick a pattern length L in [1..10]
            int L = rand_int(1, 10);
            // b) pick a repeat count R in [2..50]
            int R = rand_int(2, 50);

            // total bases for this tandem
            size_t tandemLen = static_cast<size_t>(L) * R;

            // If placing this tandem would exceed our leftover tandem budget, skip
            if (placed_tandem + tandemLen > target_tandem) {
                continue;
            }
            // If it doesn't even fit in the sequence, skip
            if (tandemLen > len) {
                continue;
            }

            // c) pick a random start position in [0..(len - tandemLen)]
            //    so we don't go out of bounds
            int startPos = rand_int(0, static_cast<int>(len - tandemLen));

            // d) Check if that region is free (not occupied)
            bool canPlace = true;
            for (int i = 0; i < static_cast<int>(tandemLen); i++) {
                if (occupied[startPos + i]) {
                    canPlace = false;
                    break;
                }
            }
            if (!canPlace) {
                continue; // skip this attempt
            }

            // e) Generate a random pattern of length L with GC bias
            auto pattern = make_pattern<T>(L, gc_level);

            // f) Place the pattern repeated R times into seq[]
            for (int r = 0; r < R; r++) {
                for (int i = 0; i < L; i++) {
                    seq[startPos + r*L + i] = pattern[i];
                    occupied[startPos + r*L + i] = true;
                }
            }

            // g) Update how many tandem bases we've placed
            placed_tandem += tandemLen;
        }

        // If we exit the loop, we either placed enough tandem bases or ran out of attempts.
        // `seq` now contains up to tandem_frac * len bases of random tandem repeats.
    }

    // Helper to print out the numeric {0,1,2,3} as letters
    char baseToChar(int b) {
        switch(b) {
            case 0: return 'A';
            case 1: return 'C';
            case 2: return 'G';
            case 3: return 'T';
            default: return '?';
        }
    }

  private:
    std::mt19937 gen;
    uint8_t alphabet_size;
    bool fix_len;
    uint32_t num_seqs;
    uint32_t seq_len;
    uint32_t group_size;
    double max_mutation_rate;
    double min_mutation_rate;
    std::string phylogeny_shape;
    uint64_t seed;
    double gc_background;
    double tandem_fraction;
};

} // namespace ts

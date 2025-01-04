#pragma once

#include "tensor.hpp"

#include "util/utils.hpp"

#include <cstddef>
#include <vector>

namespace ts {
/**
 * Computes sliding tensor sketches for a given sequence as described in
 * https://www.biorxiv.org/content/10.1101/2020.11.13.381814v1
 * @tparam seq_type the type of elements in the sequences to be sketched.
 */
template <class seq_type>
class TensorSlide : public Tensor<seq_type> {
  public:
    using sketch_type = Vec2D<double>;

    /**
     * @param alphabet_size the number of elements in the alphabet S over which sequences are
     * defined (e.g. 4 for DNA)
     * @param sketch_dim the dimension of the embedded (sketched) space, denoted by D in the paper
     * @param tup_len the length of the subsequences considered for sketching, denoted by t
     * in the paper
     * @param win_len sliding sketches are computed for substrings of size win_len
     * @param stride sliding sketches are computed every stride characters
     * @param seed the seed to initialize the random number generator used for the random hash
     * functions.
     * @param name the name of the algorithm in the output
     */
    TensorSlide(seq_type alphabet_size,
                size_t sketch_dim,
                size_t tup_len,
                size_t win_len,
                size_t stride,
                uint32_t seed,
                const std::string &name = "TSS", 
                bool alt_slide_method = false,
                bool use_optimal_hashes = false,
                size_t verbosity = 0)
        : Tensor<seq_type>(alphabet_size, sketch_dim, tup_len, seed, name, use_optimal_hashes),
          win_len(win_len),
          stride(stride),
          alt_slide_method(alt_slide_method),
          use_optimal_hashes(use_optimal_hashes),
          verbosity(verbosity)
    {
        assert(stride <= win_len && "Stride cannot be larger than the window length");
        assert(tup_len <= win_len && "Tuple length (t) cannot be larger than the window length");
    }

    /**
     * Computes sliding sketches for the given sequence.
     * A sketch is computed every #stride characters on substrings of length #window.
     * @return seq.size()/stride sketches of size #sketch_dim
     */
    Vec2D<double> compute(const std::vector<seq_type> &seq) {
        Timer timer("tensor_slide_sketch");
        Vec2D<double> sketches;
        if (seq.size() < this->subsequence_len) {
            return new2D<double>(seq.size() / this->stride, this->sketch_dim, double(0));
        }
        auto &hashes = this->hashes;
        auto &signs = this->signs;
        auto tup_len = this->subsequence_len;
        // first index: p; second index: q; third index: r
        // p,q go from 1 to tup_len; p==0 and p==tup_len+1 are sentinels for termination condition
        auto T1 = new3D<double>(tup_len + 2, tup_len + 1, this->sketch_dim, 0);
        auto T2 = new3D<double>(tup_len + 2, tup_len + 1, this->sketch_dim, 0);

        for (uint32_t p = 0; p <= tup_len; p++) {
            T1[p + 1][p][0] = 1;
        }

        if ( verbosity > 0 )
          printf("\n");

        // T[p][q] at step i represents the sketch for seq[i-w+1]...seq[i] when only using hash
        // functions 1<=p,p+1,...q<=t, where t is the sketch size
        for (uint32_t i = 0; i < seq.size(); i++) {
            if ( verbosity > 0 ) 
              printf("i=%d, seq[i] = %d\n", i, seq[i]);
            for (uint32_t p = 1; p <= tup_len; p++) {
                // q-p must be smaller than i, hence the min in the condition
                for (uint32_t q = std::min(p + i, (uint32_t)tup_len); q >= p; q--) {
                    double z = (double)(q - p + 1) / std::min(i + 1, win_len + 1);
                    auto r = hashes[q - 1][seq[i]];
                    bool s = signs[q - 1][seq[i]];
                    if ( verbosity > 5 ) 
                      printf("     p=%d, q=%d, r=%d, s=%d, z=%lf\n", p, q, r, s ? 1 : 0, z);
                    if (s) {
                        this->shift_sum_inplace(T1[p][q], T1[p][q - 1], r, z);
                        this->shift_sum_inplace(T2[p][q], T2[p][q - 1], r, z);
                    } else {
                        this->shift_sum_inplace(T1[p][q], T2[p][q - 1], r, z);
                        this->shift_sum_inplace(T2[p][q], T1[p][q - 1], r, z);
                    }
                }
            }
            if ( verbosity > 5 )
              printf("\n");

            if (i >= win_len) { // only start deleting from front after reaching #win_len
                uint32_t ws = i - win_len; // the element to be removed from the sketch
                for (uint32_t diff = 0; diff < tup_len; ++diff) {
                    for (uint32_t p = 1; p <= tup_len - diff; p++) {
                        auto r = hashes[p - 1][seq[ws]];
                        bool s = signs[p - 1][seq[ws]];
                        uint32_t q = p + diff;
                        // this computes t/(w-t); in our case t (the tuple length) is diff+1
                        double z = (double)(diff + 1) / (win_len - diff);
                        if (s) {
                            this->shift_sum_inplace(T1[p][q], T1[p + 1][q], r, -z);
                            this->shift_sum_inplace(T2[p][q], T2[p + 1][q], r, -z);
                        } else {
                            this->shift_sum_inplace(T1[p][q], T2[p + 1][q], r, -z);
                            this->shift_sum_inplace(T2[p][q], T1[p + 1][q], r, -z);
                        }
                    }
                }
            }

            if ( alt_slide_method ) {
              // RMH: The original method creates < win_len sketches at the start of the sketch.  It also
              //      does not guarantee full coverage of the sequence.
              int32_t win_start = i-win_len+1;
              if ( win_start >= 0 && ( win_start % stride == 0 || i == seq.size()-1 )) {
                if ( verbosity > 0 )
                  printf("Pushing data (alt)\n");
                sketches.push_back(diff(T1[1].back(), T2[1].back()));
              }
            } else { // Original
              if ((i + 1) % stride == 0) { // save a sketch every stride times
                if ( verbosity > 0 )
                  printf("Pushing data\n");
                sketches.push_back(diff(T1[1].back(), T2[1].back()));
              }
            }
        }

        if ( verbosity > 1 ) {
          printf("Sketch:\n");
          for (const auto& vec : sketches) {
            for (double val : vec) {
              printf("%lf, ", val);  // Use %lf for printing double values
            }
            printf("\n");
          }
        }

        return sketches;
    }

    double dist(const Vec2D<double> &a, const Vec2D<double> &b) {
        Timer timer("tensor_slide_sketch_dist");
        return l2_dist2D_minlen(a, b);
        // RMH Testing using scalar product
        //return scalar_product_dist2D_minlen(a, b);
    }


  private:
    std::vector<double> diff(const std::vector<double> &a, const std::vector<double> &b) {
        assert(a.size() == b.size());
        std::vector<double> result(a.size());
        for (uint32_t i = 0; i < result.size(); ++i) {
            result[i] = a[i] - b[i];
        }
        return result;
    }

    uint32_t win_len;
    uint32_t stride;
    bool alt_slide_method;
    bool use_optimal_hashes;
    uint32_t verbosity;
};

} // namespace ts

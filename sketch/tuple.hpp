#pragma once

namespace ts { // ts = Tensor Sketch

template <class seq_type, class size_type>
size_type subseq2ind(const Vec<seq_type> &seq, const Vec<size_type> &sub, size_type alphabet_size) {
    size_type ind = 0, coef = 1;
    for (size_type i = 0; i < sub.size(); i++) {
        ind += seq[sub[i]] * coef;
        coef *= alphabet_size;
    }
    return ind;
}

template <class seq_type, class embed_type, class size_type = std::size_t>
void tup_embed(const Seq<seq_type> &seq,
               Vec<embed_type> &embed,
               size_type alphabet_size,
               size_type tup_len) {
    size_type seq_len = seq.size();
    size_type cnt = 0, size = int_pow(alphabet_size, tup_len);
    embed = Vec<embed_type>(size, 0);
    Vec<size_type> sub(tup_len, 0);
    do {
        if (is_ascending(sub)) {
            auto ind = subseq2ind(seq, sub, alphabet_size);
            embed[ind]++;
        }
    } while (increment_sub(sub, seq_len));
}

} // namespace ts
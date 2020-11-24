#include "alphabets.hpp"

#include <algorithm>
#include <iostream>
#include <string>

namespace ts {

template <size_t n>
class log2 {
    static constexpr size_t _log2(size_t x) {
        if (x < 2) {
            return 0;
        } else {
            return _log2(x >> 1) + 1;
        }
    }

  public:
    static constexpr size_t value = _log2(n);
};

constexpr uint8_t alphabet_size_dna = 5;
constexpr char alphabet_dna[] = "ACGTN";
constexpr uint8_t bits_per_char_dna = 3;
constexpr uint8_t char2int_tab_dna[128]
        = { 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 5, 2, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5,
            0, 5, 5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 5, 2, 5, 5, 5, 3,
            5, 5, 5, 5, 5, 5, 0, 5, 5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5 };

inline uint32_t char2int_dna(uint8_t c) {
    return char2int_tab_dna[c];
}

constexpr char alphabet_protein[] = "ABCDEFGHIJKLMNOPQRSTUVWYZX";
constexpr uint8_t alphabet_size_protein = sizeof(alphabet_protein) - 1;
constexpr uint8_t bits_per_char_protein = log2<alphabet_size_protein - 1>::value + 1;
constexpr uint8_t char2int_tab_protein[128]
        = { 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25,
            25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25,
            25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 0,
            1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
            25, 23, 24, 25, 25, 25, 25, 25, 25, 0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12,
            13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 25, 23, 24, 25, 25, 25, 25, 25 };

inline uint32_t char2int_protein(uint8_t c) {
    return char2int_tab_protein[c];
}

enum class AlphabetType { DNA, Protein };

AlphabetType from_string(std::string str) {
    std::transform(str.begin(), str.end(), str.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    if (str == "dna") {
        return AlphabetType::DNA;
    } else if (str == "protein") {
        return AlphabetType::Protein;
    } else {
        throw std::logic_error("Invalid alphabet type");
    }
}

std::function<uint32_t(uint8_t c)> char2int;
const char *alphabet;
uint8_t alphabet_size;
uint8_t bits_per_char;

void init_alphabet(const std::string &alphabet_str) {
    switch (from_string(alphabet_str)) {
        case AlphabetType::DNA:
            char2int = char2int_dna;
            alphabet = alphabet_dna;
            alphabet_size = alphabet_size_dna;
            bits_per_char = bits_per_char_dna;
            return;
        case AlphabetType::Protein:
            char2int = char2int_protein;
            alphabet = alphabet_protein;
            alphabet_size = alphabet_size_protein;
            bits_per_char = bits_per_char_protein;
            return;
        default:
            std::cerr << "Invalid alphabet type: " << alphabet_str << std::endl;
            std::exit(1);
    }
}

} // namespace ts

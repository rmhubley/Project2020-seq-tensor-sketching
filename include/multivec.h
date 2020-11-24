//
// Created by Amir Joudaki on 11/7/20.
//

#ifndef SEQUENCE_SKETCHING_MULTIVEC_H
#define SEQUENCE_SKETCHING_MULTIVEC_H

#include "args.hpp"

namespace SeqSketch {
    using namespace BasicTypes;

    template<class T>
    auto new2D(int d1, int d2, T val = 0) {
        return Vec2D<T>(d1, Vec<T>(d2, val));
    }
    template<class T>
    auto new3D(int d1, int d2, int d3, T val = 0) {
        return Vec3D<T>(d1, new2D(d2, d3, val));
    }

    // utility functions for len^D-dimensional tensor
    template<class T, class = is_u_integral<T>>
    bool increment_sub(Vec<T> &sub, T len) {
        sub[0]++;
        T i = 0;
        while (sub[i] >= len) {
            sub[i++] = 0;
            if (i >= sub.size())
                return false;
            sub[i]++;
        }
        return true;
    }

    template<class T, class = is_u_integral<T>>
    T sub2ind(const Vec<T> &sub, T len) {
        T ind = 0, coef = 1;
        for (size_t i = 0; i < sub.size(); i++) {
            ind += sub[i] * coef;
            coef *= len;
        }
        return ind;
    }

    template <class T, class = is_u_integral<T>>
    T int_pow(T x, T pow) {
      int result = 1;
      for (;;) {
        if (pow & 1)
          result *= x;
        pow >>= 1;
        if (!pow)
          break;
        x *= x;
      }

      return result;
    }
}// namespace SeqSketch

#endif//SEQUENCE_SKETCHING_MULTIVEC_H

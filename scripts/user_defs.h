#ifndef MORTON_BITS_HPP
#define MORTON_BITS_HPP
#include <cstdint>

#define MORTON(a,b) encode(a,b)
#define MORTON3D(a,b,c) encode(a,b,c)
uint64_t split(const uint32_t a);

// Compute the 2d Morton code for a pair of indices
inline uint64_t encode(const uint32_t x, const uint32_t y) {
    return split(x) | split(y) << 1;
}

inline uint64_t encode(const uint32_t x, const uint32_t y, const uint32_t z) {
    return split(x) | split(y) << 1 | split(z) << 2;
}

inline uint64_t split(const uint32_t a) {
    uint64_t x = a;
    x = (x | x << 16) & 0x0000ffff0000ffffUL;
    x = (x | x <<  8) & 0x00ff00ff00ff00ffUL;
    x = (x | x <<  4) & 0x0f0f0f0f0f0f0f0fUL;
    x = (x | x <<  2) & 0x3333333333333333UL;
    x = (x | x <<  1) & 0x5555555555555555UL;
    return x;
}
#endif


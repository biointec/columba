#ifndef LARGE_INTEGER_H
#define LARGE_INTEGER_H

#include <cmath>
#include <cstdint>

/**
 * @brief A class that represents a 128-bit unsigned integer.
 *
 * This class is used to represent a 128-bit unsigned integer. It is used in the
 * BitParallelED class to store the match vectors. The class provides the basic
 * arithmetic and bitwise operators.

 */
class UInt128 {
  public:
    uint64_t high; // The high 64 bits
    uint64_t low;  // The low 64 bits

    UInt128(uint64_t high, uint64_t low) : high(high), low(low) {
    }
    UInt128(uint64_t low) : high(0), low(low) {
    }
    UInt128() : high(0), low(0) {
    }

    // left bit shift operator
    UInt128 operator<<(unsigned int shift) const {
        if (shift == 0) {
            return *this;
        } else if (shift == 64) {
            return UInt128(low, 0);
        } else if (shift < 64) {
            return UInt128((high << shift) | (low >> (64 - shift)),
                           low << shift);
        } else {
            return UInt128(low << (shift - 64), 0);
        }
    }

    // right bit shift operator
    UInt128 operator>>(unsigned int shift) const {
        if (shift == 0) {
            return *this;
        } else if (shift == 64) {
            return UInt128(0, high);
        } else if (shift < 64) {
            return UInt128(high >> shift,
                           (high << (64 - shift)) | (low >> shift));
        } else {
            return UInt128(0, high >> (shift - 64));
        }
    }

    // Right shift assignment
    UInt128& operator>>=(unsigned int shift) {
        *this = *this >> shift;
        return *this;
    }

    // Left shift assignment
    UInt128& operator<<=(unsigned int shift) {
        *this = *this << shift;
        return *this;
    }

    // Addition operator
    UInt128 operator+(const UInt128& other) const {
        uint64_t newLow = low + other.low;
        uint64_t carry = (newLow < low) ? 1 : 0; // Check for overflow
        uint64_t newHigh = high + other.high + carry;
        return UInt128(newHigh, newLow);
    }

    // Subtraction operator
    UInt128 operator-(const UInt128& other) const {
        uint64_t newLow = low - other.low;
        uint64_t borrow = (newLow > low) ? 1 : 0; // Check for underflow
        uint64_t newHigh = high - other.high - borrow;
        return UInt128(newHigh, newLow);
    }

    // Bitwise AND operator
    UInt128 operator&(const UInt128& other) const {
        return UInt128(high & other.high, low & other.low);
    }

    // Bitwise XOR operator
    UInt128 operator^(const UInt128& other) const {
        return UInt128(high ^ other.high, low ^ other.low);
    }

    // Bitwise OR operator
    UInt128 operator|(const UInt128& other) const {
        return UInt128(high | other.high, low | other.low);
    }

    // Bitwise NOT operator
    UInt128 operator~() const {
        return UInt128(~high, ~low);
    }

    // operator |=
    UInt128& operator|=(const UInt128& other) {
        high |= other.high;
        low |= other.low;
        return *this;
    }

    // operator &=
    UInt128& operator&=(const UInt128& other) {
        high &= other.high;
        low &= other.low;
        return *this;
    }

    // operator ^=
    UInt128& operator^=(const UInt128& other) {
        high ^= other.high;
        low ^= other.low;
        return *this;
    }

    // operator ==
    bool operator==(const UInt128& other) const {
        return high == other.high && low == other.low;
    }

    // conversion to bool
    operator bool() const {
        return high || low;
    }

    int popcount() const {
        return __builtin_popcountll(high) + __builtin_popcountll(low);
    }
};

// std::log2 overload for UInt128
namespace std {
inline double log2(const UInt128& x) {
    if (x.high != 0) {
        return ::log2(static_cast<double>(x.high)) + 64;
    } else {
        return ::log2(static_cast<double>(x.low));
    }
}
} // namespace std

#endif
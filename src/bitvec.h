/******************************************************************************
 *  Columba: Approximate Pattern Matching using Search Schemes                *
 *  Copyright (C) 2020-2021 - Luca Renders <luca.renders@ugent.be> and        *
 *                            Jan Fostier <jan.fostier@ugent.be>              *
 *                                                                            *
 *  This program is free software: you can redistribute it and/or modify      *
 *  it under the terms of the GNU Affero General Public License as            *
 *  published by the Free Software Foundation, either version 3 of the        *
 *  License, or (at your option) any later version.                           *
 *                                                                            *
 *  This program is distributed in the hope that it will be useful,           *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU Affero General Public License for more details.                       *
 *                                                                            *
 * You should have received a copy of the GNU Affero General Public License   *
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.     *
 ******************************************************************************/

#ifndef BITVEC_H
#define BITVEC_H

/**
 * The implementation implements the rank9 algorithm as described in
 * S. Vigna, "Broadword Implementation of Rank/Select Queries", WEA 2008
 * It relies on GCC's __builtin_popcountll, so please build this software
 * using the -mpopcnt flag to enable the SSE 4.2 POPCNT instruction.
 */

#include <cassert>
#include <cstdint>
#include <fstream>
#include <string.h>
#include <vector>

// ============================================================================
// BIT REFERENCE CLASS
// ============================================================================

class Bitref {

  private:
    size_t& wordRef; // reference to a word in the bitvector
    size_t bitmask;  // bitmask of the form (1 << bitIdx)

  public:
    /**
     * Constructor
     * @param wordRef Reference to a word in the bitvector
     * @param bitmask Bitmask of the form (1 << bitIdx)
     */
    Bitref(size_t& wordRef, size_t bitmask)
        : wordRef(wordRef), bitmask(bitmask) {
    }

    /**
     * Set bitref to particular value
     * @param val Target value
     * @return Bitref reference after modification
     */
    const Bitref& operator=(bool val) {
        if (val)
            wordRef |= bitmask;
        else
            wordRef &= ~bitmask;
        return *this;
    }

    /**
     * Set bitref to another bitref
     * @param br Another bitref
     * @return Bitref reference
     */
    const Bitref& operator=(const Bitref& br) {
        return this->operator=(bool(br));
    }

    /**
     * Bool conversion operator
     */
    operator bool() const {
        return (wordRef & bitmask) != 0;
    }
};

// ============================================================================
// BIT VECTOR CLASS
// ============================================================================

class Bitvec {

  private:
    size_t N;                   // size of the bitvector
    std::vector<size_t> bv;     // actual bitvector
    std::vector<size_t> counts; // interleaved 1st and 2nd level counts

  public:
    /**
     * Get a bit at a certain position
     * @param p Position
     * @return true or false
     */
    bool operator[](size_t p) const {
        assert(p < N);
        size_t w = p / 64;
        size_t b = p % 64;
        return (bv[w] & (1ull << b)) != 0;
    }

    /**
     * Get a bit reference at a certain position
     * @param p Position
     * @return Bit reference object
     */
    Bitref operator[](size_t p) {
        assert(p < N);
        size_t w = p / 64;
        size_t b = p % 64;
        return Bitref(bv[w], 1ull << b);
    }

    /**
     * Create an index for the bitvector to support fast rank operations
     */
    void index() {
        counts = std::vector<size_t>((bv.size() + 7) / 4, 0ull);

        size_t countL1 = 0, countL2 = 0;
        for (size_t w = 0, q = 0; w < bv.size(); w++) {
            if (w % 8 == 0) { // store the L1 counts
                countL1 += countL2;
                counts[q] = countL1;
                countL2 = __builtin_popcountll(bv[w]);
                q += 2;
            } else { // store the L2 counts
                counts[q - 1] |= (countL2 << (((w % 8) - 1) * 9));
                countL2 += __builtin_popcountll(bv[w]);
            }
        }
    }

    /**
     * Get the number of 1-bits within the range [0...p[ (preceding pos p)
     * @param p Position
     */
    size_t rank(size_t p) const {
        assert(p < N);
        size_t w = p / 64;      // word index
        size_t b = p % 64;      // bit offset
        size_t q = (w / 8) * 2; // counts index

        // add the first-level counts
        size_t rv = counts[q];

        // add the second-level counts
        int64_t t = (w % 8) - 1;
        rv += counts[q + 1] >> (t + (t >> 60 & 8)) * 9 & 0x1FF;

        // add the popcount in the final word
        return rv + __builtin_popcountll((bv[w] << 1) << (63 - b));
    }

    /**
     * Write the bitvector to an open filestream
     * @param ofs Open output filestream
     */
    void write(std::ofstream& ofs) const {
        ofs.write((char*)&N, sizeof(N));
        ofs.write((char*)bv.data(), bv.size() * sizeof(size_t));
        ofs.write((char*)counts.data(), counts.size() * sizeof(size_t));
    }

    /**
     * Read the bitvector from an open filestream
     * @param ifs Open input filestream
     */
    void read(std::ifstream& ifs) {
        ifs.read((char*)&N, sizeof(N));

        bv.resize((N + 63) / 64);
        ifs.read((char*)bv.data(), bv.size() * sizeof(size_t));

        counts.resize((bv.size() + 7) / 4);
        ifs.read((char*)counts.data(), counts.size() * sizeof(size_t));
    }

    /**
     * Return the size of the bitvector
     * @return The size of the bitvector
     */
    size_t size() const {
        return N;
    }

    /**
     * Default constructor, move constructor and move assignment operator
     */
    Bitvec() : N(0){};
    Bitvec(Bitvec&& rhs) = default;
    Bitvec& operator=(Bitvec&& rhs) = default;

    /**
     * Deleted copy constructor and copy assignment operator
     */
    Bitvec(const Bitvec&) = delete;
    Bitvec& operator=(const Bitvec&) = delete;

    /**
     * Constructor
     * @param N Number of bits in the bitvector
     */
    Bitvec(size_t N) : N(N), bv((N + 63) / 64, 0ull) {
    }
};

// ============================================================================
// TEMPLATED INTERLEAVED BIT VECTOR CLASS
// ============================================================================

template <size_t S> // S is the size of the alphabet
class BitvecIntl {

  private:
    size_t N;          // size of the bitvector
    size_t bvSize;     // number of words in the bitvector
    size_t* bv;        // interleaved bitvectors
    size_t countsSize; // number of words in the counts vector
    size_t* counts;    // interleaved 1st and 2nd level counts

    /**
     * Allocate memory for bv and counts
     */
    void allocateMem() {
        // free existing allocations
        free(bv);
        bv = NULL;
        free(counts);
        counts = NULL;

        if (N == 0) { // special case for N == 0
            bvSize = countsSize = 0;
            return;
        }

        const size_t B = 64; // memory alignment in bytes

        // allocate memory for the bitvector (and set to zero)
        bvSize = S * ((N + 63) / 64);
        // numBytes must be an integral multiple of B
        size_t numBytes = ((bvSize * sizeof(size_t) + B - 1) / B) * B;
        bv = (size_t*)aligned_alloc(B, numBytes);
        memset((void*)bv, 0, numBytes);

        // allocate memory for the counts
        countsSize = 2 * S * ((N + 511) / 512);
        // numBytes must be an integral multiple of B
        numBytes = ((countsSize * sizeof(size_t) + B - 1) / B) * B;
        counts = (size_t*)aligned_alloc(B, numBytes);
    }

    /**
     * Swap two BitvecIntl objects
     * @param lhs Left hand size
     * @param rhs Right hand size
     */
    friend void swap(BitvecIntl<S>& lhs, BitvecIntl<S>& rhs) {
        using std::swap;
        swap(lhs.N, rhs.N);
        swap(lhs.bvSize, rhs.bvSize);
        swap(lhs.bv, rhs.bv);
        swap(lhs.countsSize, rhs.countsSize);
        swap(lhs.counts, rhs.counts);
    }

  public:
    /**
     * Get a bit at position p for character c
     * @param c Character index [0,1,...,S[
     * @param p Position
     * @return true or false
     */
    bool operator()(size_t c, size_t p) const {
        assert(c < S);
        assert(p < N);
        size_t w = (p / 64) * S + c;
        size_t b = p % 64;
        return (bv[w] & (1ull << b)) != 0;
    }

    /**
     * Get a bit reference at position p for character c
     * @param c Character index [0,1,...,S[
     * @param p Position
     * @return Bit reference object
     */
    Bitref operator()(size_t c, size_t p) {
        assert(c < S);
        assert(p < N);
        size_t w = (p / 64) * S + c;
        size_t b = p % 64;
        return Bitref(bv[w], 1ull << b);
    }

    /**
     * Create an index for the bitvector to support fast rank operations
     */
    void index() {
        // reset counts to zero
        memset((void*)counts, 0, countsSize * sizeof(size_t));

        for (size_t c = 0; c < S; c++) {
            size_t countL1 = 0, countL2 = 0;
            for (size_t w = c, q = 2 * c; w < bvSize; w += S) {
                size_t numBits = __builtin_popcountll(bv[w]);
                if (w % (8 * S) == c) { // store the L1 counts
                    countL1 += countL2;
                    counts[q] = countL1;
                    countL2 = numBits;
                    q += 2 * S;
                } else { // store the L2 counts
                    size_t L2offs = 9 * ((w / S % 8) - 1);
                    counts[q + 1 - 2 * S] |= (countL2 << L2offs);
                    countL2 += numBits;
                }
            }
        }
    }

    /**
     * Get the number of 1-bits within the range [0...p[ (preceding pos p)
     * @param c Character index [0,1,...,S[
     * @param p Position
     */
    size_t rank(size_t c, size_t p) const {
        assert(c < S);
        assert(p < N);
        size_t w = (p / 64) * S + c;          // word index
        size_t b = p % 64;                    // bit offset
        size_t q = (p / 512) * 2 * S + 2 * c; // counts index

        // add the first-level counts
        size_t rv = counts[q];

        // add the second-level counts
        int64_t t = ((p / 64) % 8) - 1;
        rv += counts[q + 1] >> (t + (t >> 60 & 8)) * 9 & 0x1FF;

        // add the popcount in the final word
        return rv + __builtin_popcountll((bv[w] << 1) << (63 - b));
    }

    /**
     * Write the bitvector to an open filestream
     * @param ofs Open output filestream
     */
    void write(std::ofstream& ofs) const {
        ofs.write((char*)&N, sizeof(N));
        ofs.write((char*)bv, bvSize * sizeof(size_t));
        ofs.write((char*)counts, countsSize * sizeof(size_t));
    }

    /**
     * Read the bitvector from an open filestream
     * @param ifs Open input filestream
     */
    void read(std::ifstream& ifs) {
        ifs.read((char*)&N, sizeof(N));
        allocateMem();
        ifs.read((char*)bv, bvSize * sizeof(size_t));
        ifs.read((char*)counts, countsSize * sizeof(size_t));
    }

    /**
     * Return the size of the bitvector
     * @return The size of the bitvector
     */
    size_t size() const {
        return N;
    }

    /**
     * Constructor
     * @param N Number of bits in the interleaved bitvector per character
     */
    BitvecIntl(size_t N = 0) : N(N), bv(NULL), counts(NULL) {
        allocateMem();
    }

    /**
     * Move constructor
     * @param rhs Right hand size
     */
    BitvecIntl(BitvecIntl<S>&& rhs) : BitvecIntl() {
        swap(*this, rhs);
    }

    /**
     * Move assignment operator (shallow copy)
     * @param rhs Right hand size
     */
    BitvecIntl& operator=(BitvecIntl<S>&& rhs) {
        swap(*this, rhs);
        return *this;
    }

    /**
     * Deleted copy constructor and copy assignment operator
     */
    BitvecIntl(const BitvecIntl<S>&) = delete;
    BitvecIntl& operator=(const BitvecIntl<S>&) = delete;

    /**
     * Destructor
     */
    ~BitvecIntl() {
        free(bv);
        free(counts);
    }
};

// ============================================================================
// TEMPLATED NON-INTERLEAVED BIT VECTOR CLASS
// ============================================================================

template <size_t S> // S is the size of the alphabet
class BitvecNonIntl {

  private:
    size_t N;          // size of the bitvector
    size_t bvSize;     // number of words in the bitvector
    size_t* bv;        // interleaved bitvectors
    size_t countsSize; // number of words in the counts vector
    size_t* counts;    // interleaved 1st and 2nd level counts

    /**
     * Allocate memory for bv and counts
     */
    void allocateMem() {
        // free existing allocations
        free(bv);
        bv = NULL;
        free(counts);
        counts = NULL;

        if (N == 0) { // special case for N == 0
            bvSize = countsSize = 0;
            return;
        }

        const size_t B = 64; // memory alignment in bytes

        // allocate memory for the bitvector (and set to zero)
        bvSize = S * ((N + 63) / 64);
        // numBytes must be an integral multiple of B
        size_t numBytes = ((bvSize * sizeof(size_t) + B - 1) / B) * B;
        bv = (size_t*)aligned_alloc(B, numBytes);
        memset((void*)bv, 0, numBytes);

        // allocate memory for the counts
        countsSize = 2 * S * ((N + 511) / 512);
        // numBytes must be an integral multiple of B
        numBytes = ((countsSize * sizeof(size_t) + B - 1) / B) * B;
        counts = (size_t*)aligned_alloc(B, numBytes);
    }

    /**
     * Swap two BitvecIntl objects
     * @param lhs Left hand size
     * @param rhs Right hand size
     */
    friend void swap(BitvecNonIntl<S>& lhs, BitvecNonIntl<S>& rhs) {
        using std::swap;
        swap(lhs.N, rhs.N);
        swap(lhs.bvSize, rhs.bvSize);
        swap(lhs.bv, rhs.bv);
        swap(lhs.countsSize, rhs.countsSize);
        swap(lhs.counts, rhs.counts);
    }

  public:
    /**
     * Get a bit at position p for character c
     * @param c Character index [0,1,...,S[
     * @param p Position
     * @return true or false
     */
    bool operator()(size_t c, size_t p) const {
        assert(c < S);
        assert(p < N);
        size_t w = (bvSize / S) * c + (p / 64);
        size_t b = p % 64;
        return (bv[w] & (1ull << b)) != 0;
    }

    /**
     * Get a bit reference at position p for character c
     * @param c Character index [0,1,...,S[
     * @param p Position
     * @return Bit reference object
     */
    Bitref operator()(size_t c, size_t p) {
        assert(c < S);
        assert(p < N);
        size_t w = (bvSize / S) * c + (p / 64);
        size_t b = p % 64;
        return Bitref(bv[w], 1ull << b);
    }

    /**
     * Create an index for the bitvector to support fast rank operations
     */
    void index() {
        // reset counts to zero
        memset((void*)counts, 0, countsSize * sizeof(size_t));

        for (size_t c = 0; c < S; c++) {
            size_t countL1 = 0, countL2 = 0;
            for (size_t w = 0, q = (countsSize / S) * c; w < (bvSize / S);
                 w++) {
                if (w % 8 == 0) { // store the L1 counts
                    countL1 += countL2;
                    counts[q] = countL1;
                    countL2 = __builtin_popcountll(bv[w + (bvSize / S) * c]);
                    q += 2;
                } else { // store the L2 counts
                    counts[q - 1] |= (countL2 << (((w % 8) - 1) * 9));
                    countL2 += __builtin_popcountll(bv[w + (bvSize / S) * c]);
                }
            }
        }
    }

    /**
     * Get the number of 1-bits within the range [0...p[ (preceding pos p)
     * @param c Character index [0,1,...,S[
     * @param p Position
     */
    size_t rank(size_t c, size_t p) const {
        assert(c < S);
        assert(p < N);
        size_t w = (bvSize / S) * c + (p / 64);          // word index
        size_t b = p % 64;                               // bit offset
        size_t q = (countsSize / S) * c + 2 * (p / 512); // counts index

        // add the first-level counts
        size_t rv = counts[q];

        // add the second-level counts
        int64_t t = ((p / 64) % 8) - 1;
        rv += counts[q + 1] >> (t + (t >> 60 & 8)) * 9 & 0x1FF;

        // add the popcount in the final word
        return rv + __builtin_popcountll((bv[w] << 1) << (63 - b));
    }

    /**
     * Write the bitvector to an open filestream
     * @param ofs Open output filestream
     */
    void write(std::ofstream& ofs) const {
        ofs.write((char*)&N, sizeof(N));
        ofs.write((char*)bv, bvSize * sizeof(size_t));
        ofs.write((char*)counts, countsSize * sizeof(size_t));
    }

    /**
     * Read the bitvector from an open filestream
     * @param ifs Open input filestream
     */
    void read(std::ifstream& ifs) {
        ifs.read((char*)&N, sizeof(N));
        allocateMem();
        ifs.read((char*)bv, bvSize * sizeof(size_t));
        ifs.read((char*)counts, countsSize * sizeof(size_t));
    }

    /**
     * Return the size of the bitvector
     * @return The size of the bitvector
     */
    size_t size() const {
        return N;
    }

    /**
     * Constructor
     * @param N Number of bits in the interleaved bitvector per character
     */
    BitvecNonIntl(size_t N = 0) : N(N), bv(NULL), counts(NULL) {
        allocateMem();
    }

    /**
     * Move constructor
     * @param rhs Right hand size
     */
    BitvecNonIntl(BitvecNonIntl<S>&& rhs) : BitvecNonIntl() {
        swap(*this, rhs);
    }

    /**
     * Move assignment operator (shallow copy)
     * @param rhs Right hand size
     */
    BitvecNonIntl& operator=(BitvecNonIntl<S>&& rhs) {
        swap(*this, rhs);
        return *this;
    }

    /**
     * Deleted copy constructor and copy assignment operator
     */
    BitvecNonIntl(const BitvecNonIntl<S>&) = delete;
    BitvecNonIntl& operator=(const BitvecNonIntl<S>&) = delete;

    /**
     * Destructor
     */
    ~BitvecNonIntl() {
        free(bv);
        free(counts);
    }
};

#endif
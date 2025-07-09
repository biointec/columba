/******************************************************************************
 *  Columba: Approximate Pattern Matching using Search Schemes                *
 *  Copyright (C) 2020-2024 - Luca Renders <luca.renders@ugent.be> and        *
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

#ifndef BITPARALLELMATRIX_H
#define BITPARALLELMATRIX_H

#include "definitions.h" // for length_t, HAS_UINT128_T
#if HAS_UINT128_T
using UInt128 = __uint128_t;
#else
#include "largeinteger.h" // for UInt128
#endif                    // HAS_UINT128_T
#include "substring.h"    // for Substring

#include <algorithm> // for max, reverse, fill, min
#include <array>     // for array
#include <cassert>   // for assert
#include <cstddef>   // for size_t
#include <cstdint>   // for uint64_t
#include <memory>    // for allocator_traits<>::value_type
#include <tuple>     // for array
#include <vector>    // for vector

#include <fmt/core.h>   // for format
#include <fmt/format.h> // for to_string
#include <string>       // for string
#include <utility>      // for pair

#include "logger.h" // for Logger, logger

//  ============================================================================
//  CLASS BIT-PARALLEL-ED MATRIX
//  ============================================================================

#define N_MATCH_VECTORS (5) // Currently only ACTG and N are supported

template <typename Container> struct RefAccessor {
    const Container& ref;

    char forwardAccessor(size_t i) const {
        return ref[i];
    }
    size_t size() const {
        return ref.size();
    }

    std::string toString() const {
        std::string result;
        for (size_t i = 0; i < ref.size(); ++i) {
            result += forwardAccessor(i);
        }
        return result;
    }
};

/**
 * Interface for bit-parallel edit distance matrices, regardless of the
 * underlying word sizes.
 *
 */
class IBitParallelED {
  protected:
    const static std::vector<char>
        char2idx; ///< static dictionary where characters are mapped to index

    /**
     * Create the character to index mapping
     */
    static std::vector<char> createChar2idx() {
        std::vector<char> c2i(256, N_MATCH_VECTORS);
        c2i['A'] = 0;
        c2i['C'] = 1;
        c2i['G'] = 2;
        c2i['T'] = 3;
        c2i['N'] = 4;
        return c2i;
    }

    /**
     * Wrapper for the popcount function that counts the number of set bits in
     * the word.
     * @param x The 64-bit word to count the number of set bits in
     */
    static int popcount(uint64_t x) {
        return __builtin_popcountll(x);
    }

    /**
     * Wrapper for the popcount function that counts the number of set bits in
     * the word.
     * @param x The 128-bit word to count the number of set bits in
     */
    static int popcount(UInt128 x) {
#if HAS_UINT128_T
        uint64_t high = x >> 64; // Extract the high 64 bits
        uint64_t low =
            x; // Implicitly cast to uint64_t, so extract the low 64 bits
        return __builtin_popcountll(high) + __builtin_popcountll(low);
#else
        return x.popcount();
#endif // HAS_UINT128_T
    }

  public:
    virtual ~IBitParallelED() = default;
    /**
     * Bit-encode the horizontal sequence X. Call this routine BEFORE
     * calling initializeMatrix(). You may call initializeMatrix() multiple
     * times (with different initialization settings) with a fixed sequence
     * X.
     * @param X the substring to encode
     */
    virtual void setSequence(const Substring& X) = 0;

    /**
     * Initialize the alignment matrix
     * @param maxED Maximum edit distance allowed during alignment
     * @param initED Edit distances of column zero (default = 0, 1, ... maxED)
     */
    virtual void initializeMatrix(uint32_t maxED,
                                  const std::vector<uint32_t>& initED = {}) = 0;

    /**
     * Compute a row of the edit distance matrix in a bit-parallel manner
     * @param i row index in range [1...m]
     * @param Y character of Y-sequence at row i
     * @return false if all elements on row i exceed maxED, true otherwise
     */
    virtual bool computeRow(uint32_t i, char Y) = 0;

    /**
     * Find the minimum edit distance value and its position on a row
     * @param i Row index
     * @param jMin Column index at which minimum value is found (output)
     * @param minScore Minimum value (output)
     */
    virtual void findMinimumAtRow(uint32_t i, uint32_t& jMin,
                                  uint32_t& minScore) const = 0;

    /**
     * Check if a certain row contains the final column
     * @param i The row index
     * @return True if the final column is contained in the row
     */
    virtual bool inFinalColumn(const uint32_t i) const = 0;

#ifdef RUN_LENGTH_COMPRESSION // functions related to CIGAR strings
    /**
     * Find the CIGAR string of the alignment of a reference substring to
     * the query sequence the matrix was initialized with. The sequence
     * should be set before calling this function
     * @param ref the reference string that was aligned. Will be accessed using
     * [] operator
     * @param score the alignment score between ref and query
     * @return The CIGAR string of the alignment, must be empty
     */

    virtual void findCIGAR(const std::vector<char>& ref, const uint32_t score,
                           std::string& CIGAR) = 0;

#else  // functions related to CIGAR strings
    /**
     * Find the CIGAR string of the alignment of a reference substring to
     * the query sequence the matrix was initialized with. The sequence
     * should be set before calling this function
     * @param ref the reference string that was aligned. Will be accessed in
     * the forward direction.
     * @param score the alignment score between ref and query
     * @param CIGAR (output) the CIGAR string of the alignment, must be empty
     */
    virtual void findCIGAR(const Substring& ref, const uint32_t score,
                           std::string& CIGAR) = 0;
#endif // end not RUN_LENGTH_COMPRESSION
    /**
     * Trace the alignment and compute CIGAR string
     * @param ref reference sequence (will be accessed in forward direction),
     * pattern P should be set using setSequence(P)(...)$
     * @param refEnd End offset of the reference sequence
     * @param refBegin Begin offset of the reference sequence (output)
     * @param ED Edit distance score associated with this alignment (output)
     * @param CIGAR CIGAR string (output), must be empty
     */
    virtual void traceBack(const Substring& ref, const length_t refEnd,
                           length_t& refBegin, length_t& ED,
                           std::string& CIGAR) const = 0;

    /**
     *
     * Find cluster centers in the final column of the matrix. These centers
     * correspond to alignments that cannot be created from another alignment by
     * adding/removing gaps at the end.
     * @param lastRow the last row of the matrix that was filled in
     * @param refEnds (output) the end positions in the reference subsequence
     * that correspond to cluster centers.
     * @param maxED the maximal allowed edit distance
     * @param minED the minimal allowed edit distance
     */
    virtual void findClusterCenters(const uint32_t lastRow,
                                    std::vector<length_t>& refEnds,
                                    uint32_t maxED, uint32_t minED) const = 0;

    /**
     * Find the element at index (i, j) this procedure is O(1)
     * @param i Row index
     * @param j Column index
     * @return Score at position (i, j)
     */
    virtual uint32_t at(uint32_t i, uint32_t j) const = 0;

    /**
     * Check whether after row i, the alignment involves only vertical gaps.
     * This happens when row i includes the final column n and when all
     * values on row i decrease monotonically
     * @param i The row index
     * @return true of false
     */
    virtual bool onlyVerticalGapsLeft(uint32_t i) const = 0;

    /**
     * Retrieves the first column index that is in the band for the
     * row
     * @param i The row
     * @returns The first column of the row
     */
    virtual uint32_t getFirstColumn(uint32_t i) const = 0;

    /**
     * Get the number of columns in the matrix
     * @return The number of columns in the matrix (== X.size() + 1)
     */
    virtual uint32_t getNumberOfCols() const = 0;

    /**
     * Get the number of rows in the matrix
     * @return The number of rows in the matrix (== Y.size() + 1)
     */
    virtual uint32_t getNumberOfRows() const = 0;

    /**
     * Check whether setSequenceX() has been called
     * @return True or false
     */
    virtual bool sequenceSet() const = 0;

    /**
     * Reset the matrix
     */
    virtual void reset() = 0;

    /**
     * Get the vertical size of the final column
     * @return  the vertical size of the final column
     */
    virtual uint32_t getSizeOfFinalColumn() const = 0;

    /**
     * Get the maximum edit distance that can be computed by a matrix of this
     * type
     */
    static uint32_t getMatrixMaxED();
    /**
     * Get the maximum number of rows in the first column allowed by a matrix of
     * this type
     */
    static uint32_t getMaxFirstColRows();
};

/**
 * BitVectors struct to store the bitvectors for a row in the matrix
 */
template <typename WordType = uint64_t> struct BitVectors {
    WordType HP;    ///< bitvector to indicate which delta_H == +1
    WordType HN;    ///< bitvector to indicate which delta_H == -1
    WordType D0;    ///< bitvector to indicate which delta_D == 0
    WordType RAC;   ///< bitvector to indicate the Rightmost Active Column
                    ///< = rightmost column with a value <= maxED
    uint64_t score; ///< score at the diagonal
};

/**
 * Bit-parallel edit distance matrix implementation. This class is templated and
 * depends on the size of the word that is used to store the bitvectors.
 */
template <typename WordType = uint64_t>
class BitParallelED : public IBitParallelED {

    // check if WordType is uin64_t or UInt128
    static_assert(std::is_same<WordType, uint64_t>::value ||
                      std::is_same<WordType, UInt128>::value,
                  "WordType must be either uint64_t or UInt128");

  private:
    static const uint32_t WORD_SIZE =
        sizeof(WordType) * 8;                         ///< word size in bits
    static const uint32_t BLOCK_SIZE = WORD_SIZE / 2; ///< block size in bits
    ///< maximum supported edit distance
    static const uint32_t MATRIX_MAX_ED = (WORD_SIZE - BLOCK_SIZE - 2) / 3;
    ///< number of leftmost bits (= max size of first column)
    static const uint32_t LEFT = 2 * MATRIX_MAX_ED + 1;
    static const uint32_t DIAG_R0 = 2 * MATRIX_MAX_ED; ///< diagonal offset

    /**
     * Wrapper for the popcount function that counts the number of set bits in
     * the word.
     */

  public:
    /**
     * Constructor
     */
    BitParallelED() {
    }

    /**
     * Destructor
     */
    ~BitParallelED() {
        reset();
    }

    /**
     * @see IBitParallelED::setSequence
     */
    void setSequence(const Substring& X) override final;

    /**
     * @see IBitParallelED::initializeMatrix
     */
    void
    initializeMatrix(uint32_t maxED,
                     const std::vector<uint32_t>& initED = {}) override final;

    /**
     * @see IBitParallelED::computeRow
     */
    bool computeRow(uint32_t i, char Y) override final {
        assert(i > 0);
        assert(i < m);
        assert(char2idx[Y] < N_MATCH_VECTORS);
        assert(sequenceSet());

        // define BLOCK_SIZE as power of two to make sure this is fast:
        const uint32_t b = i / BLOCK_SIZE; // block identifier
        const uint32_t l = i % BLOCK_SIZE; // leftmost relevant bit

        // aliases to the bitvectors of the current row i (will be computed)
        WordType& HP = bv[i].HP;
        WordType& HN = bv[i].HN;

        WordType& D0 = bv[i].D0;
        WordType& RAC = bv[i].RAC;

        // select the right match vector
        const WordType& M = mv[b][char2idx[Y]];

        // copy the input vectors pertaining the previous row i-1
        HP = bv[i - 1].HP;
        HN = bv[i - 1].HN;
        RAC = bv[i - 1].RAC << 1u;

        // if we are entering a new block, shift input vectors to the right
        // so that they align with the current block
        if (i % BLOCK_SIZE == 0) {
            HP >>= BLOCK_SIZE;
            HN >>= BLOCK_SIZE;
            RAC >>= BLOCK_SIZE;
        }

        // compute the 5 bit-vectors that encode the edit distance minScore
        // (Hyyro)s
        D0 = (((M & HP) + HP) ^ HP) | M | HN;
        WordType VP = HN | ~(D0 | HP);
        WordType VN = D0 & HP;
        HP = (VN << 1u) | ~(D0 | (VP << 1u));
        HN = (D0 & (VP << 1u));

        // compute the minScore at the diagonal
        const uint32_t diagBit = l + DIAG_R0;
        bv[i].score =
            bv[i - 1].score + ((D0 & (WordType(1) << diagBit)) ? 0 : 1);

        // update the rightmost active column (Hyyro)
        // if not a match on the previous RAC, the RAC needs to be updated
        if (!(D0 & RAC)) {
            size_t val = 1u;

            while (val > 0) {
                if (HP & RAC)
                    val--;
                if (HN & RAC)
                    val++;
                if (RAC == (WordType(1) << (diagBit - Wv)))
                    return false;
                RAC >>= 1u;
            }
        }

        return true;
    }

    /**
     * @see IBitParallelED::findMinimumAtRow
     */
    void findMinimumAtRow(uint32_t i, uint32_t& jMin,
                          uint32_t& minScore) const override final {
        jMin = getFirstColumn(i);
        minScore = operator()(i, jMin);

        for (uint32_t j = getFirstColumn(i) + 1; j <= getLastColumn(i); j++) {
            uint32_t thisScore = operator()(i, j);
            if (thisScore < minScore) {
                minScore = thisScore;
                jMin = j;
            }
        }
    }

    /**
     * @see IBitParallelED::inFinalColumn
     */
    bool inFinalColumn(const uint32_t i) const override final {
        return i >= getNumberOfRows() - getSizeOfFinalColumn();
    }

    /**
     * Letters in CIGAR string M (match/mismatch), I (insertion), D (deletion)
     * or NOTHING (not initialized)
     */
    enum CIGARstate { M, I, D, NOTHING };

/**
 * @see IBitParallelED::findCIGAR
 */
#ifdef RUN_LENGTH_COMPRESSION

    void findCIGAR(const std::vector<char>& ref, uint32_t score,
                   std::string& CIGAR) override final {
        findCIGARImpl(RefAccessor<std::vector<char>>{ref}, score, CIGAR);
    }

    void findCIGARImpl(const RefAccessor<std::vector<char>>& ref,
                       uint32_t score, std::string& CIGAR) {
#else
    virtual void findCIGAR(const Substring& ref, uint32_t score,
                           std::string& CIGAR) override final {
#endif // RUN_LENGTH_COMPRESSION
        assert(CIGAR.empty());
        // initialize the matrix with the alignment score
        initializeMatrix(score);

        // compute the rows
        for (unsigned int i = 0; i < ref.size(); i++) {
            computeRow(i + 1, ref.forwardAccessor(i));
        }

        // backtrack starting from the final cell
        uint32_t i = ref.size();
        uint32_t j = n - 1;

        // make sure the score is correct
        // if the score is higher than what was found here that means there was
        // an error close to the border of a parts that caused one cluster to be
        // above the lower bound of that part, while another cluster was below
        // that lower bound. Another search will have found the occurrence at
        // this place with the correct score
        assert(operator()(i, j) <= score);

        CIGARstate state = NOTHING;
        std::vector<std::pair<char, uint32_t>> vCIGAR;

        while (j > 0 || i > 0) {

            const uint32_t b = i / BLOCK_SIZE; // block identifier
            const WordType bit = WordType(1)
                                 << ((j - b * BLOCK_SIZE) + DIAG_R0);

            if ((j > 0) && bv[i].HP & bit) { // gap in horizontal
                --j;
                if (state != CIGARstate::I) {
                    vCIGAR.emplace_back('I', 0);
                    state = CIGARstate::I;
                }

            } else if (i > 0 && j > 0 &&
                       ((mv[b][char2idx[ref.forwardAccessor(i - 1)]] |
                         ~bv[i].D0) &
                        bit)) { // diagonal
                --i;
                --j;
                if (state != CIGARstate::M) {
                    vCIGAR.emplace_back('M', 0);
                    state = CIGARstate::M;
                }

            } else { // gap in vertical
                --i;
                if (state != CIGARstate::D) {
                    vCIGAR.emplace_back('D', 0);
                    state = CIGARstate::D;
                }
            }

            vCIGAR.back().second++;
        }

        // reverse the cigar as it was created from end to beginning
        for (auto it = vCIGAR.crbegin(); it != vCIGAR.crend(); ++it) {
            CIGAR += fmt::format("{}{}", it->second, it->first);
        }
    }

    /**
     * @see IBitParallelED::traceBack
     */
    void traceBack(const Substring& ref, const length_t refEnd,
                   length_t& refBegin, length_t& ED,
                   std::string& CIGAR) const override final {

        assert(CIGAR.empty());

        std::vector<std::pair<char, uint32_t>> vCIGAR;

        uint32_t i = refEnd;
        uint32_t j = n - 1;
        ED = operator()(i, j);

        CIGARstate state = NOTHING;

        while (j > 0) {

            const uint32_t b = i / BLOCK_SIZE; // block identifier
            const WordType bit = WordType(1)
                                 << ((j - b * BLOCK_SIZE) + DIAG_R0);

            if ((bv[i].HP & bit)) { // gap in horizontal direction -> insertion
                --j;
                if (state != CIGARstate::I) {
                    vCIGAR.emplace_back('I', 0);
                    state = CIGARstate::I;
                }

            } else if ((i > 0) &&
                       ((mv[b][char2idx[ref.forwardAccessor(i - 1)]] |
                         ~bv[i].D0) &
                        bit)) { // diagonal
                --i;
                --j;
                if (state != CIGARstate::M) {
                    vCIGAR.emplace_back('M', 0);
                    state = CIGARstate::M;
                }

            } else { // gap in vertical direction
                --i;
                if (state != CIGARstate::D) {
                    vCIGAR.emplace_back('D', 0);
                    state = CIGARstate::D;
                }
            }

            vCIGAR.back().second++;
        }

        refBegin = i;

        // reverse the cigar as it was created from end to beginning
        for (auto it = vCIGAR.crbegin(); it != vCIGAR.crend(); ++it) {
            CIGAR += fmt::format("{}{}", it->second, it->first);
        }
    }

    /**
     * @see IBitParallelED::findClusterCenters
     */
    void findClusterCenters(const uint32_t lastRow,
                            std::vector<length_t>& refEnds, uint32_t maxED,
                            uint32_t minED) const override final {
        refEnds.clear();
        refEnds.reserve(getSizeOfFinalColumn());

        uint32_t firstRow = (m - 1) - getSizeOfFinalColumn();

        uint32_t col = n - 1;

        for (uint32_t i = lastRow; i > (m - 1) - getSizeOfFinalColumn(); i--) {
            uint32_t ED = operator()(i, col);
            if (ED > maxED || ED < minED) {
                continue;
            }
            bool betterThanAbove =
                (i == firstRow) || ED <= operator()(i - 1, col);
            bool betterThanBelow =
                (i == lastRow) || ED <= operator()(i + 1, col);
            if (betterThanAbove && betterThanBelow) {
                refEnds.emplace_back(i);
            }
        }
    }

    /**
     * Operator () overloading -- this procedure is O(1)
     * @param i Row index
     * @param j Column index
     * @return Score at position (i, j)
     */
    uint32_t operator()(uint32_t i, uint32_t j) const {
        // make sure i and j are within matrix bounds
        assert(i < m);
        assert(j < n);

        // we need the bits in the range [b,e[ in HN and HP
        const uint32_t bit = (i % BLOCK_SIZE) + DIAG_R0;
        uint32_t b = (i > j) ? bit - (i - j) + 1 : bit + 1;
        uint32_t e = (i > j) ? bit + 1 : bit + (j - i) + 1;

        WordType mask = ((WordType(1) << (e - b)) - WordType(1)) << b;
        int negatives = popcount(bv[i].HN & mask);
        int positives = popcount(bv[i].HP & mask);

        uint32_t score = bv[i].score;
        score += (i > j) ? (negatives - positives) : (positives - negatives);
        return score;
    }

    /**
     * @see IBitParallelED::at
     */
    uint32_t at(uint32_t i, uint32_t j) const override final {
        return operator()(i, j);
    }

    /**
     * @see IBitParallelED::onlyVerticalGapsLeft
     */
    bool onlyVerticalGapsLeft(uint32_t i) const override final {
        assert(i < m);

        if (i + LEFT < n) // if the column n is not yet reached on row i
            return false;

        const uint32_t b = i / BLOCK_SIZE;
        const uint32_t r = i % BLOCK_SIZE;

        // check if all relevant bits for HN are set to 1
        uint32_t bb = DIAG_R0 - Wv + r + 1;
        uint32_t be = DIAG_R0 + n - b * BLOCK_SIZE;

        return (((~bv[i].HN >> bb) << bb) << (WORD_SIZE - be)) == WordType(0);
    }

    /**
     * @see IBitParallelED::getFirstColumn
     */
    uint32_t getFirstColumn(uint32_t i) const override final {
        return (i <= Wv) ? 0u : i - Wv;
    }

    /**
     * @see IBitParallelED::getNumberOfCols
     */
    uint32_t getNumberOfCols() const override final {
        return n;
    }

    /**
     * @see IBitParallelED::getNumberOfRows
     */
    uint32_t getNumberOfRows() const override final {
        return m;
    }

    /**
     * @see IBitParallelED::sequenceSet
     */
    bool sequenceSet() const override final {
        return !mv.empty();
    }

    /**
     * @see IBitParallelED::reset
     */
    void reset() override final {
        mv.clear();
    }

    /**
     * @see IBitParallelED::getSizeOfFinalColumn
     */
    uint32_t getSizeOfFinalColumn() const override final {
        return Wh + Wv + 1;
    }

    /**
     * Print the banded matrix. This is a debugging function.
     * @param maxRow Last row index to print
     */
    void printMatrix(uint32_t maxRow = 500) const;

    /**
     * Get the maximum edit distance that can be computed by a matrix of this
     * type
     */
    static constexpr size_t getMatrixMaxED() {
        return MATRIX_MAX_ED;
    }

    /**
     * Get the maximum number of rows in the first column allowed by a matrix of
     * this type
     */
    static constexpr size_t getMaxFirstColRows() {
        return LEFT;
    }

  private:
    uint32_t maxED; // maximum allowed edit distance
    uint32_t m;     // number of rows
    uint32_t n;     // number of columns
    uint32_t Wv;    // vertical width of the band
    uint32_t Wh;    // horizontal width of the band

    std::vector<BitVectors<WordType>> bv;                  // bitvectors
    std::vector<std::array<WordType, N_MATCH_VECTORS>> mv; // match vectors

    /**
     * Retrieves the last column index that needs to be filled in for the
     * row
     * @param i The row to fill in
     * @returns The last column to fill in
     */
    uint32_t getLastColumn(uint32_t i) const {
        return std::min(n - 1, i + Wh);
    }
};

using BitParallelED64 = BitParallelED<uint64_t>;
using BitParallelED128 = BitParallelED<UInt128>;

#endif // BITPARALLELMATRIX_H
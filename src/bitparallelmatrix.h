/******************************************************************************
 *  Columba 1.2: Approximate Pattern Matching using Search Schemes            *
 *  Copyright (C) 2020-2023 - Luca Renders <luca.renders@ugent.be> and        *
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

#include <algorithm>

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "substring.h"

#include "wordlength.h"
// ============================================================================
// CLASS BIT-PARALLEL-ED MATRIX
// ============================================================================

typedef struct {
    uint64_t HP;    // bit vector to indicate which delta_H == +1
    uint64_t HN;    // bit vector to indicate which delta_H == -1
    uint64_t D0;    // bit vector to indicate which delta_D == 0
    uint64_t RAC;   // bit vector to indicate the Rightmost Active Column
                    // = rightmost column with a value <= maxED
    uint64_t score; // score at the diagonal
} BitVectors;

#define WORD_SIZE (64ull)
#define BLOCK_SIZE (32ull)
#define MAX_ED ((WORD_SIZE - BLOCK_SIZE - 2ull) / 3ull)
#define LEFT (2ull * MAX_ED + 1ull)
#define DIAG_R0 (2ull * MAX_ED)

class BitParallelED {
  public:
    /**
     * Constructor
     */
    BitParallelED() {
    }

    static std::vector<char> createChar2idx() {
        std::vector<char> c2i(256, 4);
        c2i['A'] = 0;
        c2i['C'] = 1;
        c2i['G'] = 2;
        c2i['T'] = 3;
        return c2i;
    }

    /**
     * Bit-encode the horizontal sequence X. Call this routine BEFORE
     * calling initializeMatrix(). You may call initializeMatrix() multiple
     * times (with different initialization settings) with a fixed sequence
     * X.
     * @param X the substring to encode
     */
    void setSequence(const Substring& X);

    /**
     * Initialize the alignment matrix
     * @param maxED Maximum edit distance allowed during alignment
     * @param initED Edit distances of column zero (default = 0, 1, ... maxED)
     */
    void initializeMatrix(uint maxED, const std::vector<uint>& initED = {});

    /**
     * Compute a row of the edit distance matrix in a bit-parallel manner
     * @param i row index in range [1...m]
     * @param Y character of Y-sequence at row i
     * @return false if all elements on row i exceed maxED, true otherwise
     */
    bool computeRow(uint i, char Y) {
        assert(i > 0);
        assert(i < m);
        assert(char2idx[Y] < 4);

        // define BLOCK_SIZE as power of two to make sure this is fast:
        const uint b = i / BLOCK_SIZE; // block identifier
        const uint l = i % BLOCK_SIZE; // leftmost relevant bit

        // aliases to the bit vectors of the current row i (will be computed)
        uint64_t& HP = bv[i].HP;
        uint64_t& HN = bv[i].HN;

        uint64_t& D0 = bv[i].D0;
        uint64_t& RAC = bv[i].RAC;

        // select the right match vector
        const uint64_t& M = mv[b][char2idx[Y]];

        // copy the input vectors pertaining the previous row i-1
        HP = bv[i - 1].HP;
        HN = bv[i - 1].HN;
        RAC = bv[i - 1].RAC << 1;

        // if we are entering a new block, shift input vectors to the right
        // so that they align with the current block
        if (i % BLOCK_SIZE == 0) {
            HP >>= BLOCK_SIZE;
            HN >>= BLOCK_SIZE;
            RAC >>= BLOCK_SIZE;
        }

        // compute the 5 bitvectors that encode the edit distance minScore
        // (Hyyro)s
        D0 = (((M & HP) + HP) ^ HP) | M | HN;
        uint64_t VP = HN | ~(D0 | HP);
        uint64_t VN = D0 & HP;
        HP = (VN << 1) | ~(D0 | (VP << 1));
        HN = (D0 & (VP << 1));

        // compute the minScore at the diagonal
        const size_t diagBit = l + DIAG_R0;
        bv[i].score = bv[i - 1].score + (D0 & (1ull << diagBit) ? 0 : 1);

        // update the rightmost active column (Hyyro)
        // if not a match on the previous RAC, the RAC needs to be updated
        if ((D0 & RAC) == 0) {
            size_t val = 1u;
            while (val > 0) {
                if (HP & RAC)
                    val--;
                if (HN & RAC)
                    val++;
                if (RAC == (1ull << (diagBit - Wv)))
                    return false;
                RAC >>= 1;
            }
        }

        return true;
    }

    /**
     * Find the minimum edit distance value and its position on a row
     * @param i Row index
     * @param jMin Column index at which minimum value is found (output)
     * @param minScore Mimumum value (output)
     */
    void findMinimumAtRow(uint i, uint& jMin, uint& minScore) const {
        jMin = getFirstColumn(i);
        minScore = operator()(i, jMin);

        for (uint j = getFirstColumn(i) + 1; j <= getLastColumn(i); j++) {
            uint thisScore = operator()(i, j);
            if (thisScore < minScore) {
                minScore = thisScore;
                jMin = j;
            }
        }
    }

    /**
     * Letters in CIGAR string M (match/mismatch), I (insertion), D (deletion)
     * or NOTHING (not initialized)
     */
    enum CIGARstate { M, I, D, NOTHING };

    /**
     * Check if a certain row contains the final column
     * @param i The row index
     * @return True if the final column is contained in the row
     */
    bool inFinalColumn(const length_t i) const {
        return i >= getNumberOfRows() - getSizeOfFinalColumn();
    }

    /**
     * Find the CIGAR string of the alignment of a reference substring to
     * the query sequence the matrix was initialized with. The sequence
     * should be set before calling this function
     * @param ref the reference string that was aligned
     * @param score the alignment score between ref and query
     * @param CIGAR (output) the CIGAR string of the alignment
     */
    void findCIGAR(const Substring& ref, const uint score,
                   std::vector<std::pair<char, uint>>& CIGAR) {

        CIGAR.clear();
        CIGAR.reserve(2 * score + 1);
        // initialize the matrix with the alignment score
        initializeMatrix(score);

        // compute the rows
        for (unsigned int i = 0; i < ref.size(); i++) {
            computeRow(i + 1, ref[i]);
        }

        // trackback starting from the final cell
        uint i = ref.size();
        uint j = n - 1;
        assert(operator()(i, j) == score);

        CIGARstate state = NOTHING;

        while (j > 0 || i > 0) {

            const uint b = i / BLOCK_SIZE; // block identifier
            uint64_t bit = 1ull << ((j - b * BLOCK_SIZE) + DIAG_R0);

            if ((j > 0) && bv[i].HP & bit) { // gap in horizontal
                j--;
                if (state != CIGARstate::I) {
                    CIGAR.emplace_back('I', 0);
                    state = CIGARstate::I;
                }

            } else {
                const uint64_t& M = mv[b][char2idx[ref[i - 1]]];
                if ((i > 0 && j > 0) && ((M | ~bv[i].D0) & bit)) { // diagonal
                    i--;
                    j--;
                    if (state != CIGARstate::M) {
                        CIGAR.emplace_back('M', 0);
                        state = CIGARstate::M;
                    }

                } else { // gap in vertical
                    i--;
                    if (state != CIGARstate::D) {
                        CIGAR.emplace_back('D', 0);
                        state = CIGARstate::D;
                    }
                }
            }

            CIGAR.back().second++;
        }

        // reverse the cigar as it was created from end to beginning
        std::reverse(CIGAR.begin(), CIGAR.end());
    }

    /**
     *
     * Find cluster centers in the final column of the matrix. These centers
     * correspond to alignments that cannot be created from anohter alignment by
     * adding/removing gaps at the end.
     * @param lastRow the last row of the matrix that was filled in
     * @param refEnds (output) the end positions in the reference subsequence
     * that correspond to cluster centers.
     * @param maxED the maximal allowed edit distance
     * @param minED the minimal allowed edit distance
     */
    void findClusterCenters(const uint lastRow, std::vector<length_t>& refEnds,
                            uint maxED, uint minED) {
        refEnds.clear();
        refEnds.reserve(getSizeOfFinalColumn());

        uint firstRow = (m - 1) - getSizeOfFinalColumn();

        uint col = n - 1;

        for (uint i = lastRow; i > (m - 1) - getSizeOfFinalColumn(); i--) {
            uint ED = operator()(i, col);
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
     * Do backtracking and compute CIGAR string
     * @param ref reference sequence, pattern P should be set using
     * setSequence(P)(...)$
     * @param refEnd End offset of the reference sequence
     * @param refBegin Begin offset of the reference sequence (output)
     * @param ED Edit distance score associated with this alignment (output)
     * @param CIGAR CIGAR string (output)
     */
    void trackBack(const Substring& ref, const length_t refEnd,
                   length_t& refBegin, length_t& ED,
                   std::vector<std::pair<char, uint>>& CIGAR) const {

        CIGAR.clear();
        CIGAR.reserve(2 * MAX_ED + 1);

        uint64_t i = refEnd;
        uint64_t j = n - 1;
        ED = operator()(i, j);

        CIGARstate state = NOTHING;

        while (j > 0) {

            const uint b = i / BLOCK_SIZE; // block identifier
            /*  const uint64_t& M = mv[b][char2idx[ref[i - 1]]]; */
            uint64_t bit = 1ull << ((j - b * BLOCK_SIZE) + DIAG_R0);

            if ((bv[i].HP & bit)) { // gap in horizontal direction -> insertion
                j--;
                if (state != CIGARstate::I) {
                    CIGAR.emplace_back('I', 0);
                    state = CIGARstate::I;
                }

            } else if ((i > 0) && ((mv[b][char2idx[ref[i - 1]]] | ~bv[i].D0) &
                                   bit)) { // diagonal
                i--;
                j--;
                if (state != CIGARstate::M) {
                    CIGAR.emplace_back('M', 0);
                    state = CIGARstate::M;
                }

            } else { // gap in vertical direction
                i--;
                if (state != CIGARstate::D) {
                    CIGAR.emplace_back('D', 0);
                    state = CIGARstate::D;
                }
            }

            CIGAR.back().second++;
        }

        std::reverse(CIGAR.begin(), CIGAR.end());
        refBegin = i;
    }

    /**
     * Operator () overloading -- this procedure is O(1)
     * @param i Row index
     * @param j Column index
     * @return Score at position (i, j)
     */
    uint operator()(uint i, uint j) const {
        // make sure i and j are within matrix bounds
        assert(i < m);
        assert(j < n);

        // we need the bits in the range [b,e[ in HN and HP
        const uint bit = (i % BLOCK_SIZE) + DIAG_R0;
        uint b = (i > j) ? bit - (i - j) + 1 : bit + 1;
        uint e = (i > j) ? bit + 1 : bit + (j - i) + 1;

        uint64_t mask = ((1ull << (e - b)) - 1ull) << b;
        int negatives = __builtin_popcountll(bv[i].HN & mask);
        int positives = __builtin_popcountll(bv[i].HP & mask);

        uint score = bv[i].score;
        score += (i > j) ? (negatives - positives) : (positives - negatives);
        return score;
    }

    /**
     * Check whether after row i, the alignment involves only vertical gaps.
     * This happens when row i includes the final column n and when all
     * values on row i decrease monotonically
     * @return true of false
     */
    bool onlyVerticalGapsLeft(uint i) const {
        assert(i < m);

        if (i + LEFT < n) // if the column n is not yet reached on row i
            return false;

        const uint b = i / BLOCK_SIZE;
        const uint r = i % BLOCK_SIZE;

        // check if all relevant bits for HN are set to 1
        size_t bb = DIAG_R0 - Wv + r + 1;
        size_t be = DIAG_R0 + n - b * BLOCK_SIZE;
        return (((~bv[i].HN >> bb) << bb) << (64 - be)) == 0ull;
    }

    /**
     * Retrieves the first column index that is in the band for the
     * row
     * @param i The row
     * @returns The first column of the row
     */
    uint getFirstColumn(uint i) const {
        return (i <= Wv) ? 0u : i - Wv;
    }

    /**
     * Retrieves the last column index that needs to be filled in for the
     * row
     * @param i The row to fill in
     * @returns The last column to fill in
     */
    uint getLastColumn(uint i) const {
        return std::min(n - 1, i + Wh);
    }

    /**
     * Get the number of columns in the matrix
     * @return The number of columns in the matrix (== X.size() + 1)
     */
    uint getNumberOfCols() const {
        return n;
    }

    /**
     * Get the number of rows in the matrix
     * @return The number of rows in the matrix (== Y.size() + 1)
     */
    uint getNumberOfRows() const {
        return m;
    }

    /**
     * Check whether setSequenceX() has been called
     * @return True or false
     */
    bool sequenceSet() const {
        return !mv.empty();
    }

    void reset() {
        mv.clear();
    }

    /**
     * Get the vertical size of the final column
     * @return  the vertical size of the final column
     */
    uint getSizeOfFinalColumn() const {
        return Wh + Wv + 1;
    }

    /**
     * Print the banded matrix
     * @param maxRow Last row index to print
     */
    void printMatrix(uint maxRow = 500) const {
        for (uint i = 0; i < std::min<uint>(maxRow + 1, m); i++) {

            uint firstCol = getFirstColumn(i);
            uint lastCol = getLastColumn(i);
            std::cout << (i < 10 ? "0" : "") << std::to_string(i);
            std::cout << " [" << getFirstColumn(i) << "," << getLastColumn(i)
                      << "]\t";
            for (uint j = 0; j < firstCol; j++)
                std::cout << "  ";
            for (uint j = firstCol; j <= lastCol; j++)
                std::cout << operator()(i, j) << " ";
            std::cout << "\tRAC:" << -(int)Wv + (int)i << "/"
                      << std::log2(bv[i].RAC) - DIAG_R0;
            std::cout << (onlyVerticalGapsLeft(i) ? " - true" : " - false");
            uint minScore, minJ;
            findMinimumAtRow(i, minJ, minScore);
            std::cout << "  Min: " << minScore << "@" << minJ;
            std::cout << " FC: " << (inFinalColumn(i) ? " true" : " false");
            std::cout << std::endl;
        }
    }

  private:
    const static std::vector<char>
        char2idx; // static dictionary where characters are mapped to index

    uint maxED; // maximum allowed edit distance
    uint m;     // number of rows
    uint n;     // number of columns
    uint Wv;    // vertical width of the band
    uint Wh;    // horizontal width of the band

    std::vector<BitVectors> bv;              // bit vectors
    std::vector<std::array<uint64_t, 4>> mv; // match vectors
};

#endif
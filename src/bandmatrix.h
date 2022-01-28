/******************************************************************************
 *  Columba: Approximate Pattern Matching using Search Schemes                *
 *  Copyright (C) 2020-2022 - Luca Renders <luca.renders@ugent.be> and        *
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
#ifndef BANDMATRIX_H
#define BANDMATRIX_H

#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>

#include "substring.h"

typedef uint32_t length_t;

// ============================================================================
// CLASS BANDED MATRIX
// ============================================================================

// The band matrix class can best be understood as m horizontal bands each
// with a width of 2W+1 elements. It is allocated as a single m*(2W+1) array.
// For example, for W = 2 and m = 6.
// XX|XXX...|
//  X|XXXX..|
//   |XXXXX.|
//   |.XXXXX|
//   |..XXXX|X
//   |...XXX|XX
// The actual matrix is between |.| The Xs left and right of the |.| are
// allocated but should in principle not be addressed.
// The storage order is row-major.

class BandMatrix {
  private:
    std::vector<length_t> matrix;
    length_t W; // The off diagonal width of this bandmatrix
    length_t m; // number of rows
    length_t n; // number of columns
    const static int Wprod = 16;

    int finalCellFirstCol; // the final row of the zeroth column
    int rowsPerColumn; // the number of rows per column (starting from the Wth
                       // column)

    void initializeMatrix(const std::vector<int>& eds, const int& increase) {

        // initialize the first column
        length_t row = 0;
        for (auto rIt = eds.begin(); rIt != eds.end(); ++rIt) {
            operator()(row++, 0) = *rIt + increase;
        }
        for (; operator()(row - 1, 0) <= eds[0] + increase + W; row++) {
            operator()(row, 0) = operator()(row - 1, 0) + 1;
        }

        // initialize the first row
        for (length_t i = 1; i <= W; i++) {
            operator()(0, i) = i + eds[0] + increase;
        }
        // initialize the cells to the right of band
        for (length_t c = W + 1; c < m; c++) {
            operator()(c - (W + 1), c) = W + eds[0] + increase + 1;
        }

        // initialize the cells to the left of band
        for (length_t r = 1; r + row <= m; r++) {
            operator()(r + row - 1, r) = W + eds[0] + increase + 1;
        }
        finalCellFirstCol = row - 2;
    }

    void initializeMatrix(length_t startValue) {
        finalCellFirstCol = W;
        for (length_t i = 0; i <= W + 1; i++) {
            operator()(0, i) = i + startValue;
            operator()(i, 0) = i + startValue;
        }
        // set max elements at sides
        // first the elements on rows [1, W]
        for (length_t i = 1; i <= W; i++) {
            // right of band
            operator()(i, i + W + 1) = W + 1 + startValue;
        }

        // then the elements on rows [W + 1, x]

        for (length_t i = W + 1; i + W + 1 < n; i++) {
            // right of band
            operator()(i, i + W + 1) = W + 1 + startValue;
            // left of band
            operator()(i, i - (W + 1)) = W + 1 + startValue;
        }

        for (length_t i = std::max<int>((int)n - (W + 1), W + 1); i < m; i++) {
            // left of band
            operator()(i, i - (W + 1)) = W + 1 + startValue;
        }
    }

  public:
    /**
     * Constructor
     * @param pieceSize, the size of the piece to match, this will initialize
     * the top row
     * @param W, the width needed for this matrix
     * @param startValue, the value found at the startmatch, this should be the
     * minimum value of the vector eds if this vector is not empty
     * @param eds, a vector to initialize the first column, if an empty vector
     * is provided the first column will be initialized starting from startvalue
     * and up to W
     */
    BandMatrix(length_t pieceSize, int W, int startValue,
               const std::vector<int>& eds)
        : W(W), n(pieceSize + 1) {

        if (eds.empty()) {
            m = pieceSize + W + 1;
            matrix.resize(m * Wprod);
            rowsPerColumn = 2 * W + 1;
            initializeMatrix(startValue);
            return;
        }
        // minimum value should be put at startmatch...

        m = pieceSize + eds.size() + (W + eds[0] - eds.back());
        matrix.resize(m * Wprod);
        rowsPerColumn = m - pieceSize + W;

        initializeMatrix(eds, startValue);
    }

    /**
     * Constructor
     * @param pieceSize, the size of the piece to match, this will initialize
     * the top row
     * @param W, the width needed for this matrix
     * @param startValue, the value found at the startmatch, this will be put at
     * the origin
     */
    BandMatrix(length_t pieceSize, int W, int startValue)
        : W(W), m(pieceSize + W + 1), n(pieceSize + 1) {
        matrix.resize(m * Wprod);
        rowsPerColumn = 2 * W + 1;
        initializeMatrix(startValue);
        return;
    }

    /**
     * Constructor
     * @param m Number of rows
     * @param W Number of off-diagonal elements (one sided)
     */
    BandMatrix(length_t m, int W) : W(W), m(m) {
        matrix.resize(m * Wprod);
        n = m - W;
        initializeMatrix(0);
    }

    /**
     * Operator () overloading
     * @param i Row index
     * @param j Column index
     * @return Element at position (i, j)
     */
    length_t operator()(length_t i, int j) const {
        return matrix[i * Wprod + j - i + W];
    }

    /**
     * Operator () overloading
     * @param i Row index
     * @param j Column index
     * @return Reference to element at position (i, j)
     */
    length_t& operator()(length_t i, int j) {
        return matrix[i * Wprod + j - i + W];
    }

    /**
     * Get the band width
     * @return The band width
     */
    const length_t getWidth() const {
        return W;
    }

    void printMatrix(length_t maxRow = 500) const {
        length_t mRow = std::min<length_t>(maxRow + 1, m);
        for (length_t i = 0; i < mRow; i++) {
            length_t firstCol = 0;
            length_t lastCol = n - 1;
            std::string rowNumber = ((i < 10) ? "0" : "") + std::to_string(i);
            std::string row = "row " + rowNumber + ": ";

            for (length_t j = firstCol; j <= lastCol; j++) {
                int number = operator()(i, j);
                row += std::to_string(number) + " ";
            }
            for (length_t j = lastCol + 1; j < 2 * W + 1; j++) {
                row += ".\t";
            }
            std::cout << row << std::endl;
        }
        std::cout << "----------------------------------------------\n";
    }

    /**
     * Update the matrix by calculating the element at position row, column.
     * @param match whether the character was a match
     * @param row the row of the element to update
     * @param column the column of the element to update
     * @returns the new value at row, column
     */
    length_t updateMatrix(bool notMatch, unsigned int row,
                          unsigned int column) {
        length_t diag = operator()(row - 1, column - 1) + notMatch;
        length_t gapX = operator()(row, column - 1) + 1;
        length_t gapY = operator()(row - 1, column) + 1;

        length_t returnValue =
            std::min<length_t>(diag, std::min<length_t>(gapX, gapY));

        operator()(row, column) = returnValue;
        return returnValue;
    }

    /**
     * Retrieves the first column that needs to be filled in for the row
     * @param row the row to fill in
     * @returns the first column to fill in
     */
    const int getFirstColumn(int row) const {
        // leftmost cell of band
        return std::max(1, row - finalCellFirstCol);
    }
    /**
     * Retrieves the last column that needs to be filled in for the row
     * @param row the row to fill in
     * @returns the last column to fill in
     */
    const int getLastColumn(int row) const {
        // rightmost cell of band
        return std::min(n - 1, W + row);
    }

    const int getNumberOfRows() const {
        return m - 1;
    }

    const int getSizeOfFinalColumn() const {
        if (n > W) {
            return rowsPerColumn;
        }
        return finalCellFirstCol + n;
    }

    const int getLastColumn() const {
        return n - 1;
    }
};

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

typedef struct {
    uint64_t A; // match vector to indicate occurrences of A in X
    uint64_t C; // match vector to indicate occurrences of C in X
    uint64_t G; // match vector to indicate occurrences of G in X
    uint64_t T; // match vector to indicate occurrences of T in X
} MatchVectors;

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
        // create the alphabet mappingS
        char2idx = std::vector<char>(256, 4);
        char2idx['A'] = 0;
        char2idx['C'] = 1;
        char2idx['G'] = 2;
        char2idx['T'] = 3;
    }

    /**
     * Bit-encode the horizontal sequence X. Call this routine BEFORE calling
     * initializeMatrix(). You may call initializeMatrix() multiple times
     * (with different initialization settings) with a fixed sequence X.
     */
    void setSequence(const Substring& X) {
        n = X.size() + 1;   // number of columns
        m = 2 * MAX_ED + n; // this is an upper bound, the exact maxED
                            // is specified during initializeMatrix()

        // allocate and initialize the match vectors
        mv.resize((m + BLOCK_SIZE - 1) / BLOCK_SIZE);

        // encode the first block
        const uint64_t init = (1ull << LEFT) - 1;
        // first left bits are set to 1 for each characters, so that
        // initialization vector can propagate to first actual column
        mv[0] = {init, init, init, init};
        uint64_t bitmask = 1ull << LEFT;
        size_t je = std::min<size_t>(X.size(), WORD_SIZE - LEFT);
        for (size_t j = 0; j < je; j++) {
            assert(char2idx[X[j]] < 4); // assert ACTG alphabet
            mv[0][char2idx[X[j]]] |= bitmask;
            bitmask <<= 1;
        }

        // encode the remaining blocks
        for (size_t b = 1; b < mv.size(); b++) {
            // first blocksize bits of block b equal last blocksize bits of
            // block b - 1
            mv[b][0] = mv[b - 1][0] >> BLOCK_SIZE;
            mv[b][1] = mv[b - 1][1] >> BLOCK_SIZE;
            mv[b][2] = mv[b - 1][2] >> BLOCK_SIZE;
            mv[b][3] = mv[b - 1][3] >> BLOCK_SIZE;

            bitmask = 1ull << (WORD_SIZE - BLOCK_SIZE);
            size_t jb = WORD_SIZE - LEFT + (b - 1) * BLOCK_SIZE;
            size_t je = std::min<size_t>(X.size(), jb + BLOCK_SIZE);
            for (size_t j = jb; j < je; j++) {
                assert(char2idx[X[j]] < 4); // assert ACTG alphabet
                mv[b][char2idx[X[j]]] |= bitmask;
                bitmask <<= 1;
            }
        }
    }

    /**
     * Initialize the alignment matrix
     * @param maxED Maximum edit distance allowed during alignment
     * @param initED Edit distances of column zero (default = 0, 1, ... maxED)
     */
    void initializeMatrix(uint maxED, const std::vector<uint>& initED = {}) {
        // make sure maxED is within supported range
        assert(maxED <= MAX_ED);

        // sanity check on the initED vector
        assert(initED.empty() || initED.front() <= maxED);
        assert(initED.empty() || initED.back() <= maxED);

        this->maxED = maxED;            // store the maximum ED
        Wv = (initED.empty()) ? maxED : // vertical width of the band
                 initED.size() - 1 + maxED - initED.back();

        // sanity check on the size of initED
        assert(Wv <= 2 * MAX_ED);

        m = Wv + n; // number of rows

        // allocate and initialize bit vectors
        bv.resize(m);
        bv[0].score = initED.empty() ? 0 : initED[0];
        Wh = maxED - bv[0].score; // horizontal width of the band

        // initialize top row as [2*MAX_ED, ..., 2, 1, 0, 1, 2, ...]
        // decrease in first LEFT bits and increase in remaining bits
        bv[0].HP = (~0ull) << LEFT;
        bv[0].HN = ~bv[0].HP;

        // correct top row if initED has been specified
        for (size_t i = 1; i < std::min<size_t>(initED.size(), LEFT + 1); i++) {
            if (initED[i] < initED[i - 1]) {
                bv[0].HP ^= 1ull << (LEFT - i); // set HP to 1
                bv[0].HN ^= 1ull << (LEFT - i); // set HN to 0
            } else if (initED[i] == initED[i - 1]) {
                bv[0].HN ^= 1ull << (LEFT - i); // set HN to 0
            }
        }

        // RAC equals the right-most active element
        bv[0].RAC = 1ull << (DIAG_R0 + Wh);
    }

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

    void
    findLocalMinimaRow(uint i, uint maxED,
                       std::vector<std::pair<uint, uint>>& posAndScore) const {

        uint jMin = getFirstColumn(i);
        uint jMax = getLastColumn(i);

        for (uint j = jMin; j <= jMax; j++) {
            uint score = operator()(i, j);

            if (score <= maxED &&
                (j == jMin || score <= operator()(i, j - 1)) &&
                (j == jMax || score <= operator()(i, j + 1))) {
                posAndScore.emplace_back(j, score);
            }
        }
    }

    /**
     * Find the minimum edit distance value and its position on a row
     * Find the minimum edit distance value and its position in a row
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

    enum CIGARstate { M, I, D, NOTHING };

    /**
     * Check if a certain row contains the final column (column n+1)
     * @param i The row index
     * @return True of false
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
            const uint64_t& M = mv[b][char2idx[ref[i - 1]]];
            uint64_t bit = 1ull << ((j - b * BLOCK_SIZE) + DIAG_R0);

            if ((j > 0) && bv[i].HP & bit) { // gap in horizontal
                j--;
                if (state != CIGARstate::I) {
                    CIGAR.emplace_back('I', 0);
                    state = CIGARstate::I;
                }
            } else if ((i > 0 && j > 0) &&
                       ((M | ~bv[i].D0) & bit)) { // diagonal
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

            CIGAR.back().second++;
        }
        std::reverse(CIGAR.begin(), CIGAR.end());
    }

    void findClusterCenters(const uint lastRow, std::vector<uint>& refEnds,
                            uint maxED, uint minED) {

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
    void trackBack(const Substring& ref, const uint refEnd, uint& refBegin,
                   uint& ED, std::vector<std::pair<char, uint>>& CIGAR) const {

        CIGAR.clear();
        CIGAR.reserve(2 * MAX_ED + 1);

        uint i = refEnd;
        uint j = n - 1;
        ED = operator()(i, j);

        CIGARstate state = NOTHING;

        while (j > 0) {
            const uint b = i / BLOCK_SIZE; // block identifier
            const uint64_t& M = mv[b][char2idx[ref[i - 1]]];
            uint64_t bit = 1ull << ((j - b * BLOCK_SIZE) + DIAG_R0);

            if (bv[i].HP & bit) { // gap in horizontal direction -> insertion
                j--;
                if (state != CIGARstate::I) {
                    CIGAR.emplace_back('I', 0);
                    state = CIGARstate::I;
                }

            } else if ((i > 0) && ((M | ~bv[i].D0) & bit)) { // diagonal
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
    std::vector<char> char2idx;

    uint maxED; // maximum allowed edit distance
    uint m;     // number of rows
    uint n;     // number of columns
    uint Wv;    // vertical width of the band
    uint Wh;    // horizontal width of the band

    std::vector<BitVectors> bv;              // bit vectors
    std::vector<std::array<uint64_t, 4>> mv; // match vectors
};

#endif

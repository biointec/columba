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

#include "bitparallelmatrix.h"
#include "logger.h"    // for Logger, logger
#include "substring.h" // for Substring

#include <cmath>   // for log2
#include <cstdint> // for uint64_t
#include <sstream> // for stringstream
#include <string>  // for char_traits, operator<<, to_string

// Make the static dictionary that maps characters to indexes
const std::vector<char> IBitParallelED::char2idx =
    IBitParallelED::createChar2idx();

template <typename WordType>
void BitParallelED<WordType>::setSequence(const Substring& X) {

    n = X.size() + 1;          // number of columns
    m = 2 * MATRIX_MAX_ED + n; // this is an upper bound, the exact maxED
                               // is specified during initializeMatrix()

    // allocate and initialize the match vectors
    mv.resize((m + BLOCK_SIZE - 1) / BLOCK_SIZE);

    // encode the first block
    const WordType init = (WordType(1) << LEFT) - WordType(1);
    // first left bits are set to 1 for each characters, so that
    // initialization vector can propagate to first actual column
    std::fill(mv[0].begin(), mv[0].end(), init);
    WordType bitmask = WordType(1) << LEFT;

    size_t je = std::min<size_t>(X.size(), WORD_SIZE - LEFT);
    for (size_t j = 0; j < je; j++) {
        assert(char2idx[X[j]] < N_MATCH_VECTORS); // assert match vector exists
        mv[0][char2idx[X[j]]] |= bitmask;
        bitmask <<= 1;
    }

    // encode the remaining blocks
    for (size_t b = 1; b < mv.size(); b++) {
        // first blocksize bits of block b equal last blocksize bits of
        // block b - 1
        for (size_t i = 0; i < N_MATCH_VECTORS; i++) {
            mv[b][i] = mv[b - 1][i] >> BLOCK_SIZE;
        }

        bitmask = WordType(1) << (WORD_SIZE - BLOCK_SIZE);
        size_t jb_b = WORD_SIZE - LEFT + (b - 1) * BLOCK_SIZE;
        size_t je_b = std::min<size_t>(X.size(), jb_b + BLOCK_SIZE);
        for (size_t j = jb_b; j < je_b; j++) {
            assert(char2idx[X[j]] < N_MATCH_VECTORS); // assert ACTG alphabet
            mv[b][char2idx[X[j]]] |= bitmask;
            bitmask <<= 1;
        }
    }
}

template <typename WordType>
void BitParallelED<WordType>::initializeMatrix(
    uint32_t maxED, const std::vector<uint32_t>& initED) {
    // make sure maxED is within supported range
    assert(maxED <= MATRIX_MAX_ED);
    assert(sequenceSet());

    // sanity check on the initED vector
    assert(initED.empty() || initED.front() <= maxED);
    assert(initED.empty() || initED.back() <= maxED);

    this->maxED = maxED;            // store the maximum ED
    Wv = (initED.empty()) ? maxED : // vertical width of the band
             initED.size() - 1 + maxED - initED.back();

    // sanity check on the size of initED
    assert(Wv <= 2 * MATRIX_MAX_ED);

    m = Wv + n; // number of rows

    // allocate and initialize bitvectors
    bv.resize(m);
    bv[0].score = initED.empty() ? 0 : initED[0];
    Wh = maxED - bv[0].score; // horizontal width of the band
    if (Wv + Wh + 1 > m) {
        m = Wv + Wh + 1;
        bv.resize(m);
    }

    // initialize top row as [2*MATRIX_MAX_ED, ..., 2, 1, 0, 1, 2, ...]
    // decrease in first LEFT bits and increase in remaining bits
    bv[0].HP = (~WordType(0)) << LEFT;
    bv[0].HN = ~bv[0].HP;

    // correct top row if initED has been specified
    const size_t n = std::min<size_t>(initED.size(), LEFT + 1);
    for (uint32_t i = 1; i < n; ++i) {
        if (initED[i] < initED[i - 1]) {
            bv[0].HP ^= WordType(1) << (LEFT - i); // set HP to 1
            bv[0].HN ^= WordType(1) << (LEFT - i); // set HN to 0
        } else if (initED[i] == initED[i - 1]) {
            bv[0].HN ^= WordType(1) << (LEFT - i); // set HN to 0
        }
    }

    // RAC equals the right-most active element
    bv[0].RAC = WordType(1) << (DIAG_R0 + Wh);
}

template <typename WordType>
void BitParallelED<WordType>::printMatrix(uint32_t maxRow) const {
    for (uint32_t i = 0; i < std::min<uint32_t>(maxRow + 1, m); i++) {

        uint32_t firstCol = getFirstColumn(i);
        uint32_t lastCol = getLastColumn(i);
        std::stringstream s;
        s << (i < 10 ? "0" : "") << std::to_string(i);
        s << " [" << getFirstColumn(i) << "," << getLastColumn(i) << "]\t";
        for (uint32_t j = 0; j < firstCol; j++)
            s << "  ";
        for (uint32_t j = firstCol; j <= lastCol; j++)
            s << operator()(i, j) << " ";
        s << "\tRAC:" << -(int)Wv + (int)i << "/"
          << std::log2(bv[i].RAC) - DIAG_R0;
        s << (onlyVerticalGapsLeft(i) ? " - true" : " - false");
        uint32_t minScore, minJ;
        findMinimumAtRow(i, minJ, minScore);
        s << "  Min: " << minScore << "@" << minJ;
        s << " FC: " << (inFinalColumn(i) ? " true" : " false");
        logger.logDeveloper(s);
    }
}

template class BitParallelED<uint64_t>;
template class BitParallelED<UInt128>;

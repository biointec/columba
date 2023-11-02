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
#include "bitparallelmatrix.h"

// Make the static dictionary that maps characters to indexes
const std::vector<char> BitParallelED::char2idx =
    BitParallelED::createChar2idx();

void BitParallelED::setSequence(const Substring& X) {
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

void BitParallelED::initializeMatrix(uint maxED,
                                     const std::vector<uint>& initED) {
    // make sure maxED is within supported range
    assert(maxED <= MAX_ED);
    assert(sequenceSet());

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

/******************************************************************************
 *  Columba 1.1: Approximate Pattern Matching using Search Schemes            *
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

#ifndef BWTREPR_H
#define BWTREPR_H

#include <cstdlib>
#include <vector>

#include "alphabet.h"
#include "bitvec.h"

// ============================================================================
// CLASS BWT REPRESENTATION (supports occ(c,k) and cumOcc(c,k) in O(1) time)
// ============================================================================

template <size_t S> // S is the size of the alphabet (including '$')
class BWTRepr {     // e.g. S = 5 for DNA (A,C,G,T + $)

  private:
    // The '$' character (cIdx == 0) is not encoded in the bitvector.
    // Hence, we use only S-1 bitvectors.

    BitvecIntl<S - 1> bv; // bitvector representation of the BWT
    size_t dollarPos;     // position of the dollar sign

  public:
    /**
     * Default constructor
     */
    BWTRepr() {
    }

    /**
     * Constructor
     * @param sigma Alphabet
     * @param BWT Burrows-Wheeler transformation
     */
    BWTRepr(const Alphabet<S>& sigma, const std::string& BWT)
        : bv(BWT.size() + 1), dollarPos(BWT.size()) {
        // The $-character (cIdx == 0) is not encoded in the bitvector.
        // Hence, use index cIdx-1 in the bitvector.

        for (size_t i = 0; i < BWT.size(); i++) {
            if (BWT[i] == '$') {
                dollarPos = i;
                continue;
            }

            for (size_t cIdx = sigma.c2i(BWT[i]); cIdx < S; cIdx++)
                bv(cIdx - 1, i) = true;
        }

        bv.index();
    }

    /**
     * Get occurrence count of character c in the range BWT[0...k[
     * @param cIdx Character index
     * @param k index
     * @return occ(c, k)
     */
    size_t occ(int cIdx, size_t k) const {
        // The $-character (cIdx == 0) is not encoded in the bitvector.
        // Hence, use index cIdx-1 in the bitvector.

        if (cIdx == 0) // special case for $-character
            return (k <= dollarPos) ? 0 : 1;

        return (cIdx == 1) ? bv.rank(cIdx - 1, k)
                           : bv.rank(cIdx - 1, k) - bv.rank(cIdx - 2, k);
    }

    /**
     * Get cumulative occurrence count of characters SMALLER than c
     * in the range BWT[0...k[
     * @param cIdx Character index
     * @param k index
     * @return cumOcc(c, k)
     */
    size_t cumOcc(int cIdx, size_t k) const {
        // The $-character (cIdx == 0) is not encoded in the bitvector.
        // Hence, use index cIdx-1 in the bitvector.

        if (cIdx == 0) // special case for $-character
            return 0;

        return (cIdx == 1) ? ((k <= dollarPos) ? 0 : 1)
                           : bv.rank(cIdx - 2, k) + ((k <= dollarPos) ? 0 : 1);
    }

    /**
     * Write table to disk
     * @param filename File name
     */
    void write(const std::string& filename) {
        std::ofstream ofs(filename);
        if (!ofs)
            throw std::runtime_error("Cannot open file: " + filename);

        ofs.write((char*)&dollarPos, sizeof(dollarPos));
        bv.write(ofs);
    }

    /**
     * Load table from disk
     * @param filename File name
     */
    bool read(const std::string& filename) {
        std::ifstream ifs(filename);
        if (!ifs)
            return false;

        ifs.read((char*)&dollarPos, sizeof(dollarPos));
        bv.read(ifs);

        return true;
    }
};

#endif
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

#ifndef ALPHABET_H
#define ALPHABET_H

#include "assert.h"
#include <array>
#include <cstdlib>
#include <vector>

#define NUM_CHAR 256

// ============================================================================
// CLASS ALPHABET (convert ASCII value <-> character index)
// ============================================================================
typedef uint32_t length_t;
template <size_t S> // S is the size of the alphabet (including '$')
class Alphabet {    // e.g. S = 5 for DNA (A,C,G,T + $)

  private:
    std::vector<int> charToIndex;    // character map (ascii -> idx)
    std::array<char, S> indexToChar; // inverse map (idx -> ascii)
    std::array<char, S - 1> redChar;
    /**
     * Initialize charMap and invCharMap vector
     * @param charCounts Vector containing the counts for each character,
     * thus indicating the presence / absence of each character
     */
    void initialize(const std::vector<length_t>& charCounts) {
        charToIndex = std::vector<int>(NUM_CHAR, -1);
        for (size_t i = 0, j = 0; i < charCounts.size(); i++) {
            if (charCounts[i] > 0) {
                charToIndex[i] = j;
                indexToChar[j++] = char(i);
            }
        }
        for (size_t i = 1; i < S; i++) {
            redChar[i - 1] = indexToChar[i];
        }
    }

  public:
    /**
     * Default constructor
     */
    Alphabet() {
    }

    /**
     * Build an alphabet for a text T
     * @param T Input text
     */
    Alphabet(const std::string& T) {
        // count the number of characters in the text
        std::vector<size_t> charCounts(NUM_CHAR, 0);
        for (char c : T)
            charCounts[(unsigned char)c]++;
        initialize(charCounts);
    }

    /**
     * Build an alphabet given the character counts
     * @param charCounts Character counts
     */
    Alphabet(const std::vector<length_t>& charCounts) {
        initialize(charCounts);
    }

    /**
     * Convert a character to a character index
     * @param c character
     * @return character index
     */
    int c2i(char c) const {
        assert(inAlphabet(c));
        return charToIndex[(unsigned char)c];
    }

    /**
     * Check whether a character exists in the extended alphabet
     * @param c character
     * @return true of false
     */
    bool inAlphabet(char c) const {
        return charToIndex[(unsigned char)c] >= 0;
    }

    /**
     * Convert a character index to a character
     * @param cIdx character index
     * @return character
     */
    char i2c(int cIdx) const {
        assert(cIdx < (int)S);
        return indexToChar[cIdx];
    }

    /**
     * Operator() overloading
     * @return an array of non sentinel characters in the alphabet
     */
    std::array<char, S - 1> operator()() const {
        return redChar;
    }

    /**
     * Return the size of the alphabet
     * @return The size of the alphabet
     */
    size_t size() const {
        return S;
    }
};

#endif
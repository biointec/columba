/******************************************************************************
 *  Columba: Approximate Pattern Matching using Search Schemes                *
 *  Copyright (C) 2020 - Luca Renders <luca.renders@ugent.be> and             *
 *                       Jan Fostier <jan.fostier@ugent.be>                   *
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
#ifndef CUSTOMTYPEDEFS_H
#define CUSTOMTYPEDEFS_H

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream> // used for printing
#include <stdexcept>
#include <string>
#include <vector>

#include "bwtrepr.h"

enum Direction { FORWARD, BACKWARD };

// compute |a-b| in a safe manner
template <typename T> T abs_diff(T a, T b) {
    return a > b ? a - b : b - a;
}
// ============================================================================
// CLASS SUBSTRING
// ============================================================================
class Substring {
  private:
    const std::string* text; // pointer to the string this is a substring of
    unsigned int startIndex; // the startIndex of this substring in the text
    unsigned int
        endIndex; // the endIndex of this substring in the text (non-inclusive)
    Direction d;

  public:
    Substring(const std::string& t, unsigned int start, unsigned int end,
              Direction dir = FORWARD)
        : startIndex(start), endIndex(std::min<int>(end, t.size())), d(dir) {
        text = &t;
        d = dir;
    }

    Substring(const std::string& t, Direction dir = FORWARD)
        : startIndex(0), endIndex(t.size()), d(dir) {
        text = &t;
    }

    Substring(const std::string* t, unsigned int start, unsigned int end,
              Direction dir = FORWARD)
        : text(t), startIndex(start), endIndex(end), d(dir) {
    }
    Substring(const Substring* s, unsigned int start, unsigned int end)
        : startIndex(start), endIndex(end) {
        text = s->text;
        d = s->d;
    }

    Substring(const Substring& s, unsigned int start, unsigned int end)
        : startIndex(start), endIndex(end) {
        text = s.text;
        d = s.d;
    }

    void setDirection(Direction nd) {
        d = nd;
    }

    const Substring getSubPiece(unsigned int start) const {
        if (d == FORWARD) {
            return Substring(this, startIndex + start, endIndex);
        } else {
            return Substring(this, startIndex, endIndex - start);
        }
    }

    /**
     * Get the character at index i of this substring
     * @param i the index to get the character from
     * @returns the character at index i
     */
    char operator[](int i) const {
        return (d == FORWARD) ? text->at(startIndex + i)
                              : text->at(endIndex - 1 - i);
    }

    /**
     * Get the size of this substring
     * @returns the size of this substring
     */
    unsigned int size() const {
        if (empty()) {
            return 0;
        }
        return endIndex - startIndex;
    }

    /**
     * Get the length of this substring (equals the size)
     * @returns the length of this substring
     */
    unsigned int length() const {
        return size();
    }

    /**
     * Check if this substring is empty
     * @returns a bool that indicates wheter the substring was empty
     */
    bool empty() const {
        return endIndex <= startIndex;
    }

    /**
     * Get the end of this substring
     * @returns the end Index of this substring (non-inclusive)
     */
    unsigned int end() const {
        return endIndex;
    }

    /**
     * Get the begin of this substring
     * @returns the begin index of the substring
     */
    unsigned int begin() const {
        return startIndex;
    }

    Substring& operator=(const Substring& other) {
        this->text = other.text;
        this->startIndex = other.begin();
        this->endIndex = other.end();
        this->d = other.d;

        return *this;
    }

    std::string getSubstring() const {
        if (empty()) {
            return "";
        }
        return text->substr(startIndex, endIndex - startIndex);
    }

    void setEnd(unsigned int newEnd) {
        endIndex = newEnd;
    }

    void setBegin(unsigned int n) {
        startIndex = n;
    }

    void incrementEnd() {
        endIndex++;
    }

    void decrementBegin() {
        startIndex--;
    }
};

// ============================================================================
// (TYPE) DEFINITIONS AND PROTOTYPES
// ============================================================================

typedef uint32_t length_t;

// ============================================================================
// HELPERS FOR BWT
// ============================================================================

// ranges to be used over any array/string/vector..
// this is a halfopen interval [begin, end[

class Range {
  public:
    length_t begin;
    length_t end;

    bool empty() const {
        return end <= begin;
    }

    length_t width() const {
        return (empty()) ? 0 : end - begin;
    }

    Range(length_t b, length_t e) {
        end = e;
        begin = b;
    }

    Range() {
        end = 0;
        begin = 0;
    }

    bool operator==(const Range& o) const {
        return o.begin == begin && o.end == end;
    }
};

class SARangePair {
  public:
    Range rangeSA;
    Range rangeSARev;

    bool empty() const {
        return rangeSA.empty();
    }

    length_t width() const {
        return rangeSA.width();
    }

    SARangePair(Range rangeSA, Range rangeSARev)
        : rangeSA(rangeSA), rangeSARev(rangeSARev) {
    }

    SARangePair() {
        rangeSA = Range();
        rangeSARev = Range();
    }

    bool operator==(const SARangePair& o) const {
        // only the first range matters as the ranges imply each other
        return o.rangeSA == rangeSA;
    }
};

// approximate matches in the text or in the SA
typedef struct {
    Range range;  // the range of the match
    int editDist; // the edit distance of this match
} AppMatch;

bool operator<(const AppMatch& l, const AppMatch& r);

// approximate matches in the suffix array, the range of the match is a range in
// the suffix array
typedef struct {
    AppMatch match; // the match with range and ED
    length_t depth; // the depth of this match
} AppMatchSA;

/**
 * Helper function to make an AppMatchSA struct
 * @param range the range in the SA of this match
 * @param ED the edit distance of this match
 * @param depth the depthe of this match (length of matched substring)
 * @returns the AppMatch with these arguments
 */
AppMatchSA makeAppMatchSA(Range range, int ED, length_t depth);

/**
 * Get the string representation of a range
 * @param range the range to get the string representation of
 */
std::string tostring(Range range);

// ============================================================================
// HELPERS FOR BIDIRECTIONAL BWT
// ============================================================================

// A bidirectional approximate match very similar to AppMatchSA
typedef struct {
    SARangePair ranges; // The ranges of this match
    length_t editDist;  // the edit distance
    length_t depth;     // the depth

} BiAppMatchSA;

/**
 * Make a biderictional approximate match in the suffix array
 * @param ranges the ranges of this approximate match (range in SA and in SA')
 * @param ED the edit distance of this approximate match
 * @param depth the depth (=length) of this approximate match
 */
const BiAppMatchSA makeBiAppMatchSA(SARangePair ranges, int ED, length_t depth);

std::ostream& operator<<(std::ostream& o, const BiAppMatchSA& m);

AppMatchSA discardReverseRange(const BiAppMatchSA& biappmatch);

bool operator<(const BiAppMatchSA& lhs, const BiAppMatchSA& rhs);

bool operator==(const BiAppMatchSA& lhs, const BiAppMatchSA& rhs);

std::vector<AppMatchSA>
discardReverseRanges(const std::vector<BiAppMatchSA>& bivector);

// Define a Search
typedef struct {
    std::vector<int> lowerBounds; // the vector with lowerbounds
    std::vector<int> upperBounds; // the vector with upperbounds
    std::vector<int> order;       // the vector with the order of the parts

    std::vector<Direction> directions;

    std::vector<bool> directionSwitch;
} Search;

/**
 * Helper function to easily create a search.
 * @param order the order of the parts
 * @param lowerBounds the lowerbounds
 * @param upperBounds the upperbounds
 */
Search makeSearch(std::vector<int> order, std::vector<int> lowerBounds,
                  std::vector<int> upperBounds);

void setPartsDirections(const Search& s, std::vector<Substring>& parts);

class PrefixOccurences {
  private:
    std::vector<Bitvec> prefixOccurences;

  public:
    PrefixOccurences(){};

    PrefixOccurences(size_t numberOfCharacters,
                     const std::vector<int>& charToIndex,
                     const std::string& bwt) {
        prefixOccurences.resize(numberOfCharacters);

        // create the bitvector for each character
        for (size_t i = 0; i < numberOfCharacters; i++) {
            prefixOccurences[i] = Bitvec(bwt.size() + 1);
        }

        // for each entry in the bwt update the bitvector of all characters
        // greater or equal
        for (size_t i = 0; i < bwt.size(); i++) {
            for (size_t idx = charToIndex[(unsigned char)bwt[i]];
                 idx < numberOfCharacters; idx++) {
                prefixOccurences[idx][i] = true;
            }
        }

        // index all bitvectors
        for (auto& bv : prefixOccurences) {
            bv.index();
        }
    }

    /**
     * Function that returns the nummber of occurences before an index of the
     * symbol at symbolindex in the alphabet
     * @param symbolIndex the index of the the symbol in the alphabet to count
     * the occrences of at index index
     * @param index the index whose entry for symbol in the occurences table is
     * asked
     * @return the number of occurences of the symbol before index in the bwt
     */
    length_t getNumberOfOccurences(length_t symbolIndex, length_t index) const {
        return (symbolIndex == 0)
                   ? prefixOccurences[symbolIndex].rank(index)
                   : prefixOccurences[symbolIndex].rank(index) -
                         prefixOccurences[symbolIndex - 1].rank(index);
    }

    /**
     * Function that returns the nummber of occurences before the index of all
     * symbols smaller than the symbol at symbolindex in the  bwt
     * @param symbolIndex the index in the alphabet whose number of prefix
     * occurences is queried.
     * @param index the index whose entry for symbol in the prefixoccurences
     * table is asked
     * @return the number of occurences of symbols smaller than symbol at
     * symbolindex before index index in the bwt
     */
    length_t getNumberOfPrefixOccurences(length_t symbolIndex,
                                         length_t index) const {

        return (symbolIndex == 0)
                   ? 0
                   : prefixOccurences[symbolIndex - 1].rank(index);
    }

    void write(const std::string& filename) {
        std::ofstream ofs(filename);
        if (!ofs)
            throw std::runtime_error("Cannot open file: " + filename);

        size_t numVec = prefixOccurences.size();
        ofs.write((char*)&numVec, sizeof(size_t));
        for (const auto& it : prefixOccurences)
            it.write(ofs);
    }

    bool read(const std::string& filename) {
        std::ifstream ifs(filename);
        if (!ifs)
            return false;

        size_t numChar;
        ifs.read((char*)&numChar, sizeof(numChar));
        prefixOccurences.resize(numChar);
        for (auto& it : prefixOccurences)
            it.read(ifs);

        return true;
    }
};

// ============================================================================
// CLASS PREFOCC for 4 characters
// ============================================================================
class PrefixOccurrences4 {

  private:
    BitvecIntl<4> prefixOccurences4;
    size_t dollarPos;

  public:
    PrefixOccurrences4(){};

    PrefixOccurrences4(size_t numberOfCharacters,
                       const std::vector<int>& charToIndex,
                       const std::string& bwt) {
        if (numberOfCharacters != 5) {
            throw std::runtime_error("Interleaved prefix occurence can only be "
                                     "used for 5 letter alphabet");
        }

        // create the interleaved bitvector for each character
        // prefixOccurences4 = BitvecIntl(bwt.size() + 1);

        // for each entry in the bwt update the bitvector of all characters
        // greater or equal
        for (size_t i = 0; i < bwt.size(); i++) {
            if (bwt[i] == '$') {
                dollarPos = i;
                continue;
            }
            for (size_t idx = charToIndex[(unsigned char)bwt[i]];
                 idx < numberOfCharacters; idx++) {

                prefixOccurences4(idx - 1, i) = true;
            }
        }

        // index
        prefixOccurences4.index();
    }

    /**
     * Function that returns the nummber of occurences before an index of the
     * symbol at symbolindex in the alphabet
     * @param symbolIndex the index of the the symbol in the alphabet to count
     * the occrences of at index index
     * @param index the index whose entry for symbol in the occurences table is
     * asked
     * @return the number of occurences of the symbol before index in the bwt
     */
    length_t getNumberOfOccurences(length_t symbolIndex, length_t index) const {
        if (symbolIndex == 0) {
            return (index <= dollarPos) ? 0 : 1;
        }
        return (symbolIndex == 1)
                   ? prefixOccurences4.rank(symbolIndex - 1, index)
                   : prefixOccurences4.rank(symbolIndex - 1, index) -
                         prefixOccurences4.rank(symbolIndex - 2, index);
    }

    /**
     * Function that returns the nummber of occurences before the index of all
     * symbols smaller than the symbol at symbolindex in the  bwt
     * @param symbolIndex the index in the alphabet whose number of prefix
     * occurences is queried.
     * @param index the index whose entry for symbol in the prefixoccurences
     * table is asked
     * @return the number of occurences of symbols smaller than symbol at
     * symbolindex before index index in the bwt
     */
    length_t getNumberOfPrefixOccurences(length_t symbolIndex,
                                         length_t index) const {
        if (symbolIndex == 0) {
            return 0;
        }
        bool notAfterDollar = index <= dollarPos;

        if (symbolIndex == 1) {
            return (notAfterDollar) ? 0 : 1;
        }

        return (symbolIndex == 1)
                   ? (notAfterDollar) ? 0 : 1
                   : prefixOccurences4.rank(symbolIndex - 2, index) +
                         ((notAfterDollar) ? 0 : 1);
    }
    void write(const std::string& filename) {
        std::ofstream ofs(filename);
        if (!ofs)
            throw std::runtime_error("Cannot open file: " + filename);

        ofs.write((char*)&dollarPos, sizeof(dollarPos));
        prefixOccurences4.write(ofs);
    }

    bool read(const std::string& filename) {
        std::ifstream ifs(filename);
        if (!ifs)
            return false;

        ifs.read((char*)&dollarPos, sizeof(dollarPos));
        prefixOccurences4.read(ifs);

        return true;
    }
};
std::ifstream getStream(const std::string& file);
std::string readString(const std::string& s);
std::vector<length_t> readArray(size_t length, std::ifstream& ifs);
void readArray(const std::string& s, size_t length, std::vector<length_t>& sa);
void readArray2(const std::string& s, size_t length, std::vector<size_t>& sa);
#endif
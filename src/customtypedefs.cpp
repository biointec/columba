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
#include "customtypedefs.h"

// ============================================================================
// HELPERS FOR BWT
// ============================================================================

Range makeRange(length_t b, length_t e) {
    Range r;
    r.begin = b;
    r.end = e;
    return r;
}

bool operator<(const AppMatch& l, const AppMatch& r) {
    if (l.range.begin != r.range.begin) {
        return l.range.begin < r.range.begin;
    } else {
        // begin is equal, better ed is smarter
        if (l.editDist != r.editDist) {
            return l.editDist < r.editDist;
        } else {
            // shorter read is smaller...
            return l.range.width() < r.range.width();
        }
    }
}

AppMatchSA makeAppMatchSA(Range range, int ED, length_t depth) {
    AppMatch match;
    match.editDist = ED;
    match.range = range;
    AppMatchSA msa;
    msa.match = match;
    msa.depth = depth;

    return msa;
}

// ============================================================================
// HELPERS FOR BIDIRECTIONAL BWT
// ============================================================================

const BiAppMatchSA makeBiAppMatchSA(SARangePair ranges, int ED,
                                    length_t depth) {
    BiAppMatchSA m;
    m.ranges = ranges;
    m.editDist = ED;
    m.depth = depth;

    return m;
}

std::ostream& operator<<(std::ostream& o, const BiAppMatchSA& m) {
    return o << "SARange: " << tostring(m.ranges.rangeSA)
             << "\tEdit distance: " << m.editDist << "\tdepth: " << m.depth;
}

std::string tostring(Range range) {
    return "[" + std::to_string(range.begin) + "," + std::to_string(range.end) +
           ")";
}

AppMatchSA discardReverseRange(const BiAppMatchSA& biappmatch) {
    return makeAppMatchSA(biappmatch.ranges.rangeSA, biappmatch.editDist,
                          biappmatch.depth);
}

std::vector<AppMatchSA>
discardReverseRanges(const std::vector<BiAppMatchSA>& bivector) {
    std::vector<AppMatchSA> vector;
    vector.reserve(bivector.size());

    std::transform(bivector.begin(), bivector.end(), std::back_inserter(vector),
                   discardReverseRange);
    return vector;
}

bool operator<(const BiAppMatchSA& lhs, const BiAppMatchSA& rhs) {
    if (lhs.ranges.rangeSA.begin != rhs.ranges.rangeSA.begin) {
        return lhs.ranges.rangeSA.begin < rhs.ranges.rangeSA.begin;
    } else {
        // begin is equal, better ed is smarter
        if (lhs.editDist != rhs.editDist) {
            return lhs.editDist < rhs.editDist;
        } else {
            // shorter read is smaller...
            return lhs.ranges.width() < rhs.ranges.width();
        }
    }
}

bool operator==(const BiAppMatchSA& lhs, const BiAppMatchSA& rhs) {
    return lhs.ranges == rhs.ranges && lhs.editDist == rhs.editDist &&
           lhs.depth == rhs.depth;
}

void setPartsDirections(const Search& s, std::vector<Substring>& parts) {

    // set the directions for the parts
    for (length_t i = 0; i < s.order.size(); i++) {
        int pieceNumber = s.order[i];
        Direction d = s.directions[i];
        Substring& piece = parts[pieceNumber];
        piece.setDirection(d);
    }
}

Search makeSearch(std::vector<int> order, std::vector<int> lowerBounds,
                  std::vector<int> upperBounds) {
    if (order.size() != lowerBounds.size() ||
        order.size() != upperBounds.size()) {
        throw std::runtime_error(
            "Could not create search, the sizes of all vectors are not equal");
    }
    Search s;
    s.upperBounds.swap(upperBounds);
    s.lowerBounds.swap(lowerBounds);
    s.order.swap(order);

    s.directions.push_back((s.order[1] > s.order[0]) ? FORWARD : BACKWARD);

    for (length_t i = 1; i < s.order.size(); i++) {
        Direction d = (s.order[i] > s.order[i - 1]) ? FORWARD : BACKWARD;
        s.directions.push_back(d);
    }
    // first partition is not a swithc
    s.directionSwitch.push_back(false);

    // second partition is never a switch
    s.directionSwitch.push_back(false);

    for (length_t i = 2; i < s.directions.size(); i++) {
        s.directionSwitch.push_back(s.directions[i] != s.directions[i - 1]);
    }

    return s;
}

std::ifstream getStream(const std::string& file) {
    std::ifstream ifs(file);
    if (!ifs) {
        throw std::runtime_error("Cannot open file: " + file);
    }

    return ifs;
}

std::string readString(const std::string& s) {
    std::ifstream ifs = getStream(s);
    std::string ret;
    ifs.seekg(0, std::ios::end);
    ret.resize(ifs.tellg());
    ifs.seekg(0, std::ios::beg);
    ifs.read((char*)&ret[0], ret.size());
    ifs.close();
    return ret;
}

std::vector<length_t> readArray(size_t length, std::ifstream& ifs) {

    std::vector<length_t> sa;
    sa.reserve(length);

    length_t element;
    while (ifs >> element) {
        sa.push_back(element);
    }
    ifs.close();
    return sa;
}

void readArray2(const std::string& s, size_t length, std::vector<size_t>& sa) {

    std::ifstream ifs = getStream(s);
    ifs.seekg(0, std::ios::end);
    size_t elementsInBinary = ifs.tellg() / sizeof(size_t);
    ifs.seekg(0, std::ios::beg);
    if (elementsInBinary == length) {
        // truely a binary file

        sa.resize(elementsInBinary);

        ifs.read((char*)&sa[0], sa.size() * sizeof(size_t));
        ifs.close();

    } else {
        // text file
        // TODO
    }
}

void readArray(const std::string& s, size_t length, std::vector<length_t>& sa) {

    std::ifstream ifs = getStream(s);
    ifs.seekg(0, std::ios::end);
    size_t elementsInBinary = ifs.tellg() / sizeof(length_t);
    ifs.seekg(0, std::ios::beg);
    if (elementsInBinary == length) {
        // truely a binary file

        sa.resize(elementsInBinary);

        ifs.read((char*)&sa[0], sa.size() * sizeof(length_t));
        ifs.close();

    } else {
        // text file
        sa = readArray(length, ifs);
    }
}

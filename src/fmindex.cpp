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
#include "fmindex.h"

using namespace std;

// ============================================================================
// CLASS FMIndex
// ============================================================================

// ----------------------------------------------------------------------------
// PRIVATE ROUTINES FOR READING IN THE FILES
// ----------------------------------------------------------------------------

void FMIndex::fromFiles(const string& prefix) {
    cout << "Reading in files with prefix " << prefix << endl;
    // read the text
    cout << "Reading " << prefix << ".txt" << endl;
    string text = readString(prefix + ".txt");
    textLength =
        (text[text.size() - 1] == '\n') ? text.size() - 1 : text.size();

    cout << "Done reading text (size: " << text.size() << ")" << endl;

    // read the counts table
    counts.assign(256, 0);
    vector<size_t> countsFromDisk(256, 0);
    cout << "Reading " << prefix << ".cct"
         << "\n";

    readArray2(prefix + ".cct", 256, counts);
    sigma = Alphabet<ALPHABET>(counts);

    cout << "Text contains " << sigma.size() << " unique characters" << endl;

    for (int i = 0; i < 256; i++) {
        if (counts[i] > 0) {
            cout << (unsigned char)i << "\t";
        }
    }
    cout << endl;

    // make the counts cumulative
    length_t cumCount = 0; // cumulative character counts
    vector<size_t> cumCounts;
    cumCounts.reserve(sigma.size());
    for (size_t i = 0; i < counts.size(); i++) {
        if (counts[i] == 0)
            continue;
        cumCounts.push_back(cumCount);
        cumCount += counts[i];
    }

    counts.swap(cumCounts);

    // read the SA
    cout << "Reading " << prefix << ".sa." << sparseFactorSA << endl;
    string saFilename = prefix + ".sa." + to_string(sparseFactorSA);
    length_t saLength = (textLength + sparseFactorSA - 1) / sparseFactorSA;
    readArray(saFilename, saLength, sa);

    cout << "Done reading suffix array (sparseness factor: " << sparseFactorSA
         << ")" << endl;

    // read the BWT
    cout << "Reading " << prefix << ".bwt" << endl;
    bwt = readString(prefix + ".bwt");

    cout << "Done reading BWT (size: " << bwt.size() << ")" << endl;

    // read the prefix occurrence table
    cout << "Reading " << prefix << ".brt" << endl;

    if (!fwdRepr.read(prefix + ".brt"))
        throw runtime_error("Cannot open file: " + prefix + ".brt");
    cout << "Done reading prefix occurrence table" << endl;
}

// ----------------------------------------------------------------------------
// ROUTINES FOR PATTERN MATCHING
// ----------------------------------------------------------------------------
length_t FMIndex::getNumberOfPrefOcc(char symbol, length_t index) const {
    return fwdRepr.cumOcc(sigma.c2i(symbol), index);
}

length_t FMIndex::getNumberOfOcc(char symbol, length_t index) const {
    return fwdRepr.occ(sigma.c2i(symbol), index);
}

length_t FMIndex::findLF(length_t k) const {
    return counts[sigma.c2i((unsigned char)bwt[k])] + getNumberOfOcc(bwt[k], k);
}

length_t FMIndex::findSA(length_t index) const {
    // if index modulo the sparsefactor is zero then this row is in the sparse
    // suffix array
    if ((index & (sparseFactorSA - 1)) == 0) {
        return sa[index >> logSparseFactorSA];
    }

    // else iterate over LF mappings untill an index is found that is in the
    // sparse suffix array.
    length_t j = 0;
    while ((index & (sparseFactorSA - 1)) != 0) {
        j++;
        index = findLF(index);
    }

    // return the entry of the new found index plus the number of iterations
    // modulo the lenght of the bwt
    length_t ret = sa[index >> logSparseFactorSA] + j;
    if (ret >= bwt.size()) {
        return ret - bwt.size();
    }
    return ret;
}

vector<length_t> FMIndex::exactMatches(const string& s) {
    // find the range in the suffix array that matches the string
    Range range = matchString(s);

    // declare the return vector
    vector<length_t> positions;
    positions.reserve(range.width());

    // fill in the vector with all values in this range in the suffix array
    for (length_t i = range.begin; i < range.end; i++) {
        positions.push_back(findSA(i));
    }

    // sort the vector and return
    sort(positions.begin(), positions.end());
    return positions;
}

Range FMIndex::matchString(const string& s) {
    // start at the end
    auto it = s.crbegin();

    // find the range for this intitial character in the BWT string
    length_t positionInAlphabet = sigma.c2i((unsigned char)*it);

    length_t start = counts[positionInAlphabet];
    length_t end;
    if (positionInAlphabet != sigma.size() - 1) {
        end = counts[positionInAlphabet + 1];
    } else {
        end = bwt.size();
    }
    nodeCounter++;

    // iterate starting from the second character over the string
    for (++it; it != s.crend(); it++) {
        // find number of occurences of this char before and after and so the
        // new range is found
        positionInAlphabet = sigma.c2i((unsigned char)*it);
        length_t startOfChar = counts[positionInAlphabet];
        start = getNumberOfOcc(*it, start) + startOfChar;
        end = getNumberOfOcc(*it, end) + startOfChar;
        nodeCounter++;
        if (start == end) {
            // no matches found
            return Range();
        }
    }

    // return this range as a pair
    return Range(start, end);
}

// ----------------------------------------------------------------------------
// ROUTINES FOR APPROXIMATE PATTERN MATCHING
// ----------------------------------------------------------------------------

vector<AppMatch> FMIndex::approxMatches(const string& pattern, int maxED) {
    resetCounters();
    vector<AppMatchSA> occurences;

    // Intialize a band-diagonal matrix. The length of the diagonal should
    // be the (P.size + maxED + 1) because (P.size + maxED) is the largest
    // possible sequence we will align. The +1 is for the extra 0th row/col.
    EditMatrix matrix(pattern.size() + maxED + 1, maxED, 0, false, false);

    // set the intial range this is the entire bwt
    Range initialRange = Range(0, bwt.size());

    recApproxMatchesNaive(initialRange, pattern, matrix, maxED, occurences, 0);

    // in the occurence vector are ranges in the suffix array that have their
    // correspondent edit distance these indexes correspond to values that are
    // the positions in the original string
    return mapOccurencesInSAToOccurencesInText(occurences, maxED);
}

vector<pair<Range, char>> FMIndex::getCharExtensions(const Range& range) {
    vector<pair<Range, char>> nextChars;
    nextChars.reserve(range.width());

    // iterate over the entire alphabet
    for (const auto& c : sigma()) {

        length_t occBef = getNumberOfOcc(c, range.begin);
        length_t occAft = getNumberOfOcc(c, range.end);

        // check if this character occurs in the specified range
        if (occAft > occBef) {
            length_t s = occBef + counts[sigma.c2i((unsigned char)c)];
            length_t e = occAft + counts[sigma.c2i((unsigned char)c)];
            Range newRange = Range(s, e);

            // push this range and character for the next iteration
            nextChars.push_back(make_pair(newRange, c));
        }
    }

    return nextChars;
}

void FMIndex::recApproxMatchesNaive(const Range& range, const string& P,
                                    EditMatrix& M, const int& maxED,
                                    vector<AppMatchSA>& occ, length_t depth) {

    // store all character extensions ((characters before this character)
    vector<pair<Range, char>> nextChar = getCharExtensions(range);
    nodeCounter += nextChar.size();

    // report if fully alligned and also check the children
    // check the last entry in the matrix at collumn corresponding to first
    // character of P
    EditDistance ED =
        (depth + maxED >= P.size()) ? M(depth, P.size()) : maxED + 1;

    if ((int)ED.getED() <= maxED) {
        // push if this is better than th allowed edit distance
        occ.push_back(makeAppMatchSA(range, ED.getED(), depth));
    }

    // checking the children so going one row down
    depth++;

    // Matrix range we're filling in: M[i, startj:endj]
    // index = vertical (index i) ; pattern P = horizontal (index j)
    length_t firstCol = max<int>(1, depth - maxED);
    length_t lastCol = min<int>(P.size(), depth + maxED);

    for (length_t ci = 0; ci < nextChar.size(); ci++) {
        char c = nextChar[ci].second;

        // compute the next edit distance matrix column
        unsigned int minChdED = maxED + 1;
        for (length_t j = firstCol; j <= lastCol; j++) {
            // search is in reverse order so the at operator needs to be
            // inversed
            M.updateMatrix(P[P.size() - j] != c, depth, j);
            matrixElementCounter++;
            minChdED = min(minChdED, M(depth, j).getED());
        }

        // edit distance threshold exceeded: backtracking
        if ((int)minChdED > maxED) {
            continue;
        }
        // continue the search on this child
        recApproxMatchesNaive(nextChar[ci].first, P, M, maxED, occ, depth);
    }
}

// ----------------------------------------------------------------------------
// POST-PROCESSING ROUTINES FOR APPROXIMAtE PATTERN MATCHING
// ----------------------------------------------------------------------------

vector<AppMatch>
FMIndex::mapOccurencesInSAToOccurencesInText(const vector<AppMatchSA>& occ,
                                             const int& maxED) {
    vector<AppMatch> occurencesInText;
    occurencesInText.reserve(bwt.size());
    // map the startposition to the best match (lowest edit distance)
    map<length_t, AppMatch> posToBestMatch;

    if (occ.size() == 0) {
        return {};
    }
    if (occ.size() == 1) {
        // all occ are distinct
        positionsInPostProcessingCounter = occ[0].match.range.width();
        auto m = convertToMatchesInText(occ[0]);
        sort(m.begin(), m.end());
        return m;
    }

    // more than 1 reported occurence, could be redundant
    for (const auto& it : occ) {

        Range range = it.match.range;
        positionsInPostProcessingCounter += range.width();

        auto matchesInTextToCheck = convertToMatchesInText(it);
        occurencesInText.insert(occurencesInText.end(),
                                matchesInTextToCheck.begin(),
                                matchesInTextToCheck.end());
    }

    sort(occurencesInText.begin(), occurencesInText.end());

    std::vector<AppMatch> nonRedundantOcc;
    nonRedundantOcc.reserve(occurencesInText.size());

    length_t maxDiff = 2 * maxED;
    length_t prevBegin = numeric_limits<length_t>::max();
    int prevED = maxED + 1;
    length_t prevDepth = numeric_limits<length_t>::max();

    for (const auto& o : occurencesInText) {
        auto diff = abs_diff<length_t>(o.range.begin, prevBegin);
        if (diff == 0) {
            continue;
        }
        if (diff <= maxDiff) {
            // check if this later occurrence is better than the previous one
            if (o.editDist > prevED) {
                continue;
            }
            if (o.editDist == prevED && o.range.width() >= prevDepth) {
                continue;
            }

            // prev was worse so pop_back
            nonRedundantOcc.pop_back();
        }

        nonRedundantOcc.push_back(o);
        prevBegin = o.range.begin;
        prevED = o.editDist;
        prevDepth = o.range.width();
    }

    return nonRedundantOcc;
}

AppMatch FMIndex::verifyInText(length_t s, length_t e, int remED,
                               const Substring& pat) const {

    string text = getText();
    // the part of the text that needs to be checked
    Substring textAt = Substring(text, s, min(e, s + remED + pat.size()));

    // need to allign startOfText to string
    // make matrix
    length_t m = textAt.size();
    length_t n = pat.size();
    EditMatrix matrix(n + remED + 1, remED, 0, s != 0, true);

    length_t currentBestRow = 0;
    length_t currentBestED = remED + 1;
    // fill in the matrix
    for (length_t i = 1; i <= m; i++) {
        int bestEDRow = remED + 1; // collumns of the matrix for this row
        length_t firstCol = max<int>(1, i - remED);
        length_t lastCol = min<int>(n, i + remED);

        for (length_t j = firstCol; j <= lastCol; j++) {
            matrix.updateMatrix(textAt[i - 1] != pat[j - 1], i, j);
            bestEDRow = min<length_t>(bestEDRow, matrix(i, j).getED());
        }

        if (bestEDRow > remED) {
            // no allignment possible -> return empty match
            AppMatch m;
            m.range = Range();
            m.editDist = remED + 1;
            return m;
        }

        if (lastCol == n && matrix(i, lastCol).getED() < currentBestED) {
            currentBestED = matrix(i, lastCol).getED();
            currentBestRow = i;
        }
    }
    // make a match in the text
    AppMatch match;
    match.range = Range(s, s + currentBestRow);
    match.editDist = currentBestED;
    return match;
}

vector<AppMatch> FMIndex::convertToMatchesInText(const AppMatchSA& saMatch) {

    vector<AppMatch> textMatches;
    textMatches.reserve(saMatch.match.range.width());

    for (length_t i = saMatch.match.range.begin; i < saMatch.match.range.end;
         i++) {
        // find the startPosition in the text by looking at the SA
        length_t startPos = findSA(i);

        length_t endPos = startPos + saMatch.depth;

        AppMatch newMatch;

        // does not go over sentinel, create simple match
        newMatch.range = Range(startPos, endPos);
        newMatch.editDist = saMatch.match.editDist;

        textMatches.push_back(newMatch);
    }
    return textMatches;
}
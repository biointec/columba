/******************************************************************************
 *  Columba: Approximate Pattern Matching using Search Schemes                *
 *  Copyright (C) 2020-2021 - Luca Renders <luca.renders@ugent.be> and        *
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
#include "fmindex.h"

using namespace std;

std::ostream& operator<<(std::ostream& os, const Range& r) {
    os << "[" << r.begin << ", " << r.end << ")";
    return os;
}

// ============================================================================
// CLASS BIFMOCC
// ============================================================================

std::ostream& operator<<(std::ostream& o, const FMOcc& m) {
    return o << "SARange: " << m.getRanges().getRangeSA()
             << "\tEdit distance: " << m.getDistance()
             << "\tdepth: " << m.getDepth();
}

std::ostream& operator<<(std::ostream& os, const Search& obj) {
    os << "{";
    int numParts = obj.getNumParts();
    for (int i = 0; i < numParts; i++) {
        os << obj.getPart(i) << ((i == numParts - 1) ? "}" : ",");
    }
    os << " {";
    for (int i = 0; i < numParts; i++) {
        os << obj.getLowerBound(i) << ((i == numParts - 1) ? "}" : ",");
    }
    os << " {";
    for (int i = 0; i < numParts; i++) {
        os << obj.getUpperBound(i) << ((i == numParts - 1) ? "}" : ",");
    }
    return os;
}

// ============================================================================
// CLASS FMIndex
// ============================================================================
thread_local std::vector<std::vector<FMPosExt>> FMIndex::stacks;
thread_local length_t FMIndex::nodeCounter;
thread_local length_t FMIndex::matrixElementCounter;
thread_local length_t FMIndex::positionsInPostProcessingCounter;
thread_local Direction FMIndex::dir = BACKWARD;
thread_local ExtraCharPtr FMIndex::extraChar;
// ----------------------------------------------------------------------------
// ROUTINES FOR ACCESSING DATA STRUCTURE
// ----------------------------------------------------------------------------

length_t FMIndex::findLF(length_t k, bool reversed) const {
    if (!reversed) {

        return counts[sigma.c2i((unsigned char)bwt[k])] +
               getNumberOfOcc(sigma.c2i((unsigned char)bwt[k]), k);
    }
    length_t posInAlphabet = sigma.c2i((unsigned char)bwt[k]);

    return counts[posInAlphabet] + getNumberOfOccRev(posInAlphabet, k);
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
        index = findLF(index, false);
    }

    // return the entry of the new found index plus the number of iterations
    // modulo the lenght of the bwt
    length_t ret = sa[index >> logSparseFactorSA] + j;
    if (ret >= bwt.size()) {
        return ret - bwt.size();
    }
    return ret;
}

// ----------------------------------------------------------------------------
// ROUTINES FOR INITIALIZATION
// ----------------------------------------------------------------------------

void FMIndex::fromFiles(const string& baseFile, bool verbose) {
    if (verbose) {
        cout << "Reading in files with baseFile " << baseFile << endl;
        // read the text
        cout << "Reading " << baseFile << ".txt" << endl;
    }
    string text;
    if (!readText(baseFile + ".txt", text)) {
        throw runtime_error("Problem reading: " + baseFile + ".txt");
    }

    textLength =
        (text[text.size() - 1] == '\n') ? text.size() - 1 : text.size();
    if (verbose) {
        cout << "Done reading text (size: " << text.size() << ")" << endl;

        // read the counts table

        cout << "Reading " << baseFile << ".cct"
             << "\n";
    }

    // read the counts table
    vector<length_t> charCounts(256, 0);
    if (!readArray(baseFile + ".cct", charCounts)) {
        throw runtime_error("Cannot open file: " + baseFile + ".cct");
    }

    length_t cumCount = 0; // cumulative character counts
    for (size_t i = 0; i < charCounts.size(); i++) {
        if (charCounts[i] == 0)
            continue;
        counts.push_back(cumCount);
        cumCount += charCounts[i];
    }
    sigma = Alphabet<ALPHABET>(charCounts);

    // read the SA
    if (verbose) {
        cout << "Reading " << baseFile << ".sa." << sparseFactorSA << endl;
    }
    string saFilename = baseFile + ".sa." + to_string(sparseFactorSA);

    if (!readArray(saFilename, sa)) {
        throw runtime_error("Cannot open file: " + saFilename);
    }
    if (verbose) {
        cout << "Done reading suffix array (sparseness factor: "
             << sparseFactorSA << ")" << endl;

        // read the BWT
        cout << "Reading " << baseFile << ".bwt" << endl;
    }
    if (!readText(baseFile + ".bwt", bwt)) {
        throw runtime_error("Cannot open file: " + baseFile + ".bwt");
    }
    if (verbose) {
        cout << "Done reading BWT (size: " << bwt.size() << ")" << endl;

        // read the baseFile occurrence table
        cout << "Reading " << baseFile << ".brt" << endl;
    }

    if (!fwdRepr.read(baseFile + ".brt"))
        throw runtime_error("Cannot open file: " + baseFile + ".brt");
    if (verbose) {
        cout << "Done reading baseFile occurrence table" << endl;
    }

    // read the reverse baseFile occurrence table
    if (!revRepr.read(baseFile + ".rev.brt"))
        throw std::runtime_error("Cannot open file: " + baseFile + ".rev.brt");
    if (verbose) {
        cout << "Done reading reverse baseFile occurrence table" << std::endl;
    }
}

void FMIndex::populateTable(bool verbose) {
    if (verbose) {
        cout << "Populating FM-range table with " << wordSize << "-mers...";
    }
    cout.flush();

    Kmer::setWordSize(wordSize);

    table.resize(1 << (2 * wordSize)); // 2 << wordSize is 4^wordSize
    setDirection(FORWARD);

    string word;
    vector<FMPosExt> stack;
    extendFMPos(getCompleteRange(), stack);
    while (!stack.empty()) {
        auto curr = stack.back();
        stack.pop_back();

        word.resize(curr.getRow());
        word[curr.getRow() - 1] = curr.getCharacter();

        if ((length_t)curr.getRow() == wordSize) { // max depth reached
            Kmer k(word);
            table.insert(make_pair(k, curr.getRanges()));

        } else // add extra characters
            extendFMPos(curr.getRanges(), stack, curr.getRow());
    }
    if (verbose) {
        cout << "done." << endl;
    }
}

// ----------------------------------------------------------------------------
// ROUTINES FOR EXACT PATTERN MATCHING
// ----------------------------------------------------------------------------
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
        start = getNumberOfOcc(positionInAlphabet, start) + startOfChar;
        end = getNumberOfOcc(positionInAlphabet, end) + startOfChar;
        nodeCounter++;
        if (start == end) {
            // no matches found
            return Range();
        }
    }

    // return this range as a pair
    return Range(start, end);
}
vector<length_t> FMIndex::exactMatches(const string& s) {
    // find the range in the suffix array that matches the string
    Range range = matchString(s);

    // declare the return vector
    vector<length_t> positions;
    positions.reserve(range.width());

    // fill in the vector with all values in this range in the suffix array
    for (length_t i = range.getBegin(); i < range.getEnd(); i++) {
        positions.push_back(findSA(i));
    }

    // sort the vector and return
    sort(positions.begin(), positions.end());
    return positions;
}

SARangePair FMIndex::matchStringBidirectionally(const Substring& pattern,
                                                SARangePair rangesOfPrev) {

    for (length_t i = 0; i < pattern.size(); i++) {

        char c = pattern[i];
        if (!addChar(c, rangesOfPrev)) {
            // rangesOfPrev was made empty
            break;
        }
    }

    return rangesOfPrev;
}
bool FMIndex::addChar(const char& c, SARangePair& startRange) {

    int posInAlphabet = sigma.c2i((unsigned char)c);
    if (posInAlphabet > -1) {

        if ((this->*extraChar)(posInAlphabet, startRange, startRange)) {
            // each character that we look at is a new node that is visited
            nodeCounter++;
            return true;
        }
    }
    // the range is now empty

    return false;
}

// ----------------------------------------------------------------------------
// ROUTINES FOR APPROXIMATE PATTERN MATCHING
// ----------------------------------------------------------------------------

std::vector<TextOccurrence>
FMIndex::approxMatchesNaive(const std::string& pattern, length_t maxED) {

    resetCounters();
    vector<FMOcc> occurences;

    BandMatrix matrix(pattern.size() + maxED + 1, maxED);

    setDirection(FORWARD);

    std::vector<FMPosExt> stack;
    stack.reserve((pattern.size() + maxED + 1) * (sigma.size() - 1));

    extendFMPos(getCompleteRange(), stack, 0);

    while (!stack.empty()) {
        const FMPosExt currentNode = stack.back();
        stack.pop_back();
        int row = currentNode.getDepth();
        length_t first = matrix.getFirstColumn(row);
        length_t last = matrix.getLastColumn(row);
        length_t minimalED = maxED + 1;
        for (length_t j = first; j <= last; j++) {
            matrix.updateMatrix(currentNode.getCharacter() != pattern[j - 1],
                                row, j);
            minimalED = min(minimalED, matrix(row, j));
        }

        if (minimalED > maxED) {
            // backtrack
            continue;
        }

        if (last == (length_t)matrix.getLastColumn()) {
            // full pattern was matched
            if (matrix(row, last) <= maxED) {
                occurences.emplace_back(currentNode, matrix(row, last));
            }
        }

        extendFMPos(currentNode, stack);
    }

    return mapOccurencesInSAToOccurencesInText(occurences, maxED);
}

bool FMIndex::findRangesWithExtraCharBackward(
    length_t positionInAlphabet, const SARangePair& rangesOfP,
    SARangePair& rangesOfChild) const {

    // first make the trivial range by  searching cP using B
    Range trivialRange = rangesOfP.getRangeSA();

    // find the new range by using the LF property
    length_t occBefore =
        getNumberOfOcc(positionInAlphabet, trivialRange.getBegin());
    length_t occAfter =
        getNumberOfOcc(positionInAlphabet, trivialRange.getEnd());

    length_t startInAlphabet = counts[positionInAlphabet];
    Range range1 =
        Range(occBefore + startInAlphabet, occAfter + startInAlphabet);

    // then make the less trivial range by counting the sizes of the ranges
    // of (dP) using B

    // first get the start of the range we are looking for of the parent
    length_t s = rangesOfP.getRangeSARev().getBegin();

    // get the start of the child within this range
    // find the number of occurences of chars smaller than c in the parent
    // range

    length_t x = getNumberOfCumOcc(positionInAlphabet, trivialRange.getEnd()) -
                 getNumberOfCumOcc(positionInAlphabet, trivialRange.getBegin());

    // make the new range with width equal to that of the trivial range
    Range range2 = Range(s + x, s + x + range1.width());

    rangesOfChild = SARangePair(range1, range2);
    return !rangesOfChild.empty();
}

bool FMIndex::findRangesWithExtraCharForward(length_t positionInAlphabet,
                                             const SARangePair& rangesOfP,
                                             SARangePair& childRanges) const {

    // first make the trivial range by searching (Pc)' using B'  if
    // searching forward we need to use B' so we need the reverse range
    Range rangeForTrivial = rangesOfP.getRangeSARev();

    // find the new range by using the LF property
    length_t occBefore =
        getNumberOfOccRev(positionInAlphabet, rangeForTrivial.getBegin());
    length_t occAfter =
        getNumberOfOccRev(positionInAlphabet, rangeForTrivial.getEnd());

    length_t startInAlphabet = counts[positionInAlphabet];
    Range range1 =
        Range(occBefore + startInAlphabet, occAfter + startInAlphabet);

    // then make the less trivial range by counting the size of the range of
    // (Pd)' using B' (forward)

    // first get the start of the range we are looking for of the parent
    length_t s = rangesOfP.getRangeSA().getBegin();

    // get the start of the child within this range
    // find the number of occurences of chars smaller than c in the parent
    // range
    Range prevRange = rangesOfP.getRangeSARev();

    length_t x = getNumberOfCumOccRev(positionInAlphabet, prevRange.getEnd()) -
                 getNumberOfCumOccRev(positionInAlphabet, prevRange.getBegin());

    // make the new range
    Range range2 = Range(s + x, s + x + range1.width());

    childRanges = SARangePair(range2, range1);
    return !childRanges.empty();
}

void FMIndex::recApproxMatchEditNaive(const Search& s, const FMOcc& startMatch,
                                      vector<FMOcc>& occ,
                                      const vector<Substring>& parts,
                                      const int& idx) {
    const Substring& p = parts[s.getPart(idx)];           // this part
    const length_t& maxED = s.getUpperBound(idx);         // maxED for this part
    const length_t& minED = s.getLowerBound(idx);         // minED for this part
    const length_t& W = maxED - startMatch.getDistance(); // Width of matrix
    const length_t& pSize = p.size();

    BandMatrix matrix = BandMatrix(pSize, W, startMatch.getDistance());

    if (matrix.getLastColumn(0) == (int)pSize && matrix(0, pSize) <= maxED) {
        // an occurrence found by gapping entire part
        if (s.isEnd(idx)) {
            occ.emplace_back(startMatch.getRanges(), matrix(0, pSize),
                             startMatch.getDepth());
        } else {
            // go to the next index
            recApproxMatchEditNaive(s,
                                    FMOcc(startMatch.getRanges(),
                                          matrix(0, pSize),
                                          startMatch.getDepth()),
                                    occ, parts, idx + 1);
            // set direction correct again
            setDirection(dir);
        }
    }

    auto& stack = stacks[idx];                  // stack for this partition
    const Direction& dir = s.getDirection(idx); // direction
    setDirection(dir);

    extendFMPos(startMatch.getRanges(), stack);

    while (!stack.empty()) {
        const auto currentNode = stack.back();
        stack.pop_back();

        const int& row = currentNode.getRow();
        const length_t firstCol = matrix.getFirstColumn(row);
        const length_t lastCol = matrix.getLastColumn(row);

        length_t minimalEDOfRow = matrix(row, firstCol - 1);
        for (length_t col = firstCol; col <= lastCol; col++) {
            length_t filledIn = matrix.updateMatrix(
                currentNode.getCharacter() != p[col - 1], row, col);
            matrixElementCounter++;
            minimalEDOfRow = std::min(minimalEDOfRow, filledIn);
        }

        if (minimalEDOfRow > maxED) {
            // backtracking
            continue;
        }

        if (lastCol == pSize && matrix(row, pSize) <= maxED &&
            matrix(row, pSize) >= minED) {
            if (s.isEnd(idx)) {
                occ.emplace_back(currentNode.getRanges(), matrix(row, pSize),
                                 startMatch.getDepth() +
                                     currentNode.getDepth());
            } else {
                // go deeper in search
                recApproxMatchEditNaive(
                    s,
                    FMOcc(currentNode.getRanges(), matrix(row, pSize),
                          startMatch.getDepth() + currentNode.getDepth()),
                    occ, parts, idx + 1);
                // set direction correct again
                setDirection(dir);
            }
        }
        if (firstCol == pSize) {
            // final cell of matrix filled in => backtracking
            continue;
        }

        // extend to the children
        extendFMPos(currentNode, stack);
    }
}

void FMIndex::recApproxMatchEditOptimized(
    const Search& s, const FMOcc& startMatch, vector<FMOcc>& occ,
    const vector<Substring>& parts, const int& idx,
    const vector<FMPosExt>& descPrevDir, const vector<int>& initPrevDir,
    const vector<FMPosExt>& descNotPrevDir, const vector<int>& initNotPrevDir) {

    // shortcut Variables
    const Substring& p = parts[s.getPart(idx)];           // this part
    const length_t& maxED = s.getUpperBound(idx);         // maxED for this part
    const length_t& W = maxED - startMatch.getDistance(); // Width of matrix
    const Direction& dir = s.getDirection(idx);           // direction
    const bool& dSwitch = s.getDirectionSwitch(idx); // has direction switched?
    auto& stack = stacks[idx];                       // stack for this partition

    // get the correct initED and descendants based on switch
    const vector<int>& initEds = dSwitch ? initNotPrevDir : initPrevDir;
    const vector<FMPosExt>& descendants =
        dSwitch ? descNotPrevDir : descPrevDir;
    const vector<int>& initOther = dSwitch ? initPrevDir : initNotPrevDir;
    const vector<FMPosExt>& descOther = dSwitch ? descPrevDir : descNotPrevDir;

    setDirection(dir);

    // calculate nescessary increase for first column of bandmatrix
    int prevED = 0;
    if (!initEds.empty()) {
        prevED = (dSwitch ? *min_element(initEds.begin(), initEds.end())
                          : initEds[0]);
    }
    int increase = startMatch.getDistance() - prevED;

    // initialize matrix and cluster
    BandMatrix matrix(p.size(), W, increase, initEds);
    Cluster clus(matrix.getSizeOfFinalColumn(), maxED, startMatch.getDepth(),
                 startMatch.getShift());

    if (matrix.getLastColumn(0) == (int)p.size()) {
        // the first row is part of the final column of the banded matrix
        // Update the first cell of the cluster with the startmatch and the
        // value found at the psizeth column of the initialization row the
        // character does not matter as the first cell of a cluster is never a
        // descendant of any of the other cells in the cluster, so this cell
        // will not be reused for the next part of the pattern
        clus.setValue(0, FMPosExt((char)0, startMatch.getRanges(), 0),
                      matrix(0, p.size()));
    }

    if (!descendants.empty()) {
        // fill in the matrix for the descendants

        int maxRow = matrix.getNumberOfRows();
        int b = 0;
        for (length_t i = 0;
             i < descendants.size() && descendants[i].getDepth() <= maxRow;
             i++) {
            // reset the depth of this node

            b = branchAndBound(
                matrix, clus, descendants[i], s, idx, parts, occ, initOther,
                descOther, {descendants.begin() + i + 1, descendants.end()});
        }
        if (b || descendants.back().getDepth() == maxRow) {
            // minimal ed exceeded no more options to get lower ed from
            // first column, or no more rows to possibly check
            return;
        }

        SARangePair pair = descendants.back().getRanges();
        if (dSwitch) {
            // after a switch the ranges of final descendant should be
            // updated
            pair = startMatch.getRanges();
        }

        // push children of final descendant
        extendFMPos(pair, stack, descendants.back().getDepth());
    } else { // get the initial nodes to check
        extendFMPos(startMatch.getRanges(), stack);
    }

    while (!stack.empty()) {

        const FMPosExt currentNode = stack.back();
        stack.pop_back();

        if (branchAndBound(matrix, clus, currentNode, s, idx, parts, occ,
                           initOther, descOther)) {

            continue;
        }

        // continue the search for children of this node
        extendFMPos(currentNode, stack);
    }
}

bool FMIndex::branchAndBound(BandMatrix& matrix, Cluster& clus,
                             const FMPosExt& currentNode, const Search& s,
                             const int& idx, const vector<Substring>& parts,
                             vector<FMOcc>& occ, const vector<int>& initOther,
                             const vector<FMPosExt>& descOther,
                             const vector<FMPosExt> remainingDesc) {

    const Substring& p = parts[s.getPart(idx)]; // The current partition
    const length_t& pSize = p.size(); // The size of the current partition

    const length_t& lowerBound =
        s.getLowerBound(idx); // lowerbound for this part
    const length_t& maxED =
        s.getUpperBound(idx); // the upperbound for this part

    const int row = currentNode.getDepth();

    // get first and last column of matrix for this row
    const length_t firstCol = matrix.getFirstColumn(row);
    const length_t lastCol = matrix.getLastColumn(row);

    // only necessary in lastCol == pSize
    bool afterwardsOnlyVerticalGap = true;

    unsigned int minimalEDOfRow = matrix(row, firstCol - 1);
    length_t lastUpdatedMatrixValue = 0;
    for (length_t j = firstCol; j <= lastCol; j++) {
        lastUpdatedMatrixValue =
            matrix.updateMatrix(p[j - 1] != currentNode.getCharacter(), row, j);
        matrixElementCounter++;

        if (minimalEDOfRow <= lastUpdatedMatrixValue) {
            afterwardsOnlyVerticalGap = false;
        }

        if (lastCol == j && minimalEDOfRow > maxED) {
            // this branch will always be executed exactly once if this is the
            // last row
            afterwardsOnlyVerticalGap = true;
        }
        minimalEDOfRow = min(minimalEDOfRow, lastUpdatedMatrixValue);
    }

    if (lastCol == pSize) {
        // update the cluster
        int clusIdx = row - matrix.getNumberOfRows() + clus.size() - 1;
        clus.setValue(clusIdx, currentNode, lastUpdatedMatrixValue);

        if (minimalEDOfRow > maxED || afterwardsOnlyVerticalGap) {
            // no need to further explore this branch for this part -> go to
            // next part
            goDeeper(clus, idx + 1, s, parts, occ, lowerBound, descOther,
                     initOther, remainingDesc);
            return true;
        }
    }

    // edit distance threshold exceeded: backtracking
    if (minimalEDOfRow > maxED) {
        return true;
    }

    return false;
}

void FMIndex::goDeeper(Cluster& cluster, const length_t& nextIdx,
                       const Search& s, const vector<Substring>& parts,
                       vector<FMOcc>& occ, const length_t& lowerBound,
                       const vector<FMPosExt>& descOtherD,
                       const vector<int>& initOtherD,
                       const vector<FMPosExt>& remainingDesc) {

    bool isEdge = s.isEdge(nextIdx - 1);

    if (isEdge) {
        // if this is final piece report highest minimum (to get shortest
        // match)
        if (nextIdx == parts.size()) {
            auto matches = cluster.reportCentersAtEnd();
            for (const auto& match : matches) {
                if (match.isValid() && match.getDistance() >= (int)lowerBound) {
                    occ.emplace_back(match);
                }
            }
        } else {
            FMOcc match = cluster.reportDeepestMinimum(this->dir);
            if (match.isValid() && match.getDistance() >= (int)lowerBound) {
                // go deeper in search
                Direction originalDir = this->dir;
                recApproxMatchEditOptimized(s, match, occ, parts, nextIdx, {},
                                            {}, descOtherD, initOtherD);
                // set direction back again
                setDirection(originalDir);
            }
        }

        return;
    }

    // one of the later stages will return to this point, so keep track of
    // the descendants and eds at this branch
    vector<FMPosExt> descendants;
    vector<int> initEds;

    FMOcc newMatch = cluster.getClusterCentra(lowerBound, descendants, initEds);
    if (!newMatch.isValid()) {
        // no centre above lowerbound found
        return;
    }

    // add the remaining descendants
    descendants.insert(descendants.end(), remainingDesc.begin(),
                       remainingDesc.end());

    // reset the depth of all the descendants
    for (int i = 0; i < (int)descendants.size(); i++) {
        descendants[i].setDepth(i + 1);
    }
    int maxEDNext = s.getUpperBound(nextIdx);

    // remove trailing initEds that are higher than maxEDNext
    while (initEds.back() > maxEDNext) {
        initEds.pop_back();
    }

    // is the next direction equal to this direction?
    bool switchAfter = s.getDirectionSwitch(nextIdx);

    if (switchAfter) {
        // switching direction as this is not the end of a search direction,
        // this means we'll get back here, thus range of newmatch should be
        // deepest point in branch
        if (!descendants.empty()) {
            newMatch.setRanges(descendants.back().getRanges());

            //  edit distance for search in other direction should be lowest
            //  value possible
            newMatch.setDistance(*min_element(initEds.begin(), initEds.end()));
        }

        Direction originalDir = this->dir;

        recApproxMatchEditOptimized(s, newMatch, occ, parts, nextIdx,
                                    descendants, initEds, descOtherD,
                                    initOtherD);

        // set direction back again
        setDirection(originalDir);
    } else {
        // go deeper on next piece
        recApproxMatchEditOptimized(s, newMatch, occ, parts, nextIdx,
                                    descendants, initEds, descOtherD,
                                    initOtherD);
    }
}

void FMIndex::recApproxMatchHamming(const Search& s, const FMOcc& startMatch,
                                    std::vector<FMOcc>& occ,
                                    const std::vector<Substring>& parts,
                                    const int& idx) {

    // shortcut variabeles
    const Substring& p = parts[s.getPart(idx)]; // the current part
    const int& pSize = p.size();                // the size of the current part
    const Direction& d = s.getDirection(idx);   // direction of current part
    const int& maxED = s.getUpperBound(idx);    // upperbound of current part
    const int& minED = s.getLowerBound(idx);    // lowerbound of currentpart
    setDirection(d);

    // create vector
    std::vector<int> vec(p.size() + 1, 0);
    // set root element of vector to the distance of the startmatch
    vec[0] = startMatch.getDistance();
    // get stack for current part
    auto& stack = stacks[idx];

    extendFMPos(startMatch.getRanges(), stack);

    while (!stack.empty()) {
        const FMPosExt node = stack.back();
        stack.pop_back();

        // update the vector
        int row = node.getRow();
        vec[row] = vec[row - 1] + (node.getCharacter() != p[row - 1]);
        if (vec[row] > maxED) {
            // backtrack
            continue;
        }

        if (row == pSize) {
            // end of part
            if (vec[row] >= minED) {
                // valid occurrence
                FMOcc match = FMOcc(node.getRanges(), vec[row],
                                    startMatch.getDepth() + pSize);
                if (s.isEnd(idx)) {
                    // end of search
                    occ.push_back(match);
                } else {
                    // contine search
                    recApproxMatchHamming(s, match, occ, parts, idx + 1);
                }
            }
            continue;
        }
        extendFMPos(node, stack);
    }
}

void FMIndex::extendFMPos(const SARangePair& parentRanges,
                          vector<FMPosExt>& stack, int row) {

    // iterate over the entire alphabet
    for (length_t i = 1; i < sigma.size(); i++) {

        SARangePair pairForNewChar;

        // check if this character occurs in the specified range
        if ((this->*extraChar)(i, parentRanges, pairForNewChar)) {
            // push this range and character for the next iteration
            stack.emplace_back(sigma.i2c(i), pairForNewChar, row + 1);

            nodeCounter++;
        }
    }
}

void FMIndex::extendFMPos(const FMPosExt& pos, vector<FMPosExt>& stack) {
    extendFMPos(pos.getRanges(), stack, pos.getDepth());
}

// ----------------------------------------------------------------------------
// POST-PROCESSING ROUTINES FOR APPROXIMAtE PATTERN MATCHING
// ----------------------------------------------------------------------------

vector<TextOccurrence> FMIndex::convertToMatchesInText(const FMOcc& saMatch) {

    vector<TextOccurrence> textMatches;
    textMatches.reserve(saMatch.getRanges().width());

    for (length_t i = saMatch.getRanges().getRangeSA().getBegin();
         i < saMatch.getRanges().getRangeSA().getEnd(); i++) {
        // find the startPosition in the text by looking at the SA
        length_t startPos = findSA(i) + saMatch.getShift();

        length_t endPos = startPos + saMatch.getDepth();

        textMatches.emplace_back(Range(startPos, endPos),
                                 saMatch.getDistance());
    }
    return textMatches;
}

vector<TextOccurrence>
FMIndex::mapOccurencesInSAToOccurencesInText(vector<FMOcc>& occ,
                                             const int& maxED) {

    sort(occ.begin(), occ.end());
    occ.erase(unique(occ.begin(), occ.end()), occ.end());
    vector<TextOccurrence> occurencesInText;
    occurencesInText.reserve(1000 * maxED);
    // map the startposition to the best match (lowest edit distance)
    map<length_t, TextOccurrence> posToBestMatch;

    if (occ.size() == 0) {
        return {};
    }
    if (occ.size() == 1) {
        // all occ are distinct
        positionsInPostProcessingCounter = occ[0].getRanges().width();
        auto m = convertToMatchesInText(occ[0]);
        sort(m.begin(), m.end());
        for (auto& occ : m) {
            occ.generateOutput();
        }
        return m;
    }

    // more than 1 reported occurence, could be redundant
    for (const auto& it : occ) {

        const Range& range = it.getRanges().getRangeSA();
        positionsInPostProcessingCounter += range.width();

        auto matchesInTextToCheck = convertToMatchesInText(it);

        occurencesInText.insert(occurencesInText.end(),
                                matchesInTextToCheck.begin(),
                                matchesInTextToCheck.end());
    }

    sort(occurencesInText.begin(), occurencesInText.end());

    std::vector<TextOccurrence> nonRedundantOcc;
    nonRedundantOcc.reserve(occurencesInText.size());

    length_t maxDiff = 2 * maxED;
    length_t prevBegin = numeric_limits<length_t>::max();
    int prevED = maxED + 1;
    length_t prevDepth = numeric_limits<length_t>::max();

    for (const auto& o : occurencesInText) {
        auto diff = abs_diff<length_t>(o.getRange().getBegin(), prevBegin);
        if (diff == 0) {
            continue;
        }
        if (diff <= maxDiff) {
            // check if this later occurrence is better than the previous
            // one
            if (o.getDistance() > prevED) {
                continue;
            }
            if (o.getDistance() == prevED &&
                o.getRange().width() >= prevDepth) {
                continue;
            }

            // prev was worse so pop_back
            nonRedundantOcc.pop_back();
        }

        prevBegin = o.getRange().getBegin();
        prevED = o.getDistance();
        prevDepth = o.getRange().width();

        nonRedundantOcc.emplace_back(o);
    }
    for (TextOccurrence& occ : nonRedundantOcc) {
        occ.generateOutput();
    }

    return nonRedundantOcc;
}

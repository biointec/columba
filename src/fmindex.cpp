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
#include "fmindex.h"

using namespace std;

ostream& operator<<(ostream& os, const Range& r) {
    os << "[" << r.begin << ", " << r.end << ")";
    return os;
}

// ============================================================================
// CLASS BIFMOCC
// ============================================================================

ostream& operator<<(ostream& o, const FMOcc& m) {
    return o << "SARange: " << m.getRanges().getRangeSA()
             << "\tEdit distance: " << m.getDistance()
             << "\tdepth: " << m.getDepth();
}

ostream& operator<<(ostream& os, const Search& obj) {
    os << "{";
    length_t numParts = obj.getNumParts();
    for (length_t i = 0; i < numParts; i++) {
        os << obj.getPart(i) << ((i == numParts - 1) ? "}" : ",");
    }
    os << " {";
    for (length_t i = 0; i < numParts; i++) {
        os << obj.getLowerBound(i) << ((i == numParts - 1) ? "}" : ",");
    }
    os << " {";
    for (length_t i = 0; i < numParts; i++) {
        os << obj.getUpperBound(i) << ((i == numParts - 1) ? "}" : ",");
    }
    return os;
}

// ============================================================================
// CLASS FMIndex
// ============================================================================
thread_local length_t Counters::nodeCounter;
thread_local length_t Counters::abortedInTextVerificationCounter;
thread_local length_t Counters::totalReportedPositions;
thread_local length_t Counters::cigarsInTextVerification;
thread_local length_t Counters::cigarsInIndex;
thread_local length_t Counters::inTextStarted;
thread_local length_t Counters::usefulCigarsInText;
thread_local length_t Counters::immediateSwitch;
thread_local length_t Counters::approximateSearchStarted;

thread_local Direction FMIndex::dir = BACKWARD;
thread_local ExtraCharPtr FMIndex::extraChar;
thread_local FindDiffPtr FMIndex::findDiff;

thread_local vector<vector<FMPosExt>> FMIndex::stacks;
thread_local vector<BitParallelED> FMIndex::matrices;

// ----------------------------------------------------------------------------
// ROUTINES FOR ACCESSING DATA STRUCTURE
// ----------------------------------------------------------------------------
length_t FMIndex::findLF(length_t k) const {

    const auto& pos = sigma.c2i((unsigned char)bwt[k]);
    return counts[pos] + getNumberOfOcc(pos, k);
}

length_t FMIndex::findSA(length_t index) const {
    length_t l = 0;
    while (!sparseSA[index]) {
        index = findLF(index);
        l++;
    }
    return sparseSA.get(index) + l;
}

// ----------------------------------------------------------------------------
// ROUTINES FOR INITIALIZATION
// ----------------------------------------------------------------------------

void FMIndex::fromFiles(const string& baseFile, bool verbose) {
    if (verbose) {

        // read the text
        cout << "Reading " << baseFile << ".txt"
             << "...";
        cout.flush();
    }

    if (!readText(baseFile + ".txt", text)) {
        throw runtime_error("Problem reading: " + baseFile + ".txt");
    }

    textLength =
        (text[text.size() - 1] == '\n') ? text.size() - 1 : text.size();
    if (verbose) {
        cout << "done (size: " << text.size() << ")" << endl;

        // read the counts table

        cout << "Reading " << baseFile << ".cct"
             << "...";
        cout.flush();
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

    if (verbose) {
        cout << "done" << endl;
        // read the BWT
        cout << "Reading " << baseFile << ".bwt"
             << "...";
        cout.flush();
    }

    if (!readText(baseFile + ".bwt", bwt)) {
        throw runtime_error("Cannot open file: " + baseFile + ".bwt");
    }
    if (verbose) {
        cout << "done (size: " << bwt.size() << ")" << endl;

        // read the baseFile occurrence table
        cout << "Reading " << baseFile << ".brt"
             << "...";
        cout.flush();
    }

    if (!fwdRepr.read(baseFile + ".brt"))
        throw runtime_error("Cannot open file: " + baseFile + ".brt");
    if (verbose) {
        cout << "done" << endl;
        cout << "Reading " << baseFile << ".rev.brt"
             << "...";
    }

    // read the reverse baseFile occurrence table
    if (!revRepr.read(baseFile + ".rev.brt"))
        throw runtime_error("Cannot open file: " + baseFile + ".rev.brt");
    if (verbose) {
        cout << "done" << endl;
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
    Counters counters;
    extendFMPos(getCompleteRange(), stack, counters);
    while (!stack.empty()) {
        auto curr = stack.back();
        stack.pop_back();

        word.resize(curr.getRow());
        word[curr.getRow() - 1] = curr.getCharacter();

        if ((length_t)curr.getRow() == wordSize) { // max depth reached
            Kmer k(word);
            table.insert(make_pair(k, curr.getRanges()));

        } else // add extra characters
            extendFMPos(curr.getRanges(), stack, counters, curr.getRow());
    }
    if (verbose) {
        cout << "done." << endl;
    }
}

// ----------------------------------------------------------------------------
// ROUTINES FOR EXACT PATTERN MATCHING
// ----------------------------------------------------------------------------
Range FMIndex::matchString(const string& s, Counters& counters) const {
    // start at the end
    auto it = s.crbegin();

    // find the range for this initial character in the BWT string
    length_t positionInAlphabet = sigma.c2i((unsigned char)*it);

    length_t start = counts[positionInAlphabet];
    length_t end;
    if (positionInAlphabet != sigma.size() - 1) {
        end = counts[positionInAlphabet + 1];
    } else {
        end = bwt.size();
    }
    counters.nodeCounter++;

    // iterate starting from the second character over the string
    for (++it; it != s.crend(); it++) {
        // find number of occurrences of this char before and after and so the
        // new range is found
        positionInAlphabet = sigma.c2i((unsigned char)*it);
        length_t startOfChar = counts[positionInAlphabet];
        start = getNumberOfOcc(positionInAlphabet, start) + startOfChar;
        end = getNumberOfOcc(positionInAlphabet, end) + startOfChar;
        counters.nodeCounter++;
        if (start == end) {
            // no matches found
            return Range();
        }
    }

    // return this range as a pair
    return Range(start, end);
}
vector<length_t> FMIndex::exactMatches(const string& s,
                                       Counters& counters) const {

    // find the range in the suffix array that matches the string
    Range range = matchString(s, counters);

    // declare the return vector
    vector<length_t> positions;
    positions.reserve(range.width());

    // fill in the vector with all values in this range in the suffix array
    for (length_t i = range.getBegin(); i < range.getEnd(); i++) {
        positions.emplace_back(findSA(i));
    }

    // sort the vector and return
    sort(positions.begin(), positions.end());
    return positions;
}

SARangePair FMIndex::matchStringBidirectionally(const Substring& pattern,
                                                SARangePair rangesOfPrev,
                                                Counters& counters) const {

    for (length_t i = 0; i < pattern.size(); i++) {

        char c = pattern[i];
        if (!addChar(c, rangesOfPrev, counters)) {
            // rangesOfPrev was made empty
            break;
        }
    }

    return rangesOfPrev;
}
bool FMIndex::addChar(const char& c, SARangePair& startRange,
                      Counters& counters) const {

    int posInAlphabet = sigma.c2i((unsigned char)c);
    if (posInAlphabet > -1) {

        if ((this->*extraChar)(posInAlphabet, startRange, startRange)) {
            // each character that we look at is a new node that is visited
            counters.nodeCounter++;
            return true;
        }
    }
    // the range is now empty

    return false;
}

// ----------------------------------------------------------------------------
// ROUTINES FOR APPROXIMATE PATTERN MATCHING
// ----------------------------------------------------------------------------

vector<TextOcc> FMIndex::approxMatchesNaive(const string& pattern,
                                            length_t maxED,
                                            Counters& counters) {

    counters.resetCounters();
    Occurrences occurrences;

    BandMatrix matrix(pattern.size() + maxED + 1, maxED);

    setDirection(FORWARD);

    vector<FMPosExt> stack;
    stack.reserve((pattern.size() + maxED + 1) * (sigma.size() - 1));

    extendFMPos(getCompleteRange(), stack, counters, 0);

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
                occurrences.addFMOcc(currentNode, matrix(row, last));
            }
        }

        extendFMPos(currentNode, stack, counters);
    }

    BitParallelED bpMatrix;
    bpMatrix.setSequence(pattern);

    return occurrences.getUniqueTextOccurrences(*this, maxED, bpMatrix,
                                                counters);
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
    // find the number of occurrences of chars smaller than c in the parent
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
    // find the number of occurrences of chars smaller than c in the parent
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
                                      Occurrences& occ,
                                      const vector<Substring>& parts,
                                      BitParallelED& inTextMatrix,
                                      Counters& counters, const int& idx) {
    const Substring& p = parts[s.getPart(idx)];           // this part
    const length_t& maxED = s.getUpperBound(idx);         // maxED for this part
    const length_t& minED = s.getLowerBound(idx);         // minED for this part
    const length_t& W = maxED - startMatch.getDistance(); // Width of matrix
    const length_t& pSize = p.size();

    BandMatrix matrix = BandMatrix(pSize, W, startMatch.getDistance());

    if (matrix.getLastColumn(0) == (int)pSize && matrix(0, pSize) <= maxED) {
        // an occurrence found by gapping entire part
        if (s.isEnd(idx)) {
            occ.addFMOcc(startMatch.getRanges(), matrix(0, pSize),
                         startMatch.getDepth());
        } else {
            // go to the next index
            recApproxMatchEditNaive(
                s,
                FMOcc(startMatch.getRanges(), matrix(0, pSize),
                      startMatch.getDepth()),
                occ, parts, inTextMatrix, counters, idx + 1);
            // set direction correct again
            setDirection(dir);
        }
    }

    auto& stack = stacks[idx];                  // stack for this partition
    const Direction& dir = s.getDirection(idx); // direction
    setDirection(dir);

    extendFMPos(startMatch.getRanges(), stack, counters);

    while (!stack.empty()) {
        const auto currentNode = stack.back();
        stack.pop_back();

        const length_t& row = currentNode.getRow();
        const length_t firstCol = matrix.getFirstColumn(row);
        const length_t lastCol = matrix.getLastColumn(row);

        length_t minimalEDOfRow = matrix(row, firstCol - 1);
        for (length_t col = firstCol; col <= lastCol; col++) {
            length_t filledIn = matrix.updateMatrix(
                currentNode.getCharacter() != p[col - 1], row, col);
            counters.abortedInTextVerificationCounter++;
            if (filledIn < minimalEDOfRow) {
                minimalEDOfRow = filledIn;
            }
        }

        if (minimalEDOfRow > maxED) {
            // backtracking
            continue;
        }

        if (lastCol == pSize && matrix(row, pSize) <= maxED &&
            matrix(row, pSize) >= minED) {
            if (s.isEnd(idx)) {
                occ.addFMOcc(currentNode.getRanges(), matrix(row, pSize),
                             startMatch.getDepth() + currentNode.getDepth());
            } else {
                // go deeper in search
                Direction originalDir = dir;
                recApproxMatchEditNaive(
                    s,
                    FMOcc(currentNode.getRanges(), matrix(row, pSize),
                          startMatch.getDepth() + currentNode.getDepth()),
                    occ, parts, inTextMatrix, counters, idx + 1);
                // set direction correct again
                setDirection(originalDir);
            }
        }
        if (firstCol == pSize) {
            // final cell of matrix filled in => backtracking
            continue;
        }

        // extend to the children
        extendFMPos(currentNode, stack, counters);
    }
}

void FMIndex::recApproxMatchEditOptimized(
    BitParallelED& intextMatrix, const Search& s, const FMOcc& startMatch,
    Occurrences& occ, const vector<Substring>& parts, Counters& counters,
    const int& idx, const vector<FMPosExt>& descPrevDir,
    const vector<uint>& initPrevDir, const vector<FMPosExt>& descNotPrevDir,
    const vector<uint>& initNotPrevDir) {

    // shortcut Variables
    const Substring& p = parts[s.getPart(idx)];      // this part
    const length_t& maxED = s.getUpperBound(idx);    // maxED for this part
    const Direction& dir = s.getDirection(idx);      // direction
    const bool& dSwitch = s.getDirectionSwitch(idx); // has direction switched?
    auto& stack = stacks[idx];                       // stack for this partition
    size_t matrixIdx = s.getPart(idx) + (dir == BACKWARD) * s.getNumParts();
    BitParallelED& bpED = matrices[matrixIdx]; // matrix for this partition

    // get the correct initED and descendants based on switch
    const vector<uint>& initEds = dSwitch ? initNotPrevDir : initPrevDir;
    const vector<FMPosExt>& descendants =
        dSwitch ? descNotPrevDir : descPrevDir;
    const vector<uint>& initOther = dSwitch ? initPrevDir : initNotPrevDir;
    const vector<FMPosExt>& descOther = dSwitch ? descPrevDir : descNotPrevDir;

    // set the direction
    setDirection(dir);

    // calculate necessary increase for first column of bandmatrix
    vector<uint> initED;
    if (initEds.empty()) {
        initED = vector<uint>(1, startMatch.getDistance());
    } else {
        uint prevED = (dSwitch ? *min_element(initEds.begin(), initEds.end())
                               : initEds[0]);
        uint increase = startMatch.getDistance() - prevED;
        initED = vector<uint>(initEds.size());
        for (size_t i = 0; i < initED.size(); i++) {
            initED[i] = initEds[i] + increase;
        }
    }

    // encode the sequence of this partition in the matrix if this has not been
    // done before
    if (!bpED.sequenceSet())
        bpED.setSequence(p);

    // initialize bit-parallel matrix
    bpED.initializeMatrix(maxED, initED);

    // initialize matrix and cluster
    Cluster clus(bpED.getSizeOfFinalColumn(), maxED, startMatch.getDepth(),
                 startMatch.getShift());

    if (bpED.inFinalColumn(0)) {
        // the first row is part of the final column of the banded matrix
        // Update the first cell of the cluster with the startmatch and the
        // value found at the psizeth column of the initialization row the
        // character does not matter as the first cell of a cluster is never
        // a descendant of any of the other cells in the cluster, so this
        // cell will not be reused for the next part of the pattern
        clus.setValue(0, FMPosExt((char)0, startMatch.getRanges(), 0),
                      bpED(0, p.size()));
    }

    if (!descendants.empty()) {
        // fill in the matrix for the descendants

        length_t maxRow = bpED.getNumberOfRows() - 1;

        for (length_t i = 0;
             i < descendants.size() && descendants[i].getDepth() <= maxRow;
             i++) {

            if (branchAndBound(
                    intextMatrix, clus, descendants[i], s, idx, parts, occ,
                    counters, initOther, descOther,
                    {descendants.begin() + i + 1, descendants.end()})) {
                return;
            }
        }
        if (descendants.back().getDepth() == maxRow) {
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
        extendFMPos(pair, stack, counters, descendants.back().getDepth());

    } else { // get the initial nodes to check
        extendFMPos(startMatch.getRanges(), stack, counters);
    }

    bool idxZero = idx == 0;

    while (!stack.empty()) {
        const FMPosExt currentNode = stack.back();
        stack.pop_back();

        if (branchAndBound(intextMatrix, clus, currentNode, s, idx, parts, occ,
                           counters, initOther, descOther)) {

            continue;
        }

        bool lastCol = bpED.inFinalColumn(currentNode.getRow());

        if (lastCol || currentNode.getRanges().width() > inTextSwitchPoint ||
            idxZero) {
            // continue the search for children of this node in-index
            extendFMPos(currentNode, stack, counters);
            continue;
        }
        // crossing-over
        // convert node to an occurrence in the text with the
        FMOcc fmocc(currentNode.getRanges(), 0,
                    startMatch.getDepth() + currentNode.getDepth(),
                    startMatch.getShift());
        const auto& textOcc = convertToMatchesInText(fmocc);

        // find the decrease and increase as compared to partialStart
        length_t lStartDec = 0, hStartDec = 0, hStartInc = 0;

        length_t startBeforeThis =
            parts[s.getLowestPartProcessedBefore(idx)].begin();

        (this->*findDiff)(lStartDec, hStartDec, hStartInc, startBeforeThis,
                          s.getMaxED(), bpED, currentNode.getRow(), maxED,
                          descOther.size(), initOther);

        inTextVerification(textOcc, s.getMaxED(), s.getMinED(), intextMatrix,
                           occ, counters, lStartDec, hStartDec, hStartInc);
    }
}

bool FMIndex::branchAndBound(BitParallelED& inTextMatrix, Cluster& clus,
                             const FMPosExt& currentNode, const Search& s,
                             const length_t& idx,
                             const vector<Substring>& parts, Occurrences& occ,
                             Counters& counters, const vector<uint>& initOther,
                             const vector<FMPosExt>& descOther,
                             const vector<FMPosExt>& remainingDesc) {
    // get the appropriate matrix
    size_t matrixIdx = s.getPart(idx) + (dir == BACKWARD) * s.getNumParts();
    BitParallelED& bpED = matrices[matrixIdx];

    // compute, in a bit-parallel manner, a single row of the ED matrix
    const length_t row = currentNode.getDepth();
    bool validED = bpED.computeRow(row, currentNode.getCharacter());

    // check if we have reached the final column of the matrix
    const length_t lastCol = bpED.getLastColumn(row);

    if ((lastCol + 1) == bpED.getNumberOfCols()) {
        // update the cluster
        length_t clusIdx = clus.size() + row - bpED.getNumberOfRows();
        clus.setValue(clusIdx, currentNode, bpED(row, lastCol));

        if (!validED || bpED.onlyVerticalGapsLeft(row)) {
            // no need to further explore this branch for this part -> go to
            // next part
            goDeeper(inTextMatrix, clus, idx + 1, s, parts, occ,
                     s.getLowerBound(idx), counters, descOther, initOther,
                     remainingDesc);
            return true;
        }
    }

    return !validED;
}

void FMIndex::goDeeper(BitParallelED& inTextMatrix, Cluster& cluster,
                       const length_t& nextIdx, const Search& s,
                       const vector<Substring>& parts, Occurrences& occ,
                       const length_t& lowerBound, Counters& counters,
                       const vector<FMPosExt>& descOtherD,
                       const vector<uint>& initOtherD,
                       const vector<FMPosExt>& remainingDesc) {

    bool isEdge = s.isEdge(nextIdx - 1);

    if (isEdge) {
        // if this is final piece report highest minimum (to get shortest
        // match)
        if (nextIdx == parts.size()) {
            auto matches = cluster.reportCentersAtEnd();
            for (const auto& match : matches) {
                if (match.isValid() && match.getDistance() >= lowerBound) {
                    occ.addFMOcc(match);
                }
            }
        } else {
            FMOcc match = cluster.reportDeepestMinimum(this->dir);
            if (match.isValid() && match.getDistance() >= lowerBound) {
                // go deeper in search
                Direction originalDir = this->dir;
                recApproxMatchEditOptimized(inTextMatrix, s, match, occ, parts,
                                            counters, nextIdx, {}, {},
                                            descOtherD, initOtherD);
                // set direction back again
                setDirection(originalDir);
            }
        }

        return;
    }

    // one of the later stages will return to this point, so keep track of
    // the descendants and eds at this branch
    vector<FMPosExt> descendants;
    vector<uint> initEds;

    FMOcc newMatch = cluster.getClusterCentra(lowerBound, descendants, initEds);
    if (!newMatch.isValid()) {
        // no centre above lower bound found
        return;
    }

    // add the remaining descendants
    descendants.insert(descendants.end(), remainingDesc.begin(),
                       remainingDesc.end());

    // reset the depth of all the descendants
    for (length_t i = 0; i < descendants.size(); i++) {
        descendants[i].setDepth(i + 1);
    }
    length_t maxEDNext = s.getUpperBound(nextIdx);

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

        recApproxMatchEditOptimized(inTextMatrix, s, newMatch, occ, parts,
                                    counters, nextIdx, descendants, initEds,
                                    descOtherD, initOtherD);

        // set direction back again
        setDirection(originalDir);
    } else {
        // go deeper on next piece
        recApproxMatchEditOptimized(inTextMatrix, s, newMatch, occ, parts,
                                    counters, nextIdx, descendants, initEds,
                                    descOtherD, initOtherD);
    }
}

void FMIndex::findDiffStartPositionBackward(
    length_t& lStartDec, length_t& hStartDec, length_t& hStartInc,
    const length_t& startBeforeThis, const length_t& maxED,
    const BitParallelED& bpED, const length_t& row, const length_t& maxEDPart,
    const length_t& descOtherSize, const vector<uint>& initOther) const {
    assert(dir == BACKWARD);
    // the cluster centers of the final row
    vector<pair<uint, uint>> centers;
    centers.reserve(2 * maxEDPart + 1);

    // find the local minima
    bpED.findLocalMinimaRow(row, maxEDPart, centers);
    lStartDec = 0, hStartDec = numeric_limits<length_t>::max();
    hStartInc = 0;
    for (const auto& c : centers) {
        length_t remPatternLeft = startBeforeThis - c.first;
        length_t remED = maxED - c.second;
        lStartDec = max(lStartDec, remPatternLeft + remED);
        if (remPatternLeft >= remED) {
            hStartDec = min(hStartDec, remPatternLeft - remED);
        } else {
            hStartInc = max(hStartInc, remED - remPatternLeft);
        }
    }
    if (hStartInc) {
        hStartDec = 0;
    }
}

void FMIndex::findDiffStartPositionForward(
    length_t& lStartDec, length_t& hStartDec, length_t& hStartInc,
    const length_t& startBeforeThis, const length_t& maxED,
    const BitParallelED& bpED, const length_t& row, const length_t& maxEDPart,
    const length_t& descOtherSize, const vector<uint>& initOther) const {
    assert(dir == FORWARD),

        lStartDec = 0, hStartDec = numeric_limits<length_t>::max();
    hStartInc = 0;

    if (startBeforeThis == 0) {
        // start position is already fixed
        hStartDec = 0;
        return;
    }
    // find the absolute minimum value of the final row
    length_t minPartialED, minIndex;
    bpED.findMinimumAtRow(row, minIndex, minPartialED);
    if (descOtherSize == 0) {
        // partial match has fixed start position
        // allow for extra characters and errors
        lStartDec = startBeforeThis + (maxED - minPartialED);
        hStartDec = startBeforeThis - (maxED - minPartialED);
        return;
    }
    // BACKWARDS already started and has influence on startpositions
    // we have started the backwards match, but not completed it
    // -> the partialStart refers to the deepest descendant

    length_t deepestMDec = startBeforeThis;

    for (length_t cc = 0; cc < initOther.size(); cc++) {
        bool betterThanAbove = (cc == 0) || initOther[cc] <= initOther[cc - 1];
        bool betterThanBelow =
            (cc == initOther.size() - 1) || initOther[cc] <= initOther[cc + 1];
        if (betterThanAbove && betterThanBelow) {
            // this is a center
            int centerDiff = initOther[cc] - initOther[0];
            length_t remED = maxED - minPartialED - centerDiff;
            length_t mDec = deepestMDec - (descOtherSize - cc);
            lStartDec = max(lStartDec, mDec + remED);
            hStartDec = min(hStartDec, mDec - remED);
        }
    }
}

void FMIndex::inTextVerification(const vector<TextOcc>& tos,
                                 const length_t& maxED, const length_t& minED,
                                 BitParallelED& intextMatrix, Occurrences& occ,
                                 Counters& counters, const length_t& lStartDec,
                                 const length_t& hStartDec,
                                 const length_t& hStartInc) const {

    // initialize matrix with correct number of zeros
    vector<uint> zeros(lStartDec - hStartDec + hStartInc + 1, 0);
    intextMatrix.initializeMatrix(maxED, zeros);

    counters.inTextStarted += tos.size();
    for (const auto& to : tos) {
        const length_t& partialStart = to.getRange().getBegin();
        // A) find the lowest  possible starts
        length_t lStart = partialStart - lStartDec;

        // B) find the highest possible end
        length_t hEnd = intextMatrix.getNumberOfRows() - 1 + lStart;
        // C) Get the reference subsequence
        Substring ref(&text, lStart, hEnd, FORWARD);

        // D) fill in the matrix row by row
        length_t i;

        for (i = 0; i < ref.size(); i++) {
            if (!intextMatrix.computeRow(i + 1, ref[i])) {
                break;
            }
        }

        // did we break before a possible match?
        if (i <= ref.size() - intextMatrix.getSizeOfFinalColumn()) {
            counters.abortedInTextVerificationCounter++;
            continue;
        }

        vector<uint> refEnds;
        intextMatrix.findClusterCenters(i, refEnds, maxED, minED);

        if (refEnds.empty()) {
            counters.abortedInTextVerificationCounter++;
            continue;
        }

        // for each valid end -> calculate CIGAR string and report
        for (const auto& refEnd : refEnds) {
            length_t bestScore = maxED + 1, bestBegin = 0;

            vector<pair<char, uint>> CIGAR;
            intextMatrix.trackBack(ref, refEnd, bestBegin, bestScore, CIGAR);
            counters.cigarsInTextVerification++;

            // make an occurrence
            occ.addTextOcc(Range(lStart + bestBegin, lStart + refEnd),
                           bestScore, CIGAR);
        }
    }
}

void FMIndex::recApproxMatchHamming(const Search& s, const FMOcc& startMatch,
                                    Occurrences& occ,
                                    const vector<Substring>& parts,
                                    Counters& counters, const int& idx) {

    // shortcut variables
    const Substring& p = parts[s.getPart(idx)]; // the current part
    const length_t& pSize = p.size();           // the size of the current part
    const Direction& d = s.getDirection(idx);   // direction of current part
    const length_t& maxED = s.getUpperBound(idx); // upper bound of current part
    const length_t& minED = s.getLowerBound(idx); // lower bound of currentpart
    setDirection(d);

    // create vector for the scores
    vector<length_t> vec(p.size() + 1, 0);
    // set root element of vector to the distance of the startmatch
    vec[0] = startMatch.getDistance();
    // get stack for current part
    auto& stack = stacks[idx];

    extendFMPos(startMatch.getRanges(), stack, counters);

    while (!stack.empty()) {
        const FMPosExt node = stack.back();
        stack.pop_back();

        if (node.getRanges().width() <= inTextSwitchPoint) {
            inTextVerificationHamming(node, s, parts, idx, occ);
            continue;
        }

        // update the vector
        length_t row = node.getRow();
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
                    occ.addFMOcc(match);
                } else {
                    // continue search
                    recApproxMatchHamming(s, match, occ, parts, counters,
                                          idx + 1);
                    setDirection(s.getDirection(idx));
                }
            }
            continue;
        }
        extendFMPos(node, stack, counters);
    }
}

void FMIndex::inTextVerificationHamming(const FMPosExt& node, const Search& s,
                                        const vector<Substring>& parts,
                                        const length_t idx,
                                        Occurrences& occ) const {

    // A) find length before and the partial occurrence
    length_t lengthBefore =
        ((idx == 0) ? 0 : parts[s.getLowestPartProcessedBefore(idx)].begin()) -
        (dir == BACKWARD) * (node.getDepth());
    length_t pSize = parts.back().end();

    length_t maxEDFull = s.getMaxED();
    length_t minEDFull = s.getMinED();

    const Range& r = node.getRanges().getRangeSA(); // SA range of current node

    for (length_t i = r.getBegin(); i < r.getEnd(); i++) {

        // look up this position in the suffix array
        length_t Tb = findSA(i);
        // subtract the length before
        Tb = (Tb > lengthBefore) ? Tb - lengthBefore : 0;

        // Calculate the end in the text  + guard that is does not go over the
        // textlength
        length_t Te = min(textLength, Tb + pSize);

        // Create the reference and pattern sequence
        Substring ref(&text, Tb, Te, FORWARD);
        Substring pattern(parts[0], 0, pSize, FORWARD);

        assert(ref.size() == pattern.size());

        length_t score = 0;

        for (length_t i = 0; i < ref.size(); i++) {
            // update the score
            score = score + (ref[i] != pattern[i]);
            if (score > maxEDFull) {
                // in text verification failed
                break;
            }
        }
        if (score <= maxEDFull && score >= minEDFull) {
            // in text verification succeeded
            vector<pair<char, uint>> CIGAR = {make_pair('M', pSize)};
            occ.addTextOcc(Range(Tb, Te), score, CIGAR);
        }
    }
}

void FMIndex::extendFMPos(const SARangePair& parentRanges,
                          vector<FMPosExt>& stack, Counters& counters,
                          length_t row) const {

    // iterate over the entire alphabet
    for (length_t i = 1; i < sigma.size(); i++) {

        SARangePair pairForNewChar;

        // check if this character occurs in the specified range
        if ((this->*extraChar)(i, parentRanges, pairForNewChar)) {
            // push this range and character for the next iteration
            stack.emplace_back(sigma.i2c(i), pairForNewChar, row + 1);

            counters.nodeCounter++;
        }
    }
}

void FMIndex::extendFMPos(const FMPosExt& pos, vector<FMPosExt>& stack,
                          Counters& counters) const {
    extendFMPos(pos.getRanges(), stack, counters, pos.getDepth());
}

// ----------------------------------------------------------------------------
// POST-PROCESSING ROUTINES FOR APPROXIMATE PATTERN MATCHING
// ----------------------------------------------------------------------------

vector<TextOcc> FMIndex::convertToMatchesInText(const FMOcc& saMatch) const {

    vector<TextOcc> textMatches;
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
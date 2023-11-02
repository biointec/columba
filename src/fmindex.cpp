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
#include "fmindex.h"

using namespace std;

// ============================================================================
// CLASS FMIndex
// ============================================================================

thread_local Direction FMIndex::dir = BACKWARD;
thread_local ExtraCharPtr FMIndex::extraChar;

thread_local vector<vector<FMPosExt>> FMIndex::stacks;
thread_local vector<BitParallelED> FMIndex::matrices;

thread_local BitParallelED FMIndex::inTextMatrix;

// ----------------------------------------------------------------------------
// ROUTINES FOR ACCESSING DATA STRUCTURE
// ----------------------------------------------------------------------------
length_t FMIndex::findLF(length_t k) const {

    const auto& pos = bwt[k];
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

    if (!bwt.read(baseFile + ".bwt")) {
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

        if (curr.getRow() == wordSize) { // max depth reached
            Kmer k(word);

            table.insert(make_pair(k, std::move(curr.getRanges())));

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
            counters.incNodeCounter();
            return true;
        }
    }
    // the range is now empty
    return false;
}

// ----------------------------------------------------------------------------
// ROUTINES FOR APPROXIMATE PATTERN MATCHING
// ----------------------------------------------------------------------------

const vector<TextOcc> FMIndex::approxMatchesNaive(const string& pattern,
                                                  length_t maxED,
                                                  Counters& counters) {

    counters.resetCounters();
    Occurrences occurrences;

    BitParallelED matrix;
    matrix.setSequence(pattern);
    matrix.initializeMatrix(maxED);

    setDirection(FORWARD);

    vector<FMPosExt> stack;
    stack.reserve((pattern.size() + maxED + 1) * (sigma.size() - 1));

    extendFMPos(getCompleteRange(), stack, counters, 0);

    length_t lastcol = pattern.size();

    while (!stack.empty()) {
        const FMPosExt currentNode = stack.back();
        stack.pop_back();
        length_t row = currentNode.getDepth();

        if (row >= matrix.getNumberOfRows()) {
            continue;
        }

        bool valid = matrix.computeRow(row, currentNode.getCharacter());

        if (!valid) {
            // backtrack
            continue;
        }

        if (matrix.inFinalColumn(row)) {

            // full pattern was matched
            if (matrix(row, lastcol) <= maxED) {
                occurrences.addFMOcc(currentNode, matrix(row, lastcol));
            }
        }

        extendFMPos(currentNode, stack, counters);
    }

    setInTextMatrixSequence(pattern);

    return getUniqueTextOccurrences(occurrences, maxED, counters);
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

                                      Counters& counters, const int& idx) {
    const Substring& p = parts[s.getPart(idx)];   // this part
    const length_t& maxED = s.getUpperBound(idx); // maxED for this part
    const length_t& minED = s.getLowerBound(idx); // minED for this part
    const length_t& pSize = p.size();
    const Direction& dir = s.getDirection(idx); // direction

    // set the direction
    setDirection(dir);

    length_t matrixID = s.getPart(idx) + (dir == BACKWARD) * s.getNumParts();
    auto& bpED = matrices[matrixID];

    if (!bpED.sequenceSet()) {
        bpED.setSequence(p);
    }

    bpED.initializeMatrix(maxED, {(uint)startMatch.getDistance()});

    if (bpED.getLastColumn(0) == pSize && bpED(0, pSize) <= maxED) {
        // an occurrence found by gapping entire part
        if (s.isEnd(idx)) {
            occ.addFMOcc(startMatch.getRanges(), bpED(0, pSize),
                         startMatch.getDepth());
        } else {
            // go to the next index
            recApproxMatchEditNaive(s,
                                    FMOcc(startMatch.getRanges(),
                                          bpED(0, pSize),
                                          startMatch.getDepth()),
                                    occ, parts, counters, idx + 1);
            // set direction correct again
            setDirection(dir);
        }
    }

    auto& stack = stacks[idx]; // stack for this partition

    extendFMPos(startMatch.getRanges(), stack, counters);

    while (!stack.empty()) {
        const auto currentNode = stack.back();
        stack.pop_back();

        const length_t& row = currentNode.getRow();
        const length_t firstCol = bpED.getFirstColumn(row);
        const length_t lastCol = bpED.getLastColumn(row);

        bool valid = bpED.computeRow(row, currentNode.getCharacter());

        if (!valid) {
            // backtracking
            continue;
        }

        if (lastCol == pSize && bpED(row, pSize) <= maxED &&
            bpED(row, pSize) >= minED) {
            if (s.isEnd(idx)) {
                occ.addFMOcc(currentNode.getRanges(), bpED(row, pSize),
                             startMatch.getDepth() + currentNode.getDepth());
            } else {
                // go deeper in search
                Direction originalDir = dir;
                recApproxMatchEditNaive(
                    s,
                    FMOcc(currentNode.getRanges(), bpED(row, pSize),
                          startMatch.getDepth() + currentNode.getDepth()),
                    occ, parts, counters, idx + 1);
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
    const Search& s, const FMOcc& startMatch, Occurrences& occ,
    const vector<Substring>& parts, Counters& counters, const int& idx,
    const vector<FMPosExt>& descPrevDir, const vector<uint>& initPrevDir,
    const vector<FMPosExt>& descNotPrevDir,
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
        length_t prevED =
            (dSwitch ? *min_element(initEds.begin(), initEds.end())
                     : initEds[0]);
        length_t increase = startMatch.getDistance() - prevED;
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

    // initialize  cluster
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
                    bpED, clus, descendants[i], s, idx, parts, occ, counters,
                    initOther, descOther,
                    {descendants.begin() + i + 1, descendants.end()})) {
                return;
            }
        }
        if (descendants.back().getDepth() == maxRow) {
            //  no more rows to possibly check
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

        if (branchAndBound(bpED, clus, currentNode, s, idx, parts, occ,
                           counters, initOther, descOther)) {

            continue;
        }

        if (currentNode.getRanges().width() > inTextSwitchPoint || idxZero) {
            // continue the search for children of this node in-index
            extendFMPos(currentNode, stack, counters);
            continue;
        }
        // crossing-over
        // convert node to an occurrence in the text with the

        length_t startBeforeThis =
            parts[s.getLowestPartProcessedBefore(idx)].begin();
        length_t globalMaxED = s.getMaxED();
        length_t startDiff = startBeforeThis + globalMaxED;

        if (startBeforeThis == 0) {
            startDiff = 0;
        } else if (dir == BACKWARD) {
            length_t row = currentNode.getRow();
            length_t col = bpED.getFirstColumn(row);
            startDiff -= col + bpED(row, col);

        } else if (!descOther.empty()) {
            // FORWARD DIRECTION
            // descOther are thus in backwards direction
            startDiff -= descOther.size() - initOther.size() + initOther.back();
        }

        const auto& startPos =
            getBeginPositions(currentNode.getRanges().getRangeSA(), startDiff,
                              startMatch.getShift());

        inTextVerification(startPos, globalMaxED, s.getMinED(), occ, counters);
    }
}

bool FMIndex::branchAndBound(BitParallelED& bpED, Cluster& clus,
                             const FMPosExt& currentNode, const Search& s,
                             const length_t& idx,
                             const vector<Substring>& parts, Occurrences& occ,
                             Counters& counters, const vector<uint>& initOther,
                             const vector<FMPosExt>& descOther,
                             const vector<FMPosExt>& remainingDesc) {

    // compute, in a bit-parallel manner, a single row of the ED matrix
    const length_t row = currentNode.getDepth();
    bool validED = bpED.computeRow(row, currentNode.getCharacter());

    // check if we have reached the final column of the matrix
    if (bpED.inFinalColumn(row)) {
        // update the cluster
        length_t clusIdx = clus.size() + row - bpED.getNumberOfRows();
        clus.setValue(clusIdx, currentNode,
                      bpED(row, bpED.getNumberOfCols() - 1));

        if (!validED || bpED.onlyVerticalGapsLeft(row)) {
            // no need to further explore this branch for this part -> go to
            // next part

            goDeeper(clus, idx + 1, s, parts, occ, counters, descOther,
                     initOther, remainingDesc);
            return true;
        }
    }

    return !validED;
}

void FMIndex::goDeeper(Cluster& cluster, const length_t& nIdx, const Search& s,
                       const vector<Substring>& parts, Occurrences& occ,
                       Counters& counters, const vector<FMPosExt>& descOtherD,
                       const vector<uint>& initOtherD,
                       const vector<FMPosExt>& remDesc) {

    bool isEdge = s.isEdge(nIdx - 1);
    const auto& lowerBound = s.getLowerBound(nIdx - 1);

    if (isEdge) {
        // if this is final piece report highest minimum (to get shortest
        // match)
        if (nIdx == parts.size()) {
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
                recApproxMatchEditOptimized(s, match, occ, parts, counters,
                                            nIdx, {}, {}, descOtherD,
                                            initOtherD);
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

    // add the remaining descendants (copy)
    descendants.insert(descendants.end(), remDesc.begin(), remDesc.end());

    // reset the depth of all the descendants
    for (length_t i = 0; i < descendants.size(); i++) {
        descendants[i].setDepth(i + 1);
    }
    length_t maxEDNext = s.getUpperBound(nIdx);

    // remove trailing initEds that are higher than maxEDNext
    while (initEds.back() > maxEDNext) {
        initEds.pop_back();
    }

    // is the next direction equal to this direction?
    bool switchAfter = s.getDirectionSwitch(nIdx);

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

        recApproxMatchEditOptimized(s, newMatch, occ, parts, counters, nIdx,
                                    descendants, initEds, descOtherD,
                                    initOtherD);

        // set direction back again
        setDirection(originalDir);
    } else {
        // go deeper on next piece
        recApproxMatchEditOptimized(s, newMatch, occ, parts, counters, nIdx,
                                    descendants, initEds, descOtherD,
                                    initOtherD);
    }
}

void FMIndex::inTextVerification(const vector<length_t>& startPos,
                                 const length_t& maxED, const length_t& minED,
                                 Occurrences& occ, Counters& counters) const {

    // initialize matrix with correct number of zeros
    vector<uint> zeros(2 * maxED + 1, 0);
    inTextMatrix.initializeMatrix(maxED, zeros);

    vector<Substring> refs;
    refs.reserve(startPos.size());

    counters.inTextStarted += startPos.size();
    for (const auto& start : startPos) {
        // B) find the highest possible end
        length_t hEnd = start + inTextMatrix.getNumberOfRows() - 1;
        // C) Get the reference subsequence
        refs.emplace_back(&text, start, hEnd, FORWARD);
    }
    IntextVerificationTask task(refs, inTextMatrix, maxED, minED);
    task.doTask(counters, occ);
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

            counters.incNodeCounter();
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

std::vector<TextOcc> FMIndex::getUniqueTextOccurrences(Occurrences& occ,
                                                       const length_t& maxED,
                                                       Counters& counters) {

    // increment reporte position counter
    counters.totalReportedPositions += occ.textOccSize();

    // erase equal occurrences from the in-index occurrences
    occ.eraseDoublesFM();

    // convert the in-index occurrences to in-text occurrences
    const auto& fmocc = occ.getFMOccurrences();
    for (const auto& f : fmocc) {

        const Range& saRange = f.getRanges().getRangeSA();

        // increment reported positions counter
        counters.totalReportedPositions += saRange.width();

        auto depth = f.getDepth(), distance = f.getDistance();
        auto end = saRange.getEnd();
        length_t shift = f.getShift();

        for (length_t i = saRange.getBegin(); i < end; i++) {
            // find the startPosition in the text by looking at the
            // SA
            length_t startPos = findSA(i) + shift;

            occ.addTextOcc(Range(startPos, startPos + depth), distance);
        }
    }
    // erase equal occurrences from the in-text occurrences, note
    // that an in-text occurrence with calculated CIGAR string takes
    // preference over an equal one without CIGAR string
    occ.eraseDoublesText();

    // find the non-redundant occurrences
    std::vector<TextOcc> nonRedundantOcc;
    nonRedundantOcc.reserve(occ.textOccSize());

    length_t maxDiff = 2 * maxED;
    length_t prevBegin = std::numeric_limits<length_t>::max();
    length_t prevDepth = std::numeric_limits<length_t>::max();
    length_t prevED = maxED + 1;

    const auto& textocc = occ.getTextOccurrences();
    for (const auto& o : textocc) {
        // find the difference between this and the previous
        // occurrence
        auto diff = abs_diff<length_t>(o.getRange().getBegin(), prevBegin);

        if (diff == 0) {
            // same location -> skip
            continue;
        }

        if (diff <= maxDiff) {
            // check if this later occurrence is better than the
            // previous one
            if (o.getDistance() > prevED ||
                (o.getDistance() == prevED &&
                 o.getRange().width() >= prevDepth)) {
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

    for (TextOcc& t : nonRedundantOcc) {

        if (!t.hasCigar()) {

            // this was an in-index occurrence which needs to
            // calculate the CIGAR string
            std::vector<std::pair<char, uint>> CIGAR;
            inTextMatrix.findCIGAR(getSubstring(t.getRange()), t.getDistance(),
                                   CIGAR);

            t.setCigar(CIGAR);

            counters.cigarsInIndex++;

        } else {
            // this was a useful in-text cigar
            counters.usefulCigarsInText++;
        }
    }

    return nonRedundantOcc;
}

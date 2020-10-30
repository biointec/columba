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
#include "bidirecfmindex.h"

using namespace std;

// ----------------------------------------------------------------------------
// ROUTINES FOR ACCESSING DATA STRUCTURE
// ----------------------------------------------------------------------------

length_t BidirecFMIndex::findLF(length_t k, bool reversed) const {
    if (!reversed) {

        return counts[sigma.c2i((unsigned char)bwt[k])] +
               getNumberOfOcc(sigma.c2i((unsigned char)bwt[k]), k);
    }
    length_t posInAlphabet = sigma.c2i((unsigned char)bwt[k]);

    return counts[posInAlphabet] + getNumberOfOccRev(posInAlphabet, k);
}

// ----------------------------------------------------------------------------
// ROUTINES FOR INITIALIZATION
// ----------------------------------------------------------------------------

void BidirecFMIndex::populateTable() {
    cout << "Populating FM-range table with " << wordSize << "-mers...";
    cout.flush();

    Kmer::setWordSize(wordSize);

    table.resize(1 << (2 * wordSize)); // 2 << wordSize is 4^wordSize
    setDirection(FORWARD);

    string word;
    vector<Node> stack;
    pushChildren(getCompleteRange(), stack);
    while (!stack.empty()) {
        auto curr = stack.back();
        stack.pop_back();

        word.resize(curr.row);
        word[curr.row - 1] = curr.c;

        if ((length_t)curr.row == wordSize) { // max depth reached
            Kmer k(word);
            table.insert(make_pair(k, curr.ranges));

        } else // add extra characters
            pushChildren(curr.getRanges(), stack, curr.getRow());
    }

    cout << "done." << endl;
}

// ----------------------------------------------------------------------------
// ROUTINES FOR EXACT PATTERN MATCHING
// ----------------------------------------------------------------------------

SARangePair
BidirecFMIndex::matchStringBidirectionally(const Substring& pattern,
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
bool BidirecFMIndex::addChar(const char& c, SARangePair& startRange) {

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

bool BidirecFMIndex::findRangesWithExtraCharBackward(
    length_t positionInAlphabet, const SARangePair& rangesOfP,
    SARangePair& rangesOfChild) const {

    // first make the trivial range by  searching cP using B
    Range trivialRange = rangesOfP.rangeSA;

    // find the new range by using the LF property
    length_t occBefore = getNumberOfOcc(positionInAlphabet, trivialRange.begin);
    length_t occAfter = getNumberOfOcc(positionInAlphabet, trivialRange.end);

    length_t startInAlphabet = counts[positionInAlphabet];
    Range range1 =
        Range(occBefore + startInAlphabet, occAfter + startInAlphabet);

    // then make the less trivial range by counting the sizes of the ranges
    // of (dP) using B

    // first get the start of the range we are looking for of the parent
    length_t s = rangesOfP.rangeSARev.begin;

    // get the start of the child within this range
    // find the number of occurences of chars smaller than c in the parent
    // range

    length_t x = getNumberOfCumOcc(positionInAlphabet, trivialRange.end) -
                 getNumberOfCumOcc(positionInAlphabet, trivialRange.begin);

    // make the new range with width equal to that of the trivial range
    Range range2 = Range(s + x, s + x + range1.width());

    rangesOfChild = SARangePair(range1, range2);
    return !rangesOfChild.empty();
}

bool BidirecFMIndex::findRangesWithExtraCharForward(
    length_t positionInAlphabet, const SARangePair& rangesOfP,
    SARangePair& childRanges) const {

    // first make the trivial range by searching (Pc)' using B'  if
    // searching forward we need to use B' so we need the reverse range
    Range rangeForTrivial = rangesOfP.rangeSARev;

    // find the new range by using the LF property
    length_t occBefore =
        getNumberOfOccRev(positionInAlphabet, rangeForTrivial.begin);
    length_t occAfter =
        getNumberOfOccRev(positionInAlphabet, rangeForTrivial.end);

    length_t startInAlphabet = counts[positionInAlphabet];
    Range range1 =
        Range(occBefore + startInAlphabet, occAfter + startInAlphabet);

    // then make the less trivial range by counting the size of the range of
    // (Pd)' using B' (forward)

    // first get the start of the range we are looking for of the parent
    length_t s = rangesOfP.rangeSA.begin;

    // get the start of the child within this range
    // find the number of occurences of chars smaller than c in the parent
    // range
    Range prevRange = rangesOfP.rangeSARev;

    length_t x = getNumberOfCumOccRev(positionInAlphabet, prevRange.end) -
                 getNumberOfCumOccRev(positionInAlphabet, prevRange.begin);

    // make the new range
    Range range2 = Range(s + x, s + x + range1.width());

    childRanges = SARangePair(range2, range1);
    return !childRanges.empty();
}

void BidirecFMIndex::recApproxMatch(
    const Search& s, const BiAppMatchSA& startMatch, vector<BiAppMatchSA>& occ,
    const vector<Substring>& parts, const int& idx,
    const vector<Node>& descPrevDir, const vector<int>& initPrevDir,
    const vector<Node>& descNotPrevDir, const vector<int>& initNotPrevDir) {

    // shortcut Variables
    const int& pNumber = s.order[idx];               // number of this part
    const length_t& maxED = s.upperBounds[idx];      // maxED for this part
    const length_t& W = maxED - startMatch.editDist; // Width of matrix
    const Substring& p = parts[pNumber];             // this part
    const Direction& dir = s.directions[idx];        // direction

    const bool& dSwitch = s.directionSwitch[idx]; // has direction switched?

    auto& stack = stacks[idx]; // stack for this partition
    // get the correct initED and descendants based on switch
    const vector<int>& initEds = dSwitch ? initNotPrevDir : initPrevDir;
    const vector<Node>& descendants = dSwitch ? descNotPrevDir : descPrevDir;
    const vector<int>& initOther = dSwitch ? initPrevDir : initNotPrevDir;
    const vector<Node>& descOther = dSwitch ? descPrevDir : descNotPrevDir;

    setDirection(dir);

    // calculate nescessary increase for first column of bandmatrix
    int prevED = 0;
    if (!initEds.empty()) {
        prevED = (dSwitch ? *min_element(initEds.begin(), initEds.end())
                          : initEds[0]);
    }
    int increase = startMatch.editDist - prevED;

    // initialize matrix and cluster
    BandMatrix matrix(p.size(), W, increase, initEds);
    Cluster clus(matrix.getSizeOfFinalColumn(), maxED, startMatch.depth);

    if (!descendants.empty()) {
        // fill in the matrix for the descendants

        int row = 1; // the current row of the matrix
        int maxRow = matrix.getNumberOfRows();
        int b = 0;
        for (length_t i = 0; i < descendants.size() && row <= maxRow;
             i++, row++) {
            b = branchAndBound(matrix, clus, row, descendants[i], s, idx, parts,
                               maxED, startMatch, occ, initOther, descOther);
        }
        if (b || row == maxRow) {
            // minimal ed exceeded no more options to get lower ed from first
            // column, or no more rows to possibly check
            return;
        }

        SARangePair pair = descendants.back().getRanges();
        if (dSwitch) {
            // after a switch the ranges of final descendant should be updated
            pair = startMatch.ranges;
        }

        // push children of final descendant
        pushChildren(pair, stack, row - 1);

    } else { // get the initial nodes to check
        pushChildren(startMatch.ranges, stack);
    }

    while (!stack.empty()) {

        const Node currentNode = stack.back();
        stack.pop_back();

        if (branchAndBound(matrix, clus, currentNode.getRow(), currentNode, s,
                           idx, parts, maxED, startMatch, occ, initOther,
                           descOther)) {
            continue;
        }

        // continue the search for children of this node
        pushChildren(currentNode.getRanges(), stack, currentNode.getRow());
    }
}

bool BidirecFMIndex::branchAndBound(
    BandMatrix& matrix, Cluster& clus, const int row, const Node& currentNode,
    const Search& s, const int& idx, const vector<Substring>& parts,
    const length_t& maxED, const BiAppMatchSA& startMatch,
    vector<BiAppMatchSA>& occ, const vector<int>& initOther,
    const vector<Node>& descOther) {
    const int& pNumber = s.order[idx];
    const Substring& p = parts[pNumber];
    const length_t& pSize = p.size();
    const int nextP = idx + 1;
    const length_t& lowerBound = s.lowerBounds[idx]; // lowerbound for this part
    // get first and last column of matrix for this row
    const length_t firstCol = matrix.getFirstColumn(row);
    const length_t lastCol = matrix.getLastColumn(row);

    // only necessary in lastCol == pSize
    bool afterwardsOnlyVerticalGap = true;

    unsigned int minimalEDOfRow = matrix(row, firstCol - 1);
    for (length_t j = firstCol; j <= lastCol; j++) {
        matrix.updateMatrix(p[j - 1] != currentNode.getCharacter(), row, j);
        matrixElementCounter++;

        if (minimalEDOfRow <= matrix(row, j)) {
            afterwardsOnlyVerticalGap = false;
        }

        if (lastCol == j && minimalEDOfRow > maxED) {
            afterwardsOnlyVerticalGap = true;
        }
        minimalEDOfRow = min(minimalEDOfRow, matrix(row, j));
    }

    if (lastCol == pSize) {
        // update the cluster
        int clusIdx = row - matrix.getNumberOfRows() + clus.size() - 1;
        clus.setValue(clusIdx, currentNode, matrix(row, lastCol));

        if (minimalEDOfRow > maxED || afterwardsOnlyVerticalGap) {
            // no need to further check nodes for this part
            // still need to check next part
            goDeeper(clus, maxED, startMatch, nextP, s, parts, occ, lowerBound,
                     descOther, initOther);
            return true;
        }
    }

    // edit distance threshold exceeded: backtracking
    if (minimalEDOfRow > maxED) {
        return -1;
    }
    // if at final row we can check if a match can be found to go to
    // next stage
    if (row == matrix.getNumberOfRows()) {
        // try to go deeper on this point

        goDeeper(clus, maxED, startMatch, nextP, s, parts, occ, lowerBound,
                 descOther, initOther);
        return true;
    }

    return false;
}

void BidirecFMIndex::goDeeper(
    Cluster& cluster, const length_t& maxED, const BiAppMatchSA& startMatch,
    const length_t& nextP, const Search& s, const vector<Substring>& parts,
    vector<BiAppMatchSA>& occ, const length_t& lowerBound,
    const vector<Node>& descOtherD, const vector<int>& initOtherD) {

    bool endOfDirection = s.order[nextP - 1] == (int)s.order.size() - 1 ||
                          s.order[nextP - 1] == 0;

    if (endOfDirection) {
        // if this is final piece report highest minimum (to get shortest
        // match) if not final piece (so we will return to here),
        // report deepest minimum

        if (nextP == parts.size()) {

            BiAppMatchSA newMatch = cluster.reportHighestMinimum();
            if (!newMatch.ranges.empty() && newMatch.editDist >= lowerBound) {
                occ.push_back(newMatch);
            }

        } else {
            BiAppMatchSA newMatch = cluster.reportDeepestMinimum();

            if (!newMatch.ranges.empty() && newMatch.editDist >= lowerBound) {

                Direction originalDir = this->dir;
                recApproxMatch(s, newMatch, occ, parts, nextP, {}, {},
                               descOtherD, initOtherD);
                // set direction back again
                setDirection(originalDir);
            }
        }
        return;
    }

    // this is not the end of the direction!
    vector<Node> descendants;
    vector<int> initEds;

    BiAppMatchSA newMatch =
        cluster.getClusterCentra(lowerBound, descendants, initEds);
    if (newMatch.ranges.empty()) {
        // no centre above lowerbound found
        return;
    }
    int maxEDNext = s.upperBounds[nextP];

    // remove trailing initEds that are higher than maxEDNext
    while (initEds.back() > maxEDNext) {
        initEds.pop_back();
    }

    // is the next direction equal to this direction?
    bool switchAfter = s.directions[nextP] != dir;

    if (switchAfter) {
        // switching direction as this is not the end of a search direction,
        // this means we'll get back here, thus range of newmatch should be
        // deepest point in branch
        if (!descendants.empty()) {
            newMatch.ranges = descendants.back().ranges;

            //  edit distance for search in other direction should be lowest
            //  value possible
            newMatch.editDist = *min_element(initEds.begin(), initEds.end());
        }

        Direction originalDir = this->dir;

        recApproxMatch(s, newMatch, occ, parts, nextP, descendants, initEds,
                       descOtherD, initOtherD);

        // set direction back again
        setDirection(originalDir);
    } else {
        // go deeper on next piece
        recApproxMatch(s, newMatch, occ, parts, nextP, descendants, initEds,
                       descOtherD, initOtherD);
    }
}

void BidirecFMIndex::pushChildren(const SARangePair& parentRanges,
                                  vector<Node>& stack, int row) {

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

vector<AppMatch> BidirecFMIndex::mapOccurencesInSAToOccurencesInText(
    vector<BiAppMatchSA>& occurences, const int& maxED) {

    sort(occurences.begin(), occurences.end());
    occurences.erase(unique(occurences.begin(), occurences.end()),
                     occurences.end());
    vector<AppMatchSA> newOccurences = discardReverseRanges(occurences);
    return FMIndex::mapOccurencesInSAToOccurencesInText(newOccurences, maxED);
}

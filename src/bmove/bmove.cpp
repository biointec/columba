/******************************************************************************
 *  Columba: Approximate Pattern Matching using Search Schemes                *
 *  Copyright (C) 2020-2024 - Lore Depuydt <lore.depuydt@ugent.be> and        *
 *                            Luca Renders <luca.renders@ugent.be> and        *
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
 *  You should have received a copy of the GNU Affero General Public License  *
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.    *
 ******************************************************************************/

#include "bmove.h"
#include "../alphabet.h" // for Alphabet
#include "../definitions.h"
#include "../indexhelpers.h"
#include "../logger.h" // for Logger
#include "../reads.h"

#include <algorithm>           // for max
#include <assert.h>            // for assert
#include <cstdint>             // for uint8_t
#include <ostream>             // for opera...
#include <sdsl/int_vector.hpp> // for int_v...
#include <stdexcept>           // for runti...
class MemoryMappedTextFile;
class Search;
class Substring;

using namespace std;

// ----------------------------------------------------------------------------
// PREPROCESSING ROUTINES
// ----------------------------------------------------------------------------

void BMove::fromFiles(const string& baseFile, bool verbose) {

    readMetaAndCounts(baseFile, verbose);
    stringstream ss;
    // Read BMove specific files
    if (verbose) {
        ss << "Reading " << baseFile << ".plcp...";
        logger.logInfo(ss);
    }
    if (!plcp.read(baseFile + ".plcp")) {
        throw runtime_error("Cannot open file: " + baseFile + ".plcp");
    }

    if (verbose) {
        ss << "Reading " << baseFile << ".move...";
        logger.logInfo(ss);
    }
    if (!move.load(baseFile, verbose)) {
        throw runtime_error("Error loading move file: " + baseFile + ".move");
    }

    if (verbose) {
        ss << "Reading " << baseFile << ".smpf...";
        logger.logInfo(ss);
    }
    if (!readIntVector(baseFile + ".smpf", samplesFirst)) {
        throw runtime_error("Cannot open file: " + baseFile + ".smpf");
    }

    if (verbose) {
        ss << "Reading " << baseFile << ".smpl...";
        logger.logInfo(ss);
    }
    if (!readIntVector(baseFile + ".smpl", samplesLast)) {
        throw runtime_error("Cannot open file: " + baseFile + ".smpl");
    }

    if (verbose) {
        ss << "Reading " << baseFile << ".prdf...";
        logger.logInfo(ss);
    }
    if (!predFirst.read(baseFile + ".prdf")) {
        throw runtime_error("Cannot open file: " + baseFile + ".prdf");
    }

    if (verbose) {
        ss << "Reading " << baseFile << ".prdl...";
        logger.logInfo(ss);
    }
    if (!predLast.read(baseFile + ".prdl")) {
        throw runtime_error("Cannot open file: " + baseFile + ".prdl");
    }

    if (verbose) {
        ss << "Reading " << baseFile << ".ftr...";
        logger.logInfo(ss);
    }
    if (!readIntVector(baseFile + ".ftr", firstToRun)) {
        throw runtime_error("Cannot open file: " + baseFile + ".ftr");
    }

    if (verbose) {
        ss << "Reading " << baseFile << ".ltr...";
        logger.logInfo(ss);
    }
    if (!readIntVector(baseFile + ".ltr", lastToRun)) {
        throw runtime_error("Cannot open file: " + baseFile + ".ltr");
    }

    if (verbose) {
        ss << "Reading " << baseFile << ".rev.move...";
        logger.logInfo(ss);
    }
    if (!moveR.load(baseFile + ".rev", verbose)) {
        throw runtime_error("Error loading reverse move file: " + baseFile +
                            ".rev.move");
    }

    if (verbose) {
        ss << "Reading " << baseFile << ".rev.smpf...";
        logger.logInfo(ss);
    }
    if (!readIntVector(baseFile + ".rev.smpf", revSamplesFirst)) {
        throw runtime_error("Cannot open file: " + baseFile + ".rev.smpf");
    }

    if (verbose) {
        ss << "Reading " << baseFile << ".rev.smpl...";
        logger.logInfo(ss);
    }
    if (!readIntVector(baseFile + ".rev.smpl", revSamplesLast)) {
        throw runtime_error("Cannot open file: " + baseFile + ".rev.smpl");
    }

    readSequenceNamesAndPositions(baseFile, verbose);
}

// ----------------------------------------------------------------------------
// ROUTINES FOR ACCESSING DATA STRUCTURE
// ----------------------------------------------------------------------------

length_t BMove::getSwitchPoint() const {
    return 0; // no in-text verification with BMove
}

void BMove::phi(length_t& pos) const {
    // Find the rank of the predecessor of pos in a circular manner.
    length_t predRank = predFirst.predecessorRankCircular(pos);
    // Select the predecessor position using its rank.
    length_t pred = predFirst.select(predRank);

    // Calculate the distance from the predecessor to the current position.
    length_t delta = pred < pos ? pos - pred : pos + 1;

    // Ensure that phi(SA[0]) is not called; the predecessor rank must be valid.
    assert(firstToRun[predRank] > 0);

    // Get the previous sample from the samplesLast array.
    length_t prev_sample = samplesLast[firstToRun[predRank] - 1];

    // Calculate and return the new position, modulo the text length.
    pos = (prev_sample + delta) % textLength;
}

void BMove::phiInverse(length_t& pos) const {
    // Find the rank of the predecessor of pos in a circular manner.
    length_t predRank = predLast.predecessorRankCircular(pos);
    // Select the predecessor position using its rank.
    length_t pred = predLast.select(predRank);

    // Calculate the distance from the predecessor to the current position.
    length_t delta = pred < pos ? pos - pred : pos + 1;

    // Ensure that phiInverse(SA[n-1]) is not called; the predecessor rank must
    // be valid.
    assert(lastToRun[predRank] < samplesFirst.size() - 1);

    // Get the next sample from the samplesFirst array.
    length_t prev_sample = samplesFirst[lastToRun[predRank] + 1];

    // Calculate and return the new position, modulo the text length.
    pos = (prev_sample + delta) % textLength;
}

length_t BMove::computeToehold(const MoveRange& range, const length_t c) const {

    length_t endRun = range.getEndRun();

    uint8_t endRunHead = move.getRunHead(endRun);

    if (endRunHead == c) {
        return samplesFirst[endRun];
    }

    length_t previousPos;
    length_t previousRun;

    move.walkToPreviousRun(range, previousPos, previousRun, c);

    return samplesLast[previousRun];
}

length_t BMove::computeToeholdRev(const MoveRange& range,
                                  const length_t c) const {

    length_t endRun = range.getEndRun();

    length_t endRunHead = moveR.getRunHead(endRun);

    if (endRunHead == c) {
        return revSamplesFirst[endRun];
    }

    length_t previousPos;
    length_t previousRun;

    moveR.walkToPreviousRun(range, previousPos, previousRun, c);

    return revSamplesLast[previousRun];
}

// ----------------------------------------------------------------------------
//  APPROXIMATE PATTERN MATCHING
// ----------------------------------------------------------------------------

void BMove::updateRangeSARuns(SARangePair& ranges) const {
    move.computeRunIndices(ranges.getRangeSAMutable());
}

void BMove::findRangeWithExtraCharBackwardAuxiliary(
    length_t positionInAlphabet, SARange& parentBackwardRange,
    SARange& childBackwardRange) const {

    if (!parentBackwardRange.getRunIndicesValid()) {
        move.computeRunIndices(parentBackwardRange);
    }
    move.addChar(parentBackwardRange, childBackwardRange, positionInAlphabet);
}

bool BMove::findRangeWithExtraCharBackward(length_t posInAlpha,
                                           const SARangeBackwards& rangeOfP,
                                           SARangeBackwards& childRange) const {

    // first make the backward range by searching cP using B
    SARange range1, trivialRange = rangeOfP;
    findRangeWithExtraCharBackwardAuxiliary(posInAlpha, trivialRange, range1);

    if (range1.empty()) {
        childRange = SARangeBackwards(range1, 0, false, 0);
        return false;
    }

    // Check if range1 has the same width as the range of the parent
    if (trivialRange.width() == range1.width()) {
        childRange = SARangeBackwards(
            range1, rangeOfP.getToehold() - !rangeOfP.getToeholdRepresentsEnd(),
            rangeOfP.getToeholdRepresentsEnd(),
            rangeOfP.getOriginalDepth() + 1);
        return true;
    }

    // aP occurs for some a != c, the toehold must be updated
    length_t newToehold = computeToehold(trivialRange, posInAlpha);

    childRange = SARangeBackwards(range1, newToehold, false,
                                  rangeOfP.getOriginalDepth() + 1);
    return true;
}

bool BMove::findRangesWithExtraCharBackward(length_t positionInAlphabet,
                                            const SARangePair& rangesOfP,
                                            SARangePair& childRanges) const {

    // first make the backward range by searching cP using B
    SARange range1, trivialRange = rangesOfP.getRangeSA();
    findRangeWithExtraCharBackwardAuxiliary(positionInAlphabet, trivialRange,
                                            range1);

    // if the range is empty, return false
    if (range1.empty()) {
        childRanges = SARangePair(range1, range1, 0, false, 0);
        return false;
    }

    // Check if range1 has the same width as the range of the parent
    const SARange& otherRange = rangesOfP.getRangeSARev();
    if (trivialRange.width() == range1.width()) {
        // otherRange remains the same
        childRanges = SARangePair(range1, otherRange,
                                  rangesOfP.getToehold() -
                                      !rangesOfP.getToeholdRepresentsEnd(),
                                  rangesOfP.getToeholdRepresentsEnd(),
                                  rangesOfP.getOriginalDepth() + 1);
        return true;
    }

    // then make the less trivial range by counting the sizes of the ranges
    // of (dP) using B

    // first get the start of the range we are looking for of the parent
    length_t s = rangesOfP.getRangeSARev().getBegin();

    // get the start of the child within this range
    // find the number of occurrences of chars smaller than c in the parent
    // range
    length_t x = move.getCumulativeCounts(trivialRange, positionInAlphabet);

    // make the new range with width equal to that of the trivial range
    SARange range2 = SARange(s + x, s + x + range1.width(),
                             otherRange.getBeginRun(), otherRange.getEndRun());
    // set the run indices to be invalid
    range2.setRunIndicesValid(false);

    // aP occurs for some a != c, the toehold must be updated
    length_t newToehold = computeToehold(trivialRange, positionInAlphabet);

    // set the final SARangePair
    childRanges = SARangePair(range1, range2, newToehold, false,
                              rangesOfP.getOriginalDepth() + 1);
    return true;
}

bool BMove::findRangesWithExtraCharForward(length_t positionInAlphabet,
                                           const SARangePair& rangesOfP,
                                           SARangePair& childRanges) const {

    // first make the trivial range by searching (Pc)' using B'  if
    // searching forward we need to use B' so we need the reverse range
    SARange trivialRange = rangesOfP.getRangeSARev();

    // If the run indices of trivialRange are not valid (e.g., after a
    // directions switch), compute them
    if (!trivialRange.getRunIndicesValid()) {
        moveR.computeRunIndices(trivialRange);
    }

    // get the range of the child by adding ona character using move
    SARange range1;
    moveR.addChar(trivialRange, range1, positionInAlphabet);

    // if the range is empty, return false
    if (range1.empty()) {
        childRanges = SARangePair(range1, range1, 0, false, 0);
        return false;
    }

    // Check if range1 has the same width as the range of the parent
    const SARange& otherRange = rangesOfP.getRangeSA();
    if (trivialRange.width() == range1.width()) {
        // otherRange remains the same
        childRanges = SARangePair(otherRange, range1,
                                  rangesOfP.getToehold() +
                                      rangesOfP.getToeholdRepresentsEnd(),
                                  rangesOfP.getToeholdRepresentsEnd(),
                                  rangesOfP.getOriginalDepth() + 1);
        return true;
    }

    // then make the less trivial range by counting the size of the range of
    // (Pd)' using B' (forward)

    // first get the start of the range we are looking for of the parent
    length_t s = rangesOfP.getRangeSA().getBegin();

    // get the start of the child within this range
    // find the number of occurrences of chars smaller than c in the parent
    // range

    length_t x = moveR.getCumulativeCounts(trivialRange, positionInAlphabet);

    // make the new range with width equal to that of the trivial range
    SARange range2 = SARange(s + x, s + x + range1.width(),
                             otherRange.getBeginRun(), otherRange.getEndRun());
    // set the run indices to be invalid
    range2.setRunIndicesValid(false);

    // aP occurs for some a != c, the toehold must be updated
    length_t newToehold =
        textLength - 1 - computeToeholdRev(trivialRange, positionInAlphabet);

    // set the final SARangePair
    childRanges = SARangePair(range2, range1, newToehold, true,
                              rangesOfP.getOriginalDepth() + 1);
    return true;
}

bool BMove::findRangesWithExtraCharBackwardUniDirectional(
    length_t positionInAlphabet, const SARangePair& rangesOfP,
    SARangePair& rangesOfChild) const {

    // first make the trivial range by  searching cP using B
    SARange range1, trivialRange = rangesOfP.getRangeSA();
    findRangeWithExtraCharBackwardAuxiliary(positionInAlphabet, trivialRange,
                                            range1);

    // if the range is empty, return false
    if (range1.empty()) {
        rangesOfChild = SARangePair(range1, range1, 0, false, 0);
        return false;
    }

    // Check if range1 has the same width as the range of the parent
    if (trivialRange.width() == range1.width()) {
        // otherRange remains the same
        rangesOfChild = SARangePair(range1, SARange(),
                                    rangesOfP.getToehold() -
                                        !rangesOfP.getToeholdRepresentsEnd(),
                                    rangesOfP.getToeholdRepresentsEnd(),
                                    rangesOfP.getOriginalDepth() + 1);
        return true;
    }

    // aP occurs for some a != c, the toehold must be updated
    length_t newToehold = computeToehold(trivialRange, positionInAlphabet);

    // set the final SARangePair
    rangesOfChild = SARangePair(range1, SARange(), newToehold, false,
                                rangesOfP.getOriginalDepth() + 1);
    return true;
}

// ----------------------------------------------------------------------------
// ROUTINES FOR ACCESSING THE DATA STRUCTURE
// ----------------------------------------------------------------------------

SARangePair BMove::getRangeOfSingleChar(char c) const {
    int i = sigma.c2i(c);
    if (i < 0) {
        return SARangePair();
    }

    SARangePair pair = getCompleteRange();
    if ((unsigned int)i < sigma.size() - 1) {
        pair.getRangeSAMutable().setBegin(counts[i]);
        pair.getRangeSAMutable().setEnd(counts[i + 1]);
        pair.getRangeSARevMutable().setBegin(counts[i]);
        pair.getRangeSARevMutable().setEnd(counts[i + 1]);
        pair.getRangeSAMutable().setRunIndicesValid(false);
        pair.getRangeSARevMutable().setRunIndicesValid(false);
        return pair;
    }
    pair.getRangeSAMutable().setBegin(counts[i]);
    pair.getRangeSAMutable().setEnd(textLength);
    pair.getRangeSARevMutable().setBegin(counts[i]);
    pair.getRangeSARevMutable().setEnd(textLength);
    pair.getRangeSAMutable().setRunIndicesValid(false);
    pair.getRangeSARevMutable().setRunIndicesValid(false);
    return pair;
}

// ----------------------------------------------------------------------------
// POST-PROCESSING ROUTINES FOR APPROXIMATE MATCHING
// ----------------------------------------------------------------------------

void BMove::collectTextPositions(length_t firstPos, length_t originalDepth,
                                 std::vector<length_t>& positions) const {

    length_t currentPos = firstPos;

    // Collect positions for the first while loop
    assert(currentPos < textLength);
    positions.push_back(currentPos);

    while (plcp[currentPos] >= originalDepth) {
        phi(currentPos);
        assert(currentPos < textLength);
        positions.push_back(currentPos);
    }

    // Reverse the collected positions from the first while loop
    std::reverse(positions.begin(), positions.end());

    currentPos = firstPos;

    // Collect positions for the second while loop
    while (currentPos != getInitialToehold()) {
        phiInverse(currentPos);
        if (plcp[currentPos] < originalDepth)
            break;
        assert(currentPos < textLength);
        positions.push_back(currentPos);
    }
}

void BMove::getTextPositionsFromSARange(
    const SARangePair& ranges, std::vector<length_t>& positions) const {

    assert(ranges.width() > 0);
    positions.reserve(ranges.width());

    assert(!ranges.getToeholdRepresentsEnd() ||
           ranges.getToehold() >= ranges.getOriginalDepth() - 1);

    length_t firstPos =
        ranges.getToehold() - (ranges.getToeholdRepresentsEnd()
                                   ? (ranges.getOriginalDepth() - 1)
                                   : 0);

    collectTextPositions(firstPos, ranges.getOriginalDepth(), positions);

    assert(ranges.width() == positions.size());
}

std::vector<length_t> BMove::getBeginPositions(const SARangeBackwards& rangeSA,
                                               length_t startDiff,
                                               length_t shift) const {

    assert(rangeSA.width() > 0);
    std::vector<length_t> positions;
    positions.reserve(rangeSA.width());

    assert(!rangeSA.getToeholdRepresentsEnd() ||
           rangeSA.getToehold() >= rangeSA.getOriginalDepth() - 1);

    length_t firstPos =
        rangeSA.getToehold() - (rangeSA.getToeholdRepresentsEnd()
                                    ? (rangeSA.getOriginalDepth() - 1)
                                    : 0);

    collectTextPositions(firstPos, rangeSA.getOriginalDepth(), positions);

    assert(rangeSA.width() == positions.size());
    return positions;
}

// ----------------------------------------------------------------------------
// NOT TO BE USED FUNCTIONS
// ----------------------------------------------------------------------------

void BMove::inTextVerification(const vector<length_t>& startPos,
                               const length_t& maxED, const length_t& minED,
                               Occurrences& occ, Counters& counters,
                               const Substring& pattern,
                               bool fixedStartPos) const {
    throw runtime_error("Function " + string(__func__) +
                        " should not be called for this index.");
}

const std::string& BMove::getText() const {
    throw runtime_error("Function " + string(__func__) +
                        " should not be called for this index.");
}

bool BMove::getNoCIGAR() const {
    throw runtime_error("Function " + string(__func__) +
                        " should not be called for this index.");
}

void BMove::resetInTextMatrices() {
    throw runtime_error("Function " + string(__func__) +
                        " should not be called for this index.");
}

void BMove::inTextVerificationOneString(const length_t startPos,
                                        const length_t endPos,
                                        const length_t& maxED,
                                        const length_t& minED, Occurrences& occ,
                                        Counters& counters,
                                        const std::string& pattern) const {
    throw runtime_error("Function " + string(__func__) +
                        " should not be called for this index.");
}

void BMove::inTextVerificationHamming(const FMPosExt& node, const Search& s,
                                      const std::vector<Substring>& parts,
                                      const length_t idx, Occurrences& occ,
                                      Counters& counters) const {
    throw runtime_error("Function " + string(__func__) +
                        " should not be called for this index.");
}

void BMove::verifyExactPartialMatchInText(const FMOcc& startMatch,
                                          const length_t& beginInPattern,
                                          const length_t& maxED,
                                          Occurrences& occ, Counters& counters,
                                          length_t minED,
                                          const Substring& pattern) {
    throw runtime_error("Function " + string(__func__) +
                        " should not be called for this index.");
}

void BMove::verifyExactPartialMatchInTextHamming(
    FMOcc& startMatch, length_t beginInPattern, length_t maxD,
    const std::vector<Substring>& parts, Occurrences& occ, Counters& counters,
    length_t minD) const {
    throw runtime_error("Function " + string(__func__) +
                        " should not be called for this index.");
}

length_t BMove::findSA(length_t index) const {
    throw runtime_error("Function " + string(__func__) +
                        " should not be called for this index.");
}

void BMove::generateCIGARS(std::vector<TextOcc>& occs, Counters& counters,
                           const ReadBundle& bundle) {
    throw runtime_error("Function " + string(__func__) +
                        " should not be called for this index.");
}

void BMove::generateCIGAR(TextOcc& t, Counters& counters,
                          const Substring& read) const {
    throw runtime_error("Function " + string(__func__) +
                        " should not be called for this index.");
}
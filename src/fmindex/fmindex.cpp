/******************************************************************************
 *  Columba: Approximate Pattern Matching using Search Schemes                *
 *  Copyright (C) 2020-2024 - Luca Renders <luca.renders@ugent.be> and        *
 *                            Lore Depuydt <lore.depuydt@ugent.be> and        *
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
#include "../alphabet.h"
#include "../bitparallelmatrix.h"
#include "../definitions.h"
#include "../indexhelpers.h"
#include "../indexinterface.h"
#include "../logger.h"
#include "../search.h"
#include "../substring.h"

#include <assert.h>
#include <memory>
#include <stdexcept>
#include <utility>

using namespace std;

// ============================================================================
// CLASS FMIndex
// ============================================================================

thread_local vector<uint32_t> FMIndex::zerosBuffer(0, 2 * MAX_K + 1);

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

void FMIndex::getTextPositionsFromSARange(
    const SARangePair& ranges, std::vector<length_t>& positions) const {
    const Range& saRange = ranges.getRangeSA();
    positions.reserve(saRange.width());
    for (length_t i = saRange.getBegin(); i < saRange.getEnd(); i++) {
        positions.push_back(findSA(i));
    }
}

// ----------------------------------------------------------------------------
// ROUTINES FOR INITIALIZATION
// ----------------------------------------------------------------------------

void FMIndex::fromFiles(const string& baseFile, bool verbose) {

    string textFile = baseFile + ".txt.bin";
    if (verbose) {

        logger.logInfo("Reading " + textFile + "...");
    }
    std::ifstream inFile(textFile, std::ios::binary);
    if (!inFile) {
        throw runtime_error("Error opening file for reading: " + textFile);
    }
    // First, read the length of the text
    inFile.read(reinterpret_cast<char*>(&textLength), sizeof(length_t));
    logger.logDeveloper("Text length: " + to_string(textLength) +
                        " characters");
    // Then, read the text data itself
    text.resize(textLength);
    inFile.read(&text[0], textLength);

    inFile.close();

    // read FMIndex specific files
    stringstream ss;
    if (verbose) {

        // read the BWT
        ss << "Reading " << baseFile << ".bwt"
           << "...";
        logger.logInfo(ss);
    }

    if (!bwt.read(baseFile + ".bwt")) {
        throw runtime_error("Cannot open file: " + baseFile + ".bwt");
    }
    if (verbose) {

        // read the baseFile occurrence table
        ss << "Reading " << baseFile << ".brt"
           << "...";
        logger.logInfo(ss);
    }

    if (!forwardRepresentation.read(baseFile + ".brt"))
        throw runtime_error("Cannot open file: " + baseFile + ".brt");
    if (verbose) {

        ss << "Reading " << baseFile << ".rev.brt"
           << "...";
        logger.logInfo(ss);
    }

    // read the reverse baseFile occurrence table
    if (!reverseRepresentation.read(baseFile + ".rev.brt"))
        throw runtime_error("Cannot open file: " + baseFile + ".rev.brt");

    readSequenceNamesAndPositions(baseFile, verbose);
}

// ----------------------------------------------------------------------------
// ROUTINES FOR APPROXIMATE PATTERN MATCHING
// ----------------------------------------------------------------------------

bool FMIndex::findRangesWithExtraCharBackward(
    length_t positionInAlphabet, const SARangePair& rangesOfP,
    SARangePair& rangesOfChild) const {

    // first make the trivial range by  searching cP using B
    const Range& trivialRange = rangesOfP.getRangeSA();

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
    const Range& rangeForTrivial = rangesOfP.getRangeSARev();

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

    length_t x =
        getNumberOfCumOccRev(positionInAlphabet, rangeForTrivial.getEnd()) -
        getNumberOfCumOccRev(positionInAlphabet, rangeForTrivial.getBegin());

    // make the new range
    Range range2 = Range(s + x, s + x + range1.width());

    childRanges = SARangePair(range2, range1);
    return !childRanges.empty();
}

bool FMIndex::findRangeWithExtraCharBackward(length_t posInAlpha,
                                             const Range& rangeOfP,
                                             Range& childRange) const {

    // find the new range by using the LF property
    length_t occBefore = getNumberOfOcc(posInAlpha, rangeOfP.getBegin());
    length_t occAfter = getNumberOfOcc(posInAlpha, rangeOfP.getEnd());

    length_t startInAlphabet = counts[posInAlpha];
    childRange = Range(occBefore + startInAlphabet, occAfter + startInAlphabet);
    return !childRange.empty();
}

bool FMIndex::findRangesWithExtraCharBackwardUniDirectional(
    length_t positionInAlphabet, const SARangePair& rangesOfP,
    SARangePair& rangesOfChild) const {

    // first make the trivial range by  searching cP using B
    const Range& bRange = rangesOfP.getRangeSA();

    // find the new range by using the LF property
    length_t occBefore = getNumberOfOcc(positionInAlphabet, bRange.getBegin());
    length_t occAfter = getNumberOfOcc(positionInAlphabet, bRange.getEnd());

    length_t startInAlphabet = counts[positionInAlphabet];
    Range range1 =
        Range(occBefore + startInAlphabet, occAfter + startInAlphabet);

    rangesOfChild = SARangePair(range1, Range());
    return !rangesOfChild.empty();
}

void FMIndex::verifyExactPartialMatchInText(const FMOcc& startMatch,
                                            const length_t& beginInPattern,
                                            const length_t& maxED,
                                            Occurrences& occ,
                                            Counters& counters, length_t minED,
                                            const Substring& pattern) {
    assert(pattern.getDirection() == FORWARD);
    // Immediately switch to in-text verification
    counters.inc(Counters::IMMEDIATE_SWITCH);

    // A) find out highest possible difference
    length_t startDiff = (beginInPattern == 0) ? 0 : beginInPattern + maxED;

    // B) get the possible begin positions
    const auto& partialStarts =
        getBeginPositions(startMatch.getRanges().getRangeSA(), startDiff);

    // C) verify these positions in the text
    inTextVerification(partialStarts, maxED, minED, occ, counters, pattern,
                       beginInPattern == 0);
}

void FMIndex::inTextVerification(const vector<length_t>& startPos,
                                 const length_t& maxED, const length_t& minED,
                                 Occurrences& occ, Counters& counters,
                                 const Substring& pattern,
                                 bool fixedStartPos) const {

    assert(pattern.getDirection() == FORWARD);

    // find the number of required zeros
    length_t nZeros = (fixedStartPos) ? 1 : 2 * maxED + 1;
    bool matrix64 = use64Matrix(nZeros, maxED);

    IBitParallelED* matrix;
    if (matrix64)
        matrix = fullReadMatrix;
    else
        matrix = fullReadMatrix128;

    // initialize the matrix
    initializeMatrix(matrix, pattern, nZeros, maxED);
    length_t nRows = matrix->getNumberOfRows();

    // get the substrings and their start positions
    vector<Substring> refs;
    refs.reserve(startPos.size());
    for (const auto& start : startPos) {

        // A) find the highest possible end
        length_t maxEnd = textLength - 1; // -1 to remove the $ character
        length_t hEnd = min(maxEnd, start + nRows - 1);
        // B) Get the reference subsequence
        refs.emplace_back(getSubstring(start, hEnd));
    }

    if (matrix64) {
        InTextVerificationTask<uint64_t> task(
            refs, fullReadMatrix, maxED, minED, strand, pairStatus, noCIGAR);
        task.doTask(counters, occ);
    } else {
        InTextVerificationTask<UInt128> task(
            refs, fullReadMatrix128, maxED, minED, strand, pairStatus, noCIGAR);
        task.doTask(counters, occ);
    }
}

void FMIndex::inTextVerificationOneString(const length_t startPos,
                                          const length_t endPos,
                                          const length_t& maxED,
                                          const length_t& minED,
                                          Occurrences& occ, Counters& counters,
                                          const string& pattern) const {
    // with one string we need a single zero in our start column
    bool matrix64 = use64Matrix(1, maxED);

    IBitParallelED* matrix;
    if (matrix64)
        matrix = fullReadMatrix;
    else
        matrix = fullReadMatrix128;

    // initialize the matrix
    initializeMatrix(matrix, pattern, 1, maxED);

    // get the reference substring
    std::vector<Substring> refs = {getSubstring(startPos, endPos)};

    if (matrix64) {
        InTextVerificationTask<uint64_t> task(
            refs, fullReadMatrix, maxED, minED, strand, pairStatus, noCIGAR);
        task.doTask(counters, occ);
    } else {
        InTextVerificationTask<UInt128> task(
            refs, fullReadMatrix128, maxED, minED, strand, pairStatus, noCIGAR);
        task.doTask(counters, occ);
    }
}

void FMIndex::verifyExactPartialMatchInTextHamming(
    const FMOcc& startMatch, length_t beginInPattern, length_t maxD,
    const std::vector<Substring>& parts, Occurrences& occ, Counters& counters,
    length_t minD) const {
    const auto& r = startMatch.getRanges().getRangeSA();

    // get length of pattern (= sum length of parts)
    length_t pSize = parts.back().end();

    Substring pattern(parts[0], 0, pSize, FORWARD);
    inTextVerificationHamming(r, pattern, maxD, minD, beginInPattern, occ,
                              counters);
}

void FMIndex::inTextVerificationHamming(
    const SARange& r, const Substring& pattern, const length_t maxEDFull,
    const length_t minEDFull, const length_t lengthBefore, Occurrences& occ,
    Counters& counters) const {
    std::string CIGAR = "*";
    const auto& pSize = pattern.size();

    if (!noCIGAR) {
        // CIGAR for hamming distance
        CIGAR = fmt::format("{}M", pSize);
    }

    for (length_t i = r.getBegin(); i < r.getEnd(); i++) {

        // look up this position in the suffix array
        length_t Tb = findSA(i);
        counters.inc(Counters::IN_TEXT_STARTED);
        // subtract the length before
        Tb = (Tb > lengthBefore) ? Tb - lengthBefore : 0;

        // Calculate the end in the text  + guard that is does not go over the
        // text length
        length_t Te = Tb + pSize;
        if (Te > textLength) {
            // No match with hamming distance possible
            continue;
        }

        // Create the reference and pattern sequence
        Substring ref = getSubstring(Tb, Te);

        assert(ref.size() == pattern.size());

        length_t score = 0;

        for (length_t j = 0; j < ref.size(); j++) {
            // update the score
            score =
                score + (ref.forwardAccessor(j) != pattern.forwardAccessor(j));
            if (score > maxEDFull) {
                // in text verification failed
                break;
            }
        }
        if (score <= maxEDFull && score >= minEDFull) {
            // in text verification succeeded
            occ.addTextOcc(Range(Tb, Te), score, CIGAR, strand, pairStatus);
        }
    }
}

void FMIndex::inTextVerificationHamming(const FMPosExt& node, const Search& s,
                                        const vector<Substring>& parts,
                                        const length_t idx, Occurrences& occ,
                                        Counters& counters) const {
    // A) find length before and the partial occurrence
    length_t lengthBefore =
        ((idx == 0) ? 0 : parts[s.getLowestPartProcessedBefore(idx)].begin()) -
        (dir == BACKWARD) * (node.getDepth());
    length_t pSize = parts.back().end();

    length_t maxEDFull = s.getMaxED();
    length_t minEDFull = s.getMinED();

    const SARange& r =
        node.getRanges().getRangeSA(); // SA range of current node

    Substring pattern(parts[0], 0, pSize, FORWARD);
    inTextVerificationHamming(r, pattern, maxEDFull, minEDFull, lengthBefore,
                              occ, counters);
}

// ----------------------------------------------------------------------------
// ROUTINES FOR ACCESSING THE DATA STRUCTURE
// ----------------------------------------------------------------------------

SARangePair FMIndex::getRangeOfSingleChar(char c) const {
    int i = sigma.c2i(c);
    if (i < 0) {
        return SARangePair();
    }
    if ((unsigned int)i < sigma.size() - 1) {
        return SARangePair(SARange(counts[i], counts[i + 1]),
                           SARange(counts[i], counts[i + 1]));
    }
    return SARangePair(SARange(counts[i], bwt.size()),
                       SARange(counts[i], bwt.size()));
}

// ----------------------------------------------------------------------------
// NOT TO BE USED FUNCTIONS
// ----------------------------------------------------------------------------

void FMIndex::updateRangeSARuns(SARangePair& ranges) const {
    throw runtime_error("Function " + string(__func__) +
                        " should not be called for this index.");
}
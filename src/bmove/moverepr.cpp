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
 * You should have received a copy of the GNU Affero General Public License   *
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.     *
 ******************************************************************************/

#include "moverepr.h"
#include "../definitions.h"
#include "../indexhelpers.h"
#include "moverow.h"

#include <cassert>
#include <istream>

using namespace std;

void MoveRepr::getRunIndex(const length_t position, length_t& runIndex,
                           pair<length_t, length_t>& possibleRange) const {
    // Iteratively narrow down the possible range using binary search.
    while (possibleRange.second > possibleRange.first) {
        // Use the middle of the possible range as a test value.
        length_t testIndex =
            (possibleRange.first + possibleRange.second + 1) / 2;

        // Eliminate half of the possible range based on the comparison.
        if (this->rows[testIndex].inputStartPos <= position) {
            possibleRange.first = testIndex;
        } else {
            possibleRange.second = testIndex - 1;
        }
    }

    runIndex = possibleRange.first;

    assert(position >= rows[possibleRange.first].inputStartPos);
    assert(possibleRange.first == nrOfRuns - 1 ||
           position < rows[possibleRange.first + 1].inputStartPos);
}

void MoveRepr::computeRunIndices(MoveRange& range) const {
    length_t begin = range.getBegin();
    length_t end = range.getEnd() - 1;
    length_t beginRun = range.getBeginRun();
    length_t endRun = range.getEndRun();
    pair<length_t, length_t> possibleRange(beginRun, endRun);

    getRunIndex(begin, beginRun, possibleRange);
    range.setBeginRun(beginRun);
    possibleRange.second = endRun;
    getRunIndex(end, endRun, possibleRange);
    range.setEndRun(endRun);
    range.setRunIndicesValid(true);
}

bool MoveRepr::walkToNextRun(const MoveRange& startRange, length_t& nextPos,
                             length_t& nextRun, const length_t c) const {
    nextPos = startRange.getBegin();
    nextRun = startRange.getBeginRun();

    length_t endRun = startRange.getEndRun();

    // Iterate through the runs to find the next run containing character c
    while (rows[nextRun].c != c && nextRun <= endRun) {
        nextRun++;
        nextPos = rows[nextRun].inputStartPos;
    }

    // Return true if a run containing character c was found
    return nextRun <= endRun;
}

void MoveRepr::walkToPreviousRun(const MoveRange& startRange,
                                 length_t& previousPos, length_t& previousRun,
                                 const length_t c) const {
    previousPos = startRange.getEnd() - 1;
    previousRun = startRange.getEndRun();

    // Iterate backwards through the runs to find the previous run containing
    // character c
    while (rows[previousRun].c != c) {
        previousPos = rows[previousRun].inputStartPos - 1;
        previousRun--;
    }
}

void MoveRepr::fastForward(const length_t& positionIndex,
                           length_t& runIndex) const {
    // Fast forward the runIndex until it contains
    // the run that contains the positionIndex.
    while (rows[runIndex].inputStartPos <= positionIndex) {
        runIndex++;
        // Ensure runIndex stays within bounds
        assert(runIndex < nrOfRuns + 1);
    }
    runIndex--;
}

void MoveRepr::findLF(length_t& positionIndex, length_t& runIndex) const {
    const auto& row = rows[runIndex];
    length_t offset = positionIndex - row.inputStartPos;
    positionIndex = row.outputStartPos + offset;
    runIndex = row.outputStartRun;
    // Fast forward to the correct runIndex after LF operation
    fastForward(positionIndex, runIndex);
}

void MoveRepr::findLFWithoutFastForward(length_t& positionIndex,
                                        const length_t& runIndex) const {
    const auto& row = rows[runIndex];
    length_t offset = positionIndex - row.inputStartPos;
    positionIndex = row.outputStartPos + offset;
}

void MoveRepr::addChar(const SARange& parentRange, SARange& childRange,
                       const length_t& c) const {
    length_t nextPos;
    length_t nextRun;
    if (!walkToNextRun(parentRange, nextPos, nextRun, c)) {
        childRange.setEmpty();
        return;
    }

    length_t previousPos;
    length_t previousRun;
    walkToPreviousRun(parentRange, previousPos, previousRun, c);

    findLF(nextPos, nextRun);
    findLF(previousPos, previousRun);

    childRange = SARange(nextPos, previousPos + 1, nextRun, previousRun);
}

length_t MoveRepr::countChar(const SARange& parentRange,
                             const length_t& c) const {
    length_t nextPos;
    length_t nextRun;
    if (!walkToNextRun(parentRange, nextPos, nextRun, c)) {
        return 0;
    }

    length_t previousPos;
    length_t previousRun;
    walkToPreviousRun(parentRange, previousPos, previousRun, c);

    findLFWithoutFastForward(nextPos, nextRun);
    findLFWithoutFastForward(previousPos, previousRun);

    return previousPos + 1 - nextPos;
}

length_t MoveRepr::getCumulativeCounts(const SARange& range,
                                       length_t positionInAlphabet) const {

    assert(positionInAlphabet < ALPHABET);
    assert(positionInAlphabet > 0);

    length_t cumulativeCount = 0;
    // zero character, check inputRange to zeroCharPos
    if (range.getBegin() <= zeroCharPos && range.getEnd() > zeroCharPos) {
        // range contains the zeroCharPos
        cumulativeCount++;
    }

    for (length_t c = 1; c < positionInAlphabet; c++) {
        cumulativeCount += countChar(range, c);
    }

    return cumulativeCount;
}

bool MoveRepr::load(const string& baseFile, bool verbose) {
    string fileName = baseFile + ".move";

    ifstream ifs(fileName, ios::binary);
    if (!ifs) {
        return false;
    }

    // Load the bwtSize, amount of input intervals and the alphabet size.
    ifs.read((char*)&bwtSize, sizeof(bwtSize));
    ifs.read((char*)&nrOfRuns, sizeof(nrOfRuns));
    ifs.read((char*)&zeroCharPos, sizeof(zeroCharPos));

    // Load rows: allocate memory and read.
    rows.reserve(nrOfRuns + 1);
    for (length_t i = 0; i < nrOfRuns + 1; i++) {
        rows[i] = MoveRow(ifs);
    }

    ifs.close();
    return true;
}

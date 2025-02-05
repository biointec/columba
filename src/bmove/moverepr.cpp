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
#include "../logger.h"

#include <cassert>
#include <cstdint>
#include <istream>
#include <sys/types.h>

using namespace std;

// ----------------------------------------------------------------------------
//  BitPackedRepresentation
// ----------------------------------------------------------------------------

length_t BitPackedRepresentation::getRowValue(length_t rowIndex, uint16_t bitOffset,
                                   uint8_t numBits) const {
    length_t mask = (1ULL << numBits) - 1;
    uint64_t byteIndex = ((uint64_t)rowIndex) * totalBytes + (bitOffset / 8);
    uint8_t bitIndex = bitOffset % 8;

    const __uint128_t* block =
        reinterpret_cast<const __uint128_t*>(&buffer[byteIndex]);
    __uint128_t currentValue = *block;

    return (currentValue >> bitIndex) & mask;
}

void BitPackedRepresentation::setRowValue(length_t rowIndex, length_t value,
                               uint16_t bitOffset, uint8_t numBits) {
    length_t mask = (1ULL << numBits) - 1;
    value &= mask; // Ensure value fits in the specified number of bits

    uint64_t byteIndex = ((uint64_t)rowIndex) * totalBytes + (bitOffset / 8);
    uint8_t bitIndex = bitOffset % 8;

    __uint128_t* block = reinterpret_cast<__uint128_t*>(&buffer[byteIndex]);
    __uint128_t currentValue = *block;

    __uint128_t clearMask = ~(__uint128_t(mask) << bitIndex);
    currentValue &= clearMask;
    currentValue |= (__uint128_t(value) << bitIndex);
    *block = currentValue;
}

// ----------------------------------------------------------------------------
//  MoveLFReprBP
// ----------------------------------------------------------------------------

bool MoveLFReprBP::initialize(length_t nrOfRuns, length_t textSize) {
    this->nrOfRuns = nrOfRuns;
    this->textSize = textSize;

    // Calculate the number of bits required for n and r
    bitsForC = static_cast<uint8_t>(std::ceil(std::log2(ALPHABET)));
    bitsForN = static_cast<uint8_t>(std::ceil(std::log2(textSize)));
    bitsForR = static_cast<uint8_t>(std::ceil(std::log2(nrOfRuns)));

    logger.logInfo("\tBit-packing move structure with " +
                   std::to_string(nrOfRuns) + " runs and text size " +
                   std::to_string(textSize) + "...");

    // Calculate total bits for each row and convert to bytes
    totalBits = bitsForC + 2 * bitsForN +
                bitsForR; // One char, two n values and one r value
    totalBytes =
        (totalBits + 7) / 8; // Convert total bits to bytes (rounded up)

    logger.logDeveloper("\tBits for c: " + std::to_string(bitsForC));
    logger.logDeveloper("\tBits for n: " + std::to_string(bitsForN));
    logger.logDeveloper("\tBits for r: " + std::to_string(bitsForR));
    logger.logDeveloper("\tTotal bits per row: " + std::to_string(totalBits));
    logger.logDeveloper("\tTotal bytes per row: " + std::to_string(totalBytes));

    // Manually allocate buffer memory
    buffer = new uint8_t[static_cast<uint64_t>(totalBytes) *
                         static_cast<uint64_t>(nrOfRuns + 1)]();

    return true;
}

bool MoveLFReprBP::load(const std::string& baseFile) {
    string fileName = baseFile + ".LFBP";
    std::ifstream ifs(fileName, std::ios::binary);
    if (!ifs) {
        return false; // File could not be opened
    }

    // Load the textSize and nrOfRuns from the file
    ifs.read(reinterpret_cast<char*>(&textSize), sizeof(textSize));
    ifs.read(reinterpret_cast<char*>(&nrOfRuns), sizeof(nrOfRuns));
    ifs.read((char*)&zeroCharPos, sizeof(zeroCharPos));

    // Now calculate the number of bits required for n and r
    bitsForN = static_cast<uint8_t>(std::ceil(std::log2(textSize)));
    bitsForR = static_cast<uint8_t>(std::ceil(std::log2(nrOfRuns)));
    bitsForC = static_cast<uint8_t>(std::ceil(std::log2(ALPHABET)));

    // Calculate total bits for each row and convert to bytes
    totalBits =
        bitsForC + 2 * bitsForN + bitsForR; // Two n values and one r value
    totalBytes =
        (totalBits + 7) / 8; // Convert total bits to bytes (rounded up)

    // Manually allocate buffer memory
    buffer = new uint8_t[static_cast<uint64_t>(totalBytes) *
                         static_cast<uint64_t>(nrOfRuns + 1)]();

    // Load the packed rows into the buffer
    ifs.read(reinterpret_cast<char*>(buffer),
             static_cast<uint64_t>(totalBytes) *
                 static_cast<uint64_t>(nrOfRuns + 1));

    // Ensure the file was properly read
    if (!ifs) {
        return false; // Read failed
    }

    ifs.close();
    return true;
}

bool MoveLFReprBP::write(const std::string& fileName) const {
    std::ofstream ofs(fileName, std::ios::binary);
    if (!ofs) {
        return false; // File could not be opened for writing
    }

    // Write the textSize and nrOfRuns to the file
    ofs.write(reinterpret_cast<const char*>(&textSize), sizeof(textSize));
    ofs.write(reinterpret_cast<const char*>(&nrOfRuns), sizeof(nrOfRuns));
    ofs.write(reinterpret_cast<const char*>(&zeroCharPos), sizeof(zeroCharPos));

    // Write the buffer (packed rows) to the file
    ofs.write(reinterpret_cast<const char*>(buffer),
              static_cast<uint64_t>(totalBytes) *
                  static_cast<uint64_t>(nrOfRuns + 1));

    // Ensure the file was properly written
    if (!ofs) {
        return false; // Write failed
    }

    ofs.close();
    return true;
}

void MoveLFReprBP::setRowValues(length_t rowIndex, uint8_t runChar,
                                length_t inputStartPos, length_t outputStartPos,
                                length_t outputStartRun) {
    setRowValue(rowIndex, runChar, 0, bitsForC); // Pack runChar at bit offset 0
    setRowValue(rowIndex, inputStartPos, bitsForC,
                bitsForN); // Pack inputStartPos at bit offset 0
    setRowValue(rowIndex, outputStartPos, bitsForC + bitsForN,
                bitsForN); // Pack outputStartPos at bit offset bitsForN
    setRowValue(rowIndex, outputStartRun, bitsForC + 2 * bitsForN,
                bitsForR); // Pack outputStartRun at bit offset 2 * bitsForN
}

uint8_t MoveLFReprBP::getRunHead(length_t runIndex) const {
    assert(runIndex < nrOfRuns + 1); // Ensure the index is within bounds
    return getRowValue(runIndex, 0,
                       bitsForC); // Use getRowValue to get the run head
}

length_t MoveLFReprBP::getInputStartPos(length_t i) const {
    assert(i < nrOfRuns + 1); // Ensure the index is within bounds
    return getRowValue(i, bitsForC,
                       bitsForN); // Use getRowValue to get the inputStartPos
}

length_t MoveLFReprBP::getOutputStartPos(length_t i) const {
    assert(i < nrOfRuns + 1); // Ensure the index is within bounds
    return getRowValue(i, bitsForC + bitsForN,
                       bitsForN); // Use getRowValue to get the outputStartPos
}

length_t MoveLFReprBP::getOutputStartRun(length_t i) const {
    assert(i < nrOfRuns + 1); // Ensure the index is within bounds
    return getRowValue(i, bitsForC + 2 * bitsForN,
                       bitsForR); // Use getRowValue to get the outputStartRun
}

void MoveLFReprBP::setOutputStartRun(length_t rowIndex, length_t value) {
    assert(rowIndex < nrOfRuns); // Ensure the index is within bounds
    setRowValue(rowIndex, value, bitsForC + 2 * bitsForN,
                bitsForR); // Pack value into the correct bit offset
}

void MoveLFReprBP::getRunIndex(const length_t position, length_t& runIndex,
                               pair<length_t, length_t>& possibleRange) const {
    // Iteratively narrow down the possible range using binary search.
    while (possibleRange.second > possibleRange.first) {
        // Use the middle of the possible range as a test value.
        length_t testIndex =
            (possibleRange.first + possibleRange.second + 1) / 2;

        // Eliminate half of the possible range based on the comparison.
        if (getInputStartPos(testIndex) <= position) {
            possibleRange.first = testIndex;
        } else {
            possibleRange.second = testIndex - 1;
        }
    }

    runIndex = possibleRange.first;

    assert(position >= getInputStartPos(possibleRange.first));
    assert(possibleRange.first == nrOfRuns - 1 ||
           position < getInputStartPos(possibleRange.first + 1));
}

void MoveLFReprBP::computeRunIndices(MoveRange& range) const {
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

bool MoveLFReprBP::walkToNextRun(const MoveRange& startRange, length_t& nextPos,
                                 length_t& nextRun, const length_t c) const {
    nextPos = startRange.getBegin();
    nextRun = startRange.getBeginRun();

    length_t endRun = startRange.getEndRun();

    // Iterate through the runs to find the next run containing character c
    while (getRunHead(nextRun) != c && nextRun <= endRun) {
        nextRun++;
        nextPos = getInputStartPos(nextRun);
    }

    // Return true if a run containing character c was found
    return nextRun <= endRun;
}

void MoveLFReprBP::walkToPreviousRun(const MoveRange& startRange,
                                     length_t& previousPos,
                                     length_t& previousRun,
                                     const length_t c) const {
    previousPos = startRange.getEnd() - 1;
    previousRun = startRange.getEndRun();

    // Iterate backwards through the runs to find the previous run containing
    // character c
    while (getRunHead(previousRun) != c) {
        previousPos = getInputStartPos(previousRun) - 1;
        previousRun--;
    }
}

void MoveLFReprBP::fastForward(const length_t& positionIndex,
                               length_t& runIndex) const {
    // Fast forward the runIndex until it contains
    // the run that contains the positionIndex.
    while (getInputStartPos(runIndex) <= positionIndex) {
        runIndex++;
        // Ensure runIndex stays within bounds
        assert(runIndex < nrOfRuns + 1);
    }
    runIndex--;
}

void MoveLFReprBP::findLF(length_t& positionIndex, length_t& runIndex) const {
    length_t offset = positionIndex - getInputStartPos(runIndex);
    positionIndex = getOutputStartPos(runIndex) + offset;
    runIndex = getOutputStartRun(runIndex);
    // Fast forward to the correct runIndex after LF operation
    fastForward(positionIndex, runIndex);
}

void MoveLFReprBP::findLFWithoutFastForward(length_t& positionIndex,
                                            const length_t& runIndex) const {
    length_t offset = positionIndex - getInputStartPos(runIndex);
    positionIndex = getOutputStartPos(runIndex) + offset;
}

void MoveLFReprBP::addChar(const SARange& parentRange, SARange& childRange,
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

length_t MoveLFReprBP::countChar(const SARange& parentRange,
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

length_t MoveLFReprBP::getCumulativeCounts(const SARange& range,
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

// ----------------------------------------------------------------------------
//  MovePhiReprBP
// ----------------------------------------------------------------------------

bool MovePhiReprBP::initialize(length_t nrOfRuns, length_t textSize) {
    this->nrOfRuns = nrOfRuns;
    this->textSize = textSize;

    // Calculate the number of bits required for n and r
    bitsForN = static_cast<uint8_t>(std::ceil(std::log2(textSize)));
    bitsForR = static_cast<uint8_t>(std::ceil(std::log2(nrOfRuns)));

    logger.logInfo("\tInitializing MovePhiReprBP with " +
                   std::to_string(nrOfRuns) + " runs and text size " +
                   std::to_string(textSize) + "...");

    // Calculate total bits for each row and convert to bytes
    totalBits = 2 * bitsForN + bitsForR; // Two n values and one r value
    totalBytes =
        (totalBits + 7) / 8; // Convert total bits to bytes (rounded up)

    logger.logDeveloper("\tBits for n: " + std::to_string(bitsForN));
    logger.logDeveloper("\tBits for r: " + std::to_string(bitsForR));
    logger.logDeveloper("\tTotal bits per row: " + std::to_string(totalBits));
    logger.logDeveloper("\tTotal bytes per row: " + std::to_string(totalBytes));

    // Manually allocate buffer memory
    buffer = new uint8_t[static_cast<uint64_t>(totalBytes) * static_cast<uint64_t>(nrOfRuns + 1)]();

    return true;
}

bool MovePhiReprBP::load(const std::string& fileName) {
    std::ifstream ifs(fileName, std::ios::binary);
    if (!ifs) {
        return false; // File could not be opened
    }

    // Load the textSize and nrOfRuns from the file
    ifs.read(reinterpret_cast<char*>(&textSize), sizeof(textSize));
    ifs.read(reinterpret_cast<char*>(&nrOfRuns), sizeof(nrOfRuns));

    // Now calculate the number of bits required for n and r
    bitsForN = static_cast<uint8_t>(std::ceil(std::log2(textSize)));
    bitsForR = static_cast<uint8_t>(std::ceil(std::log2(nrOfRuns)));

    // Calculate total bits for each row and convert to bytes
    totalBits = 2 * bitsForN + bitsForR; // Two n values and one r value
    totalBytes =
        (totalBits + 7) / 8; // Convert total bits to bytes (rounded up)

    // Manually allocate buffer memory
    buffer = new uint8_t[static_cast<uint64_t>(totalBytes) * static_cast<uint64_t>(nrOfRuns + 1)]();

    // Load the packed rows into the buffer
    ifs.read(reinterpret_cast<char*>(buffer), static_cast<uint64_t>(totalBytes) * static_cast<uint64_t>(nrOfRuns + 1));

    // Ensure the file was properly read
    if (!ifs) {
        return false; // Read failed
    }

    ifs.close();
    return true;
}

bool MovePhiReprBP::write(const std::string& fileName) const {
    std::ofstream ofs(fileName, std::ios::binary);
    if (!ofs) {
        return false; // File could not be opened for writing
    }

    // Write the textSize and nrOfRuns to the file
    ofs.write(reinterpret_cast<const char*>(&textSize), sizeof(textSize));
    ofs.write(reinterpret_cast<const char*>(&nrOfRuns), sizeof(nrOfRuns));

    // Write the buffer (packed rows) to the file
    ofs.write(reinterpret_cast<const char*>(buffer),
              static_cast<uint64_t>(totalBytes) * static_cast<uint64_t>(nrOfRuns + 1));

    // Ensure the file was properly written
    if (!ofs) {
        return false; // Write failed
    }

    ofs.close();
    return true;
}

void MovePhiReprBP::setRowValues(length_t rowIndex, length_t inputStartPos,
                                 length_t outputStartPos,
                                 length_t outputStartRun) {
    setRowValue(rowIndex, inputStartPos, 0,
                bitsForN); // Pack inputStartPos at bit offset 0
    setRowValue(rowIndex, outputStartPos, bitsForN,
                bitsForN); // Pack outputStartPos at bit offset bitsForN
    setRowValue(rowIndex, outputStartRun, 2 * bitsForN,
                bitsForR); // Pack outputStartRun at bit offset 2 * bitsForN
}

length_t MovePhiReprBP::getInputStartPos(length_t i) const {
    assert(i < nrOfRuns + 1); // Ensure the index is within bounds
    return getRowValue(i, 0,
                       bitsForN); // Use getRowValue to get the inputStartPos
}

length_t MovePhiReprBP::getOutputStartPos(length_t i) const {
    assert(i < nrOfRuns + 1); // Ensure the index is within bounds
    return getRowValue(i, bitsForN,
                       bitsForN); // Use getRowValue to get the outputStartPos
}

length_t MovePhiReprBP::getOutputStartRun(length_t i) const {
    assert(i < nrOfRuns + 1); // Ensure the index is within bounds
    return getRowValue(i, 2 * bitsForN,
                       bitsForR); // Use getRowValue to get the outputStartRun
}

void MovePhiReprBP::setOutputStartRun(length_t rowIndex, length_t value) {
    assert(rowIndex < nrOfRuns); // Ensure the index is within bounds
    setRowValue(rowIndex, value, 2 * bitsForN,
                bitsForR); // Pack value into the correct bit offset
}

void MovePhiReprBP::fastForward(const length_t& positionIndex,
                                length_t& runIndex) const {
    // Fast forward the runIndex until it contains
    // the run that contains the positionIndex.
    while (getInputStartPos(runIndex) <= positionIndex) {
        runIndex++;
        // Ensure runIndex stays within bounds
        assert(runIndex < nrOfRuns + 1);
    }
    runIndex--;
}

void MovePhiReprBP::phi(length_t& positionIndex, length_t& runIndex) const {
    assert(runIndex < nrOfRuns);
    assert(positionIndex < textSize);
    length_t offset = positionIndex - getInputStartPos(runIndex);
    positionIndex = getOutputStartPos(runIndex) + offset;
    runIndex = getOutputStartRun(runIndex);
    // Fast forward to the correct runIndex after LF operation
    fastForward(positionIndex, runIndex);
}

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
#ifndef MOVELFREPR_H
#define MOVELFREPR_H

#include "../definitions.h"  // for length_t
#include "../indexhelpers.h" // for SARange
#include <memory>            // for allocator_traits<>::value_type
#include <stdint.h>          // for uint8_t
#include <string>            // for string
#include <utility>           // for pair
#include <vector>            // for vector

// ----------------------------------------------------------------------------
//  BitPackedRepresentation
// ----------------------------------------------------------------------------

class BitPackedRepresentation {
  protected:
    uint8_t* buffer; // Buffer for bit-packed data
    length_t nrOfRuns;
    length_t textSize;
    uint8_t bitsForN, bitsForR; // Common bits for position and run indices
    uint16_t totalBits, totalBytes;

    // Helper function to calculate bits needed for maximum value
    uint8_t calculateBits(length_t maxValue) const {
        return static_cast<uint8_t>(std::ceil(std::log2(maxValue)));
    }

  public:
    BitPackedRepresentation()
        : buffer(nullptr), nrOfRuns(0), textSize(0), bitsForN(0), bitsForR(0),
          totalBits(0), totalBytes(0) {
    }
    virtual ~BitPackedRepresentation() {
        delete[] buffer; // Clean up buffer
    }

    // return the number of runs
    length_t size() const {
        return nrOfRuns;
    }

    // Pure virtual methods to be implemented in each derived class
    virtual bool initialize(length_t nrOfRuns, length_t textSize) = 0;
    virtual bool load(const std::string& fileName) = 0;
    virtual bool write(const std::string& fileName) const = 0;

    // Shared helpers for reading/writing bits in buffer
    length_t getRowValue(length_t rowIndex, uint16_t bitOffset,
                         uint8_t numBits) const;
    void setRowValue(length_t rowIndex, length_t value, uint16_t bitOffset,
                     uint8_t numBits);
};

// ----------------------------------------------------------------------------
//  MoveLFReprBP
// ----------------------------------------------------------------------------

class MoveLFReprBP : public BitPackedRepresentation {
  private:
    uint8_t bitsForC; // Bits for the 'c' field, unique to MoveLFReprBP

    length_t zeroCharPos; // The position index of the zero character

  public:
    // Default constructor
    MoveLFReprBP() : bitsForC(0), zeroCharPos(0) {
    }

    // Initialize function to set up bit packing
    bool initialize(length_t nrOfRuns, length_t bwtSize) override;

    // Load the move representation from a file
    bool load(const std::string& fileName) override;

    // Write the move representation to a file
    bool write(const std::string& fileName) const override;

    // Initialize values for each row
    void setRowValues(length_t rowIndex, uint8_t runChar,
                      length_t inputStartPos, length_t outputStartPos,
                      length_t outputStartRun);

    // Get the run head (character) for a given row
    uint8_t getRunHead(length_t runIndex) const;

    // Get the inputStartPos for a given row
    length_t getInputStartPos(length_t runIndex) const;

    // Get the outputStartPos for a given row
    length_t getOutputStartPos(length_t runIndex) const;

    // Get the outputStartRun for a given row
    length_t getOutputStartRun(length_t runIndex) const;

    // Set the outputStartRun value for a specific row
    void setOutputStartRun(length_t rowIndex, length_t value);

    // Set the zeroCharPos value
    void setZeroCharPos(length_t value) {
        zeroCharPos = value;
    }

    // Perform the LF-mapping step for positionIndex, updating runIndex
    void findLF(length_t& positionIndex, length_t& runIndex) const;

    /**
     * @brief Perform the LF operation on the given positionIndex and runIndex
     * without fast forwarding.
     *
     * @param positionIndex The position index to perform the LF operation on.
     * @param runIndex The run index to perform the LF operation on.
     */
    void findLFWithoutFastForward(length_t& positionIndex,
                                  const length_t& runIndex) const;

    // Fast-forward to the run containing a given position index
    void fastForward(const length_t& positionIndex, length_t& runIndex) const;

    /**
     * Find the run index for a given position.
     * @param position The position for which to find the run index.
     * @param runIndex The run index to be updated.
     * @param possibleRange The possible range for binary search.
     *
     * @return void
     */
    void getRunIndex(const length_t position, length_t& runIndex,
                     std::pair<length_t, length_t>& possibleRange) const;

    /**
     * Compute run indices of the positions of the given range.
     * @param range The range of which we will compute the run indices.
     *
     * @return void
     */
    void computeRunIndices(MoveRange& range) const;

    /**
     * @brief Let startRange be an SA interval, which we want to extend with a
     * character c. This function will update SA index nextPos and run index
     * nextRun to indicate the next position in the SA interval corresponding to
     * an occurrence of c in the BWT (starting from the beginning of the
     * interval).
     *
     * @param startRange The SA interval to extend with character c.
     * @param nextPos The position in the SA interval that will contain the next
     * occurrence of c
     * @param nextRun The run index that will contain the next occurrence of c
     * @param c The character to extend the SA interval with
     * @return true - if the SA interval is not empty, false otherwise
     */
    bool walkToNextRun(const MoveRange& startRange, length_t& nextPos,
                       length_t& nextRun, const length_t c) const;

    /**
     * @brief Let startRange be an SA interval, which we want to extend with a
     * character c. This function will update SA index previousPos and run index
     * previousRun to indicate the previous position in the SA interval
     * corresponding to an occurrence of c in the BWT (starting from the end of
     * the interval).
     *
     * @param startRange The SA interval to extend with character c.
     * @param previousPos The position in the SA interval that will contain the
     * previous occurrence of c
     * @param previousRun The run index that will contain the previous
     * occurrence of c
     * @param c The character to extend the SA interval with
     * @return true - if the SA interval is not empty, false otherwise
     */
    void walkToPreviousRun(const MoveRange& startRange, length_t& previousPos,
                           length_t& previousRun, const length_t c) const;

    /**
     * @brief Extend the match corresponding to parentRange by prepending
     * character c to the match.
     *
     * @param parentRange The match to extend.
     * @param childRange The extended match.
     * @param c The character to prepend to the match.
     */
    void addChar(const SARange& parentRange, SARange& childRange,
                 const length_t& c) const;

    /**
     * @brief Extend the match corresponding to parentRange by appending
     * character c to the match. Only the count is necessary so no fast
     * forwarding.
     *
     * @param parentRange The match to extend.
     * @param c The character to append to the match.
     * @return length_t - The count of the character c in the extended match.
     */
    length_t countChar(const SARange& parentRange, const length_t& c) const;

    /**
     * @brief Extend the match corresponding to range by appending all
     * characters lexicographically smaller than c to the match. Only the count
     * is necessary so no fast forwarding.
     *
     * @param range The match to extend.
     * @param c The character to append to the match.
     */
    length_t getCumulativeCounts(const SARange& range,
                                 length_t positionInAlphabet) const;
};

// ----------------------------------------------------------------------------
//  MovePhiReprBP
// ----------------------------------------------------------------------------

class MovePhiReprBP : public BitPackedRepresentation {
  public:
    // Default constructor
    MovePhiReprBP() {
    }

    // Initialize function to handle memory allocation
    bool initialize(length_t nrOfRuns, length_t textSize) override;

    // Load function to handle initialization
    bool load(const std::string& fileName) override;

    // Write function to handle serialization
    bool write(const std::string& fileName) const override;

    // Initialize values for each row
    void setRowValues(length_t rowIndex, length_t inputStartPos,
                      length_t outputStartPos, length_t outputStartRun);

    // Get the inputStartPos value for a specific row
    length_t getInputStartPos(length_t i) const;

    // Get the outputStartPos value for a specific row
    length_t getOutputStartPos(length_t i) const;

    // Get the outputStartRun value for a specific row
    length_t getOutputStartRun(length_t i) const;

    // Set the outputStartRun value for a specific row
    void setOutputStartRun(length_t rowIndex, length_t value);

    // Fast-forward runIndex to the row containing positionIndex
    void fastForward(const length_t& positionIndex, length_t& runIndex) const;

    // Perform the LF-mapping step for positionIndex, updating the runIndex
    void phi(length_t& positionIndex, length_t& runIndex) const;
};

#endif
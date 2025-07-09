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
#ifndef BMOVE_H
#define BMOVE_H

#include "../definitions.h"
#include "../indexhelpers.h"
#include "../indexinterface.h"
#include "../reads.h"
#include "moverepr.h"
#include "plcp.h" // for PLCP
#include "sparsebitvec.h"

#include <ios>
#include <sdsl/int_vector.hpp>
#include <string>
#include <vector>
class MemoryMappedTextFile;
class Search;
class Substring;

class BMove : public IndexInterface {
  private:
    // PLCP (permuted longest common prefix) array
    PLCP plcp;

    // First and last sa samples of each input interval, needed for keeping the
    // toehold
    sdsl::int_vector<> samplesFirst;
    sdsl::int_vector<> samplesLast;
    sdsl::int_vector<> revSamplesFirst;
    sdsl::int_vector<> revSamplesLast;

    // Predecessor structures.
    SparseBitvec predFirst;
    SparseBitvec predLast;

#ifndef PHI_MOVE
    // Contain predecessor to run mappings.
    sdsl::int_vector<> firstToRun;
    sdsl::int_vector<> lastToRun;
#else
    // Contains bit packed phi and phi inverse move structures
    MovePhiReprBP phiMove;
    MovePhiReprBP phiInvMove;
#endif

    // The Move data structure of the text and the reverse text
    MoveLFReprBP move;
    MoveLFReprBP moveR;

    // ----------------------------------------------------------------------------
    // PREPROCESSING ROUTINES
    // ----------------------------------------------------------------------------

    /**
     * Private helper function that reads in all the necessary files
     * @param baseFile the baseFile of the files that will be read in
     * @param verbose if true the steps will we written to cout
     */
    void fromFiles(const std::string& baseFile, bool verbose);

    /**
     * Read a binary file and stores content in sdsl int_vector
     * @param filename File name
     * @param intVector output int_vector (contents will be overwritten)
     * @returns True if successful, false otherwise
     */
    bool readIntVector(const std::string& filename,
                       sdsl::int_vector<>& intVector) {
        std::ifstream ifs(filename, std::ios::binary);
        if (!ifs) {
            return false;
        }
        intVector.load(ifs);
        ifs.close();
        return true;
    }

    // ----------------------------------------------------------------------------
    // ROUTINES FOR ACCESSING DATA STRUCTURE
    // ----------------------------------------------------------------------------

#ifndef PHI_MOVE
    /**
     * Finds the text position given the position that succeeds it in the SA,
     * in other words, finds SA[i-1] given SA[i]
     * @param pos Position SA[i] in the text, to be updated to SA[i-1]
     */
    void phi(length_t& pos) const;

    /**
     * Finds the text position given the position that precedes it in the SA,
     * in other words, finds SA[i+1] given SA[i]
     * @param pos Position SA[i] in the text, to be updated to SA[i+1]
     */
    void phiInverse(length_t& pos) const;
#endif

    /**
     * Get toehold for the given character to match.
     * @param range Range in the sa array in which to match char and get
     * toehold.
     * @param c Character to match.
     * @return Index in the SA where the toehold starts.
     */
    length_t computeToehold(const SARange& range, const length_t c) const;

    /**
     * Get toehold for the given character to match in the reverse text.
     * @param range Range in the sa array in which to match char and get
     * toehold.
     * @param c Character to match.
     * @return Index in the SA where the toehold starts.
     */
    length_t computeToeholdRev(const SARange& range, const length_t c) const;

    /**
     * Get an initial toehold for the full range in the BWT.
     * @returns A toehold for the full range in the BWT.
     */
    length_t getInitialToehold() const {
        assert(samplesLast[samplesLast.size() - 1]);
        return samplesLast[samplesLast.size() - 1] - 1;
    }

    // ----------------------------------------------------------------------------
    //  APPROXIMATE PATTERN MATCHING
    // ----------------------------------------------------------------------------

    /**
     * @brief Update the runs corresponding to the SA range
     *
     * @param ranges SARangePair for which the SA range runs must be updated
     */
    virtual void updateRangeSARuns(SARangePair& ranges) const override;

    /**
     * Helper function for backwards adding of character (both uni- and
     * bidirectional). Finds the backwards (trivial) range for extending the
     * parent range with a character.
     * @param positionInAlphabet the position in the alphabet of the character.
     * @param parentBackwardRange the backward range of the parent
     * @param childBackwardRange the backward range of the child (output)
     *
     */
    void
    findRangeWithExtraCharBackwardAuxiliary(length_t positionInAlphabet,
                                            SARange& parentBackwardRange,
                                            SARange& childBackwardRange) const;

    /**
     * Finds the ranges of cP using the principle explained in the paper of
     * Lam.
     * @param positionInAlphabet the position in alphabet of the character
     * that is added in the front.
     * @param rangesOfP the ranges of pattern P.
     * @param childRanges the ranges corresponding to string cP, this
     * will be set during execution
     * @returns true if the range is not empty, false otherwise
     */
    virtual bool
    findRangesWithExtraCharBackward(length_t positionInAlphabet,
                                    const SARangePair& rangesOfP,
                                    SARangePair& childRanges) const override;

    /**
     * Finds the ranges of string Pc using the principle explained in the
     * paper of Lam
     * @param positionInAlphabet the position in alphabet of the character c
     * that is added in the back
     * @param rangesOfP the ranges of pattern P
     * @param childRanges the ranges corresponding to string Pc, this will
     * be set during execution
     * @returns true if the range is not empty, false otherwise
     */
    virtual bool
    findRangesWithExtraCharForward(length_t positionInAlphabet,
                                   const SARangePair& rangesOfP,
                                   SARangePair& childRanges) const override;

    /**
     * Finds the range of Pc using unidirectional backward matching
     * @param positionInAlphabet the position in alphabet of the character c
     * that is added in the back
     * @param rangeOfP the range over the suffix array of the text
     * corresponding to pattern P
     * @param childRange the range over the suffix array of text
     * corresponding to pattern Pc this will be set during execution
     * @returns true if the range is not empty, false otherwise
     */
    virtual bool
    findRangeWithExtraCharBackward(length_t positionInAlphabet,
                                   const SARangeBackwards& rangeOfP,
                                   SARangeBackwards& childRange) const override;

    /**
     * Finds the range of cP and stores a 'dummy' range over the suffix
     * array of the reversed text
     * @param positionInAlphabet the position in alphabet of the character
     * that is added in the front.
     * @param rangesOfP the ranges of pattern P, only the backwards range is
     * used
     * @param childRanges the ranges corresponding to string cP, this
     * will be set during execution. The getRangeSARev() function will
     * return a dummy!
     * @returns true if the range is not empty, false otherwise
     */
    virtual bool findRangesWithExtraCharBackwardUniDirectional(
        length_t positionInAlphabet, const SARangePair& rangesOfP,
        SARangePair& rangesOfChild) const override;

    // ----------------------------------------------------------------------------
    // IN TEXT VERIFICATION ROUTINES
    // ----------------------------------------------------------------------------

    /**
     * Verifies the text occurrences in the text and adds them to the
     * occurrences for the edit distance metric. WARNING: This function
     * makes use of the fullReadMatrix pointer. Before calling this function
     * make sure the pointer points to the correct bit-parallel matrix and
     * that the sequence for this matrix has been set.
     * @param sPos The start positions to be verified.
     * @param maxED The maximal edit distance that is allowed for the
     * search.
     * @param minED The minimal edit distance that is allowed for the
     * search.
     * @param occ Data structure with the in-index and in-text occurrences,
     * if an occurrence, either in-index or in-text, is found it will be
     * added to this data structure.
     * @param counters The performance counters.
     * @param pattern The pattern to verify.
     * @param fixedStartPos Whether the start positions are fixed or not.
     */
    virtual void inTextVerification(const std::vector<length_t>& sPos,
                                    const length_t& maxED,
                                    const length_t& minED, Occurrences& occ,
                                    Counters& counters,
                                    const Substring& pattern,
                                    bool fixedStartPos) const override;

    /**
     * Verifies an in-text occurrence and adds it to the occurrences for the
     * edit distance metric. WARNING: This function makes use of the
     * fullReadMatrix pointer. Before calling this function make sure the
     * pointer points to the correct bit-parallel matrix and that the
     * sequence for this matrix has been set.
     * @param startPos the start position of the in-text occurrence to
     * verify
     * @param endPos the end position of the in-text occurrence to verify
     * @param maxED the maximal allowed edit distance
     * @param minED the minimal allowed edit distance
     * @param occ the Occurrences data structure to add a valid in-text
     * occurrence to
     * @param counters the performance counters
     * @param pattern the pattern to verify
     */
    virtual void
    inTextVerificationOneString(const length_t startPos, const length_t endPos,
                                const length_t& maxED, const length_t& minED,
                                Occurrences& occ, Counters& counters,
                                const std::string& pattern) const override;
    /**
     * In text verification for the Hamming distance.
     * @param node The node in the index at which the switch to in-text
     * verification will happen.
     * @param s The current search.
     * @param parts The parts of the pattern for this search.
     * @param idx The current index in the search.
     * @param occ Data structure with all occurrences, both in FM Index and
     * in text. If in-text verification leads to valid text occurrences,
     * these will be added to occ
     */
    virtual void inTextVerificationHamming(const FMPosExt& node,
                                           const Search& s,
                                           const std::vector<Substring>& parts,
                                           const length_t idx, Occurrences& occ,
                                           Counters& counters) const override;

    void inTextVerificationHamming(const Range& r, const Substring& pattern,
                                   const length_t maxEDFull,
                                   const length_t minEDFull,
                                   const length_t lengthBefore,
                                   Occurrences& occ,
                                   Counters& counters) const override;

    // ----------------------------------------------------------------------------
    // LOCATION ROUTINES
    // ----------------------------------------------------------------------------

    /**
     * Find the begin positions of the range that needs to be verified
     * in-text
     * @param rangeSA The range over the SA of the partial in-index match
     * that needs to be verified in text.
     * @param startDiff The highest possible difference between the start of
     * the partial match and that of a valid full match in the text.
     * @param shift The shift of the node. Necessary in case the search in
     * the backwards direction has ended. (default = 0).
     * @returns The lowest possible begin position for each of the partial
     * matches in the range.
     */
    virtual std::vector<length_t>
    getBeginPositions(const SARangeBackwards& rangeSA, length_t startDiff,
                      length_t shift = 0) const override;

  public:
    // ----------------------------------------------------------------------------
    // INITIALIZATION ROUTINES
    // ----------------------------------------------------------------------------

    /**
     * Constructor
     * @param baseFile The base name of the files that contain the info
     * about the index.
     * @param verbose If true, the steps will be written to cout. [default =
     * true]
     * @param noCIGAR If true, the CIGAR string will not be calculated.
     * @param wordSize The size of the mers to be stored in the hashtable.
     * Used for quick look-ups of exact seeds. [default = 10]
     */
    BMove(const std::string& baseFile, bool verbose = true, bool noCIGAR = true,
          length_t wordSize = 10)
        : IndexInterface(baseFile, verbose, noCIGAR, wordSize) {

        // read in files
        fromFiles(baseFile, verbose);

        // populate sparse hash table
        populateTable(verbose);

        // set the index in FORWARD_STRAND mode
        BMove::setIndexInMode(FORWARD_STRAND);
    }

    /**
     * @brief Destructor
     *
     */
    virtual ~BMove() override = default;

    // ----------------------------------------------------------------------------
    // ROUTINES FOR ACCESSING THE DATA STRUCTURE
    // ----------------------------------------------------------------------------

    /**
     * Get the complete range of this index corresponding to an exact match
     * of the empty string.
     * @returns A SARangePair with both ranges the complete range of the
     * index.
     */
    virtual SARangePair getCompleteRange() const override {
        return SARangePair(SARange(0, textLength, 0, move.size() - 1),
                           SARange(0, textLength, 0, moveR.size() - 1),
                           getInitialToehold(), false, 0);
    }

    virtual FMPos getEmptyStringFMPos() const override {
        return FMPos(getCompleteRange(), 0);
    }

    /**
     * Get the original text
     */
    virtual const std::string& getText() const override;

    /**
     * Get the cross-over point form in-index to in-text verification
     */
    virtual length_t getSwitchPoint() const override;

    /**
     * Find the range of an exact match of single character in this index.
     * Warning: this function assumes that the character is in the alphabet.
     * @returns the ranges of a single character in this index
     */
    virtual SARangePair getRangeOfSingleChar(char c) const override;

    // ----------------------------------------------------------------------------
    // ROUTINES FOR APPROXIMATE MATCHING
    // ----------------------------------------------------------------------------

    /**
     * Sets the index in the correct mode. All found occurrences will now be
     * be labeled as found along this strand and as an occurrence of this
     * read in the pair.
     * @param strand On which strand the current read lies
     * @param pairStatus The status of the current read in its pair.
     * [default = FIRST_IN_PAIR]
     */
    virtual void setIndexInMode(Strand reverseComplement,
                                PairStatus firstRead = FIRST_IN_PAIR) override {
        setIndexInModeSubRoutine(reverseComplement, firstRead);

        // TODO extra compile option to disable this for if CIGAR not required
        // point to the correct fullReadMatrix based on the two parameters
        fullReadMatrix = &getFullReadMatrix(strand, pairStatus);
        fullReadMatrix128 = &getFullReadMatrix128(strand, pairStatus);
    }

    /**
     * Verifies an exact partial match in the text for all occurrences of
     * that exact partial match using the edit distance. If a valid match is
     * found it will be added to the occurrences.
     * @param startMatch The match containing the SA ranges corresponding to
     * the exact partial match  to start from and depth of the exact match.
     * @param beginInPattern The position in the pattern where the partial
     * match starts.
     * @param maxED The maximal allowed edit distance.
     * @param occ A data structure with approximate occurrences of the
     * complete pattern. Both in-index and in-text occurrences are stored.
     * If a new approximate occurrence is found, either in-index or in-text
     * it will be added to this data structure.
     * @param counters The performance counters.
     * @param minED the minimal allowed distance for found occurrences
     * @param pattern The pattern to verify.
     */
    virtual void verifyExactPartialMatchInText(
        const FMOcc& startMatch, const length_t& beginInPattern,
        const length_t& maxED, Occurrences& occ, Counters& counters,
        length_t minED, const Substring& pattern) override;

    /**
     * Verifies an exact partial match in the text for all occurrences of
     * that exact partial match using the Hamming distance. If a valid match
     * is found it will be added to the occurrences.
     * @param startMatch The match containing the SA ranges corresponding to
     * the exact partial match  to start from and depth of the exact match.
     * @param beginInPattern The position in the pattern where the partial
     * match starts.
     * @param maxD The maximal allowed hamming distance.
     * @param occ A data structure with approximate occurrences of the
     * complete pattern. Both in-index and in-text occurrences are stored.
     * If a new approximate occurrence is found, either in-index or in-text
     * it will be added to this data structure.
     * @param counters The performance counters.
     * @param minD the minimal allowed Hamming distance for found
     * occurrences
     */
    virtual void verifyExactPartialMatchInTextHamming(
        const FMOcc& startMatch, length_t beginInPattern, length_t maxD,
        const std::vector<Substring>& parts, Occurrences& occ,
        Counters& counters, length_t minD) const override;

    // ----------------------------------------------------------------------------
    // POST-PROCESSING ROUTINES FOR APPROXIMATE MATCHING
    // ----------------------------------------------------------------------------

    /**
     * Finds the entry in the suffix array of this index. This is
     * computed from the sparse suffix array and the bwt.
     * @param index The index for which the entry in the uncompressed suffix
     * array is computed.
     * @returns the entry in the SA of the index
     */
    virtual length_t findSA(length_t index) const override;

    /**
     * @brief Get the text positions based on one known position in the text,
     * and the depth of the match
     *
     * @param firstPos The first known text position
     * @param originalDepth The depth of the match
     * @param positions The vector to which the text positions will be added
     */
    void collectTextPositions(length_t firstPos, length_t originalDepth,
                              std::vector<length_t>& positions) const;

    /**
     * @brief Get the text positions corresponding to a suffix array range
     *
     * @param range The range in the suffix array
     * @param positions The vector to which the text positions will be added
     */
    virtual void getTextPositionsFromSARange(
        const SARangePair& ranges,
        std::vector<length_t>& positions) const override;
};

#endif // BMOVE_H
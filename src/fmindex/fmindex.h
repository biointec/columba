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

#ifndef FMINDEX_H
#define FMINDEX_H

#include "../bitparallelmatrix.h" // for IBitParallelED
#include "../definitions.h"       // for length_t
#include "../indexhelpers.h"      // for SARange
#include "../indexinterface.h"    // for IndexInterface
#include "../reads.h"             // for ReadBundle
#include "bwtrepr.h"              // for BWTRepresentation
#include "encodedtext.h"          // for EncodedText
#include "suffixArray.h"          // for SparseSuffixArray
#include <algorithm>              // for max
#include <assert.h>               // for assert
#include <cmath>                  // for log2
#include <cstdint>                // for uint16_t
#include <fstream>                // for ifstream
#include <string>                 // for string
#include <vector>                 // for vector

class Search;
class Substring;

class FMIndex : public IndexInterface {
  private:
    // info about the text
    std::string text; // The text

    // info about the suffix array
    int logSparseFactorSA = 5; // the log of the sparse factor

    // bidirectional fm index data structures
    EncodedText<ALPHABET> bwt;  // the bwt string of the reference genome
    SparseSuffixArray sparseSA; // the suffix array of the reference genome
    BWTRepresentation<ALPHABET>
        forwardRepresentation; // the baseFile occurrences table
    BWTRepresentation<ALPHABET>
        reverseRepresentation; // the baseFile occurrences table of the BWT of
                               // the reversed text

    // in-text verification
    length_t inTextSwitchPoint = 50;

    thread_local static std::vector<uint32_t> zerosBuffer;

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
     * Read a text file (e.g. input text, BWT, ...)
     * @param filename File name
     * @param buf Buffer (contents will be overwritten)
     * @returns True if successful, false otherwise
     */
    static bool readText(const std::string& filename, std::string& buf) {
        std::ifstream ifs(filename);
        if (!ifs)
            return false;

        ifs.seekg(0, std::ios::end);
        buf.resize(ifs.tellg());
        ifs.seekg(0, std::ios::beg);
        ifs.read((char*)&buf[0], buf.size());

        return true;
    }

    // ----------------------------------------------------------------------------
    // ROUTINES FOR ACCESSING DATA STRUCTURE
    // ----------------------------------------------------------------------------
    /**
     * Finds the LF mapping of the character at index k in the bwt string
     * @param k the index to find the LF mapping off
     * @returns the row that is the LF mapping of k. It is so that the entry
     * in the suffix array of this return value is one less than the entry
     * in the suffix array at index k
     */
    length_t findLF(length_t k) const;

    /**
     * Function that returns the number of occurrences before an index of
     * the symbol at symbolIndex in the alphabet
     * @param symbolIndex the index of the the symbol in the alphabet to
     * count the occurrences of at index index
     * @param index the index whose entry for symbol in the occurrences table
     * is asked
     * @return the number of occurrences of the symbol before index in the
     * bwt
     */
    length_t getNumberOfOcc(length_t symbolIndex, length_t index) const {
        return forwardRepresentation.occ(symbolIndex, index);
    }
    /**
     * Same as getNumberOfOccurrences, but now in the bwt of the
     * reversed text
     * @param symbolIndex the index in the alphabet to get the number of
     * occurrences of
     * @param index the index in the occurrences table, for which the number
     * of occurrences of the symbol at symbolIndex is asked.
     * @return the number of occurrences at index index in the occurrences
     * table of the bwt of the reversed text for the symbol at symbolIndex
     * in the alphabet
     */
    length_t getNumberOfOccRev(length_t symbolIndex, length_t index) const {
        return reverseRepresentation.occ(symbolIndex, index);
    }

    /**
     * Function that returns the number of occurrences before the index of
     * all symbols smaller than the symbol at symbolIndex in the bwt
     * @param symbolIndex the index in the alphabet whose number of baseFile
     * occurrences is queried.
     * @param index the index whose entry for symbol in the prefixOccurrences
     * table is asked
     * @return the number of occurrences of symbols smaller than symbol at
     * symbolIndex before index index in the bwt
     */
    length_t getNumberOfCumOcc(length_t symbolIndex, length_t index) const {
        return forwardRepresentation.cumOcc(symbolIndex, index);
    }

    /**
     * Function that returns the number of occurrences before the index of
     * all symbols smaller than the symbol at symbolIndex in the bwt of the
     * reversed text
     * @param symbolIndex the index in the alphabet whose number of baseFile
     * occurrences is queried.
     * @param index the index whose entry for symbol in the
     * prefixOccurrences table of the reverse text is asked
     * @return the number of occurrences of symbols smaller than symbol at
     * symbolIndex before index index in the bwt
     */
    length_t getNumberOfCumOccRev(length_t symbolIndex, length_t index) const {
        return reverseRepresentation.cumOcc(symbolIndex, index);
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
     * Finds the ranges of string Pc using the principle explained in the paper
     * of Lam
     * @param positionInAlphabet the position in alphabet of the character c
     * that is added in the back
     * @param rangesOfP the ranges of pattern P
     * @param childRanges the ranges corresponding to string Pc, this will be
     * set during execution
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
                                   const Range& rangeOfP,
                                   Range& childRange) const override;

    /**
     * Finds the range of cP and stores a 'dummy' range over the suffix array of
     * the reversed text
     * @param positionInAlphabet the position in alphabet of the character
     * that is added in the front.
     * @param rangesOfP the ranges of pattern P, only the backwards range is
     * used
     * @param childRanges the ranges corresponding to string cP, this
     * will be set during execution. The getRangeSARev() function will return a
     * dummy!
     * @returns true if the range is not empty, false otherwise
     */
    virtual bool findRangesWithExtraCharBackwardUniDirectional(
        length_t positionInAlphabet, const SARangePair& rangesOfP,
        SARangePair& rangesOfChild) const override;

    // ----------------------------------------------------------------------------
    // IN TEXT VERIFICATION ROUTINES
    // ----------------------------------------------------------------------------

    static bool use64Matrix(length_t nZeros, length_t maxED) {
        // find the  required length of the first column of the matrix
        // this column is initialized with zeros, followed by 1.. maxED

        return BitParallelED64::getMaxFirstColRows() >= nZeros + maxED &&
               BitParallelED64::getMatrixMaxED() >= maxED;
    }
    void initializeMatrix(IBitParallelED* matrix, const Substring& pattern,
                          length_t nZeros, length_t maxED) const {
        zerosBuffer.resize(nZeros);

        if (!matrix->sequenceSet()) {
            assert(pattern.getDirection() == FORWARD);
            matrix->setSequence(pattern);
        }
        matrix->initializeMatrix(maxED, zerosBuffer);
    }
    /**
     * Initializes and selects the correct bit-parallel matrix for the in-text
     * verification. The first column of this matrix will be initialized with
     * nZeros zeros followed by 1..maXEd.
     * @param nZeros The number of zeros to initialize the first column with.
     * @param maxED The maximal edit distance that is allowed for the search.
     * @param pattern The pattern to verify.
     * @returns a pointer to the initialized and selected bit-parallel matrix.
     */
    IBitParallelED*
    initializeAndSelectInTextMatrix(length_t nZeros, length_t maxED,
                                    const Substring& pattern) const {

        bool matrix64 = use64Matrix(nZeros, maxED);

        IBitParallelED* matrix;
        if (matrix64)
            matrix = fullReadMatrix;
        else
            matrix = fullReadMatrix128;

        // initialize the matrix
        initializeMatrix(matrix, pattern, nZeros, maxED);

        return matrix;
    }

    /**
     * Verifies the text occurrences in the text and adds them to the
     * occurrences for the edit distance metric. WARNING: This function makes
     * use of the fullReadMatrix pointer. Before calling this function make sure
     * the pointer points to the correct bit-parallel matrix and that the
     * sequence for this matrix has been set.
     * @param sPos The start positions to be verified.
     * @param maxED The maximal edit distance that is allowed for the search.
     * @param minED The minimal edit distance that is allowed for the search.
     * @param occ Data structure with the in-index and in-text occurrences, if
     * an occurrence, either in-index or in-text, is found it will be added to
     * this data structure.
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
     * pointer points to the correct bit-parallel matrix and that the sequence
     * for this matrix has been set.
     * @param startPos the start position of the in-text occurrence to verify
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
     * @param occ Data structure with all occurrences, both in FM Index and in
     * text. If in-text verification leads to valid text occurrences, these will
     * be added to occ
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
     * Find the begin positions of the range that needs to be verified in-text
     * @param rangeSA The range over the SA of the partial in-index match that
     * needs to be verified in text.
     * @param startDiff The highest possible difference between the start of the
     * partial match and that of a valid full match in the text.
     * @param shift The shift of the node. Necessary in case the search in the
     * backwards direction has ended. (default = 0).
     * @returns The lowest possible begin position for each of the partial
     * matches in the range.
     */
    virtual std::vector<length_t>
    getBeginPositions(const Range& rangeSA, length_t startDiff,
                      length_t shift = 0) const override {

        const auto &saBegin = rangeSA.getBegin(), saEnd = rangeSA.getEnd();

        std::vector<length_t> positions(rangeSA.width(), 0);

        for (length_t i = saBegin; i < saEnd; i++) {
            // Calculate the sum of findSA(i) and shift
            length_t sum = findSA(i) + shift;
            // subtract the startDiff while handling underflow
            positions[i - saBegin] = sum >= startDiff ? sum - startDiff : 0;
        }

        return positions;
    }

  public:
    // ----------------------------------------------------------------------------
    // INITIALIZATION ROUTINES
    // ----------------------------------------------------------------------------

    /**
     * Constructor
     * @param baseFile The base name of the files that contain the info about
     * the index.
     * @param inTextSwitch The switch point from in-index to in-text
     * verification.
     * @param noCIGAR If true, the CIGAR string will not be calculated.
     * @param sa_sparse The sparseness factor of suffix array. It is assumed
     * this is a power of two. [default = 1]
     * @param verbose If true, the steps will be written to cout. [default =
     * true]
     * @param wordSize The size of the mers to be stored in the hashtable. Used
     * for quick look-ups of exact seeds. [default = 10]
     */
    FMIndex(const std::string& baseFile, length_t inTextSwitch, bool noCIGAR,
            int sa_sparse = 1, bool verbose = true, length_t wordSize = 10)
        : IndexInterface(baseFile, verbose, noCIGAR, wordSize),
          logSparseFactorSA(log2(sa_sparse)), sparseSA(baseFile, sa_sparse),
          inTextSwitchPoint(inTextSwitch) {

        // read in files
        fromFiles(baseFile, verbose);

        // populate sparse hash table
        populateTable(verbose);

        // set the index in FORWARD_STRAND mode
        FMIndex::setIndexInMode(FORWARD_STRAND);
    }

    /**
     * @brief Destructor
     *
     */
    virtual ~FMIndex() override = default;

    /**
     * Get the complete range of this index corresponding to an exact match of
     * the empty string.
     * @returns A SARangePair with both ranges the complete range of the
     * index.
     */
    virtual SARangePair getCompleteRange() const override {
        return SARangePair(SARange(0, textLength), SARange(0, textLength));
    }

    virtual FMPos getEmptyStringFMPos() const override {
        return FMPos(getCompleteRange(), 0);
    }

    // ----------------------------------------------------------------------------
    // ROUTINES FOR ACCESSING THE DATA STRUCTURE
    // ----------------------------------------------------------------------------

    /**
     * Get the original text
     */
    virtual const std::string& getText() const override {
        return text;
    }

    /**
     * Get a reference subsequence.
     * @param b The start position of the subsequence.
     * @param e The end position of the subsequence.
     * @returns the reference subsequence.
     */
    Substring getSubstring(length_t start, length_t end) const {
        return Substring(text, start, end);
    }

    /**
     * Get a reference subsequence.
     * @param range The range of the subsequence in the concatenated reference
     * text.
     * @returns the reference subsequence.
     */
    Substring getSubstring(const Range& range) const {
        return getSubstring(range.getBegin(), range.getEnd());
    }

    /**
     * Get the cross-over point form in-index to in-text verification
     */
    virtual length_t getSwitchPoint() const override {
        return inTextSwitchPoint;
    }

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
     * read in the pair. The matrix used for in-text verification will now
     * also use the correct matrix (assuming that this matrix has been set
     * correctly with setSequenceInTextMatrix()).
     * @param strand On which strand the current read lies
     * @param pairStatus The status of the current read in its pair. [default =
     * FIRST_IN_PAIR]
     */
    virtual void setIndexInMode(Strand reverseComplement,
                                PairStatus firstRead = FIRST_IN_PAIR) override {
        setIndexInModeSubRoutine(reverseComplement, firstRead);
        // point to the correct fullReadMatrix based on the two parameters
        fullReadMatrix = &getFullReadMatrix(strand, pairStatus);
        fullReadMatrix128 = &getFullReadMatrix128(strand, pairStatus);
    }

    /**
     * Verifies an exact partial match in the text for all occurrences of that
     * exact partial match using the edit distance. If a valid match is found it
     * will be added to the occurrences.
     * @param startMatch The match containing the SA ranges corresponding to
     * the exact partial match  to start from and depth of the exact match.
     * @param beginInPattern The position in the pattern where the partial match
     * starts.
     * @param maxED The maximal allowed edit distance.
     * @param occ A data structure with approximate occurrences of the complete
     * pattern. Both in-index and in-text occurrences are stored. If a new
     * approximate occurrence is found, either in-index or in-text it will be
     * added to this data structure.
     * @param counters The performance counters.
     * @param minED the minimal allowed distance for found occurrences
     * @param pattern The pattern to verify.
     */
    virtual void verifyExactPartialMatchInText(
        const FMOcc& startMatch, const length_t& beginInPattern,
        const length_t& maxED, Occurrences& occ, Counters& counters,
        length_t minED, const Substring& pattern) override;

    /**
     * Verifies an exact partial match in the text for all occurrences of that
     * exact partial match using the Hamming distance. If a valid match is found
     * it will be added to the occurrences.
     * @param startMatch The match containing the SA ranges corresponding to
     * the exact partial match  to start from and depth of the exact match.
     * @param beginInPattern The position in the pattern where the partial match
     * starts.
     * @param maxD The maximal allowed hamming distance.
     * @param occ A data structure with approximate occurrences of the complete
     * pattern. Both in-index and in-text occurrences are stored. If a new
     * approximate occurrence is found, either in-index or in-text it will be
     * added to this data structure.
     * @param counters The performance counters.
     * @param minD the minimal allowed Hamming distance for found occurrences
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
     * @brief Get the text positions corresponding to a suffix array range
     *
     * @param range The range in the suffix array
     * @param positions The vector to which the text positions will be added
     */
    virtual void getTextPositionsFromSARange(
        const SARangePair& ranges,
        std::vector<length_t>& positions) const override;
};

#endif

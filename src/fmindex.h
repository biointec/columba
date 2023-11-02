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

#ifndef FMINDEX_H
#define FMINDEX_H

#include <google/sparse_hash_map>

#include "alphabet.h"
#include "bwtrepr.h"
#include "encodedtext.h"
#include "fmindexhelpers.h"
#include "search.h"
#include "suffixArray.h"
#include "tkmer.h"

#include <fstream>  // used for reading in files
#include <iostream> // used for printing
#include <math.h>   // for taking the log
#include <numeric>  // for summing over vector
#include <string>   // strings
#include <vector>   // vectors

class FMIndex;
typedef bool (FMIndex::*ExtraCharPtr)(length_t, const SARangePair&,
                                      SARangePair&) const;

class FMIndex {
  private:
    // info about the text
    const std::string baseFile; //  The basefile of the reference text
    length_t textLength;        // the length of the text
    std::string text;           // The text

    Alphabet<ALPHABET> sigma; // the alphabet
    length_t sigmaSize;       // Size of the alfabet

    // info about the suffix array
    length_t sparseFactorSA =
        32; // the sparseness factor of the suffix array, defaults to 32
    int logSparseFactorSA = 5; // the log of the sparse factor

    // bidirectional fm index data structures
    EncodedText<ALPHABET> bwt;    // the bwt string of the reference genome
    std::vector<length_t> counts; // the counts array of the reference genome
    SparseSuffixArray sparseSA;   // the suffix array of the reference genome
    BWTRepr<ALPHABET> fwdRepr;    // the baseFile occurrences table
    BWTRepr<ALPHABET> revRepr;    // the baseFile occurrences of the rev BWT

    // in-text verification
    length_t inTextSwitchPoint = 50;

    // direction variables
    thread_local static Direction dir; // the direction of the index
    thread_local static ExtraCharPtr
        extraChar; // pointer to extra char method (for direction)

    // stacks for search schemes
    thread_local static std::vector<std::vector<FMPosExt>>
        stacks; // stacks of nodes for the different partitions
    thread_local static std::vector<BitParallelED>
        matrices; // alignment matrices for the different partitions

    thread_local static BitParallelED
        inTextMatrix; // alignment matrix for the whole read

    // sparse hash info
    const size_t wordSize = 10; // the size of the mers to be stored in a table
    google::sparse_hash_map<Kmer, SARangePair, KmerHash>
        table; // hashtable that contains all wordSize-mers

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
     * Read a binary file and stores content in array
     * @param filename File name
     * @param array Suffix array (contents will be overwritten)
     * @returns True if successful, false otherwise
     */
    static bool readArray(const std::string& filename,
                          std::vector<length_t>& array) {
        std::ifstream ifs(filename, std::ios::binary);
        if (!ifs)
            return false;

        ifs.seekg(0, std::ios::end);
        array.resize(ifs.tellg() / sizeof(length_t));
        ifs.seekg(0, std::ios::beg);
        ifs.read((char*)&array[0], array.size() * sizeof(length_t));

        return true;
    }

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
    /**
     * Populate the hash table
     * @param verbose if steps are written to cout
     */
    void populateTable(bool verbose);

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
     * the symbol at symbolindex in the alphabet
     * @param symbolIndex the index of the the symbol in the alphabet to
     * count the occrences of at index index
     * @param index the index whose entry for symbol in the occurrences table
     * is asked
     * @return the number of occurrences of the symbol before index in the
     * bwt
     */
    length_t getNumberOfOcc(length_t symbolIndex, length_t index) const {
        return fwdRepr.occ(symbolIndex, index);
    }
    /**
     * Same as BiBWT::getNumberOfoccurrences, but now in the bwt of the
     * reversed text
     * @param symbolIndex the index in the alphabet to get the number of
     * occurrences of
     * @param index the index in the occurrences table, for which the number
     * of occurrences of the symbol at symbolindex is asked.
     * @return the number of occurrences at index index in the occurrences
     * table of the bwt of the reversed text for the symbol at symbolIndex
     * in the alphabet
     */
    length_t getNumberOfOccRev(length_t symbolIndex, length_t index) const {
        return revRepr.occ(symbolIndex, index);
    }

    /**
     * Function that returns the number of occurrences before the index of
     * all symbols smaller than the symbol at symbolindex in the bwt
     * @param symbolIndex the index in the alphabet whose number of baseFile
     * occurrences is queried.
     * @param index the index whose entry for symbol in the prefixoccurrences
     * table is asked
     * @return the number of occurrences of symbols smaller than symbol at
     * symbolindex before index index in the bwt
     */
    length_t getNumberOfCumOcc(length_t symbolIndex, length_t index) const {
        return fwdRepr.cumOcc(symbolIndex, index);
    }

    /**
     * Function that returns the number of occurrences before the index of
     * all symbols smaller than the symbol at symbolindex in the bwt of the
     * reversed text
     * @param symbolIndex the index in the alphabet whose number of baseFile
     * occurrences is queried.
     * @param index the index whose entry for symbol in the
     * rprefixoccurrences table of the reverse text is asked
     * @return the number of occurrences of symbols smaller than symbol at
     * symbolindex before index index in the bwt
     */
    length_t getNumberOfCumOccRev(length_t symbolIndex, length_t index) const {
        return revRepr.cumOcc(symbolIndex, index);
    }

    // ----------------------------------------------------------------------------
    //  APPROXIMATE PATTERN MATCHING
    // ----------------------------------------------------------------------------

    /**
     * Matches a search recursively with a depth first approach (each branch
     * of the tree is fully examined until the backtracking condition is
     * met) using the edit distance metric. This function uses all
     * optimizations for eliminating redundancy in the edit distance metric

     * @param search, the search to follow
     * @param startMatch, the approximate match found for all previous
     * partitions of the search
     * @param occ, a datastructure with matches of the complete search, if  such
     * a match is found it will be added to this datastructure
     * @param parts, the parts of the pattern
     * @param counters the performance counters
     * @param idx, the index of the partition to match, defaults to 1 as an
     * exact search for the zeroth partition is assumed
     * @param descPrevDir, the descendants of the previous direction,
     * defaults to empty vector
     * @param initPrevDir, the initialization eds of the previous direction,
     * defaults to empty vector
     * @param descNotPrevDir, the descendants of the other direction,
     * defaults to empty vector
     * @param initNotPrevDir, the initialization eds of the other direction,
     * defaults to empty vector
     */
    void recApproxMatchEditOptimized(
        const Search& search, const FMOcc& startMatch, Occurrences& occ,
        const std::vector<Substring>& parts, Counters& counters,
        const int& idx = 1,
        const std::vector<FMPosExt>& descPrevDir = std::vector<FMPosExt>(),
        const std::vector<uint>& initPrevDir = std::vector<uint>(),
        const std::vector<FMPosExt>& descNotPrevDir = std::vector<FMPosExt>(),
        const std::vector<uint>& initNotPrevDir = std::vector<uint>());

    /**
     * Finds the ranges of cP using the principle explained in the paper of
     * Lahm
     * @param positionInAlphabet the position in alphabet of the character
     * that is added in the front
     * @param rangesOfP the ranges of pattern P
     * @returns the ranges cP
     */
    bool findRangesWithExtraCharBackward(length_t positionInAlphabet,
                                         const SARangePair& rangesOfP,
                                         SARangePair& childRanges) const;

    /**
     * Finds the ranges of Pc using the principle explained in the paper of
     * Lahm
     * @param positionInAlphabet the position in alphabet of the character c
     * that is added in the back
     * @param rangesOfP the ranges of pattern P
     * @returns the ranges of Pc
     */
    bool findRangesWithExtraCharForward(length_t positionInAlphabet,
                                        const SARangePair& rangesOfP,
                                        SARangePair& childRanges) const;

    /**
  * Helper function for the approximate matching. This function fills in
  * the matrix for the current node at the current row and goes deeper
  * for the next part is that is necessary The function true if the search
  * should backtrack.

  * @param clus, the cluster corresponding to the final column of the
  * matrix
  * @param currentnode, the node for which the matrix is filled in
  * @param s, the current search
  * @param idx, the idx of the current part
  * @param parts, the parts of the pattern
  * @param bpEDv Vector containing an alignment matrix per part
  * @param occ Datastructure with the in-index and in-text occurrences, if an
  * occurrence is found it will be added to this datastructure
  * @param counters the performance counters
  * @param initOther, eds of descendants in other direction
  * @param descOther, descendants in other direction
  * @param remainingDesc, the remaining descendants on the current
  * branch, that are already created but aren't checked yet and might
  * need to be checked for the next part, defaults to an empty vector
  * @return false if the search can continue along this branch for the
  * current part, true if the search can backtrack
  */
    bool branchAndBound(BitParallelED& bpED, Cluster& clus,
                        const FMPosExt& currentNode, const Search& s,
                        const length_t& idx,
                        const std::vector<Substring>& parts, Occurrences& occ,
                        Counters& counters, const std::vector<uint>& initOther,
                        const std::vector<FMPosExt>& descOther,
                        const std::vector<FMPosExt>& remainingDesc = {});

    /**
    * Goes deeper in a search if a valid approximate match is found in the
    * cluster

    * @param cluster, the cluster to search for a valid approximate match
    * @param nIdx, the idx of next part to research
    * @param s, the search
    * @param parts, the parts of the pattern
    * @param occ Datastructure with the in-index and in-text occurrences, if an
    * occurrence is found it will be added to this datastructure
    * @param cnts the performance counters
    * @param descsOtherD, the descendants of the other direction,
    * defaults to empty vector
    * @param initOhterD, the initialization eds of the other direction,
    * defaults to empty vector
    * @param remDesc, the remaining descendants on the current
    * branch, that are already created but aren't checked yet and need to
    * be checked for the next part, defaults to an empty vector
    */
    void goDeeper(Cluster& cluster, const length_t& nIdx, const Search& s,
                  const std::vector<Substring>& parts, Occurrences& occ,
                  Counters& cnts,
                  const std::vector<FMPosExt>& descndansndansOtherD,
                  const std::vector<uint>& intitOtherD,
                  const std::vector<FMPosExt>& remDesc);

    // ----------------------------------------------------------------------------
    // IN TEXT VERIFICATION ROUTINES
    // ----------------------------------------------------------------------------

    /**
     * Verifies the text occurrences in the text and adds them to the
     * occurrences for the edit distance metric.
     * @param startPos the start positions to be verified
     * @param maxED the maximal edit distance that is allowed for the search
     * @param minED the minimal edit distance that is allowed for the search
     * @param occ the occurrences, both those found in the FM Index and in the
     * text, if a valid text occurrence is found it will be added
     * @param counters the performance counters
     */
    void inTextVerification(const std::vector<length_t>& startPos,
                            const length_t& maxED, const length_t& minED,
                            Occurrences& occ, Counters& counters) const;

    /**
     * In text verification for the Hamming distance.
     * @param node the node in the index at which the switch to in-text
     * verification will happen
     * @param s the current search
     * @param parts the parts of the pattern for this search
     * @param idx the current index in the search
     * @param occ Datastructure with all occurrences, both in FM Index and in
     * text. If in-text verification leads to valid text occurrences, these will
     * be added to occ
     */
    void inTextVerificationHamming(const FMPosExt& node, const Search& s,
                                   const std::vector<Substring>& parts,
                                   const length_t idx, Occurrences& occ) const;

    /**
     * Converts a match in the suffix array to matches in the text.
     * @param matchInSA the match that will be converted
     * @returns a vector with the corresponding text occurrences
     */
    std::vector<TextOcc> convertToMatchesInText(const FMOcc& matchInSA) const;

    /**
     * Find the begin positions of the range that needs to be verified in-text
     * @param rangeSA the range over the SA of the partial in-index match that
     * needs to be verified in text
     * @param startDiff the highest possible difference between the start of the
     * partial match and that of a valid full match in the text
     * @param shift the shift of the current node (default = 0)
     * @returns the lowest possible begin position for each of the values in the
     * range
     */
    std::vector<length_t> getBeginPositions(const Range& rangeSA,
                                            length_t startDiff,
                                            length_t shift = 0) const {

        std::vector<length_t> positions(rangeSA.width(), 0);
        const auto &saBegin = rangeSA.getBegin(), saEnd = rangeSA.getEnd();

        for (length_t i = saBegin; i < saEnd; i++) {
            positions[i - saBegin] = findSA(i) + shift - startDiff;
        }

        return positions;
    }

    // ----------------------------------------------------------------------------
    // EXTEND ROUTINES
    // ----------------------------------------------------------------------------

    /**
     * Pushes all the children corresponding to the node with ranges equal
     * to ranges
     * @param ranges the ranges to get the children of
     * @param stack, the stack to push the children on
     * @param counters the performance counters
     * @param row, the row of the parentNode (defaults to 0)
     */
    void extendFMPos(const SARangePair& ranges, std::vector<FMPosExt>& stack,
                     Counters& counters, length_t row = 0) const;

    /**
     * Pushes all the children corresponding to the this position
     * @param pos, the position to get the children of
     * @param stack, the stack to push the children on
     * @param counters the performance counters
     */
    void extendFMPos(const FMPosExt& pos, std::vector<FMPosExt>& stack,
                     Counters& counters) const;

  public:
    // ----------------------------------------------------------------------------
    // INITIALIZATION ROUTINES
    // ----------------------------------------------------------------------------

    /**
     * Constructor
     * @param baseFile baseFile of the files that contain the info
     * @param sa_spase sparseness factor of suffix array. It is assumed this
     * is a power of two
     * @param verbose, will write to cout
     */
    FMIndex(const std::string& baseFile, length_t inTextSwitch,
            int sa_sparse = 1, bool verbose = true, length_t wordSize = 10)
        : baseFile(baseFile), sparseFactorSA(sa_sparse),
          logSparseFactorSA(log2(sa_sparse)), sparseSA(baseFile, sa_sparse),
          inTextSwitchPoint(inTextSwitch), wordSize(wordSize) {
        // read in files
        fromFiles(baseFile, verbose);

        // populate table
        populateTable(verbose);
    }

    /**
     * Get the complete range of this index
     * @returns an SARangePair with both ranges the complete range of the
     * index
     */
    SARangePair getCompleteRange() const {
        return SARangePair(Range(0, bwt.size()), Range(0, bwt.size()));
    }

    // ----------------------------------------------------------------------------
    // ROUTINES FOR ACCESSING THE DATASTRUCTURE
    // ----------------------------------------------------------------------------

    /**
     * Get the original text
     */
    const std::string& getText() const {
        return text;
    }

    const length_t getTextSize() const {
        return text.size();
    }

    /**
     * Get the cross-over point form in-index to in-text verification
     */
    length_t getSwitchPoint() const {
        return inTextSwitchPoint;
    }

    /**
     * @returns the wordsize of the mers stored in the table
     */
    length_t getWordSize() const {
        return wordSize;
    }

    /**
     * @returns the ranges of a single character in this index
     */
    SARangePair getRangeOfSingleChar(char c) const {
        assert(sigma.inAlphabet(c));
        unsigned int i = sigma.c2i(c);
        if (i < sigma.size() - 1) {
            return SARangePair(Range(counts[i], counts[i + 1]),
                               Range(counts[i], counts[i + 1]));
        }
        return SARangePair(Range(counts[i], bwt.size()),
                           Range(counts[i], bwt.size()));
    }

    /**
     * Looks up the SARangePair corresponding to p in the hashtable. Assumes
     * p is of size wordSize
     * @param p, the substring to find the ranges of
     * @returns the ranges corresponding to substring p, if no pair can be
     * found returns empty ranges
     */
    SARangePair lookUpInKmerTable(const Substring& p) const {
        Kmer k(p.tostring());

        auto it = table.find(k);
        if (it != table.end()) {
            return it->second;
        }

        return SARangePair();
    }

    // ----------------------------------------------------------------------------
    // ROUTINES FOR EXACT MATCHING
    // ----------------------------------------------------------------------------

    /**
     * Calculates the exact matches to the string in the index and returns them
     * with their CIGAR string
     * @param s the string to match in the reference genome
     * @param counters the performance counters
     * @returns a sorted vector containing the start positions of all exact
     * substring matches of s in the reference sequence
     */
    std::vector<TextOcc> exactMatchesOutput(const std::string& s,
                                            Counters& counters) const {
        if (s.size() == 0) {
            return {};
        }
        SARangePair r = getCompleteRange();
        assert(dir == FORWARD);

        std::vector<TextOcc> tOcc;

        bool broken = false;
        for (const char& c : s) {
            findRangesWithExtraCharForward(sigma.c2i(c), r, r);
            counters.incNodeCounter();
            if (r.empty() || r.width() <= inTextSwitchPoint) {
                broken = true;
                break;
            }
        }

        // THE CIGAR for the matches
        std::vector<std::pair<char, uint>> CIGAR = {
            std::make_pair('M', s.size())};

        const auto& p = getBeginPositions(r.getRangeSA(), 0);
        if (broken && !r.empty()) {
            // for each of the positions check if the substring in the text
            // equals the string

            for (const auto& pos : p) {
                if (pos + s.size() <= text.size() &&
                    text.compare(pos, s.size(), s) == 0) {
                    tOcc.emplace_back(Range(pos, pos + s.size()), 0, CIGAR);
                }
            }
            counters.inTextStarted += p.size();
            counters.usefulCigarsInText += tOcc.size();
            counters.abortedInTextVerificationCounter += p.size() - tOcc.size();
            counters.cigarsInTextVerification += tOcc.size();
        } else {
            // No verification needed, just add the cigar
            for (const auto& pos : p) {
                tOcc.emplace_back(Range(pos, pos + s.size()), 0, CIGAR);
            }
            counters.cigarsInIndex += p.size();
        }

        counters.totalReportedPositions += tOcc.size();
        // sort the vector and return
        std::sort(tOcc.begin(), tOcc.end());
        return tOcc;
    }

    /**
     * This function matches a string exactly and returns the ranges in the
     * sa and saRev
     * @param pattern the string to match
     * @param counters the performance counters
     * @returns the pair of ranges of this pattern
     */
    SARangePair matchStringBidirectionally(const Substring& pattern,
                                           Counters& counters) const {
        return matchStringBidirectionally(pattern, getCompleteRange(),
                                          counters);
    }

    /**
     * This function matches a string exactly starting form startRange
     * @param pattern the string to match
     * @param startRange, the range to search in
     * @param counters the performance counters
     * @returns the pair of ranges
     */
    SARangePair matchStringBidirectionally(const Substring& pattern,
                                           SARangePair startRange,
                                           Counters& counters) const;

    /**
     * Adds one character and updates the range. If the character can't be
     * added the range will be set to an empty range
     * @param c the character to be added (in the current direction of the
     * index)
     * @param range the range to extend
     * @param counters the performance counters
     */
    bool addChar(const char& c, SARangePair& range, Counters& counters) const;

    // ----------------------------------------------------------------------------
    // PREPROCESSING ROUTINES FOR APPROXIMATE MATCHING
    // ----------------------------------------------------------------------------

    /**
     * Resizes to the required number of stacks and reserves space on each
     * stack, so that each stack can match the entire pattern
     * @param number, the number of stacks required
     * @param size, the size of the pattern
     */
    void reserveStacks(const length_t number, const length_t size) {
        stacks.resize(number);
        length_t stackSize = size * sigma.size();
        for (auto& stack : stacks) {
            stack.reserve(stackSize);
        }
    }
    /**
     * Reset the in-text matrices to be empty matrices
     * @param number the number of partitions
     */
    void resetMatrices(const length_t number) {
        matrices.resize(2 * number);

        for (auto& matrix : matrices) {
            matrix.reset();
        }
    }

    /**
     * Set the sequence for the in-text verification matrix
     * @param s the sequence to be set
     */
    void setInTextMatrixSequence(const Substring& s) {
        inTextMatrix.setSequence(s);
    }

    // ----------------------------------------------------------------------------
    // ROUTINES FOR APPROXIMATE MATCHING
    // ----------------------------------------------------------------------------

    /**
     * Matches the pattern approximately. All matches are at most a certain
     * edit distance away from the pattern
     * @param pattern the pattern to match
     * @param maxED the maximum edit distance
     * @param counters the performance counters
     * @returns a vector with matches which contain a range (the range of
     * the text that matched) and the edit distance this substring is away
     * from the pattern
     */
    const std::vector<TextOcc> approxMatchesNaive(const std::string& pattern,
                                                  length_t maxED,
                                                  Counters& counters);

    /**
     * Sets the search direction of the fm-index
     * @param d the direction to search in, either FORWARD or BACKWARD
     */
    void setDirection(Direction d) {
        dir = d;
        extraChar = (d == FORWARD) ? &FMIndex::findRangesWithExtraCharForward
                                   : &FMIndex::findRangesWithExtraCharBackward;
    }
    /**
     * Matches a search recursively with a depth first approach (each branch
     * of the tree is fully examined until the backtracking condition is
     * met) using hamming distance metric
     * @param search, the search to follow
     * @param startMatch, the approximate match found for all previous parts
     * of the search
     * @param occ, a datastructure with matches of the complete search, if  such
     * a match is found it will be added to this datastructure
     * @param counters the performance counters
     * @param idx, the index of the partition to match, defaults to 1 as an
     * exact search for the zeroth part is assumed
     */
    void recApproxMatchHamming(const Search& s, const FMOcc& startMatch,
                               Occurrences& occ,
                               const std::vector<Substring>& parts,
                               Counters& counters, const int& idx = 1);

    /**
     * Entry to the recusive approximate matching procedure for the edit
     * distance

     * @param search, the search to follow
     * @param startMatch the match containing the SA ranges corresponding to
     * the match of the first part of the search
     * @param occ, a datastructure with matches of the complete search, if  such
     * a match is found it will be added to this datastructure
     * @param parts, the parts of the pattern
     * @param counters the performance counters
     * @param idx, the index of the partition to match, defaults to 1 as an
     * exact search for the zeroth partition is assumed
     */
    void recApproxMatchEditOptimizedEntry(const Search& search,
                                          const FMOcc& startMatch,
                                          Occurrences& occ,
                                          const std::vector<Substring>& parts,
                                          Counters& counters,
                                          const int& idx = 1) {

        if (startMatch.getRanges().width() > inTextSwitchPoint) {
            counters.approximateSearchStarted++;
            recApproxMatchEditOptimized(search, startMatch, occ, parts,
                                        counters, idx);
            return;
        }
        // verify the partial match in text
        verifyExactPartialMatchInText(
            startMatch, parts[search.getLowestPartProcessedBefore(idx)].begin(),
            search.getMaxED(), occ, counters);
    }

    /**
     * Verifies an exact partial match in the text for all occurrences of that
     * exact partial match

     * @param startMatch the match containing the SA ranges corresponding to
     * this exact match and depth of the exact match
     * @param beginInPattern the begin position of the exact part in the pattern
     * to be searched
     * @param maxED the maximal allowed edit distance
     * @param occ the occurrences, to this vector new in-text verified
     * occurrences can be added during the function
     * @param counters the performace counters
     */
    void verifyExactPartialMatchInText(const FMOcc& startMatch,
                                       const length_t& beginInPattern,
                                       const length_t& maxED, Occurrences& occ,
                                       Counters& counters) {

        // Immediately switch to in-text verification
        counters.immediateSwitch++;

        // A) find out highest possible difference
        length_t startDiff = beginInPattern + maxED;

        // B) get the possible begin positions
        const auto& partialStarts =
            getBeginPositions(startMatch.getRanges().getRangeSA(), startDiff);

        // C) verify these positions in the text
        inTextVerification(partialStarts, maxED, 0, occ, counters);
    }

    /**
     * Matches a search recursively with a depth first approach (each branch
     * of the tree is fully examined until the backtracking condition is
     * met) using the edit distance metric. This function does not use any
     * optimizations for eliminating redundancy in the edit distance metric.
     * It simply matches the current part starting from startrange and each
     * node found that has an edit distance between the lower and upper
     * bound is used to start a search for the next part
     * @param s, the search to follow
     * @param startMatch, the approximate match found for all previous
     * partitions of the search
     * @param occ, a datastructure with matches of the complete search, if
     * such a match is found it will be added to this datastructure
     * @param parts, the parts of the pattern
     * @param counters the performance counters
     * @param idx, the index of the partition to match, defaults to 1 as an
     * exact search for the zeroth partition is assumed
     */
    void recApproxMatchEditNaive(const Search& s, const FMOcc& startMatch,
                                 Occurrences& occ,
                                 const std::vector<Substring>& parts,
                                 Counters& counters, const int& idx);
    // ----------------------------------------------------------------------------
    // POST-PROCESSING ROUTINES FOR APPROXIMATE MATCHING
    // ----------------------------------------------------------------------------

    /**
     * Finds the entry in the suffix array of this index. This is
     * computed from the sparse suffix array and the bwt
     * @param index the index to find the entry in the SA off
     * @returns the entry in the SA of the index
     */
    length_t findSA(length_t index) const;

    /**
     * Get a reference subsequence
     * @param b the start position of the subsequence
     * @param e the end position of the subsequence
     * @returns the reference subsequence
     */
    Substring getSubstring(const length_t b, const length_t e) const {
        return Substring(text, b, e);
    }

    /**
     * Get a reference subsequence
     * @param range the range of the subsequence
     * @returns the reference subsequence
     */
    Substring getSubstring(const Range& range) const {
        return getSubstring(range.getBegin(), range.getEnd());
    }

    /**
     * Find unqiue text occurrences for the hamming distance
     */
    std::vector<TextOcc> getTextOccHamming(Occurrences& occ,
                                           Counters& counters) const {
        counters.totalReportedPositions += occ.textOccSize();
        // erase in-index doubles
        occ.eraseDoublesFM();

        // Get the in-index occurrences
        const auto& fmoccs = occ.getFMOccurrences();

        // The size of the match
        length_t size = (fmoccs.empty()) ? 0 : fmoccs[1].getDepth();

        for (const auto& fmOcc : fmoccs) {
            // Get the range
            const Range& saRange = fmOcc.getRanges().getRangeSA();
            counters.totalReportedPositions += saRange.width();

            for (length_t i = saRange.getBegin(); i < saRange.getEnd(); i++) {
                length_t b = findSA(i);
                std::vector<std::pair<char, uint>> CIGAR = {
                    std::make_pair('M', size)};
                occ.addTextOcc(Range(b, b + size), fmOcc.getDistance(), CIGAR);
            }
        }

        // remove doubles
        occ.eraseDoublesText();
        // Generate string output
        occ.generateOutput();

        return occ.getTextOccurrences();
    }

    /**
     * Find the unique text occurrences
     * @param occ occurrences, both in-text and in-index
     * @param maxED the maximal allowed edit distance
     * @param counters performance counters
     */
    std::vector<TextOcc> getUniqueTextOccurrences(Occurrences& occ,
                                                  const length_t& maxED,
                                                  Counters& counters);
};

#endif
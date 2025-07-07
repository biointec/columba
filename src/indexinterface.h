/******************************************************************************
 *  Columba: Approximate Pattern Matching using Search Schemes                *
 *  Copyright (C) 2020-2022 - Luca Renders <luca.renders@ugent.be> and        *
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
#ifndef INDEXINTERFACE_H
#define INDEXINTERFACE_H

#include "alphabet.h"
#include "bitparallelmatrix.h" // for BitParallelED64
#include "definitions.h"       // for length_t, PairStatus, Strand
#include "indexhelpers.h"      // for SARangePair, Range, TextOcc
#include "reads.h"             // for ReadBundle
#include "substring.h"         // for Substring
#include "tkmer.h"             // for Kmer, KmerHash

#include <algorithm>                // for max
#include <parallel_hashmap/phmap.h> // phmap::parallel_flat_hash_map
#include <stddef.h>                 // for size_t
#include <stdint.h>                 // for uint16_t
#include <string>                   // for string, char_traits<>::pos_type
#include <utility>                  // for pair
#include <vector>                   // for vector
                                    // RUN_LENGTH_COMPRESSION

class MemoryMappedTextFile;
class Search;

class IndexInterface;
typedef bool (IndexInterface::*ExtraCharPtr)(length_t, const SARangePair&,
                                             SARangePair&) const;

/**
 * Class representing the bidirectional index and all operations on that
 * Index
 */
class IndexInterface {
  protected:
    // info about the text
    const std::string baseFile; // The base file of the reference text
    length_t textLength;        // the length of the text

    std::vector<length_t> counts; // the counts array of the reference genome

    Alphabet<ALPHABET> sigma; // the alphabet

    // direction variables
    thread_local static Direction dir; // the direction of the index
    thread_local static ExtraCharPtr
        extraChar; // pointer to extra char method (for direction)
                   // RUN_LENGTH_COMPRESSION

    // stacks for search schemes
    thread_local static std::vector<std::vector<FMPosExt>>
        stacks; // stacks of nodes for the different searches
    thread_local static std::vector<BitParallelED64>
        matrices; // alignment matrices for the different parts

    thread_local static std::vector<BitParallelED128>
        matrices128; // big alignment matrices for the different parts

    // info about the current read/search
    thread_local static Strand
        strand; // The strand on which the current search is executed. If strand
                // is REVERSE_C_STRAND, then the read was reverse complemented
                // before starting the search and all found occurrences should
                // be labeled as reverse complemented.
    thread_local static PairStatus
        pairStatus; // The status of the read that is currently being searched
                    // in the pair.

    // sparse hash info
    const size_t wordSize = 10; // the size of the mers to be stored in a table
    phmap::parallel_flat_hash_map<Kmer, SARangePair, KmerHash>
        table; // hashtable that contains all wordSize-mers

    // assignment/SAM info
    std::vector<length_t> startPos;    // the start positions of the sequences
    std::vector<std::string> seqNames; // the names of the sequences
    std::vector<length_t> firstSeqIDPerFile; // the first seqID per FASTA file

    bool noCIGAR = false; // do not calculate the CIGAR string

    // in-text verification matrices
    thread_local static std::vector<BitParallelED64>
        fullReadMatrices; // the in-text verification matrices
    thread_local static BitParallelED64*
        fullReadMatrix; // pointer to the current fullReadMatrix

    thread_local static std::vector<BitParallelED128>
        fullReadMatrices128; // the in-text verification matrices for many
                             // errors/start positions
    thread_local static BitParallelED128*
        fullReadMatrix128; // pointer to the current fullReadMatrix for many
                           // errors/start positions

    // ----------------------------------------------------------------------------
    // PREPROCESSING ROUTINES
    // ----------------------------------------------------------------------------

    /**
     * Read the counts array, the build_tag info, the flavor and the
     * compilation-bits info of the index.
     * @param baseFile The base name of the files that contain the info about
     * the index.
     * @param verbose If true, the steps will be written to cout.
     */
    void readMetaAndCounts(const std::string& baseFile, bool verbose);

    /**
     * Read the sequence names and start positions of the sequences in the
     * index.
     * @param baseFile The base name of the files that contain the info about
     * the index.
     * @param verbose If true, the steps will be written to cout.
     */
    void readSequenceNamesAndPositions(const std::string& baseFile,
                                       bool verbose);

    /**
     * Read a binary file and stores content in array
     * @param filename File name
     * @param array Suffix array (contents will be overwritten)
     * @returns True if successful, false otherwise
     */
    static bool readArray(const std::string& filename,
                          std::vector<length_t>& array);

    /**
     * Populate the hash table
     * @param verbose if steps are written to cout
     */
    void populateTable(bool verbose);

    // ----------------------------------------------------------------------------
    // ROUTINES FOR ACCESSING DATA STRUCTURE
    // ----------------------------------------------------------------------------

    /**
     * Get the start positions of the sequences
     * @return The start positions vector
     */
    const std::vector<length_t>& getStartPositions() const {
        return startPos;
    }

    // ----------------------------------------------------------------------------
    //  APPROXIMATE PATTERN MATCHING
    // ----------------------------------------------------------------------------
    /**
     * Helper function for approximate matching with the edit distance metric.
     * This function crosses over from in-index matching to in-text
     * verification.
     * @param node The current node in the index.
     * @param s The current search.
     * @param parts The parts of the pattern.
     * @param occ The data structure with the in-index and in-text occurrences,
     * if an in-text occurrence is found it will be added to this data
     * structure.
     * @param counters The performance counters.
     * @param pattern The pattern to verify.
     * @param idx The index in the search
     * @param bpED The bit-parallel alignment matrix of the current part.
     * @param sMatch The approximate match found for the previous part
     * @param dOther The descendants of the other direction.
     * @param iOther The initialization eds of the other direction.
     */
    void goToInTextVerificationEdit(const FMPosExt& node, const Search& s,
                                    const std::vector<Substring>& parts,
                                    Occurrences& occ, Counters& counters,
                                    const Substring& pattern, length_t idx,
                                    IBitParallelED* bpED, const FMOcc& sMatch,
                                    const std::vector<FMPosExt>& dOther,
                                    const std::vector<uint16_t>& iOther);
    /**
     * Matches a search recursively with a depth first approach (each branch
     * of the tree is fully examined until the backtracking condition is
     * met) using the edit distance metric. This function uses all
     * optimizations for eliminating redundancy in the edit distance metric
     * @param search the search to follow
     * @param startMatch the approximate match found for all previous
     * partitions of the search
     * @param occ a data structure with matches of the complete search, if
     * such a match is found it will be added to this data structure
     * @param parts the parts of the pattern
     * @param counters the performance counters
     * @param idx the index of the partition to match, defaults to 1 as an
     * exact search for the zeroth partition is assumed
     * @param descPrevDir the descendants of the previous direction,
     * defaults to empty vector
     * @param initPrevDir the initialization eds of the previous direction,
     * defaults to empty vector
     * @param descNotPrevDir the descendants of the other direction,
     * defaults to empty vector
     * @param initNotPrevDir the initialization eds of the other direction,
     * defaults to empty vector
     */
    void recApproxMatchEdit(
        const Search& search, const FMOcc& startMatch, Occurrences& occ,
        const std::vector<Substring>& parts, Counters& counters,
        const int& idx = 1,
        const std::vector<FMPosExt>& descPrevDir = std::vector<FMPosExt>(),
        const std::vector<uint16_t>& initPrevDir = std::vector<uint16_t>(),
        const std::vector<FMPosExt>& descNotPrevDir = std::vector<FMPosExt>(),
        const std::vector<uint16_t>& initNotPrevDir = std::vector<uint16_t>());

    /**
     * @brief Update the runs corresponding to the SA range
     *
     * @param ranges SARangePair for which the SA range runs must be updated
     */
    virtual void updateRangeSARuns(SARangePair& ranges) const = 0;

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
                                    SARangePair& childRanges) const = 0;

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
                                   SARangePair& childRanges) const = 0;

    /**
     * Finds the range of Pc using unidirectional backward matching. With
     * toeholds in the case of run-length compression
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
                                   SARangeBackwards& childRange

    ) const = 0;

    /**
     * Finds the range of cP and stores a 'dummy' range over the suffix array of
     * the reversed text
     * @param positionInAlphabet the position in alphabet of the character
     * that is added in the front.
     * @param rangesOfP the ranges of pattern P, only the backwards range is
     * used
     * @param rangesOfChild the ranges corresponding to string cP, this
     * will be set during execution. The getRangeSARev() function will return a
     * dummy!
     * @returns true if the range is not empty, false otherwise
     */
    virtual bool findRangesWithExtraCharBackwardUniDirectional(
        length_t positionInAlphabet, const SARangePair& rangesOfP,
        SARangePair& rangesOfChild) const = 0;

#ifdef RUN_LENGTH_COMPRESSION

    void updateMatchStr(std::vector<char>& matchStr, char c, length_t row,
                        length_t startSize) const {
        assert(c != 0 && "Cannot update match string with a null character.");
        size_t numToKeep = startSize + (row - 1);
        assert(matchStr.size() >= numToKeep);
        matchStr.resize(numToKeep);
        matchStr.push_back(c);

        return;
    }
#endif // RUN_LENGTH_COMPRESSION

    /**
     * Helper function for the approximate matching. This function fills in
     * the matrix for the current node at the current row and goes deeper
     * for the next part if that is necessary. The function returns true if
     * the search should backtrack.
     * @param bpED The pointer to the bit-parallel alignment matrix of the
     * current part.
     * @param cluster The cluster corresponding to the final column of the
     * matrix.
     * @param currentNode The node for which the matrix is filled in.
     * @param s The current search.
     * @param idx The index of the current part in the search.
     * @param parts The parts of the pattern.
     * @param occ Data structure with the in-index and in-text occurrences,
     * if an occurrence, either in-index or in-text, is found it will be
     * added to this data structure.
     * @param counters The performance counters.
     * @param initOther Edit distance values of known descendants in the
     * other direction.
     * @param descOther The known descendants in the other direction.
     * @param remainingDesc The remaining descendants on the current branch,
     * that are already created but aren't checked yet and might need to be
     * checked for the next part, defaults to an empty vector.
     * @return false if the search can continue along this branch for the
     * current part, true if the search should backtrack.
     */
    bool branchAndBound(IBitParallelED* bpED, MatrixMetaInfo& cluster,
                        const FMPosExt& currentNode, const Search& s,
                        const length_t& idx,
                        const std::vector<Substring>& parts, Occurrences& occ,
                        Counters& counters,
                        const std::vector<uint16_t>& initOther,
                        const std::vector<FMPosExt>& descOther,
                        const std::vector<FMPosExt>& remainingDesc = {});

    /**
     * Goes deeper in a search if a valid partial approximate match is found in
     * the cluster. This function is called by the recursive function's helper
     * branchAndBound function..
     * @param cluster The cluster corresponding to the final column of the
     * matrix.
     * @param nIdx The index of next part to in search s.
     * @param s The search.
     * @param parts The parts of the pattern.
     * @param occ Data structure with the in-index and in-text occurrences, if
     * an occurrence, either in-index or in-text, is found it will be added to
     * this data structure.
     * @param counters The performance counters.
     * @param descOtherD The known descendants in the other direction.
     * @param initOtherD Edit distance values of known descendants in the other
     * direction.
     * @param remDesc The remaining descendants on the current branch,
     * that are already created but aren't checked yet and might need to be
     * checked for the next part.
     */
    void goDeeper(MatrixMetaInfo& cluster, const length_t& nIdx,
                  const Search& s, const std::vector<Substring>& parts,
                  Occurrences& occ, Counters& counters,
                  const std::vector<FMPosExt>& descOtherD,
                  const std::vector<uint16_t>& initOtherD,
                  const std::vector<FMPosExt>& remDesc);

    // ----------------------------------------------------------------------------
    // IN TEXT VERIFICATION ROUTINES
    // ----------------------------------------------------------------------------

    /**
     * Verifies the text occurrences in the text and adds them to the
     * occurrences for the edit distance metric. WARNING: This function makes
     * use of the fullReadMatrix pointer. Before calling this function make sure
     * the pointer has been set correctly (using the setIndexInMode function).
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
                                    bool fixedStartPos) const = 0;

    /**
     * Verifies an in-text occurrence and adds it to the occurrences for the
     * edit distance metric. WARNING: This function makes use of the
     * fullReadMatrix pointer. Before calling this function make sure the
     * pointer has been set correctly (using the setIndexInMode function).
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
                                const std::string& pattern) const = 0;
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
     * @param counters The performance counters.
     */
    virtual void inTextVerificationHamming(const FMPosExt& node,
                                           const Search& s,
                                           const std::vector<Substring>& parts,
                                           const length_t idx, Occurrences& occ,
                                           Counters& counters) const = 0;

    virtual void inTextVerificationHamming(
        const Range& r, const Substring& pattern, const length_t maxEDFull,
        const length_t minEDFull, const length_t lengthBefore, Occurrences& occ,
        Counters& counters) const = 0;

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
    getBeginPositions(const SARangeBackwards& rangeSA, length_t startDiff,
                      length_t shift = 0) const = 0;

    // ----------------------------------------------------------------------------
    // EXTEND ROUTINES
    // ----------------------------------------------------------------------------

    /**
     * Creates the children nodes of the node with the given ranges and pushes
     * them to the stack.
     * @param ranges The ranges of the parent node to get the children of.
     * @param stack The stack to push the children on.
     * @param counters The performance counters.
     * @param row The row in the current matrix of the parentNode (defaults to
     * 0).
     */
    void extendFMPos(const SARangePair& ranges, std::vector<FMPosExt>& stack,
                     Counters& counters, length_t row = 0) const;

    /**
     * Creates the children nodes of the node with the given position and pushes
     * them to the stack.
     * @param pos The position of the parent node to get the children of.
     * @param stack The stack to push the children on.
     * @param counters The performance counters.
     */
    void extendFMPos(const FMPos& pos, std::vector<FMPosExt>& stack,
                     Counters& counters) const;

    // ----------------------------------------------------------------------------
    // SUBROUTINES FOR PATTERN MATCHING
    // ----------------------------------------------------------------------------

    /**
     * Sets the index in the correct mode. All found occurrences will now be
     * be labeled as found along this strand and as an occurrence of this
     * read in the pair.
     * @param strand On which strand the current read lies
     * @param pairStatus The status of the current read in its pair. [default =
     * FIRST_IN_PAIR]
     */
    void setIndexInModeSubRoutine(Strand strand,
                                  PairStatus pairStatus = FIRST_IN_PAIR) {
        this->strand = strand;
        this->pairStatus = pairStatus;
    }

  public:
    // ----------------------------------------------------------------------------
    // INITIALIZATION ROUTINES
    // ----------------------------------------------------------------------------

    /**
     * Constructor
     * @param baseFile The base name of the files that contain the info about
     * the index.
     * @param verbose If true, the steps will be written to cout. [default =
     * true]
     * @param noCIGAR If true, the CIGAR string will not be calculated.
     * @param wordSize The size of the mers to be stored in the hashtable. Used
     * for quick look-ups of exact seeds. [default = 10]
     */
    IndexInterface(const std::string& baseFile, bool verbose = true,
                   bool noCIGAR = false, length_t wordSize = 10)
        : baseFile(baseFile), wordSize(wordSize), noCIGAR(noCIGAR) {
        readMetaAndCounts(baseFile, verbose);
    }

    /**
     * @brief Destructor
     *
     */
    virtual ~IndexInterface() = default;

    /**
     * Get the complete range of this index corresponding to an exact match of
     * the empty string.
     * @returns A SARangePair with both ranges the complete range of the
     * index.
     */
    virtual SARangePair getCompleteRange() const = 0;

    /**
     * @brief Get the Empty String FMPos object
     *
     * @return FMPosExt
     */
    virtual FMPos getEmptyStringFMPos() const = 0;

    // ----------------------------------------------------------------------------
    // ROUTINES FOR ACCESSING THE DATA STRUCTURE
    // ----------------------------------------------------------------------------

    /**
     * Get the original text
     */
    virtual const std::string& getText() const = 0;

    /**
     * Get the length of the text
     */
    const length_t getTextSize() const {
        return textLength;
    }

    /**
     * @returns the wordsize of the mers stored in the table
     */
    length_t getWordSize() const {
        return wordSize;
    }

    /**
     * Get the cross-over point form in-index to in-text verification
     */
    virtual length_t getSwitchPoint() const = 0;

    /**
     * @brief Check if the CIGAR string must be reported
     */
    bool getNoCIGAR() const {
        return noCIGAR;
    }

    /**
     * Find the range of an exact match of single character in this index.
     * Warning: this function assumes that the character is in the alphabet.
     * @returns the ranges of a single character in this index
     */
    virtual SARangePair getRangeOfSingleChar(char c) const = 0;

    /**
     * Looks up the SARangePair of the exact match of c p in the hashtable.
     * Warning: This function assumes p is of size wordSize.
     * @param p The substring to find the ranges of.
     * @returns the ranges corresponding to substring p, if no pair can be
     * found it returns empty ranges.
     */
    SARangePair lookUpInKmerTable(const Substring& p) const {

        auto it = table.find(Kmer(p.getText(), p.begin()));
        return (it == table.end() || p.containsN()) ? SARangePair()
                                                    : it->second;
    }

    const std::vector<std::string>& getSeqNames() const {
        return seqNames;
    }

    /**
     * Finds the sequence name to which the given text-occurrence belongs. Also
     * changes the range of the text-occurrence to be relative to this sequence.
     * @param t The text-occurrence to find the sequence name of.
     * @param seqID the unique ID of the sequence to which the text-occurrence
     * belongs (if a sequence could be found)
     * @param counters The performance counters.
     * @param largestStratum the largest stratum of distance that is allowed.
     * @param metric the used distance metric. In case of hamming distance no
     * trimming will be performed.
     * @param pattern The pattern that was matched, needed for trimming
     * @returns SeqNameFound::FOUND if the name was found without
     * trimming, SeqNameFound::FOUND_WITH_TRIMMING if the name was found with
     * trimming and SeqNameFound::NOT_FOUND if no name could be found (i.e. t
     * spans two references and trimming cannot fix it).
     */
    SeqNameFound findSeqName(TextOcc& t, length_t& seqID, Counters& counters,
                             length_t largestStratum,
                             const DistanceMetric& metric,
                             const std::string& pattern) const;

    inline bool isInFirstFile(const TextOcc& occ) const {
        return firstSeqIDPerFile.size() == 1 ||
               occ.getEnd() < startPos[firstSeqIDPerFile[1]];
    }
    // ----------------------------------------------------------------------------
    // ROUTINES FOR EXACT MATCHING
    // ----------------------------------------------------------------------------

    /**
     * Calculates the exact matches to the string in the index. Creates
     * corresponding text-occurrences and adds them to the output vector.
     * Make sure to set the index in the correct mode before calling this
     * function to correctly label the found occurrences.
     * @param s The string to match in the reference genome.
     * @param counters The performance counters.
     * @param tOcc The vector to which the found text occurrences will be added.
     */
    void exactMatchesOutput(const std::string& s, Counters& counters,

                            std::vector<TextOcc>& tOcc) const;

    /**
     * This function matches a string exactly and returns the ranges in the
     * suffix array of the text and of the reverse text. If the string cannot be
     * matched empty ranges will be returned.
     * Warning: this function assumes the direction of the index
     * corresponds to the direction of pattern!
     * @param pattern The string to match.
     * @param counters The performance counters.
     * @returns The pair of ranges corresponding to an exact match of this
     * pattern.
     */
    SARangePair matchStringBidirectionally(const Substring& pattern,
                                           Counters& counters) const {
        return matchStringBidirectionally(pattern, getCompleteRange(),
                                          counters);
    }

    /**
     * This function matches a string exactly starting form a previous start
     * range.
     * Warning: this function assumes the direction of the index corresponds to
     * the direction of pattern!
     * @param pattern The string to match.
     * @param startRange The range to start the match from.
     * @param counters The performance counters.
     * @returns The pair of ranges corresponding to an exact match of this
     * pattern.
     */
    SARangePair matchStringBidirectionally(const Substring& pattern,
                                           SARangePair startRange,
                                           Counters& counters) const;

    /**
     * Adds one character and updates the range. If the character cannot be
     * added the range will be set to an empty range.
     * Warning: this function assumes the direction of the index is set
     * correctly!
     * @param c The character to be added (in the current direction of the
     * index).
     * @param range The range to update.
     * @param counters The performance counters.
     */
    bool addChar(const char& c, SARangePair& range, Counters& counters) const;

    // ----------------------------------------------------------------------------
    // PREPROCESSING ROUTINES FOR APPROXIMATE MATCHING
    // ----------------------------------------------------------------------------

    /**
     * Resizes the stacks vector to the required number of stacks and reserves
     * space on each stack, so that each stack can match the entire pattern
     * @param number The number of stacks required. This is equal to the number
     * of parts.
     * @param size The size of the pattern.
     */
    void reserveStacks(const length_t number, const length_t size) {
        stacks.resize(number);
        length_t stackSize = size * sigma.size();
        for (auto& stack : stacks) {
            stack.reserve(stackSize);
        }
    }
    /**
     * Reset the alignment matrices for the different partitions.
     * @param number The number of partitions.
     */
    void resetMatrices(const length_t number) {
        matrices.resize(2 * number);
        matrices128.resize(2 * number);

        for (auto& matrix : matrices) {
            matrix.reset();
        }
        for (auto& matrix : matrices128) {
            matrix.reset();
        }
    }

    // ----------------------------------------------------------------------------
    // ROUTINES FOR APPROXIMATE MATCHING
    // ----------------------------------------------------------------------------

    /**
     * Matches the pattern approximately via a very naive backtracking approach.
     * All found matches are at most a certain edit distance away from the
     * pattern and will be added to the return vector.
     * @param pattern the pattern to match
     * @param maxED the maximum edit distance
     * @param counters the performance counters
     * @param occ The vector to which the matches will be added.
     */
    void approxMatchesNaive(const std::string& pattern, length_t maxED,
                            Counters& counters, std::vector<TextOcc>& occ);

    void approxMatchesNaiveHamming(const std::string& pattern, length_t maxED,
                                   Counters& counters,
                                   std::vector<TextOcc>& occ);

    /**
     * Get the in-Text verification matrix for the read on the given strand with
     * the given pair status.
     * @param strand The strand to align to.
     * @param pairStatus The pair status of the read (first or second).
     */
    BitParallelED64& getFullReadMatrix(Strand strand,
                                       PairStatus pairStatus) const {
        return fullReadMatrices[strand * 2 + pairStatus];
    }

    BitParallelED128& getFullReadMatrix128(Strand strand,
                                           PairStatus pairStatus) const {
        return fullReadMatrices128[strand * 2 + pairStatus];
    }

    void resetFullReadMatrices() {
        for (auto& matrix : fullReadMatrices) {
            matrix.reset();
        }
        for (auto& matrix : fullReadMatrices128) {
            matrix.reset();
        }
    }

    /**
     * Sets the search direction of the fm-index
     * @param d the direction to search in, either FORWARD or BACKWARD
     * @param isUniDirectionalBackward if the search is unidirectional backward
     */
    void setDirection(Direction d, bool isUniDirectionalBackward) {
        dir = d;
        extraChar =
            (isUniDirectionalBackward)
                ? &IndexInterface::findRangesWithExtraCharBackwardUniDirectional
                : ((d == FORWARD)
                       ? &IndexInterface::findRangesWithExtraCharForward
                       : &IndexInterface::findRangesWithExtraCharBackward);
    }

    /**
     * Sets the index in the correct mode. All found occurrences will now be
     * be labeled as found along this strand and as an occurrence of this
     * read in the pair.
     * @param strand On which strand the current read lies
     * @param pairStatus The status of the current read in its pair. [default =
     * FIRST_IN_PAIR]
     */
    virtual void setIndexInMode(Strand strand,
                                PairStatus pairStatus = FIRST_IN_PAIR) = 0;

    /**
     * Matches a search recursively with a depth first approach (each branch
     * of the tree is fully examined until the backtracking condition is
     * met) using the hamming distance metric.
     * @param s The search to follow.
     * @param startMatch The partial approximate match found for the previous
     * parts of the search
     * @param occ A data structure with approximate occurrences of the complete
     * pattern. Both in-index and in-text occurrences are stored. If a new
     * approximate occurrence is found, either in-index or in-text it will be
     * added to this data structure.
     * @param parts The parts of the pattern.
     * @param counters The performance counters.
     * @param idx The index of current part in the search, defaults to 1 as an
     * exact search for the zeroth part is assumed.
     */
    void recApproxMatchHamming(const Search& s, const FMOcc& startMatch,
                               Occurrences& occ,
                               const std::vector<Substring>& parts,
                               Counters& counters, const int& idx = 1);

    /**
     * Entry to the recursive approximate matching procedure for the edit
     * distance metric.
     * @param search The search to follow.
     * @param startMatch The match containing the SA ranges corresponding to
     * the match of the first part of the search
     * @param occ A data structure with approximate occurrences of the complete
     * pattern. Both in-index and in-text occurrences are stored. If a new
     * approximate occurrence is found, either in-index or in-text it will be
     * added to this data structure.
     * @param parts The parts of the pattern.
     * @param counters The performance counters.
     * @param idx The index of current part in the search, defaults to 1 as an
     * exact search for the zeroth part is assumed.
     */
    void recApproxMatchEditEntry(const Search& search, const FMOcc& startMatch,
                                 Occurrences& occ,
                                 const std::vector<Substring>& parts,
                                 Counters& counters, const int& idx = 1);

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
        length_t minED, const Substring& pattern) = 0;

    /**
     * Verifies an exact partial match in the text for all occurrences of that
     * exact partial match using the Hamming distance. If a valid match is found
     * it will be added to the occurrences.
     * @param startMatch The match containing the SA ranges corresponding to
     * the exact partial match  to start from and depth of the exact match.
     * @param beginInPattern The position in the pattern where the partial match
     * starts.
     * @param maxD The maximal allowed hamming distance.
     * @param parts The parts of the pattern.
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
        Counters& counters, length_t minD) const = 0;

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
    virtual length_t findSA(length_t index) const = 0;

    /**
     * @brief Get the text positions corresponding to a suffix array range
     *
     * @param ranges The ranges in the suffix array
     * @param positions The vector to which the text positions will be added
     */
    virtual void
    getTextPositionsFromSARange(const SARangePair& ranges,
                                std::vector<length_t>& positions) const = 0;

    /**
     * Find unique text occurrences for the hamming distance. The in-index
     * occurrences are converted to in-text occurrences with CIGAR string.
     * Afterwards redundant occurrences are removed. Eventually all
     * non-redundant text-occurrences are returned.
     * @param occ The occurrences to find the unique text occurrences for. Will
     * be invalidated during execution.
     * @param counters The performance counters.
     * @returns the non-redundant text occurrences.
     */
    std::vector<TextOcc> getTextOccHamming(Occurrences& occ,
                                           Counters& counters) const;

    /**
     * Find unique text occurrences for the edit distance. The in-index
     * occurrences are converted to in-text occurrences.
     * Afterwards redundant occurrences are removed. For all non-redundant
     * in-text occurrences that do not have a CIGAR string a CIGAR string is
     * calculated. Eventually all non-redundant text-occurrences are returned
     * sorted.
     * @param occ The occurrences to find the unique text occurrences for.
     * @param maxED The maximal edit distance that is allowed for the search.
     * @param counters The performance counters.
     * @param bundle The read bundle with info about the read/reverse complement
     * @returns the non-redundant text occurrences.
     */
    std::vector<TextOcc>
    getUniqueTextOccurrences(Occurrences& occ, const length_t& maxED,
                             Counters& counters,
                             const ReadBundle& bundle) const;

    /**
     * Generate the CIGAR strings for the in-index occurrences that do not have
     * one yet. This function assumes that the in-text verification matrices
     * have correctly been set.
     * @param occs The occurrences to generate the CIGAR strings for.
     * @param counters The performance counters.
     * @param bundle The read bundle with info about the read
     */
    void generateCIGARS(std::vector<TextOcc>& occs, Counters& counters,
                        const ReadBundle& bundle) {
        if (noCIGAR) {
            return;
        }
        for (TextOcc& t : occs) {
            if (!t.hasCigar()) {
                generateCIGAR(t, counters, bundle.getSequence(t.getStrand()));
            }
        }
    }

    /**
     * Generate the CIGAR string for the in-index occurrence. This function
     * assumes that the in-text verification matrix associated with the strand
     * and pairStatus of the occurrence has been correctly set.
     * @param t The occurrence to generate the CIGAR string for.
     * @param counters The performance counters.
     * @param read The read to which t is a match
     */
    void generateCIGAR(TextOcc& t, Counters& counters,
                       const Substring& read) const {
        if (noCIGAR || t.hasCigar()) {
            return;
        }

// aliases
#ifdef RUN_LENGTH_COMPRESSION
        const std::vector<char>& ref = t.getMatchedStr();
#else
        const auto& ref = Substring(getText(), t.getRange().getBegin(),
                                    t.getRange().getEnd());
#endif
        const auto& score = t.getDistance();
        const auto& strand = t.getStrand();
        const auto& status = t.getPairStatus();

        // select the correct matrix
        IBitParallelED* matrix;
        if (BitParallelED64::getMatrixMaxED() >= score) {
            matrix = &getFullReadMatrix(strand, status);
        } else {
            matrix = &getFullReadMatrix128(strand, status);
        }

        // set the sequence if needed
        if (!matrix->sequenceSet()) {
            assert(read.getDirection() == FORWARD);
            matrix->setSequence(read);
        }

        // find the CIGAR string
        std::string strCIGAR;
        matrix->findCIGAR(ref, score, strCIGAR);
        t.setCigar(strCIGAR);

        // increase the counter
        counters.inc(Counters::CIGARS_IN_INDEX);
    }
};

#endif // INDEXINTERFACE_H
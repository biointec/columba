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

#ifndef SEARCHSTRATEGY_H
#define SEARCHSTRATEGY_H

#include "definitions.h"    // for length_t, MAX_K, MappingMode, Distance...
#include "indexhelpers.h"   // for TextOcc, PairedTextOccs, Occurrences
#include "indexinterface.h" // for the IndexInterface
#include "reads.h"          // for ReadPair, ReadBundle
#include "search.h"         // for Search, SearchScheme
#include "substring.h"      // for Substring

#include <algorithm> // for max, move, min
#include <cassert>   // for assert
#include <cstdint>   // for uint16_t
#include <ios>       // for ifstream, basic_ios
#include <iterator>  // for move_iterator, back_insert_iterator
#include <memory>    // for allocator, allocator_traits<>::value_type
#include <numeric>   // for accumulate
#include <set>       // for set
#include <stdexcept> // for runtime_error, invalid_argument
#include <string>    // for string, operator+, char_traits, to_string
#include <utility>   // for pair, move
#include <vector>    // for vector, _Bit_iterator

#ifdef _WIN32
#include <Windows.h>
#include <direct.h> // For _mkdir on Windows
#define PATH_SEPARATOR "\\"
#else
#include <sys/stat.h> // for stat on POSIX systems
#include <unistd.h>   // For access on POSIX systems
#define PATH_SEPARATOR "/"
#endif

#define Distribution std::vector<int>

// ============================================================================
// CLASS SEARCHSTRATEGY
// ============================================================================

// This is an abstract class. Every derived class should be able to create
// searches for a given value of k. This abstract base class handles the
// partitioning (either with values provided in the derived class or default
// uniform values) and approximate matching (either hamming or edit distance)
class SearchStrategy;

// Pointer to a partition function
typedef void (SearchStrategy::*PartitionPtr)(const std::string&,
                                             std::vector<Substring>&,
                                             const int&, const int&,
                                             std::vector<SARangePair>&,
                                             Counters&) const;

// Pointer to function that starts the index on a particular search
typedef void (SearchStrategy::*StartIdxPtr)(const Search&, const FMOcc&,
                                            Occurrences&,
                                            const std::vector<Substring>&,
                                            Counters&, const int&) const;

// Pointer to match function (ALL or BEST scenario)
typedef void (SearchStrategy::*MatchPtrSE)(ReadBundle&, const length_t,
                                           Counters&, std::vector<TextOcc>&);

typedef std::vector<PairedTextOccs> (SearchStrategy::*MatchPtrPE)(
    ReadPair&, length_t, length_t, length_t, Counters&, std::vector<TextOcc>&);

// Pointer to function for pairing SE alignments (in ALL or BEST scenario's)
typedef std::vector<PairedTextOccs> (SearchStrategy::*PairSEPtr)(
    std::vector<TextOcc>&, std::vector<TextOcc>&, length_t, length_t, Counters&,
    ReadPair&, length_t, std::vector<TextOcc>&, bool);

// pointers to functions for processing combinations according to orientation
typedef void (SearchStrategy::*ProcessOriBestPtr)(
    ReadPair&, std::vector<std::pair<bool, std::vector<TextOcc>>>&,
    std::vector<std::pair<bool, std::vector<TextOcc>>>&,
    std::vector<std::pair<bool, std::vector<TextOcc>>>&,
    std::vector<std::pair<bool, std::vector<TextOcc>>>&,
    std::vector<PairedTextOccs>&, int, int, length_t, length_t, Counters&);

typedef void (SearchStrategy::*ProcessOriAllPtr)(
    const ReadPair& pair, std::pair<bool, std::vector<TextOcc>>& m1,
    std::pair<bool, std::vector<TextOcc>>& mRC1,
    std::pair<bool, std::vector<TextOcc>>& m2,
    std::pair<bool, std::vector<TextOcc>>& mRC2,
    std::vector<PairedTextOccs>& pairs, int maxFragSize, int minFragSize,
    length_t maxD, Counters& counters);

// pointer to filter function (hamming or edit distance)
typedef std::vector<TextOcc> (SearchStrategy::*FilterPtr)(
    Occurrences&, length_t, Counters&, const ReadBundle& bundle) const;

#ifndef RUN_LENGTH_COMPRESSION
// pointer to in-text verification function (hamming or edit distance)
typedef void (SearchStrategy::*InTextVerificationPtr)(
    const FMOcc&, length_t, length_t, const std::vector<Substring>&,
    Occurrences&, Counters&, length_t) const;
#endif

// Pointer to function that generates SAM lines for single-end
typedef void (SearchStrategy::*GenerateOutputSEPtr)(ReadBundle&, size_t, int,
                                                    std::vector<TextOcc>&,
                                                    Counters& counters) const;
// Pointer to the naive approximate matching procedure (hamming or edit)
typedef void (IndexInterface::*NaiveMatchPtr)(const std::string& pattern,
                                              length_t maxED,
                                              Counters& counters,
                                              std::vector<TextOcc>& occ);

typedef TextOcc (*CreateUnmappedSEPtr)(const ReadBundle&);

/**
 * Abstract class for searching a pattern in an index. Supports both all- and
 * best-alignment, as well as single- and paired-end mapping. Is based on search
 * schemes for different values of k
 */
class SearchStrategy {

    // ATTRIBUTES
    //-------------------------------------------------------------------------------
  protected:
    IndexInterface&
        index;        // reference to the index of the text that is searched
    std::string name; // the name of the particular search strategy
    DistanceMetric distanceMetric; // which distance metric to use

    uint32_t x = 0; // the value of x for best+X mode (0 by default)
  private:
    // helper type for best-matching in paired end alignment. The bool
    // indicates whether the distance (= index in OccVector) has been
    // processed and the second element is the approximate matches for this
    // distance
    typedef std::pair<bool, std::vector<TextOcc>> BoolAndVector;
    typedef std::vector<BoolAndVector> OccVector;

    // variables for getting info about strategy used
    PartitionStrategy partitionStrategy; // the partitioning strategy
    MappingMode mode;                    // whether to find all or best matches
    Orientation orientation; // the expected orientation of the reads

    // pointers for correct partitioning and correct distance metric
    PartitionPtr partitionPtr;   // pointer to the partition method
    StartIdxPtr startIdxPtr;     // pointer to start method (hamming or edit)
    NaiveMatchPtr naiveMatchPtr; // pointer to the correct naive approximate
                                 // matching function

    // pointers to functions for mapping mode (ALL or BEST)
    MatchPtrSE matchPtrSE; // pointer to the correct match function SE
    MatchPtrPE matchPtrPE; // pointer to the correct match function PE
    PairSEPtr pairSEPtr;   // pointer to the correct pairing function for SE
                           // alignments

    // pointers to functions for processing combinations according to
    // orientation
    ProcessOriBestPtr processOriBestPtr; // pointer to the correct function for
                                         // processing combinations in BEST mode
    ProcessOriAllPtr processOriAllPtr;   // pointer to the correct function for
                                         // processing combinations in ALL mode

    // pointers for correct filtering and in-text verification
    FilterPtr
        filterPtr; // pointer to the correct filter (edit or hamming based)
#ifndef RUN_LENGTH_COMPRESSION
    InTextVerificationPtr
        inTextVerificationPtr; // pointer to the correct in-text verification
    // (edit or hamming based)
#endif

    // pointer for generating Single End output lines (SAM or RHS)
    GenerateOutputSEPtr generateOutputSEPtr;
    CreateUnmappedSEPtr createUnmappedSEPtr;

    bool samOutput = false;
    bool generateSAMLines = true; // whether to generate SAM lines, currently
                                  // only false when mapping single-end to infer
                                  // the fragment size and orientation
    bool unmappedSAM = true; // whether to generate SAM lines for unmapped reads
    bool xaTag = false; // whether to use the XA tag for secondary alignments
    bool discordantAllowed = false; // whether discordant alignments are allowed
    length_t nDiscordantPairs = 100000; // number of discordant pairs allowed

    length_t count = 0;
    // ----------------------------------------------------------------------------
    // PROTECTED FUNCTIONS
    // ----------------------------------------------------------------------------

  protected:
    /**
     * Constructor
     * @param index the index to be used
     * @param p the partitioning strategy to be used
     * @param distanceMetric the distance metric to be used
     * @param mode whether all or best matches should be found
     * @param sm whether single or paired end mode should be used
     */
    SearchStrategy(IndexInterface& index, PartitionStrategy p,
                   DistanceMetric distanceMetric, MappingMode mode,
                   SequencingMode sm);

    /**
     * Copy constructor
     */
    SearchStrategy(const SearchStrategy& other) = default;

    // PROTECTED FUNCTIONS: PARTITIONING
    // ----------------------------------------------------------------------------

    /**
     * Calculates the number of parts for a certain max edit distance. This
     * calculation is strategy dependent
     * @param maxED the maximal allowed edit distance for the aligning
     */
    virtual uint32_t calculateNumParts(unsigned int maxED) const = 0;

    // PROTECTED FUNCTIONS: PARTITIONING (STATIC)
    // ----------------------------------------------------------------------------

    /**
     * Default function that retrieves the begin positions for  static
     * partitioning. If derived class does not implement this function then
     * uniform positions are given.
     * @param numParts  how many parts are needed
     * @param maxScore the maximal allowed score, (irrelevant for this
     * function)
     * @returns vector with doubles indicating the position (relative to the
     * length of the pattern) where a seed should be placed.
     */
    virtual const std::vector<double> getBegins(const int& numParts,
                                                const int& maxScore) const {
        std::vector<double> b;
        double u = 1.0 / numParts;
        for (int i = 1; i < numParts; i++) {
            b.push_back(i * u);
        }
        return b;
    }
    // ----------------------------------------------------------------------------

    // PROTECTED FUNCTIONS: PARTITIONING (DYNAMIC)
    // ----------------------------------------------------------------------------

    /**
     * Default function that retrieves the seeding positions for dynamic
     * partitioning. If derived class does not implement this function then
     * uniform seeds are given.
     * @param numParts how many parts are needed
     * @param maxScore the maximal allowed score, (irrelevant for this default
     * function)
     * @returns vector with doubles indicating the position (relative to the
     * length of the pattern) where a seed should be placed.
     */
    virtual const std::vector<double>
    getSeedingPositions(const int& numParts, const int& maxScore) const {

        double u = 1.0 / (numParts - 1);
        std::vector<double> s;
        for (int i = 1; i < numParts - 1; i++) {
            s.push_back(i * u);
        }
        return s;
    }

    /**
     * Helper function for dynamic partitioning. Seeds the parts.
     * @param pattern the pattern to partition
     * @param parts empty vector tro which the seeds are added
     * @param numParts how many parts are needed
     * @param maxScore the maximal allowed edit or hamming distance
     * @param exactMatchRanges a vector corresponding to the ranges for the
     * exact matches of the parts, will be cleared and filled during the
     * execution
     * @returns the number of characters used by the seeding operation
     */
    int seed(const std::string& pattern, std::vector<Substring>& parts,
             const int& numParts, const int& maxScore,
             std::vector<SARangePair>& exactMatchRanges) const;
    /**
     * Function that retrieves the weights for dynamic partitioning.
     * If derived class does not implement this function then uniform weights
     * are given with a preference for the two edges.
     * @param numParts how many parts are needed
     * @param maxScore the maximal allowed score, (irrelevant for this default
     * function)
     * @returns vector with weights
     */
    virtual const std::vector<uint64_t> getWeights(const int& numParts,
                                                   const int& maxScore) const {
        std::vector<uint64_t> w(numParts, 1);
        w.front() = 2;
        w.back() = 2;
        return w;
    }
    // ----------------------------------------------------------------------------

    // PROTECTED FUNCTIONS: ACCESS
    // ----------------------------------------------------------------------------

    /**
     * Find the maximum supported distance score of this strategy.
     */
    virtual length_t getMaxSupportedDistance() const = 0;

    /**
     * Static function that calculates the number of elements in an OccVector.
     */
    static length_t numElements(const OccVector& ov) {
        length_t num = 0;
        for (const auto& p : ov) {
            num += p.second.size();
        }
        return num;
    }

    /**
     * Static helper function for matchApproxBestPlusX. Combines the forward and
     * reverse complemented matches, removes duplicates and only keeps matches
     * that have a distance score of at most max.
     * @param ovFW the forward occurrences
     * @param ovRC the reverse complement occurrences
     * @param best the best found distance score
     * @param max the maximal allowed distance score
     * @returns a vector with all unique occurrences, sorted on distance score
     */
    static std::vector<TextOcc> combineOccVectors(OccVector& ovFW,
                                                  OccVector& ovRC,
                                                  length_t best, length_t max);
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------

    // PRIVATE FUNCTIONS
    // ----------------------------------------------------------------------------

  private:
    // PRIVATE FUNCTIONS: PARTITIONING
    // ----------------------------------------------------------------------------

    /**
     * Splits the pattern into numParts parts, either by uniform range or
     * uniform size
     * @param pattern the pattern to be split
     * @param parts the vector containing the substrings of this pattern,
     * will be cleared and filled during the execution of this method. If
     * the splitting fails for some reason, the vector will be empty
     * @param numParts  how many parts are needed
     * @param maxScore the maximum allowed edit distance
     * @param exactMatchRanges a vector corresponding to the ranges for the
     * exact matches of the parts, will be cleared and filled during the
     * execution
     * @param counters the performance counters
     */
    void partition(const std::string& pattern, std::vector<Substring>& parts,
                   const int& numParts, const int& maxScore,
                   std::vector<SARangePair>& exactMatchRanges,
                   Counters& counters) const;

    /**
     * Helper function for uniform and static partitioning. This function
     * calculates the exact match ranges of the parts
     * @param parts the parts of the pattern (they must be set before calling)
     * @param exactMatchRanges a vector corresponding to the ranges for the
     * exact matches of each part
     * @param counters the performance counters
     */
    void calculateExactMatchRanges(std::vector<Substring>& parts,
                                   std::vector<SARangePair>& exactMatchRanges,
                                   Counters& counters) const;

    // PRIVATE FUNCTIONS: PARTITIONING (UNIFORM)
    // ----------------------------------------------------------------------------

    /**
     * Splits the pattern into numParts parts, such that each part has the
     * same size
     * @param pattern the pattern to be split
     * @param parts an empty vector which will be filled with the different
     * parts
     * @param numParts, how many parts are needed
     * @param maxScore, the maximum allowed edit distance
     * @param exactMatchRanges, a vector corresponding to the ranges for the
     * exact matches of the parts, will be cleared and filled during the
     * execution
     * @param counters the performance counters
     */
    void partitionUniform(const std::string& pattern,
                          std::vector<Substring>& parts, const int& numParts,
                          const int& maxScore,
                          std::vector<SARangePair>& exactMatchRanges,
                          Counters& counters) const;

    // ----------------------------------------------------------------------------

    // PRIVATE FUNCTIONS: PARTITIONING (STATIC)
    // ----------------------------------------------------------------------------

    /**
     * Splits the pattern into numParts parts, such that each search carries
     * the same weight (on average)
     * @param pattern the pattern to be split
     * @param parts an empty vector which will be filled with the different
     * parts
     * @param numParts, how many parts are needed
     * @param exactMatchRanges, a vector corresponding to the ranges for the
     * exact matches of the parts, will be cleared and filled during the
     * execution
     * @param counters the performance counters
     */
    void partitionOptimalStatic(const std::string& pattern,
                                std::vector<Substring>& parts,
                                const int& numParts, const int& maxScore,
                                std::vector<SARangePair>& exactMatchRanges,
                                Counters& counters) const;

    /**
     * Helper function for optimal static partitioning. This function
     * creates the optimal static parts
     * @param pattern the pattern to partition
     * @param part  empty vector to which the parts are added
     * @param numParts how many parts there need to be in the partition
     */
    void setParts(const std::string& pattern, std::vector<Substring>& parts,
                  const int& numParts, const int& maxScore) const;

    // ----------------------------------------------------------------------------

    // PRIVATE FUNCTIONS: PARTITIONING (DYNAMIC)
    // ----------------------------------------------------------------------------

    /**
     * Splits the pattern into numParts parts, such that each part has
     * (approximately) the same range. The exactMatchRanges are also
     * calculated.
     * @param pattern the pattern to be split
     * @param parts an empty vector which will be filled with the different
     * parts
     * @param numParts how many parts are needed
     * @param exactMatchRanges a vector corresponding to the ranges for the
     * exact matches of the parts, will be cleared and filled during the
     * execution
     */
    void partitionDynamic(const std::string& pattern,
                          std::vector<Substring>& parts, const int& numParts,
                          const int& maxScore,
                          std::vector<SARangePair>& exactMatchRanges,
                          Counters& counters) const;

    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------

    // PRIVATE FUNCTIONS: MATCHING
    // ----------------------------------------------------------------------------
    /**
     * Executes the search recursively. If U[0] != 1, then the search will
     * start at pi[0], else the search will start with idx i and U[i]!=0 and
     * U[j]=0 with j < i. The search will not be executed if the U[0]== 0 and
     * exactMatchRanges[pi[0]].width() <= switch point of the index
     * @param s the search to follow
     * @param parts the parts of the pattern
     * @param occ Data structure containing in-index and in-text occurrences.
     * If during the search occurrences are found they will be added to this
     * data structure.
     * @param exactMatchRanges a vector corresponding to the ranges for the
     * exact matches of the parts
     * @param counters the performance counters
     */
    void doRecSearch(const Search& s, std::vector<Substring>& parts,
                     Occurrences& occ,
                     const std::vector<SARangePair>& exactMatchRanges,
                     Counters& counters) const;

    /**
     * Matches a sequence using searches with a maximal allowed distance of k
     * Found occurrences either in-text or in-index are added to the occs
     * data structure.
     * Warning: this function assumes that the index is in the correct mode
     * (revCompl or forward strand) and that the in-text verification matrix for
     * this mode has been set.
     * @param seq the sequence to match
     * @param k the maximal allowed edit distance
     * @param counters the performance counters
     * @param occs the data structure containing the occurrences.
     * @param minED the minimal  distance of all considered matches
     */
    void matchWithSearches(const std::string& seq, const length_t k,
                           Counters& counters, Occurrences& occs,
                           const length_t minED = 0);

    /**
     * Maps a read to the index using this strategy, along the correct strand.
     * @param read the read to map (in the correct orientation)
     * @param maxED the maximal allowed edit distance
     * @param counters the performance counters
     * @param pairStatus whether this is the first or second read of a pair
     * @param strand the strand of the read
     * @param minD the minimum allowed distance
     * @returns a vector with approximate occurrences in the text of this read,
     * flagged with the correct strand and pairStatus.
     */
    std::vector<TextOcc> mapRead(const std::string& read, length_t maxED,
                                 Counters& counters, PairStatus pairStatus,
                                 Strand strand, length_t minD = 0) {
        Occurrences occurrences;
        // set the index in the correct mode
        index.setIndexInMode(strand, pairStatus);
        if (maxED == 0) {
            std::vector<TextOcc> exactMatches;
            index.exactMatchesOutput(read, counters, exactMatches);
            // sort the exact matches needed for correct pairing
            std::sort(exactMatches.begin(), exactMatches.end());
            return exactMatches;
        }

        // match the sequence and return the filtered occurrences
        matchWithSearches(read, maxED, counters, occurrences, 0);
        // filter out redundant matches
        const auto& bundle = ReadBundle::createBundle(read, strand);
        auto textOcc = (this->*filterPtr)(occurrences, maxED, counters, bundle);

        // remove elements under minD
        textOcc.erase(std::remove_if(textOcc.begin(), textOcc.end(),
                                     [minD](const TextOcc& elem) {
                                         return elem.getDistance() < minD;
                                     }),
                      textOcc.end());
        return textOcc;
    }

    // PRIVATE FUNCTIONS: MATCHING (HAMMING)
    // ----------------------------------------------------------------------------

    /**
     * Starts the index with hamming distance
     * @param s the search to follow
     * @param startMatch the partial match to start from.
     * @param occ Data structure containing in-index and in-text occurrences. If
     * during the search occurrences are found they will be added to this
     * @param parts the parts of the pattern
     * @param idx the index in the search to match next
     */
    void startIndexHamming(const Search& s, const FMOcc& startMatch,
                           Occurrences& occ,
                           const std::vector<Substring>& parts,
                           Counters& counters, const int& idx) const {
        index.recApproxMatchHamming(s, startMatch, occ, parts, counters, idx);
    }

#ifndef RUN_LENGTH_COMPRESSION
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
    void inTextVerificationHamming(const FMOcc& startMatch,
                                   length_t beginInPattern, length_t maxD,
                                   const std::vector<Substring>& parts,
                                   Occurrences& occ, Counters& counters,
                                   length_t minD) const {
        index.verifyExactPartialMatchInTextHamming(
            startMatch, beginInPattern, maxD, parts, occ, counters, minD);
    }
#endif

    // ----------------------------------------------------------------------------

    // PRIVATE FUNCTIONS: MATCHING (EDIT)
    // ----------------------------------------------------------------------------

    /**
     * Starts the index with edit distance and optimized alignment for the
     * edit distance metric
     * @param s the search to follow
     * @param startMatch the partial match to start from.
     * @param occ Data structure containing in-index and in-text occurrences. If
     * during the search occurrences are found they will be added to this
     * @param parts the parts of the pattern
     * @param idx the index in the search to match next
     */
    void startIndexEdit(const Search& s, const FMOcc& startMatch,
                        Occurrences& occ, const std::vector<Substring>& parts,
                        Counters& counters, const int& idx) const {
        index.recApproxMatchEditEntry(s, startMatch, occ, parts, counters, idx);
    }

#ifndef RUN_LENGTH_COMPRESSION
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
     */
    void inTextVerificationEdit(const FMOcc& startMatch,
                                length_t beginInPattern, length_t maxD,
                                const std::vector<Substring>& parts,
                                Occurrences& occ, Counters& counters,
                                length_t minD) const {
        const auto pattern =
            Substring(parts.back(), 0, parts.back().end(), FORWARD);
        index.verifyExactPartialMatchInText(startMatch, beginInPattern, maxD,
                                            occ, counters, minD, pattern);
    }
#endif
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------

    // PRIVATE FUNCTIONS: SINGLE END
    // ----------------------------------------------------------------------------

    // PRIVATE FUNCTIONS: SINGLE END - ALL
    // ----------------------------------------------------------------------------
    /**
     * Matches a pattern approximately using this strategy in, finding all
     * approximate occurrences.
     * @param bundle the read bundle with info about the curren read
     * @param maxED the maximal allowed edit distance (or  hamming
     * distance)
     * @param counters the performance counters
     *
     * @param result a vector with approximate occurrences in the text of this
     * pattern (output)
     *
     */
    virtual void matchApproxAllMap(ReadBundle& bundle, length_t maxED,
                                   Counters& counters,
                                   std::vector<TextOcc>& result);
    // ----------------------------------------------------------------------------

    // PRIVATE FUNCTIONS: SINGLE END - BEST (+X)
    // ----------------------------------------------------------------------------

    /**
     * Helper function for matchApproxBestPlusX. Checks if the occurrences are
     * valid and assigns sequences to them if possible
     * @param occVector the vector of occurrences per edit distance
     * @param best the best edit distance found so far, can be updated
     * @param l the current edit distance
     * @param counters the performance counters
     * @param cutOff the maximal allowed edit distance
     * @param seq the sequence to which the occurrences should match, needed for
     * CIGAR string generation and trimming
     */
    void checkAlignments(OccVector& occVector, uint32_t& best, uint32_t l,
                         Counters& counters, uint32_t cutOff,
                         const std::string& seq) const;

    bool findBestAlignments(const ReadBundle& bundle, OccVector& fw,
                            OccVector& rc, Counters& counters, uint32_t x,
                            uint32_t& best, PairStatus status = FIRST_IN_PAIR);
    /**
     * Matches a pattern approximately using this strategy in, finding all
     * approximate occurrences.
     * @param readBundle the read bundle with info about the curren read
     * @param maxED the maximal allowed edit distance (or  hamming
     * distance)
     * @param counters the performance counters
     * @param minIdentity the minimal identity of all considered matches
     * @param result a vector with approximate occurrences in the text of this
     * pattern (output)
     */
    virtual void matchApproxBestPlusX(ReadBundle& readBundle, length_t x,
                                      Counters& counters,
                                      const length_t minIdentity,
                                      std::vector<TextOcc>& result);

    /**
     * Matches a pattern approximately using this strategy, finding all
     * best approximate occurrences.
     * @param readBundle the read bundle with info about the curren read
     * @param maxED the maximal allowed edit distance (or  hamming
     * distance)
     * @param counters the performance counters
     * @param minIdentity the minimal identity of all considered matches
     * @param result a vector with approximate occurrences in the text of this
     * pattern (output)
     */
    virtual void matchApproxBestMap(ReadBundle& readBundle,
                                    const length_t minIdentity,
                                    Counters& counters,
                                    std::vector<TextOcc>& result) {

        return matchApproxBestPlusX(readBundle, x, counters, minIdentity,
                                    result);
    }

    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------

    // PRIVATE FUNCTIONS: PAIRED END
    // ----------------------------------------------------------------------------

    // PRIVATE FUNCTIONS: PAIRED END - ALL
    // ----------------------------------------------------------------------------
    /**
     * @brief Finds approximate paired-end matches between two sets of
     * reads.
     *
     * This function identifies paired matches between two sets of reads,
     * considering various parameters for the matching process. The matches
     * are represented as instances of the PairedTextOccs struct and stored
     * in a vector.
     *
     * @param pair A ReadPair object containing the two bundles with info
     * @param maxEDOrIdentity The maximum allowed edit distance or min
     * identity for a match.
     * @param maxFragSize The maximum DNA fragment size.
     * @param minFragSize The minimum DNA fragment size.
     * @param counters A reference to the Counters object for tracking
     * statistics.
     * @param unpairedRead1 A vector of TextOcc instances representing
     * the occurrences that could not be paired if no valid pair was found
     * (output)
     * @return A vector of PairedTextOccs instances representing the
     * identified paired matches.
     * @see PairedTextOccs
     */
    std::vector<PairedTextOccs>
    matchApproxPairedEndAll(ReadPair& pair, length_t maxEDOrIdentity,
                            length_t maxFragSize, length_t minFragSize,
                            Counters& counters,
                            std::vector<TextOcc>& unpairedOcc);

    /**
     * Helper function for paired-end ALL matching. Processes the combination
     * of upstream sequence and downstream sequence for given fragment size
     * bounds and total allowed distance.
     * @param uSeq the upstream sequence
     * @param uStatus the status of the upstream read (first or second)
     * @param uStrand the strand of the upstream read
     * @param uVector the occurrences of the upstream read
     * (can be updated)
     * @param dSeq the downstream sequence
     * @param dStatus the status of the downstream read (first or second)
     * @param dStrand the strand of the downstream read
     * @param dVector the occurrences of the downstream read
     * (can be updated)
     * @param pairs the vector with the paired occurrences (can be updated)
     * @param maxFragSize the maximum fragment size
     * @param minFragSize the minimum fragment size
     * @param maxD the maximal distance allowed
     * @param counters the performance counters
     */
    void processComb(const std::string& uSeq, PairStatus uStatus,
                     Strand uStrand, BoolAndVector& uVector,
                     const std::string& dSeq, PairStatus dStatus,
                     Strand dStrand, BoolAndVector& dVector,
                     std::vector<PairedTextOccs>& pairs, int maxFragSize,
                     int minFragSize, length_t maxD, Counters& counters) {

        if (!uVector.first) { // first try to align the upstream sequence
            uVector = {true, mapRead(uSeq, maxD, counters, uStatus, uStrand)};
        }
        if (uVector.second.empty()) {
            return; // no possible pairs
        }
        if (!dVector.first) { // then align downstream read
            dVector = {true, mapRead(dSeq, maxD, counters, dStatus, dStrand)};
        }
        if (dVector.second.empty()) {
            return; // no possible pairs
        }
        pairOccurrences(uVector.second, dVector.second, pairs, maxFragSize,
                        minFragSize, maxD, counters, uSeq, dSeq);
    }

    /**
     * Processes the possible combinations in case of FR orientation for All
     * mapping.
     * @param pair the read pair
     * @param m1 the matches of the first read (can be updated)
     * @param mRC1 the matches of the first read revComp (can be updated)
     * @param m2 the matches of the second read (can be updated)
     * @param mRC2 the matches of the second read revComp (can be updated)
     * @param pairs the vector with the paired occurrences (can be updated)
     * @param maxFragSize the maximum fragment size
     * @param minFragSize the minimum fragment size
     * @param maxD the maximal distance allowed
     * @param counters the performance counters
     */
    void processCombFRAll(const ReadPair& pair, BoolAndVector& m1,
                          BoolAndVector& mRC1, BoolAndVector& m2,
                          BoolAndVector& mRC2,
                          std::vector<PairedTextOccs>& pairs, int maxFragSize,
                          int minFragSize, length_t maxD, Counters& counters) {
        // first 1 forward and upstream/ 2 reverse complement and downstream
        processComb(pair.getRead1(), FIRST_IN_PAIR, FORWARD_STRAND, m1,
                    pair.getRevComp2(), SECOND_IN_PAIR, REVERSE_C_STRAND, mRC2,
                    pairs, maxFragSize, minFragSize, maxD, counters);
        // then 2 forward and upstream/ 1 reverse complement and downstream
        processComb(pair.getRead2(), SECOND_IN_PAIR, FORWARD_STRAND, m2,
                    pair.getRevComp1(), FIRST_IN_PAIR, REVERSE_C_STRAND, mRC1,
                    pairs, maxFragSize, minFragSize, maxD, counters);
    }

    /**
     * Processes the possible combinations in case of FF orientation for All
     * mapping.
     * @param pair the read pair
     * @param m1 the matches of the first read (can be updated)
     * @param mRC1 the matches of the first read revComp (can be updated)
     * @param m2 the matches of the second read (can be updated)
     * @param mRC2 the matches of the second read revComp (can be updated)
     * @param pairs the vector with the paired occurrences (can be updated)
     * @param maxFragSize the maximum fragment size
     * @param minFragSize the minimum fragment size
     * @param maxD the maximal distance allowed
     * @param counters the performance counters
     */
    void processCombFFAll(const ReadPair& pair, BoolAndVector& m1,
                          BoolAndVector& mRC1, BoolAndVector& m2,
                          BoolAndVector& mRC2,
                          std::vector<PairedTextOccs>& pairs, int maxFragSize,
                          int minFragSize, length_t maxD, Counters& counters) {
        // first 1 forward and upstream/ 2 forward and downstream
        processComb(pair.getRead1(), FIRST_IN_PAIR, FORWARD_STRAND, m1,
                    pair.getRead2(), SECOND_IN_PAIR, FORWARD_STRAND, m2, pairs,
                    maxFragSize, minFragSize, maxD, counters);
        // then 2 reverse comp and upstream/ 1 reverse complement and downstream
        processComb(pair.getRevComp2(), SECOND_IN_PAIR, REVERSE_C_STRAND, mRC2,
                    pair.getRevComp1(), FIRST_IN_PAIR, REVERSE_C_STRAND, mRC1,
                    pairs, maxFragSize, minFragSize, maxD, counters);
    }

    /**
     * Processes the possible combinations in case of RF orientation for All
     * mapping.
     * @param pair the read pair
     * @param m1 the matches of the first read (can be updated)
     * @param mRC1 the matches of the first read revComp (can be updated)
     * @param m2 the matches of the second read (can be updated)
     * @param mRC2 the matches of the second read revComp (can be updated)
     * @param pairs the vector with the paired occurrences (can be updated)
     * @param maxFragSize the maximum fragment size
     * @param minFragSize the minimum fragment size
     * @param maxD the maximal distance allowed
     * @param counters the performance counters
     */
    void processCombRFAll(const ReadPair& pair, BoolAndVector& m1,
                          BoolAndVector& mRC1, BoolAndVector& m2,
                          BoolAndVector& mRC2,
                          std::vector<PairedTextOccs>& pairs, int maxFragSize,
                          int minFragSize, length_t maxD, Counters& counters) {
        // first 1 rev comp and upstream/ 2 forward and downstream
        processComb(pair.getRevComp1(), FIRST_IN_PAIR, REVERSE_C_STRAND, mRC1,
                    pair.getRead2(), SECOND_IN_PAIR, FORWARD_STRAND, m2, pairs,
                    maxFragSize, minFragSize, maxD, counters);
        // then 2 reverse comp and upstream/ 1 forward and downstream
        processComb(pair.getRevComp2(), SECOND_IN_PAIR, REVERSE_C_STRAND, mRC2,
                    pair.getRead1(), FIRST_IN_PAIR, FORWARD_STRAND, m1, pairs,
                    maxFragSize, minFragSize, maxD, counters);
    }

    // ----------------------------------------------------------------------------

    // PRIVATE FUNCTIONS: PAIRED END - BEST
    // ----------------------------------------------------------------------------
    /**
     * Process the sequence.
     * @param seq the sequence to process
     * @param status the status of the read (first or second)
     * @param strand the strand of the read
     * @param maxDist the maximum distance allowed
     * @param vector the vector with the occurrences
     * @param counters the performance counters
     * @returns true if any of the distances between 0 and maxDist have an
     * occurrence
     */
    bool processSeq(const std::string& seq, const PairStatus status,
                    const Strand strand, const length_t maxDist,
                    OccVector& vector, Counters& counters);

    /**
     * Handles trimmed occurrences (whose distance is now larger than before)
     * @param trimmedIdsSet the list with the unique ids of the trimmed
     * occurrences in vector[oDist]
     * @param oDist the original distance of the occurrences
     * @param vector the vector with the occurrences for each distance
     */
    void handleTrimmedOccs(std::vector<length_t>& trimmedIdsSet,
                           const length_t oDist, OccVector& vector);

    /**
     * Helper function for paired-end BEST matching. Processes the combination
     * of upstream sequence and downstream sequence for given fragment size
     * bounds and total allowed distance.
     * @param uSeq the upstream sequence
     * @param dSeq the downstream sequence
     * @param counters the performance counters
     * @param uStatus the status of the upstream read (first or second)
     * @param uStrand the strand of the upstream read
     * @param dStatus the status of the downstream read (first or second)
     * @param dStrand the strand of the downstream read
     * @param uVector the vector with the occurrences of the upstream read
     * (can be updated)
     * @param dVector the vector with the occurrences of the downstream read
     * (can be updated)
     * @param pairs the vector with the paired occurrences (can be updated)
     * @param maxFragSize the maximum fragment size
     * @param minFragSize the minimum fragment size
     * @param totDist the total distance allowed
     * @param minTotDist the minimum total distance allowed
     */
    void processComb(const std::string& uSeq, const std::string& dSeq,
                     Counters& counters, PairStatus uStatus, Strand uStrand,
                     PairStatus dStatus, Strand dStrand, OccVector& uVector,
                     OccVector& dVector, std::vector<PairedTextOccs>& pairs,
                     int maxFragSize, int minFragSize, length_t totDist,
                     length_t minTotDist);

    /**
     * Processes the possible combinations in case of FF orientation for BEST
     * mapping.
     * @param pair the read pair
     * @param matches1 the matches of the first read (can be updated)
     * @param matches2 the matches of the second read (can be updated)
     * @param matchesRC1 the matches of the first read revComp (can be
     * updated)
     * @param matchesRC2 the matches of the second read revComp (can be
     * updated)
     * @param pairs the vector with the paired occurrences (can be updated)
     * @param maxFragSize the maximum fragment size
     * @param minFragSize the minimum fragment size
     * @param totDist the total distance allowed
     * @param minTotDist the minimum total distance allowed
     * @param counters the performance counters
     */
    void processCombFF(ReadPair& pair, OccVector& matches1, OccVector& matches2,
                       OccVector& matchesRC1, OccVector& matchesRC2,
                       std::vector<PairedTextOccs>& pairs, int maxFragSize,
                       int minFragSize, length_t totDist, length_t minTotDist,
                       Counters& counters);

    /**
     * Processes the possible combinations in case of RF orientation for BEST
     * mapping.
     * @param pair the read pair
     * @param matches1 the matches of the first read (can be updated)
     * @param matches2 the matches of the second read (can be updated)
     * @param matchesRC1 the matches of the first read revComp (can be
     * updated)
     * @param matchesRC2 the matches of the second read revComp (can be
     * updated)
     * @param pairs the vector with the paired occurrences (can be updated)
     * @param maxFragSize the maximum fragment size
     * @param minFragSize the minimum fragment size
     * @param totDist the total distance allowed
     * @param minTotDist the minimum total distance allowed
     * @param counters the performance counters
     */
    void processCombRF(ReadPair& pair, OccVector& matches1, OccVector& matches2,
                       OccVector& matchesRC1, OccVector& matchesRC2,
                       std::vector<PairedTextOccs>& pairs, int maxFragSize,
                       int minFragSize, length_t totDist, length_t minTotDist,
                       Counters& counters);
    /**
     * Processes the possible combinations in case of FR orientation for BEST
     * mapping.
     * @param pair the read pair
     * @param matches1 the matches of the first read (can be updated)
     * @param matches2 the matches of the second read (can be updated)
     * @param matchesRC1 the matches of the first read revComp (can be
     * updated)
     * @param matchesRC2 the matches of the second read revComp (can be
     * updated)
     * @param pairs the vector with the paired occurrences (can be updated)
     * @param maxFragSize the maximum fragment size
     * @param minFragSize the minimum fragment size
     * @param totDist the total distance allowed
     * @param minTotDist the minimum total distance allowed
     * @param counters the performance counters
     */
    void processCombFR(ReadPair& pair, OccVector& matches1, OccVector& matches2,
                       OccVector& matchesRC1, OccVector& matchesRC2,
                       std::vector<PairedTextOccs>& pairs, int maxFragSize,
                       int minFragSize, length_t totDist, length_t minTotDist,
                       Counters& counters);

    /**
     * Matches the pair approximately and reports the best possible matches
     * (+ x strata). A stratum is considered to contain the pairs where the
     * total distance of that pair equals the stratum. If no concordant pair
     * is found, and discordant pairing is allowed a best discordant pair is
     * found. Again the best discordant pair is the pair with the smallest
     * total distance. If no discordant pair is found, the matches are added
     * as unpaired.
     * @param pair the read pair
     * @param minIdentity the minimal identity of all considered matches (in
     * percent, for individual reads)
     * @param maxFragSize the maximum DNA fragment size
     * @param minFragSize the minimum DNA fragment size
     * @param counters a reference to the Counters object for tracking
     * statistics
     * @param unpairedOcc a vector of TextOcc instances representing the
     * occurrences that could not be paired if no valid pair was found
     * (output)
     * @param x the number of strata above the best match
     * @param matches1 the single ended matches of the first read  if read
     * has been processed as single ended
     * @param matches2 the single ended matches of the second read if read
     * has been processed as single ended
     * @param singleEnd whether the reads have been processed as single
     * ended [default true]
     * @param read2done whether the second read has been processed [default =
     * true]
     * @returns a vector of PairedTextOccs instances representing the best
     * paired matches of the pair
     */
    std::vector<PairedTextOccs> matchApproxPairedEndBestPlusX(
        ReadPair& pair, length_t minIdentity, length_t maxFragSize,
        length_t minFragSize, Counters& counters,
        std::vector<TextOcc>& unpairedOcc, length_t x,
        std::vector<TextOcc>& matches1, std::vector<TextOcc>& matches2,
        bool singleEnd = true, bool read2done = true);

    /**
     * Matches the pair approximately and reports the best possible matches
     * (+ x strata). A stratum is considered to contain the pairs where the
     * total distance of that pair equals the stratum. If no concordant pair
     * is found, and discordant pairing is allowed a best discordant pair is
     * found. Again the best discordant pair is the pair with the smallest
     * total distance. If no discordant pair is found, the matches are added
     * as unpaired.
     * @param pair the read pair
     * @param minIdentity the minimal identity of all considered matches (in
     * percent, for individual reads)
     * @param maxFragSize the maximum DNA fragment size
     * @param minFragSize the minimum DNA fragment size

     * @param counters a reference to the Counters object for tracking
     * statistics
     * @param unpairedOcc a vector of TextOcc instances representing the
     * occurrences that could not be paired if no valid pair was found
     * (output)
     * @param x the number of strata above the best match
     * @returns a vector of PairedTextOccs instances representing the best
     * paired matches of the pair
     */
    std::vector<PairedTextOccs> matchApproxPairedEndBestPlusX(
        ReadPair& pair, length_t minIdentity, length_t maxFragSize,
        length_t minFragSize, Counters& counters,
        std::vector<TextOcc>& unpairedOcc, length_t x) {

        std::vector<TextOcc> matches1 = {}, matches2 = {};
        return matchApproxPairedEndBestPlusX(pair, minIdentity, maxFragSize,
                                             minFragSize, counters, unpairedOcc,
                                             x, matches1, matches2, false);
    }

    /**
     * Matches the pair approximately and reports all matches in the best
     * possible stratum. A stratum is considered to contain the pairs where the
     * total distance of that pair equals the stratum. If no concordant pair
     * is found, and discordant pairing is allowed a best discordant pair is
     * found. Again the best discordant pair is the pair with the smallest
     * total distance. If no discordant pair is found, the matches are added
     * as unpaired.
     * @param pair the read pair
     * @param minIdentity the minimal identity of all considered matches (in
     * percent, for individual reads)
     * @param maxFragSize the maximum DNA fragment size
     * @param minFragSize the minimum DNA fragment size
     * @param counters a reference to the Counters object for tracking
     * statistics
     * @param unpairedOcc a vector of TextOcc instances representing the
     * occurrences that could not be paired if no valid pair was found
     * (output)
     * @param x the number of strata above the best match
     * @returns a vector of PairedTextOccs instances representing the best
     * paired matches of the pair
     */
    std::vector<PairedTextOccs>
    matchApproxPairedEndBest(ReadPair& pair, length_t minIdentity,
                             length_t maxFragSize, length_t minFragSize,
                             Counters& counters,
                             std::vector<TextOcc>& unpairedOcc) {
        return matchApproxPairedEndBestPlusX(pair, minIdentity, maxFragSize,
                                             minFragSize, counters, unpairedOcc,
                                             x);
    }

    // ----------------------------------------------------------------------------

    // PRIVATE FUNCTIONS: PAIRED END - PAIRING
    // ----------------------------------------------------------------------------

    /**
     * Finds approximate paired-end matches between two sets of occurrences:
     * one set for each read in the pair. The matches are represented as
     * instances of the PairedTextOccs struct and stored in a vector.
     * @param upstreamOccs A vector of TextOcc instances representing the
     * matches of the upstream read
     * @param downStreamOccs A vector of TextOcc instances representing the
     * matches of the downstream read
     * @param pairs A vector of PairedTextOccs instances representing the
     * identified paired matches (output)
     * @param maxFragSize The maximum DNA fragment size
     * @param minFragSize The minimum DNA fragment size
     * @param maxED The maximum allowed edit distance for one read
     * @param counters A reference to the Counters object for tracking
     * statistics
     * @param uSeq The sequence of the upstream read along the correct strand
     * @param dSeq The sequence of the downstream read along the correct strand
     * @see PairedTextOccs
     */
    void pairOccurrences(std::vector<TextOcc>& upstreamOccs,
                         std::vector<TextOcc>& downStreamOccs,
                         std::vector<PairedTextOccs>& pairs,
                         length_t maxFragSize, length_t minFragSize,
                         length_t maxED, Counters& counters,
                         const std::string& uSeq,
                         const std::string& dSeq) const;

    /**
     * Finds the best approximate paired-end matches between two sets of
     * occurrences: one set for each read in the pair.
     * @param upstreamOccs The occurrences for the upstream read
     * @param downStreamOccs The occurrences for the downstream read
     * @param pairs The vector to store the best paired occurrences
     * @param maxFragSize The maximum DNA fragment size
     * @param minFragSize The minimum DNA fragment size
     * @param uMax The maximum allowed edit distance for the upstream read
     * @param dMax The maximum allowed edit distance for the downstream read
     * @param counters A reference to the Counters object for tracking
     * statistics
     * @param uTrimmedIds A set to store the IDs of the upstream occurrences
     * that might have been changed due to trimming
     * @param dTrimmedIds A set to store the IDs of the downstream occurrences
     * that might have been changed due to trimming
     * @param uSeq The sequence of the upstream read along the correct strand
     * @param dSeq The sequence of the downstream read along the correct strand
     */
    void pairOccurrencesForBestMapping(
        std::vector<TextOcc>& upstreamOccs,
        std::vector<TextOcc>& downStreamOccs,
        std::vector<PairedTextOccs>& pairs, length_t maxFragSize,
        length_t minFragSize, length_t uMax, length_t dMax, Counters& counters,
        std::set<length_t>& uTrimmedIds, std::set<length_t>& dTrimmedIds,
        const std::string& uSeq, const std::string& dSeq) const;

    /**
     * Creates unpaired matches for the current read bundle and adds them to the
     * unpairedOcc vector.
     * @param matches the matches of the read
     * @param bundle the bundle with info of the read
     * @param status the status within the pair (first or second)
     * @param maxED the maximum allowed edit distance
     * @param counters the performance counters
     * @param unpairedOcc vector with the unpaired occurrences, this will be
     * updated
     */
    void addUnpairedMatches(std::vector<TextOcc>& matches, ReadBundle& bundle,
                            PairStatus status, length_t maxED,
                            Counters& counters,
                            std::vector<TextOcc>& unpairedOcc);

    /**
     * Create the unpaired occurrences for the current read and add them to the
     * unpairedOcc vector.
     * @param forward the forward matches of the read
     * @param revComp the reverse complemented matches of the read
     * @param bundle the bundle with info of the read
     * @param status the status within the pair (first or second)
     * @param maxED the maximum allowed edit distance
     * @param counters the performance counters
     * @param unpairedOcc vector with the unpaired occurrences, this will be
     * written to
     */
    void addUnpairedMatches(std::vector<TextOcc>& forward,
                            std::vector<TextOcc>& revComp, ReadBundle& bundle,
                            PairStatus status, length_t maxED,
                            Counters& counters,
                            std::vector<TextOcc>& unpairedOcc) {
        // combine forward and revComp in a single vector without unnecessary
        // copying
        std::vector<TextOcc> allMatches;
        allMatches.reserve(forward.size() + revComp.size());

        // Move elements from forward and revComp vectors into allMatches
        allMatches.insert(allMatches.end(),
                          std::make_move_iterator(forward.begin()),
                          std::make_move_iterator(forward.end()));
        allMatches.insert(allMatches.end(),
                          std::make_move_iterator(revComp.begin()),
                          std::make_move_iterator(revComp.end()));

        // Clear forward and revComp vectors to release resources
        forward.clear(), revComp.clear();

        addUnpairedMatches(allMatches, bundle, status, maxED, counters,
                           unpairedOcc);
    }

    /**
     * Adds unpaired occurrences for both reads. This is needed if no concordant
     * pair is found and discordant pairing is not allowed or if there are too
     * many discordant candidates.
     * @param fw1 the forward matches of the first read
     * @param rc1 the reverse complemented matches of the first read
     * @param fw2 the forward matches of the second read
     * @param rc2 the reverse complemented matches of the second read
     * @param maxED the maximum allowed edit distance
     * @param counters the performance counters
     * @param unpairedOcc vector with the unpaired occurrences, this will be
     * written to
     */
    void addUnpairedMatches(std::vector<TextOcc>& fw1,
                            std::vector<TextOcc>& rc1,
                            std::vector<TextOcc>& fw2,
                            std::vector<TextOcc>& rc2, ReadPair& reads,
                            length_t maxED, Counters& counters,
                            std::vector<TextOcc>& unpairedOcc) {
        unpairedOcc.clear();
        unpairedOcc.reserve(fw1.size() + rc1.size() + fw2.size() + rc2.size());
        addUnpairedMatches(fw1, rc1, reads.getBundle1(), FIRST_IN_PAIR, maxED,
                           counters, unpairedOcc);
        addUnpairedMatches(fw2, rc2, reads.getBundle2(), SECOND_IN_PAIR, maxED,
                           counters, unpairedOcc);
    }

    /**
     * Add an unmapped record for both reads
     * @param reads the info about the reads
     * @param pairs the pairs vector, a record with a pair of unmapped
     * occurrences will be added
     */
    void addBothUnmapped(const ReadPair& reads,
                         std::vector<PairedTextOccs>& pairs) const {
        if (!unmappedSAM)
            return;
        // if both reads unmapped -> report both as unmapped
        auto m1 = TextOcc::createUnmappedSAMOccurrencePE(reads.getBundle1(),
                                                         FIRST_IN_PAIR);
        auto m2 = TextOcc::createUnmappedSAMOccurrencePE(reads.getBundle2(),
                                                         SECOND_IN_PAIR);
        PairedTextOccs pair = {std::move(m1), std::move(m2), 0};
        pairs.emplace_back(std::move(pair));
    }

    /**
     * Adds pairs to the vector where one of the reads has occurrences and
     * the other has not.
     * @param reads the info about the reads
     * @param fw1 the forward matches of the first read
     * @param rc1 the reverse complemented matches of the first read
     * @param fw2 the forward matches of the second read
     * @param rc2 the reverse complemented matches of the second read
     * @param pairs the pairs vector, records with one mapped and one
     * unmapped occurrence will be added
     * @param maxED the maximum allowed edit distance
     * @param counters the performance counters
     */
    void addOneUnmapped(ReadPair& reads, std::vector<TextOcc>& fw1,
                        std::vector<TextOcc>& rc1, std::vector<TextOcc>& fw2,
                        std::vector<TextOcc>& rc2,
                        std::vector<PairedTextOccs>& pairs, length_t maxED,
                        Counters& counters) {
        // Combine fw1 and rc1 into matches1 using move semantics
        std::vector<TextOcc> matches1;
        matches1.reserve(fw1.size() + rc1.size());
        std::move(fw1.begin(), fw1.end(), std::back_inserter(matches1));
        std::move(rc1.begin(), rc1.end(), std::back_inserter(matches1));

        // Combine fw2 and rc2 into matches2 using move semantics
        std::vector<TextOcc> matches2;
        matches2.reserve(fw2.size() + rc2.size());
        std::move(fw2.begin(), fw2.end(), std::back_inserter(matches2));
        std::move(rc2.begin(), rc2.end(), std::back_inserter(matches2));

        // Call addOneUnmapped function with matches1 and matches2
        addOneUnmapped(reads, matches1, matches2, pairs, maxED, counters);
    }

    /**
     * Adds pairs to the vector where one of the reads has occurrences and
     * the other has not.
     * @param reads the info about the reads
     * @param matches1 the matches of the first read
     * @param matches2 the matches of the second read
     * @param pairs the pairs vector, records with one mapped and one
     * unmapped occurrence will be added
     * @param maxED the maximum allowed edit distance
     * @param counters the performance counters
     */
    void addOneUnmapped(ReadPair& reads, std::vector<TextOcc>& matches1,
                        std::vector<TextOcc>& matches2,
                        std::vector<PairedTextOccs>& pairs, length_t maxED,
                        Counters& counters);
    /**
     * Adds discordant pairs to the vector
     * @param fw1 the forward matches of the first read
     * @param rc1 the reverse complemented matches of the first read
     * @param fw2 the forward matches of the second read
     * @param rc2 the reverse complemented matches of the second read
     * @param pairs the pairs vector, records with discordant pairs will be
     * added
     * @param maxED the maximum allowed edit distance
     * @param counters the performance counters
     * @param reads the info about the reads
     * @param unpairedOcc vector with the unpaired occurrences, this will be
     * written to on the off-chance one of the reads has no valid matches
     * after assigning sequences
     */
    void addDiscPairs(std::vector<TextOcc>& fw1, std::vector<TextOcc>& rc1,
                      std::vector<TextOcc>& fw2, std::vector<TextOcc>& rc2,
                      std::vector<PairedTextOccs>& pairs, length_t maxED,
                      Counters& counters, ReadPair& reads,
                      std::vector<TextOcc>& unpairedOcc);

    /**
     * Pair the occurrences of the reads discordantly. If no discordant
     * pairs could be found (or if it is not allowed) the mapped and valid
     * occurrences will be added, with SAM lines, to the unpaired vector.
     * @param fw1 the forward matches of the first read
     * @param rc1 the reverse complemented matches of the first read
     * @param fw2 the forward matches of the second read
     * @param rc2 the reverse complemented matches of the second read
     * @param pairs the pairs vector, records with discordant pairs will be
     * added here
     * @param reads the info about the reads
     * @param maxED the maximum allowed edit distance
     * @param counters the performance counters
     * @param unpairedOcc vector with the unpaired occurrences, this will be
     * written to on the off-chance one of the reads has no valid matches
     * after assigning sequences
     */
    void pairDiscordantly(BoolAndVector& fw1, BoolAndVector& rc1,
                          BoolAndVector& fw2, BoolAndVector& rc2,
                          std::vector<PairedTextOccs>& pairs, ReadPair& reads,
                          length_t maxED, Counters& counters,
                          std::vector<TextOcc>& unpairedOcc);

    /**
     * Helper function for pairDiscordantlyBest. Maps the
     * stratum if it has not been mapped yet. The function only maps the current
     * stratum.
     * @param v the vector with occurrences and a boolean flag indicating if the
     * stratum has been mapped.
     * @param maxD the maximal distance of this stratum
     * @param counters the performance counters
     * @param seq the sequence to map
     * @param status the status within the pair (first or second)
     * @param strand the strand of the the sequence
     */
    void mapStratum(OccVector& v, length_t maxD, Counters& counters,
                    const std::string& seq, PairStatus status, Strand strand) {
        if (!v[maxD].first) {

            v[maxD] = {true,
                       mapRead(seq, maxD, counters, status, strand, maxD)};
        }
    }

    /**
     * Helper function for pairDiscordantlyBest if no discordant pair could be
     * found. Finds all occurrences in the best +x strata for the given read
     * bundle.
     * @param fw the forward occurrences (can be updated). If some strata have
     * already been processed, this should be reflected in this  OccVector.
     * @param rc the reverse complemented occurrences (can be updated). If some
     * strata have already been processed, this should be reflected in this
     * OccVector.
     * @param bundle the read bundle
     * @param counters the performance counters
     * @param status the status within the pair (first or second)
     * @param x the number of strata above the best match to explore
     */
    std::vector<TextOcc> findBestMapping(OccVector& fw, OccVector& rc,
                                         const ReadBundle& bundle,
                                         Counters& counters, PairStatus status,
                                         uint32_t x);

    /**
     * Finds the best possible pairs, without regards for orientation or
     * fragment size.
     * @param fw1 the forward matches of the first read and a flag indicating if
     * the strata have been mapped
     * @param rc1 the reverse complemented matches of the first read and a flag
     * indicating if the strata have been mapped
     * @param fw2 the forward matches of the second read and a flag indicating
     * if the strata have been mapped
     * @param rc2 the reverse complemented matches of the second read and a flag
     * indicating if the strata have been mapped
     * @param pairs the pairs vector, records with the best paired occurrences
     * will be added here
     * @param reads the info about the reads
     * @param counters the performance counters
     * @param unpairedOcc vector with the unpaired occurrences, this will be
     * written to on the off-chance one of the reads has no valid matches after
     * assigning sequences
     * @param x the number of strata above the best match to explore
     */
    void pairDiscordantlyBest(OccVector& fw1, OccVector& rc1, OccVector& fw2,
                              OccVector& rc2,
                              std::vector<PairedTextOccs>& pairs,
                              ReadPair& reads, Counters& counters,
                              std::vector<TextOcc>& unpairedOcc, length_t x);

    // ----------------------------------------------------------------------------

    // PRIVATE FUNCTIONS: PAIRED END - PROCESSING SINGLE ENDED OCCS FROM
    // INFERENCE
    // ----------------------------------------------------------------------------
    /**
     * Pair single-ended matches. This function assumes that the occurrences
     * have a CIGAR string and an assigned sequence. The matches in the
     * second vector will also be flagged as matches for the second read.
     * @param matches1 the matches of the first read
     * @param matches2 the matches of the second read
     * @param maxFragSize the maximum DNA fragment size
     * @param minFragSize the minimum DNA fragment size
     * @param counters a reference to the Counters object for tracking
     * statistics
     * @param readPair the read pair with info about the current read
     * @param maxED the maximum allowed edit distance for one read
     * @param unpairedOcc a vector of TextOcc instances representing the
     * occurrences that could not be paired if no valid pair was found
     * (output)
     * @param read2done whether the second read has been processed
     * @returns a vector of PairedTextOccs instances representing the
     * identified paired matches
     * @see PairedTextOccs
     */
    std::vector<PairedTextOccs> pairSingleEndedMatchesAll(
        std::vector<TextOcc>& matches1, std::vector<TextOcc>& matches2,
        length_t maxFragSize, length_t minFragSize, Counters& counters,
        ReadPair& readPair, length_t maxED, std::vector<TextOcc>& unpairedOcc,
        bool read2done);

    /**
     * Pair single-ended matches in BEST mapping mode. This function assumes
     * that the occurrences have a CIGAR string and an assigned sequence.
     * The matches in the second vector will also be flagged as matches for the
     * second read.
     * @param matches1 the matches of the first read
     * @param matches2 the matches of the second read
     * @param maxFragSize the maximum DNA fragment size
     * @param minFragSize the minimum DNA fragment size
     * @param counters a reference to the Counters object for tracking
     * statistics
     * @param readPair the read pair with info about the current pair
     * @param minIdentity the minimal identity of all considered matches
     * @param read2done whether the second read has been processed
     */
    std::vector<PairedTextOccs> pairSingleEndedMatchesBest(
        std::vector<TextOcc>& matches1, std::vector<TextOcc>& matches2,
        length_t maxFragSize, length_t minFragSize, Counters& counters,
        ReadPair& readPair, length_t minIdentity,
        std::vector<TextOcc>& unpairedOcc, bool read2done) {
        // x = 0, singleEnd = true
        return matchApproxPairedEndBestPlusX(
            readPair, minIdentity, maxFragSize, minFragSize, counters,
            unpairedOcc, 0, matches1, matches2, true, read2done);
    }

    /**
     * Static function to process the single ended matches in a best-map
     * scenario. This will add the occurrences to the correct OccVector for
     * further pairing. After the function the elements of matches1 and matches2
     * will be in an undefined state.
     * @param matchesSE1 the matches of the first read
     * @param matchesSE2 the matches of the second read
     * @param fw1 the forward matches of the first read (can be updated)
     * @param rc1 the reverse complemented matches of the first read (can be
     * updated)
     * @param fw2 the forward matches of the second read (can be updated)
     * @param rc2 the reverse complemented matches of the second read (can be
     * updated)
     * @param read2done whether the second read has been processed
     */
    static void addSingleEndedForBest(std::vector<TextOcc>& matchesSE1,
                                      std::vector<TextOcc>& matchesSE2,
                                      OccVector& fw1, OccVector& rc1,
                                      OccVector& fw2, OccVector& rc2,
                                      bool read2done);

    // ----------------------------------------------------------------------------

    // PRIVATE FUNCTIONS: POST PROCESSING
    // ----------------------------------------------------------------------------
    /**
     * Filter the occurrences (both in-index and in-text) based on hamming
     * distance. The in-index occurrences are first converted to in-text
     * occurrences. The function also ensures that each in-text occurrence has a
     * CIGAR string.
     * @param occ the occurrences
     * @param maxED the maximal allowed edit distance
     * @param counters performance counters
     * @param bundle the read bundle with info about the current read (dummy
     * because filter for edit distance needs it).
     * @returns the filtered in-text occurrences with CIGAR strings
     */
    std::vector<TextOcc> filterHamming(Occurrences& occ, length_t maxED,
                                       Counters& counters,
                                       const ReadBundle& bundle) const {
        return index.getTextOccHamming(occ, counters);
    }

    /**
     * Filter the occurrences (both in-index and in-text) based on edit
     * distance. The in-index occurrences are first converted to in-text
     * occurrences.
     * @param occ the occurrences
     * @param maxED the maximal allowed edit distance
     * @param counters performance counters
     * @param bundle the read bundle with info about the current read
     * @returns the filtered in-text occurrences ()
     */
    std::vector<TextOcc>
    filterEditWithCIGARCalculation(Occurrences& occ, length_t maxED,
                                   Counters& counters,
                                   const ReadBundle& bundle) const {
        auto occs =
            index.getUniqueTextOccurrences(occ, maxED, counters, bundle);
        index.generateCIGARS(occs, counters, bundle);
        return occs;
    }

    /**
     * Filter the occurrences (both in-index and in-text) based on edit
     * distance. The in-index occurrences are first converted to in-text
     * occurrences. The function also ensures that each in-text occurrence has a
     * CIGAR string.
     * @param occ the occurrences
     * @param maxED the maximal allowed edit distance
     * @param counters performance counters
     * @param bundle the read bundle with info about the current read
     * @returns the filtered in-text occurrences ()
     */
    std::vector<TextOcc>
    filterEditWithoutCIGARCalculation(Occurrences& occ, length_t maxED,
                                      Counters& counters,
                                      const ReadBundle& bundle) const {
        return index.getUniqueTextOccurrences(occ, maxED, counters, bundle);
    }

    /**
     * Assigns a sequence (and the CIGAR string if not RLC) to an occurrence.
     * Returns how and if the sequence was found.
     * @param occ the occurrence to assign a sequence to
     * @param counters the performance counters
     * @param maxED the maximal allowed edit distance
     * @param pattern the pattern that was matched
     * @returns how and if the sequence was found
     */
    SeqNameFound assignSequenceAndCIGAR(TextOcc& occ, Counters& counters,
                                        length_t maxED,
                                        const std::string& pattern) const {
        // first handle CIGAR string

        if (!occ.hasCigar()) {
            index.generateCIGAR(occ, counters, pattern);
        }

        return assignSequence(occ, counters, maxED, pattern);
    }

    /**
     * Assigns a reference sequence to an occurrence if it has not been
     * assigned. The range of the occurrence is also changed to show its
     * position within the assigned sequence
     * @param occ the occurrence to assign a sequence to
     * @param counters the performance counters
     * @param maxED the maximal allowed edit distance
     * @param pattern the pattern that was matched
     * @returns true if a sequence was assigned, false if no sequence could
     * be assigned
     */
    SeqNameFound assignSequence(TextOcc& occ, Counters& counters,
                                length_t maxED,
                                const std::string& pattern) const {

        if (!occ.isSeqNameChecked()) {
            index.setIndexInMode(occ.getStrand(), occ.getPairStatus());

            length_t seqID;
            // if a sequence name is found the coordinates of the occurrence
            // will change
            auto found = index.findSeqName(occ, seqID, counters, maxED,
                                           distanceMetric, pattern);
            if (found != SeqNameFound::NOT_FOUND) {
                occ.setAssignedSequence(found, seqID);
            } else {
                occ.setAssignedSequenceNotFound();
            }
        }

        return occ.getSeqNameFound();
    }

    TextOcc createUnmappedRecordSE(ReadBundle& bundle) const {
        return createUnmappedSEPtr(bundle);
    }

    /**
     * Generates the SAM line for the first occurrence and put the other
     * occurrences in the XA tag.
     * @param bundle the read bundle with info about the current read
     * @param nHits the number of hits with the best score
     * @param minScore the best score of the occurrences
     * @param occs all the occurrences in the text
     * @param counters the performance counters
     */
    void generateSE_SAM_XATag(ReadBundle& bundle, size_t nHits, int minScore,
                              std::vector<TextOcc>& occs,
                              Counters& counters) const {

        occs.front().generateSAMSingleEndXA(bundle, nHits, minScore,
                                            occs.begin() + 1, occs.end(),
                                            index.getSeqNames());
        counters.inc(Counters::DROPPED_UNIQUE_MATCHES, occs.size() - 1);
        occs.resize(1);
    }

    /**
     * Generates the SAM lines for all occurrences. The first occurrence is the
     * primary occurrence and the rest are secondary.
     * @param bundle the read bundle with info about the current read
     * @param nHits the number of hits with the best score
     * @param minScore the best score of the occurrences
     * @param occs all the occurrences in the text
     * @param counters the performance counters
     */
    void generateSE_SAM(ReadBundle& bundle, size_t nHits, int minScore,
                        std::vector<TextOcc>& occs, Counters& counters) const {
        occs.front().generateSAMSingleEndFirst(bundle, nHits, minScore,
                                               index.getSeqNames());
        for (length_t i = 1; i < occs.size(); i++) {
            occs[i].generateSAMSingleEndNotFirst(bundle.getSeqID(), nHits,
                                                 minScore, index.getSeqNames());
        }
    }

    /**
     * Generates the RHS line for the occurrence.
     * @param bundle the read bundle with info about the current read
     * @param nHits the number of hits with the best score (dummy parameter)
     * @param minScore the best score of the occurrences (dummy parameter)
     * @param occs all the occurrences in the text
     * @param counters the performance counters
     */
    void generateSE_RHS(ReadBundle& bundle, size_t nHits, int minScore,
                        std::vector<TextOcc>& occs, Counters& counters) const {
        // INFO: parameter bundle cannot be const because it needs same
        // signature as generateSE_SAM, which might update the revQual in the
        // bundle

        // for each distinct assigned sequence + distance combination keep only
        // the first ones

        // Sort the vector using a lambda comparator
        std::sort(occs.begin(), occs.end(),
                  [](const TextOcc& lhs, const TextOcc& rhs) {
                      if (lhs.getDistance() != rhs.getDistance())
                          return lhs.getDistance() < rhs.getDistance();
                      return lhs.getAssignedSequenceID() <
                             rhs.getAssignedSequenceID();
                  });

        // Use std::unique with a lambda to remove consecutive duplicates
        auto uniqueEnd =
            std::unique(occs.begin(), occs.end(),
                        [](const TextOcc& lhs, const TextOcc& rhs) {
                            return lhs.getDistance() == rhs.getDistance() &&
                                   lhs.getAssignedSequenceID() ==
                                       rhs.getAssignedSequenceID();
                        });

        // Resize the vector to remove the redundant elements
        occs.erase(uniqueEnd, occs.end());

        occs.front().generateRHSSingleEnd(bundle, occs.begin() + 1, occs.end(),
                                          index.getSeqNames());
        counters.inc(Counters::DROPPED_UNIQUE_MATCHES, occs.size() - 1);
        occs.resize(1);
    }

    /**
     * Generates the SAM lines for the occurrences. This function will assign
     * sequences to the occurrences sometimes with the use of trimming (if
     * possible).
     * @param occs the occurrences in the text
     * @param bundle the read bundle with info about the curren read
     * @param counters the performance counters
     * @param cutOff the maximum allowed distance for an occurrence
     */
    void generateOutputSingleEnd(std::vector<TextOcc>& occs, ReadBundle& bundle,
                                 Counters& counters, length_t cutOff) const;
    /**
     * Generates the SAM lines for paired matches. This function assumes
     * that the occurrences have a CIGAR string and an assigned sequence.
     * @param pairedMatches the paired matches. The SAM line of each match
     * will be set.
     * @param readPair the read pair with info about the current read
     */

    void generateSAMPairedEnd(std::vector<PairedTextOccs>& pairedMatches,
                              ReadPair& readPair);

    /**
     * Generates the RHS lines for the occurrences. This function will assign
     * sequences to the occurrences sometimes with the use of trimming (if
     * possible).
     * @param occs the occurrences in the text
     * @param bundle the read bundle with info about the curren read
     * @param counters the performance counters
     * @param cutOff the maximum allowed distance for an occurrence
     */
    void generateRHSSingleEnd(std::vector<TextOcc>& occs, ReadBundle& bundle,
                              Counters& counters, length_t cutOff) const;

    //-------------------------------------------------------------------------
    //-------------------------------------------------------------------------

  public:
    // PUBLIC FUNCTIONS
    // ----------------------------------------------------------------------------

    virtual ~SearchStrategy() {
    }

    // PUBLIC FUNCTIONS: ACCESS
    // ----------------------------------------------------------------------------

    /**
     * Retrieves the name of this strategy, derived classes should set a
     * meaningful name
     */
    std::string getName() const {
        return name;
    }
    /**
     * Retrieve the partitioning strategy in string format
     */
    std::string getPartitioningStrategy() const;

    /**
     * Retrieve the distance metric in string format
     */
    std::string getDistanceMetric() const;

    /**
     * Retrieve the mapping in string format
     */
    std::string getMappingModeString() const;

    /**
     * Retrieve the mapping mode
     */
    MappingMode getMappingMode() const {
        return mode;
    }

    /**
     * Check whether this occurrence is in the reference text that originates
     * from the first FASTA file
     * @param occ the occurrence to check
     * @returns true if the occurrence is in the first file
     */
    inline bool isInFirstFile(const TextOcc& occ) const {
        return index.isInFirstFile(occ);
    }

#ifndef RUN_LENGTH_COMPRESSION

    /**
     * Retrieves the text of the index (for debugging purposes)
     */
    std::string getText() const {
        return index.getText();
    }
#endif

    /**
     * Get the tipping point for in-text verification
     */
    length_t getSwitchPoint() const {
        return index.getSwitchPoint();
    }

    /**
     * Calculate the maximal allowed distance for a given minimal identity
     * and sequence size
     * @param minIdentity the minimal identity
     * @param seqSize the size of the sequence
     * @returns the maximal allowed edit distance
     */
    length_t getMaxED(length_t minIdentity, length_t seqSize) const {
        assert(minIdentity >= 50 && minIdentity <= 100);

        length_t cutOff = (seqSize * (100 - minIdentity)) / 100;
        const length_t distanceCutoff = BEST_CUTOFF_COLUMBA;

        return std::min(std::min(distanceCutoff, getMaxSupportedDistance()),
                        cutOff);
    }

    /**
     * Turn on or off the generation of SAM lines
     */
    void setGenerateSAMLines(bool generateSAMLines) {
        this->generateSAMLines = generateSAMLines;
    }

    /**
     * Turn on or off the allowing of discordant pairs
     */
    void setDiscordantAllowed(bool discordantAllowed) {
        this->discordantAllowed = discordantAllowed;
    }

    /**
     * Set the number of allowed discordant pairs
     */
    void setNumDiscordantAllowed(length_t numDiscordantAllowed) {
        this->nDiscordantPairs = numDiscordantAllowed;
    }

    /**
     * Set whether the output should be SAM format or the custom Read Hit
     * Summary Format
     */
    void setSamOutput(bool samOutput) {
        this->samOutput = samOutput;
        if (samOutput) {
            generateOutputSEPtr = &SearchStrategy::generateSE_SAM;
            createUnmappedSEPtr = &TextOcc::createUnmappedSAMOccurrenceSE;
        } else {
            generateOutputSEPtr = &SearchStrategy::generateSE_RHS;
            createUnmappedSEPtr = &TextOcc::createUnmappedRHSOccurrenceSE;
        }
    }

    /**
     * Set the value of x for best+x mapping
     * @param strataAfterBest the number of strata to search after the best
     */
    void setStrataAfterBest(uint32_t strataAfterBest) {
        this->x = strataAfterBest;
    }

    /**
     * Set the mapping mode to use
     * @param mode the mapping mode to use
     */
    void setMappingMode(MappingMode mode) {
        this->mode = mode;
        if (mode == ALL) {
            matchPtrSE = &SearchStrategy::matchApproxAllMap;
            matchPtrPE = &SearchStrategy::matchApproxPairedEndAll;
            pairSEPtr = &SearchStrategy::pairSingleEndedMatchesAll;
        } else {
            // BEST
            matchPtrSE = &SearchStrategy::matchApproxBestMap;
            matchPtrPE = &SearchStrategy::matchApproxPairedEndBest;
            pairSEPtr = &SearchStrategy::pairSingleEndedMatchesBest;
        }
    }

    /**
     * Set whether to use the XA tag for secondary alignments or use separate
     * SAM lines for them.
     * @param xaTag whether to use the XA tag for secondary alignments
     */
    void setXATag(bool xaTag) {
        this->xaTag = xaTag;
        if (xaTag) {
            generateOutputSEPtr = &SearchStrategy::generateSE_SAM_XATag;
        } else {
            generateOutputSEPtr = &SearchStrategy::generateSE_SAM;
        }
    }

    /**
     * Set whether unmapped SAM records should not be generated
     * @param unmappedSam whether unmapped SAM records should not be generated
     */
    void setUnmappedSam(bool unmappedSam) {
        this->unmappedSAM = unmappedSam;
    }

    /**
     * Sets the orientation to use for paired-end alignment
     * @param orientation the orientation of the paired-end reads (FF, FR,
     * RF)
     */
    void setOrientation(Orientation orientation) {
        this->orientation = orientation;
        if (orientation == FF) {
            processOriBestPtr = &SearchStrategy::processCombFF;
            processOriAllPtr = &SearchStrategy::processCombFFAll;
        } else if (orientation == FR) {
            processOriBestPtr = &SearchStrategy::processCombFR;
            processOriAllPtr = &SearchStrategy::processCombFRAll;
        } else {
            processOriBestPtr = &SearchStrategy::processCombRF;
            processOriAllPtr = &SearchStrategy::processCombRFAll;
        }
    }

    /**
     * Whether the search scheme supports this distance score
     * @param maxScore the maximal score
     * @returns true if the search scheme supports this distance score
     */
    virtual bool supportsDistanceScore(const int& maxScore) const = 0;

    /**
     * Whether the search scheme supports the best mapping mode
     * @param max the maximal score for best mapping mode (output)
     * @returns true if the search scheme supports the best mapping mode
     */
    bool supportsBestMapping(int& max) const {
        max = 0;
        for (length_t i = 0; i < MAX_K; i++) {
            if (supportsDistanceScore(i)) {
                return (max == 0) ? false : true;
            }
            max = i;
        }
        return (max == 0) ? false : true;
    }
    // ----------------------------------------------------------------------------

    // PUBLIC FUNCTIONS: MATCHING
    // ----------------------------------------------------------------------------

    /**
     * Creates all searches for this specific strategy. This is strategy
     * dependent
     * @param maxED the maximal allowed edit distance for the aligning
     * @param exactMatchRanges the ranges for the exact matches of the parts
     */
    virtual const std::vector<Search>&
    createSearches(unsigned int maxED,
                   const std::vector<SARangePair>& exactMatchRanges) const = 0;

    /**
     * Matches a pattern approximately using this strategy in the current
     * mode. WARNING: make sure that the orientation is set correctly before
     * calling. (@see setOrientation)
     * @param readBundle the read bundle with info about the curren read
     * @param maxEDOrIdentity the maximal allowed edit distance (or  hamming
     * distance) in case of all-map mode and the minimal identity in case of
     * best-map mode
     * @param counters the performance counters
     * @param result a vector with approximate occurrences in the text of this
     * pattern (output)
     *
     */
    virtual void matchApprox(ReadBundle& readBundle, length_t maxEDOrIdentity,
                             Counters& counters, std::vector<TextOcc>& result) {
        (this->*matchPtrSE)(readBundle, maxEDOrIdentity, counters, result);
    }

    /**
     * Matches a pair of reads approximately using this strategy in the current
     * mode. WARNING: make sure that the orientation is set correctly before
     * calling.
     * @param pair the read pair with info about the reads
     * @param maxEDOrIdentity the maximal allowed distance or minimal required
     * identity
     * @param maxFragSize the maximum fragment size
     * @param minFragSize the minimum fragment size
     * @param counters the performance counters
     * @param unpairedOccurrences the unpaired occurrences in case no valid pair
     * could be found
     * @returns a vector with approximate paired occurrences
     */
    virtual std::vector<PairedTextOccs>
    matchApproxPE(ReadPair& pair, length_t maxEDOrIdentity,
                  length_t maxFragSize, length_t minFragSize,
                  Counters& counters,
                  std::vector<TextOcc>& unpairedOccurrences) {
        return (this->*matchPtrPE)(pair, maxEDOrIdentity, maxFragSize,
                                   minFragSize, counters, unpairedOccurrences);
    }
    // ----------------------------------------------------------------------------

    // PUBLIC FUNCTIONS: PAIRING
    // ----------------------------------------------------------------------------
    /**
     * Pairs the single ended matches in the current mapping mode (ALL
     * or BEST) and with the current metric
     * @param readPair the read pair with info about the reads
     * @param matches1 the matches of the first read
     * @param matches2 the matches of the second read
     * @param maxFragSize the maximum fragment size
     * @param minFragSize the minimum fragment size
     * @param maxEDOrMinIdentity the maximal allowed distance or minimal
     * identity
     * @param unpairedOccurrences the unpaired occurrences in case no
     * valid pair could be found
     * @param counters the performance counters
     * @param read2done whether the second read has been processed
     * @returns a vector with approximate occurrences as pairs
     *
     */
    virtual std::vector<PairedTextOccs>
    pairSingleEnd(ReadPair& readPair, std::vector<TextOcc>& matches1,
                  std::vector<TextOcc>& matches2, length_t maxFragSize,
                  length_t minFragSize, length_t maxEDOrMinIdentity,
                  std::vector<TextOcc>& unpairedOccurrences, Counters& counters,
                  bool read2done) {
        return (this->*pairSEPtr)(matches1, matches2, maxFragSize, minFragSize,
                                  counters, readPair, maxEDOrMinIdentity,
                                  unpairedOccurrences, read2done);
    }
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
};

// ============================================================================
// CLASS CUSTOM SEARCHSTRATEGY
// ============================================================================

class CustomSearchStrategy;

// Pointer to the correct getBegins() function for static partitioning
typedef const std::vector<double> (CustomSearchStrategy::*GetBeginsPtr)(
    const int& numParts, const int& maxScore) const;

// Pointer to the correct getSeedingPositions() function for dynamic
// partitioning
typedef const std::vector<double> (
    CustomSearchStrategy::*GetSeedingPositionsPtr)(const int& numParts,
                                                   const int& maxScore) const;
// Pointer to the correct getWeights() function for dynamic partitioning
typedef const std::vector<uint64_t> (CustomSearchStrategy::*GetWeightsPtr)(
    const int& numParts, const int& maxScore) const;

/**
 * This is a derived class of SearchStrategy. It creates a custom scheme
 * using files provided by the user. It takes a specified folder in which a
 * file "name.txt", containing the name of the scheme on the first line, and
 * for each supported distance score a subfolder exists. Such a subfolder
 * has as name the distance score. Each subfolder must contain at least a
 * file "searches.txt". The different searches of the scheme for this
 * distance score should be written on separate lines of this file. Each
 * search consists out of three arrays, pi, L and U, the arrays are
 * separated by a single space. Each array is written between curly braces
 * {} and the different values are separated. The pi array must be
 * zero-based.
 *
 * The different sub-folders can also contain files for static and dynamic
 * partitioning. The file "static_partitioning.txt" should consist out of
 * one line of space- separated percentages (between 0 and 1 - exclusive).
 * These percentages point to start positions of the second to the last part
 * of the partitioned pattern (relative to the size of the pattern). Hence,
 * if a search scheme partitions a pattern in k parts, then k+1 percentages
 * should be provided. The file "dynamic_partitioning.txt" should consist
 * out of two lines. The first line contains k-1 space- separated
 * percentages. These percentages are the seeding positions of the middle
 * parts (the first and last part are seeded at the begin and end and thus
 * do not need a percentage). Note that this line can be empty if the
 * pattern is partitioned into 2 parts. The second line contains k integers,
 * where k is the number of parts. Each integer corresponds to the weight
 * given to that part in the dynamic partitioning process.
 */
class CustomSearchStrategy : public SearchStrategy {
  private:
    std::vector<std::vector<Search>>
        schemePerED; // the search schemes for each distance score,
    std::array<bool, MAX_K>
        supportsMaxScore; // if a particular distance score is supported

    // helpers for static partitioning
    std::vector<std::vector<double>>
        staticPositions;                     // the static positions per k
    std::vector<GetBeginsPtr> beginsPointer; // pointer to the correct
                                             // getBegins() function,
                                             // either default or custom

    // helpers for dynamic partitioning
    std::vector<std::vector<double>>
        seedingPositions; // the seeds for dynamic partitioning per

    std::vector<std::vector<uint64_t>>
        weights; // the weights for dynamic partitioning per score

    std::vector<GetSeedingPositionsPtr>
        seedingPointer; // pointer to the correct
                        // getSeedingPositions() function,
                        // either default or custom

    std::vector<GetWeightsPtr>
        weightsPointers; // pointer to the correct getWeights() function,
                         // either default or custom

    /**
     * Retrieves the search scheme from a folder, also checks if the scheme
     * is valid
     * @param pathToFolder the path to the folder containing the search
     * scheme
     * @param verbose if the sanity check should be verbose
     */
    void getSearchSchemeFromFolder(const std::string& pathToFolder,
                                   bool verbose);

    /**
     * If the values provided for dynamic partitioning for the given max
     * score are valid (i.e. strictly increasing and between 0 and 1). Will
     * throw a runtime error if this is not the case
     * @param maxScore, the score to check
     */
    void sanityCheckDynamicPartitioning(const int& maxScore) const;

    /**
     * If the values provided for static partitioning for the given max
     * score are valid (i.e. strictly increasing and between 0 and 1). Will
     * throw a runtime error if this is not the case
     * @param maxScore, the score to check
     */
    void sanityCheckStaticPartitioning(const int& maxScore) const;

    /**
     * Parse the search from a line.
     * @param line the line to parse
     * @param idx the index of the search
     * @returns the parsed line as a search, if the line is not valid a
     * runtime error will be thrown.
     */
    Search makeSearch(const std::string& line, length_t idx) const;

    /**
     * Checks whether the connectivity property is satisfied for all
     * searches. Will throw a runtime error if one of these is not satisfied
     * @param verbose if the information about which search covers which
     * error distribution should be written to standard out
     */
    void sanityCheck(bool verbose) const;

    // static partitioning
    /**
     * Gets the static positions in the default manner.
     * @param numParts the number of parts of the pattern
     * @param maxScore the maximal allowed score
     */
    const std::vector<double> getBeginsDefault(const int& numParts,
                                               const int& maxScore) const {
        return SearchStrategy::getBegins(numParts, maxScore);
    }

    /**
     * Gets the static positions in the custom manner (i.e. with values
     * provided by user in "static_partitioning.txt").
     * @param numParts the number of parts of the pattern
     * @param maxScore the maximal allowed score
     */
    const std::vector<double> getBeginsCustom(const int& numParts,
                                              const int& maxScore) const {
        return staticPositions[maxScore - 1];
    }
    /**
     * Overridden function of the base class. Retrieves the begin positions
     * in the custom way if values were provided in a
     * "static_partitioning.txt" file, otherwise the base class function
     * will be called;
     */
    const std::vector<double> getBegins(const int& numParts,
                                        const int& maxScore) const override {
        assert(supportsMaxScore[maxScore - 1]);
        return (this->*beginsPointer[maxScore - 1])(numParts, maxScore);
    }

    // dynamic partitioning
    /**
     * Gets the seeding positions in the default manner.
     * @param numParts the number of parts of the pattern
     * @param maxScore the maximal allowed score
     */
    const std::vector<double>
    getSeedingPositionsDefault(const int& numParts, const int& maxScore) const {
        return SearchStrategy::getSeedingPositions(numParts, maxScore);
    }

    /**
     * Gets the seeding positions in the custom manner (i.e. with values
     * provided by user in "dynamic_partitioning.txt").
     * @param numParts the number of parts of the pattern
     * @param maxScore the maximal allowed score
     */
    const std::vector<double>
    getSeedingPositionsCustom(const int& numParts, const int& maxScore) const {
        return seedingPositions[maxScore - 1];
    }

    /**
     * Overridden function of the base class. Retrieves the seeding
     * positions in the custom way if values were provided in a
     * "dynamic_partitioning.txt" file, otherwise the base class function
     * will be called;
     */
    const std::vector<double>
    getSeedingPositions(const int& numParts,
                        const int& maxScore) const override {
        assert(supportsMaxScore[maxScore - 1]);
        return (this->*seedingPointer[maxScore - 1])(numParts, maxScore);
    }

    /**
     * Gets the seeding positions in the default manner.
     * @param numParts the number of parts of the pattern
     * @param maxScore the maximal allowed score
     */
    const std::vector<uint64_t> getWeightsDefault(const int& numParts,
                                                  const int& maxScore) const {
        return SearchStrategy::getWeights(numParts, maxScore);
    }

    /**
     * Gets the weights in the custom manner (i.e. with values
     * provided by user in "dynamic_partitioning.txt").
     * @param numParts the number of parts of the pattern
     * @param maxScore the maximal allowed score
     */
    const std::vector<uint64_t> getWeightsCustom(const int& numParts,
                                                 const int& maxScore) const {
        return weights[maxScore - 1];
    }
    /**
     * Overridden function of the base class. Retrieves the weights
     * positions in the custom way if values were provided in a
     * "dynamic_partitioning.txt" file, otherwise the base class function
     * will be called;
     */
    const std::vector<uint64_t> getWeights(const int& numParts,
                                           const int& maxScore) const override {
        assert(supportsMaxScore[maxScore - 1]);
        return (this->*weightsPointers[maxScore - 1])(numParts, maxScore);
    }

  public:
    CustomSearchStrategy(IndexInterface& index, const std::string& pathToFolder,
                         PartitionStrategy p, DistanceMetric metric,
                         MappingMode mode, SequencingMode sMode,
                         bool verbose = false)
        : SearchStrategy(index, p, metric, mode, sMode) {

        // resize and fill the vectors
        schemePerED.resize(MAX_K);
        staticPositions.resize(MAX_K);
        beginsPointer.resize(MAX_K, &CustomSearchStrategy::getBeginsDefault);
        seedingPositions.resize(MAX_K);
        weights.resize(MAX_K);
        seedingPointer.resize(
            MAX_K, &CustomSearchStrategy::getSeedingPositionsDefault);
        weightsPointers.resize(MAX_K, &CustomSearchStrategy::getWeightsDefault);

        // read the schemes
        getSearchSchemeFromFolder(pathToFolder, verbose);
    }

    uint32_t calculateNumParts(unsigned int maxED) const override {
        assert(supportsMaxScore[maxED - 1]);
        return schemePerED[maxED - 1][0].getNumParts();
    }
    const std::vector<Search>&
    createSearches(unsigned int maxED,
                   const std::vector<SARangePair>& ranges = {}) const override {
        assert(supportsMaxScore[maxED - 1]);
        return schemePerED[maxED - 1];
    }
    length_t getMaxSupportedDistance() const override {
        // find biggest value for which supportsMaxScore[0: value -1] is
        // true
        length_t max = 0;
        for (length_t i = 0; i < MAX_K; i++) {
            if (supportsMaxScore[i]) {
                max++;
            } else {
                break;
            }
        }
        return max;
    }

    bool supportsDistanceScore(const int& maxScore) const override {
        return supportsMaxScore[maxScore - 1];
    }
};

// ============================================================================
// CLASS MULTIPLE SCHEMES
// ============================================================================

/**
 * @class MultipleSchemes
 *
 * @brief A class representing a collection of search schemes for
 * approximate string matching.
 *
 * This class manages multiple search schemes. It calculates the number of parts
 * and provides a method to choose the search scheme to use based on the number
 * of exact matches of the parts.
 */
class MultipleSchemes {
  private:
    unsigned int k; // the maximum number of errors for these schemes
    std::vector<SearchScheme> schemes; // the list of schemes
    unsigned int numParts = 0; // the number of parts (p) for this scheme

    /**
     * @brief Private helper function to read search schemes from a folder.
     *
     * This function reads search schemes from text files in a specified
     * folder. Each file should be named "scheme1.txt," "scheme2.txt," and
     * so on.
     *
     * @param pathToFolder The path to the folder containing the scheme
     * files.
     * @param verbose A flag indicating whether to display verbose output
     * during reading.
     */
    void getSchemesFromFolder(const std::string& pathToFolder, bool verbose) {
        // in this folder there should be several files:
        // scheme1.txt, scheme2.txt, scheme3.txt, ...

        std::string base = "scheme";
        int x = 1; // Starting value of x

        while (true) {
            std::string filename = base + std::to_string(x) + ".txt";
            std::string filePath = pathToFolder + PATH_SEPARATOR + filename;

            std::ifstream file(filePath);
            if (!file.is_open()) {
                // File does not exist, break the loop
                break;
            }

            // Read in the scheme
            schemes.emplace_back(SearchScheme::readScheme(file, filename, k));

            // Increment x for the next iteration
            ++x;
        }

        if (!schemes.empty()) {

            numParts = schemes.front().getNumParts();

            if (std::any_of(schemes.begin(), schemes.end(),
                            [this](const SearchScheme& scheme) {
                                return scheme.getNumParts() != numParts;
                            })) {
                throw std::runtime_error(
                    "Not all schemes have same number of parts in: " +
                    pathToFolder);
            }
        }
    }

  public:
    /**
     * @brief Constructor for the MultipleSchemes class.
     *
     * Initializes the MultipleSchemes object by reading search schemes from
     * a folder.
     *
     * @param pathToFolder The path to the folder containing the scheme
     * files.
     * @param k The maximum number of errors for these schemes.
     * @param verbose A flag indicating whether to display verbose output
     * during reading.
     */
    MultipleSchemes(const std::string& pathToFolder, unsigned int k,
                    bool verbose = false)
        : k(k) {
        getSchemesFromFolder(pathToFolder, verbose);
    }

    /**
     * @brief Constructor for an MultipleSchemes instance.
     *
     * Initializes the MultipleSchemes object by creating an empty list of
     * schemes.
     *
     * @param k The maximum number of errors for these schemes.
     * @param verbose A flag indicating whether to display verbose output
     * during reading.
     */
    MultipleSchemes(unsigned int k, bool verbose = false) : k(k) {
    }

    /**
     * Create the multiple schemes from a single strategy by mirroring the
     * pi-strings
     * @param strategy the strategy to take the scheme from
     * @param k the maximum number of errors
     */
    MultipleSchemes(const SearchStrategy& strategy, unsigned int k) : k(k) {
        std::vector<SARangePair> dummyRanges;
        const auto& searches = strategy.createSearches(k, dummyRanges);
        SearchScheme scheme(searches, k);
        numParts = scheme.getNumParts();
        SearchScheme reversedScheme = scheme.mirrorPiStrings();
        schemes.push_back(scheme);
        schemes.push_back(reversedScheme);
    }

    /**
     * @brief Calculate the number of parts for a given maximum edit
     * distance. This function should only be called with maxED == k.
     *
     * @param maxED The maximum edit distance for which to calculate the
     * number of parts.
     * @return The number of parts (p) corresponding to the maximum edit
     * distance.
     */
    uint32_t calculateNumParts(unsigned int maxED) const {
        assert(maxED == k);
        return numParts;
    }

    /**
     * @brief Create searches based on a maximum edit distance and specified
     * ranges (=number of exact matches of the searches).
     *
     * This method creates searches based on a given maximum edit distance
     * and a set of SARangePair ranges. It iterates over the schemes and
     * chooses the scheme for which the starting part of the critical search
     * has the least amount of exact matches.
     *
     * @param maxED The maximum edit distance for creating searches.
     * @param ranges A vector of SARangePair objects representing search
     * ranges.
     * @return A vector of Search objects based on the specified parameters.
     */
    const std::vector<Search>&
    createSearches(unsigned int maxED,
                   const std::vector<SARangePair>& ranges) const {
        assert(maxED == k);
        assert(!isEmpty());
        // find the scheme for which the critical search has the smallest
        // exact range to start from

        // if the total number of exact matches of the parts is smaller than
        // the number of parts dynamic selection has too much overhead
        unsigned int total =
            std::accumulate(ranges.begin(), ranges.end(), 0,
                            [](unsigned int sum, const SARangePair& range) {
                                return sum + range.width();
                            });

        if (total <= numParts) {
            return schemes[0].getSearches();
        }

        int minIndex = 0;
        unsigned int minValue =
            ranges[schemes[0].getCriticalPartIndex()].width();

        for (unsigned int i = 1; i < schemes.size(); ++i) {
            unsigned int criticalPartIndex = schemes[i].getCriticalPartIndex();
            if (ranges[criticalPartIndex].width() < minValue) {
                minValue = ranges[criticalPartIndex].width();
                minIndex = i;
            }
        }
        return schemes[minIndex].getSearches();
    }

    /**
     * @brief find out whether the list of schemes is empty
     *
     * @returns true if the list of schemes is empty.
     */
    bool isEmpty() const {
        return schemes.empty();
    }

    /**
     * Add a single scheme to the list of schemes. The scheme to add must be
     * defined for the same number of maximal errors and the same number of
     * parts.
     * @param scheme the scheme to add
     */
    void addScheme(const SearchScheme& scheme) {
        if (scheme.getNumParts() != numParts) {
            throw std::runtime_error(
                "The number of parts of the scheme does not match the number "
                "of parts of the other schemes");
        }
        if (scheme.getK() != k) {
            throw std::runtime_error(
                "The number of errors of the scheme does not match the number "
                "of errors of the other schemes");
        }
        schemes.push_back(scheme);
    }
};
// ============================================================================
// CLASS MULTIPLE SCHEME STRATEGY
// ============================================================================

/**
 * @class MultipleSchemesStrategy
 *
 * @brief A class representing a strategy for handling multiple search
 * schemes with different edit distances.
 *
 * This class is designed to manage multiple search schemes with varying
 * edit distances. It inherits from the SearchStrategy class and provides
 * methods for calculating the number of parts and creating searches based
 * on a specified edit distance and search ranges. The Multiple Schemes will
 * dynamically choose what scheme to use.
 */
class MultipleSchemesStrategy : public SearchStrategy {
  private:
    std::vector<MultipleSchemes>
        schemesPerED; // A vector of MultipleSchemes objects, one for each
                      // edit distance.

    static bool directoryExists(const std::string& path) {
#ifdef _WIN32
        DWORD ftyp = GetFileAttributesA(path.c_str());
        if (ftyp == INVALID_FILE_ATTRIBUTES)
            return false; // something is wrong with your path!

        if (ftyp & FILE_ATTRIBUTE_DIRECTORY)
            return true; // this is a directory!

        return false; // this is not a directory!
#else
        struct stat info;
        return (stat(path.c_str(), &info) == 0 && (info.st_mode & S_IFDIR));
#endif
    }

    static bool fileExists(const std::string& path) {
        std::ifstream file(path);
        return file.good();
    }

    /**
     * @brief Private helper function to read and initialize multiple search
     * schemes.
     *
     * This function reads and initializes multiple search schemes from a
     * specified folder for various edit distances. If an edit distance is
     * not provided the function stops reading.
     *
     * @param pathToFolder The path to the folder containing the search
     * schemes.
     * @param verbose A flag indicating whether to display verbose output
     * during reading.
     */
    void readSchemes(const std::string& pathToFolder, bool verbose) {

        // get the name of the scheme
        std::string line;
        {
            std::ifstream ifs(pathToFolder + "name.txt");
            if (!ifs) {
                throw std::runtime_error(
                    "Problem reading: " + pathToFolder +
                    "name.txt\nDid you provide a directory to "
                    "a search scheme without a name file?");
            }
            std::getline(ifs, line);
            name = line;
            ifs.close();
        }

        // read in the schemes
        for (uint16_t k = 1; k <= MAX_K; k++) {
            std::string dir = pathToFolder + PATH_SEPARATOR + std::to_string(k);

            if (directoryExists(pathToFolder)) {
                std::string file = dir + PATH_SEPARATOR + "scheme1.txt";
                if (fileExists(file)) {
                    schemesPerED.emplace_back(MultipleSchemes(dir, k, verbose));
                    continue;
                }
            }
            // Add an empty search scheme
            schemesPerED.emplace_back(MultipleSchemes(k, verbose));
        }
    }

  public:
    /**
     * @brief Constructor for the MultipleSchemesStrategy class.
     *
     * Initializes the MultipleSchemesStrategy object by reading and
     * initializing search schemes for various edit distances.
     *
     * @param index The IndexInterface object used for searching.
     * @param pathToFolder The path to the folder containing the search
     * schemes.
     * @param p The partition strategy to use
     * @param metric The distance metric to use
     * @param mode The mapping mode to use
     * @param sMode The sequencing mode to use
     * @param verbose A flag indicating whether to display verbose output
     * during reading (default: false).
     */
    MultipleSchemesStrategy(IndexInterface& index,
                            const std::string& pathToFolder,
                            PartitionStrategy p, DistanceMetric metric,
                            MappingMode mode, SequencingMode sMode,
                            bool verbose = false)
        : SearchStrategy(index, p, metric, mode, sMode) {

        readSchemes(pathToFolder, verbose);
    }

    /**
     * Construct a multiple schemes strategy from a single strategy. This will
     * take the schemes in the strategy and mirror them and add both versions to
     * the multiple schemes.
     * @param strategy the strategy to take the scheme from
     * @param verbose if the constructor should be verbose
     */
    MultipleSchemesStrategy(SearchStrategy* strategy, bool verbose = false)
        : SearchStrategy(*strategy) {
        schemesPerED.reserve(MAX_K);
        for (uint16_t k = 1; k <= MAX_K; k++) {
            if (strategy->supportsDistanceScore(k)) {
                // create the multiple schemes
                schemesPerED.emplace_back(MultipleSchemes(*strategy, k));
            } else {
                // add an empty scheme
                schemesPerED.emplace_back(MultipleSchemes(k, verbose));
            }
        }
    }

    void addScheme(const SearchScheme& scheme) {
        unsigned int k = scheme.getK();
        if (k > MAX_K) {
            throw std::runtime_error("The maximum number of errors is " +
                                     std::to_string(MAX_K));
        }
        if (schemesPerED.size() < k) {
            for (uint16_t i = schemesPerED.size(); i < k; i++) {
                schemesPerED.emplace_back(MultipleSchemes(i + 1));
            }
        }
        schemesPerED[k - 1].addScheme(scheme);
    }

    /**
     * @brief Calculate the number of parts for a given maximum edit
     * distance.
     *
     * @param maxED The maximum edit distance for which to calculate the
     * number of parts.
     * @return The number of parts (p) corresponding to the specified
     * maximum edit distance.
     */
    uint32_t calculateNumParts(unsigned int maxED) const override {
        assert(maxED - 1 < schemesPerED.size());
        assert(!schemesPerED[maxED - 1].isEmpty());
        return schemesPerED[maxED - 1].calculateNumParts(maxED);
    }

    length_t getMaxSupportedDistance() const override {
        return schemesPerED.size();
    }
    bool supportsDistanceScore(const int& maxScore) const override {
        return (length_t)maxScore - 1 < schemesPerED.size() &&
               !schemesPerED[maxScore - 1].isEmpty();
    }

    /**
     * @brief Create searches based on a maximum edit distance and specified
     * ranges.
     *
     * This method creates searches based on a given maximum edit distance
     * and a set of SARangePair ranges.
     *
     * @param maxED The maximum edit distance for creating searches.
     * @param ranges A vector of SARangePair objects representing search
     * ranges.
     * @return A vector of Search objects based on the specified parameters.
     */
    const std::vector<Search>&
    createSearches(unsigned int maxED,
                   const std::vector<SARangePair>& ranges) const override {
        assert(maxED - 1 < schemesPerED.size());
        assert(!schemesPerED[maxED - 1].isEmpty());
        return schemesPerED[maxED - 1].createSearches(maxED, ranges);
    }
};

// ============================================================================
// CLASS NaiveBackTrackingStrategy
// ============================================================================

/**
 * Search strategy that uses naive backtracking for approximate string matching.
 */
class NaiveBackTrackingStrategy : public SearchStrategy {
  private:
    const std::vector<std::vector<Search>> searches = createSearchesVector();
    uint32_t calculateNumParts(unsigned int maxED) const override {
        return 1;
    }
    const std::vector<Search>&
    createSearches(unsigned int maxED,
                   const std::vector<SARangePair>& ranges) const override {

        return searches[maxED - 1];
    }

    length_t getMaxSupportedDistance() const override {
        return MAX_K;
    }
    bool supportsDistanceScore(const int& maxScore) const override {
        return maxScore <= (int)getMaxSupportedDistance();
    }

    static std::vector<std::vector<Search>> createSearchesVector() {
        std::vector<std::vector<Search>> searches;
        searches.reserve(MAX_K - 1);
        for (length_t k = 1; k < MAX_K; k++)
            searches.push_back({Search::makeSearch({0}, {0}, {k}, 0)});
        return searches;
    }

  public:
    NaiveBackTrackingStrategy(IndexInterface& index, PartitionStrategy p,
                              DistanceMetric metric, MappingMode mode,
                              SequencingMode sMode)
        : SearchStrategy(index, p, metric, mode, sMode) {
        name = "Naive backtracking";
    }
};

// ============================================================================
// HARDCODED CUSTOM CLASSES
// ============================================================================

/**
 * Hardcoded class for the Kucherov  k+1 strategy.
 */
class KucherovKPlus1 : public SearchStrategy {
  private:
    const std::vector<Search> ED1 = {
        Search::makeSearch({0, 1}, {0, 1}, {0, 1}, 0),
        Search::makeSearch({1, 0}, {0, 0}, {0, 1}, 1)};
    const std::vector<Search> ED2 = {
        Search::makeSearch({0, 1, 2}, {0, 0, 0}, {0, 2, 2}, 0),
        Search::makeSearch({2, 1, 0}, {0, 0, 0}, {0, 1, 2}, 1),
        Search::makeSearch({1, 0, 2}, {0, 0, 1}, {0, 1, 2}, 2)};

    const std::vector<Search> ED3 = {
        Search::makeSearch({0, 1, 2, 3}, {0, 0, 0, 0}, {0, 1, 3, 3}, 0),
        Search::makeSearch({1, 0, 2, 3}, {0, 0, 1, 1}, {0, 1, 3, 3}, 1),
        Search::makeSearch({2, 3, 1, 0}, {0, 0, 0, 0}, {0, 1, 3, 3}, 2),
        Search::makeSearch({3, 2, 1, 0}, {0, 0, 1, 1}, {0, 1, 3, 3}, 3)};

    const std::vector<Search> ED4 = {
        Search::makeSearch({0, 1, 2, 3, 4}, {0, 0, 0, 0, 0}, {0, 2, 2, 4, 4},
                           0),
        Search::makeSearch({4, 3, 2, 1, 0}, {0, 0, 0, 0, 0}, {0, 1, 3, 4, 4},
                           1),
        Search::makeSearch({1, 0, 2, 3, 4}, {0, 0, 1, 3, 3}, {0, 1, 3, 3, 4},
                           2),
        Search::makeSearch({0, 1, 2, 3, 4}, {0, 0, 1, 3, 3}, {0, 1, 3, 3, 4},
                           3),
        Search::makeSearch({3, 2, 4, 1, 0}, {0, 0, 0, 1, 1}, {0, 1, 2, 4, 4},
                           4),
        Search::makeSearch({2, 1, 0, 3, 4}, {0, 0, 0, 1, 3}, {0, 1, 2, 4, 4},
                           5),
        Search::makeSearch({1, 0, 2, 3, 4}, {0, 0, 1, 2, 4}, {0, 1, 2, 4, 4},
                           6),
        Search::makeSearch({0, 1, 2, 3, 4}, {0, 0, 0, 3, 4}, {0, 0, 4, 4, 4},
                           7)};

    const std::vector<std::vector<Search>> schemePerED = {ED1, ED2, ED3, ED4};

    const std::vector<std::vector<double>> seedingPositions = {
        {}, {0.57}, {0.38, 0.65}, {0.38, 0.55, 0.73}};

    const std::vector<std::vector<uint64_t>> weights = {
        {1, 1}, {39, 10, 40}, {400, 4, 5, 400}, {100, 5, 1, 6, 105}};

    const std::vector<std::vector<double>> staticPositions = {
        {0.5}, {0.41, 0.7}, {0.25, 0.50, 0.75}, {0.27, 0.47, 0.62, 0.81}};
    uint32_t calculateNumParts(unsigned int maxED) const override {
        return maxED + 1;
    }
    const std::vector<Search>&
    createSearches(unsigned int maxED,
                   const std::vector<SARangePair>& ranges) const override {
        assert(maxED >= 1);
        assert(maxED <= 4);

        return schemePerED[maxED - 1];
    }
    const std::vector<double> getBegins(const int& numParts,
                                        const int& maxScore) const override {
        return staticPositions[maxScore - 1];
    }
    const std::vector<uint64_t> getWeights(const int& numParts,
                                           const int& maxScore) const override {
        return weights[maxScore - 1];
    }

    const std::vector<double>
    getSeedingPositions(const int& numParts,
                        const int& maxScore) const override {
        return seedingPositions[maxScore - 1];
    }

    length_t getMaxSupportedDistance() const override {
        return 4;
    }

  public:
    KucherovKPlus1(IndexInterface& index, PartitionStrategy p,
                   DistanceMetric metric, MappingMode mode,
                   SequencingMode sMode)
        : SearchStrategy(index, p, metric, mode, sMode) {
        name = "KUCHEROV K + 1";
    };
    bool supportsDistanceScore(const int& maxScore) const override {
        return maxScore >= 1 && maxScore <= 4;
    }
};

/**
 * Hardcoded class for the Kucherov  k+2 strategy.
 */
class KucherovKPlus2 : public SearchStrategy {
  private:
    const std::vector<Search> ED1 = {
        Search::makeSearch({0, 1, 2}, {0, 0, 0}, {0, 1, 1}, 0),
        Search::makeSearch({1, 2, 0}, {0, 0, 0}, {0, 0, 1}, 1)};
    const std::vector<Search> ED2 = {
        Search::makeSearch({0, 1, 2, 3}, {0, 0, 0, 0}, {0, 1, 1, 2}, 0),
        Search::makeSearch({3, 2, 1, 0}, {0, 0, 0, 0}, {0, 1, 2, 2}, 1),
        Search::makeSearch({1, 2, 3, 0}, {0, 0, 0, 1}, {0, 0, 1, 2}, 2),
        Search::makeSearch({0, 1, 2, 3}, {0, 0, 0, 2}, {0, 0, 2, 2}, 3)};

    const std::vector<Search> ED3 = {
        Search::makeSearch({0, 1, 2, 3, 4}, {0, 0, 0, 0, 0}, {0, 1, 2, 3, 3},
                           0),
        Search::makeSearch({1, 2, 3, 4, 0}, {0, 0, 0, 0, 0}, {0, 1, 2, 2, 3},
                           1),
        Search::makeSearch({2, 3, 4, 1, 0}, {0, 0, 0, 0, 1}, {0, 1, 1, 3, 3},
                           2),
        Search::makeSearch({3, 4, 2, 1, 0}, {0, 0, 0, 1, 2}, {0, 0, 3, 3, 3},
                           3)};

    const std::vector<Search> ED4 = {
        Search::makeSearch({0, 1, 2, 3, 4, 5}, {0, 0, 0, 0, 0, 0},
                           {0, 1, 2, 3, 4, 4}, 0),
        Search::makeSearch({1, 2, 3, 4, 5, 0}, {0, 0, 0, 0, 0, 0},
                           {0, 1, 2, 3, 4, 4}, 1),
        Search::makeSearch({5, 4, 3, 2, 1, 0}, {0, 0, 0, 0, 0, 1},
                           {0, 1, 2, 2, 4, 4}, 2),
        Search::makeSearch({3, 4, 5, 2, 1, 0}, {0, 0, 0, 0, 1, 2},
                           {0, 1, 1, 3, 4, 4}, 3),
        Search::makeSearch({2, 3, 4, 5, 1, 0}, {0, 0, 0, 0, 2, 3},
                           {0, 1, 1, 2, 4, 4}, 4),
        Search::makeSearch({4, 5, 3, 2, 1, 0}, {0, 0, 0, 1, 3, 3},
                           {0, 0, 3, 3, 4, 4}, 5),
        Search::makeSearch({0, 1, 2, 3, 4, 5}, {0, 0, 0, 3, 3, 3},
                           {0, 0, 3, 3, 4, 4}, 6),
        Search::makeSearch({0, 1, 2, 3, 4, 5}, {0, 0, 0, 0, 4, 4},
                           {0, 0, 2, 4, 4, 4}, 7),
        Search::makeSearch({2, 3, 1, 0, 4, 5}, {0, 0, 0, 1, 2, 4},
                           {0, 0, 2, 2, 4, 4}, 8),
        Search::makeSearch({4, 5, 3, 2, 1, 0}, {0, 0, 0, 0, 4, 4},
                           {0, 0, 1, 4, 4, 4}, 9)};

    const std::vector<std::vector<Search>> schemePerED = {ED1, ED2, ED3, ED4};

    uint32_t calculateNumParts(unsigned int maxED) const override {
        return maxED + 2;
    }
    const std::vector<Search>&
    createSearches(unsigned int maxED,
                   const std::vector<SARangePair>& ranges) const override {
        return schemePerED[maxED - 1];
    }

    const std::vector<std::vector<double>> seedingPositions = {
        {0.94}, {0.48, 0.55}, {0.4, 0.63, 0.9}, {0.34, 0.5, 0.65, 0.7}};

    const std::vector<std::vector<uint64_t>> weights = {
        {11, 10, 1},
        {400, 4, 1, 800},
        {6, 3, 2, 1, 1},
        {52, 42, 16, 14, 1, 800}};

    const std::vector<std::vector<double>> staticPositions = {
        {0.47, 0.94},
        {0.35, 0.50, 0.65},
        {0.22, 0.44, 0.66, 0.88},
        {0.18, 0.37, 0.53, 0.69, 0.83}};

    const std::vector<double> getBegins(const int& numParts,
                                        const int& maxScore) const override {
        return staticPositions[maxScore - 1];
    }
    const std::vector<uint64_t> getWeights(const int& numParts,
                                           const int& maxScore) const override {
        return weights[maxScore - 1];
    }

    virtual const std::vector<double>
    getSeedingPositions(const int& numParts,
                        const int& maxScore) const override {
        return seedingPositions[maxScore - 1];
    }

    length_t getMaxSupportedDistance() const override {
        return 4;
    }

  public:
    KucherovKPlus2(IndexInterface& index, PartitionStrategy p,
                   DistanceMetric metric, MappingMode mode,
                   SequencingMode sMode)
        : SearchStrategy(index, p, metric, mode, sMode) {
        name = "KUCHEROV K + 2";
    };

    bool supportsDistanceScore(const int& maxScore) const override {
        return maxScore >= 1 && maxScore <= 4;
    }
};

/**
 * Hardcoded class for the optimal strategy by Kianfar et al.
 */
class OptimalKianfar : public SearchStrategy {
  private:
    const std::vector<Search> ED1 = {
        Search::makeSearch({0, 1}, {0, 0}, {0, 1}, 0),
        Search::makeSearch({1, 0}, {0, 1}, {0, 1}, 1)};
    const std::vector<Search> ED2 = {
        Search::makeSearch({0, 1, 2}, {0, 0, 2}, {0, 1, 2}, 0),
        Search::makeSearch({2, 1, 0}, {0, 0, 0}, {0, 2, 2}, 1),
        Search::makeSearch({1, 2, 0}, {0, 1, 1}, {0, 1, 2}, 2)};

    const std::vector<Search> ED3 = {
        Search::makeSearch({0, 1, 2, 3}, {0, 0, 0, 3}, {0, 2, 3, 3}, 0),
        Search::makeSearch({1, 2, 3, 0}, {0, 0, 0, 0}, {1, 2, 3, 3}, 1),
        Search::makeSearch({2, 3, 1, 0}, {0, 0, 2, 2}, {0, 0, 3, 3}, 2)};

    const std::vector<Search> ED4 = {
        Search::makeSearch({0, 1, 2, 3, 4}, {0, 0, 0, 0, 4}, {0, 3, 3, 4, 4},
                           0),
        Search::makeSearch({1, 2, 3, 4, 0}, {0, 0, 0, 0, 0}, {2, 2, 3, 3, 4},
                           1),
        Search::makeSearch({4, 3, 2, 1, 0}, {0, 0, 0, 3, 3}, {0, 0, 4, 4, 4},
                           2)};

    const std::vector<std::vector<Search>> schemePerED = {ED1, ED2, ED3, ED4};
    const std::vector<std::vector<double>> seedingPositions = {
        {}, {0.50}, {0.34, 0.66}, {0.42, 0.56, 0.67}};

    const std::vector<std::vector<uint64_t>> weights = {
        {1, 1}, {10, 1, 5}, {1, 1, 1, 1}, {7, 2, 1, 3, 5}};

    const std::vector<std::vector<double>> staticPositions = {
        {0.5}, {0.30, 0.60}, {0.17, 0.69, 0.96}, {0.2, 0.5, 0.6, 0.8}};
    uint32_t calculateNumParts(unsigned int maxED) const override {
        return maxED + 1;
    }
    const std::vector<Search>&
    createSearches(unsigned int maxED,
                   const std::vector<SARangePair>& ranges) const override {
        if (maxED < 1 || maxED > 5) {
            throw std::invalid_argument("max ED should be between 1 and 4");
        }
        return schemePerED[maxED - 1];
    }

    const std::vector<double> getBegins(const int& numParts,
                                        const int& maxScore) const override {
        return staticPositions[maxScore - 1];
    }
    const std::vector<uint64_t> getWeights(const int& numParts,
                                           const int& maxScore) const override {
        return weights[maxScore - 1];
    }

    const std::vector<double>
    getSeedingPositions(const int& numParts,
                        const int& maxScore) const override {
        return seedingPositions[maxScore - 1];
    }

    length_t getMaxSupportedDistance() const override {
        return 4;
    }

  public:
    OptimalKianfar(IndexInterface& index, PartitionStrategy p,
                   DistanceMetric metric, MappingMode m, SequencingMode sMode)
        : SearchStrategy(index, p, metric, m, sMode) {
        name = "OPTIMAL KIANFAR";
    };

    bool supportsDistanceScore(const int& maxScore) const override {
        return maxScore >= 1 && maxScore <= 4;
    }
};

// ============================================================================
// CLASS O1StarSearchStrategy
// ============================================================================

/**
 * A concrete derived class of SearchStrategy.The strategy here is founded
 * on this observation: if x errors are allowed and the pattern is divided
 * up in (x + 2) parts then every match with max x errors contains a seed
 * consisting of n parts, where the first and last part of the seed contain no
 * errors and all parts between these contain exactly one error.(2 <= n <= x +
 * 2) */
class O1StarSearchStrategy : public SearchStrategy {
  private:
    const std::vector<Search> ED1 = {
        Search::makeSearch({0, 1, 2}, {0, 0, 0}, {0, 1, 1}, 0),
        Search::makeSearch({1, 2, 0}, {0, 0, 0}, {0, 0, 1}, 1)};
    const std::vector<Search> ED2 = {
        Search::makeSearch({0, 1, 2, 3}, {0, 0, 0, 0}, {0, 1, 2, 2}, 0),
        Search::makeSearch({1, 2, 3, 0}, {0, 0, 0, 0}, {0, 1, 2, 2}, 1),
        Search::makeSearch({2, 3, 1, 0}, {0, 0, 0, 0}, {0, 0, 2, 2}, 2)};

    const std::vector<Search> ED3 = {
        Search::makeSearch({0, 1, 2, 3, 4}, {0, 0, 0, 0, 0}, {0, 1, 3, 3, 3},
                           0),
        Search::makeSearch({1, 2, 3, 4, 0}, {0, 0, 0, 0, 0}, {0, 1, 3, 3, 3},
                           1),
        Search::makeSearch({2, 3, 4, 1, 0}, {0, 0, 0, 0, 0}, {0, 1, 3, 3, 3},
                           2),
        Search::makeSearch({3, 4, 2, 1, 0}, {0, 0, 0, 0, 0}, {0, 0, 3, 3, 3},
                           3)};

    const std::vector<Search> ED4 = {
        Search::makeSearch({0, 1, 2, 3, 4, 5}, {0, 0, 0, 0, 0, 0},
                           {0, 1, 4, 4, 4, 4}, 0),
        Search::makeSearch({1, 2, 3, 4, 5, 0}, {0, 0, 0, 0, 0, 0},
                           {0, 1, 4, 4, 4, 4}, 1),
        Search::makeSearch({2, 3, 4, 5, 1, 0}, {0, 0, 0, 0, 0, 0},
                           {0, 1, 4, 4, 4, 4}, 2),
        Search::makeSearch({3, 4, 5, 2, 1, 0}, {0, 0, 0, 0, 0, 0},
                           {0, 1, 4, 4, 4, 4}, 3),
        Search::makeSearch({4, 5, 3, 2, 1, 0}, {0, 0, 0, 0, 0, 0},
                           {0, 0, 4, 4, 4, 4}, 4),
    };

    const std::vector<std::vector<Search>> schemePerED = {ED1, ED2, ED3, ED4};

    uint32_t calculateNumParts(unsigned int maxED) const override {
        return maxED + 2;
    }
    const std::vector<Search>&
    createSearches(unsigned int maxED,
                   const std::vector<SARangePair>& ranges) const override {
        return schemePerED[maxED - 1];
    }

    const std::vector<std::vector<double>> seedingPositions = {
        {0.94}, {0.51, 0.93}, {0.34, 0.64, 0.88}, {0.28, 0.48, 0.63, 0.94}};

    const std::vector<std::vector<uint64_t>> weights = {
        {11, 10, 1}, {20, 11, 11, 10}, {3, 2, 2, 1, 1}, {1, 2, 2, 1, 2, 1}};

    const std::vector<std::vector<double>> staticPositions = {
        {0.50, 0.96},
        {0.26, 0.64, 0.83},
        {0.22, 0.46, 0.67, 0.95},
        {0.19, 0.37, 0.57, 0.74, 0.96}};

    const std::vector<double> getBegins(const int& numParts,
                                        const int& maxScore) const override {
        return staticPositions[maxScore - 1];
    }
    const std::vector<uint64_t> getWeights(const int& numParts,
                                           const int& maxScore) const override {
        return weights[maxScore - 1];
    }

    virtual const std::vector<double>
    getSeedingPositions(const int& numParts,
                        const int& maxScore) const override {
        return seedingPositions[maxScore - 1];
    }

    length_t getMaxSupportedDistance() const override {
        return 4;
    }

  public:
    O1StarSearchStrategy(IndexInterface& index, PartitionStrategy p,
                         DistanceMetric metric, MappingMode mode,
                         SequencingMode sMode)
        : SearchStrategy(index, p, metric, mode, sMode) {
        name = "01*0";
    };

    bool supportsDistanceScore(const int& maxScore) const override {
        return maxScore >= 1 && maxScore <= 4;
    }
};

// ============================================================================
// CLASS PIGEONHOLE SEARCHSTRATEGY
// ============================================================================

/**
 * A concrete derived class of SearchStrategy. The strategy here is founded
 * on this observation: if x errors are allowed and the pattern is divided
 * up in (x
 * + 1) sections then every approximate match has an exact match with at
 * least one of the sections. The strategy iterates over the sections, it
 * tries to exactly match the current section, then approximately match the
 * pattern before this section and after the the pattern after this section
 * with the remaining edit distance.
 */
class PigeonHoleSearchStrategy : public SearchStrategy {
  private:
    const std::vector<Search> ED1 = {
        Search::makeSearch({0, 1}, {0, 0}, {0, 1}, 0),
        Search::makeSearch({1, 0}, {0, 0}, {0, 1}, 1)};
    const std::vector<Search> ED2 = {
        Search::makeSearch({0, 1, 2}, {0, 0, 0}, {0, 2, 2}, 0),
        Search::makeSearch({1, 2, 0}, {0, 0, 0}, {0, 2, 2}, 1),
        Search::makeSearch({2, 1, 0}, {0, 0, 0}, {0, 2, 2}, 2)};

    const std::vector<Search> ED3 = {
        Search::makeSearch({0, 1, 2, 3}, {0, 0, 0, 0}, {0, 3, 3, 3}, 0),
        Search::makeSearch({1, 0, 2, 3}, {0, 0, 0, 0}, {0, 3, 3, 3}, 1),
        Search::makeSearch({2, 3, 1, 0}, {0, 0, 0, 0}, {0, 3, 3, 3}, 2),
        Search::makeSearch({3, 2, 1, 0}, {0, 0, 0, 0}, {0, 3, 3, 3}, 3)};

    const std::vector<Search> ED4 = {
        Search::makeSearch({0, 1, 2, 3, 4}, {0, 0, 0, 0, 0}, {0, 4, 4, 4, 4},
                           0),
        Search::makeSearch({1, 2, 3, 4, 0}, {0, 0, 0, 0, 0}, {0, 4, 4, 4, 4},
                           1),
        Search::makeSearch({2, 3, 4, 1, 0}, {0, 0, 0, 0, 0}, {0, 4, 4, 4, 4},
                           2),
        Search::makeSearch({3, 4, 2, 1, 0}, {0, 0, 0, 0, 0}, {0, 4, 4, 4, 4},
                           3),
        Search::makeSearch({4, 3, 2, 1, 0}, {0, 0, 0, 0, 0}, {0, 4, 4, 4, 4},
                           4)};

    const std::vector<std::vector<Search>> schemePerED = {ED1, ED2, ED3, ED4};
    uint32_t calculateNumParts(unsigned int maxED) const override {
        return maxED + 1;
    }
    const std::vector<Search>&
    createSearches(unsigned int maxED,
                   const std::vector<SARangePair>& ranges) const override {
        return schemePerED[maxED - 1];
    }

    length_t getMaxSupportedDistance() const override {
        return 4;
    }

  public:
    PigeonHoleSearchStrategy(IndexInterface& index, PartitionStrategy p,
                             DistanceMetric metric, MappingMode mode,
                             SequencingMode sMode)
        : SearchStrategy(index, p, metric, mode, sMode) {
        name = "PIGEON HOLE";
    };

    bool supportsDistanceScore(const int& maxScore) const override {
        return maxScore >= 1 && maxScore <= 4;
    }
};

// ============================================================================
// CLASS MINU SEARCHSTRATEGY
// ============================================================================

/**
 * A concrete derived class of SearchStrategy. This search strategy contains
 * the search schemes created by Hato.
 */
class MinUSearchStrategy : public SearchStrategy {
  private:
    const std::vector<Search> ED1 = {
        Search::makeSearch({0, 1}, {0, 0}, {0, 1}, 0),
        Search::makeSearch({1, 0}, {0, 0}, {0, 1}, 1)};
    const std::vector<Search> ED2 = {
        Search::makeSearch({0, 1, 2}, {0, 1, 1}, {0, 2, 2}, 0),
        Search::makeSearch({1, 0, 2}, {0, 0, 0}, {0, 1, 2}, 1),
        Search::makeSearch({2, 1, 0}, {0, 0, 2}, {0, 1, 2}, 2)};

    const std::vector<Search> ED3 = {
        Search::makeSearch({0, 1, 2, 3}, {0, 0, 0, 0}, {0, 1, 3, 3}, 0),
        Search::makeSearch({1, 0, 2, 3}, {0, 1, 1, 1}, {0, 1, 3, 3}, 1),
        Search::makeSearch({2, 3, 1, 0}, {0, 0, 0, 2}, {0, 1, 3, 3}, 2),
        Search::makeSearch({3, 2, 1, 0}, {0, 1, 1, 3}, {0, 1, 3, 3}, 3)};

    const std::vector<Search> ED4 = {
        Search::makeSearch({0, 1, 2, 3, 4}, {0, 0, 2, 2, 2}, {0, 2, 2, 4, 4},
                           0),
        Search::makeSearch({1, 2, 0, 3, 4}, {0, 0, 0, 0, 0}, {0, 1, 2, 4, 4},
                           1),
        Search::makeSearch({2, 1, 0, 3, 4}, {0, 1, 1, 1, 1}, {0, 1, 2, 4, 4},
                           2),
        Search::makeSearch({3, 4, 2, 1, 0}, {0, 0, 0, 0, 3}, {0, 1, 4, 4, 4},
                           3),
        Search::makeSearch({4, 3, 2, 1, 0}, {0, 1, 1, 1, 4}, {0, 1, 4, 4, 4},
                           4)};

    const std::vector<Search> ED5 = {
        Search::makeSearch({0, 1, 2, 3, 4, 5}, {0, 0, 0, 2, 2, 2},
                           {0, 1, 3, 5, 5, 5}, 0),
        Search::makeSearch({1, 0, 2, 3, 4, 5}, {0, 1, 1, 3, 3, 3},
                           {0, 1, 3, 5, 5, 5}, 1),
        Search::makeSearch({2, 3, 1, 0, 4, 5}, {0, 0, 0, 0, 0, 0},
                           {0, 1, 3, 3, 5, 5}, 2),
        Search::makeSearch({3, 2, 1, 0, 4, 5}, {0, 1, 1, 1, 1, 1},
                           {0, 1, 3, 3, 5, 5}, 3),
        Search::makeSearch({4, 5, 3, 2, 1, 0}, {0, 0, 0, 0, 0, 4},
                           {0, 1, 3, 5, 5, 5}, 4),
        Search::makeSearch({5, 4, 3, 2, 1, 0}, {0, 1, 1, 1, 1, 5},
                           {0, 1, 3, 5, 5, 5}, 5)};

    const std::vector<Search> ED6 = {
        Search::makeSearch({0, 1, 2, 3, 4, 5, 6}, {0, 0, 2, 2, 2, 2, 6},
                           {0, 2, 2, 6, 6, 6, 6}, 0),
        Search::makeSearch({1, 2, 0, 3, 4, 5, 6}, {0, 1, 1, 1, 1, 1, 5},
                           {0, 1, 2, 6, 6, 6, 6}, 1),
        Search::makeSearch({2, 1, 0, 3, 4, 5, 6}, {0, 0, 0, 0, 0, 0, 4},
                           {0, 1, 2, 6, 6, 6, 6}, 2),
        Search::makeSearch({3, 4, 5, 6, 2, 1, 0}, {0, 0, 0, 0, 0, 0, 0},
                           {0, 1, 3, 3, 6, 6, 6}, 3),
        Search::makeSearch({4, 3, 5, 6, 2, 1, 0}, {0, 1, 1, 1, 1, 1, 1},
                           {0, 1, 3, 3, 6, 6, 6}, 4),
        Search::makeSearch({5, 6, 4, 3, 2, 1, 0}, {0, 0, 0, 2, 2, 2, 2},
                           {0, 1, 3, 3, 6, 6, 6}, 5),
        Search::makeSearch({6, 5, 4, 3, 2, 1, 0}, {0, 1, 1, 3, 3, 3, 3},
                           {0, 1, 3, 3, 6, 6, 6}, 6)};
    const std::vector<Search> ED7 = {
        Search::makeSearch({0, 1, 2, 3, 4, 5, 6, 7}, {0, 0, 0, 0, 0, 0, 0, 0},
                           {0, 1, 3, 3, 7, 7, 7, 7}, 0),
        Search::makeSearch({1, 0, 2, 3, 4, 5, 6, 7}, {0, 1, 1, 1, 1, 1, 1, 1},
                           {0, 1, 3, 3, 7, 7, 7, 7}, 1),
        Search::makeSearch({2, 3, 1, 0, 4, 5, 6, 7}, {0, 0, 0, 2, 2, 2, 2, 2},
                           {0, 1, 3, 3, 7, 7, 7, 7}, 2),
        Search::makeSearch({3, 2, 1, 0, 4, 5, 6, 7}, {0, 1, 1, 3, 3, 3, 3, 3},
                           {0, 1, 3, 3, 7, 7, 7, 7}, 3),
        Search::makeSearch({4, 5, 6, 7, 3, 2, 1, 0}, {0, 0, 0, 0, 0, 0, 0, 4},
                           {0, 1, 3, 3, 7, 7, 7, 7}, 4),
        Search::makeSearch({5, 4, 6, 7, 3, 2, 1, 0}, {0, 1, 1, 1, 1, 1, 1, 5},
                           {0, 1, 3, 3, 7, 7, 7, 7}, 5),
        Search::makeSearch({6, 7, 5, 4, 3, 2, 1, 0}, {0, 0, 0, 2, 2, 2, 2, 6},
                           {0, 1, 3, 3, 7, 7, 7, 7}, 6),
        Search::makeSearch({7, 6, 5, 4, 3, 2, 1, 0}, {0, 1, 1, 3, 3, 3, 3, 7},
                           {0, 1, 3, 3, 7, 7, 7, 7}, 7)};

    const std::vector<std::vector<Search>> schemePerED = {ED1, ED2, ED3, ED4,
                                                          ED5, ED6, ED7};
    uint32_t calculateNumParts(unsigned int maxED) const override {
        return maxED + 1;
    }

    length_t getMaxSupportedDistance() const override {
        return 7;
    }

  protected:
    virtual const std::vector<Search>&
    createSearches(unsigned int maxED,
                   const std::vector<SARangePair>& ranges) const override {
        assert(maxED <= 7);
        assert(maxED >= 1);
        return schemePerED[maxED - 1];
    }

  public:
    MinUSearchStrategy(IndexInterface& index, PartitionStrategy p,
                       DistanceMetric metric, MappingMode mode,
                       SequencingMode sMode)
        : SearchStrategy(index, p, metric, mode, sMode) {
        name = "minU";
    };

    bool supportsDistanceScore(const int& maxScore) const override {
        return maxScore >= 1 && maxScore <= 7;
    }
};

/**
 * The default search strategy for Columba without dynamic selection. It uses
 * the MinU search schemes for k <=7 and the greedy search schemes for k > 7 and
 * k <= 13.
 */
class ColumbaSearchStrategy : public MinUSearchStrategy {
  public:
    ColumbaSearchStrategy(IndexInterface& index, PartitionStrategy p,
                          DistanceMetric metric, MappingMode mode,
                          SequencingMode sMode)
        : MinUSearchStrategy(index, p, metric, mode, sMode) {
        name = "Columba";
    };

    const std::vector<Search>&
    createSearches(unsigned int maxED,
                   const std::vector<SARangePair>& ranges) const override {
        if (maxED <= 7) {
            return MinUSearchStrategy::createSearches(maxED, ranges);
        } else {
            return schemePerEDHigh[maxED - 8];
        }
    }

    bool supportsDistanceScore(const int& maxScore) const override {
        return maxScore >= 1 && maxScore <= 13;
    }

  private:
    const std::vector<Search> ED8 = {
        Search::makeSearch({0, 1, 2, 3, 4, 5, 6, 7, 8},
                           {0, 1, 2, 3, 4, 5, 6, 7, 8},
                           {0, 2, 3, 5, 7, 8, 8, 8, 8}, 0),
        Search::makeSearch({1, 0, 2, 3, 4, 5, 6, 7, 8},
                           {0, 0, 1, 2, 3, 4, 5, 6, 7},
                           {0, 1, 4, 5, 7, 8, 8, 8, 8}, 1),
        Search::makeSearch({2, 1, 0, 3, 4, 5, 6, 7, 8},
                           {0, 0, 0, 0, 0, 1, 2, 3, 4},
                           {0, 1, 2, 6, 7, 8, 8, 8, 8}, 2),
        Search::makeSearch({3, 2, 1, 0, 4, 5, 6, 7, 8},
                           {0, 0, 0, 1, 1, 2, 3, 4, 5},
                           {0, 1, 3, 3, 7, 8, 8, 8, 8}, 3),
        Search::makeSearch({4, 3, 2, 1, 0, 5, 6, 7, 8},
                           {0, 0, 0, 0, 2, 3, 4, 5, 6},
                           {0, 1, 3, 4, 4, 8, 8, 8, 8}, 4),
        Search::makeSearch({5, 4, 3, 2, 1, 0, 6, 7, 8},
                           {0, 0, 0, 0, 0, 0, 0, 1, 2},
                           {0, 1, 2, 5, 5, 5, 8, 8, 8}, 5),
        Search::makeSearch({6, 5, 4, 3, 2, 1, 0, 7, 8},
                           {0, 0, 0, 0, 0, 0, 1, 2, 3},
                           {0, 1, 2, 5, 6, 6, 6, 8, 8}, 6),
        Search::makeSearch({7, 6, 5, 4, 3, 2, 1, 0, 8},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0},
                           {0, 1, 2, 3, 7, 7, 7, 7, 8}, 7),
        Search::makeSearch({8, 7, 6, 5, 4, 3, 2, 1, 0},
                           {0, 0, 0, 0, 0, 0, 0, 0, 1},
                           {0, 1, 2, 3, 8, 8, 8, 8, 8}, 8),
    };

    const std::vector<Search> ED9 = {
        Search::makeSearch({0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                           {0, 1, 2, 3, 4, 9, 9, 9, 9, 9}, 0),
        Search::makeSearch({1, 2, 3, 4, 5, 6, 7, 8, 9, 0},
                           {0, 0, 0, 0, 0, 0, 0, 0, 1, 2},
                           {0, 1, 2, 3, 8, 8, 8, 8, 8, 9}, 1),
        Search::makeSearch({2, 3, 4, 5, 6, 7, 8, 9, 1, 0},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
                           {0, 1, 2, 3, 7, 7, 7, 7, 9, 9}, 2),
        Search::makeSearch({3, 4, 5, 6, 7, 8, 9, 2, 1, 0},
                           {0, 0, 0, 0, 0, 0, 1, 2, 3, 4},
                           {0, 1, 2, 6, 6, 6, 6, 9, 9, 9}, 3),
        Search::makeSearch({4, 5, 6, 7, 8, 9, 3, 2, 1, 0},
                           {0, 0, 0, 0, 0, 0, 0, 1, 2, 3},
                           {0, 1, 2, 5, 5, 5, 9, 9, 9, 9}, 4),
        Search::makeSearch({5, 6, 7, 8, 9, 4, 3, 2, 1, 0},
                           {0, 0, 0, 0, 3, 4, 5, 6, 7, 8},
                           {0, 1, 4, 4, 4, 9, 9, 9, 9, 9}, 5),
        Search::makeSearch({6, 7, 8, 9, 5, 4, 3, 2, 1, 0},
                           {0, 0, 0, 0, 0, 1, 2, 3, 4, 5},
                           {0, 1, 3, 3, 7, 9, 9, 9, 9, 9}, 6),
        Search::makeSearch({7, 8, 9, 6, 5, 4, 3, 2, 1, 0},
                           {0, 0, 2, 3, 4, 5, 6, 7, 8, 9},
                           {0, 2, 2, 5, 7, 9, 9, 9, 9, 9}, 7),
        Search::makeSearch({8, 9, 7, 6, 5, 4, 3, 2, 1, 0},
                           {0, 0, 0, 1, 1, 2, 3, 4, 5, 6},
                           {0, 1, 3, 5, 7, 9, 9, 9, 9, 9}, 8),
        Search::makeSearch({9, 8, 7, 6, 5, 4, 3, 2, 1, 0},
                           {0, 1, 1, 2, 2, 3, 4, 5, 6, 7},
                           {0, 1, 3, 5, 7, 9, 9, 9, 9, 9}, 9),
    };

    const std::vector<Search> ED10 = {
        Search::makeSearch({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10},
                           {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10},
                           {0, 2, 3, 5, 7, 9, 10, 10, 10, 10, 10}, 0),
        Search::makeSearch({1, 0, 2, 3, 4, 5, 6, 7, 8, 9, 10},
                           {0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
                           {0, 1, 4, 5, 7, 9, 10, 10, 10, 10, 10}, 1),
        Search::makeSearch({2, 1, 0, 3, 4, 5, 6, 7, 8, 9, 10},
                           {0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6},
                           {0, 1, 2, 6, 7, 9, 10, 10, 10, 10, 10}, 2),
        Search::makeSearch({3, 2, 1, 0, 4, 5, 6, 7, 8, 9, 10},
                           {0, 0, 0, 1, 1, 2, 3, 4, 5, 6, 7},
                           {0, 1, 3, 3, 7, 9, 10, 10, 10, 10, 10}, 3),
        Search::makeSearch({4, 3, 2, 1, 0, 5, 6, 7, 8, 9, 10},
                           {0, 0, 0, 0, 2, 3, 4, 5, 6, 7, 8},
                           {0, 1, 3, 4, 4, 9, 10, 10, 10, 10, 10}, 4),
        Search::makeSearch({5, 4, 3, 2, 1, 0, 6, 7, 8, 9, 10},
                           {0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4},
                           {0, 1, 2, 5, 5, 5, 10, 10, 10, 10, 10}, 5),
        Search::makeSearch({6, 5, 4, 3, 2, 1, 0, 7, 8, 9, 10},
                           {0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5},
                           {0, 1, 2, 5, 6, 6, 6, 10, 10, 10, 10}, 6),
        Search::makeSearch({7, 6, 5, 4, 3, 2, 1, 0, 8, 9, 10},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2},
                           {0, 1, 2, 3, 7, 7, 7, 7, 10, 10, 10}, 7),
        Search::makeSearch({8, 7, 6, 5, 4, 3, 2, 1, 0, 9, 10},
                           {0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3},
                           {0, 1, 2, 3, 8, 8, 8, 8, 8, 10, 10}, 8),
        Search::makeSearch({9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 10},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                           {0, 1, 2, 3, 4, 9, 9, 9, 9, 9, 10}, 9),
        Search::makeSearch({10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
                           {0, 1, 2, 3, 4, 10, 10, 10, 10, 10, 10}, 10),
    };

    const std::vector<Search> ED11 = {
        Search::makeSearch({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                           {0, 1, 2, 3, 4, 5, 11, 11, 11, 11, 11, 11}, 0),
        Search::makeSearch({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2},
                           {0, 1, 2, 3, 4, 10, 10, 10, 10, 10, 10, 11}, 1),
        Search::makeSearch({2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 1, 0},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
                           {0, 1, 2, 3, 4, 9, 9, 9, 9, 9, 11, 11}, 2),
        Search::makeSearch({3, 4, 5, 6, 7, 8, 9, 10, 11, 2, 1, 0},
                           {0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4},
                           {0, 1, 2, 3, 8, 8, 8, 8, 8, 11, 11, 11}, 3),
        Search::makeSearch({4, 5, 6, 7, 8, 9, 10, 11, 3, 2, 1, 0},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3},
                           {0, 1, 2, 3, 7, 7, 7, 7, 11, 11, 11, 11}, 4),
        Search::makeSearch({5, 6, 7, 8, 9, 10, 11, 4, 3, 2, 1, 0},
                           {0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6},
                           {0, 1, 2, 5, 6, 6, 6, 11, 11, 11, 11, 11}, 5),
        Search::makeSearch({6, 7, 8, 9, 10, 11, 5, 4, 3, 2, 1, 0},
                           {0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5},
                           {0, 1, 2, 5, 5, 5, 11, 11, 11, 11, 11, 11}, 6),
        Search::makeSearch({7, 8, 9, 10, 11, 6, 5, 4, 3, 2, 1, 0},
                           {0, 0, 0, 0, 2, 3, 4, 5, 6, 7, 8, 9},
                           {0, 1, 3, 4, 4, 9, 11, 11, 11, 11, 11, 11}, 7),
        Search::makeSearch({8, 9, 10, 11, 7, 6, 5, 4, 3, 2, 1, 0},
                           {0, 0, 0, 1, 1, 2, 3, 4, 5, 6, 7, 8},
                           {0, 1, 3, 3, 7, 9, 11, 11, 11, 11, 11, 11}, 8),
        Search::makeSearch({9, 10, 11, 8, 7, 6, 5, 4, 3, 2, 1, 0},
                           {0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7},
                           {0, 1, 2, 6, 7, 9, 11, 11, 11, 11, 11, 11}, 9),
        Search::makeSearch({10, 11, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0},
                           {0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10},
                           {0, 1, 4, 5, 7, 9, 11, 11, 11, 11, 11, 11}, 10),
        Search::makeSearch({11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0},
                           {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11},
                           {0, 2, 3, 5, 7, 9, 11, 11, 11, 11, 11, 11}, 11)};
    const std::vector<Search> ED12 = {
        Search::makeSearch({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12},
                           {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12},
                           {0, 2, 3, 5, 7, 9, 11, 12, 12, 12, 12, 12, 12}, 0),
        Search::makeSearch({1, 0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12},
                           {0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11},
                           {0, 1, 4, 5, 7, 9, 11, 12, 12, 12, 12, 12, 12}, 1),
        Search::makeSearch({2, 1, 0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12},
                           {0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8},
                           {0, 1, 2, 6, 7, 9, 11, 12, 12, 12, 12, 12, 12}, 2),
        Search::makeSearch({3, 2, 1, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12},
                           {0, 0, 0, 1, 1, 2, 3, 4, 5, 6, 7, 8, 9},
                           {0, 1, 3, 3, 7, 9, 11, 12, 12, 12, 12, 12, 12}, 3),
        Search::makeSearch({4, 3, 2, 1, 0, 5, 6, 7, 8, 9, 10, 11, 12},
                           {0, 0, 0, 0, 2, 3, 4, 5, 6, 7, 8, 9, 10},
                           {0, 1, 3, 4, 4, 9, 11, 12, 12, 12, 12, 12, 12}, 4),
        Search::makeSearch({5, 4, 3, 2, 1, 0, 6, 7, 8, 9, 10, 11, 12},
                           {0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6},
                           {0, 1, 2, 5, 5, 5, 11, 12, 12, 12, 12, 12, 12}, 5),
        Search::makeSearch({6, 5, 4, 3, 2, 1, 0, 7, 8, 9, 10, 11, 12},
                           {0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7},
                           {0, 1, 2, 5, 6, 6, 6, 12, 12, 12, 12, 12, 12}, 6),
        Search::makeSearch({7, 6, 5, 4, 3, 2, 1, 0, 8, 9, 10, 11, 12},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4},
                           {0, 1, 2, 3, 7, 7, 7, 7, 12, 12, 12, 12, 12}, 7),
        Search::makeSearch({8, 7, 6, 5, 4, 3, 2, 1, 0, 9, 10, 11, 12},
                           {0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5},
                           {0, 1, 2, 3, 8, 8, 8, 8, 8, 12, 12, 12, 12}, 8),
        Search::makeSearch({9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 10, 11, 12},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2},
                           {0, 1, 2, 3, 4, 9, 9, 9, 9, 9, 12, 12, 12}, 9),
        Search::makeSearch({10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 11, 12},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3},
                           {0, 1, 2, 3, 4, 10, 10, 10, 10, 10, 10, 12, 12}, 10),
        Search::makeSearch({11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 12},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                           {0, 1, 2, 3, 4, 5, 11, 11, 11, 11, 11, 11, 12}, 11),
        Search::makeSearch({12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
                           {0, 1, 2, 3, 4, 5, 12, 12, 12, 12, 12, 12, 12}, 12)};

    const std::vector<Search> ED13 = {
        Search::makeSearch({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                           {0, 1, 2, 3, 4, 5, 6, 13, 13, 13, 13, 13, 13, 13},
                           0),
        Search::makeSearch({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 0},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2},
                           {0, 1, 2, 3, 4, 5, 12, 12, 12, 12, 12, 12, 12, 13},
                           1),
        Search::makeSearch({2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 1, 0},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
                           {0, 1, 2, 3, 4, 5, 11, 11, 11, 11, 11, 11, 13, 13},
                           2),
        Search::makeSearch({3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 2, 1, 0},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4},
                           {0, 1, 2, 3, 4, 10, 10, 10, 10, 10, 10, 13, 13, 13},
                           3),
        Search::makeSearch({4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 3, 2, 1, 0},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3},
                           {0, 1, 2, 3, 4, 9, 9, 9, 9, 9, 13, 13, 13, 13}, 4),
        Search::makeSearch({5, 6, 7, 8, 9, 10, 11, 12, 13, 4, 3, 2, 1, 0},
                           {0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6},
                           {0, 1, 2, 3, 8, 8, 8, 8, 8, 13, 13, 13, 13, 13}, 5),
        Search::makeSearch({6, 7, 8, 9, 10, 11, 12, 13, 5, 4, 3, 2, 1, 0},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5},
                           {0, 1, 2, 3, 7, 7, 7, 7, 13, 13, 13, 13, 13, 13}, 6),
        Search::makeSearch({7, 8, 9, 10, 11, 12, 13, 6, 5, 4, 3, 2, 1, 0},
                           {0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8},
                           {0, 1, 2, 5, 6, 6, 6, 13, 13, 13, 13, 13, 13, 13},
                           7),
        Search::makeSearch({8, 9, 10, 11, 12, 13, 7, 6, 5, 4, 3, 2, 1, 0},
                           {0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7},
                           {0, 1, 2, 5, 5, 5, 11, 13, 13, 13, 13, 13, 13, 13},
                           8),
        Search::makeSearch({9, 10, 11, 12, 13, 8, 7, 6, 5, 4, 3, 2, 1, 0},
                           {0, 0, 0, 0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11},
                           {0, 1, 3, 4, 4, 9, 11, 13, 13, 13, 13, 13, 13, 13},
                           9),
        Search::makeSearch({10, 11, 12, 13, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0},
                           {0, 0, 0, 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10},
                           {0, 1, 3, 3, 7, 9, 11, 13, 13, 13, 13, 13, 13, 13},
                           10),
        Search::makeSearch({11, 12, 13, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0},
                           {0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
                           {0, 1, 2, 6, 7, 9, 11, 13, 13, 13, 13, 13, 13, 13},
                           11),
        Search::makeSearch({12, 13, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0},
                           {0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12},
                           {0, 1, 4, 5, 7, 9, 11, 13, 13, 13, 13, 13, 13, 13},
                           12),
        Search::makeSearch({13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0},
                           {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13},
                           {0, 2, 3, 5, 7, 9, 11, 13, 13, 13, 13, 13, 13, 13},
                           13)};

    const std::vector<std::vector<Search>> schemePerEDHigh = {ED8,  ED9,  ED10,
                                                              ED11, ED12, ED13};

    length_t getMaxSupportedDistance() const override {
        return 13;
    }
};

/**
 * A search strategy based on the minU seach schemes, and the greedy solution
 * seach schemes that uses dynamic selection. This is used by default in
 * Columba. It is based on ColumbaSearchStrategy but adds the search schemes
 * that have the critical search start with a different part.
 */
class DynamicColumbaStrategy : public MultipleSchemesStrategy {

  private:
    const static std::vector<Search> getMidSearch2() {
        return {Search::makeSearch({2, 1, 0}, {0, 1, 1}, {0, 2, 2}, 0),
                Search::makeSearch({1, 2, 0}, {0, 0, 0}, {0, 1, 2}, 1),
                Search::makeSearch({0, 1, 2}, {0, 0, 2}, {0, 1, 2}, 2)};
    }
    // keep the middle search schemes for even k

    const static std::vector<Search> getMidSearch4() {
        return {Search::makeSearch({0, 1, 2, 3, 4}, {0, 1, 1, 1, 4},
                                   {0, 1, 4, 4, 4}, 0),
                Search::makeSearch({1, 0, 2, 3, 4}, {0, 0, 0, 0, 3},
                                   {0, 1, 4, 4, 4}, 1),
                Search::makeSearch({2, 3, 4, 1, 0}, {0, 1, 1, 1, 1},
                                   {0, 2, 2, 4, 4}, 2),
                Search::makeSearch({3, 2, 4, 1, 0}, {0, 0, 0, 0, 0},
                                   {0, 1, 2, 4, 4}, 3),
                Search::makeSearch({4, 3, 2, 1, 0}, {0, 0, 2, 2, 2},
                                   {0, 1, 2, 4, 4}, 4)};
    }

    const static std::vector<Search> getMidSearch6() {
        return {Search::makeSearch({0, 1, 2, 3, 4, 5, 6}, {0, 1, 1, 1, 1, 1, 5},
                                   {0, 1, 2, 6, 6, 6, 6}, 0),
                Search::makeSearch({1, 0, 2, 3, 4, 5, 6}, {0, 0, 0, 0, 0, 0, 4},
                                   {0, 1, 2, 6, 6, 6, 6}, 1),
                Search::makeSearch({2, 1, 0, 3, 4, 5, 6}, {0, 0, 2, 2, 2, 2, 6},
                                   {0, 2, 2, 6, 6, 6, 6}, 2),
                Search::makeSearch({3, 4, 5, 6, 2, 1, 0}, {0, 0, 0, 2, 2, 2, 2},
                                   {0, 1, 3, 3, 6, 6, 6}, 3),
                Search::makeSearch({4, 3, 5, 6, 2, 1, 0}, {0, 1, 1, 3, 3, 3, 3},
                                   {0, 1, 3, 3, 6, 6, 6}, 4),
                Search::makeSearch({5, 6, 4, 3, 2, 1, 0}, {0, 0, 0, 0, 0, 0, 0},
                                   {0, 1, 3, 3, 6, 6, 6}, 5),
                Search::makeSearch({6, 5, 4, 3, 2, 1, 0}, {0, 1, 1, 1, 1, 1, 1},
                                   {0, 1, 3, 3, 6, 6, 6}, 6)};
    }

    /**
     * @brief Create a new dynamic Columba strategy.
     */
    DynamicColumbaStrategy(SearchStrategy* strategy)
        : MultipleSchemesStrategy(strategy) {
    }

  public:
    /**
     * Statically create a pointer to dynamic columba strategy.
     * @param index the index to use
     * @param p the partition strategy to use
     * @param metric the distance metric to use
     * @param mode the mapping mode to use
     * @param sMode the sequencing mode to use
     */
    static DynamicColumbaStrategy*
    createDynamicColumbaStrategy(IndexInterface& index, PartitionStrategy p,
                                 DistanceMetric metric, MappingMode mode,
                                 SequencingMode sMode) {
        ColumbaSearchStrategy columbaStrategy(index, p, metric, mode, sMode);
        auto instance = new DynamicColumbaStrategy(&columbaStrategy);

        instance->addScheme(SearchScheme(getMidSearch2(), 2));
        instance->addScheme(SearchScheme(getMidSearch4(), 4));
        SearchScheme scheme6 = SearchScheme(getMidSearch6(), 6);
        instance->addScheme(scheme6);
        instance->addScheme(scheme6.mirrorPiStrings());
        return instance;
    }
};

/**
 * A search strategy based on a custom search strategy provided by the user but
 * with dynamic selection between the schemes and their mirrored variants.
 * @see CustomSearchStrategy for the requirements of the definition of the
 * custom search strategy.
 */
class DynamicCustomStrategy : public MultipleSchemesStrategy {
  private:
    /**
     * @brief Create a new dynamic custom strategy.
     */
    DynamicCustomStrategy(SearchStrategy* strategy)
        : MultipleSchemesStrategy(strategy) {
    }

  public:
    /**
     * Statically create a pointer to the dynamic variant of a custom search
     * strategy. @see CustomSearchStrategy for the requirements in the folder.
     * @param index the index to use
     * @param pathToFolder the path to the folder where the custom strategy is
     * defined
     * @param p the partition strategy to use
     * @param metric the distance metric to use
     * @param mode the mapping mode to use
     * @param sMode the sequencing mode to use
     * @param verbose whether to print debug information
     */
    static DynamicCustomStrategy* createDynamicCustomSearchStrategy(
        IndexInterface& index, const std::string& pathToFolder,
        PartitionStrategy p, DistanceMetric metric, MappingMode mode,
        SequencingMode sMode, bool verbose = false) {
        CustomSearchStrategy customStrategy(index, pathToFolder, p, metric,
                                            mode, sMode, verbose);
        auto instance = new DynamicCustomStrategy(&customStrategy);
        return instance;
    }
};

#endif
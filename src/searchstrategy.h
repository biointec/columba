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

#ifndef SEARCHSTRATEGY_H
#define SEARCHSTRATEGY_H

#define MAX_K 13 // the maximum number of allowed errors

#include "fmindex.h"
#include "fmindexhelpers.h"

#include <sys/stat.h>

#define Distribution std::vector<int>

// An enum for partion strategy
enum PartitionStrategy { UNIFORM, STATIC, DYNAMIC };
// An enum for which distance metric to use
enum DistanceMetric { HAMMING, EDITNAIVE, EDITOPTIMIZED };

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
                                            std::vector<Substring>&, Counters&,
                                            const int&) const;

typedef std::vector<TextOcc> (SearchStrategy::*FilterPtr)(Occurrences&,
                                                          length_t,
                                                          Counters&) const;

class SearchStrategy {
  protected:
    FMIndex& index; // pointer to the index of the text that is searched

    // variables for getting info about strategy used
    PartitionStrategy partitionStrategy; // the partitioning strategy
    DistanceMetric distanceMetric;       // which distance metric to use
    std::string name; // the name of the particular search strategy

    // pointers for correct partitioning and correct distance metric
    PartitionPtr partitionPtr; // pointer to the partition method
    StartIdxPtr startIdxPtr;   // pointer to start method (hamming or
                               // (naive/optimized) edit distance)

    FilterPtr
        filterPtr; // pointer to the correct filter (edit or hamming based)

    // ----------------------------------------------------------------------------
    // CONSTRUCTOR
    // ----------------------------------------------------------------------------

    /**
     * Constructor
     * @param argument, pointer to the bidirectional FM index to use
     * @param p, partition strategy
     * @param edit, true if edit distance should be used, false if hamming
     * distance should be used
     */
    SearchStrategy(FMIndex& argument, PartitionStrategy p,
                   DistanceMetric distanceMetric);

    // ----------------------------------------------------------------------------
    // SANITY CHECKS
    // ----------------------------------------------------------------------------

    /**
     * Static function which generates all error distributions with P parts and
     * K errors.
     * @param P the number of parts
     * @param K the number of allowed errors
     * @param distributions vector to store the error distributions in
     */
    static void genErrorDistributions(int P, int K,
                                      std::vector<Distribution>& distributions);

    /**
     * Static function to check if a particular search scheme covers all
     * distributions.
     * @param distributions the error distributions to check
     * @param scheme the search scheme to check
     * @param verbose if true the details about which search covers which
     * distribution will be written to stdout
     */
    static bool
    coversDistributions(const std::vector<Distribution>& distributions,
                        const std::vector<Search>& scheme, bool verbose);

    // ----------------------------------------------------------------------------
    // PARTITIONING
    // ----------------------------------------------------------------------------

    /**
     * Splits the pattern into numParts parts, either by uniform range or
     * uniform size
     * @param pattern the pattern to be split
     * @param parts the vector containing the substrings of this pattern,
     * will be cleared and filled during the execution of this method. If
     * the splitting fails for some reason, the vector will be empty
     * @param numparts, how many parts are needed
     * @param maxScore the maximum allowed edit distance
     * @param exactMatchRanges, a vector corresponding to the ranges for the
     * exact matches of the parts, will be cleared and filled during the
     * execution
     */
    void partition(const std::string& pattern, std::vector<Substring>& parts,
                   const int& numParts, const int& maxScore,
                   std::vector<SARangePair>& exactMatchRanges,
                   Counters& counters) const;
    /**
     * Calculates the number of parts for a certain max edit distance. This
     * calculation is strategy dependent
     * @param maxED the maximal allowed edit distance for the alligning
     */
    virtual uint calculateNumParts(unsigned int maxED) const = 0;

    // Uniform Partitioning

    /**
     * Splits the pattern into numParts parts, such that each part has the
     * same size
     * @param pattern the pattern to be split
     * @param parts an empty vector which will be filled with the different
     * parts
     * @param numparts, how many parts are needed
     * @param maxScore, the maximum allowed edit distance
     * @param exactMatchRanges, a vector corresponding to the ranges for the
     * exact matches of the parts, will be cleared and filled during the
     * execution
     */
    void partitionUniform(const std::string& pattern,
                          std::vector<Substring>& parts, const int& numParts,
                          const int& maxScore,
                          std::vector<SARangePair>& exactMatchRanges,
                          Counters& counters) const;

    // Optimal static partitioning

    /**
     * Splits the pattern into numParts parts, such that each search carries
     * the same weight (on average)
     * @param pattern the pattern to be split
     * @param parts an empty vector which will be filled with the different
     * parts
     * @param numparts, how many parts are needed
     * @param exactMatchRanges, a vector corresponding to the ranges for the
     * exact matches of the parts, will be cleared and filled during the
     * execution
     */
    void partitionOptimalStatic(const std::string& pattern,
                                std::vector<Substring>& parts,
                                const int& numParts, const int& maxScore,
                                std::vector<SARangePair>& exactMatchRanges,
                                Counters& counters) const;

    /**
     * Helper function for optimal static partitioning. This function
     * creates the optimal static parts
     * @param pattern, the pattern to partition
     * @param parts, empty vector to which the parts are added
     * @param numParts, how many parts there need to be in the partition
     */
    void setParts(const std::string& pattern, std::vector<Substring>& parts,
                  const int& numParts, const int& maxScore) const;

    /**
     * Function that retrieves the begin positions for optimal static
     * partitioning. If derived class does not implement this function then
     * uniform positions are given.
     * @param numparts, how many parts are needed
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

    // Dynamic Partitioning

    /**
     * Splits the pattern into numParts parts, such that each part has
     * (approximately) the same range. The exactMatchRanges are also
     * calculated.
     * @param pattern the pattern to be split
     * @param parts an empty vector which will be filled with the different
     * parts
     * @param numparts, how many parts are needed
     * @param exactMatchRanges, a vector corresponding to the ranges for the
     * exact matches of the parts, will be cleared and filled during the
     * execution
     */
    void partitionDynamic(const std::string& pattern,
                          std::vector<Substring>& parts, const int& numParts,
                          const int& maxScore,
                          std::vector<SARangePair>& exactMatchRanges,
                          Counters& counters) const;

    /**
     * Function that retrieves the seeding positions for dynamic partitioning.
     * If derived class does not implement this function then uniform seeds are
     * given.
     * @param numParts, how many parts are needed
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
     * @param pattern, the pattern to partition
     * @param parts, empty vector tro which the seeds are added
     * @param numparts, how many parts are needed
     * @param exactMatchRanges, a vector corresponding to the ranges for the
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
     * @param numparts how many parts are needed
     * @returns vector with weights
     */
    virtual const std::vector<int> getWeights(const int& numParts,
                                              const int& maxScore) const {
        std::vector<int> w(numParts, 1);
        w.front() = 2;
        w.back() = 2;
        return w;
    }

    /**
     * Helper function for dynamic partitioning. This function extends the parts
     * so that nothing of the pattern is not allocated to any part. This does
     * not keep track of the ranges over the suffix array, so should only be
     * called if this does not matter (e.g. when the parts that can be extended
     * all correspond to empty ranges)
     * @param pattern, the pattern that is split
     * @param parts, the current parts, they are updated so that all characters
     * of pattern are part of exactly one part
     */
    void extendParts(const std::string& pattern,
                     std::vector<Substring>& parts) const;

    // ----------------------------------------------------------------------------
    // (APPROXIMATE) MATCHING
    // ----------------------------------------------------------------------------

    /**
     * Creates all searches for this specific strategy. This is strategy
     * dependent
     * @param maxED the maximal allowed edit distance for the aligning
     */
    virtual const std::vector<Search>&
    createSearches(unsigned int maxED,
                   const std::vector<SARangePair>& exactMatchRanges) const = 0;

    /**
     * Executes the search recursively. If U[0] != 1, then the search will
     * start at pi[0], else the search will start with idx i and U[i]!=0 and
     * U[j]=0 with j < i
     * @param s, the search to follow
     * @param parts, the parts of the pattern
     * @param allMatches, vector to add occurrences to
     * @param exactMatchRanges, a vector corresponding to the ranges for the
     * exact matches of the parts
     */
    void doRecSearch(const Search& s, std::vector<Substring>& parts,
                     Occurrences& occ,
                     const std::vector<SARangePair>& exactMatchRanges,
                     Counters& counters) const;

    /**
     * Filter the occurrences (both in-index and in-text) based on hamming
     * distance. The in-index occurrences are first converted to in-text
     * occurrences
     * @param occ the occurrences
     * @param maxED the maximal allowed edit distance
     * @param counters perfomance counters
     * @returns the filterd in-text occurrences ()
     */
    std::vector<TextOcc> filterHamming(Occurrences& occ, length_t maxED,
                                       Counters& counters) const {
        return index.getTextOccHamming(occ, counters);
    }
    /**
     * Filter the occurrences (both in-index and in-text) based on edit
     * distance. The in-index occurrences are first converted to in-text
     * occurrences
     * @param occ the occurrences
     * @param maxED the maximal allowed edit distance
     * @param counters perfomance counters
     * @returns the filterd in-text occurrences ()
     */
    std::vector<TextOcc> filterEdit(Occurrences& occ, length_t maxED,
                                    Counters& counters) const {
        return index.getUniqueTextOccurrences(occ, maxED, counters);
    }

    /**
     * Starts the index with hamming distance
     * @param s, the search to follow
     * @param startMatch, the startMatch that corresponds to the first piece
     * @param occ, vector to add occurrences to
     * @param parts, the parts of the pattern
     * @param idx, the index in the search to match next
     */
    void startIndexHamming(const Search& s, const FMOcc& startMatch,
                           Occurrences& occ, std::vector<Substring>& parts,
                           Counters& counters, const int& idx) const {
        index.recApproxMatchHamming(s, startMatch, occ, parts, counters, idx);
    }

    /**
     * Starts the index with edit distance and optimized alignment for the
     * edit distance metric
     * @param s, the search to follow
     * @param startMatch, the startMatch that corresponds to the first piece
     * @param occ, vector to add occurrences to
     * @param parts, the parts of the pattern
     * @param idx, the index in the search to match next
     */
    void startIndexEditOptimized(const Search& s, const FMOcc& startMatch,
                                 Occurrences& occ,
                                 std::vector<Substring>& parts,
                                 Counters& counters, const int& idx) const {
        index.recApproxMatchEditOptimizedEntry(s, startMatch, occ, parts,
                                               counters, idx);
    }

    /**
     * Starts the index with naive edit distance (= redundancy between parts
     * of a search)
     * @param s, the search to follow
     * @param startMatch, the startMatch that corresponds to the first piece
     * @param occ, vector to add occurrences to
     * @param parts, the parts of the pattern
     * @param idx, the index in the search to match next
     */
    void startIndexEditNaive(const Search& s, const FMOcc& startMatch,
                             Occurrences& occ, std::vector<Substring>& parts,
                             Counters& counters, const int& idx) const {
        index.recApproxMatchEditNaive(s, startMatch, occ, parts, counters, idx);
    }

    void generateSAM(std::vector<TextOcc>& occs, const std::string& pattern,
                     const std::string& ID, const std::string& qual,
                     bool revCompl) const {
        // TODO move logic to after both forward and revcompl have finished

        if (occs.empty()) {
            return;
        }

        // sort the occurrences on score
        std::sort(occs.begin(), occs.end(),
                  [](const TextOcc& a, const TextOcc& b) {
                      if (a.getDistance() != b.getDistance()) {
                          return a.getDistance() < b.getDistance();
                      }
                      return a.getRange().getBegin() < b.getRange().getBegin();
                  });

        // generate the sam lines
        bool notPrimary = false;

        length_t minscore = occs.front().getDistance();
        length_t nHits = 1;

        for (length_t i = 1; i < occs.size(); i++) {
            if (occs[i].getDistance() == minscore) {
                nHits++;
            } else {
                break;
            }
        }

        for (auto& t : occs) {
            t.generateSamSE(revCompl, notPrimary, ID, pattern, qual, nHits,
                            minscore);

            notPrimary = true;
        }
    }

  public:
    // ----------------------------------------------------------------------------
    //  Destructor
    // ----------------------------------------------------------------------------
    virtual ~SearchStrategy() {
    }
    // ----------------------------------------------------------------------------
    // INFORMATION
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
     * Retrieves the text of the index (for debugging purposes)
     */
    std::string getText() const {
        return index.getText();
    }

    /**
     * Get the tipping point for in-text verification
     */
    length_t getSwitchPoint() const {
        return index.getSwitchPoint();
    }

    /**
     * Matches a pattern approximately using this strategy
     * @param pattern the pattern to match
     * @param maxED the maximal allowed edit distance (or  hamming
     * distance)
     * @param counters the performance counters
     * @param ID the id of the read
     * @param qual the quality string of the read
     * @param revComp boolean to indicate if this the read or the reverse
     * complement of the read
     * @returns a vector with approximate occurrences in the text of this
     * pattern
     *
     */
    virtual std::vector<TextOcc>
    matchApprox(const std::string& pattern, length_t maxED, Counters& counters,
                const std::string& ID, const std::string& qual, bool revComp);
};

// ============================================================================
// CLASS CUSTOM SEARCHSTRATEGY
// ============================================================================

// This is a derived class of SearchStrategy. It creates a custom scheme using
// files provided by the user. It takes a specified folder in which a file
// "name.txt", containing the name of the scheme on the first line, and for each
// supported distance score a subfolder exists. Such a subfolder has as name the
// distance score. Each subfolder must contain at least a file "searches.txt".
// The different searches of the scheme for this distance score should be
// written on separate lines of this file. Each search consists out of three
// arrays, pi, L and U, the arrays are separated by a single space. Each array
// is written between curly braces {} and the different values are separated.
// The pi array must be zero-based.
//
// The different subfolders can also contain files for static and dynamic
// partitioning. The file "static_partitioning.txt" should consist out of one
// line of space- separated percentages (between 0 and 1 - exclusive). These
// percentages point to start positions of the second to the last part of the
// partitioned pattern (relative to the size of the pattern). Hence, if a search
// scheme partitions a pattern in k parts, then k+1 percentages should be
// provided. The file "dynamic_partitioning.txt" should consist out of two
// lines. The first line contains k-1 space- separated percentages. These
// percentages are the seeding positions of the middle parts (the first and last
// part are seeded at the begin and end and thus do not need a percentage). Note
// that this line can be empty if the pattern is partitioned into 2 parts. The
// second line contains k integers, where k is the number of parts. Each integer
// corresponds to the weight given to that part in the dynamic partitioning
// process.
class CustomSearchStrategy;

// Pointer to the correct getBegins() function for static partitioning
typedef const std::vector<double> (CustomSearchStrategy::*GetBeginsPtr)(
    const int& numParts, const int& maxScore) const;

// Pointer to the correct getSeedingPositions() function for dynamic
// partitioning
typedef const std::vector<double> (
    CustomSearchStrategy::*GetSeedingPostitionsPtr)(const int& numParts,
                                                    const int& maxScore) const;
// Pointer to the correct getWeights() function for dynamic partitioning
typedef const std::vector<int> (CustomSearchStrategy::*GetWeightsPtr)(
    const int& numParts, const int& maxScore) const;

class CustomSearchStrategy : public SearchStrategy {

  private:
    std::vector<std::vector<Search>>
        schemePerED; // the search schemes for each distance score,
    std::vector<bool>
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

    std::vector<std::vector<int>>
        weights; // the weights for dynamic partitioning per score

    std::vector<GetSeedingPostitionsPtr>
        seedingPointer; // pointer to the correct
                        // getSeedingPositions() function,
                        // either default or custom

    std::vector<GetWeightsPtr>
        weightsPointers; // pointer to the correct getWeigths() function, either
                         // default or custom

    /**
     * Retrieves the search scheme from a folder, also checks if the scheme
     * is valid
     * @param pathToFolder the path to the folder containing the search
     * scheme
     * @param verbose if the sanity check should be verbose
     */
    void getSearchSchemeFromFolder(std::string pathToFolder, bool verbose);

    /**
     * If the values provided for dynamic partitioning for the given max score
     * are valid (i.e. strictly increasing and between 0 and 1). Will throw a
     * runtime error if this is not the case
     * @param maxScore, the score to check
     */
    void sanityCheckDynamicPartitioning(const int& maxScore) const;

    /**
     * If the values provided for static partitioning for the given max score
     * are valid (i.e. strictly increasing and between 0 and 1). Will throw a
     * runtime error if this is not the case
     * @param maxScore, the score to check
     */
    void sanityCheckStaticPartitioning(const int& maxScore) const;

    /**
     * Parse the search from a line.
     * @param line the line to parse
     * @param idx the index of the search
     * @returns the parsed line as a search, if the line is not valid a runtime
     * error will be thrown.
     */
    Search makeSearch(const std::string& line, length_t idx) const;

    /**
     * Parses an array from a string.
     * @param vectorString the string to parse
     * @param vector the vector with the parsed array as values
     */
    void getVector(const std::string& vectorString,
                   std::vector<length_t>& vector) const;

    /**
     * Checks whether the connectivity property is satisfied for all searches.
     * Will throw a runtime error if one of these is not satisfied
     * @param verbose if the information about which search covers which error
     * distribution should be written to standard out
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
     * Gets the static positions in the custom manner (i.e. with values provided
     * by user in "static_partitioning.txt").
     * @param numParts the number of parts of the pattern
     * @param maxScore the maximal allowed score
     */
    const std::vector<double> getBeginsCustom(const int& numParts,
                                              const int& maxScore) const {
        return staticPositions[maxScore - 1];
    }
    /**
     * Overridden function of the base class. Retrieves the begin positions in
     * the custom way if values were provided in a "static_partitioning.txt"
     * file, otherwise the base class function will be called;
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
     * Overridden function of the base class. Retrieves the seeding positions in
     * the custom way if values were provided in a "dynamic_partitioning.txt"
     * file, otherwise the base class function will be called;
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
    const std::vector<int> getWeightsDefault(const int& numParts,
                                             const int& maxScore) const {
        return SearchStrategy::getWeights(numParts, maxScore);
    }

    /**
     * Gets the weights in the custom manner (i.e. with values
     * provided by user in "dynamic_partitioning.txt").
     * @param numParts the number of parts of the pattern
     * @param maxScore the maximal allowed score
     */
    const std::vector<int> getWeightsCustom(const int& numParts,
                                            const int& maxScore) const {
        return weights[maxScore - 1];
    }
    /**
     * Overridden function of the base class. Retrieves the weights positions in
     * the custom way if values were provided in a "dynamic_partitioning.txt"
     * file, otherwise the base class function will be called;
     */
    const std::vector<int> getWeights(const int& numParts,
                                      const int& maxScore) const override {
        assert(supportsMaxScore[maxScore - 1]);
        return (this->*weightsPointers[maxScore - 1])(numParts, maxScore);
    }

  public:
    CustomSearchStrategy(FMIndex& index, const std::string& pathToFolder,
                         PartitionStrategy p = DYNAMIC,
                         DistanceMetric metric = EDITOPTIMIZED,
                         bool verbose = false)
        : SearchStrategy(index, p, metric) {

        // resize and fill the vectors
        schemePerED.resize(MAX_K);
        supportsMaxScore.resize(MAX_K, false);
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

    uint calculateNumParts(unsigned int maxED) const {
        assert(supportsMaxScore[maxED - 1]);
        return schemePerED[maxED - 1][0].getNumParts();
    }
    const std::vector<Search>&
    createSearches(unsigned int maxED,
                   const std::vector<SARangePair>& ranges) const {
        assert(supportsMaxScore[maxED - 1]);
        return schemePerED[maxED - 1];
    }
};

// ============================================================================
// CLASS MULTIPLE SCHEMES
// ============================================================================

/**
 * @class MultipleSchemes
 *
 * @brief A class representing a collection of search schemes for approximate
 * string matching.
 *
 * This class manages multiple search schemes, each represented by a text file.
 * It calculates the number of parts and provides a method to choose the search
 * scheme to use based on the number of exact matches of the parts.
 */
class MultipleSchemes {

  private:
    unsigned int k; // the maximum number of errors for these schemes
    std::vector<SearchScheme> schemes; // the list of schemes
    unsigned int numParts = 0; // the number of parts (p) for this scheme

    /**
     * @brief Private helper function to read search schemes from a folder.
     *
     * This function reads search schemes from text files in a specified folder.
     * Each file should be named "scheme1.txt," "scheme2.txt," and so on.
     *
     * @param pathToFolder The path to the folder containing the scheme files.
     * @param verbose A flag indicating whether to display verbose output during
     * reading.
     */
    void getSchemesFromFolder(const std::string& pathToFolder, bool verbose) {
        // in this folder there should be several files:
        // scheme1.txt, scheme2.txt, scheme3.txt, ...

        std::string base = "scheme";
        int x = 1; // Starting value of x

        while (true) {
            std::string filename = base + std::to_string(x) + ".txt";
            std::string filePath = pathToFolder + "/" + filename;

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

            for (const auto& scheme : schemes) {
                if (scheme.getNumParts() != numParts) {
                    throw std::runtime_error(
                        "Not all schemes have same amount of parts in: " +
                        pathToFolder);
                }
            }
        }
    }

  public:
    /**
     * @brief Constructor for the MultipleSchemes class.
     *
     * Initializes the MultipleSchemes object by reading search schemes from a
     * folder.
     *
     * @param pathToFolder The path to the folder containing the scheme files.
     * @param k The maximum number of errors for these schemes.
     * @param verbose A flag indicating whether to display verbose output during
     * reading.
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
     * @param verbose A flag indicating whether to display verbose output during
     * reading.
     */
    MultipleSchemes(unsigned int k, bool verbose = false) : k(k) {
    }

    /**
     * @brief Calculate the number of parts for a given maximum edit distance.
     * This function should only be calles with maxED == k.
     *
     * @param maxED The maximum edit distance for which to calculate the number
     * of parts.
     * @return The number of parts (p) corresponding to the maximum edit
     * distance.
     */
    uint calculateNumParts(unsigned int maxED) const {
        assert(maxED == k);
        return numParts;
    }

    /**
     * @brief Create searches based on a maximum edit distance and specified
     * ranges (=number of exact matches of the searches).
     *
     * This method creates searches based on a given maximum edit distance and a
     * set of SARangePair ranges. It iterates over the schemes and chooses the
     * scheme for which the starting part of the critical search has the least
     * amount of exact matches.
     *
     * @param maxED The maximum edit distance for creating searches.
     * @param ranges A vector of SARangePair objects representing search ranges.
     * @return A vector of Search objects based on the specified parameters.
     */
    const std::vector<Search>&
    createSearches(unsigned int maxED,
                   const std::vector<SARangePair>& ranges) const {
        assert(maxED == k);
        assert(!isEmpty());
        // find the scheme for which the critical search has the smallest exact
        // range to start from

        // if the total number of exact matches of the parts is smaller than the
        // number of parts dynamic selection has too much overhead
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
};
// ============================================================================
// CLASS MULTIPLE SCHEME STRATEGY
// ============================================================================

/**
 * @class MultipleSchemesStrategy
 *
 * @brief A class representing a strategy for handling multiple search schemes
 * with different edit distances.
 *
 * This class is designed to manage multiple search schemes with varying edit
 * distances. It inherits from the SearchStrategy class and provides methods for
 * calculating the number of parts and creating searches based on a specified
 * edit distance and search ranges.
 * The Multiple Schemes will dynamically choose what scheme to use.
 */
class MultipleSchemesStrategy : public SearchStrategy {
  private:
    std::vector<MultipleSchemes>
        schemesPerED; // A vector of MultipleSchemes objects, one for each edit
                      // distance.

    /**
     * @brief Private helper function to read and initialize multiple search
     * schemes.
     *
     * This function reads and initializes multiple search schemes from a
     * specified folder for various edit distances. If an edit distance is not
     * provided the function stops reading.
     *
     * @param pathToFolder The path to the folder containing the search schemes.
     * @param verbose A flag indicating whether to display verbose output during
     * reading.
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
            std::string dir = pathToFolder + "/" + std::to_string(k);
            struct stat info;
            if (!(stat(dir.c_str(), &info) == 0 && S_ISDIR(info.st_mode))) {
                // no schemes present, create an empty instance
                schemesPerED.emplace_back(MultipleSchemes(k, verbose));
                continue;
            };

            struct stat buffer;
            if (stat((dir + "/scheme1.txt").c_str(), &buffer) != 0) {
                // no schemes present, create an empty instance
                schemesPerED.emplace_back(MultipleSchemes(k, verbose));
                continue;
            }

            schemesPerED.emplace_back(MultipleSchemes(dir, k, verbose));
        }
    }

  public:
    /**
     * @brief Constructor for the MultipleSchemesStrategy class.
     *
     * Initializes the MultipleSchemesStrategy object by reading and
     * initializing search schemes for various edit distances.
     *
     * @param index The FMIndex object used for searching.
     * @param pathToFolder The path to the folder containing the search schemes.
     * @param p The partition strategy to use (default: DYNAMIC).
     * @param metric The distance metric to use (default: EDITOPTIMIZED).
     * @param verbose A flag indicating whether to display verbose output during
     * reading (default: false).
     */
    MultipleSchemesStrategy(FMIndex& index, const std::string& pathToFolder,
                            PartitionStrategy p = DYNAMIC,
                            DistanceMetric metric = EDITOPTIMIZED,
                            bool verbose = false)
        : SearchStrategy(index, p, metric) {

        readSchemes(pathToFolder, verbose);
    }

    /**
     * @brief Calculate the number of parts for a given maximum edit distance.
     *
     * @param maxED The maximum edit distance for which to calculate the number
     * of parts.
     * @return The number of parts (p) corresponding to the specified maximum
     * edit distance.
     */
    uint calculateNumParts(unsigned int maxED) const {
        assert(maxED - 1 < schemesPerED.size());
        assert(!schemesPerED[maxED - 1].isEmpty());
        return schemesPerED[maxED - 1].calculateNumParts(maxED);
    }

    /**
     * @brief Create searches based on a maximum edit distance and specified
     * ranges.
     *
     * This method creates searches based on a given maximum edit distance and a
     * set of SARangePair ranges.
     *
     * @param maxED The maximum edit distance for creating searches.
     * @param ranges A vector of SARangePair objects representing search ranges.
     * @return A vector of Search objects based on the specified parameters.
     */
    const std::vector<Search>&
    createSearches(unsigned int maxED,
                   const std::vector<SARangePair>& ranges) const {
        assert(maxED - 1 < schemesPerED.size());
        assert(!schemesPerED[maxED - 1].isEmpty());
        return schemesPerED[maxED - 1].createSearches(maxED, ranges);
    }
};

// ============================================================================
// CLASS NaiveBackTrackingStrategy
// ============================================================================

// Matches a pattern using the naive backtracking strategy.
class NaiveBackTrackingStrategy : public SearchStrategy {
  private:
    std::vector<Search> searches = {};
    uint calculateNumParts(unsigned int maxED) const {
        return 1;
    }
    const std::vector<Search>&
    createSearches(unsigned int maxED,
                   const std::vector<SARangePair>& ranges) const {
        return searches;
    }

  public:
    virtual std::vector<TextOcc> matchApprox(const std::string& pattern,
                                             length_t maxED, Counters& counters,
                                             const std::string& ID,
                                             const std::string& qual,
                                             bool revComp) const {

        counters.resetCounters();
        if (maxED == 0) {

            auto m = index.exactMatchesOutput(pattern, counters);
            generateSAM(m, pattern, ID, qual, revComp);
            return m;
        }
        auto m = index.approxMatchesNaive(pattern, maxED, counters);
        generateSAM(m, pattern, ID, qual, revComp);
        return m;
    }

    NaiveBackTrackingStrategy(FMIndex& index, PartitionStrategy p = DYNAMIC,
                              DistanceMetric metric = EDITOPTIMIZED)
        : SearchStrategy(index, p, metric) {
        name = "Naive backtracking";
    };
};

// ============================================================================
// HARDCODED CUSTOM CLASSES (TODO remove these)
// ============================================================================

class KucherovKplus1 : public SearchStrategy {
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

    const std::vector<std::vector<int>> weights = {
        {1, 1}, {39, 10, 40}, {400, 4, 5, 400}, {100, 5, 1, 6, 105}};

    const std::vector<std::vector<double>> staticPositions = {
        {0.5}, {0.41, 0.7}, {0.25, 0.50, 0.75}, {0.27, 0.47, 0.62, 0.81}};
    uint calculateNumParts(unsigned int maxED) const {
        return maxED + 1;
    }
    const std::vector<Search>&
    createSearches(unsigned int maxED,
                   const std::vector<SARangePair>& ranges) const {
        assert(maxED >= 1);
        assert(maxED <= 4);

        return schemePerED[maxED - 1];
    }
    const std::vector<double> getBegins(const int& numParts,
                                        const int& maxScore) const override {
        return staticPositions[maxScore - 1];
    }
    const std::vector<int> getWeights(const int& numParts,
                                      const int& maxScore) const override {
        return weights[maxScore - 1];
    }

    const std::vector<double>
    getSeedingPositions(const int& numParts,
                        const int& maxScore) const override {
        return seedingPositions[maxScore - 1];
    }

  public:
    KucherovKplus1(FMIndex& index, PartitionStrategy p = DYNAMIC,
                   DistanceMetric metric = EDITOPTIMIZED)
        : SearchStrategy(index, p, metric) {
        name = "KUCHEROV K + 1";
    };
};

class KucherovKplus2 : public SearchStrategy {
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

    uint calculateNumParts(unsigned int maxED) const {
        return maxED + 2;
    }
    const std::vector<Search>&
    createSearches(unsigned int maxED,
                   const std::vector<SARangePair>& ranges) const {
        return schemePerED[maxED - 1];
    }

    const std::vector<std::vector<double>> seedingPositions = {
        {0.94}, {0.48, 0.55}, {0.4, 0.63, 0.9}, {0.34, 0.5, 0.65, 0.7}};

    const std::vector<std::vector<int>> weights = {{11, 10, 1},
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
    const std::vector<int> getWeights(const int& numParts,
                                      const int& maxScore) const override {
        return weights[maxScore - 1];
    }

    virtual const std::vector<double>
    getSeedingPositions(const int& numParts,
                        const int& maxScore) const override {
        return seedingPositions[maxScore - 1];
    }

  public:
    KucherovKplus2(FMIndex& index, PartitionStrategy p = DYNAMIC,
                   DistanceMetric metric = EDITOPTIMIZED)
        : SearchStrategy(index, p, metric) {
        name = "KUCHEROV K + 2";
    };
};

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

    const std::vector<std::vector<int>> weights = {
        {1, 1}, {10, 1, 5}, {1, 1, 1, 1}, {7, 2, 1, 3, 5}};

    const std::vector<std::vector<double>> staticPositions = {
        {0.5}, {0.30, 0.60}, {0.17, 0.69, 0.96}, {0.2, 0.5, 0.6, 0.8}};
    uint calculateNumParts(unsigned int maxED) const {
        return maxED + 1;
    }
    const std::vector<Search>&
    createSearches(unsigned int maxED,
                   const std::vector<SARangePair>& ranges) const {
        if (maxED < 1 || maxED > 5) {
            throw std::invalid_argument("max ED should be between 1 and 4");
        }
        return schemePerED[maxED - 1];
    }

    const std::vector<double> getBegins(const int& numParts,
                                        const int& maxScore) const override {
        return staticPositions[maxScore - 1];
    }
    const std::vector<int> getWeights(const int& numParts,
                                      const int& maxScore) const override {
        return weights[maxScore - 1];
    }

    const std::vector<double>
    getSeedingPositions(const int& numParts,
                        const int& maxScore) const override {
        return seedingPositions[maxScore - 1];
    }

  public:
    OptimalKianfar(FMIndex& index, PartitionStrategy p = DYNAMIC,
                   DistanceMetric metric = EDITOPTIMIZED)
        : SearchStrategy(index, p, metric) {
        name = "OPTIMAL KIANFAR";
    };
};

// ============================================================================
// CLASS O1StarSearchStrategy
// ============================================================================

// A concrete derived class of SearchStrategy. The strategy here is founded
// on this observation: if x errors are allowed and the pattern is divided
// up in (x
// + 2) parts then every match with max x errors contains a seed consisting
// of n parts, where the first and last part of the seed contain no errors
// and all parts inbetween these contain exacly one error. (2 <= n <= x +
// 2)
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

    uint calculateNumParts(unsigned int maxED) const {
        return maxED + 2;
    }
    const std::vector<Search>&
    createSearches(unsigned int maxED,
                   const std::vector<SARangePair>& ranges) const {
        return schemePerED[maxED - 1];
    }

    const std::vector<std::vector<double>> seedingPositions = {
        {0.94}, {0.51, 0.93}, {0.34, 0.64, 0.88}, {0.28, 0.48, 0.63, 0.94}};

    const std::vector<std::vector<int>> weights = {
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
    const std::vector<int> getWeights(const int& numParts,
                                      const int& maxScore) const override {
        return weights[maxScore - 1];
    }

    virtual const std::vector<double>
    getSeedingPositions(const int& numParts,
                        const int& maxScore) const override {
        return seedingPositions[maxScore - 1];
    }

  public:
    O1StarSearchStrategy(FMIndex& index, PartitionStrategy p = DYNAMIC,
                         DistanceMetric metric = EDITOPTIMIZED)
        : SearchStrategy(index, p, metric) {
        name = "01*0";
    };
};

class ManBestStrategy : public SearchStrategy {
  private:
    const std::vector<Search> ED4 = {
        Search::makeSearch({0, 1, 2, 3, 4, 5}, {0, 0, 0, 0, 0, 4},
                           {0, 3, 3, 3, 4, 4}, 0),
        Search::makeSearch({1, 2, 3, 4, 5, 0}, {0, 0, 0, 0, 0, 0},
                           {0, 2, 2, 3, 3, 4}, 1),
        Search::makeSearch({2, 1, 3, 4, 5, 0}, {0, 1, 1, 1, 1, 1},
                           {0, 2, 2, 3, 3, 4}, 2),
        Search::makeSearch({3, 2, 1, 4, 5, 0}, {0, 1, 2, 2, 2, 2},
                           {0, 1, 2, 3, 3, 4}, 3),
        Search::makeSearch({5, 4, 3, 2, 1, 0}, {0, 0, 0, 0, 3, 3},
                           {0, 0, 4, 4, 4, 4}, 4)};

    uint calculateNumParts(unsigned int maxED) const {
        return maxED + 2;
    }
    const std::vector<Search>&
    createSearches(unsigned int maxED,
                   const std::vector<SARangePair>& ranges) const {
        assert(maxED == 4);
        return ED4;
    }

    const std::vector<double> seedingPositions = {0.35, 0.59, 0.67, 0.9};

    const std::vector<int> weights = {89, 15, 90, 1, 48, 84};

    const std::vector<double> staticPositions = {0.24, 0.43, 0.62, 0.73, 0.77};

    const std::vector<double> getBegins(const int& numParts,
                                        const int& maxScore) const override {
        assert(maxScore == 4);
        return staticPositions;
    }
    const std::vector<int> getWeights(const int& numParts,
                                      const int& maxScore) const override {
        assert(maxScore == 4);
        return weights;
    }

    virtual const std::vector<double>
    getSeedingPositions(const int& numParts,
                        const int& maxScore) const override {
        assert(maxScore == 4);
        return seedingPositions;
    }

  public:
    ManBestStrategy(FMIndex& index, PartitionStrategy p = DYNAMIC,
                    DistanceMetric metric = EDITOPTIMIZED)
        : SearchStrategy(index, p, metric) {
        name = "MANBEST";
    };
};
// ============================================================================
// CLASS PIGEONHOLESEARCHSTRATEGY
// ============================================================================

// A concrete derived class of SearchStrategy. The strategy here is founded
// on this observation: if x errors are allowed and the pattern is divided
// up in (x
// + 1) sections then every approximate match has an exact match with at
// least one of the sections. The strategy iterates over the sections, it
// tries to exactly match the current section, then approximately match the
// pattern before this section and after the the pattern after this section
// with the remaining edit distance.
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
    uint calculateNumParts(unsigned int maxED) const {
        return maxED + 1;
    }
    const std::vector<Search>&
    createSearches(unsigned int maxED,
                   const std::vector<SARangePair>& ranges) const {
        return schemePerED[maxED - 1];
    }

  public:
    PigeonHoleSearchStrategy(FMIndex& index, PartitionStrategy p = DYNAMIC,
                             DistanceMetric metric = EDITOPTIMIZED)
        : SearchStrategy(index, p, metric) {
        name = "PIGEON HOLE";
    };
};

#endif

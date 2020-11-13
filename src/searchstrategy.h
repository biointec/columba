/******************************************************************************
 *  Columba: Approximate Pattern Matching using Search Schemes                *
 *  Copyright (C) 2020 - Luca Renders <luca.renders@ugent.be> and             *
 *                       Jan Fostier <jan.fostier@ugent.be>                   *
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

#include "bidirecfmindex.h"

// ============================================================================
// CLASS SEARCHSTRATEGY
// ============================================================================

// This is an abstract class. Every derived class should be able to create
// searches for a given value of k. This abstract base class handles the
// partitioning (either with values provided in the derived class or default
// uniform values) and approximate matching (either hamming or edit distance)

enum PartitionStrategy { UNIFORM, STATIC, DYNAMIC };

class SearchStrategy;
typedef void (SearchStrategy::*PartitionPtr)(const std::string&,
                                             std::vector<Substring>&);
typedef void (SearchStrategy::*StartIdxPtr)(const Search&, const BiAppMatchSA&,
                                            std::vector<BiAppMatchSA>&,
                                            std::vector<Substring>&,
                                            const int&) const;

class SearchStrategy {
  protected:
    int numParts; // the number of parts the pattern will be split into
    std::vector<Search> searches; // the searches of this strategy

    BidirecFMIndex* index; // pointer to the index of the text that is searched

    std::string name; // the name of the partical search strategy

    std::vector<SARangePair> exactMatchRanges; // ranges corresponding to
                                               // exact matches of the parts

    PartitionPtr partitionPtr; // pointer to the partition method
    StartIdxPtr
        startIdxPtr;    // pointer to start method (hamming or edit distance)
    int lastMaxED = -1; // the last maxED so that searches do not need to be
                        // recalculated for next read

    PartitionStrategy partitionStrategy;
    bool isEdit;

    /**
     * Constructor
     * @param argument, pointer to the bidirectional FM index to use
     * @param p, partition strategy
     * @param edit, true if edit distance should be used, false if hamming
     * distance should be used
     */
    SearchStrategy(BidirecFMIndex* argument, PartitionStrategy p, bool edit) {
        index = argument;
        partitionStrategy = p;
        isEdit = edit;
        switch (p) {
        case UNIFORM:
            partitionPtr = &SearchStrategy::partitionUniform;
            break;
        case DYNAMIC:
            partitionPtr = &SearchStrategy::partitionUniformRange;
            break;
        case STATIC:
            partitionPtr = &SearchStrategy::partitionOptimalStatic;
            break;
        default:
            break;
        }

        if (edit) {
            startIdxPtr = &SearchStrategy::startIndexEdit;
        } else {
            startIdxPtr = &SearchStrategy::startIndexHamming;
        }
    }

    /**
     * Calculates the number of parts for a certain max edit distance. This
     * calculation is strategy dependent
     * @param maxED the maximal allowed edit distance for the alligning
     */
    virtual void calculateNumParts(unsigned int maxED) = 0;

    /**
     * Creates all searches for this specific strategy
     * @param maxED the maximal allowed edit distance for the aligning
     */
    virtual void createSearches(unsigned int maxED) = 0;

    /**
     * Splits the pattern into numParts parts, either by uniform range or
     * uniform size
     * @param pattern the pattern to be split
     * @param parts the vector containing the substrings of this pattern,
     * will be cleared and filled during the execution of this method. If
     * the splitting fails for some reason, the vector will be empty
     * @param maxED, the maximum allowed edit distance
     */
    void partition(const std::string& pattern, std::vector<Substring>& parts,
                   int maxED) {

        parts.clear();

        if (lastMaxED != maxED) {

            // calculate how many parts there will be
            calculateNumParts(maxED);
            // create the searches
            createSearches(maxED);
            lastMaxED = maxED;
        }

        if (numParts >= (int)pattern.size() || numParts == 1) {
            // no need of splitting up since all parts would be one character
            // or less or there is only one part
            return;
        }

        (this->*partitionPtr)(pattern, parts);
    }
    /**
     * Helper function for optimal static partitioning. This function creates
     * the optimal static parts
     * @param pattern, the pattern to partition
     * @param parts, empty vector to which the parts are added
     */
    void setParts(const std::string& pattern, std::vector<Substring>& parts) {
        const std::vector<double>& begins = getBegins();

        int pSize = pattern.size();
        // set the first part
        parts.emplace_back(pattern, 0, begins[0] * pSize);

        for (unsigned int i = 0; i < begins.size() - 1; i++) {
            parts.emplace_back(pattern, begins[i] * pSize,
                               begins[i + 1] * pSize);
        }
        parts.emplace_back(pattern, begins.back() * pSize, pattern.size());
    }

    /**
     * Splits the pattern into numParts parts, such that each part has the same
     * size
     * @param pattern the pattern to be split
     * @param parts an empty vector which will be filled with the different
     * parts
     * @param maxED, the maximum allowed edit distance
     */
    void partitionUniform(const std::string& pattern,
                          std::vector<Substring>& parts) {

        for (int i = 0; i < numParts; i++) {
            parts.emplace_back(pattern, i * pattern.size() / numParts,
                               (i + 1) * pattern.size() / numParts);
        }
        // set end of final part correct
        parts.back().setEnd(pattern.size());

        // match the exactRanges for each part
        index->setDirection(FORWARD);
        SARangePair initialRanges = index->getCompleteRange();

        exactMatchRanges.resize(numParts);
        std::vector<bool> partNumberSeen(numParts, false);

        for (int i = 0; i < numParts; i++) {
            exactMatchRanges[i] =
                index->matchStringBidirectionally(parts[i], initialRanges);
        }
    }

    /**
     * Splits the pattern into numParts parts, such that each search carries the
     * same weight (on average)
     * @param pattern the pattern to be split
     * @param parts an empty vector which will be filled with the different
     * parts
     * @param maxED, the maximum allowed edit distance
     */
    void partitionOptimalStatic(const std::string& pattern,
                                std::vector<Substring>& parts) {

        setParts(pattern, parts);

        // match the exactRanges for each part
        index->setDirection(FORWARD);
        SARangePair initialRanges = index->getCompleteRange();

        exactMatchRanges.resize(numParts);
        std::vector<bool> partNumberSeen(numParts, false);

        for (int i = 0; i < numParts; i++) {
            exactMatchRanges[i] =
                index->matchStringBidirectionally(parts[i], initialRanges);
        }
    }

    /**
     * Function that retrieves the seeding positions for dynamic partitioning.
     * If derived class does not implement this function then uniform seeds are
     * given.
     * @returns vector with doubles indicating the position (relative to the
     * length of the pattern) where a seed should be placed.
     */
    virtual std::vector<double> getSeedingPositions() const {

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
     * @returns the number of characters used by the seeding operation
     */
    int seed(const std::string& pattern, std::vector<Substring>& parts) {
        int pSize = pattern.size();
        int wSize = index->getWordSize();

        const auto& seedPercent = getSeedingPositions();

        std::vector<int> seeds;
        // push the seed for the first part
        seeds.push_back(0);

        // push the optimal seeds for the middle parts
        for (int i = 1; i < numParts - 1; i++) {
            seeds.push_back((seedPercent[i - 1] * pSize) - (wSize / 2));
        }

        for (int i = 0; i < numParts - 1; i++) {
            parts.emplace_back(pattern, seeds[i], seeds[i] + wSize);
        }

        // push the seeds for the final parts
        parts.emplace_back(pattern, pSize - wSize, pSize);

        exactMatchRanges.resize(numParts);
        for (int i = 0; i < numParts; i++) {
            exactMatchRanges[i] = index->lookUpInKmerTable(parts[i]);
        }
        return numParts * wSize;
    }

    /**
     * Function that retrieves the weights for dynamic partitioning.
     * If derived class does not implement this function then uniform weights
     * are given.
     * @returns vector with weights
     */
    virtual std::vector<int> getWeights() const {
        std::vector<int> w(numParts, 1);
        return w;
    }

    /**
     * Splits the pattern into numParts parts, such that each part has
     * (approximately) the same range. The exactMatchRanges are also
     * calculated.
     * @param pattern the pattern to be split
     * @param parts an empty vector which will be filled with the different
     * parts
     * @param maxED, the maximum allowed edit distance
     */
    void partitionUniformRange(const std::string& pattern,
                               std::vector<Substring>& parts) {

        int matchedChars = seed(pattern, parts);
        int pSize = pattern.size();
        std::vector<int> weights = getWeights();

        Direction dir = FORWARD;
        int partToExtend = 0;

        // extend the part with the largest range, as to minimize the range for
        // each part do this untill all characters are assigned to a part
        for (int j = matchedChars; j < pSize; j++) {

            // find the part with the largest range
            length_t maxRange = 0;
            for (int i = 0; i < numParts; i++) {
                bool noLeftExtension =
                    (i == 0) || parts[i].begin() == parts[i - 1].end();
                bool noRightExtension = (i == numParts - 1) ||
                                        parts[i].end() == parts[i + 1].begin();
                if (noLeftExtension && noRightExtension) {
                    continue;
                }
                if (exactMatchRanges[i].width() * weights[i] >= maxRange) {
                    maxRange = exactMatchRanges[i].width() * weights[i];
                    partToExtend = i;
                    if (noLeftExtension) {
                        // only right extension
                        dir = FORWARD;
                    } else if (noRightExtension) {
                        // only left extension
                        dir = BACKWARD;
                    } else {
                        // both directions possible, choose direction of
                        // smallest neighbour
                        dir = (exactMatchRanges[i - 1].width() <
                               exactMatchRanges[i + 1].width())
                                  ? BACKWARD
                                  : FORWARD;
                    }
                }
            }

            if (maxRange == 0) {
                // no need to keep calculating new range, just extend the parts
                extendParts(pattern, parts);
                return;
            }

            // extend partToExtend in direction
            char c; // the new character
            if (dir == FORWARD) {
                parts[partToExtend].incrementEnd();
                c = pattern[parts[partToExtend].end() - 1];
            } else {
                parts[partToExtend].decrementBegin();
                c = pattern[parts[partToExtend].begin()];
            }

            // match the new character
            index->setDirection(dir);
            index->addChar(c, exactMatchRanges.at(partToExtend));
        }
    }

    /**
     * Helper function for dynamic partitioning. This function extends the parts
     * so that nothing of the pattern is not allocated to any part. This does
     * not keep track of the ranges over the suffix array, so should only be
     * called if this does not matter (e.g. when the parts that can be extended
     * all correspond to empty ranges)
     */
    void extendParts(const std::string& pattern,
                     std::vector<Substring>& parts) const {
        for (length_t i = 0; i < parts.size(); i++) {

            if ((i != (length_t)numParts - 1) &&
                (parts[i].end() != parts[i + 1].begin())) {
                // extend completely to the right
                // it is known that the range will stay [0,0)
                parts[i].setEnd(parts[i + 1].begin());
            }
            if ((i != 0) && (parts[i].begin() != parts[i - 1].end())) {
                // extend completly to the left
                parts[i].setBegin(parts[i - 1].end());
            }
        }
    }

    /**
     * Function that retrieves the begin postiiions for optimal static
     * partitioning. If derived class does not implement this function then
     * uniform positions are given.
     * @returns vector with doubles indicating the position (relative to the
     * length of the pattern) where a seed should be placed.
     */
    virtual const std::vector<double> getBegins() const {
        std::vector<double> b;
        double u = 1.0 / numParts;
        for (int i = 1; i < numParts; i++) {
            b.push_back(i * u);
        }
        return b;
    }

    /**
     * Starts the index with hamming distance
     * @param s, the search to follow
     * @param startMatch, the startMatch that corresponds to the first piece
     * @param occ, vector to add occurrences to
     * @param parts, the parts of the pattern
     * @param idx, the index in the search to match next
     */
    void startIndexHamming(const Search& s, const BiAppMatchSA& startMatch,
                           std::vector<BiAppMatchSA>& occ,
                           std::vector<Substring>& parts,
                           const int& idx) const {
        index->recApproxMatchHamming(s, startMatch, occ, parts, idx);
    }

    /**
     * Starts the index with edit distance
     * @param s, the search to follow
     * @param startMatch, the startMatch that corresponds to the first piece
     * @param occ, vector to add occurrences to
     * @param parts, the parts of the pattern
     * @param idx, the index in the search to match next
     */

    void startIndexEdit(const Search& s, const BiAppMatchSA& startMatch,
                        std::vector<BiAppMatchSA>& occ,
                        std::vector<Substring>& parts, const int& idx) const {
        index->recApproxMatch(s, startMatch, occ, parts, idx);
    }
    /**
     * Executes the search recursively. If U[0] != 1, then the search will start
     * at pi[0], else the search will start with idx i and U[i]!=0 and U[j]=0
     * with j < i
     * @param s, the search to follow
     * @param parts, the parts of the pattern
     * @param allMatches, vector to add occurrences to
     */
    void doRecSearch(const Search& s, std::vector<Substring>& parts,
                     std::vector<BiAppMatchSA>& allMatches) const {

        if (s.upperBounds[0] > 0) {
            // first part is allowed an error so start with an empty match
            setPartsDirections(s, parts);
            SARangePair startRange = index->getCompleteRange();
            BiAppMatchSA startMatch = makeBiAppMatchSA(startRange, 0, 0);
            (this->*startIdxPtr)(s, startMatch, allMatches, parts, 0);
            return;
        }

        // first get the bidirectional match of first part
        int first = s.order[0];
        SARangePair startRange = exactMatchRanges[first];

        if (!startRange.empty()) {

            setPartsDirections(s, parts);

            int partInSearch = 1;
            length_t exactLength = parts[first].size();

            while (s.upperBounds[partInSearch] == 0) {
                // extend the exact match
                index->setDirection(s.directions[partInSearch - 1]);
                startRange = index->matchStringBidirectionally(
                    parts[s.order[partInSearch]], startRange);
                if (startRange.empty()) {
                    return;
                }
                exactLength += parts[s.order[partInSearch]].size();
                partInSearch++;
            }
            BiAppMatchSA startMatch =
                makeBiAppMatchSA(startRange, 0, exactLength);

            (this->*startIdxPtr)(s, startMatch, allMatches, parts,
                                 partInSearch);
        }
    }

  public:
    /**
     * Retrieves the name of this strategy, derived classes should set a
     * meaningful name
     */
    std::string getName() const {
        return name;
    }

    std::string getPartitioningStrategy() const {
        switch (partitionStrategy) {
        case UNIFORM:
            return "UNIFORM";
            break;
        case DYNAMIC:
            return "DYNAMIC";
            break;
        case STATIC:
            return "STATIC";
            break;

        default:
            // should not get here
            return "";
        }
    }

    std::string getEditOrHamming() const {
        return isEdit ? "EDIT" : "HAMMING";
    }

    /**
     * Mathces a pattern approximately using this strategy
     * @param pattern, the pattern to match
     * @param maxED, the maximal allowed edit distance (or  hamming distance)
     */
    std::vector<AppMatch> matchApprox(const std::string& pattern,
                                      length_t maxED) {
        index->resetCounters();

        if (maxED == 0) {
            index->setDirection(BACKWARD);
            auto result = index->exactMatches(pattern);
            std::vector<AppMatch> returnvalue;
            for (length_t startpos : result) {
                AppMatch m;
                m.editDist = 0;
                m.range = Range(startpos, startpos + pattern.size());
                returnvalue.push_back(m);
            }
            return returnvalue;
        }
        // create the parts of the pattern
        std::vector<Substring> parts;

        partition(pattern, parts, maxED);

        if (parts.empty() || numParts * maxED >= pattern.size()) {
            // splitting up was not viable just search the entire pattern
            std::cerr << "Warning: Normal bidirectional search was used as "
                         "entered pattern is too short"
                      << std::endl;

            return index->approxMatches(pattern, maxED);
        }

        // the vector containing all matches in the sufffix array
        std::vector<BiAppMatchSA> allMatches;

        index->reserveStacks(numParts, pattern.length());
        // do all searches

        for (const Search& s : searches) {

            doRecSearch(s, parts, allMatches);
        }

        // return all matches mapped to the text
        return index->mapOccurencesInSAToOccurencesInText(allMatches, maxED);
    }

    virtual ~SearchStrategy() {
    }
};

class KucherovKplus1 : public SearchStrategy {
  private:
    const std::vector<Search> ED1 = {makeSearch({0, 1}, {0, 0}, {0, 1}),
                                     makeSearch({1, 0}, {0, 0}, {0, 1})};
    const std::vector<Search> ED2 = {
        makeSearch({0, 1, 2}, {0, 0, 0}, {0, 2, 2}),
        makeSearch({2, 1, 0}, {0, 0, 0}, {0, 1, 2}),
        makeSearch({1, 0, 2}, {0, 0, 1}, {0, 1, 2})};

    const std::vector<Search> ED3 = {
        makeSearch({0, 1, 2, 3}, {0, 0, 0, 0}, {0, 1, 3, 3}),
        makeSearch({1, 0, 2, 3}, {0, 0, 1, 1}, {0, 1, 3, 3}),
        makeSearch({2, 3, 1, 0}, {0, 0, 0, 0}, {0, 1, 3, 3}),
        makeSearch({3, 2, 1, 0}, {0, 0, 1, 1}, {0, 1, 3, 3})};

    const std::vector<Search> ED4 = {
        makeSearch({0, 1, 2, 3, 4}, {0, 0, 0, 0, 0}, {0, 2, 2, 4, 4}),
        makeSearch({4, 3, 2, 1, 0}, {0, 0, 0, 0, 0}, {0, 1, 3, 4, 4}),
        makeSearch({1, 0, 2, 3, 4}, {0, 0, 1, 3, 3}, {0, 1, 3, 4, 4}),
        makeSearch({0, 1, 2, 3, 4}, {0, 0, 1, 3, 3}, {0, 1, 3, 4, 4}),
        makeSearch({3, 2, 4, 1, 0}, {0, 0, 0, 1, 1}, {0, 1, 2, 4, 4}),
        makeSearch({2, 1, 0, 3, 4}, {0, 0, 0, 1, 3}, {0, 1, 2, 4, 4}),
        makeSearch({1, 0, 2, 3, 4}, {0, 0, 1, 2, 4}, {0, 1, 2, 4, 4}),
        makeSearch({0, 1, 2, 3, 4}, {0, 0, 0, 3, 4}, {0, 0, 4, 4, 4})};

    const std::vector<std::vector<Search>> schemePerED = {ED1, ED2, ED3, ED4};

    const std::vector<std::vector<double>> seedingPositions = {
        {}, {0.59}, {0.34, 0.66}, {0.42, 0.56, 0.67}};

    const std::vector<std::vector<int>> weights = {
        {1, 1}, {2, 1, 2}, {1, 1, 1, 1}, {7, 2, 1, 3, 5}};

    const std::vector<std::vector<double>> staticPositions = {
        {0.5}, {0.40, 0.66}, {0.26, 0.50, 0.75}, {0.26, 0.46, 0.61, 0.80}};
    void calculateNumParts(unsigned int maxED) {
        numParts = maxED + 1;
    }
    void createSearches(unsigned int maxED) {
        assert(maxED >= 1);
        assert(maxED <= 4);

        searches = schemePerED[maxED - 1];
    }
    const std::vector<double> getBegins() const override {
        return staticPositions[numParts - 2];
    }
    std::vector<int> getWeights() const override {
        return weights[numParts - 2];
    }

    virtual std::vector<double> getSeedingPositions() const override {
        return seedingPositions[numParts - 2];
    }

  public:
    KucherovKplus1(BidirecFMIndex* index, PartitionStrategy p = DYNAMIC,
                   bool edit = true)
        : SearchStrategy(index, p, edit) {
        name = "KUCHEROV K + 1";
    };
};

class KucherovKplus2 : public SearchStrategy {
  private:
    const std::vector<Search> ED1 = {
        makeSearch({0, 1, 2}, {0, 0, 0}, {0, 1, 1}),
        makeSearch({1, 2, 0}, {0, 0, 0}, {0, 0, 1})};
    const std::vector<Search> ED2 = {
        makeSearch({0, 1, 2, 3}, {0, 0, 0, 0}, {0, 1, 1, 2}),
        makeSearch({3, 2, 1, 0}, {0, 0, 0, 0}, {0, 1, 2, 2}),
        makeSearch({1, 2, 3, 0}, {0, 0, 0, 1}, {0, 0, 1, 2}),
        makeSearch({0, 1, 2, 3}, {0, 0, 0, 2}, {0, 0, 2, 2})};

    const std::vector<Search> ED3 = {
        makeSearch({0, 1, 2, 3, 4}, {0, 0, 0, 0, 0}, {0, 1, 2, 3, 3}),
        makeSearch({1, 2, 3, 4, 0}, {0, 0, 0, 0, 0}, {0, 1, 2, 2, 3}),
        makeSearch({2, 3, 4, 1, 0}, {0, 0, 0, 0, 1}, {0, 1, 1, 3, 3}),
        makeSearch({3, 4, 2, 1, 0}, {0, 0, 0, 1, 2}, {0, 0, 3, 3, 3})};

    const std::vector<Search> ED4 = {
        makeSearch({0, 1, 2, 3, 4, 5}, {0, 0, 0, 0, 0, 0}, {0, 1, 2, 3, 4, 4}),
        makeSearch({1, 2, 3, 4, 5, 0}, {0, 0, 0, 0, 0, 0}, {0, 1, 2, 3, 4, 4}),
        makeSearch({5, 4, 3, 2, 1, 0}, {0, 0, 0, 0, 0, 1}, {0, 1, 2, 2, 4, 4}),
        makeSearch({3, 4, 5, 2, 1, 0}, {0, 0, 0, 0, 1, 2}, {0, 1, 1, 3, 4, 4}),
        makeSearch({2, 3, 4, 5, 1, 0}, {0, 0, 0, 0, 2, 3}, {0, 1, 1, 2, 4, 4}),
        makeSearch({4, 5, 3, 2, 1, 0}, {0, 0, 0, 1, 3, 3}, {0, 0, 3, 3, 4, 4}),
        makeSearch({0, 1, 2, 3, 4, 5}, {0, 0, 0, 3, 3, 3}, {0, 0, 3, 3, 4, 4}),
        makeSearch({0, 1, 2, 3, 4, 5}, {0, 0, 0, 0, 4, 4}, {0, 0, 2, 4, 4, 4}),
        makeSearch({2, 3, 1, 0, 4, 5}, {0, 0, 0, 1, 2, 4}, {0, 0, 2, 2, 4, 4}),
        makeSearch({4, 5, 3, 2, 1, 0}, {0, 0, 0, 0, 4, 4}, {0, 0, 1, 4, 4, 4})};

    const std::vector<std::vector<Search>> schemePerED = {ED1, ED2, ED3, ED4};

    void calculateNumParts(unsigned int maxED) {
        numParts = maxED + 2;
    }
    void createSearches(unsigned int maxED) {
        searches = schemePerED[maxED - 1];
    }

    const std::vector<std::vector<double>> seedingPositions = {
        {0.94}, {0.5, 0.55}, {0.4, 0.63, 0.9}, {0.35, 0.5, 0.65, 0.69}};

    const std::vector<std::vector<int>> weights = {
        {11, 10, 1}, {8, 3, 1, 8}, {6, 3, 2, 1, 1}, {26, 21, 8, 2, 7, 21}};

    const std::vector<std::vector<double>> staticPositions = {
        {0.33, 0.66},
        {0.25, 0.50, 0.75},
        {0.2, 0.4, 0.6, 0.8},
        {0.16, 0.33, 0.49, 0.66, 0.82}};

    const std::vector<double> getBegins() const override {
        return staticPositions[numParts - 3];
    }
    std::vector<int> getWeights() const override {
        return weights[numParts - 3];
    }

    virtual std::vector<double> getSeedingPositions() const override {
        return seedingPositions[numParts - 3];
    }

  public:
    KucherovKplus2(BidirecFMIndex* index, PartitionStrategy p = DYNAMIC,
                   bool edit = true)
        : SearchStrategy(index, p, true) {
        name = "KUCHEROV K + 2";
    };
};

class OptimalKianfar : public SearchStrategy {
  private:
    const std::vector<Search> ED1 = {makeSearch({0, 1}, {0, 0}, {0, 1}),
                                     makeSearch({1, 0}, {0, 1}, {0, 1})};
    const std::vector<Search> ED2 = {
        makeSearch({0, 1, 2}, {0, 0, 2}, {0, 1, 2}),
        makeSearch({2, 1, 0}, {0, 0, 0}, {0, 2, 2}),
        makeSearch({1, 2, 0}, {0, 1, 1}, {0, 1, 2})};

    const std::vector<Search> ED3 = {
        makeSearch({0, 1, 2, 3}, {0, 0, 0, 3}, {0, 2, 3, 3}),
        makeSearch({1, 2, 3, 0}, {0, 0, 0, 0}, {1, 2, 3, 3}),
        makeSearch({2, 3, 1, 0}, {0, 0, 2, 2}, {0, 0, 3, 3})};

    const std::vector<Search> ED4 = {
        makeSearch({0, 1, 2, 3, 4}, {0, 0, 0, 0, 4}, {0, 3, 3, 4, 4}),
        makeSearch({1, 2, 3, 4, 0}, {0, 0, 0, 0, 0}, {2, 2, 3, 3, 4}),
        makeSearch({4, 3, 2, 1, 0}, {0, 0, 0, 3, 3}, {0, 0, 4, 4, 4})};

    const std::vector<std::vector<Search>> schemePerED = {ED1, ED2, ED3, ED4};
    const std::vector<std::vector<double>> seedingPositions = {
        {}, {0.59}, {0.34, 0.66}, {0.42, 0.56, 0.67}};

    const std::vector<std::vector<int>> weights = {
        {1, 1}, {2, 1, 2}, {1, 1, 1, 1}, {7, 2, 1, 3, 5}};

    const std::vector<std::vector<double>> staticPositions = {
        {0.5}, {0.33, 0.66}, {0.25, 0.50, 0.75}, {0.2, 0.4, 0.6, 0.80}};
    void calculateNumParts(unsigned int maxED) {
        numParts = maxED + 1;
    }
    void createSearches(unsigned int maxED) {
        if (maxED < 1 || maxED > 4) {
            throw std::invalid_argument("max ED should be between 1 and 4");
        }
        searches = schemePerED[maxED - 1];
    }

    const std::vector<double> getBegins() const override {
        return staticPositions[numParts - 2];
    }
    std::vector<int> getWeights() const override {
        return weights[numParts - 2];
    }

    std::vector<double> getSeedingPositions() const override {
        return seedingPositions[numParts - 2];
    }

  public:
    OptimalKianfar(BidirecFMIndex* index, PartitionStrategy p = DYNAMIC,
                   bool edit = true)
        : SearchStrategy(index, p, edit) {
        name = "OPTIMAL KIANFAR";
    };
};

// ============================================================================
// CLASS O1StarSearchStrategy
// ============================================================================

// A concrete derived class of SearchStrategy. The strategy here is founded
// on this observation: if x errors are allowed and the pattern is divided
// up in (x
// + 2) parts then every match with max x erros contains a seed consisting
// of n parts, where the first and last part of the seed contain no errors
// and all parts inbetween these contain exacly one error. (2 <= n <= x +
// 2)

class O1StarSearchStrategy : public SearchStrategy {
  private:
    const std::vector<Search> ED1 = {
        makeSearch({0, 1, 2}, {0, 0, 0}, {0, 1, 1}),
        makeSearch({1, 2, 0}, {0, 0, 0}, {0, 0, 1})};
    const std::vector<Search> ED2 = {
        makeSearch({0, 1, 2, 3}, {0, 0, 0, 0}, {0, 1, 2, 2}),
        makeSearch({1, 2, 3, 0}, {0, 0, 0, 0}, {0, 1, 2, 2}),
        makeSearch({2, 3, 1, 0}, {0, 0, 0, 0}, {0, 0, 2, 2})};

    const std::vector<Search> ED3 = {
        makeSearch({0, 1, 2, 3, 4}, {0, 0, 0, 0, 0}, {0, 1, 3, 3, 3}),
        makeSearch({1, 2, 3, 4, 0}, {0, 0, 0, 0, 0}, {0, 1, 3, 3, 3}),
        makeSearch({2, 3, 4, 1, 0}, {0, 0, 0, 0, 0}, {0, 1, 3, 3, 3}),
        makeSearch({3, 4, 2, 1, 0}, {0, 0, 0, 0, 0}, {0, 0, 3, 3, 3})};

    const std::vector<Search> ED4 = {
        makeSearch({0, 1, 2, 3, 4, 5}, {0, 0, 0, 0, 0, 0}, {0, 1, 4, 4, 4, 4}),
        makeSearch({1, 2, 3, 4, 5, 0}, {0, 0, 0, 0, 0, 0}, {0, 1, 4, 4, 4, 4}),
        makeSearch({2, 3, 4, 5, 1, 0}, {0, 0, 0, 0, 0, 0}, {0, 1, 4, 4, 4, 4}),
        makeSearch({3, 4, 5, 2, 1, 0}, {0, 0, 0, 0, 0, 0}, {0, 1, 4, 4, 4, 4}),
        makeSearch({4, 5, 3, 2, 1, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 4, 4, 4, 4}),
    };

    const std::vector<std::vector<Search>> schemePerED = {ED1, ED2, ED3, ED4};

    void calculateNumParts(unsigned int maxED) {
        numParts = maxED + 2;
    }
    void createSearches(unsigned int maxED) {
        searches = schemePerED[maxED - 1];
    }

    const std::vector<std::vector<double>> seedingPositions = {
        {0.94}, {0.5, 0.75}, {0.33, 0.66, 0.88}, {0.29, 0.45, 0.62, 0.9}};

    const std::vector<std::vector<int>> weights = {
        {11, 10, 1}, {1, 1, 1, 1}, {1, 1, 1, 1, 1}, {1, 2, 2, 2, 2, 1}};

    const std::vector<std::vector<double>> staticPositions = {
        {0.50, 0.96},
        {0.32, 0.69, 0.9},
        {0.25, 0.5, 0.75, 0.96},
        {0.16, 0.33, 0.49, 0.66, 0.82}};

    const std::vector<double> getBegins() const override {
        return staticPositions[numParts - 3];
    }
    std::vector<int> getWeights() const override {
        return weights[numParts - 3];
    }

    virtual std::vector<double> getSeedingPositions() const override {
        return seedingPositions[numParts - 3];
    }

  public:
    O1StarSearchStrategy(BidirecFMIndex* index, PartitionStrategy p = DYNAMIC,
                         bool edit = true)
        : SearchStrategy(index, p, edit) {
        name = "01*0";
    };
};

class ManBestStrategy : public SearchStrategy {
  private:
    const std::vector<Search> ED4 = {
        makeSearch({0, 1, 2, 3, 4, 5}, {0, 0, 0, 0, 0, 4}, {0, 3, 3, 3, 4, 4}),
        makeSearch({1, 2, 3, 4, 5, 0}, {0, 0, 0, 0, 0, 0}, {0, 2, 2, 3, 3, 4}),
        makeSearch({2, 1, 3, 4, 5, 0}, {0, 1, 1, 1, 1, 1}, {0, 2, 2, 3, 3, 4}),
        makeSearch({3, 2, 1, 4, 5, 0}, {0, 1, 2, 2, 2, 2}, {0, 1, 2, 3, 3, 4}),
        makeSearch({5, 4, 3, 2, 1, 0}, {0, 0, 0, 0, 3, 3}, {0, 0, 4, 4, 4, 4})};

    void calculateNumParts(unsigned int maxED) {
        numParts = maxED + 2;
    }
    void createSearches(unsigned int maxED) {
        assert(maxED == 4);
        searches = ED4;
    }

    const std::vector<double> seedingPositions = {0.39, 0.6, 0.68, 0.9};

    const std::vector<int> weights = {15, 4, 3, 1, 2, 1};

    const std::vector<double> staticPositions = {0.26, 0.48, 0.66, 0.76, 0.96};

    const std::vector<double> getBegins() const override {
        return staticPositions;
    }
    std::vector<int> getWeights() const override {
        return weights;
    }

    virtual std::vector<double> getSeedingPositions() const override {
        return seedingPositions;
    }

  public:
    ManBestStrategy(BidirecFMIndex* index, PartitionStrategy p = DYNAMIC,
                    bool edit = true)
        : SearchStrategy(index, p, edit) {
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
    const std::vector<Search> ED1 = {makeSearch({0, 1}, {0, 0}, {0, 1}),
                                     makeSearch({1, 0}, {0, 0}, {0, 1})};
    const std::vector<Search> ED2 = {
        makeSearch({0, 1, 2}, {0, 0, 0}, {0, 2, 2}),
        makeSearch({1, 2, 0}, {0, 0, 0}, {0, 2, 2}),
        makeSearch({2, 1, 0}, {0, 0, 0}, {0, 2, 2})};

    const std::vector<Search> ED3 = {
        makeSearch({0, 1, 2, 3}, {0, 0, 0, 0}, {0, 3, 3, 3}),
        makeSearch({1, 0, 2, 3}, {0, 0, 0, 0}, {0, 3, 3, 3}),
        makeSearch({2, 3, 1, 0}, {0, 0, 0, 0}, {0, 3, 3, 3}),
        makeSearch({3, 2, 1, 0}, {0, 0, 0, 0}, {0, 3, 3, 3})};

    const std::vector<Search> ED4 = {
        makeSearch({0, 1, 2, 3, 4}, {0, 0, 0, 0, 0}, {0, 4, 4, 4, 4}),
        makeSearch({1, 2, 3, 4, 0}, {0, 0, 0, 0, 0}, {0, 4, 4, 4, 4}),
        makeSearch({2, 3, 4, 1, 0}, {0, 0, 0, 0, 0}, {0, 4, 4, 4, 4}),
        makeSearch({3, 4, 2, 1, 0}, {0, 0, 0, 0, 0}, {0, 4, 4, 4, 4}),
        makeSearch({4, 3, 2, 1, 0}, {0, 0, 0, 0, 0}, {0, 4, 4, 4, 4})};

    const std::vector<std::vector<Search>> schemePerED = {ED1, ED2, ED3, ED4};
    void calculateNumParts(unsigned int maxED) {
        numParts = maxED + 1;
    }
    void createSearches(unsigned int maxED) {
        searches = schemePerED[maxED - 1];
    }

  public:
    PigeonHoleSearchStrategy(BidirecFMIndex* index,
                             PartitionStrategy p = DYNAMIC, bool edit = true)
        : SearchStrategy(index, p, edit) {
        name = "PIGEON HOLE";
    };
};

class BackTrackStrategyNaive : public SearchStrategy {
  private:
    FMIndex* bwt;
    void calculateNumParts(unsigned int maxED) {
        numParts = 1;
    }
    void createSearches(unsigned int maxED) {
    }

  public:
    virtual std::vector<AppMatch> matchApprox(const std::string& pattern,
                                              length_t maxED) {

        bwt->resetCounters();
        if (maxED == 0) {

            auto result = bwt->exactMatches(pattern);
            std::vector<AppMatch> returnvalue;
            for (length_t startpos : result) {
                AppMatch m;
                m.editDist = 0;
                m.range = Range(startpos, startpos + pattern.size());
                returnvalue.push_back(m);
            }
            return returnvalue;
        }

        return bwt->approxMatches(pattern, maxED);
    }

    BackTrackStrategyNaive(BidirecFMIndex* index, PartitionStrategy p = DYNAMIC,
                           bool edit = true)
        : SearchStrategy(index, p, edit) {
        name = "Naive backtracking";
        bwt = (FMIndex*)index;
    };
};

#endif

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

#ifndef SEARCH_H
#define SEARCH_H

#include "definitions.h" // for length_t, Direction, BACKWARD, FORWARD
#include "substring.h"   // for Substring

#include <algorithm> // for copy, max, any_of, min_element
#include <cassert>   // for assert
#include <cmath>     // for pow
#include <cstdint>   // for uint16_t
#include <fstream>   // for basic_istream, stringstream, ifstream, ostream
#include <iterator>  // for distance
#include <memory>    // for allocator, allocator_traits<>::value_type
#include <numeric>   // for accumulate
#include <sstream>   // for stringstream
#include <stdexcept> // for runtime_error
#include <string>    // for string, operator+, char_traits, basic_string
#include <utility>   // for pair, make_pair
#include <vector>    // for vector, _Bit_iterator

#define Distribution std::vector<int>

// ============================================================================
// CLASS SEARCH
// ============================================================================
/**
 * @class Search
 *
 * @brief A class representing a search  for approximate string
 * matching.
 *
 * This class defines a search with lower and upper bounds and a fixed order of
 * parts,
 */
class Search {
  protected:
    std::vector<length_t> L;     // the vector with lower bounds
    std::vector<length_t> U;     // the vector with upper bounds
    std::vector<length_t> order; // the vector with the order of the parts
    length_t sIdx;               // the index of the search

  private:
    std::vector<Direction> directions; // the directions of each phase
    std::vector<bool>
        directionSwitch; // has the direction switched for each phase

    std::vector<std::pair<length_t, length_t>>
        lowestAndHighestPartsProcessedBefore; // The lowest and highest parts
                                              // processed before each phase
    bool unidirectionalBackwards = false;     // flag to indicate unidirectional
                                              // backwards search
    length_t uniDirectionalBackwardsIndex; // the index from which the search is
                                           // only bidirectional backwards

    /**
     * @brief Constructor for the Search class.
     *
     * Creates a Search object with the provided parameters.
     *
     * @param order The order of parts in the search.
     * @param lowerBounds The lower bounds for each part.
     * @param upperBounds The upper bounds for each part.
     * @param directions The directions of each phase.
     * @param dSwitch Indicates whether the direction switches for each phase.
     * @param lowestAndHighestPartsProcessedBefore The lowest and highest parts
     * processed before each phase.
     * @param idx The unique index of the search.
     */
    Search(std::vector<length_t>& order, std::vector<length_t>& lowerBounds,
           std::vector<length_t>& upperBounds,
           std::vector<Direction>& directions, std::vector<bool>& dSwitch,
           std::vector<std::pair<length_t, length_t>>
               lowestAndHighestPartsProcessedBefore,
           length_t idx, bool uniBackwards, length_t uniBackwardsIndex)
        : L(lowerBounds), U(upperBounds), order(order), sIdx(idx),
          directions(directions), directionSwitch(dSwitch),
          lowestAndHighestPartsProcessedBefore(
              lowestAndHighestPartsProcessedBefore),
          unidirectionalBackwards(uniBackwards),
          uniDirectionalBackwardsIndex(uniBackwardsIndex) {
    }

  public:
    /**
     * @brief Static method to create a Search object.
     *
     * Constructs a Search object with the provided order, lower bounds, upper
     * bounds, and search index.
     *
     * @param order The order of parts in the search.
     * @param lowerBounds The lower bounds for each part.
     * @param upperBounds The upper bounds for each part.
     * @param sIdx The index of the search.
     * @return A Search object created with the given parameters.
     */
    static Search makeSearch(std::vector<length_t> order,
                             std::vector<length_t> lowerBounds,
                             std::vector<length_t> upperBounds, length_t sIdx) {
        // check correctness of sizes
        if (order.size() != lowerBounds.size() ||
            order.size() != upperBounds.size()) {
            throw std::runtime_error("Could not create search, the sizes of "
                                     "all vectors are not equal");
        }

        // compute the directions
        std::vector<Direction> directions;
        directions.reserve(order.size());
        // set the direction of the first part (copy the direction of second
        // part as first part's direction could be either way)
        directions.push_back((order[1] > order[0]) ? FORWARD : BACKWARD);

        for (length_t i = 1; i < order.size(); i++) {
            Direction d = (order[i] > order[i - 1]) ? FORWARD : BACKWARD;
            directions.push_back(d);
        }

        // compute the direction switches
        std::vector<bool> directionSwitch;
        directionSwitch.reserve(order.size());
        // first partition is not a switch
        directionSwitch.push_back(false);

        // second partition is never a switch
        directionSwitch.push_back(false);

        // add the other partitions
        for (length_t i = 2; i < directions.size(); i++) {
            directionSwitch.push_back(directions[i] != directions[i - 1]);
        }

        // compute the lowest and highest parts processed before
        std::vector<std::pair<length_t, length_t>>
            lowestAndHighestPartsProcessedBefore;
        lowestAndHighestPartsProcessedBefore.reserve(order.size());
        auto firstPair = std::make_pair(order[0], order[0]);
        lowestAndHighestPartsProcessedBefore.emplace_back(firstPair);

        for (length_t i = 1; i < order.size(); i++) {
            auto& pairBefore = lowestAndHighestPartsProcessedBefore.back();
            auto current = order[i];
            if (current < pairBefore.first) {
                lowestAndHighestPartsProcessedBefore.emplace_back(
                    current, pairBefore.second);
            } else {
                lowestAndHighestPartsProcessedBefore.emplace_back(
                    pairBefore.first, current);
            }
        }

        // check if unidirectional backwards search
        bool uniBackwards = (order[0] == order.size() - 1);
        length_t uniBackwardsIndex = order.size(); // initialize on never

        if (order.back() != 0) {
            uniBackwardsIndex = order.size(); // never uni-directional backwards
        } else if (uniBackwards) {
            uniBackwardsIndex = 0; // always unidirectional backwards
        } else {
            // find the index from which the search is only bidirectional
            // backwards this is the part just after the last part has been
            // processed
            for (length_t idx = 0; idx < order.size(); idx++) {
                if (order[idx] == order.size() - 1) {
                    uniBackwardsIndex = idx + 1;
                    break;
                }
            }
        }

        return Search(order, lowerBounds, upperBounds, directions,
                      directionSwitch, lowestAndHighestPartsProcessedBefore,
                      sIdx, uniBackwards, uniBackwardsIndex);
    }

    /**
     * @brief Static method to create a Search object with a lower bound.
     * @param originalSearch the original search
     * @param minD the minimal allowed distance (the last lower bound of the
     * search will be set to at least this value)
     */
    static Search makeSearchWithLowerBound(const Search& originalSearch,
                                           length_t minD) {
        Search s(originalSearch);
        s.setMinED(minD);
        return s;
    }

    /**
     * @brief Set the directions of the parts to the directions of the search.
     *
     * This method sets the directions for the parts based on the search's
     * directions.
     *
     * @param parts The parts to set the direction for.
     */
    void setDirectionsInParts(std::vector<Substring>& parts) const {
        // set the directions for the parts
        for (length_t i = 0; i < order.size(); i++) {
            parts[order[i]].setDirection(directions[i]);
        }
    }
    /**
     * @brief Get the lower bound for the specified part.
     *
     * @param idx The index of the part.
     * @return The lower bound for the part.
     */
    length_t getLowerBound(length_t idx) const {
        assert(idx < L.size());
        return L[idx];
    }

    /**
     * @brief Get the upper bound for the specified part.
     *
     * @param idx The index of the part.
     * @return The upper bound for the part.
     */
    length_t getUpperBound(length_t idx) const {
        assert(idx < U.size());
        return U[idx];
    }
    /**
     * @brief Get the the index of the part that is processed in phase idx
     *
     * @param idx The index of the phase.
     * @return The part that is processed in the specified phase.
     */
    length_t getPart(length_t idx) const {
        assert(idx < order.size());
        return order[idx];
    }

    /**
     * @brief Get the lowest part processed before the specified phase.
     *
     * @param idx The index of the part.
     * @return The lowest part processed before the specified phase.
     */
    length_t getLowestPartProcessedBefore(length_t idx) const {
        assert(idx >= 1);
        assert(idx < order.size());
        return lowestAndHighestPartsProcessedBefore[idx - 1].first;
    }

    /**
     * @brief Get the highest part processed before the specified phase.
     *
     * @param idx The index of the part.
     * @return The highest part processed before the specified phase.
     */
    length_t getHighestPartProcessedBefore(length_t idx) const {
        assert(idx >= 1);
        assert(idx < order.size());
        return lowestAndHighestPartsProcessedBefore[idx - 1].second;
    }

    /**
     * @brief Get the maximum allowed edit distance at the end of the search.
     *
     * @return The maximum allowed edit distance (upper bound).
     */
    length_t getMaxED() const {
        return U.back();
    }

    /**
     * @brief Get the minimum allowed edit distance at the end of the search.
     *
     * @return The minimum allowed edit distance (lower bound).
     */
    length_t getMinED() const {
        return L.back();
    }

    /**
     * @brief Get the direction for the specified phase.
     *
     * @param idx The index of the phase.
     * @return The direction for the phase (FORWARD or BACKWARD).
     */
    Direction getDirection(length_t idx) const {
        assert(idx < directions.size());
        return directions[idx];
    }
    /**
     * @brief Check if the direction switches at the specified phase.
     *
     * @param idx The index of the phase.
     * @return True if the direction switches at the phase, false otherwise.
     */
    bool getDirectionSwitch(length_t idx) const {
        assert(idx < directionSwitch.size());
        return directionSwitch[idx];
    }

    /**
     * @brief Get the number of parts in this search.
     *
     * @return The number of parts in the search.
     */
    length_t getNumParts() const {
        return order.size();
    }

    /**
     * @brief Get the index of the search.
     *
     * @return The index of the search.
     */
    length_t getIndex() const {
        return sIdx;
    }

    /**
     * @brief Check if the specified phase is the first or last part of the
     * pattern.
     *
     * @param idx The index of the phase.
     * @return True if this phase process the  the first or last part of the
     * pattern, false otherwise.
     */
    bool isEdge(length_t idx) const {
        return order[idx] == 0 || order[idx] == getNumParts() - 1;
    }

    /**
     * @brief Check if the specified phase is the final phase of the search.
     *
     * @param idx The index of the phase.
     * @return True if the phase is the final phase of the search, false
     * otherwise.
     */
    bool isEnd(length_t idx) const {
        return idx == order.size() - 1;
    }

    /**
     * @brief Check if the connectivity property is satisfied.
     *
     * The connectivity property is satisfied if the permutation is connected.
     *
     * @return True if the connectivity property is satisfied, false otherwise.
     */
    bool connectivitySatisfied() const {
        length_t highestSeen = order[0];
        length_t lowestSeen = order[0];
        for (length_t i = 1; i < order.size(); i++) {
            if (order[i] == highestSeen + 1) {
                highestSeen++;
            } else if (order[i] == lowestSeen - 1) {
                lowestSeen--;
            } else {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Check if the upper and lower bounds are valid.
     *
     * The bounds are valid if they are non-decreasing and L[i] <= U[i].
     *
     * @return True if the bounds are valid, false otherwise.
     */
    bool validBounds() const {
        // check if bounds at index 0 make sense
        if (L[0] > U[0]) {
            return false;
        }
        for (length_t i = 1; i < order.size(); i++) {
            if (L[i] > U[i] || L[i] < L[i - 1] || U[i] < U[i - 1]) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Check if the search is zero based.
     *
     * The search is zero-based if its order contains zero.
     *
     * @return True if the search is zero-based, false otherwise.
     */

    bool zeroBased() const {
        return *std::min_element(order.begin(), order.end()) == 0;
    }

    /**
     * @brief Compare two Search objects to determine their order.
     *
     * Searches are compared based on their U-strings, L-strings, and index.
     * A higher U-string corresponds to a lower search.
     *
     * @param rhs The Search object to compare with.
     * @return True if this Search is less than rhs, false otherwise.
     */
    bool operator<(const Search& rhs) const {
        assert(rhs.getNumParts() == getNumParts());

        // A) compare U-string
        for (length_t i = 0; i < getNumParts(); i++) {
            if (U[i] != rhs.U[i]) {
                return U[i] > rhs.U[i];
            }
        }

        // B) compare L-string
        for (length_t i = 0; i < getNumParts(); i++) {
            if (L[i] != rhs.L[i]) {
                return L[i] < rhs.L[i];
            }
        }

        // C) sort on index
        return sIdx < rhs.sIdx;
    }

    /**
     * @brief Check if the search covers a given distribution of errors.
     *
     * This method checks if the search covers a specified distribution of
     * errors based on the lower and upper bounds.
     *
     * @param dist The distribution to check.
     * @return True if the search covers the distribution, false otherwise.
     */
    bool coversDistribution(const Distribution& dist) const {
        unsigned int total = 0;
        for (unsigned int i = 0; i < getNumParts(); i++) {
            unsigned int part = getPart(i);
            total += dist[part];
            if (total > getUpperBound(i) || total < getLowerBound(i)) {
                return false;
            }
        }

        return true;
    }

    /**
     * @brief Sets the minimal allowed distance at the end of the search.
     * @param minED the minimal allowed distance
     */
    void setMinED(length_t minED) {
        assert(minED <= getMaxED());
        L.back() = std::max(L.back(), minED);
    }

    /**
     * Is the search unidirectional backwards from a certain index until the
     * end?
     * @param idx the index from which to check
     */
    bool isUnidirectionalBackwards(length_t idx) const {
        return unidirectionalBackwards || idx >= uniDirectionalBackwardsIndex;
    }

    /**
     * Create a search that is a copy of this search but with mirrored
     * pi-strings
     * @return the search with mirrored pi-strings
     */
    Search mirrorPiStrings() const {
        std::vector<length_t> mirroredOrder = order;
        for (length_t i = 0; i < order.size(); i++) {
            mirroredOrder[i] = order.size() - 1 - order[i];
        }
        return Search::makeSearch(mirroredOrder, L, U, sIdx);
    }
};
/**
 * Operator overloading. Outputs the search to the output stream
 * @param os, the output stream
 * @param obj, the search to print
 */
std::ostream& operator<<(std::ostream& os, const Search& obj);

/**
 * @class SearchScheme
 *
 * @brief A class representing a search scheme for approximate string matching.
 *
 */
class SearchScheme {
  private:
    const std::vector<Search>
        searches; // A vector of search operations representing the scheme.
    uint16_t k;   // The maximum number of errors allowed in the scheme.
    uint16_t
        criticalPartIndex; // The index of the critical search in the  scheme.

    /**
     * @brief Private method to determine the critical search index in the
     * scheme.
     *
     * The critical search is the one with the smallest U-string according to
     * the operator<. This corresponds to the lexicographically highest
     * U-string.
     */
    void setCriticalSearchIndex() {
        // the index of the search with the heaviest U-string
        // heaviest U-string is smallest search according to operator<
        // Find the iterator pointing to the minimum element
        auto minElementIterator =
            std::min_element(searches.begin(), searches.end());

        // Calculate the index by subtracting the begin iterator from the
        // minimum element iterator
        auto indexInSearches =
            std::distance(searches.begin(), minElementIterator);

        const auto& criticalSearch = searches[indexInSearches];
        criticalPartIndex = criticalSearch.getPart(0);
    }

    /**
     * @brief Private method to perform a sanity check on the validity of the
     * searches in the scheme.
     *
     * This method checks if the searches in the scheme satisfy various
     * properties, including the number of parts, zero-based, connectivity, and
     * valid bounds.
     *
     * @param verbose A flag indicating whether to display verbose information
     * during the check.
     */
    void sanityCheck(bool verbose) const {
        const Search& firstSearch = searches.front();
        length_t P = firstSearch.getNumParts();
        // check if all searches have same number of parts
        if (std::any_of(searches.begin(), searches.end(), [P](const Search& s) {
                return s.getNumParts() != P;
            })) {
            throw std::runtime_error("Not all searches for distance " +
                                     std::to_string(k) +
                                     " have the same number of parts");
        }

        // check if zero based
        if (std::any_of(searches.begin(), searches.end(),
                        [](const Search& s) { return !s.zeroBased(); })) {
            throw std::runtime_error(
                "Not all searches are zero based for distance " +
                std::to_string(k) + "!");
        }

        // check if connectivity satisfied
        if (std::any_of(searches.begin(), searches.end(), [](const Search& s) {
                return !s.connectivitySatisfied();
            })) {
            throw std::runtime_error("Connectivity property not satisfied "
                                     "for all searches with distance " +
                                     std::to_string(k) + "!");
        }

        // check if U and L std::string are valid
        if (std::any_of(searches.begin(), searches.end(),
                        [](const Search& s) { return !s.validBounds(); })) {
            throw std::runtime_error("Decreasing lower or upper bounds "
                                     "for a search for K  = " +
                                     std::to_string(k));
        }
    }

    /**
     * @brief Static method to create a Search object from a line of text.
     *
     * This method parses a line of text to create a Search object.
     *
     * @param line The line of text containing search information.
     * @param idx The index of the search within the scheme.
     * @return A Search object created from the provided information.
     */
    static Search makeSearch(const std::string& line, length_t idx) {
        std::stringstream ss(line);

        std::vector<std::string> tokens;
        std::string token;
        while (ss >> token) {
            tokens.push_back(token);
        }

        if (tokens.size() != 3) {
            throw std::runtime_error("A search should have 3 vectors: order, "
                                     "lower bound and upper bound!");
        }

        std::vector<length_t> order;
        getVector(tokens[0], order);

        std::vector<length_t> lower_bound;
        getVector(tokens[1], lower_bound);

        std::vector<length_t> upper_bound;
        getVector(tokens[2], upper_bound);

        return Search::makeSearch(order, lower_bound, upper_bound, idx);
    }

    /**
     * @brief Static method to extract a vector of integers from a string
     * representation between curly brackets and separated by comma's.
     *
     * This method parses a string containing a vector of integers and converts
     * it into a vector of integers.
     *
     * @param vectorString The string representation of the vector.
     * @param vector The vector to store the parsed integers.
     */
    static void getVector(const std::string& vectorString,
                          std::vector<length_t>& vector) {

        if (vectorString.size() < 2) {
            throw std::runtime_error(vectorString +
                                     " is not a valid vector for a search");
        }
        std::string bracketsRemoved =
            vectorString.substr(1, vectorString.size() - 2);

        std::stringstream ss(bracketsRemoved);
        std::string token;
        while (std::getline(ss, token, ',')) {
            vector.emplace_back(stoull(token));
        }
    }

    /**
     * @brief Static method to generate all error distributions with a given
     * number of parts and errors.
     *
     * This method generates error distributions with a specified number of
     * parts and allowed errors.
     *
     * @param P The number of parts in the error distributions.
     * @param K The number of allowed errors.
     * @param distributions A vector to store the generated error distributions.
     */
    static void
    genErrorDistributions(int P, int K,
                          std::vector<Distribution>& distributions) {
        Distribution distribution(P, 0);

        for (int i = 0; i < pow(K + 1, P); i++) {
            int sum =
                std::accumulate(distribution.begin(), distribution.end(), 0);
            if (sum <= K)
                distributions.push_back(distribution);

            for (int j = 0; j < P; j++) {
                distribution[j]++;
                if (distribution[j] != K + 1)
                    break;
                distribution[j] = 0;
            }
        }
    }

  public:
    /**
     * @brief Constructor for the SearchScheme class.
     *
     * Initializes a SearchScheme object with a vector of searches and a maximum
     * error constraint. Also finds the index of the critical search.
     *
     * @param searches A vector of Search objects representing the search
     * operations.
     * @param k The maximum number of errors allowed in the scheme.
     * @param verbose A flag indicating whether to display verbose information
     * during the initialization.
     */
    SearchScheme(const std::vector<Search>& searches, unsigned int k,
                 bool verbose = false)
        : searches(searches), k(k) {
        sanityCheck(verbose);
        setCriticalSearchIndex();
    }

    /**
     * @brief Static method to read and create a SearchScheme from a file.
     *
     * This method reads search information from a file and creates a
     * SearchScheme object.
     *
     * @param stream_searches An input file stream containing search
     * information.
     * @param fileName The name of the file being processed.
     * @param k The maximum number of errors allowed in the scheme.
     * @return A SearchScheme object created from the file data.
     */
    static SearchScheme readScheme(std::ifstream& stream_searches,
                                   const std::string& fileName,
                                   unsigned int k) {

        length_t sIdx = 0; // the index of the current search

        // read the searches line by line
        std::vector<Search> rSearches;
        std::string line;
        while (std::getline(stream_searches, line)) {
            try {
                rSearches.push_back(makeSearch(line, sIdx++));
            } catch (const std::runtime_error& e) {
                throw std::runtime_error(
                    "Something went wrong with processing line: " + line +
                    "\nin file: " + fileName + "\n" + e.what());
            }
        }

        if (rSearches.empty()) {
            throw std::runtime_error("Empty scheme in: " + fileName);
        }

        return SearchScheme(rSearches, k);
    }

    /**
     * @brief Get the vector of search operations in the scheme.
     *
     * @return A reference to the vector of search operations in the scheme.
     */
    const std::vector<Search>& getSearches() const {
        return searches;
    }

    /**
     * @brief Get the index of the critical search in the scheme.
     *
     *
     * @return The index of the critical search.
     */
    uint16_t getCriticalPartIndex() const {
        return criticalPartIndex;
    }

    /**
     * @brief Get the number of parts in the search scheme.
     *
     * @return The number of parts in the search scheme.
     */
    uint16_t getNumParts() const {
        return searches.front().getNumParts();
    }

    /**
     * @brief Mirror the pi-strings of all searches in the scheme
     * @return the scheme with mirrored pi-strings
     */
    SearchScheme mirrorPiStrings() const {
        std::vector<Search> reversedSearches;
        for (const Search& s : searches) {
            reversedSearches.push_back(s.mirrorPiStrings());
        }
        return SearchScheme(reversedSearches, k);
    }

    unsigned int getK() const {
        return k;
    }
};
#endif
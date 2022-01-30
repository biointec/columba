/******************************************************************************
 *  Columba 1.1: Approximate Pattern Matching using Search Schemes            *
 *  Copyright (C) 2020-2022 - Luca Renders <luca.renders@ugent.be> and        *
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
#include "searchstrategy.h"

using namespace std;

// ============================================================================
// CLASS SEARCHSTRATEGY
// ============================================================================

// ----------------------------------------------------------------------------
// CONSTRUCTOR
// ----------------------------------------------------------------------------

SearchStrategy::SearchStrategy(FMIndex& argument, PartitionStrategy p,
                               DistanceMetric distanceMetric)
    : index(argument), partitionStrategy(p), distanceMetric(distanceMetric) {

    // set the partition strategy
    switch (p) {
    case UNIFORM:
        partitionPtr = &SearchStrategy::partitionUniform;
        break;
    case DYNAMIC:
        partitionPtr = &SearchStrategy::partitionDynamic;
        break;
    case STATIC:
        partitionPtr = &SearchStrategy::partitionOptimalStatic;
        break;
    default:
        break;
    }

    // set the distance metric
    switch (distanceMetric) {
    case HAMMING:
        startIdxPtr = &SearchStrategy::startIndexHamming;
        break;
    case EDITNAIVE:
        startIdxPtr = &SearchStrategy::startIndexEditNaive;
        break;
    case EDITOPTIMIZED:
        startIdxPtr = &SearchStrategy::startIndexEditOptimized;
    default:
        break;
    }
}

// ----------------------------------------------------------------------------
// INFORMATION
// ----------------------------------------------------------------------------
string SearchStrategy::getPartitioningStrategy() const {
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

string SearchStrategy::getDistanceMetric() const {
    switch (distanceMetric) {
    case HAMMING:
        return "HAMMING";
        break;
    case EDITOPTIMIZED:
        return "(OPTIMIZED) EDIT";
        break;
    case EDITNAIVE:
        return "(NAIVE) EDIT";
        break;

    default:
        // should not get here
        return "";
    }
}
// ----------------------------------------------------------------------------
// SANITY CHECKS
// ----------------------------------------------------------------------------

void SearchStrategy::genErrorPatterns(int P, int K, vector<Pattern>& patterns) {
    Pattern pattern(P, 0);

    for (int i = 0; i < pow(K + 1, P); i++) {
        int sum = accumulate(pattern.begin(), pattern.end(), 0);
        if (sum <= K)
            patterns.push_back(pattern);

        for (int j = 0; j < P; j++) {
            pattern[j]++;
            if (pattern[j] != K + 1)
                break;
            pattern[j] = 0;
        }
    }
}

bool SearchStrategy::coversPatterns(const vector<Pattern>& patterns,
                                    const vector<Search>& scheme,
                                    bool verbose) {
    vector<int> numCover(scheme.size(), 0);

    bool ret = true;

    // and error pattern is simply a vector containing
    // the number of errors in each partition P
    for (const Pattern& pattern : patterns) {
        int numberOfCovers = 0;

        vector<int> searchesThatCover;
        // check if a search covers the pattern
        bool patternCovered = false;
        for (size_t si = 0; si < scheme.size(); si++) {
            const Search& s = scheme[si];
            bool thisCover = true;
            length_t numErrors = 0;

            // check if search covers pattern
            for (length_t i = 0; i < s.getNumParts(); i++) {
                length_t p = s.getPart(i);
                numErrors += pattern[p];
                if ((numErrors > s.getUpperBound(i)) ||
                    (numErrors < s.getLowerBound(i)))
                    thisCover = false;
            }

            // print the pattern and the search that covers it
            if (thisCover) {
                searchesThatCover.push_back(si + 1);
                numberOfCovers++;
                numCover[si]++;

                if (verbose) {
                    cout << "Pattern: ";
                    for (length_t i = 0; i < pattern.size(); i++)
                        cout << pattern[i];
                    cout << " is covered by search:" << si + 1 << " (";
                    for (length_t i = 0; i < pattern.size(); i++)
                        cout << s.getPart(i);
                    cout << "), (";
                    for (length_t i = 0; i < s.getNumParts(); i++)
                        cout << s.getLowerBound(i);
                    cout << "), (";
                    for (length_t i = 0; i < s.getNumParts(); i++)
                        cout << s.getUpperBound(i);
                    cout << ")";

                    if (patternCovered)
                        cout << " * "; // * means redundant
                    cout << endl;
                }
                patternCovered = true;
            }
        }

        if (!patternCovered && verbose) {
            cout << "Pattern is not covered: ";
            for (length_t i = 0; i < pattern.size(); i++)
                cout << pattern[i];
            cout << endl;
        }

        ret &= patternCovered;

        if (!ret && !verbose)
            return ret;

        if (verbose) {
            cout << "Pattern ";
            for (length_t i = 0; i < pattern.size(); i++)
                cout << pattern[i];
            cout << " is covered by searches: ";
            for (auto si : searchesThatCover)
                cout << si << ",";
            cout << "\n";
        }
    }

    if (verbose)
        cout << patterns.size() << " patterns covered" << endl;

    // check if all searches cover at least one pattern
    for (size_t i = 0; i < numCover.size(); i++) {
        if (verbose)
            cout << "Search " << i << " is used " << numCover[i] << " times\n";

        if (numCover[i] == 0)
            cout << "Warning: search " << scheme[i]
                 << " covers no error patterns!\n";
    }

    return ret;
}

// ----------------------------------------------------------------------------
// PARTITIONING
// ----------------------------------------------------------------------------
void SearchStrategy::partition(const string& pattern, vector<Substring>& parts,
                               const int& numParts, const int& maxScore,
                               vector<SARangePair>& exactMatchRanges,
                               Counters& counters) const {

    parts.clear();

    if (numParts >= (int)pattern.size() || numParts == 1) {
        // no need of splitting up since all parts would be one
        // character or less or there is only one part
        return;
    }
    parts.reserve(numParts);
    (this->*partitionPtr)(pattern, parts, numParts, maxScore, exactMatchRanges,
                          counters);
}

// Uniform Partitioning

void SearchStrategy::partitionUniform(const string& pattern,
                                      vector<Substring>& parts,
                                      const int& numParts, const int& maxScore,
                                      vector<SARangePair>& exactMatchRanges,
                                      Counters& counters) const {

    for (int i = 0; i < numParts; i++) {
        parts.emplace_back(pattern, (i * 1.0 / numParts) * pattern.size(),
                           ((i + 1) * 1.0 / numParts) * pattern.size());
    }
    // set end of final part correct
    parts.back().setEnd(pattern.size());

    // match the exactRanges for each part
    index.setDirection(FORWARD);
    SARangePair initialRanges = index.getCompleteRange();

    vector<bool> partNumberSeen(numParts, false);

    for (int i = 0; i < numParts; i++) {
        exactMatchRanges[i] =
            index.matchStringBidirectionally(parts[i], initialRanges, counters);
    }
}

// Static Partitioning
void SearchStrategy::partitionOptimalStatic(
    const string& pattern, vector<Substring>& parts, const int& numParts,
    const int& maxScore, vector<SARangePair>& exactMatchRanges,
    Counters& counters) const {

    setParts(pattern, parts, numParts, maxScore);

    // match the exactRanges for each part
    index.setDirection(FORWARD);
    SARangePair initialRanges = index.getCompleteRange();

    for (int i = 0; i < numParts; i++) {
        exactMatchRanges[i] =
            index.matchStringBidirectionally(parts[i], initialRanges, counters);
    }
}

void SearchStrategy::setParts(const string& pattern, vector<Substring>& parts,
                              const int& numParts, const int& maxScore) const {
    const vector<double>& begins = getBegins(numParts, maxScore);

    int pSize = pattern.size();
    // set the first part
    parts.emplace_back(pattern, 0, begins[0] * pSize);

    for (unsigned int i = 0; i < begins.size() - 1; i++) {
        parts.emplace_back(pattern, begins[i] * pSize, begins[i + 1] * pSize);
    }
    parts.emplace_back(pattern, begins.back() * pSize, pattern.size());
}

// Dynamic Partitioning

void SearchStrategy::partitionDynamic(const string& pattern,
                                      vector<Substring>& parts,
                                      const int& numParts, const int& maxScore,
                                      vector<SARangePair>& exactMatchRanges,
                                      Counters& counters) const {

    int matchedChars =
        seed(pattern, parts, numParts, maxScore, exactMatchRanges);
    int pSize = pattern.size();
    vector<int> weights = getWeights(numParts, maxScore);

    Direction dir = FORWARD;
    int partToExtend = 0;

    // extend the part with the largest range, as to minimize the range
    // for each part do this until all characters are assigned to a
    // part
    for (int j = matchedChars; j < pSize; j++) {

        // find the part with the largest range
        length_t maxRangeWeighted = 0;

        for (int i = 0; i < numParts; i++) {
            bool noLeftExtension =
                (i == 0) || parts[i].begin() == parts[i - 1].end();
            bool noRightExtension =
                (i == numParts - 1) || parts[i].end() == parts[i + 1].begin();
            if (noLeftExtension && noRightExtension) {
                continue;
            }
            if (exactMatchRanges[i].width() * weights[i] > maxRangeWeighted) {
                maxRangeWeighted = exactMatchRanges[i].width() * weights[i];
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

        if (maxRangeWeighted == 0) {
            // no need to keep calculating new range, just extend the
            // parts
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
        index.setDirection(dir);
        index.addChar(c, exactMatchRanges.at(partToExtend), counters);
    }
}

int SearchStrategy::seed(const string& pattern, vector<Substring>& parts,
                         const int& numParts, const int& maxScore,
                         vector<SARangePair>& exactMatchRanges) const {
    int pSize = pattern.size();
    bool useKmerTable = (pSize >= 100);
    int wSize = (useKmerTable) ? index.getWordSize() : 1;

    const auto& seedPercent = getSeedingPositions(numParts, maxScore);

    vector<int> seeds;
    // push the seed for the first part
    seeds.emplace_back(0);

    // push the optimal seeds for the middle parts
    for (int i = 1; i < numParts - 1; i++) {
        seeds.emplace_back((seedPercent[i - 1] * pSize) - (wSize / 2));
    }

    for (int i = 0; i < numParts - 1; i++) {
        parts.emplace_back(pattern, seeds[i], seeds[i] + wSize);
    }

    // push the seeds for the final parts
    parts.emplace_back(pattern, pSize - wSize, pSize);

    exactMatchRanges.resize(numParts);
    for (int i = 0; i < numParts; i++) {
        exactMatchRanges[i] = (useKmerTable)
                                  ? index.lookUpInKmerTable(parts[i])
                                  : index.getRangeOfSingleChar(parts[i][0]);
    }
    return numParts * wSize;
}

void SearchStrategy::extendParts(const string& pattern,
                                 vector<Substring>& parts) const {
    for (length_t i = 0; i < parts.size(); i++) {

        if ((i != parts.size() - 1) &&
            (parts[i].end() != parts[i + 1].begin())) {
            // extend completely to the right
            // it is known that the range will stay [0,0)
            parts[i].setEnd(parts[i + 1].begin());
        }
        if ((i != 0) && (parts[i].begin() != parts[i - 1].end())) {
            // extend completely to the left
            parts[i].setBegin(parts[i - 1].end());
        }
    }
}

// ----------------------------------------------------------------------------
// (APPROXIMATE) MATCHING
// ----------------------------------------------------------------------------
vector<TextOcc> SearchStrategy::matchApprox(const string& pattern,
                                            length_t maxED,
                                            Counters& counters) const {

    if (maxED == 0) {
        return index.exactMatchesOutput(pattern, counters);
    }
    // create the parts of the pattern
    vector<Substring> parts;

    // calculate how many parts there will be
    uint numParts = calculateNumParts(maxED);

    // the ranges corresponding to the exact match for each part
    vector<SARangePair> exactMatchRanges(numParts);

    // partition the read
    partition(pattern, parts, numParts, maxED, exactMatchRanges, counters);

    if (parts.empty() || numParts * maxED >= pattern.size()) {
        // splitting up was not viable -> just search the entire pattern
        cerr << "Warning: Normal bidirectional search was used as "
                "entered pattern is too short "
             << pattern.size() << endl;

        return index.approxMatchesNaive(pattern, maxED, counters);
    }

    // The occurrences in the text and index
    Occurrences occ;

    // create the bit-parallel matrix for in-text verification
    BitParallelED intextMatrix;
    // set sequence to this matrix
    intextMatrix.setSequence(pattern);

    // END of preprocessing

    // A) do in-text verification for the parts that occur less than the index's
    // switch point
    for (uint i = 0; i < numParts; i++) {
        size_t width = exactMatchRanges[i].width();
        if (width != 0 && width <= index.getSwitchPoint()) {
            // do in-text verification on this part
            const auto& part = parts[i];
            FMOcc startMatch(exactMatchRanges[i], 0, part.size());
            index.verifyExactPartialMatchInText(
                intextMatrix, startMatch, part.begin(), maxED, occ, counters);
        }
    }

    // B) do the searches for which the first part occurs more than the switch
    // point

    // create the searches
    const vector<Search>& searches = createSearches(maxED);
    index.reserveStacks(numParts,
                        pattern.length()); // reserve stacks for each part

    // create the bit-parallel alignment matrices
    index.resetMatrices(parts.size()); // reset the alignment matrix that will
                                       // be (possibly) used for each part

    for (const Search& s : searches) {
        doRecSearch(intextMatrix, s, parts, occ, exactMatchRanges, counters);
    }

    // return all matches mapped to the text
    return (distanceMetric != HAMMING)
               ? occ.getUniqueTextOccurrences(index, maxED, intextMatrix,
                                              counters)
               : occ.getTextOccurrencesHamming(index, pattern.size(), counters);
}

void SearchStrategy::doRecSearch(BitParallelED& intextMatrix, const Search& s,
                                 vector<Substring>& parts, Occurrences& occ,
                                 const vector<SARangePair>& exactMatchRanges,
                                 Counters& counters) const {

    if (s.getUpperBound(0) > 0) {
        // first part is allowed an error so start with an empty match
        s.setDirectionsInParts(parts);

        SARangePair startRange = index.getCompleteRange();
        FMOcc startMatch = FMOcc(startRange, 0, 0);
        (this->*startIdxPtr)(intextMatrix, s, startMatch, occ, parts, counters,
                             0);
        return;
    }

    // first get the bidirectional match of first part
    int first = s.getPart(0);
    SARangePair startRange = exactMatchRanges[first];

    // if this range is bigger than the switch point
    if (startRange.width() > index.getSwitchPoint()) {

        // prepare the parts for this search
        s.setDirectionsInParts(parts);

        // can we continue exact matching according to this search?
        uint partInSearch = 1;
        length_t exactLength = parts[first].size();

        while (s.getUpperBound(partInSearch) == 0) {
            // extend the exact match
            index.setDirection(s.getDirection(partInSearch - 1));
            const auto& part = parts[s.getPart(partInSearch)];

            startRange =
                index.matchStringBidirectionally(part, startRange, counters);
            if (startRange.empty()) {
                return;
            }

            exactLength += parts[s.getPart(partInSearch)].size();
            partInSearch++;
        }

        // Create a match corresponding to the exact match
        FMOcc startMatch = FMOcc(startRange, 0, exactLength);
        // Start approximate matching in the index
        (this->*startIdxPtr)(intextMatrix, s, startMatch, occ, parts, counters,
                             partInSearch);
    }
}
// ============================================================================
// CLASS CUSTOMSEARCHSTRATEGY
// ============================================================================

// ----------------------------------------------------------------------------
// CONSTRUCTION
// ----------------------------------------------------------------------------

void CustomSearchStrategy::getSearchSchemeFromFolder(string pathToFolder,
                                                     bool verbose) {

    // get the name of the file
    string line;
    {
        ifstream ifs(pathToFolder + "name.txt");
        if (!ifs) {
            throw runtime_error("Problem reading: " + pathToFolder +
                                "name.txt\nDid you provide a directory to "
                                "a search scheme without a name file?");
        }
        getline(ifs, line);
        name = line;
        ifs.close();
    }

    // get the info per distance score (scores between 1 and 4 are looked
    // at)
    for (int i = 1; i <= 4; i++) {

        ifstream stream_searches(pathToFolder + to_string(i) + "/searches.txt");
        if (!stream_searches) {
            // this score is not supported
            supportsMaxScore[i - 1] = false;
            continue;
        }
        // max score is supported
        // read the searches line by line
        while (getline(stream_searches, line)) {
            try {
                schemePerED[i - 1].push_back(makeSearch(line));
            } catch (const runtime_error& e) {
                throw runtime_error(
                    "Something went wrong with processing line: " + line +
                    "\nin file: " + pathToFolder + to_string(i) +
                    "/searches.txt\n" + e.what());
            }
        }
        if (schemePerED[i - 1].size() > 0) {
            supportsMaxScore[i - 1] = true;
        }
        stream_searches.close();
    }

    // check if the searches are valid
    sanityCheck(verbose);

    // get static positions (if they exist)
    for (int i = 1; i <= 4; i++) {
        ifstream stream_static(pathToFolder + to_string(i) +
                               "/static_partitioning.txt");

        if (stream_static) {
            // a file with static partitioning positions exists
            getline(stream_static, line);
            vector<string> postionsAsString = {};
            stringstream ss(line);
            string token;
            while (ss >> token) {
                postionsAsString.push_back(token);
            }

            if (postionsAsString.size() != calculateNumParts(i) - 1) {
                throw runtime_error(
                    "Not enough static positions provided in " + pathToFolder +
                    to_string(i) + "/static_partitioning.txt\nExpected: " +
                    to_string(calculateNumParts(i) - 1) + " parts\nProvided: " +
                    to_string(postionsAsString.size()) + " parts");
            }

            for (auto str : postionsAsString) {
                staticPositions[i - 1].push_back(stod(str));
            }

            // check if these positions are valid
            sanityCheckStaticPartitioning(i);
            // if valid set getBeginsPointer to custom
            beginsPointer[i - 1] = &CustomSearchStrategy::getBeginsCustom;

            stream_static.close();
        }
    }

    // get dynamic seeds and weights (if file exists)
    for (int i = 1; i <= 4; i++) {
        ifstream stream_dynamic(pathToFolder + to_string(i) +
                                "/dynamic_partitioning.txt");

        if (stream_dynamic) {
            // a file with dynamic partitioning positions exists
            getline(stream_dynamic, line);
            vector<string> seedsAsString = {};
            stringstream ss(line);
            string token;
            while (ss >> token) {
                seedsAsString.push_back(token);
            }

            if (seedsAsString.size() != calculateNumParts(i) - 2) {
                throw runtime_error(
                    "Not enough seeding positions provided in " + pathToFolder +
                    to_string(i) + "/dynamic_partitioning.txt\nExpected: " +
                    to_string(calculateNumParts(i) - 1) + " seeds\nProvided: " +
                    to_string(seedsAsString.size()) + " seeds");
            }

            for (auto str : seedsAsString) {
                seedingPositions[i - 1].push_back(stod(str));
            }

            // check if these seeds are valid
            sanityCheckDynamicPartitioning(i);

            // get the weights
            getline(stream_dynamic, line);
            stringstream ss_w(line);
            string stringWeight;
            while (ss_w >> stringWeight) {
                weights[i - 1].push_back(stoi(stringWeight));
            }

            if (weights[i - 1].size() != calculateNumParts(i)) {
                throw runtime_error(
                    "Not enough weights provided for max score " +
                    to_string(i));
            }

            // set the pointers to custom
            seedingPointer[i - 1] =
                &CustomSearchStrategy::getSeedingPositionsCustom;
            weightsPointers[i - 1] = &CustomSearchStrategy::getWeightsCustom;
        }
    }
}

Search CustomSearchStrategy::makeSearch(const string& line) const {
    stringstream ss(line);

    vector<string> tokens;
    string token;
    while (ss >> token) {
        tokens.push_back(token);
    }

    if (tokens.size() != 3) {
        throw runtime_error("A search should have 3 vectors: order, "
                            "lower bound and upper bound!");
    }

    vector<length_t> order;
    getVector(tokens[0], order);

    vector<length_t> lower_bound;
    getVector(tokens[1], lower_bound);

    vector<length_t> upper_bound;
    getVector(tokens[2], upper_bound);

    return Search::makeSearch(order, lower_bound, upper_bound);
}

void CustomSearchStrategy::getVector(const string& vectorString,
                                     vector<length_t>& vector) const {

    if (vectorString.size() < 2) {
        throw runtime_error(vectorString +
                            " is not a valid vector for a search");
    }
    string bracketsRemoved = vectorString.substr(1, vectorString.size() - 2);

    stringstream ss(bracketsRemoved);
    string token;
    while (getline(ss, token, ',')) {
        vector.emplace_back(stoull(token));
    }
}

// ----------------------------------------------------------------------------
// SANITY CHECKS
// ----------------------------------------------------------------------------

void CustomSearchStrategy::sanityCheckStaticPartitioning(
    const int& maxScore) const {
    const auto& positions = staticPositions[maxScore - 1];

    // no length zero + increasing + all smaller than 1 and greater than 0
    for (unsigned int i = 0; i < positions.size(); i++) {
        if (positions[i] <= 0 || positions[i] >= 1) {
            throw runtime_error("One of the provided static positions for " +
                                to_string(maxScore) +
                                " is not between 0 and 1 (exclusive)");
        }
        if (i < positions.size() - 1 && positions[i] - positions[i + 1] >= 0) {
            throw runtime_error("Provided static positions for " +
                                to_string(maxScore) +
                                " are not strictly increasing");
        }
    }
}
void CustomSearchStrategy::sanityCheckDynamicPartitioning(
    const int& maxScore) const {

    const auto& seeds = seedingPositions[maxScore - 1];

    // no length zero + increasing + all smaller than 1 and greater than 0
    for (unsigned int i = 0; i < seeds.size(); i++) {
        if (seeds[i] <= 0 || seeds[i] >= 1) {
            throw runtime_error("One of the provided static positions for " +
                                to_string(maxScore) +
                                " is not between 0 and 1 (exclusive)!");
        }
        if (i < seeds.size() - 1 && seeds[i] - seeds[i + 1] >= 0) {
            throw runtime_error("Provided seeding positions for " +
                                to_string(maxScore) +
                                " are not strictly increasing");
        }
    }
}

void CustomSearchStrategy::sanityCheck(bool verbose) const {

    // check if for each supported edit distance all error patterns  are
    // covered
    for (int K = 1; K <= 4; K++) {

        const auto& scheme = schemePerED[K - 1];
        if (!supportsMaxScore[K - 1]) {
            continue;
        }
        const Search& firstSearch = scheme.front();
        length_t P = firstSearch.getNumParts();
        // check if all searches have same number of parts
        if (any_of(scheme.begin(), scheme.end(),
                   [P](const Search& s) { return s.getNumParts() != P; })) {
            throw runtime_error("Not all searches for distance " +
                                to_string(K) +
                                " have the same number of parts");
        }

        // check if zero based
        if (any_of(scheme.begin(), scheme.end(),
                   [](const Search& s) { return !s.zeroBased(); })) {
            throw runtime_error(
                "Not all searches are zero based for distance " + to_string(K) +
                "!");
        }

        // check if connectivity satisfied
        if (any_of(scheme.begin(), scheme.end(), [](const Search& s) {
                return !s.connectivitySatisfied();
            })) {
            throw runtime_error("Connectivity property not satisfied "
                                "for all searches with distance " +
                                to_string(K) + "!");
        }

        // check if U and L string are valid
        if (any_of(scheme.begin(), scheme.end(),
                   [](const Search& s) { return !s.noDecreasingInBounds(); })) {
            throw runtime_error("Decreasing lower or upper bounds "
                                "for a search for K  = " +
                                to_string(K));
        }
        vector<Pattern> errorPatterns;
        genErrorPatterns(P, K, errorPatterns);
        if (!coversPatterns(errorPatterns, scheme, verbose)) {
            throw runtime_error("Search scheme does not cover all "
                                "error patterns for K = " +
                                to_string(K) + "!");
        }
    }
}

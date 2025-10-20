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

#include "searchstrategy.h"
#include "logger.h"    // for Logger, logger
#include "substring.h" // for Substring

#include <sstream> // for char_traits, basic_istream, ifstream, strings...
using namespace std;

// ============================================================================
// CLASS SEARCH STRATEGY
// ============================================================================

// ----------------------------------------------------------------------------
// CONSTRUCTOR
// ----------------------------------------------------------------------------

SearchStrategy::SearchStrategy(IndexInterface& argument, PartitionStrategy p,
                               DistanceMetric distanceMetric, MappingMode m,
                               SequencingMode sequencingMode, length_t useKmerCutOff)
    : index(argument), distanceMetric(distanceMetric), partitionStrategy(p),
      mode(m), useKmerCutOff(useKmerCutOff) {

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
        filterPtr = &SearchStrategy::filterHamming;
        naiveMatchPtr = &IndexInterface::approxMatchesNaiveHamming;
#ifndef RUN_LENGTH_COMPRESSION
        inTextVerificationPtr = &SearchStrategy::inTextVerificationHamming;
#endif
        break;
    case EDIT:
        startIdxPtr = &SearchStrategy::startIndexEdit;
        naiveMatchPtr = &IndexInterface::approxMatchesNaive;
        filterPtr = (sequencingMode == SINGLE_END)
                        ? &SearchStrategy::filterEditWithCIGARCalculation
                        : &SearchStrategy::filterEditWithoutCIGARCalculation;

#ifndef RUN_LENGTH_COMPRESSION
        inTextVerificationPtr = &SearchStrategy::inTextVerificationEdit;

#endif

    default:
        break;
    }

    // set the mode
    setMappingMode(m);
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
    case EDIT:
        return "EDIT";
        break;

    default:
        // should not get here
        return "";
    }
}

string SearchStrategy::getMappingModeString() const {
    switch (mode) {
    case BEST:
        return "BEST";
        break;
    case ALL:
        return "ALL";
        break;

    default:
        // should not get here
        return "";
    }
}

// ----------------------------------------------------------------------------
// PARTITIONING
// ----------------------------------------------------------------------------
void SearchStrategy::partition(const string& pattern, vector<Substring>& parts,
                               const int& numParts, const int& maxScore,
                               vector<SARangePair>& exactMatchRanges,
                               Counters& counters) const {
    assert(exactMatchRanges.size() == (uint32_t)numParts);
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

void SearchStrategy::calculateExactMatchRanges(
    vector<Substring>& parts, vector<SARangePair>& exactMatchRanges,
    Counters& counters) const {

    index.setDirection(FORWARD, false); // Bidirectional ranges
    length_t wordSize = index.getWordSize();

    // Match exact ranges for each part bidirectionally except the last
    for (length_t i = 0; i < parts.size() - 0; ++i) {
        auto& current = parts[i];
        auto size = current.size();
        auto start = current.begin() + ((size >= wordSize) ? wordSize : 0);
        auto initRanges =
            (size >= wordSize)
                ? index.lookUpInKmerTable({current, current.begin(), start})
                : index.getCompleteRange();
        exactMatchRanges[i] = index.matchStringBidirectionally(
            {current, start, current.end()}, initRanges, counters);
    }

    // Last part matched unidirectional backwards
    index.setDirection(BACKWARD, true);
    auto& lastPart = parts.back();
    lastPart.setDirection(BACKWARD);
    auto size = lastPart.size();
    auto end = (size >= wordSize) ? lastPart.end() - wordSize : lastPart.end();
    auto initRanges =
        (size >= wordSize)
            ? index.lookUpInKmerTable({lastPart, end, lastPart.end()})
            : index.getCompleteRange();
    exactMatchRanges.back() = index.matchStringBidirectionally(
        {lastPart, lastPart.begin(), end}, initRanges, counters);
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

    // calculate the exact match ranges
    calculateExactMatchRanges(parts, exactMatchRanges, counters);
}

// Static Partitioning
void SearchStrategy::partitionOptimalStatic(
    const string& pattern, vector<Substring>& parts, const int& numParts,
    const int& maxScore, vector<SARangePair>& exactMatchRanges,
    Counters& counters) const {

    setParts(pattern, parts, numParts, maxScore);
    calculateExactMatchRanges(parts, exactMatchRanges, counters);
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

    // assert part i ends where part i+1 begins
    for (unsigned int i = 0; i < parts.size() - 1; i++) {
        assert(parts[i].end() == parts[i + 1].begin());
    }
}

// Dynamic Partitioning

/**
 * @brief Seed the parts of the pattern with no k-mer look-up table, parts will
 * be four characters long (matched in the index) if possible. Helper function
 * for SearchStrategy::partitionDynamic.
 * @param pSize The size of the pattern.
 * @param parts The parts of the pattern.
 * @param numParts The number of parts.
 * @param exactMatchRanges The exact match ranges of the parts.
 * @param counters The performance counters.
 * @param index The index.
 * @param pattern The pattern.
 * @param matchedChars The number of characters matched in the index. (output)
 */
void seedIfNoKmers(const length_t& pSize, vector<Substring>& parts,
                   const int& numParts, vector<SARangePair>& exactMatchRanges,
                   Counters& counters, IndexInterface& index,
                   const string& pattern, int& matchedChars) {
    // wordsize was 0, all parts are of size 0
    // extend all parts but the last in forwards direction to size 4 (or 1 if
    // too short)
    length_t size = (pSize >= 4 * parts.size()) ? 4 : 1;
    for (int i = 0; i < numParts - 1; i++) {
        parts[i].setEnd(min(parts[i].end() + size, pSize));
        index.setDirection(FORWARD, false);
        for (length_t j = parts[i].begin(); j < parts[i].end(); j++) {
            index.addChar(pattern[j], exactMatchRanges.at(i), counters);
            matchedChars++;
        }
    }
}

/**
 * Helper function for dynamic partitioning. This function extends the parts
 * so that nothing of the pattern is not allocated to any part. This does
 * not keep track of the ranges over the suffix array, so should only be
 * called if this does not matter (e.g. when the parts that can be extended
 * all correspond to empty ranges)
 * @param pattern the pattern that is split
 * @param parts  the current parts, they are updated so that all characters
 * of pattern are part of exactly one part
 */
void extendParts(const string& pattern, vector<Substring>& parts) {
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

void SearchStrategy::partitionDynamic(const string& pattern,
                                      vector<Substring>& parts,
                                      const int& numParts, const int& maxScore,
                                      vector<SARangePair>& exactMatchRanges,
                                      Counters& counters) const {

    int matchedChars =
        seed(pattern, parts, numParts, maxScore, exactMatchRanges);

    length_t pSize = pattern.size();
    vector<uint64_t> weights = getWeights(numParts, maxScore);

    Direction dir = FORWARD;
    int partToExtend = 0;

    if (matchedChars == 0) {
        // no characters were matched, seed the parts with no k-mer look-up
        // table
        seedIfNoKmers(pSize, parts, numParts, exactMatchRanges, counters, index,
                      pattern, matchedChars);
    }

    // extend the part with the largest range, as to minimize the range
    // for each part do this until all characters are assigned to a
    // part
    for (length_t j = matchedChars; j < pSize; j++) {

        // find the part with the largest range
        uint64_t maxRangeWeighted = 0;

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

        // match the new character, if this is the last part use unidirectional
        // backward matching
        index.setDirection(dir, partToExtend == numParts - 1);
        index.addChar(c, exactMatchRanges.at(partToExtend), counters);
    }
}

int SearchStrategy::seed(const string& pattern, vector<Substring>& parts,
                         const int& numParts, const int& maxScore,
                         vector<SARangePair>& exactMatchRanges) const {
    auto pSize = pattern.size();
    bool useKmerTable = (numParts * index.getWordSize() < (pSize * 2) / 3) && useKmerTableInSeeding(pSize);
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

    // assert part i ends before part i+1 starts
    for (int i = 0; i < numParts - 1; i++) {
        assert(parts[i].end() <= parts[i + 1].begin());
    }

    exactMatchRanges.resize(numParts);
    for (int i = 0; i < numParts; i++) {
        exactMatchRanges[i] = (useKmerTable)
                                  ? index.lookUpInKmerTable(parts[i])
                                  : index.getRangeOfSingleChar(parts[i][0]);
    }

    return numParts * wSize;
}

// ----------------------------------------------------------------------------
// (APPROXIMATE) MATCHING
// ----------------------------------------------------------------------------

void SearchStrategy::matchWithSearches(const string& seq, const length_t k,
                                       Counters& counters, Occurrences& occs,
                                       const length_t minDistance) {

    // calculate number of parts
    length_t numParts = calculateNumParts(k);

    // create ranges for exact matches of parts
    vector<SARangePair> exactMatchRanges(numParts);

    // create a vector for the parts
    vector<Substring> parts;

    // partition the read
    count++;
    partition(seq, parts, numParts, k, exactMatchRanges, counters);

    if (parts.empty()) {
        // splitting up was not viable -> just search the entire pattern
        if (name != "Naive backtracking" && numParts > 1) {
            stringstream ss;
            ss << "Warning: The pattern size (" << seq.size()
               << ") is too short for the current number of parts (" << numParts
               << "). Using normal bidirectional search instead.\n"
               << "Pattern: " << seq;
            logger.logWarning(ss);
        }

        vector<TextOcc> textOccs;

        (index.*naiveMatchPtr)(seq, k, counters, textOccs);

        occs.addTextOccs(textOccs);
        return;
    }

#ifndef RUN_LENGTH_COMPRESSION // No in-text verification in B-move
    // A) In-text verification for parts that have low number of exact
    // matches
    for (uint16_t i = 0; i < numParts; i++) {
        size_t width = exactMatchRanges[i].width();
        if (width != 0 && width <= index.getSwitchPoint()) {
            // do in-text verification on this part
            const auto& part = parts[i];
            FMOcc startMatch(exactMatchRanges[i], 0, part.size(),
                             FORWARD_STRAND,
                             FIRST_IN_PAIR); // strand and pairStatus of
                                             // startMatch do not matter
            (this->*inTextVerificationPtr)(startMatch, part.begin(), k, parts,
                                           occs, counters, minDistance);
        }
    }
#endif

    // B) do the searches from search schemes
    // create the searches
    const vector<Search>& searches = createSearches(k, exactMatchRanges);

    // reserve stacks and matrices for each part
    index.reserveStacks(numParts,
                        seq.length()); // reserve stacks for each part
    // create the bit-parallel alignment matrices
    index.resetMatrices(parts.size()); // reset the alignment matrix that will
                                       // be (possibly) used for each part

    for (const Search& s : searches) {
        doRecSearch(s, parts, occs, exactMatchRanges, counters);
    }
}

void SearchStrategy::matchApproxAllMap(ReadBundle& bundle, length_t maxED,
                                       Counters& counters,
                                       vector<TextOcc>& result) {

    if (maxED == 0) {

        index.setIndexInMode(FORWARD_STRAND, FIRST_IN_PAIR);
        index.exactMatchesOutput(bundle.getRead(), counters, result);
        index.setIndexInMode(REVERSE_C_STRAND, FIRST_IN_PAIR);
        index.exactMatchesOutput(bundle.getRevComp(), counters, result);

        if (generateSAMLines) {
            generateOutputSingleEnd(result, bundle, counters, 0);
        }
        return;
    }

    // The occurrences in the text and index
    Occurrences occ;

    // reset all in-text verification matrices
    index.resetFullReadMatrices();

    // set the index in forward mode with the correct sequence
    index.setIndexInMode(FORWARD_STRAND);
    // map with the searches
    matchWithSearches(bundle.getRead(), maxED, counters, occ);

    // set the index in reverse complement mode with the correct sequence
    index.setIndexInMode(REVERSE_C_STRAND);
    // map with the searches
    matchWithSearches(bundle.getRevComp(), maxED, counters, occ);

    // C) get all non-redundant matches mapped to the text
    result = (this->*filterPtr)(occ, maxED, counters, bundle);

    // D) generate sam output
    if (generateSAMLines) {
        generateOutputSingleEnd(result, bundle, counters, maxED);
    }
}

void SearchStrategy::checkAlignments(OccVector& occVector, uint32_t& best,
                                     uint32_t l, Counters& counters,
                                     uint32_t cutOff, const string& seq) const {
    std::vector<TextOcc> trimmedOccs = {};
    std::vector<TextOcc> assignedOccs = {};
    for (length_t i = 0; i < occVector[l].second.size(); i++) {
        auto& occ = occVector[l].second[i];

        if (!occ.hasCigar()) {
            index.generateCIGAR(occ, counters, seq);
        }

        auto seqFound = assignSequence(occ, counters, cutOff, seq);
        if (seqFound != FOUND) {
            if (seqFound == FOUND_WITH_TRIMMING && occ.getDistance() > l) {
                trimmedOccs.emplace_back(std::move(occ));
            }
        } else {
            assignedOccs.emplace_back(std::move(occ));
            if (l < best) {
                best = l;
            }
        }
    }
    occVector[l].second = std::move(assignedOccs);
    // the occurrences found with trimming should be added
    // to the correct vector according to their new distance
    for (auto& occ : trimmedOccs) {
        occ.removeTrimmingLabel();
        occVector[occ.getDistance()].second.emplace_back(std::move(occ));
    }
}

vector<TextOcc> SearchStrategy::combineOccVectors(OccVector& ovFW,
                                                  OccVector& ovRC,
                                                  length_t best, length_t max) {
    assert(ovRC.size() == ovFW.size());
    assert(max < ovFW.size());

    auto compare = [](const TextOcc& a, const TextOcc& b) {
        return a.getAssignedSequenceID() < b.getAssignedSequenceID() ||
               (a.getAssignedSequenceID() == b.getAssignedSequenceID() &&
                a.getRange().getBegin() < b.getRange().getBegin());
    };
    auto equal = [](const TextOcc& a, const TextOcc& b) {
        return a.getAssignedSequenceID() == b.getAssignedSequenceID() &&
               a.getRange().getBegin() == b.getRange().getBegin();
    };

    std::vector<TextOcc> matches;
    matches.reserve(numElements(ovFW) + numElements(ovRC));

    for (length_t i = best; i <= max; i++) {
// some occurrences might be represented more than once
// sort the occurrences first on assigned sequence name and then on
// begin position in the text
#ifdef DEVELOPER_MODE
        std::stable_sort(ovFW[i].second.begin(), ovFW[i].second.end(), compare);
        std::stable_sort(ovRC[i].second.begin(), ovRC[i].second.end(), compare);
#else
        std::sort(ovFW[i].second.begin(), ovFW[i].second.end(), compare);
        std::sort(ovRC[i].second.begin(), ovRC[i].second.end(), compare);
#endif

        // remove duplicates
        // Remove duplicates from ovFW[i].second
        auto fw_end =
            std::unique(ovFW[i].second.begin(), ovFW[i].second.end(), equal);
        ovFW[i].second.erase(fw_end, ovFW[i].second.end());

        // Remove duplicates from ovRC[i].second
        auto rc_end =
            std::unique(ovRC[i].second.begin(), ovRC[i].second.end(), equal);
        ovRC[i].second.erase(rc_end, ovRC[i].second.end());

        matches.insert(matches.end(),
                       make_move_iterator(ovFW[i].second.begin()),
                       make_move_iterator(ovFW[i].second.end()));
        matches.insert(matches.end(),
                       make_move_iterator(ovRC[i].second.begin()),
                       make_move_iterator(ovRC[i].second.end()));
    }

    return matches;
}

bool SearchStrategy::findBestAlignments(const ReadBundle& bundle,
                                        OccVector& ovFW, OccVector& ovRC,
                                        Counters& counters, uint32_t x,
                                        uint32_t& best, PairStatus status) {

    length_t cutOff = ovFW.size() - 1;
    best = cutOff + 1;
    bool bestFound = false;

    // get the read and the revC from the bundle
    const auto& read = bundle.getRead();
    const auto& revC = bundle.getRevComp();

    if (x == 0) {
        // start with exact match, this does not need an in-text
        // verification matrix as an exact match has no insertions or
        // deletions
        if (!ovFW[0].first) {

            index.setIndexInMode(FORWARD_STRAND, status);
            index.exactMatchesOutput(bundle.getRead(), counters,
                                     ovFW[0].second);
            ovFW[0].first = true;
        }
        if (!ovRC[0].first) {
            index.setIndexInMode(REVERSE_C_STRAND, status);
            index.exactMatchesOutput(bundle.getRevComp(), counters,
                                     ovRC[0].second);
            ovRC[0].first = true; // processing for 0 finished
        }

        if (!ovFW[0].second.empty() || !ovRC[0].second.empty()) {
            checkAlignments(ovFW, best, 0, counters, cutOff, read);
            checkAlignments(ovRC, best, 0, counters, cutOff, revC);
            if (best == 0) {
                bestFound = true;
            }
        }
    }

    // moving to approximate pattern matching if necessary
    uint32_t maxED = (best == 0) ? x : cutOff;
    uint32_t prevK = 0;

    // define a function that checks whether a match was found in this stratum
    auto hasUpdate = [&](const string& seq, PairStatus status, Strand strand,
                         uint32_t k, OccVector& vector, Counters& counters) {
        if (vector[k].first) {
            return !vector[k].second.empty();
        }
        return processSeq(seq, status, strand, k, vector, counters);
    };

    for (uint32_t k = std::max(x, (uint32_t)1); k <= maxED;) {
        bool update = false;
        // start with forward direction
        update |= hasUpdate(read, status, FORWARD_STRAND, k, ovFW, counters);
        // then reverse complement
        update |= hasUpdate(revC, status, REVERSE_C_STRAND, k, ovRC, counters);

        if (update) {
            // new matches were found, we need to process them to see if
            // they need trimming
            for (length_t l = prevK + 1; l <= std::min(k, best + x); l++) {
                checkAlignments(ovFW, best, l, counters, maxED, read);
                checkAlignments(ovRC, best, l, counters, maxED, revC);
            }
        }
        if (bestFound) {
            break; // this was the last iteration
        }
        if (update && best < cutOff + 1) {
            bestFound = true;
            if (x == 0) {
                break;
            }
            prevK = k, k = std::min(best + x,
                                    maxED); // check the final x strata
        } else {
            // continue to next stratum
            if (k == maxED)
                break;
            uint32_t step = (k < 5) ? 2 : 4;
            prevK = k;
            k = std::min(k + x + step, maxED);
        }
    }

    return bestFound;
}

void SearchStrategy::matchApproxBestPlusX(ReadBundle& bundle, length_t x,
                                          Counters& counters,
                                          const length_t minIdentity,
                                          vector<TextOcc>& result) {

    // reset the in-text verification matrices from previous read
    index.resetFullReadMatrices();

    length_t cutOff = getMaxED(minIdentity, bundle.size());
    OccVector occVectorFW(cutOff + 1), occVectorRC(cutOff + 1);

    uint32_t best;
    bool bestFound =
        findBestAlignments(bundle, occVectorFW, occVectorRC, counters, x, best);

    if (!bestFound) {
        if (unmappedSAM) {
            result.emplace_back(createUnmappedRecordSE(bundle));
        }
        return;
    }

    // find out the number of hits for best
    uint32_t nHits =
        occVectorFW[best].second.size() + occVectorRC[best].second.size();

    result = combineOccVectors(occVectorFW, occVectorRC, best,
                               std::min(best + x, cutOff));

    (this->*generateOutputSEPtr)(bundle, nHits, best, result, counters);
}

vector<PairedTextOccs> SearchStrategy::matchApproxPairedEndAll(
    ReadPair& pair, length_t maxEDOrIdentity, length_t maxFragSize,
    length_t minFragSize, Counters& counters,
    vector<TextOcc>& unpairedOccurrences) {

    length_t maxED = maxEDOrIdentity;     // maximal distance
    vector<PairedTextOccs> pairedMatches; // output vector

    // vector of occurrences per distance per strand and read
    BoolAndVector matches1forward, matches2forward, matches2revCompl,
        matches1revCompl;

    // reset all in-text verification matrices from the previous pair
    index.resetFullReadMatrices();

    (this->*processOriAllPtr)(pair, matches1forward, matches1revCompl,
                              matches2forward, matches2revCompl, pairedMatches,
                              maxFragSize, minFragSize, maxED, counters);

    if (pairedMatches.empty()) {
        pairDiscordantly(matches1forward, matches1revCompl, matches2forward,
                         matches2revCompl, pairedMatches, pair, maxED, counters,
                         unpairedOccurrences);
    }

    // now we need to generate SAM lines for the paired matches
    generateSAMPairedEnd(pairedMatches, pair);

    return pairedMatches;
}

bool SearchStrategy::processSeq(const string& seq, const PairStatus status,
                                const Strand strand, const length_t maxDist,
                                OccVector& vector, Counters& counters) {
    assert(maxDist < vector.size());
    // if maxDist has been processed then all distances below it have
    // also been processed
    if (!vector[maxDist].first) {

        // find the minimum distance (first distance d for which
        // vector[d].first == false)
        auto it =
            find_if(vector.begin(), vector.end(),
                    [](const BoolAndVector& pair) { return !pair.first; });
        length_t minD = min((length_t)distance(vector.begin(), it), maxDist);

        auto occs = mapRead(seq, maxDist, counters, status, strand, minD);

        // the mapped occurrences can be between the min and max
        // distance they need to be split up according to their
        // distances and added to the correct vector
        for (auto& occ : occs) {
            vector[occ.getDistance()].second.emplace_back(std::move(occ));
        }
        // set distance processed to true
        for (length_t i = minD; i <= maxDist; i++) {
            vector[i].first = true;
        }
    }
    // return true if any of the distances between 0 and maxDist have an
    // occurrence
    bool found =
        any_of(vector.begin(), vector.begin() + maxDist + 1,
               [](const BoolAndVector& pair) { return !pair.second.empty(); });
    return found;
}

void SearchStrategy::handleTrimmedOccs(vector<length_t>& trimmedIds,
                                       const length_t oDist, OccVector& v) {

// sort the vector by id
#ifdef DEVELOPER_MODE
    stable_sort(trimmedIds.begin(), trimmedIds.end());
#else
    sort(trimmedIds.begin(), trimmedIds.end());
#endif
    for (auto idIt = trimmedIds.rbegin(); idIt != trimmedIds.rend(); ++idIt) {
        auto id = *idIt;
        auto& occ = v[oDist].second[id];
        occ.removeTrimmingLabel();
        // insert occ in the correct vector and remove from old vector
        auto& tVector = v[occ.getDistance()].second;
        auto point = lower_bound(tVector.begin(), tVector.end(), occ);
        tVector.insert(point, std::move(occ));
        v[oDist].second.erase(v[oDist].second.begin() + id);
    }
}

void SearchStrategy::processComb(
    const string& uSeq, const string& dSeq, Counters& counters,
    PairStatus uStatus, Strand uStrand, PairStatus dStatus, Strand dStrand,
    OccVector& uVector, OccVector& dVector, vector<PairedTextOccs>& pairs,
    int maxFragSize, int minFragSize, length_t totDist, length_t minTotDist) {
    assert(!uVector.empty());
    assert(!dVector.empty());

    auto findFirstPosDist = [](const OccVector& vector) {
        return find_if(vector.begin(), vector.end(),
                       [](const BoolAndVector& pair) {
                           return !pair.second.empty() || !pair.first;
                       });
    };

    auto processRead = [&](const string& seq, PairStatus status, Strand strand,
                           OccVector& vector, const length_t& max,
                           length_t& maxOther) {
        if (!processSeq(seq, status, strand, max, vector, counters)) {
            return false; // no occurrences found for read
        }
        // Recalculate minDist and maxOther
        length_t minD = distance(vector.cbegin(), findFirstPosDist(vector));
        maxOther = min(totDist - minD, maxOther);
        return true;
    };

    // find the max possible distance for up and downstream
    // the max distance is the total distance minus the minimum distance
    // for the other stream

    length_t minDDist = distance(dVector.cbegin(), findFirstPosDist(dVector));
    length_t minUDist = distance(uVector.cbegin(), findFirstPosDist(uVector));

    length_t maxUp = min(totDist - minDDist, (length_t)uVector.size() - 1);
    length_t maxDown = min(totDist - minUDist, (length_t)dVector.size() - 1);

    if (maxUp <= maxDown) {
        if (!(processRead(uSeq, uStatus, uStrand, uVector, maxUp, maxDown) &&
              processRead(dSeq, dStatus, dStrand, dVector, maxDown, maxUp)))
            return;

    } else {
        if (!(processRead(dSeq, dStatus, dStrand, dVector, maxDown, maxUp) &&
              processRead(uSeq, uStatus, uStrand, uVector, maxUp, maxDown)))
            return;
    }

    for (length_t dist = minUDist + minDDist; dist <= totDist; dist++) {
        for (length_t uDist = minUDist; uDist <= min(maxUp, dist); uDist++) {

            length_t dDist = dist - uDist;
            if (dDist > maxDown || dDist < minDDist) {
                continue;
            }
            set<length_t> uTrimmedIds;
            set<length_t> dTrimmedIds;
            // pair the occurrences that have a distance of uDist +
            // dDist
            pairOccurrencesForBestMapping(
                uVector[uDist].second, dVector[dDist].second, pairs,
                maxFragSize, minFragSize, maxUp, maxDown, counters, uTrimmedIds,
                dTrimmedIds, uSeq, dSeq);

            // convert the sets to vectors
            vector<length_t> uTrimmedIdsVec(uTrimmedIds.begin(),
                                            uTrimmedIds.end());
            vector<length_t> dTrimmedIdsVec(dTrimmedIds.begin(),
                                            dTrimmedIds.end());

            // the occurrences found with trimming should be added to
            // the correct vector according to their new distance and
            // removed from the old vector
            handleTrimmedOccs(uTrimmedIdsVec, uDist, uVector);
            handleTrimmedOccs(dTrimmedIdsVec, dDist, dVector);
        }
        if (!pairs.empty()) {
            // no need to further explore higher strata
            return;
        }
    }
}

/**
 * Add the best pairs for the 1-2 and 2-1 combination to the pairs vector.
 * @param pairs12 the pairs for the 1-2 combination
 * @param pairs21 the pairs for the 2-1 combination
 * @param pairs the vector with the paired occurrences (can be updated)
 */
void mergeOrMovePairs(vector<PairedTextOccs>& pairs12,
                      vector<PairedTextOccs>& pairs21,
                      vector<PairedTextOccs>& pairs) {
    if (pairs12.empty() || pairs21.empty()) {
        pairs = std::move(pairs12.empty() ? pairs21 : pairs12);
        return;
    }
    length_t distance12 = pairs12.front().getDistance(),
             distance21 = pairs21.front().getDistance();
    if (distance12 <= distance21) {
        pairs = std::move(pairs12); // pairs12 is best or ex-aequo
        if (distance12 == distance21) {
            pairs.insert(pairs.end(), make_move_iterator(pairs21.begin()),
                         make_move_iterator(pairs21.end()));
        }
        return;
    }
    pairs = std::move(pairs21);
    return;
}

void SearchStrategy::processCombFF(ReadPair& pair, OccVector& matches1,
                                   OccVector& matches2, OccVector& matchesRC1,
                                   OccVector& matchesRC2,
                                   vector<PairedTextOccs>& pairs,
                                   int maxFragSize, int minFragSize,
                                   length_t totDist, length_t minTotDist,
                                   Counters& counters) {

    vector<PairedTextOccs> pairsFF; // the pairs that are forward-forward
    vector<PairedTextOccs> pairsRR; // the pairs that are reverseC-reverseC

    auto processFF = [&]() {
        processComb(pair.getBundle1().getRead(), pair.getBundle2().getRead(),
                    counters, FIRST_IN_PAIR, FORWARD_STRAND, SECOND_IN_PAIR,
                    FORWARD_STRAND, matches1, matches2, pairsFF, maxFragSize,
                    minFragSize, totDist, minTotDist);
    };
    auto processRR = [&]() {
        processComb(pair.getBundle2().getRevComp(),
                    pair.getBundle1().getRevComp(), counters, SECOND_IN_PAIR,
                    REVERSE_C_STRAND, FIRST_IN_PAIR, REVERSE_C_STRAND,
                    matchesRC2, matchesRC1, pairsRR, maxFragSize, minFragSize,
                    totDist, minTotDist);
    };

    // choose ordering of what combination to process first based on of
    // the comb already has a match for dist < minTotDist
    if (any_of(matches1.begin(), matches1.end(),
               [](BoolAndVector& bv) { return !bv.second.empty(); }) ||
        any_of(matches2.begin(), matches2.end(),
               [](BoolAndVector& bv) { return !bv.second.empty(); })) {

        processFF();
        totDist = (pairsFF.empty()) ? totDist : pairsFF.front().getDistance();
        processRR();
    } else {
        processRR();
        totDist = (pairsRR.empty()) ? totDist : pairsRR.front().getDistance();
        processFF();
    }
    mergeOrMovePairs(pairsFF, pairsRR, pairs);
}

void SearchStrategy::processCombRF(ReadPair& pair, OccVector& matches1,
                                   OccVector& matches2, OccVector& matchesRC1,
                                   OccVector& matchesRC2,
                                   vector<PairedTextOccs>& pairs,
                                   int maxFragSize, int minFragSize,
                                   length_t totDist, length_t minTotDist,
                                   Counters& counters) {
    vector<PairedTextOccs> pairs12; // the pairs that are 1-2
    vector<PairedTextOccs> pairs21; // the pairs that are 2-1
    auto process12 = [&]() {
        processComb(pair.getBundle1().getRevComp(), pair.getBundle2().getRead(),
                    counters, FIRST_IN_PAIR, REVERSE_C_STRAND, SECOND_IN_PAIR,
                    FORWARD_STRAND, matchesRC1, matches2, pairs12, maxFragSize,
                    minFragSize, totDist, minTotDist);
    };

    auto process21 = [&]() {
        processComb(pair.getBundle2().getRevComp(), pair.getBundle1().getRead(),
                    counters, SECOND_IN_PAIR, REVERSE_C_STRAND, FIRST_IN_PAIR,
                    FORWARD_STRAND, matchesRC2, matches1, pairs21, maxFragSize,
                    minFragSize, totDist, minTotDist);
    };

    if (any_of(matchesRC1.begin(), matchesRC1.begin() + minTotDist,
               [](BoolAndVector& bv) { return !bv.second.empty(); }) ||
        any_of(matches2.begin(), matches2.begin() + minTotDist,
               [](BoolAndVector& bv) { return !bv.second.empty(); })) {
        process12();
        totDist = (pairs12.empty()) ? totDist : pairs12.front().getDistance();
        process21();
    } else {
        process21();
        totDist = (pairs21.empty()) ? totDist : pairs21.front().getDistance();
        process12();
    }

    mergeOrMovePairs(pairs12, pairs21, pairs);
}

void SearchStrategy::processCombFR(ReadPair& pair, OccVector& matches1,
                                   OccVector& matches2, OccVector& matchesRC1,
                                   OccVector& matchesRC2,
                                   vector<PairedTextOccs>& pairs,
                                   int maxFragSize, int minFragSize,
                                   length_t totDist, length_t minTotDist,
                                   Counters& counters) {
    vector<PairedTextOccs> pairs12; // the pairs that are 1-2
    vector<PairedTextOccs> pairs21; // the pairs that are 2-1

    auto process12 = [&]() {
        processComb(pair.getBundle1().getRead(), pair.getBundle2().getRevComp(),
                    counters, FIRST_IN_PAIR, FORWARD_STRAND, SECOND_IN_PAIR,
                    REVERSE_C_STRAND, matches1, matchesRC2, pairs12,
                    maxFragSize, minFragSize, totDist, minTotDist);
    };

    auto process21 = [&]() {
        processComb(pair.getBundle2().getRead(), pair.getBundle1().getRevComp(),
                    counters, SECOND_IN_PAIR, FORWARD_STRAND, FIRST_IN_PAIR,
                    REVERSE_C_STRAND, matches2, matchesRC1, pairs21,
                    maxFragSize, minFragSize, totDist, minTotDist);
    };

    if (any_of(matches1.begin(), matches1.begin() + minTotDist,
               [](BoolAndVector& bv) { return !bv.second.empty(); }) ||
        any_of(matchesRC2.begin(), matchesRC2.begin() + minTotDist,
               [](BoolAndVector& bv) { return !bv.second.empty(); })) {
        process12();
        totDist = (pairs12.empty()) ? totDist : pairs12.front().getDistance();
        process21();
    } else {
        process21();
        totDist = (pairs21.empty()) ? totDist : pairs21.front().getDistance();
        process12();
    }

    mergeOrMovePairs(pairs12, pairs21, pairs);
}

void SearchStrategy::addSingleEndedForBest(vector<TextOcc>& matchesSE1,
                                           vector<TextOcc>& matchesSE2,
                                           OccVector& fw1, OccVector& rc1,
                                           OccVector& fw2, OccVector& rc2,
                                           bool read2done) {

    for (length_t i = 0; i < matchesSE1.size(); i++) {
        auto occ = std::move(matchesSE1[i]);
        auto& vector = (occ.getStrand() == FORWARD_STRAND) ? fw1 : rc1;
        vector[occ.getDistance()].second.emplace_back(std::move(occ));
    }
    for (length_t i = 0; i < matchesSE2.size(); i++) {
        auto occ = std::move(matchesSE2[i]);
        occ.setPairStatus(SECOND_IN_PAIR); // set flag for  second read
        auto& vector = (occ.getStrand() == FORWARD_STRAND) ? fw2 : rc2;
        vector[occ.getDistance()].second.emplace_back(std::move(occ));
    }
    for (length_t i = 0; i < fw1.size(); i++) {
        fw1[i].first = rc1[i].first = true;
    }
    for (length_t i = 0; i < fw2.size(); i++) {
        fw2[i].first = rc2[i].first = read2done;
    }
}

vector<PairedTextOccs> SearchStrategy::matchApproxPairedEndBestPlusX(
    ReadPair& pair, length_t minIdentity, length_t maxFragSize,
    length_t minFragSize, Counters& counters, vector<TextOcc>& unpairedOcc,
    length_t x, vector<TextOcc>& matchesSE1, vector<TextOcc>& matchesSE2,
    bool singleEnd, bool read2done) {

    length_t cutOff1 = getMaxED(minIdentity, pair.getRead1().size());
    length_t cutOff2 = getMaxED(minIdentity, pair.getRead2().size());

    length_t best = cutOff1 + cutOff2 + 1;

    vector<PairedTextOccs> pairs; // return value

    // create for the read1 and revComp1 a vector of vector of TextOcc
    // of size cutOff1 + 1 (one index per valid distance),  do the same
    // for read2 and revComp2
    OccVector fw1(cutOff1 + 1), rc1(cutOff1 + 1), fw2(cutOff2 + 1),
        rc2(cutOff2 + 1);
    if (singleEnd) {
        addSingleEndedForBest(matchesSE1, matchesSE2, fw1, rc1, fw2, rc2,
                              read2done);
    }

    // reset all in-text verification matrices from previous pair
    index.resetFullReadMatrices();

    // the minimal distance that has not been explored
    length_t minDistNotExplored = 0;

    // try to find a pair that matches exactly according to the
    // orientation
    if (x == 0) {

        (this->*processOriBestPtr)(pair, fw1, fw2, rc1, rc2, pairs, maxFragSize,
                                   minFragSize, 0, 0, counters);
        minDistNotExplored = 1; // all distances below 1 have been explored
    }
    bool bestFound = false;
    if (!pairs.empty()) {
        best = 0;
        bestFound = true;
    }

    length_t maxStratum = (best == 0) ? x : cutOff1 + cutOff2;

    for (length_t k = max(x, (length_t)1); k <= maxStratum;) {

        (this->*processOriBestPtr)(pair, fw1, fw2, rc1, rc2, pairs, maxFragSize,
                                   minFragSize, k, minDistNotExplored,
                                   counters);

        if (!bestFound) {
            if (!pairs.empty()) {
                // A best pair has been found
                best = k;
                bestFound = true;
                maxStratum = min(best + x, cutOff1 + cutOff2);
                minDistNotExplored = k + 1;
                if (x == 0) {
                    break;
                }

                k = maxStratum;

            } else {
                if (k == maxStratum) {
                    // no best found and at the end of the loop
                    break;
                }
                // increase the stratum
                length_t step = (k < 6) ? 2 : 4;
                k = min(maxStratum, k + x + step);
            }
        } else {
            // best was found already, this was best+x so break the loop
            break;
        }
    }

    if (pairs.empty()) {

        // no pairs found, try discordant pairs
        pairDiscordantlyBest(fw1, rc1, fw2, rc2, pairs, pair, counters,
                             unpairedOcc, x);
    }

    generateSAMPairedEnd(pairs, pair);
    return pairs;
}

void SearchStrategy::doRecSearch(const Search& s, vector<Substring>& parts,
                                 Occurrences& occ,
                                 const vector<SARangePair>& exactMatchRanges,
                                 Counters& counters) const {

    if (s.getUpperBound(0) > 0) {
        // first part is allowed an error so start with an empty match
        s.setDirectionsInParts(parts);

        FMOcc startMatch = FMOcc(index.getEmptyStringFMPos(), 0);
        (this->*startIdxPtr)(s, startMatch, occ, parts, counters, 0);
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
        uint16_t partInSearch = 1;
        length_t exactLength = parts[first].size();

        while (s.getUpperBound(partInSearch) == 0) {
            // extend the exact match
            index.setDirection(s.getDirection(partInSearch),
                               s.isUnidirectionalBackwards(partInSearch));
            const auto& part = parts[s.getPart(partInSearch)];

            startRange =
                index.matchStringBidirectionally(part, startRange, counters);
            if (startRange.empty()) {
                return;
            }

            exactLength += part.size();
            partInSearch++;
        }

        // Create a match corresponding to the exact match
        FMOcc startMatch = FMOcc(startRange, 0, exactLength);

#ifdef RUN_LENGTH_COMPRESSION
        // TODO: do this once for all parts in the beginning
        // add the exact match string to the start match

        length_t lPartIdx = s.getLowestPartProcessedBefore(partInSearch);
        length_t hPartIdx = s.getHighestPartProcessedBefore(partInSearch);

        std::string str =
            Substring(parts[lPartIdx].getText(), parts[lPartIdx].begin(),
                      parts[hPartIdx].end())
                .tostring();

        std::vector<char> strVector(str.begin(), str.end());

        if (s.getDirection(partInSearch) == BACKWARD) {
            // reverse the string if the direction is backward
            // so that we can append to the vector in the phase
            std::reverse(strVector.begin(), strVector.end());
        }

        startMatch.setMatchedStr(strVector);
#endif

        // Start approximate matching in the index
        (this->*startIdxPtr)(s, startMatch, occ, parts, counters, partInSearch);
    }
}

// ----------------------------------------------------------------------------
// PAIRING
// ----------------------------------------------------------------------------

/**
 * @brief Find the first occurrence in a sorted vector of TextOcc that is after
 * a specified position
 * @param occs The vector of TextOcc
 * @param position The position to search for
 * @return An iterator to the first occurrence that is after the specified
 * position
 */
vector<TextOcc>::iterator findFirstOccAfter(vector<TextOcc>& occs,
                                            length_t position) {

    assert(is_sorted(occs.begin(), occs.end(),
                     [](const TextOcc& a, const TextOcc& b) {
                         return a.getIndexBegin() < b.getIndexBegin();
                     }));
    return lower_bound(occs.begin(), occs.end(), position,
                       [](const TextOcc& occ, length_t pos) {
                           return occ.getIndexBegin() < pos;
                       });
}

void SearchStrategy::pairOccurrences(vector<TextOcc>& upstreamOccs,
                                     vector<TextOcc>& downStreamOccs,
                                     vector<PairedTextOccs>& pairs,
                                     length_t maxFragSize, length_t minFragSize,
                                     length_t maxED, Counters& counters,
                                     const std::string& uSeq,
                                     const std::string& dSeq) const {
    // TODO (maybe): if upStreamOccs.size() >> downStreamOccs.size()
    // then loop over downStreamOccs and find the upstream matches for
    // each

    if (upstreamOccs.size() == 0 || downStreamOccs.size() == 0) {
        return;
    }

    // loop over all matches of the upstream read
    for (auto& upStreamOcc : upstreamOccs) {
        length_t upstreamPos = upStreamOcc.getIndexBegin();

        // find the index of first match of the downstream read that is
        // after the upstream read
        auto it = findFirstOccAfter(downStreamOccs, upstreamPos);

        // loop over all matches of downstream read that are after the
        // upstream read
        for (; it != downStreamOccs.end(); it++) {
            // check if the matches are within the fragment size
            // and have the correct orientation

            length_t fragSize = (it->getIndexEnd()) - upstreamPos;

            if (fragSize <= maxFragSize && fragSize >= minFragSize) {
                // found a pair
                // assign sequence if not assigned yet

                if (assignSequenceAndCIGAR(upStreamOcc, counters, maxED,
                                           uSeq) == NOT_FOUND) {
                    break; // no sequence could be assigned, break inner loop
                }

                if (assignSequenceAndCIGAR(*it, counters, maxED, dSeq) ==
                    NOT_FOUND) {
                    continue; // no sequence could be assigned
                }

                if (upStreamOcc.getAssignedSequenceID() !=
                    it->getAssignedSequenceID()) {
                    // the assigned sequences are not the same, so no
                    // pair can be formed (discordant pairs are handled
                    // later on)
                    continue;
                }

                pairs.push_back({upStreamOcc, *it});

            } else if (fragSize > maxFragSize) {
                // since the matches are sorted, we can stop the
                // loop
                break;
            }
        }
    }
}

vector<PairedTextOccs> SearchStrategy::pairSingleEndedMatchesAll(
    vector<TextOcc>& matches1, vector<TextOcc>& matches2, length_t maxFragSize,
    length_t minFragSize, Counters& counters, ReadPair& readPair,
    length_t maxED, vector<TextOcc>& unpairedOccurrences, bool read2done) {

    vector<PairedTextOccs> pairedMatches; // return value

    // we need to split up matches1 and matches2 in forward and revCompl
    vector<TextOcc> fw1, rc1, fw2, rc2;
    // iterate over matches1 and move into correct vector
    for (auto& match : matches1) {
        auto& vector = match.isRevCompl() ? rc1 : fw1;
        vector.emplace_back(std::move(match));
    }
    // iterate over matches2 and move into correct vector
    for (auto& match : matches2) {
        // set flag for matches2 as second read in pair
        match.setPairStatus(SECOND_IN_PAIR);
        auto& vector = match.isRevCompl() ? rc2 : fw2;
        vector.emplace_back(std::move(match));
    }

// sort the matches on their begin position (needed for pairing)
#ifdef DEVELOPER_MODE
    std::stable_sort(fw1.begin(), fw1.end());
    std::stable_sort(fw2.begin(), fw2.end());
    std::stable_sort(rc1.begin(), rc1.end());
    std::stable_sort(rc2.begin(), rc2.end());
#else
    sort(fw1.begin(), fw1.end());
    sort(fw2.begin(), fw2.end());
    sort(rc1.begin(), rc1.end());
    sort(rc2.begin(), rc2.end());
#endif

    // reset the in-text verification matrices from a previous pair
    index.resetFullReadMatrices();

    // create BoolAndVector for each vector with bool set to true
    BoolAndVector fw1bv = {true, std::move(fw1)},
                  rc1bv = {true, std::move(rc1)},
                  fw2bv = {read2done, std::move(fw2)},
                  rc2bv = {read2done, std::move(rc2)};

    (this->*processOriAllPtr)(readPair, fw1bv, rc1bv, fw2bv, rc2bv,
                              pairedMatches, maxFragSize, minFragSize, maxED,
                              counters);

    if (pairedMatches.empty()) {
        pairDiscordantly(fw1bv, rc1bv, fw2bv, rc2bv, pairedMatches, readPair,
                         maxED, counters, unpairedOccurrences);
    }
    generateSAMPairedEnd(pairedMatches, readPair);
    return pairedMatches;
}

void SearchStrategy::addUnpairedMatches(vector<TextOcc>& allMatches,
                                        ReadBundle& bundle, PairStatus status,
                                        length_t maxED, Counters& counters,
                                        vector<TextOcc>& unpairedOcc) {
    vector<TextOcc> temp; // a temp vector with the unpaired matches

    // move the elements of forward and revComp to temp if a
    // sequence can be assigned
    for (auto& occ : allMatches) {

        if (!occ.hasCigar()) {
            const auto& seq = bundle.getSequence(occ.getStrand());
            index.generateCIGAR(occ, counters, seq);
        }

        const auto& seq = bundle.getSequence(occ.getStrand());
        if (assignSequence(occ, counters, maxED, seq) != NOT_FOUND) {
            temp.emplace_back(std::move(occ));
        }
    }

    if (temp.empty()) {
        // probably won't happen but you never know, just create an
        // unmapped occurrence
        if (unmappedSAM) {
            unpairedOcc.emplace_back(
                TextOcc::createUnmappedSAMOccurrencePE(bundle, status));
        }
        return;
    }

// sort on distance score
#ifdef DEVELOPER_MODE
    std::stable_sort(
#else
    sort(
#endif
        temp.begin(), temp.end(), [](const TextOcc& a, const TextOcc& b) {
            return a.getDistance() < b.getDistance();
        });

    // Get best score and number of times it occurs
    uint32_t best = temp.front().getDistance();
    // count the number of occurrences with best score
    auto bestCount =
        count_if(temp.begin(), temp.end(), [best](const TextOcc& occ) {
            return occ.getDistance() == best;
        });

    // generate SAM lines
    bool firstOcc = true;
    for (auto& occ : temp) {
        occ.generateSAMUnpaired(bundle, bestCount, best, firstOcc,
                                index.getSeqNames());
        firstOcc = false;
    }

    // add the temp occurrences to the unpairedOcc vector
    unpairedOcc.insert(unpairedOcc.end(), std::make_move_iterator(temp.begin()),
                       std::make_move_iterator(temp.end()));
}

void SearchStrategy::addOneUnmapped(ReadPair& reads, vector<TextOcc>& matches1,
                                    vector<TextOcc>& matches2,
                                    vector<PairedTextOccs>& pairs,
                                    length_t maxED, Counters& counters) {
    // Define lambda function to create paired occurrences with one
    // mapped read
    auto createPairedOccurrences = [&](vector<TextOcc>& occurrences,
                                       const ReadBundle& bundle1,
                                       const ReadBundle& bundle2,
                                       const bool isFirstMapped,
                                       vector<PairedTextOccs>& pairs) {
        const auto& mappedBundle = (isFirstMapped) ? bundle1 : bundle2;

        for (auto& occurrence : occurrences) {
            // assign sequence and generate CIGAR string
            const auto& seq = mappedBundle.getSequence(occurrence.getStrand());

            if (!occurrence.hasCigar()) {
                index.generateCIGAR(occurrence, counters, seq);
            }

            if (assignSequence(occurrence, counters, maxED, seq) == NOT_FOUND) {
                continue;
            }

            PairStatus unmappedReadStatus =
                isFirstMapped ? SECOND_IN_PAIR : FIRST_IN_PAIR;

            auto unMapped = TextOcc::createUnmappedSAMOccurrencePE(
                isFirstMapped ? bundle2 : bundle1, unmappedReadStatus, true,
                occurrence.getStrand());
            PairedTextOccs pair = {occurrence, std::move(unMapped), 0};
            pairs.emplace_back(std::move(pair));
        }
    };

    bool firstMapped = matches1.size() > 0;
    assert(firstMapped == (matches2.size() == 0));

    // report pairs where the mapped one is upstream and
    // mate is empty (unmapped) downstream
    if (firstMapped) {
        createPairedOccurrences(matches1, reads.getBundle1(),
                                reads.getBundle2(), firstMapped, pairs);
    } else {
        createPairedOccurrences(matches2, reads.getBundle1(),
                                reads.getBundle2(), firstMapped, pairs);
    }
    if (pairs.empty()) {
        // unlikely but possible that no match with assigned sequence
        // could be found
        addBothUnmapped(reads, pairs);
    }
}

void SearchStrategy::addDiscPairs(vector<TextOcc>& fw1, vector<TextOcc>& rc1,
                                  vector<TextOcc>& fw2, vector<TextOcc>& rc2,
                                  vector<PairedTextOccs>& pairs, length_t maxED,
                                  Counters& counters, ReadPair& reads,
                                  vector<TextOcc>& unpairedOcc) {
    if ((fw1.empty() && rc1.empty()) || (fw2.empty() && rc2.empty())) {
        // no pairs can be formed
        return;
    }
    auto pairOccs = [&](TextOcc& first, TextOcc& second,
                        vector<PairedTextOccs>& pairs) {
        // generate CIGAR strings if necessary
        const auto& firstSeq =
            reads.getSequence(first.getStrand(), first.getPairStatus());
        const auto& secondSeq =
            reads.getSequence(second.getStrand(), second.getPairStatus());

        if (!first.hasCigar()) {
            index.generateCIGAR(first, counters, firstSeq);
        }
        if (!second.hasCigar()) {
            index.generateCIGAR(second, counters, secondSeq);
        }

        // try to assign a sequence to both occurrences
        if (assignSequence(first, counters, maxED, firstSeq) != NOT_FOUND &&
            assignSequence(second, counters, maxED, secondSeq) != NOT_FOUND) {

            // if they are not aligned to the same reference, the
            // fragment size is 0 to mark that template length does not
            // make sense here
            bool sameRef =
                first.getAssignedSequenceID() == second.getAssignedSequenceID();

            bool firstUpstream = first.getBegin() < second.getBegin();
            uint32_t fragSize =
                (sameRef) ? (firstUpstream ? second.getEnd() - first.getBegin()
                                           : first.getEnd() - second.getBegin())
                          : 0;

            auto& upstream = firstUpstream ? first : second;
            auto& downstream = firstUpstream ? second : first;

            // create the discordant paired match
            PairedTextOccs pair = {upstream, downstream, fragSize};
            pair.setDiscordant();
            pairs.emplace_back(std::move(pair));
        }
    };

    for (auto& first : fw1) {
        for (auto& second : fw2) {
            pairOccs(first, second, pairs);
        }
        for (auto& second : rc2) {
            pairOccs(first, second, pairs);
        }
    }
    for (auto& first : rc1) {
        for (auto& second : fw2) {
            pairOccs(first, second, pairs);
        }
        for (auto& second : rc2) {
            pairOccs(first, second, pairs);
        }
    }
}

void SearchStrategy::pairDiscordantly(BoolAndVector& fw1, BoolAndVector& rc1,
                                      BoolAndVector& fw2, BoolAndVector& rc2,
                                      vector<PairedTextOccs>& pairs,
                                      ReadPair& reads, length_t maxED,
                                      Counters& counters,
                                      vector<TextOcc>& unpairedOcc) {

    // step 1: since some reads/rev_complements might not have been
    // matched yet, we need to map them

    auto map = [&](BoolAndVector& vector, const string& seq, PairStatus status,
                   Strand strand) {
        if (!vector.first) {
            vector = {true, mapRead(seq, maxED, counters, status, strand)};
        }
    };

    map(fw1, reads.getRead1(), FIRST_IN_PAIR, FORWARD_STRAND);
    map(rc1, reads.getRevComp1(), FIRST_IN_PAIR, REVERSE_C_STRAND);
    map(fw2, reads.getRead2(), SECOND_IN_PAIR, FORWARD_STRAND);
    map(rc2, reads.getRevComp2(), SECOND_IN_PAIR, REVERSE_C_STRAND);

    length_t matchesFirst = fw1.second.size() + rc1.second.size();
    length_t matchesSecond = fw2.second.size() + rc2.second.size();

    bool firstMapped = matchesFirst != 0;
    bool secondMapped = matchesSecond != 0;
    // step 2: pair discordantly if both reads are mapped
    if (discordantAllowed && firstMapped && secondMapped) {

        if (matchesFirst * matchesSecond > 10000) {
            logger.logWarning("Too many possible discordant pairs to pair "
                              "discordantly for reads " +
                              reads.getBundle1().getSeqID() + " and " +
                              reads.getBundle2().getSeqID() + ": " +
                              to_string(matchesFirst * matchesSecond));
            addUnpairedMatches(fw1.second, rc1.second, fw2.second, rc2.second,
                               reads, maxED, counters, unpairedOcc);

        } else {
            addDiscPairs(fw1.second, rc1.second, fw2.second, rc2.second, pairs,
                         maxED, counters, reads, unpairedOcc);
            if (!pairs.empty()) {
                return;
            }
        }
    }

    // step2: no discordant pairs were found or allowed, so we need to
    // report the temp mappings
    if (firstMapped && secondMapped) {
        addUnpairedMatches(fw1.second, rc1.second, fw2.second, rc2.second,
                           reads, maxED, counters, unpairedOcc);
    } else if (!firstMapped && !secondMapped) {
        addBothUnmapped(reads, pairs);
    } else {
        addOneUnmapped(reads, fw1.second, rc1.second, fw2.second, rc2.second,
                       pairs, maxED, counters);
    }
}

vector<TextOcc> SearchStrategy::findBestMapping(OccVector& fw, OccVector& rc,
                                                const ReadBundle& bundle,
                                                Counters& counters,
                                                PairStatus status, uint32_t x) {
    assert(fw.size() == rc.size());

    uint32_t best;
    bool bestFound =
        findBestAlignments(bundle, fw, rc, counters, x, best, status);

    if (bestFound) {
        uint32_t max = fw.size() - 1;
        return combineOccVectors(fw, rc, best, std::min(best + x, max));
    }
    return {};
}

void SearchStrategy::pairDiscordantlyBest(OccVector& fw1, OccVector& rc1,
                                          OccVector& fw2, OccVector& rc2,
                                          vector<PairedTextOccs>& pairs,
                                          ReadPair& reads, Counters& counters,
                                          vector<TextOcc>& unpairedOcc,
                                          length_t x) {

    assert(fw1.size() == rc1.size());
    assert(fw2.size() == rc2.size());
    length_t max1 = fw1.size() - 1;
    length_t max2 = fw2.size() - 1;
    if (discordantAllowed) {

        length_t bestStratum = fw1.size() + fw2.size() + 1;
        bool bestFound = false;

        // we need to find the pair with the best total distance score

        for (length_t i = 0; i < fw1.size() + fw2.size(); i++) {
            if (i <= max1) {
                mapStratum(fw1, i, counters, reads.getRead1(), FIRST_IN_PAIR,
                           FORWARD_STRAND);
                mapStratum(rc1, i, counters, reads.getRevComp1(), FIRST_IN_PAIR,
                           REVERSE_C_STRAND);
            }
            if (i <= max2) {
                mapStratum(fw2, i, counters, reads.getRead2(), SECOND_IN_PAIR,
                           FORWARD_STRAND);
                mapStratum(rc2, i, counters, reads.getRevComp2(),
                           SECOND_IN_PAIR, REVERSE_C_STRAND);
            }
            length_t min1 = (i > max2) ? i - max2 : 0;
            for (length_t e1 = min1; e1 <= min(i, max1); e1++) {
                length_t e2 = i - e1;

                // get all possible discordant pairs for this stratum
                addDiscPairs(fw1[e1].second, rc1[e1].second, fw2[e2].second,
                             rc2[e2].second, pairs, i, counters, reads,
                             unpairedOcc);
            }
            if (!pairs.empty()) {
                // found discordant pairs
                if (!bestFound) {
                    bestStratum = i;
                    bestFound = true;
                }
                if (i == bestStratum + x) {
                    return;
                }
            }
        }
    }
    // discordant was not allowed or no discordant pairs were found
    // we need to add the best mappings for each read as unpaired
    // step 1: find the best mappings
    auto best1 = findBestMapping(fw1, rc1, reads.getBundle1(), counters,
                                 FIRST_IN_PAIR, x);
    auto best2 = findBestMapping(fw2, rc2, reads.getBundle2(), counters,
                                 SECOND_IN_PAIR, x);

    if (best1.empty() && best2.empty()) {
        // no mappings found, add both reads as unmapped
        addBothUnmapped(reads, pairs);
    } else if (best1.empty()) {
        // only read 1 has no mappings, add read 1 as unmapped
        addOneUnmapped(reads, best1, best2, pairs, max2, counters);
    } else if (best2.empty()) {
        // only read 2 has no mappings, add read 2 as unmapped
        addOneUnmapped(reads, best1, best2, pairs, max1, counters);

    } else {
        // step 2: add the best mappings as unpaired
        addUnpairedMatches(best1, reads.getBundle1(), FIRST_IN_PAIR,
                           fw1.size() - 1, counters, unpairedOcc);
        addUnpairedMatches(best2, reads.getBundle2(), SECOND_IN_PAIR,
                           fw2.size() - 1, counters, unpairedOcc);
    }
}

void SearchStrategy::pairOccurrencesForBestMapping(
    vector<TextOcc>& upstreamOccs, vector<TextOcc>& downStreamOccs,
    vector<PairedTextOccs>& pairs, length_t maxFragSize, length_t minFragSize,
    length_t uMax, length_t dMax, Counters& counters,
    set<length_t>& uTrimmedIds, set<length_t>& dTrimmedIds, const string& uSeq,
    const string& dSeq) const {

    if (upstreamOccs.size() == 0 || downStreamOccs.size() == 0) {
        return;
    }

    // sort the downstream occurrences on their begin postion of the range
    sort(downStreamOccs.begin(), downStreamOccs.end(),
         [](const TextOcc& a, const TextOcc& b) {
             return a.getIndexBegin() < b.getIndexBegin();
         });

    // loop over all matches of the upstream read
    for (length_t i = 0; i < upstreamOccs.size(); i++) {
        auto& upStreamOcc = upstreamOccs[i];
        length_t upstreamPos = upStreamOcc.getIndexBegin();

        // find the index of first match of the downstream read that is
        // after the upstream read
        auto it = findFirstOccAfter(downStreamOccs, upstreamPos);

        // loop over all matches of downstream read that are after the
        // upstream read
        for (; it != downStreamOccs.end(); it++) {
            // check if the matches are within the fragment size
            // and have the correct orientation

            length_t fragSize = (it->getIndexEnd()) - upstreamPos;

            if (fragSize <= maxFragSize && fragSize >= minFragSize) {
                // found a pair
                // assign sequence if not assigned yet

                auto uFound =
                    assignSequenceAndCIGAR(upStreamOcc, counters, uMax, uSeq);
                if (uFound != FOUND) {
                    // no pair with this upStreamOcc will be valid so
                    // break the loop on downstream to continue to next
                    if (uFound == FOUND_WITH_TRIMMING) {
                        uTrimmedIds.insert(i);
                    }
                    break;
                }

                auto dFound = assignSequenceAndCIGAR(*it, counters, dMax, dSeq);

                if (dFound != FOUND) {
                    // no pair with this downStreamOcc will be valid so
                    // continue to next downstream occ
                    if (dFound == FOUND_WITH_TRIMMING) {
                        dTrimmedIds.insert(it - downStreamOccs.begin());
                    }
                    continue;
                }
                if (upStreamOcc.getAssignedSequenceID() !=
                    it->getAssignedSequenceID()) {
                    // the assigned sequences are not the same, so no
                    // pair can be formed (discordant pairs are handled
                    // later on)
                    continue;
                }

                pairs.push_back({upStreamOcc, *it});

            } else if (fragSize > maxFragSize) {
                // since the matches are sorted, we can stop the
                // loop
                break;
            }
        }
    }
}
// ----------------------------------------------------------------------------
// POST PROCESSING
// ----------------------------------------------------------------------------

void SearchStrategy::generateOutputSingleEnd(vector<TextOcc>& occs,
                                             ReadBundle& bundle,
                                             Counters& counters,
                                             length_t cutOff) const {

    if (occs.empty()) {
        if (unmappedSAM) {
            occs.emplace_back(createUnmappedRecordSE(bundle));
        }
        return;
    }

    vector<TextOcc> assignedOccs;
    assignedOccs.reserve(occs.size());

    // for each occurrence find the assigned sequence
    for (auto& t : occs) {
        index.setIndexInMode(t.getStrand());
        const auto& seq = bundle.getSequence(t.getStrand());
        length_t seqID;
        SeqNameFound found =
            index.findSeqName(t, seqID, counters, cutOff, distanceMetric, seq);
        if (found != SeqNameFound::NOT_FOUND) {
            t.setAssignedSequence(found, seqID);
            assignedOccs.emplace_back(std::move(t));
        }
    }

    // keep only the occurrences that have an assigned sequence
    occs = std::move(assignedOccs);

    if (occs.empty()) {
        if (unmappedSAM) {
            occs.emplace_back(createUnmappedRecordSE(bundle));
        }
        return;
    }
#ifdef DEVELOPER_MODE
    // sort the occurrences on score
    // Use stable sort because we expect many elements with same key
    stable_sort(occs.begin(), occs.end(),
                [](const TextOcc& a, const TextOcc& b) {
                    return a.getDistance() < b.getDistance();
                });

    // find the minimal score
    uint32_t minScore = occs.front().getDistance();

    uint32_t nHits = 1;

    // count the number of hits with the minimal score
    // take into account that the occurrences are sorted on score
    // find the index of the first occurrence with a score higher than
    // minScore with a binary search
    auto it = lower_bound(occs.begin(), occs.end(), minScore + 1);
    nHits = distance(occs.begin(), it);

#else
    // Find minimal distance
    auto minElement = std::min_element(
        occs.begin(), occs.end(), [](const TextOcc& a, const TextOcc& b) {
            return a.getDistance() < b.getDistance();
        });

    length_t minScore = minElement->getDistance();

    // Count number of elements with minimal distance
    length_t nHits =
        std::count_if(occs.begin(), occs.end(), [&](const TextOcc& elem) {
            return elem.getDistance() == minScore;
        });

    // Swap element with minimal distance to primary position
    if (minElement != occs.begin()) {
        std::iter_swap(occs.begin(), minElement);
    }
#endif
    (this->*generateOutputSEPtr)(bundle, nHits, minScore, occs, counters);
}

void SearchStrategy::generateSAMPairedEnd(vector<PairedTextOccs>& pairedMatches,
                                          ReadPair& rp) {

    if (pairedMatches.empty()) {
        return;
    }

// sort the paired matches on their
// cumulative distance score
// This gets expensive if the number of pairs is large/ We can solve
// this by using a vector of reference_wrappers or pointers, but
// this will need some refactoring, which might not be worth it
#ifdef DEVELOPER_MODE

    std::stable_sort(pairedMatches.begin(), pairedMatches.end(),
                     [](const PairedTextOccs& a, const PairedTextOccs& b) {
                         return a.getDistance() < b.getDistance();
                     });

    uint32_t bestScore = pairedMatches.front().getDistance();
    // find the number of pairs that have the best score
    // using a binary search
    auto it =
        lower_bound(pairedMatches.begin(), pairedMatches.end(), bestScore + 1);
    uint32_t nPairs = distance(pairedMatches.begin(), it);

#else

    // Find minimal distance
    auto minElement =
        std::min_element(pairedMatches.begin(), pairedMatches.end(),
                         [](const PairedTextOccs& a, const PairedTextOccs& b) {
                             return a.getDistance() < b.getDistance();
                         });

    uint32_t bestScore = minElement->getDistance();

    // Count number of elements with minimal distance
    length_t nPairs = std::count_if(pairedMatches.begin(), pairedMatches.end(),
                                    [&](const PairedTextOccs& elem) {
                                        return elem.getDistance() == bestScore;
                                    });

    // Swap element with minimal distance to primary position
    if (minElement != pairedMatches.begin()) {
        std::iter_swap(pairedMatches.begin(), minElement);
    }

#endif

    bool primaryAlignment = true;
    // generate the SAM lines for the paired matches
    for (auto& match : pairedMatches) {

        // upstream is first in pair
        bool upFirst = match.getUpStream().isFirstReadInPair();

        // generate SAM for upstream
        if (match.getUpStream().isValid()) {
            auto& bundle = upFirst ? rp.getBundle1() : rp.getBundle2();
            match.getUpStream().generateSAMPairedEnd(
                bundle, nPairs, bestScore, match.getDownStream(),
                match.getFragSize(), match.isDiscordant(), primaryAlignment,
                index.getSeqNames());
        }

        // generate SAM for downstream
        if (match.getDownStream().isValid()) {
            auto& bundle = upFirst ? rp.getBundle2() : rp.getBundle1();
            match.getDownStream().generateSAMPairedEnd(
                bundle, nPairs, bestScore, match.getUpStream(),
                match.getFragSize(), match.isDiscordant(), primaryAlignment,
                index.getSeqNames());
        }
        primaryAlignment = false;
    }
}

// ============================================================================
// CLASS CUSTOM SEARCH STRATEGY
// ============================================================================

// ----------------------------------------------------------------------------
// CONSTRUCTION
// ----------------------------------------------------------------------------

void CustomSearchStrategy::getSearchSchemeFromFolder(const string& pathToFolder,
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

    // get the info per distance score (scores between 1 and MAX_K are
    // looked at)
    for (int i = 1; i <= MAX_K; i++) {
        string name =
            pathToFolder + to_string(i) + PATH_SEPARATOR + "searches.txt";
        ifstream stream_searches(name);
        if (!stream_searches) {
            // this score is not supported
            supportsMaxScore[i - 1] = false;
            continue;
        }

        schemePerED[i - 1] =
            SearchScheme::readScheme(stream_searches, name, i).getSearches();

        if (schemePerED[i - 1].size() > 0) {
            supportsMaxScore[i - 1] = true;
        }
        stream_searches.close();
    }

    // check if the searches are valid
    sanityCheck(verbose);

    // get static positions (if they exist)
    for (int i = 1; i <= MAX_K; i++) {
        ifstream stream_static(pathToFolder + to_string(i) + PATH_SEPARATOR +
                               "static_partitioning.txt");

        if (stream_static) {
            // a file with static partitioning positions exists
            getline(stream_static, line);
            vector<string> positionsAsString = {};
            stringstream ss(line);
            string token;
            while (ss >> token) {
                positionsAsString.push_back(token);
            }

            if (positionsAsString.size() != calculateNumParts(i) - 1) {
                throw runtime_error(
                    "Not enough static positions provided in " + pathToFolder +
                    to_string(i) + PATH_SEPARATOR +
                    "static_partitioning.txt\nExpected: " +
                    to_string(calculateNumParts(i) - 1) + " parts\nProvided: " +
                    to_string(positionsAsString.size()) + " parts");
            }

            for (auto str : positionsAsString) {
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
    for (int i = 1; i <= MAX_K; i++) {
        ifstream stream_dynamic(pathToFolder + to_string(i) + PATH_SEPARATOR +
                                "dynamic_partitioning.txt");

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
                    to_string(i) + PATH_SEPARATOR +
                    "dynamic_partitioning.txt\nExpected: " +
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

/**
 * Helper function for CustomSearchStrategy::makeSearch.
 * Extracts the vector from a string and stores it in the vector.
 * @param vectorString The string representation of the vector.
 * @param vector The vector to store the values in.
 */
void getVector(const string& vectorString, vector<length_t>& vector) {

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

Search CustomSearchStrategy::makeSearch(const string& line,
                                        length_t idx) const {
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

    return Search::makeSearch(order, lower_bound, upper_bound, idx);
}

// ----------------------------------------------------------------------------
// SANITY CHECKS
// ----------------------------------------------------------------------------

void CustomSearchStrategy::sanityCheckStaticPartitioning(
    const int& maxScore) const {
    const auto& positions = staticPositions[maxScore - 1];

    // no length zero + increasing + all smaller than 1 and greater than
    // 0
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

    // no length zero + increasing + all smaller than 1 and greater than
    // 0
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

    // check if for each supported edit distance all error distributions
    // are covered

    for (int K = 1; K <= MAX_K; K++) {

        const auto& scheme = schemePerED[K - 1];
        if (!supportsMaxScore[K - 1]) {
            continue;
        }
        const Search& firstSearch = scheme.front();
        length_t P = firstSearch.getNumParts();
        // check if all searches have same number of parts
        if (any_of(scheme.begin(), scheme.end(),
                   [P](const Search& s) { return s.getNumParts() != P; })) {
            logger.logError("Not all searches for distance " + to_string(K) +
                            " have the same number of "
                            "parts");
            throw runtime_error("Not all searches for "
                                "distance " +
                                to_string(K) +
                                " have the same "
                                "number of parts");
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
                   [](const Search& s) { return !s.validBounds(); })) {
            throw runtime_error("Decreasing lower or upper bounds "
                                "for a search for K  = " +
                                to_string(K));
        }
    }
}

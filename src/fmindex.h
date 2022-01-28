/******************************************************************************
 *  Columba: Approximate Pattern Matching using Search Schemes                *
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
#ifndef FMINDEX_H
#define FMINDEX_H

#include <google/sparse_hash_map>

#include "alphabet.h"
#include "bandmatrix.h"
#include "bwtrepr.h"
#include "suffixArray.h"
#include "tkmer.h"

#include <algorithm> //used for sorting
#include <fstream>   // used for reading in files
#include <iostream>  // used for printing
#include <math.h>    //for taking the log
#include <numeric>   // for summing over vector
#include <sstream>   // used for splitting strings
#include <string>    // strings
#include <vector>    // vectors

// compute |a-b| in a safe manner
template <typename T> T abs_diff(T a, T b) {
    return a > b ? a - b : b - a;
}

// ============================================================================
// (TYPE) DEFINITIONS AND PROTOTYPES
// ============================================================================

typedef uint32_t length_t;

// ============================================================================
// CLASS RANGE
// ============================================================================

class Range {
  private:
    length_t begin; // beginning of the range
    length_t end;   // end of the range (non-inclusive)

  public:
    /**
     * Constructor
     * @param b, the beginning of the range
     * @param e, the end of the range (non-inclusive)
     */
    Range(length_t b, length_t e) : begin(b), end(e) {
    }

    /**
     * Default constructor, initializes an empty range
     */
    Range() : begin(0), end(0) {
    }

    length_t getBegin() const {
        return begin;
    }
    length_t getEnd() const {
        return end;
    }
    /**
     * Check if this range is empty
     * @returns true if the range is empty, false otherwise
     */
    bool empty() const {
        return end <= begin;
    }

    /**
     * Gets the width of the range (end - begin)
     * @returns the width of this range
     */
    length_t width() const {
        return (empty()) ? 0 : end - begin;
    }

    /**
     * Operator overloading, two ranges are equal if their begin and end field
     * are equal
     */
    bool operator==(const Range& o) const {
        return o.getBegin() == begin && o.getEnd() == end;
    }
    friend std::ostream& operator<<(std::ostream& os, const Range& r);
};

/**
 * Operator overloading. Outputs the range as [begin, end) of the output
 * stream
 * @param output, the output stream
 * @param r, the range to print
 */
std::ostream& operator<<(std::ostream& output, const Range& r);

// ============================================================================
// CLASS TextOccurrence
// ============================================================================
class TextOcc {
  private:
    Range range;       // the range in the text
    length_t distance; // the distance to this range (edit or hamming)
    std::vector<std::pair<char, uint>> CIGAR; // The CIGAR string of the match

    std::string output; // the corresponding output for this occurrence (for
                        // now a custom format)

  public:
    /**
     * Constructor
     * @param range, the range of this occurrence in the text
     * @param distance, the (edit or hamming) distance to the mapped read of
     * this occurrence
     */
    TextOcc(Range range, length_t distance)
        : range(range), distance(distance), CIGAR(), output() {
    }

    /**
     * Constructor
     * @param range, the range of this occurrence in the text
     * @param distance, the (edit or hamming) distance to the mapped read of
     * this occurrence
     * @param CIGAR the CIGAR string of the match
     */
    TextOcc(Range range, length_t distance,
            std::vector<std::pair<char, uint>>& CIGAR)
        : range(range), distance(distance), CIGAR(CIGAR), output() {
    }

    /**
     * Constructor for an invalid text occurrence
     */
    TextOcc() : range(0, 0) {
    }

    /**
     * Generates the output of this occurrence, for now in format:
     * startposition\twidth\tdistance\tCIGAR, where startposition is the
     * beginning of thetext occurrence, width is the length of this occurrence,
     * distance is the (edit or hamming) distance to the mapped read and CIGAR
     * is the CIGAR string of the match
     */
    void generateOutput() {
        output = std::to_string(range.getBegin()) + "\t" +
                 std::to_string(range.width()) + "\t" +
                 std::to_string(distance) + "\t";
        for (const auto& p : CIGAR) {
            output += std::to_string(p.second) + p.first;
        }
    }

    const Range getRange() const {
        return range;
    }
    const length_t getDistance() const {
        return distance;
    }

    const std::string& getOutput() const {
        return output;
    }

    /**
     * Operator overloading for sorting the occurrences.
     * Occurrences are first sorted on their begin position, then on their
     * distance , their width  and finally on the existence of the CIGAR string
     */
    bool operator<(const TextOcc& r) {

        if (range.getBegin() != r.getRange().getBegin()) {
            return range.getBegin() < r.getRange().getBegin();
        } else {
            // begin is equal, better ed is smarter
            if (distance != r.getDistance()) {
                return distance < r.getDistance();
            } else if (range.width() != r.getRange().width()) {
                // shorter read is smaller...
                return range.width() < r.getRange().width();
            } else {
                return hasCigar() && !r.hasCigar();
            }
        }
    }

    bool operator==(const TextOcc& r) {
        return r.getRange() == range && r.getDistance() == distance;
    }

    bool isValid() const {
        return !range.empty();
    }

    length_t width() const {
        return range.width();
    }

    bool hasCigar() const {
        return !CIGAR.empty();
    }

    void setCigar(std::vector<std::pair<char, uint>>& cigar) {
        CIGAR = cigar;
    }
};

// ============================================================================
// CLASS SARANGEPAIR
// ============================================================================

/**
 * A pair of ranges. The first range is range over the suffix array of the
 * text. The second range is the corresponding range over the suffix array
 * of the reversed text
 */
class SARangePair {
  private:
    Range rangeSA;    // the range over the suffix array
    Range rangeSARev; // the range over the suffix array of the reversed text

  public:
    /**
     * Default constructor, creates two empty ranges
     */
    SARangePair() : rangeSA(Range()), rangeSARev(Range()) {
    }
    SARangePair(Range rangeSA, Range rangeSARev)
        : rangeSA(rangeSA), rangeSARev(rangeSARev) {
    }

    const Range& getRangeSA() const {
        return rangeSA;
    }

    const Range& getRangeSARev() const {
        return rangeSARev;
    }
    /**
     * @returns true if the ranges are empty, false otherwise
     */
    bool empty() const {
        return rangeSA.empty();
    }

    length_t width() const {
        return rangeSA.width();
    }
    /**
     * Operator overloading
     * @returns true if this is equal to rhs
     */
    bool operator==(const SARangePair& o) const {
        // only the first range matters as the ranges imply each other
        return o.getRangeSA() == rangeSA;
    }
};
// ============================================================================
// CLASS FMPos
// ============================================================================
/**
 * A position in the bidirectional FM-index.
 */
class FMPos {
  protected:
    SARangePair ranges; // the ranges over the suffix arrays
    length_t depth; // the depth of the prefix of the suffixes of this position

  public:
    /**
     * Default constructor for empty position (= empty ranges and depth of
     * zero)
     */
    FMPos() : ranges(SARangePair()), depth(0) {
    }

    FMPos(SARangePair& ranges, length_t depth) : ranges(ranges), depth(depth) {
    }

    const SARangePair& getRanges() const {
        return ranges;
    }

    const length_t& getDepth() const {
        return depth;
    }

    void setRanges(SARangePair ranges) {
        this->ranges = ranges;
    }

    void setDepth(length_t depth) {
        this->depth = depth;
    }

    /**
     * Operator overloading, two FMPos are equal if their ranges and depth
     * are equal
     * @param rhs the FMPos to compare to this
     * @returns true if this is equal to rhs
     */
    bool operator==(const FMPos& rhs) const {
        return ranges == rhs.getRanges() && depth == rhs.getDepth();
    }
    /**
     * @returns true if the ranges are not empty, false otherwise
     */
    bool isValid() const {
        return !ranges.empty();
    }
};
// ============================================================================
// CLASS FMOcc
// ============================================================================

/**
 * An occurrence in the bidirectional FM-index
 */
class FMOcc {
  private:
    FMPos pos;         // The FM position of this occurrence
    length_t distance; // the distance (hamming or edit)
    length_t shift; // A right-sift to the corresponding positions in the text

  public:
    FMOcc() : pos(), distance(0), shift(0) {
    }
    /**
     * Make a bidirectional approximate match in the suffix array
     * @param ranges the ranges of this approximate match (range in SA and
     * in SA')
     * @param distance the (edit or hamming) distance of this approximate
     * match
     * @param depth the depth (=length) of this approximate match
     * @param shift The right shift to the corresponding positions in the
     * text, defaults to zero
     */
    FMOcc(SARangePair ranges, length_t distance, length_t depth,
          length_t shift = 0)
        : pos(ranges, depth), distance(distance), shift(shift) {
    }
    /**
     * Make a bidirectional approximate match in the suffix array
     * @param pos, the position in the FMIndex of this approximate match
     * @param distance the (edit or hamming) distance of this approximate
     * match
     * @param shift The right shift to the corresponding positions in the
     * text, defaults to zero
     */
    FMOcc(FMPos pos, length_t distance, length_t shift = 0)
        : pos(pos), distance(distance), shift(shift) {
    }

    const SARangePair& getRanges() const {
        return pos.getRanges();
    }
    const length_t& getDistance() const {
        return distance;
    }

    const length_t& getDepth() const {
        return pos.getDepth();
    }

    length_t getWidth() const {
        return pos.getRanges().width();
    }

    const length_t& getShift() const {
        return shift;
    }

    void setRanges(SARangePair ranges) {
        pos.setRanges(ranges);
    }

    void setDistance(length_t distance) {
        this->distance = distance;
    }

    void setDepth(length_t depth) {
        pos.setDepth(depth);
    }
    /**
     * @returns true if the position is valid, false otherwise
     */
    bool isValid() const {
        return pos.isValid();
    }
    /**
     * Operator overloading to sort FMOcc
     * First the FMOcc are sorted on the begin of the range over the suffix
     * array of their position Then they are sorted on their distance
     * Lastly they are sorted on the width of their ranges
     * @param rhs the FMOcc to compare to this
     * @returns true if this is smaller than rhs
     */
    bool operator<(const FMOcc& rhs) const {
        if (pos.getRanges().getRangeSA().getBegin() !=
            rhs.getRanges().getRangeSA().getBegin()) {
            return getRanges().getRangeSA().getBegin() <
                   rhs.getRanges().getRangeSA().getBegin();
        } else if (distance != rhs.getDistance()) {
            // begin is equal, better ed is smarter
            return distance < rhs.getDistance();
        } else if (getRanges().width() != rhs.getRanges().width()) {
            // shorter read is smaller...
            return getRanges().width() < rhs.getRanges().width();
        } else {
            // prefer no shift
            return getShift() < rhs.getShift();
        }
    }
    /**
     * Operator overloading
     * Two FMocc are equal if their ranges, distance and depth are all equal
     * @param returns true if this is equal to rhs
     */
    bool operator==(const FMOcc& rhs) {
        return getRanges() == rhs.getRanges() &&
               distance == rhs.getDistance() && getDepth() == rhs.getDepth() &&
               getShift() == rhs.getShift();
    }
    friend std::ostream& operator<<(std::ostream& os, const FMOcc& biocc);
};

// ============================================================================
// CLASS FMPosExt
// ============================================================================
/**
 * A single node in the bidirectional FM-index. Its depth is the depth from
 * the startmatch for a particular phase of a search
 */
class FMPosExt : public FMPos {
  private:
    char c;                // the character of this node
    bool reported = false; // has this particular node already reported?
  public:
    /**
     * Create a node of the search tree
     * @param character the character of this node
     * @param ranges the ranges over the suffix and reversed suffix array
     * that go to this node
     * @param row the row of this node in the alignment matrix = depth of
     * this node
     */
    FMPosExt(char character, SARangePair ranges, length_t row)
        : FMPos(ranges, row), c(character), reported(false) {
    }

    /**
     * Default constructor, this Node will have empty ranges
     */
    FMPosExt() : FMPos(), c(char(0)) {
    }

    /**
     * Sets the report flag to true
     */
    void report() {
        reported = true;
    }

    /**
     * Reports the match (with added depth) at this node,
     * @param occ the match will be stored here
     * @param startDepth the depth to add to the match
     * @param EDFound the found edit distance for this node
     * @param noDoubleReports false if this node is allowed to report more
     * than once, defaults to false
     * @param shift, right shift of the matche, defaults to zero
     */
    void report(FMOcc& occ, const length_t& startDepth, const length_t& EDFound,
                const bool& noDoubleReports = false, length_t shift = 0) {
        if (!reported) {
            occ = FMOcc(getRanges(), EDFound, depth + startDepth, shift);

            // if final part, report only once
            if (noDoubleReports) {
                report();
            }
        }
    }

    /**
     * Gets the ranges of this node
     * @returns the ranges of this node
     */
    const SARangePair& getRanges() const {
        return ranges;
    }

    /**
     * Get the character of this node
     * @returns the character of this node
     */
    const char getCharacter() const {
        return c;
    }

    /**
     * Get the row of this node
     * @returns the row of this node
     */
    length_t getRow() const {
        return depth;
    }
};

// ============================================================================
// CLASS CLUSTER
// ============================================================================

class Cluster {
  private:
    std::vector<uint> eds;       // the edit distances of this cluster
    std::vector<FMPosExt> nodes; // the nodes of this cluster

    length_t lastCell;   // the lastCell of the cluster that was filled in
    uint maxED;          // the maxEd for this cluster
    length_t startDepth; // the startdepth for this cluster (= depth of match
                         // before matrix of this cluster)

    length_t shift; // the right shift of the occurrences in the text
  public:
    /**
     * Constructor
     * @param size, the size of the cluster
     * @param maxED, the maximal allowed edit distance
     * @param startDepth, the depth before this cluster
     * @param shift, the right shift of the occurrences in the text
     */
    Cluster(length_t size, length_t maxED, length_t startDepth, length_t shift)
        : eds(size, maxED + 1), nodes(size), lastCell(-1), maxED(maxED),
          startDepth(startDepth), shift(shift) {
    }

    /**
     * Sets the ed and node at index idx to ed and node. Also updates
     * lastCell to be idx
     * @param idx, the idx to change
     * @param node, the node to set at index idx
     * @param ed, the ed to set at index idx
     */
    void setValue(length_t idx, const FMPosExt& node, const length_t& ed) {
        eds[idx] = ed;
        nodes[idx] = node;
        lastCell = idx;
    }

    /**
     * Returns the size of this cluster
     */
    const length_t size() const {
        return eds.size();
    }

    /**
     * @returns vector with all nodes in the cluster that are a centre and
     * under the maximal allowed distance, be aware that if there are
     * multiple centers in the cluster it is very likely that only one of
     * them will be redundant, but the others might eliminate another
     * occurrence
     */
    std::vector<FMOcc> reportCentersAtEnd() {

        std::vector<FMOcc> centers;
        centers.reserve(lastCell + 1);

        for (length_t i = 0; i <= lastCell; i++) {
            if (eds[i] <= maxED && (i == 0 || eds[i] <= eds[i - 1]) &&
                (i == lastCell || eds[i] <= eds[i + 1])) {
                FMOcc m;
                nodes[i].report(m, startDepth, eds[i], true, shift);
                centers.emplace_back(m);
            }
        }

        return centers;
    }

    /**
     * @returns approximate match that corresponds to the ranges of the
     * deepest global minimum of this cluster, but with the depth of the
     * highest global minimum of this cluster. If the direction is backward
     * a shift will be set such that the occurrence in the text will be as
     * short as possible
     */
    FMOcc reportDeepestMinimum(Direction dir) {
        uint minED = maxED + 1;
        length_t highestBestIdx = -1;
        length_t deepestBestIdx = -1;

        for (length_t i = 0; i <= lastCell; i++) {
            if (eds[i] < minED) {
                minED = eds[i];
                highestBestIdx = i;
                deepestBestIdx = i;
            }
            if (eds[i] == minED) {
                deepestBestIdx = i;
            }
        }
        FMOcc m;
        if (minED <= maxED) {
            nodes[deepestBestIdx].report(
                m, startDepth - (deepestBestIdx - highestBestIdx), minED, true,
                ((dir == BACKWARD) ? (deepestBestIdx - highestBestIdx) : 0) +
                    shift);
        }
        return m;
    }

    /**
     * This method returns a match that corresponds to the highest cluster
     * centre. Its descendants and the corresponding initializationeds are
     * updated. Eds of descendants that are part of a cluster centre which
     * is lower than the lower bound will be updated in the initEds vector
     * @param lowerBound, the lower bound for this iteration
     * @param desc, the descendants of the highest cluster centre, these
     * will be inserted during the method
     * @param initEds, the initialization eds for the next iteration, these
     * correspond to the eds of the highest centre and its descendants,
     * where eds part of a cluster of which the centre is below the
     * lower bound are updated. These values will be inserted during the
     * method
     * @returns The occurrence corresponding to the upper cluster centre
     * which has a valid distance
     */
    FMOcc getClusterCentra(uint lowerBound, std::vector<FMPosExt>& desc,
                           std::vector<uint>& initEds) {
        desc.reserve(eds.size());
        initEds.reserve(eds.size());
        FMOcc m;
        for (length_t i = 0; i <= lastCell; i++) {
            if (eds[i] > maxED || eds[i] < lowerBound) {
                continue;
            }
            bool betterThanParent = (i == 0) || eds[i] <= eds[i - 1];
            bool betterThanChild = (i == lastCell) || eds[i] <= eds[i + 1];

            if (betterThanParent && betterThanChild) {
                // this is a valid centre
                nodes[i].report(m, startDepth, eds[i], false, shift);

                // get all the descendants
                initEds.emplace_back(eds[i]);
                for (length_t j = i + 1; j <= lastCell; j++) {
                    desc.emplace_back(nodes[j]);
                    initEds.emplace_back(eds[j]);
                }

                // replace the clusters under the lower bound
                for (length_t k = 1; k < initEds.size(); k++) {
                    if (initEds[k] < lowerBound &&
                        initEds[k] <= initEds[k - 1] &&
                        (k == initEds.size() - 1 ||
                         initEds[k] <= initEds[k + 1])) {
                        // k is a centre under the lower bound

                        length_t highestPoint = 0;
                        length_t lowestPoint = initEds.size() - 1;
                        // find highest point of this cluster
                        for (length_t l = k; l-- > 0;) {
                            if (initEds[l] != initEds[l + 1] + 1) {
                                highestPoint = l + 1;
                                break;
                            }
                        }
                        // find lowest point of this cluster
                        for (length_t l = k + 1; l < initEds.size(); l++) {
                            if (initEds[l] != initEds[l - 1] + 1) {
                                lowestPoint = l - 1;
                                break;
                            }
                        }

                        // highest and lowest cannot span entire
                        // initEds.size(), otherwise there would not be a
                        // valid cluster centre above the lower bound
                        if (highestPoint != 0 &&
                            lowestPoint != initEds.size() - 1) {
                            // Make /\ with ed values of this cluster
                            // do iE[hp] = ie[hp - 1] + 1 and iE[lp] = iE[lp
                            // + 1] +1 until entire cluster has been
                            // replaced
                            length_t lC = lowestPoint;
                            length_t hC = highestPoint;
                            bool highest = true;
                            // do not go over maxED + 1, to ensure
                            // continuity at the other end
                            while (lC > hC) {
                                if (highest) {
                                    initEds[hC] = std::min(maxED + 1,
                                                           initEds[hC - 1] + 1);
                                    hC++;
                                } else {
                                    initEds[lC] = std::min(maxED + 1,
                                                           initEds[lC + 1] + 1);
                                    lC--;
                                }
                                highest = !highest;
                            }
                            if (lC == hC) {
                                // change middle element of cluster
                                initEds[lC] = std::min(initEds[lC + 1] + 1,
                                                       initEds[lC - 1] + 1);
                            }

                        } else if (highestPoint == 0 &&
                                   lowestPoint != initEds.size() - 1) {
                            // monotonous rise from lowestPoint to
                            // highestPoint
                            for (length_t l = lowestPoint; l-- > 0;) {
                                initEds[l] = initEds[l + 1] + 1;
                            }
                        } else if (highestPoint != 0 &&
                                   lowestPoint == initEds.size() - 1) {
                            // monotonous rise from highestPoint to
                            // lowestPoint
                            for (length_t l = highestPoint; l < initEds.size();
                                 l++) {
                                initEds[l] = initEds[l - 1] + 1;
                            }
                        }
                    }
                }
                // stop searching
                break;
            }
        }

        return m;
    }
};

// ============================================================================
// CLASS SEARCH
// ============================================================================
class Search {
  private:
    std::vector<length_t> lowerBounds; // the vector with lower bounds
    std::vector<length_t> upperBounds; // the vector with upper bounds
    std::vector<length_t> order;       // the vector with the order of the parts
    std::vector<Direction> directions; // the directions of each phase
    std::vector<bool>
        directionSwitch; // has the direction switched for each phase

    std::vector<std::pair<length_t, length_t>>
        lowestAndHighestPartsProcessedBefore;

    Search(std::vector<length_t>& order, std::vector<length_t>& lowerBounds,
           std::vector<length_t>& upperBounds,
           std::vector<Direction>& directions, std::vector<bool>& dSwitch,
           std::vector<std::pair<length_t, length_t>>
               lowestAndHighestPartsProcessedBefore)
        : lowerBounds(lowerBounds), upperBounds(upperBounds), order(order),
          directions(directions), directionSwitch(dSwitch),
          lowestAndHighestPartsProcessedBefore(
              lowestAndHighestPartsProcessedBefore) {
    }

  public:
    /**
     * Static function to construct a search. The directions and switches of
     * the search are calculated
     * @param order, the order of the search
     * @param lowerBounds, the lower bounds of the search
     * @param upperBounds, the upper bounds of the search
     */
    static Search makeSearch(std::vector<length_t> order,
                             std::vector<length_t> lowerBounds,
                             std::vector<length_t> upperBounds) {
        // check correctness of sizes
        if (order.size() != lowerBounds.size() ||
            order.size() != upperBounds.size()) {
            throw std::runtime_error("Could not create search, the sizes of "
                                     "all vectors are not equal");
        }

        // compute the directions
        std::vector<Direction> directions;
        directions.reserve(order.size());
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
        // TODO check in case of Kianfar
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

        return Search(order, lowerBounds, upperBounds, directions,
                      directionSwitch, lowestAndHighestPartsProcessedBefore);
    }

    /**
     * Sets the directions of the parts to the directions of the search
     * @param parts, the parts to set the direction of
     */
    void setDirectionsInParts(std::vector<Substring>& parts) const {
        // set the directions for the parts
        for (length_t i = 0; i < order.size(); i++) {
            parts[order[i]].setDirection(directions[i]);
        }
    }

    /**
     * @returns the lower bound for the idx'th part
     */
    length_t getLowerBound(length_t idx) const {
        assert(idx < lowerBounds.size());
        return lowerBounds[idx];
    }

    /**
     * @returns the upper bound for the idx'th part
     */
    length_t getUpperBound(length_t idx) const {
        assert(idx < upperBounds.size());
        return upperBounds[idx];
    }
    /**
     * @returns  the idx'th part
     */
    length_t getPart(length_t idx) const {
        assert(idx < order.size());
        return order[idx];
    }

    length_t getLowestPartProcessedBefore(length_t idx) const {
        assert(idx >= 1);
        assert(idx < order.size());
        return lowestAndHighestPartsProcessedBefore[idx - 1].first;
    }
    length_t getHighestPartProcessedBefore(length_t idx) const {
        assert(idx >= 1);
        assert(idx < order.size());
        return lowestAndHighestPartsProcessedBefore[idx - 1].second;
    }

    length_t getMaxED() const {
        return upperBounds.back();
    }

    length_t getMinED() const {
        return lowerBounds.back();
    }

    /**
     * @returns the direction for the idxith part
     */
    Direction getDirection(length_t idx) const {
        assert(idx < directions.size());
        return directions[idx];
    }
    /**
     * @returns if the direction switches at the idxth part
     */
    bool getDirectionSwitch(length_t idx) const {
        assert(idx < directionSwitch.size());
        return directionSwitch[idx];
    }

    /**
     * Get the number of parts in this search
     * @return the number of parts
     */
    length_t getNumParts() const {
        return order.size();
    }

    /**
     * Checks if the idxth part is the first or last part of pattern
     * @returns true if the idxth part is the first or last part, false
     * otherwise
     */
    bool isEdge(length_t idx) const {
        return order[idx] == 0 || order[idx] == getNumParts() - 1;
    }

    /**
     * Checks if the idxth part is the final part of the search
     * @returns true if idx is the final part of the search, false otherwise
     */
    bool isEnd(length_t idx) const {
        return idx == order.size() - 1;
    }

    /**
     * Checks if the connectivity property is satisfied
     * @returns true if the property is satisfied
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
     * Checks if the upper and lower bounds are not decreasing
     * @returns true if the bounds are valid
     */
    bool noDecreasingInBounds() const {
        for (length_t i = 1; i < order.size(); i++) {
            if (lowerBounds[i] < lowerBounds[i - 1]) {
                return false;
            }
            if (upperBounds[i] < upperBounds[i - 1]) {
                return false;
            }
        }
        return true;
    }

    /**
     * Check if the search is zero based (order must contain zero)
     * @returns  true if the search is zero based else false
     */
    bool zeroBased() const {
        return *std::min_element(order.begin(), order.end()) == 0;
    }
};
/**
 * Operator overloading. Outputs the search to theoutput stream
 * @param os, the output stream
 * @param obj, the search to print
 */
std::ostream& operator<<(std::ostream& os, const Search& obj);

// ============================================================================
// STRUCT COUNTERS
// ============================================================================
// A struct of performance counters
struct Counters {
    // performance counters
    thread_local static length_t
        nodeCounter; // counts the number of nodes visited in the index

    thread_local static length_t
        totalReportedPositions; // counts the number of matches reported (either
                                // via in-text verification or in-index
                                // matching)

    thread_local static length_t
        cigarsInIndex; // counts the number of cigar strings calculated for
                       // matches in the index, note that this is only
                       // calculated if the match is non-redundant

    thread_local static length_t
        inTextStarted; // counts the number of times in-text verification was
                       // started, this equals the number of look-ups in the
                       // suffix array for in-text verification
    thread_local static length_t
        abortedInTextVerificationCounter; // counts the number of unsuccesful
                                          // in-text verifications

    thread_local static length_t
        cigarsInTextVerification; // counts the number of cigars strings
                                  // calculated for matches in the text, note
                                  // that this is done for each match in the
                                  // text as at the point of calculation it is
                                  // not known if this match will turn out to be
                                  // redundant
    thread_local static length_t
        usefulCigarsInText; // counts the number of cigar strings calculated for
                            // non-redundant matches in the text

    thread_local static length_t
        immediateSwitch; // Counts the number of times the partial matches after
                         // the first part has been matched are immediately
                         // in-text verified
    thread_local static length_t
        approximateSearchStarted; // Counts the number of times a search does
                                  // start

    /**
     * Reset all counters to 0
     */
    void resetCounters() {
        nodeCounter = 0, abortedInTextVerificationCounter = 0,
        totalReportedPositions = 0, cigarsInIndex = 0,
        cigarsInTextVerification = 0, inTextStarted = 0, usefulCigarsInText = 0,
        immediateSwitch = 0, approximateSearchStarted = 0;
    }
};

// ============================================================================
// CLASS FMIndex
// ============================================================================
class Occurrences;
class FMIndex;
typedef bool (FMIndex::*ExtraCharPtr)(length_t, const SARangePair&,
                                      SARangePair&) const;

typedef void (FMIndex::*FindDiffPtr)(length_t&, length_t&, length_t&,
                                     const length_t&, const length_t&,
                                     const BitParallelED&, const length_t&,
                                     const length_t&, const length_t&,
                                     const std::vector<uint>&) const;

class FMIndex {
  private:
    // info about the text
    const std::string baseFile; //  The basefile of the reference text
    length_t textLength;        // the length of the text
    std::string text;

    Alphabet<ALPHABET> sigma; // the alphabet

    // info about the suffix array
    length_t sparseFactorSA =
        32; // the sparseness factor of the suffix array, defaults to 32
    int logSparseFactorSA = 5; // the log of the sparse factor

    // bidirectional fm index data structures
    std::string bwt;              // the bwt string of the reference genome
    std::vector<length_t> counts; // the counts array of the reference genome
    SparseSuffixArray sparseSA;   // the suffix array of the reference genome
    BWTRepr<ALPHABET> fwdRepr;    // the baseFile occurrences table
    BWTRepr<ALPHABET> revRepr;    // the baseFile occurrences of the rev BWT

    // in-text verification
    length_t inTextSwitchPoint = 5;

    // direction variables
    thread_local static Direction dir; // the direction of the index
    thread_local static ExtraCharPtr
        extraChar; // pointer to extra char method (for direction)

    thread_local static FindDiffPtr findDiff; // pointer to correct finddiffmode

    // stacks for search schemes
    thread_local static std::vector<std::vector<FMPosExt>>
        stacks; // stacks of nodes for the different partitions
    thread_local static std::vector<BitParallelED>
        matrices; // alignment matrices for the different partitions

    // sparse hash info
    static const size_t wordSize =
        4; // the size the mers to be stored in a table
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
     * @param symbolIndex the index of the symbol in the alphabet to
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
     * Same as FMIndex::getNumberOfOccurrences, but now in the bwt of the
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
     * @param intextMatrix the bit-parallel matrix for in-text verification
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
        BitParallelED& intextMatrix, const Search& search,
        const FMOcc& startMatch, Occurrences& occ,
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
     * Goes deeper in a search if a valid approximate match is found in the
     * cluster
     * @param intextMatrix the matrix used for in-text verifications of the
     * current pattern
     * @param cluster, the cluster to search for a valid approximate match
     * @param nextIdx, the idx of next part to research
     * @param s, the search
     * @param parts, the parts of the pattern
     * @param occ Datastructure with the in-index and in-text occurrences, if an
     * occurrence is found it will be added to this datastructure
     * @param lowerbound, the lower bound for this partition
     * @param counters the performance counters
     * @param descendantsOtherD, the descendants of the other direction,
     * defaults to empty vector
     * @param ininEdsOtherD, the initialization eds of the other direction,
     * defaults to empty vector
     * @param remainingDesc, the remaining descendants on the current
     * branch, that are already created but aren't checked yet and need to
     * be checked for the next part, defaults to an empty vector
     */
    void goDeeper(BitParallelED& intextMatrix, Cluster& cluster,
                  const length_t& nextIdx, const Search& s,
                  const std::vector<Substring>& parts, Occurrences& occ,
                  const length_t& lowerBound, Counters& counters,
                  const std::vector<FMPosExt>& descendantsOtherD = {},
                  const std::vector<uint>& initEdsOtherD = {},
                  const std::vector<FMPosExt>& remainingDesc = {});

    /**
     * Helper function for the approximate matching. This function fills in
     * the matrix for the current node at the current row and goes deeper
     * for the next part is that is necessary The function true if the search
     * should backtrack.
     * @param intextMatrix the matrix used for in-text verifications of the
     * current pattern
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
    bool branchAndBound(BitParallelED& intextMatrix, Cluster& clus,
                        const FMPosExt& currentNode, const Search& s,
                        const length_t& idx,
                        const std::vector<Substring>& parts, Occurrences& occ,
                        Counters& counters, const std::vector<uint>& initOther,
                        const std::vector<FMPosExt>& descOther,
                        const std::vector<FMPosExt>& remainingDesc = {});

    // ----------------------------------------------------------------------------
    // IN TEXT VERIFICATION ROUTINES
    // ----------------------------------------------------------------------------

    /**
     * Verifies the text occurrences in the text and adds them to the
     * occurrences for the edit distance metric.
     * Note: hStartDec and hStartInc cannot both be higher than 0
     * @param to the partial text occurrences to be checked, these are converted
     * from FM occurrences
     * @param maxED the maximal edit distance that is allowed for the search
     * @param minED the minimal edit distance that is allowed for the search
     * @param intextMatrix the matrix for in text verification, the match
     * vectors should be set according to the pattern
     * @param occ the occurrences, both those found in the FM Index and in the
     * text, if a valid text occurrence is found it will be added
     * @param counters the performance counters
     * @param lStartDec the decrease,as compared to the start of the partial
     * match, which leads to the lowest possible start position
     * @param hStartDec the decrease, as compared to the start of the partial
     * match, which leads to the highest possible start position
     * @param hStartInc the increase, as compared to the start of the
     * partial match, which leads to the highest possible start position
     */
    void inTextVerification(const std::vector<TextOcc>& to,
                            const length_t& maxED, const length_t& minED,
                            BitParallelED& intextMatrix, Occurrences& occ,
                            Counters& counters, const length_t& lStartDec,
                            const length_t& hStartDec,
                            const length_t& hStartInc) const;
    /**
     * Helper function, needed for in-text verification for the edit distance.
     * Finds the differences with the partial start, which lead to the lowest
     * and highest possible start positions
     * Corresponds to case 1: the search is backwards
     * @param lStartDec the decrease,as compared to the start of the partial
     * match, which leads to the lowest possible start position (=output)
     * @param hStartDec the decrease, as compared to the start of the partial
     * match, which leads to the highest possible start position (=output)
     * @param hStartInc the increase, as compared to the start of the
     * partial match, which leads to the highest possible start position
     * (=output)
     * @param startBeforeThis the startposition in the pattern of the part with
     * lowest index that has already been processed
     * @param maxED the maximum allowed edit distance in the entire search
     * @param bpED the bit parall edit distance that was used for the current
     * part
     * @param maxEDPart the maximum allowed edit distance for the current part
     * @param descOtherSize the number of descendants in the other direction
     * (irrelevant for backwards)
     * @param initOther the final column in the other direction (irrelevant for
     * backwards)
     */
    void findDiffStartPositionBackward(
        length_t& lStartDec, length_t& hStartDec, length_t& hStartInc,
        const length_t& startBeforeThis, const length_t& maxED,
        const BitParallelED& bpED, const length_t& row,
        const length_t& maxEDPart, const length_t& descOtherSize,
        const std::vector<uint>& initOther) const;

    /**
     * Helper function, needed for in-text verification for the edit distance.
     * Finds the differences with the partial start, which lead to the lowest
     * and highest possible start positions
     * Corresponds to case 2, 3, 4: the search is currently forwards
     * @param lStartDec the decrease,as compared to the start of the partial
     * match, which leads to the lowest possible start position (=output)
     * @param hStartDec the decrease, as compared to the start of the partial
     * match, which leads to the highest possible start position (=output)
     * @param hStartInc the increase, as compared to the start of the
     * partial match, which leads to the highest possible start position
     * (=output)
     * @param startBeforeThis the startposition in the pattern of the part with
     * lowest index that has already been processed
     * @param maxED the maximum allowed edit distance in the entire search
     * @param bpED the bit parall edit distance that was used for the current
     * part
     * @param maxEDPart the maximum allowed edit distance for the current part
     * @param descOtherSize the number of descendants in the other direction
     * @param initOther the final column in the other direction
     */
    void findDiffStartPositionForward(
        length_t& lStartDec, length_t& hStartDec, length_t& hStartInc,
        const length_t& startBeforeThis, const length_t& maxED,
        const BitParallelED& bpED, const length_t& row,
        const length_t& maxEDPart, const length_t& descOtherSize,
        const std::vector<uint>& initOther) const;

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

    /**
     * Converts a match in the suffix array to matches in the text.
     * @param matchInSA the match that will be converted
     * @returns a vector with the corresponding text occurrences
     */
    std::vector<TextOcc> convertToMatchesInText(const FMOcc& matchInSA) const;

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
            int sa_sparse = 1, bool verbose = true)
        : baseFile(baseFile), sparseFactorSA(sa_sparse),
          logSparseFactorSA(log2(sa_sparse)), sparseSA(baseFile, sa_sparse),
          inTextSwitchPoint(inTextSwitch) {
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
     * Private helper function for exactMatches. This finds the
     * range in the suffix array that corresponds to matches in the
     * reference genome
     * @param s the string to be matched in the reference genome
     * @param counters the performance counters
     * @returns a Range containing the start and end values of the range
     * ([start, end[)
     */
    Range matchString(const std::string& s, Counters& counters) const;
    /**
     * Calculates the positions in the reference genome where exact matches
     * to the argument string start.
     * @param s the string to match in the reference genome
     * @param counters the performance counters
     * @returns a sorted vector containing the start positions of all exact
     * substring matches of s in the reference sequence
     */
    std::vector<length_t> exactMatches(const std::string& s,
                                       Counters& counters) const;

    /**
     * Calculates the exact matches to the string in the index and returns them
     * with their CIGAR string
     * @param s the string to match in the reference genome
     * @param counters the performance counters
     * @returns a sorted vector containing the start positions of all exact
     * substring matches of s in the reference sequence
     */
    std::vector<TextOcc> exactMatchesOutput(const std::string& s,
                                            Counters& counters) {
        setDirection(BACKWARD);
        const auto& positions = exactMatches(s, counters);

        length_t length = s.size();
        // create CIGAR string of all matches
        std::vector<std::pair<char, uint>> CIGAR(1, {'M', length});
        std::vector<TextOcc> textOccurrences;
        // Create text occurrences and set CIGAR string
        for (const auto p : positions) {
            textOccurrences.emplace_back(Range(p, p + length), 0, CIGAR);
            textOccurrences.back().generateOutput();
        }
        return textOccurrences;
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
    std::vector<TextOcc> approxMatchesNaive(const std::string& pattern,
                                            length_t maxED, Counters& counters);

    /**
     * Sets the search direction of the fm-index
     * @param d the direction to search in, either FORWARD or BACKWARD
     */
    void setDirection(Direction d) {
        dir = d;
        extraChar = (d == FORWARD) ? &FMIndex::findRangesWithExtraCharForward
                                   : &FMIndex::findRangesWithExtraCharBackward;
        findDiff = (d == FORWARD) ? &FMIndex::findDiffStartPositionForward
                                  : &FMIndex::findDiffStartPositionBackward;
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
     * @param intextMatrix the matrix for in-text verification, setSequence MUST
     * have been called before this function
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
    void recApproxMatchEditOptimizedEntry(BitParallelED& intextMatrix,
                                          const Search& search,
                                          const FMOcc& startMatch,
                                          Occurrences& occ,
                                          const std::vector<Substring>& parts,
                                          Counters& counters,
                                          const int& idx = 1) {

        if (startMatch.getRanges().width() > inTextSwitchPoint) {
            counters.approximateSearchStarted++;
            recApproxMatchEditOptimized(intextMatrix, search, startMatch, occ,
                                        parts, counters, idx);
            return;
        }
        // verify the partial match in text
        verifyExactPartialMatchInText(
            intextMatrix, startMatch,
            parts[search.getLowestPartProcessedBefore(idx)].begin(),
            search.getMaxED(), occ, counters);
    }

    /**
     * Verifies an exact partial match in the text for all occurrences of that
     * exact partial match
     * @param intextMatrix the matrix for in-text verification, setSequence MUST
     * have been called before this function
     * @param startMatch the match containing the SA ranges corresponding to
     * this exact match and depth of the exact match
     * @param beginInPattern the begin position of the exact part in the pattern
     * to be searched
     * @param maxED the maximal allowed edit distance
     * @param occ the occurrences, to this vector new in-text verified
     * occurrences can be added during the function
     * @param counters the performace counters
     */
    void verifyExactPartialMatchInText(BitParallelED& intextMatrix,
                                       const FMOcc& startMatch,
                                       const length_t& beginInPattern,
                                       const length_t& maxED, Occurrences& occ,
                                       Counters& counters) {
        // Immediately switch to in-text verification
        counters.immediateSwitch++;

        // A) convert startmatch to occurrences in the text
        const auto& to = convertToMatchesInText(startMatch);

        // B) find out lowest and highest startpostion
        length_t remBefore = beginInPattern;
        length_t lStartDec = remBefore + maxED;
        length_t hStartDec = (remBefore > maxED) ? remBefore - maxED : 0;
        length_t hStartInc = 0;
        inTextVerification(to, maxED, 0, intextMatrix, occ, counters, lStartDec,
                           hStartDec, hStartInc);
    }
    /**
     * Matches a search recursively with a depth first approach (each branch
     * of the tree is fully examined until the backtracking condition is
     * met) using the edit distance metric. This function does not use any
     * optimizations for eliminating redundancy in the edit distance metric.
     * It simply matches the current part starting from startrange and each
     * node found that has an edit distance between the lower and upper
     * bound is used to start a search for the next part
     * @param search, the search to follow
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
                                 BitParallelED& inTextMatrix,
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
};

// ============================================================================
// CLASS Occurrences
// ============================================================================
// This class combines the occurrences in the text and the occurrences
// in the fm index into one datastructure
class Occurrences {
  private:
    std::vector<TextOcc> inTextOcc; // the in-text occurrences
    std::vector<FMOcc> inFMOcc;     // the in-index occurrences

    /**
     * Erase all double in-index occurrences and sorts the occurrences
     */
    void eraseDoublesFM() {
        sort(inFMOcc.begin(), inFMOcc.end());
        inFMOcc.erase(unique(inFMOcc.begin(), inFMOcc.end()), inFMOcc.end());
    }

    /**
     * Erase all double in-text occurrences and sorts the occurrences
     */
    void eraseDoublesText() {
        sort(inTextOcc.begin(), inTextOcc.end());
        inTextOcc.erase(unique(inTextOcc.begin(), inTextOcc.end()),
                        inTextOcc.end());
    }

  public:
    Occurrences(const length_t reserve = 200) {
        inTextOcc.reserve(reserve);
        inFMOcc.reserve(reserve);
    }

    /**
     * Add an FM occurrence
     */
    void addFMOcc(const FMOcc& match) {
        inFMOcc.emplace_back(match);
    }

    /**
     * Add an FM occurrence
     */
    void addFMOcc(const SARangePair& ranges, const length_t& score,
                  const length_t& depth) {
        inFMOcc.emplace_back(ranges, score, depth);
    }

    /**
     * Add an FM occurrence
     */
    void addFMOcc(const FMPosExt& currentNode, const length_t& score) {
        inFMOcc.emplace_back(currentNode, score);
    }

    /**
     * Add an i-text occurrence
     */
    void addTextOcc(const Range& range, const length_t& score,
                    std::vector<std::pair<char, uint>>& CIGAR) {
        inTextOcc.emplace_back(range, score, CIGAR);
    }

    /**
     * Get all unique text occurrences for the hamming distance
     * @param index the FM index to use
     * @param patternSize the size of the pattern
     * @param counters performance counters
     */
    std::vector<TextOcc> getTextOccurrencesHamming(const FMIndex& index,
                                                   length_t patternSize,
                                                   Counters& counters) {

        eraseDoublesFM();

        counters.totalReportedPositions += inTextOcc.size();

        for (const auto& fmOcc : inFMOcc) {
            const Range saRange = fmOcc.getRanges().getRangeSA();
            counters.totalReportedPositions += saRange.width();
            for (length_t i = saRange.getBegin(); i < saRange.getEnd(); i++) {
                length_t b = index.findSA(i);

                std::vector<std::pair<char, uint>> CIGAR = {
                    std::make_pair('M', patternSize)};
                inTextOcc.emplace_back(Range(b, b + patternSize),
                                       fmOcc.getDistance(), CIGAR);
            }
        }

        // remove doubles
        eraseDoublesText();
        for (auto& t : inTextOcc) {
            t.generateOutput();
        }

        return inTextOcc;
    }

    /**
     * Get all unique text occurrences for the edit distance
     * @param index the FM index to use
     * @param maxED the maximum allowed edit distance, needed for filtering
     * redundant occurrences
     * @param patternMatrix the bit parallel matrix needed to find the CIGAR
     * strings of the in-index occurrences
     * @param counters performance counters
     */
    std::vector<TextOcc> getUniqueTextOccurrences(const FMIndex& index,
                                                  const length_t& maxED,
                                                  BitParallelED& patternMatrix,
                                                  Counters& counters) {

        // increment reporte position counter
        counters.totalReportedPositions += inTextOcc.size();

        // erase equal occurrences from the in-index occurrences
        eraseDoublesFM();

        // convert the in-index occurrences to in-text occurrences
        for (const auto& fmocc : inFMOcc) {

            const Range saRange = fmocc.getRanges().getRangeSA();

            // increment reported positions counter
            counters.totalReportedPositions += saRange.width();

            for (length_t i = saRange.getBegin(); i < saRange.getEnd(); i++) {
                // find the startPosition in the text by looking at the
                // SA
                length_t startPos = index.findSA(i) + fmocc.getShift();

                inTextOcc.emplace_back(
                    Range(startPos, startPos + fmocc.getDepth()),
                    fmocc.getDistance());
            }
        }

        // erase equal occurrences from the in-text occurrences, note
        // that an in-text occurrence with calculated CIGAR string takes
        // preference over an equal one without CIGAR string
        eraseDoublesText();

        // find the non-redundant occurrences
        std::vector<TextOcc> nonRedundantOcc;
        nonRedundantOcc.reserve(inTextOcc.size());

        length_t maxDiff = 2 * maxED;
        length_t prevBegin = std::numeric_limits<length_t>::max();
        length_t prevDepth = std::numeric_limits<length_t>::max();

        length_t prevED = maxED + 1;

        for (const auto& o : inTextOcc) {
            // find the difference between this and the previous
            // occurrence
            auto diff = abs_diff<length_t>(o.getRange().getBegin(), prevBegin);

            if (diff == 0) {
                // same location -> skip
                continue;
            }

            if (diff <= maxDiff) {
                // check if this later occurrence is better than the
                // previous one
                if (o.getDistance() > prevED) {
                    continue;
                }
                if (o.getDistance() == prevED &&
                    o.getRange().width() >= prevDepth) {
                    continue;
                }

                // prev was worse so pop_back
                nonRedundantOcc.pop_back();
            }

            prevBegin = o.getRange().getBegin();
            prevED = o.getDistance();
            prevDepth = o.getRange().width();

            nonRedundantOcc.emplace_back(o);
        }

        for (TextOcc& occ : nonRedundantOcc) {
            if (!occ.hasCigar()) {
                // this was an in-index occurrence which needs to
                // calculate the CIGAR string find the reference
                // sequence
                Substring ref = index.getSubstring(occ.getRange());
                // calculate the cigar string
                std::vector<std::pair<char, uint>> CIGAR;
                patternMatrix.findCIGAR(ref, occ.getDistance(), CIGAR);
                occ.setCigar(CIGAR);

                counters.cigarsInIndex++;

            } else {
                // this was a useful in-text cigar
                counters.usefulCigarsInText++;
            }
            // convert the information to text-format
            occ.generateOutput();
        }

        return nonRedundantOcc;
    }

    length_t getMaxSize() const {
        return std::max(inTextOcc.size(), inFMOcc.size());
    }
};

#endif

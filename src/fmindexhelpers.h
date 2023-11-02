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

#ifndef FMINDEXHELPERS_H
#define FMINDEXHELPERS_H

#include "bitparallelmatrix.h"
#include "substring.h"
#include <sstream>

// ============================================================================
// (TYPE) DEFINITIONS AND PROTOTYPES
// ============================================================================

#include "wordlength.h"

// ============================================================================
// HELPERS
// ============================================================================

// compute |a-b| in a safe manner
template <typename T> T abs_diff(T a, T b) {
    return a > b ? a - b : b - a;
}

#define MAX_MAPQ 60

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
 * Operator overloading. Outputs the range as [begin, end) of theoutput stream
 * @param output, the output stream
 * @param r, the range toprint
 */
std::ostream& operator<<(std::ostream& output, const Range& r);

// ============================================================================
// CLASS CIGAR
// ============================================================================

class CIGARString : public std::vector<std::pair<char, uint>> {
  public:
    friend std::ostream& operator<<(std::ostream& o, const CIGARString& c) {
        for (const auto& p : c) {
            o << p.second << p.first;
        }
        return o;
    }
    CIGARString(std::vector<std::pair<char, uint>>& super)
        : std::vector<std::pair<char, uint>>(super) {
    }
    CIGARString() : std::vector<std::pair<char, uint>>() {
    }
};

// ============================================================================
// CLASS TextOccurrence
// ============================================================================
class TextOcc {
  private:
    Range range;       // the range in the text
    length_t distance; // the distance to this range (edit or hamming)
    CIGARString CIGAR; // The CIGAR string of the match

    std::string samLine; // the corresponding output for this occurrence (for
    // now a custom format)

    uint16_t getFlagsSE(bool reverseComplement, bool notPrimary) const {
        uint16_t result = 0;

        // set rev complement flag
        result |= (reverseComplement << 4);
        // set not primary flag
        result |= (notPrimary << 8);

        return result;
    }

    length_t getMapQ(length_t score, length_t hitNumbers,
                     length_t minScore) const {
        if (score != minScore) {
            return 0;
        }
        if (hitNumbers == 1) {
            return MAX_MAPQ;
        } else {
            return round(-10.0 * log10(1 - 1.0 / hitNumbers));
        }
    }

  public:
    /**
     * Constructor
     * @param range, the range of this occurrence in the text
     * @param distance, the (edit or hamming) distance to the mapped read of
     * this occurrence
     */
    TextOcc(Range range, length_t distance)
        : range(range), distance(distance), CIGAR(), samLine() {
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
        : range(range), distance(distance), CIGAR(CIGAR), samLine() {
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
        samLine = std::to_string(range.getBegin()) + "\t" +
                  std::to_string(range.width()) + "\t" +
                  std::to_string(distance) + "\t";
        for (const auto& p : CIGAR) {
            samLine += std::to_string(p.second) + p.first;
        }
    }

    void generateSamSE(bool revCompl, bool notPrimary, const std::string& seqID,
                       const std::string& seq, const std::string& qual,
                       length_t nHits, length_t minscore) {
        std::stringstream s;

        s << seqID << "\t";                            // read name
        s << getFlagsSE(revCompl, notPrimary) << "\t"; // flags
        s << "*\t";                        // reference sequence name
        s << range.getBegin() + 1 << "\t"; // 1-based pos in ref seq
        s << getMapQ(distance, nHits, minscore) << "\t"; // mapping quality
        s << CIGAR << "\t";                              // CIGAR string
        s << "*\t";       // mate ref seq name, always set to *
        s << "0\t";       // mate pos, always set to 0
        s << "0\t";       // inferred insert size, always set to 0
        s << seq << "\t"; // read sequence
        s << (qual.empty() ? "*" : qual) << "\t"; // read quality
        s << "AS:i:" << distance << "\t";         // alignment score
        s << "PG:Z:Columba";                      // program name
        samLine = s.str();
    }

    const Range getRange() const {
        return range;
    }
    const length_t getDistance() const {
        return distance;
    }

    const std::string& getSamLine() const {
        return samLine;
    }

    /**
     * Operator overloading for sorting the occurrences.
     * Occurrences are first sorted on their begin position, then on their
     * distance , their length and finally on the existence of the CIGAR string
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
            // begins of ranges are unequal, return smallest begin
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
     * Operatoroverloading
     * Two FMocc are equal if their ranges, distance, depth and shift are all
     * equal
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
     * @param row the row of this node in thealignment matrix = depth of
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

            // if finalPiece, report only once
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
                           std::vector<uint>& initEds);
};

// ============================================================================
// STRUCT COUNTERS
// ============================================================================
// A struct of performance counters
struct Counters {
    // performance counters
    uint64_t nodeCounter; // counts the number of nodes visited in the index

    uint64_t totalReportedPositions; // counts the number of matches
                                     // reported (either via in-text
                                     // verification or in-index matching)

    uint64_t cigarsInIndex; // counts the number of cigar strings calculated
                            // for matches in the index, note that this is
                            // only calculated if the match is non-redundant

    uint64_t inTextStarted; // counts the number of times in-text verification
                            // was started, this equals the number of look-ups
                            // in the suffix array for in-text verification
    uint64_t abortedInTextVerificationCounter; // counts the number of
                                               // unsuccesful in-text
                                               // verifications

    uint64_t cigarsInTextVerification; // counts the number of cigars strings
                                       // calculated for matches in the text,
                                       // note that this is done for each match
                                       // in the text as at the point of
                                       // calculation it is not known if this
                                       // match will turn out to be redundant
    uint64_t
        usefulCigarsInText; // counts the number of cigar strings calculated
                            // for non-redundant matches in the text

    uint64_t immediateSwitch; // Counts the number of times the partial
                              // matches after the first part has been
                              // matched are immediately in-text verified
    uint64_t approximateSearchStarted; // Counts the number of times a
                                       // search does start

    /**
     * Reset all counters to 0
     */
    void resetCounters() {
        nodeCounter = 0, abortedInTextVerificationCounter = 0,
        totalReportedPositions = 0, cigarsInIndex = 0,
        cigarsInTextVerification = 0, inTextStarted = 0, usefulCigarsInText = 0,
        immediateSwitch = 0, approximateSearchStarted = 0;
    }

    Counters() {
        resetCounters();
    }

    void addCounters(const Counters& o) {
        nodeCounter += o.nodeCounter,
            abortedInTextVerificationCounter +=
            o.abortedInTextVerificationCounter,
            totalReportedPositions += o.totalReportedPositions,
            cigarsInIndex += o.cigarsInIndex,
            cigarsInTextVerification += o.cigarsInTextVerification,
            inTextStarted += o.inTextStarted,
            usefulCigarsInText += o.usefulCigarsInText,
            immediateSwitch += o.immediateSwitch,
            approximateSearchStarted += o.approximateSearchStarted;
    }
    void incNodeCounter() {
        nodeCounter++;
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
     * Add an in-text occurrence
     */
    void addTextOcc(const Range& range, const length_t& score,
                    std::vector<std::pair<char, uint>>& CIGAR) {
        inTextOcc.emplace_back(range, score, CIGAR);
    }

    /**
     * Add an in-text occurrence
     */
    void addTextOcc(const Range& range, const length_t score) {
        inTextOcc.emplace_back(range, score);
    }

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

    const std::vector<FMOcc>& getFMOccurrences() const {
        return inFMOcc;
    };

    size_t textOccSize() const {
        return inTextOcc.size();
    }

    void generateOutput() {
        for (auto& t : inTextOcc) {
            t.generateOutput();
        }
    }

    const std::vector<TextOcc>& getTextOccurrences() const {
        return inTextOcc;
    }

    length_t getMaxSize() const {
        return std::max(inTextOcc.size(), inFMOcc.size());
    }
};

// ============================================================================
// CLASS INTEXTVERIFICATIONTASK
// ============================================================================

class IntextVerificationTask {

  private:
    const std::vector<Substring>&
        refs;              // the reference subsequence to be checked
    BitParallelED& matrix; // the matrix to use
    const length_t maxED;  // maximum allowed edit distance
    const length_t minED;  // minimum allowed edit distance

  public:
    IntextVerificationTask(const std::vector<Substring>& refs,
                           BitParallelED& matrix, const length_t maxED,
                           const length_t minED)
        : refs(refs), matrix(matrix), maxED(maxED), minED(minED) {
    }

    /**
     * Verify the ref subsequence to the matrix
     * @param counters performance counters
     * @param occ the occurrences, new in-text occurrences may be added
     */
    void doTask(Counters& counters, Occurrences& occ);
};

#endif
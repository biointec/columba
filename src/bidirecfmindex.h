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
#ifndef BIDIRECFMINDEX_H
#define BIDIRECFMINDEX_H

#include "tkmer.h"
#include <google/sparse_hash_map>

#include "fmindex.h"
#include <chrono>

class Node {
  public:
    char c;                // the character of this node
    SARangePair ranges;    // the ranges of this node
    int row;               // the row of this node in this search
    bool reported = false; // has this particular node already reported?

    /**
     * Create a node of the search tree
     * @param character the character of this node
     * @param ranges the ranges over the suffix and reversed suffix array
     * that go to this node
     * @param row the row of this node in the allignment matrix = depth of
     * this node
     */
    Node(const char& character, const SARangePair& ranges, const length_t& row)
        : c(character), ranges(ranges), row(row), reported(false) {
    }

    /**
     * Default constructor, this Node will have empty ranges
     */
    Node() : c(char(0)), ranges(SARangePair()), row(0), reported(false) {
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
     * @param noDoubleReports false if this node is allowed to report more than
     once, defaults to false
     */
    void report(BiAppMatchSA& occ, const length_t& startDepth,
                const length_t& EDFound, const bool& noDoubleReports = false) {
        if (!reported) {
            occ = makeBiAppMatchSA(ranges, EDFound, row + startDepth);
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
    unsigned int getRow() const {
        return row;
    }
};

class Cluster {
  private:
    std::vector<int> eds;    // the edit distances of this cluster
    std::vector<Node> nodes; // the nodes of this cluster

    int lastCell;   // the lastCell of the cluster that was filled in
    int maxED;      // the maxEd for this cluster
    int startDepth; // the startdepth for this cluster

  public:
    /**
     * Constructor
     * @param size, the size of the cluster
     * @param maxED, the maximal allowed edit distance
     * @param startDepth, the depth before this cluster
     */
    Cluster(int size, int maxED, int startDepth)
        : eds(size, maxED + 1), nodes(size), lastCell(-1), maxED(maxED),
          startDepth(startDepth) {
    }

    /**
     * Sets the ed and node at index idx to ed and node. Also updates lastCell
     * to be idx
     * @param idx, the idx to change
     * @param node, the node to set at index idx
     * @param ed, the ed to set at index idx
     */
    void setValue(int idx, const Node& node, const int& ed) {
        eds[idx] = ed;
        nodes[idx] = node;
        lastCell = idx;
    }

    /**
     * Returns the size of this cluster
     */
    const unsigned int size() const {
        return eds.size();
    }

    /**
     * @returns approximate match that corresponds to the highest minimum of
     * this cluster (the highest node in the tree that has the minimal ed of
     * this cluster). If no node corresponds to an ed lower than or eqaul to
     * maxED, an empty match is returned
     *
     */
    BiAppMatchSA reportHighestMinimum() {
        int minED = maxED + 1;
        int bestIdx = -1;

        for (int i = 0; i <= lastCell; i++) {
            if (eds[i] < minED) {
                minED = eds[i];
                bestIdx = i;
            }
        }
        BiAppMatchSA m;
        if (bestIdx > -1) {
            nodes[bestIdx].report(m, startDepth, eds[bestIdx], true);
        }
        return m;
    }
    /**
     * @returns approximate match that corresponds to the deepest minimum of
     * this cluster (the deepest node in the tree that has the minimal ed of
     * this cluster). If no node corresponds to an ed lower than or eqaul to
     * maxED, an empty match is returned
     *
     */
    BiAppMatchSA reportDeepestMinimum() {
        int minED = maxED;
        int bestIdx = -1;

        for (int i = 0; i <= lastCell; i++) {
            if (eds[i] <= minED) {
                minED = eds[i];
                bestIdx = i;
            }
        }
        BiAppMatchSA m;
        if (bestIdx > -1) {
            nodes[bestIdx].report(m, startDepth, eds[bestIdx], true);
        }
        return m;
    }

    /**
     * This method returns a match that corresponds to the highest cluster
     * centre. Its descendants and the corresponding initalization eds are
     * updated. Eds of descendants that are part of a cluster centre which is
     * lower than the lowerbound will be updated in the initEds vector
     * @param lowerBound, the lowerbound for this iteration
     * @param desc, the descendants of the highest cluster centre, these will be
     * inserted during the method
     * @param initEds, the initialization eds for the next iteration, these
     * correspond to the eds of the highest centre and its descendants, where
     * eds part of a cluster of which the centre is below the lowerbound are
     * updated. These values will be inserted during the method
     */
    BiAppMatchSA getClusterCentra(int lowerBound, std::vector<Node>& desc,
                                  std::vector<int>& initEds) {
        desc.reserve(eds.size());
        initEds.reserve(eds.size());
        BiAppMatchSA m;
        for (int i = 0; i <= lastCell; i++) {
            if (eds[i] > maxED || eds[i] < lowerBound) {
                continue;
            }
            bool betterThanParent = (i == 0) || eds[i] <= eds[i - 1];
            bool betterThanChild = (i == lastCell) || eds[i] <= eds[i + 1];

            if (betterThanParent && betterThanChild) {
                // this is a valid centre
                nodes[i].report(m, startDepth, eds[i]);

                // get all the descendants
                initEds.push_back(eds[i]);
                for (int j = i + 1; j <= lastCell; j++) {
                    desc.push_back(nodes[j]);
                    initEds.push_back(eds[j]);
                }

                // replace the clusters under the lowerbound
                for (unsigned int k = 1; k < initEds.size(); k++) {
                    if (initEds[k] < lowerBound &&
                        initEds[k] <= initEds[k - 1] &&
                        (k == initEds.size() - 1 ||
                         initEds[k] <= initEds[k + 1])) {
                        // this is a centre under the lowerbound

                        unsigned int highestPoint = k;
                        unsigned int lowestPoint = initEds.size() - 1;
                        // find highest point of this cluster
                        for (unsigned int l = k - 1; l > 0; l--) {
                            if (initEds[l] != initEds[l - 1] - 1) {
                                highestPoint = l;
                                break;
                            }
                        }
                        // replace values higher than centre
                        for (unsigned int l = highestPoint; l < k; l++) {
                            initEds[l] = initEds[l - 1] + 1;
                        }

                        // find lowest point of this cluster
                        for (unsigned int l = k; l < initEds.size() - 1; l++) {
                            if (initEds[l] != initEds[l + 1] - 1) {
                                lowestPoint = l;
                                break;
                            }
                        }

                        // replace values below centre
                        int increase = -1;
                        if (lowestPoint == initEds.size() - 1) {
                            increase = 1;
                        }
                        for (unsigned int l = k; l <= lowestPoint; l++) {
                            initEds[l] = initEds[l - 1] + increase;
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
class BidirecFMIndex;
typedef bool (BidirecFMIndex::*ExtraCharPtr)(length_t, const SARangePair&,
                                             SARangePair&) const;

class BidirecFMIndex : public FMIndex {
  private:
    BWTRepr<ALPHABET> revRepr; // the prefix occurences of the rev BWT
    Direction dir = BACKWARD;  // the direction of the index
    ExtraCharPtr extraChar;    // pointer to extra char method (for direction)

    std::vector<std::vector<Node>>
        stacks; // stacks of nodes for the different partitions

    static const size_t wordSize =
        4; // the size the mers to be stored in a table
    google::sparse_hash_map<Kmer, SARangePair, KmerHash>
        table; // hashtable that contains all wordSize-mers

    // ----------------------------------------------------------------------------
    // PREPROCESSING ROUTINES
    // ----------------------------------------------------------------------------

    /**
     * Populate the hash table
     */
    void populateTable();

    // ----------------------------------------------------------------------------
    // ROUTINES FOR ACCESSING DATA STRUCTURE
    // ----------------------------------------------------------------------------

    /**
     * Finds the LF mapping of the character at index k in the bwt string
     * @param k the index to find the LF mapping off
     * @param reversed boolean to indicate if the reversed bwt should be
     * used or the normal bwt
     * @returns the row that is the LF mapping of k. It is so that the entry
     * in the suffix array of this return value is one less than the entry
     * in the sufffix array at index k
     */
    length_t findLF(length_t k, bool reversed) const;

    /**
     * Function that returns the nummber of occurences before an index of
     * the symbol at symbolindex in the alphabet
     * @param symbolIndex the index of the the symbol in the alphabet to
     * count the occrences of at index index
     * @param index the index whose entry for symbol in the occurences table
     * is asked
     * @return the number of occurences of the symbol before index in the
     * bwt
     */
    length_t getNumberOfOcc(length_t symbolIndex, length_t index) const {
        return fwdRepr.occ(symbolIndex, index);
    }
    /**
     * Same as BiBWT::getNumberOfOccurences, but now in the bwt of the
     * reversed text
     * @param symbolIndex the index in the alphabet to get the number of
     * occurences of
     * @param index the index in the occurences table, for which the number
     * of occurences of the symbol at symbolindex is asked.
     * @return the number of occurences at index index in the occurences
     * table of the bwt of the reversed text for the symbol at symbolIndex
     * in the alphabet
     */
    length_t getNumberOfOccRev(length_t symbolIndex, length_t index) const {
        return revRepr.occ(symbolIndex, index);
    }

    /**
     * Function that returns the nummber of occurences before the index of
     * all symbols smaller than the symbol at symbolindex in the  bwt
     * @param symbolIndex the index in the alphabet whose number of prefix
     * occurences is queried.
     * @param index the index whose entry for symbol in the prefixoccurences
     * table is asked
     * @return the number of occurences of symbols smaller than symbol at
     * symbolindex before index index in the bwt
     */
    length_t getNumberOfCumOcc(length_t symbolIndex, length_t index) const {
        return fwdRepr.cumOcc(symbolIndex, index);
    }

    /**
     * Function that returns the nummber of occurences before the index of
     * all symbols smaller than the symbol at symbolindex in the  bwt of the
     * reversed text
     * @param symbolIndex the index in the alphabet whose number of prefix
     * occurences is queried.
     * @param index the index whose entry for symbol in the
     * rprefixoccurences table of the reverse text is asked
     * @return the number of occurences of symbols smaller than symbol at
     * symbolindex before index index in the bwt
     */
    length_t getNumberOfCumOccRev(length_t symbolIndex, length_t index) const {
        return revRepr.cumOcc(symbolIndex, index);
    }

    // ----------------------------------------------------------------------------
    // HELP ROUTINES FOR APPROXIMATE PATTERN MATCHING
    // ----------------------------------------------------------------------------

    /**
     * Private helper function that gets all possible next characters and
     * their ranges given the ranges of the current pattern
     * @param rangesOfParent the ranges of the parent in the bwt
     * @returns a vector with pairs of ranges and charachters, these are the
     * ranges of all possible previous/next characters
     */
    std::vector<std::pair<SARangePair, char>>
    getCharExtensions(const SARangePair& rangesOfParent) const;

    /**
     * Finds the ranges of cP using the principle explained in the paper of Lahm
     * @param positionInAlphabet the postition in alphabet of the character
     * that is added in the front
     * @param rangesOfP the ranges of pattern P
     * @returns the ranges cP
     */
    bool findRangesWithExtraCharBackward(length_t positionInAlphabet,
                                         const SARangePair& rangesOfP,
                                         SARangePair& childRanges) const;

    /**
     * Finds the ranges of Pc using the principle explained in the paper of Lahm
     * @param positionInAlphabet the postition in alhabet of the character c
     * that is added in the back
     * @param rangesOfP the ranges of pattern P
     * @returns the ranges of Pc
     */
    bool findRangesWithExtraCharForward(length_t positionInAlphabet,
                                        const SARangePair& rangesOfP,
                                        SARangePair& childRanges) const;

    // ----------------------------------------------------------------------------
    // HELPER ROUTINES FOR APPROXIMATE MATCHING (ITERATIVELY)
    // ----------------------------------------------------------------------------
    /**
     * Goes deeper in a search if a valid approximate match is found in the
     * cluster
     * @param cluster, the cluster to search for a valid approximate
     * @param maxED, the maximal allowed edit distance for this partition
     * @param startMatch, the original match for this partition
     * @param nextP, the next partion to research
     * @param s, the search
     * @param occ, the vector with all occurrences of the entrie pattern, if the
     * current partition is the final partition of the search then the match (if
     * one found) will be pushed onto this vector
     * @param lowerbound, the lowerbound for this partition
     * @param descendantsOtherD, the desencdants of the other direction,
     * defaults to empty vector
     * @param ininEdsOtherD, the initialization eds of the other direction,
     * defaults to empty vector
     */
    void goDeeper(Cluster& cluster, const length_t& maxED,
                  const BiAppMatchSA& startMatch, const length_t& nextp,
                  const Search& s, const std::vector<Substring>& parts,
                  std::vector<BiAppMatchSA>& occ, const length_t& lowerBound,
                  const std::vector<Node>& descendantsOtherD = {},
                  const std::vector<int>& initEdsOtherD = {});

    /**
     * Helper function for the approximate matching. This function fills in the
     * matrix for the current node at the current row and goes deeper for the
     * next part is that is necessary
     * The function returns zero if no branching is needed, 1 if the search went
     * deeper and -1 if maxED was exceeded
     * @param matrix, the matrix to fill in
     * @param clust, the cluster corresponding to the final column of the matrix
     * @param row, the row of the matrix to fill in
     * @param node, the node for which the matrix is filled in
     * @param idx, the idx of the current part
     * @param parts, the parts of the pattern
     * @param maxED, the curren allowed max ed
     * @param startMatch, the match for which this part started
     * @param occ, vector to push occurrences to
     * @param initOther, eds of descendants in other direction
     * @param descOther, descendants in other direction
     * @return false if the search can continue along this branch, true if the
     * search can backtrack
     */
    bool branchAndBound(BandMatrix& matrix, Cluster& clus, const int row,
                        const Node& currentNode, const Search& s,
                        const int& idx, const std::vector<Substring>& parts,
                        const length_t& maxED, const BiAppMatchSA& startMatch,
                        std::vector<BiAppMatchSA>& occ,
                        const std::vector<int>& initOther,
                        const std::vector<Node>& descOther);

    /**
     * Pushes all the children corresponding to the node with ranges equal to
     * ranges
     * @param ranges the ranges to get the children of
     * @param stack, the stack to push the children on
     * @param row, the row of the parentNode (defaults to 0)
     */
    void pushChildren(const SARangePair& ranges, std::vector<Node>& stack,
                      int row = 0);

  public:
    // ----------------------------------------------------------------------------
    // INITIALIZATION ROUTINES
    // ----------------------------------------------------------------------------

    /**
     * Constructor
     * @param prefix prefix of the files that contain the info
     * @param sa_spase sparseness factor of suffix array. It is assumed this
     * is a power of two
     */
    BidirecFMIndex(const std::string& prefix, int sa_sparse = 1)
        : FMIndex(prefix, sa_sparse) {

        // read the reverse prefix occurrence table
        if (!revRepr.read(prefix + ".rev.brt"))
            throw std::runtime_error("Cannot open file: " + prefix +
                                     ".rev.brt");
        std::cout << "Done reading reverse prefix occurrence table"
                  << std::endl;
        extraChar = &BidirecFMIndex::findRangesWithExtraCharBackward;
        populateTable();
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

    length_t getNodes() const {
        return nodeCounter;
    }

    /**
     * @returns the wordsize of the mers stored in the table
     */
    int getWordSize() const {
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
     * Looks up the SARangePair corresponding to p in the hashtable. Assumes p
     * is of size wordSize
     * @param p, the substring to find the ranges of
     * @returns the ranges corresponding to substring p, if no pair can be found
     * returns empty ranges
     */
    SARangePair lookUpInKmerTable(const Substring& p) const {
        Kmer k(p.getSubstring());

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
     * This function matches a string exactly and returns the ranges in the
     * sa and saRev
     * @param pattern the string to match
     * @returns the pair of ranges of this pattern
     */
    SARangePair matchStringBidirectionally(const Substring& pattern) {
        return matchStringBidirectionally(pattern, getCompleteRange());
    }

    /**
     * This function matches a string exactly starting form startRange
     * @param pattern the string to match
     * @param startRange, the range to search in
     * @returns the pair of ranges
     */
    SARangePair matchStringBidirectionally(const Substring& pattern,
                                           SARangePair startRange);
    /**
     * Adds one character and updates the range. If the character can't be
     * added the range will be set to an empty range
     * @param c the character to be added (in the current direction of the
     * index)
     * @param range the range to extend
     */
    bool addChar(const char& c, SARangePair& range);

    // ----------------------------------------------------------------------------
    // PREPROCESSING ROUTINES FOR APPROXIMATE MATCHING
    // ----------------------------------------------------------------------------

    /**
     * Resizes to the required number of stacks and reserves space on each
     * stack, so that each stack can match the entire pattern
     * @param number, the number of stacks required
     * @param size, the size of the pattern
     */
    void reserveStacks(const length_t number, const int size) {
        stacks.resize(number);
        for (auto& stack : stacks) {
            stack.reserve(size * sigma.size());
        }
    }

    // ----------------------------------------------------------------------------
    // ROUTINES FOR APPROXIMATE MATCHING
    // ----------------------------------------------------------------------------

    /**
     * Sets the search direction of the fm-index
     * @param d the direction to search in, either FORWARD or BACKWARD
     */
    void setDirection(Direction d) {
        dir = d;
        extraChar = (d == FORWARD)
                        ? &BidirecFMIndex::findRangesWithExtraCharForward
                        : &BidirecFMIndex::findRangesWithExtraCharBackward;
    }

    /**
     * Matches a search recursively with a depth first approach (each branch of
     * the tree is fully examined untill the backtracking condition is met)
     * @param search, the search to folloow
     * @param startMatch, the approximate match found for all previous partions
     * of the search
     * @param occ, a vector with matches of the complete search, if a such a
     * match is found is a pushed upon this vector
     * @param idx, the index of the partition to match, defaults to 1 as an
     * exact search for the zeroth partition is assumed
     * @param descPrevDir, the descendants of the previous direction, defaults
     * to empty vector
     * @param initPrevDir, the initialization eds of the previous direction,
     * defaults to empty vector
     * @param descNotPrevDir, the descendants of the other direction, defaults
     * to empty vector
     * @param initNotPrevDir, the initialization eds of the other direciton,
     * defaults to empty vector
     */
    void recApproxMatch(
        const Search& search, const BiAppMatchSA& startMatch,
        std::vector<BiAppMatchSA>& occ, const std::vector<Substring>& parts,
        const int& idx = 1,
        const std::vector<Node>& descPrevDir = std::vector<Node>(),
        const std::vector<int>& initPrevDir = std::vector<int>(),
        const std::vector<Node>& descNotPrevDir = std::vector<Node>(),
        const std::vector<int>& initNotPrevDir = std::vector<int>());

    // ----------------------------------------------------------------------------
    // POST-PROCESSING ROUTINES FOR APPROXIMATE MATCHING
    // ----------------------------------------------------------------------------

    /**
     * This function maps matches in the sa to matches in the text. It takes
     * the ranges of the matches and together with the depth this is matched
     * to a range in the text (this new range has a width of depth). The
     * edit distance is also mapped to this new range in the text. This function
     * also filters out the redundant matches using the maximal allowed edit
     * distance
     * @param occurrences, the vector with all approximate matches and their
     * ranges in the SA and revSA
     * @param maxED, the maximal allowed edit distance
     * @returns a vector of matches in the text containing a range and the
     * edit distance
     */
    std::vector<AppMatch>
    mapOccurencesInSAToOccurencesInText(std::vector<BiAppMatchSA>& occurences,
                                        const int& maxED);
};

#endif

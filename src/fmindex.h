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
#ifndef FMINDEX_H
#define FMINDEX_H

#include "bandmatrix.h"

#include "customtypedefs.h"
#include <algorithm> //used for sorting
#include <cstdint>
#include <fstream>  // used for reading in files
#include <iostream> // used for printing
#include <map>      // used to store the static alphabet
#include <math.h>   //for taking the log
#include <numeric>  // for summing over vector
#include <set>
#include <sstream> // used for splitting strings
#include <string>
#include <vector>

// ============================================================================
// CLASS FMIndex
// ============================================================================

//

class FMIndex {
  protected:
    // define the special end character
    static const char DOLLAR = '$';

    // alphabet that is used (ASCII or DNA)

    Alphabet<ALPHABET> sigma;
    length_t sparseFactorSA =
        32; // the sparseness factor of the suffix array, defaults to 32

    int logSparseFactorSA = 5; // the log of the sparse factor

    std::string bwt;            // the bwt string of the reference genome
    std::vector<size_t> counts; // the counts array of the reference genome
    std::vector<length_t>
        sa; // the (sparse) suffix array of the reference genome
    BWTRepr<ALPHABET> fwdRepr; // the prefix occurences table

    const std::string
        prefix; // remember the prefix for retrieving the text again
    length_t textLength = 0;

    length_t nodeCounter = 0; // a nodeCounter for counting how many nodes are
                              // visited during a specific matching procedure
    length_t matrixElementCounter =
        0; // a counter for counting how many matrix elements are filled in
    length_t positionsInPostProcessingCounter =
        0; // a counter that checks how many start positions were reported
           // (= length of all reported ranges over the suffix array)

    // ----------------------------------------------------------------------------
    // INITIALIZATION ROUTINES
    // ----------------------------------------------------------------------------
    /**
     * Private helper function that reads in all the necessary files
     * @param prefix the prefix of the files that will be read in
     */
    void fromFiles(const std::string& prefix);

    // ----------------------------------------------------------------------------
    // ROUTINES FOR ACCESSING DATA STRUCTURE
    // ----------------------------------------------------------------------------

    /**
     * Function that returns the nummber of occurences before the index of
     * all symbols smaller than this symbol in the bwt
     * @param symbol the symbol to count the occrences of smaller symbols of
     * at index index
     * @param index the index whose entry for symbol in the prefixoccurences
     * table is asked
     * @return the number of occurences of symbols smaller than symbol
     * before index in the bwt
     */
    length_t getNumberOfPrefOcc(char symbol, length_t index) const;

    /**
     * Function that returns the nummber of occurences before the index of
     * this symbol in the bwt
     * @param symbol the symbol to count the occrences of at index index
     * @param index the index whose entry for symbol in the occurences table
     * is asked
     * @return the number of occurences of the symbol before index in the
     * bwt
     */
    length_t getNumberOfOcc(char symbol, length_t index) const;

    /**
     * Finds the LF mapping of the character at index k in the bwt string
     * @param k the index to find the LF mapping off
     * @returns the row that is the LF mapping of k. It is so that the entry
     * in the suffix array of this return value is one less than the entry
     * in the sufffix array at index k
     */
    length_t findLF(length_t k) const;

    // ----------------------------------------------------------------------------
    // ROUTINES FOR PATTERN MATCHING
    // ----------------------------------------------------------------------------

    /**
     * Private helper function for ReadMapper::exactMatches. This finds the
     * range in the suffix array that corresponds to matches in the
     * reference genome
     * @param s the string to be matched in the reference genome
     * @returns a Range containing the start and end values of the range
     * ([start, end[)
     */
    Range matchString(const std::string& s);

    // ----------------------------------------------------------------------------
    // ROUTINES FOR APPROXIMAE PATTERN MATCHING
    // ----------------------------------------------------------------------------

    /**
     * Private helper function that gets all possible next characters and
     * their ranges given the ranges of the current pattern
     * @param rangesOfParent the ranges of the parent in the bwt
     * @returns a vector with pairs of ranges and charachters, these are the
     * ranges of all possible previous/next characters
     */
    std::vector<std::pair<Range, char>> getCharExtensions(const Range& range);

    /**
     * Private helper function that recursively calculates which ranges in
     * the suffix array correspond to approximate matches to a pattern.
     * @param range the initial range we are considering, this is a range in
     * the sorted text
     * @param P the pattern to match
     * @param M the banded matrix for keeping track of the edit distance
     * @param bestED the bestED untill this position
     * @param occ a vector of found matches, these contain ranges in the
     * suffix array and their ED and depth
     * @param depth how many iterations (= length of matched substring
     * untill now)
     */
    void recApproxMatchesNaive(const Range& range, const std::string& P,
                               EditMatrix& M, const int& maxED,
                               std::vector<AppMatchSA>& occ, length_t depth);

    /**
     * Private helper function for the recApproxMatchesNaive, this functions
     * goes over all children (all possible prefixes) of the current matched
     * pattern. It will update the chdED and chdReported vectors during the
     * execution. This function will call the recAppproxMatches for the
     * children
     * @param nextChar a vector with all the children that need to be check
     * (a child has a range in the alphabetic ordering and a char)
     * @param chED a vector containing the edit distances found at the
     * children, this will be update during this function
     * @param chdRepoted a vector containg bools that tell if a child found
     * a match
     * @param P the pattern to match
     * @param grandParentED the bestED found at the grandparent
     * @param M the BandMatrix filled in for all rows previous to these
     * children
     * @param row the rownumber of these children in the matrix (this is
     * also the length of the matched substring )
     * @param occ the vector with all matches that have been found
     */
    void checkChildren(std::vector<std::pair<Range, char>>& nextChar,
                       std::vector<EditDistance>& chdED,
                       std::vector<bool>& chdReported, const std::string& P,
                       const EditDistance& grandparentED, EditMatrix& M,
                       length_t row, std::vector<AppMatchSA>& occ);

    // ----------------------------------------------------------------------------
    // POST-PROCESSING ROUTINES FOR APPROXIMATE PATTERN MATCHING
    // ----------------------------------------------------------------------------

    /**
     * Checks the edit distance between a given string and a piece of the
     * original text.
     * @param startPosition the position in the text where the check needs
     * to happen
     * @param endPosition the position just after the piece of the text that
     * needs to be checked, if no restrictions are wanted, set this to the
     * size of the text
     * @param remainingED the maximal allowed ED for this verification
     * @param stringToCheck the string that will be checked against the text
     * from startPosition to endPosition
     * @returns a match in the text with found ed not greater to remaininged
     * or a match with an empty range and ed higher than remaininged
     */
    AppMatch verifyInText(length_t startPosition, length_t endPosition,
                          int remainingED,
                          const Substring& stringToCheck) const;

    /**
    * Converts a match in the suffix array to matches in the text. Also
    * verifies wheter these matches go over a sentinel character. If that
    is
    * the case, it is checked if without the sentinel (and characters
    before or
    * after it) it still forms a match. This is the match that will be
    added to the list of matches.
    *
    * REMARK: if two sentinels/word delimiting chracter are present in
    match this is incorrect however for genomes the distance between
    sentinels is magnitudes larger than the size of the reads, since entrire
    chromosomes cannot be read in at once
    *
    * @param matchInSA the match that will be converted
    * @param pattern the patter for which this was a match
    * @param ED the max allowed edit distance
    */
    std::vector<AppMatch> convertToMatchesInText(const AppMatchSA& matchInSA);

  public:
    /**
     * Get the counters
     * @returns a tuple (nodeCounter, matrixElementCounter,
     * positionsInPostProcessingCounter)
     */
    std::tuple<length_t, length_t, length_t> getCounters() {
        return std::make_tuple(nodeCounter, matrixElementCounter,
                               positionsInPostProcessingCounter);
    }

    /**
     * Resets the counters
     */
    void resetCounters() {
        nodeCounter = 0;
        matrixElementCounter = 0;
        positionsInPostProcessingCounter = 0;
    }

    /**
     * Constructor with default values for the sparseness factors
     * @param prefix prefix of the files that contain the info
     */
    FMIndex(const std::string& prefix) : prefix(prefix) {
        // read in all the files
        fromFiles(prefix);
    }

    /**
     * Constructor
     * @param prefix prefix of the files that contain the info
     * @param sa_spase sparseness factor of suffix array. It is assumed this
     * is a power of two
     */

    FMIndex(const std::string& prefix, int sa_sparse) : prefix(prefix) {
        // set the correct sparsefactors
        this->sparseFactorSA = sa_sparse;

        this->logSparseFactorSA =
            log2(sa_sparse); // assumed sa_sparse is power of two

        // read in all the files
        fromFiles(prefix);
    }

    /**
     * Get the original text
     */
    std::string getText() const {
        return readString(prefix + ".txt");
    }

    // --------------------------------------------------------------------
    // ROUTINES FOR PATTERN MATCHING
    // --------------------------------------------------------------------

    /**
     * Calculates the positions in the reference genome where exact matches
     * to the argument string start.
     * @param s the string to match in the reference genome
     * @returns a sorted vector containing the start positions of all exact
     * substring matches of s in the reference sequence
     */
    std::vector<length_t> exactMatches(const std::string& s);

    // --------------------------------------------------------------------
    // ROUTINES FOR APPROXIMATE PATTERN MATCHING
    // --------------------------------------------------------------------

    /**
     * Matches the pattern approximately. All matches are at most a certain
     * edit distance away from the pattern
     * @param pattern the pattern to match
     * @param maxED the maximum edit distance
     * @returns a vector with matches which contain a range (the range of
     * the text that matched) and the edit distance this substring is away
     * from the pattern
     */
    virtual std::vector<AppMatch> approxMatches(const std::string& pattern,
                                                int maxED);
    // --------------------------------------------------------------------
    // ROUTINES FOR APPROXIMATE PATTERN MATCHING
    // --------------------------------------------------------------------

    /**
     * Maps the found occurences in the suffix array to actual occurences in
     * the text. It also filters out redundant matches using the maximal allowed
     * edit distance
     * @param occurences the found occurences in the suffix array
     * @param maxED the maximal allowed edit distance
     */
    std::vector<AppMatch> mapOccurencesInSAToOccurencesInText(
        const std::vector<AppMatchSA>& occurences, const int& maxED);

    /**
     * Finds the entry in the suffix array of this index. This is computed
     * from the sparse suffix array and the bwt
     * @param index the index to find the entry in the SA off
     * @returns the entry in the SA of the index
     */
    length_t findSA(length_t index) const;
};

#endif

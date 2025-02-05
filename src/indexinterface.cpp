/******************************************************************************
 *  Columba: Approximate Pattern Matching using Search Schemes                *
 *  Copyright (C) 2020-2022 - Luca Renders <luca.renders@ugent.be> and        *
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
#include "indexinterface.h"
#include "definitions.h"
#include "indexhelpers.h"
#include "logger.h"

#include "search.h" // for Search

#include <assert.h>                 // for assert
#include <cstdint>                  // for uint16_t
#include <fmt/format.h>             // for fmt::to_string
#include <fstream>                  // for operator<<, ifstream, basic_ostream
#include <iterator>                 // for distance
#include <limits>                   // for numeric_limits
#include <memory>                   // for allocator_traits<>::value_type
#include <parallel_hashmap/phmap.h> // phmap::parallel_flat_hash_map
#include <stdexcept>                // for runtime_error
#include <type_traits>              // for __strip_reference_wrapper<>::__type

#ifndef RUN_LENGTH_COMPRESSION
#include <fmt/core.h> // for fmt::format
#endif                // not RUN_LENGTH_COMPRESSION

#ifdef _WIN32 // Windows specific headers
#include <windows.h>
#else // Unix-like system headers
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

using namespace std;

// ============================================================================
// CLASS IndexInterface
// ============================================================================

thread_local Direction IndexInterface::dir = BACKWARD;
thread_local ExtraCharPtr IndexInterface::extraChar;

thread_local vector<vector<FMPosExt>> IndexInterface::stacks;
thread_local vector<BitParallelED64> IndexInterface::matrices;
thread_local vector<BitParallelED128> IndexInterface::matrices128;

thread_local Strand IndexInterface::strand = FORWARD_STRAND;
thread_local PairStatus IndexInterface::pairStatus = FIRST_IN_PAIR;

// ----------------------------------------------------------------------------
// PREPROCESSING ROUTINES
// ----------------------------------------------------------------------------

void IndexInterface::readMetaAndCounts(const string& baseFile, bool verbose) {
    logger.logDeveloper("Reading meta and counts");
    // Open the meta file
    ifstream ifs(baseFile + ".meta");
    if (ifs) {
        // Read the build tag
        length_t tag;
        ifs >> tag;
        if (tag < COLUMBA_BUILD_INDEX_TAG) {
            logger.logWarning(
                "The index was built with an older version of columba_build. "
                "Proceed with caution or re-build the index.");
        } else if (tag > COLUMBA_BUILD_INDEX_TAG) {
            logger.logDeveloper(
                "The index was built with a newer version of columba_build. "
                "Proceed with caution or re-build the index.");
        }
        // read the compiled info
        size_t build_length_t_size;
        ifs >> build_length_t_size;
        if (build_length_t_size != sizeof(length_t)) {
            size_t thisSize = sizeof(length_t) * 8;
            build_length_t_size *= 8;
            throw runtime_error(
                "The index was built with a compiled version that uses " +
                to_string(build_length_t_size) +
                "-bit numbers, while the current programme was compiled "
                "using " +
                to_string(thisSize) +
                "-bit numbers. Recompile the programme with the correct "
                "THIRTY_TWO flag set or rebuild the index.");
        }

        // read the flavor
        string flavour;
        ifs >> flavour;
        if (flavour != COLUMBA_FLAVOUR) {
            throw runtime_error(
                "The index was built with a different flavor of Columba. "
                "\n\tCurrent flavor: " +
                fmt::to_string(COLUMBA_FLAVOUR) +
                "\n\tIndex flavor: " + fmt::to_string(flavour) +
                "\n\tRecompile the programme with the correct flavor or "
                "re-build the index.");
        }

    } else {
        logger.logWarning(
            "The index was built with an older version of columba_build or "
            "the meta data is missing! "
            "Proceed with caution or re-build the index.");
    }

    stringstream ss;
    if (verbose) {
        ss << "Reading " << baseFile << ".cct" << "...";
        logger.logInfo(ss);
    }

    // Read the counts table
    vector<length_t> charCounts(256, 0);
    if (!readArray(baseFile + ".cct", charCounts)) {
        throw runtime_error("Cannot open file: " + baseFile + ".cct");
    }

    length_t cumCount = 0; // Cumulative character counts
    for (size_t i = 0; i < charCounts.size(); i++) {
        if (charCounts[i] == 0)
            continue;
        counts.push_back(cumCount);
        cumCount += charCounts[i];
    }
    textLength = cumCount;
    sigma = Alphabet<ALPHABET>(charCounts);
}

void IndexInterface::readSequenceNamesAndPositions(const string& baseFN,
                                                   bool verbose) {
    stringstream ss;
    if (verbose) {
        ss << "Reading " << baseFile << ".pos...";
        logger.logInfo(ss);
    }
    // read the start positions
    if (!readArray(baseFile + ".pos", startPos)) {
        throw runtime_error("Cannot open file: " + baseFile +
                            ".pos\nIs the "
                            "reference index outdated?");
    }

    if (verbose) {
        ss << "Reading " << baseFile << ".sna...";
        logger.logInfo(ss);
    }
    {
        seqNames.reserve(startPos.size());
        ifstream ifs(baseFile + ".sna", ios::binary);
        if (!ifs) {
            throw runtime_error("Cannot open file: " + baseFile +
                                ".sna\nIs the reference index outdated?");
        }

        while (ifs.good()) {
            size_t len;
            ifs.read(reinterpret_cast<char*>(&len), sizeof(len));

            string str(len, ' ');
            ifs.read(&str[0], len);

            seqNames.push_back(str);
        }

        ifs.close();
    }

    // read the first seqID per file (.fsid)
    if (verbose) {
        ss << "Reading " << baseFile << ".fsid...";
        logger.logInfo(ss);
    }
    {
        if (!readArray(baseFile + ".fsid", firstSeqIDPerFile)) {
            throw runtime_error("Cannot open file: " + baseFile +
                                ".fsid\nIs the reference index outdated?");
        }
    }
}

bool IndexInterface::readArray(const string& filename,
                               vector<length_t>& array) {
    // Platform-specific file handling and memory mapping
#ifdef _WIN32
    HANDLE hFile =
        CreateFileA(filename.c_str(), GENERIC_READ, FILE_SHARE_READ, NULL,
                    OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
    if (hFile == INVALID_HANDLE_VALUE) {
        // Handle error
        return false;
    }

    LARGE_INTEGER fileSize;
    if (!GetFileSizeEx(hFile, &fileSize)) {
        CloseHandle(hFile);
        return false;
    }

    if (fileSize.QuadPart % sizeof(length_t) != 0) {
        CloseHandle(hFile);
        return false;
    }

    size_t numElements = fileSize.QuadPart / sizeof(length_t);

    HANDLE hMapFile = CreateFileMapping(hFile, NULL, PAGE_READONLY, 0, 0, NULL);
    if (hMapFile == NULL) {
        CloseHandle(hFile);
        return false;
    }

    length_t* data = static_cast<length_t*>(
        MapViewOfFile(hMapFile, FILE_MAP_READ, 0, 0, fileSize.LowPart));
    if (data == NULL) {
        CloseHandle(hMapFile);
        CloseHandle(hFile);
        return false;
    }

    array.assign(data, data + numElements);

    UnmapViewOfFile(data);
    CloseHandle(hMapFile);
    CloseHandle(hFile);

    return true;
#else // Unix-like systems
    int fd = open(filename.c_str(), O_RDONLY);
    if (fd == -1) {
        // Handle error
        return false;
    }

    struct stat sb;
    if (fstat(fd, &sb) == -1) {
        close(fd);
        return false;
    }

    size_t fileSize = sb.st_size;
    if (fileSize % sizeof(length_t) != 0) {
        close(fd);
        return false;
    }

    size_t numElements = fileSize / sizeof(length_t);

    void* mapped = mmap(NULL, fileSize, PROT_READ, MAP_PRIVATE, fd, 0);
    if (mapped == MAP_FAILED) {
        close(fd);
        return false;
    }

    length_t* data = static_cast<length_t*>(mapped);
    try {
        array.assign(data, data + numElements);
    } catch (...) {
        munmap(mapped, fileSize);
        throw;
    }

    munmap(mapped, fileSize);
    close(fd);

    return true;
#endif
}

void IndexInterface::populateTable(bool verbose) {
    if (wordSize == 0) {
        // table only needs to store full range for empty string
        Kmer::setWordSize(0);
        Kmer kmer("");
        table.insert(make_pair(kmer, getCompleteRange()));
        return;
    }
    if (verbose) {
        stringstream ss;
        ss << "Populating FM-range table with " << wordSize << "-mers...";
        logger.logInfo(ss);
    }

    Kmer::setWordSize(wordSize);

    table.reserve(1 << (2 * wordSize)); // 2 << wordSize is 4^wordSize (DNA)
    setDirection(FORWARD, false);       // table must have bidirectional ranges

    string word;
    vector<FMPosExt> stack;
    Counters counters;
    extendFMPos(getCompleteRange(), stack, counters);
    while (!stack.empty()) {
        auto curr = stack.back();
        stack.pop_back();

        word.resize(curr.getRow());
        word[curr.getRow() - 1] = curr.getCharacter();

        if (curr.getRow() == wordSize) { // max depth reached
            Kmer k(word);
#ifdef RUN_LENGTH_COMPRESSION
            updateRangeSARuns(curr.getRangesMutable());
#endif
            table.insert(make_pair(k, std::move(curr.getRanges())));

        } else // add extra characters
            extendFMPos(curr.getRanges(), stack, counters, curr.getRow());
    }
}

// ----------------------------------------------------------------------------
// ROUTINES FOR APPROXIMATE PATTERN MATCHING
// ----------------------------------------------------------------------------
void IndexInterface::goToInTextVerificationEdit(
    const FMPosExt& node, const Search& s, const vector<Substring>& parts,
    Occurrences& occ, Counters& counters, const Substring& pattern,
    length_t idx, IBitParallelED* bpED, const FMOcc& sMatch,
    const vector<FMPosExt>& dOther, const vector<uint16_t>& iOther) {
#ifdef RUN_LENGTH_COMPRESSION
    // function does nothing in case of run-length compression
    return;
#else

    // we need to find the earliest possible position in the text
    // where the pattern could occur. To do this we need to find the highest
    // startDifference with the position in the current node
    length_t st = parts[s.getLowestPartProcessedBefore(idx)].begin();
    length_t startDiff = st + s.getMaxED();

    if (st == 0) {
        startDiff = 0;
    } else if (dir == BACKWARD) {
        length_t row = node.getRow();
        length_t col = bpED->getFirstColumn(row);
        startDiff -= col + bpED->at(row, col);

    } else if (!dOther.empty()) {
        // FORWARD DIRECTION
        // descOther are thus in backwards direction
        startDiff -= dOther.size() - iOther.size() + iOther.back();
    }

    const auto& pos = getBeginPositions(node.getRanges().getRangeSA(),
                                        startDiff, sMatch.getShift());

    inTextVerification(pos, s.getMaxED(), s.getMinED(), occ, counters, pattern,
                       st == 0);
#endif // RUN_LENGTH_COMPRESSION
}

void IndexInterface::recApproxMatchEdit(
    const Search& s, const FMOcc& startMatch, Occurrences& occ,
    const vector<Substring>& parts, Counters& counters, const int& idx,
    const vector<FMPosExt>& descPrevDir, const vector<uint16_t>& initPrevDir,
    const vector<FMPosExt>& descNotPrevDir,
    const vector<uint16_t>& initNotPrevDir) {
    // shortcut Variables
    const Substring& p = parts[s.getPart(idx)];      // this part
    const length_t& maxED = s.getUpperBound(idx);    // maxED for this part
    const Direction& dir = s.getDirection(idx);      // direction
    const bool& dSwitch = s.getDirectionSwitch(idx); // has direction switched?
    auto& stack = stacks[idx];                       // stack for this partition
    size_t matrixIdx = s.getPart(idx) + (dir == BACKWARD) * s.getNumParts();

    // select the correct bit parallel matrix
    IBitParallelED* bpED;
    if (maxED <= BitParallelED64::getMatrixMaxED()) {
        bpED = &matrices[matrixIdx]; // Assuming matrices is of a type
                                     // compatible with IBitParallelED
    } else {
        bpED = &matrices128[matrixIdx]; // Assuming matrices128 is of a type
                                        // compatible with IBitParallelED
    }

    // get the correct initED and descendants based on switch
    const vector<uint16_t>& initEds = dSwitch ? initNotPrevDir : initPrevDir;
    const vector<FMPosExt>& descendants =
        dSwitch ? descNotPrevDir : descPrevDir;
    const vector<uint16_t>& initOther = dSwitch ? initPrevDir : initNotPrevDir;
    const vector<FMPosExt>& descOther = dSwitch ? descPrevDir : descNotPrevDir;

    // set the direction
    setDirection(dir, s.isUnidirectionalBackwards(idx));

    // calculate necessary increase for first column of bandMatrix
    vector<uint32_t> initED;
    if (initEds.empty()) {
        initED = vector<uint32_t>(1, startMatch.getDistance());
    } else {
        uint16_t prevED =
            (dSwitch ? *min_element(initEds.begin(), initEds.end())
                     : initEds[0]);
        uint32_t increase = startMatch.getDistance() - prevED;
        initED = vector<uint32_t>(initEds.size());
        for (size_t i = 0; i < initED.size(); i++) {
            initED[i] = initEds[i] + increase;
        }
    }

    // encode the sequence of this partition in the matrix if this has not been
    // done before
    if (!bpED->sequenceSet())
        bpED->setSequence(p);

    // initialize bit-parallel matrix
    bpED->initializeMatrix(maxED, initED);

    // initialize  cluster
    Cluster cluster(bpED->getSizeOfFinalColumn(), maxED, startMatch.getDepth(),
                    startMatch.getShift());

    if (bpED->inFinalColumn(0)) {
        // the first row is part of the final column of the banded matrix
        // Update the first cell of the cluster with the startMatch and the
        // value found at the pSize'th column of the initialization row the
        // character does not matter as the first cell of a cluster is never
        // a descendant of any of the other cells in the cluster, so this
        // cell will not be reused for the next part of the pattern
        cluster.setValue(0, FMPosExt((char)0, startMatch.getRanges(), 0),
                         bpED->at(0, p.size()));
    }

    if (!descendants.empty()) {
        // fill in the matrix for the descendants

        length_t maxRow = bpED->getNumberOfRows() - 1;

        for (length_t i = 0;
             i < descendants.size() && descendants[i].getDepth() <= maxRow;
             i++) {

            if (branchAndBound(
                    bpED, cluster, descendants[i], s, idx, parts, occ, counters,
                    initOther, descOther,
                    {descendants.begin() + i + 1, descendants.end()})) {
                return;
            }
        }
        if (descendants.back().getDepth() == maxRow) {
            //  no more rows to possibly check
            return;
        }

        SARangePair pair = descendants.back().getRanges();
        if (dSwitch) {
            // after a switch the ranges of final descendant should be
            // updated
            pair = startMatch.getRanges();
        }

        // push children of final descendant
        extendFMPos(pair, stack, counters, descendants.back().getDepth());

    } else { // get the initial nodes to check
        extendFMPos(startMatch.getRanges(), stack, counters);
    }

#ifndef RUN_LENGTH_COMPRESSION
    // Variables for in-text verification
    bool idxZero = idx == 0;
    Substring pattern(parts.back(), 0, parts.back().end(), FORWARD);
    length_t inTextSwitchPoint = getSwitchPoint();
#endif // not RUN_LENGTH_COMPRESSION

    while (!stack.empty()) {
        const FMPosExt currentNode = stack.back();
        stack.pop_back();

        if (branchAndBound(bpED, cluster, currentNode, s, idx, parts, occ,
                           counters, initOther, descOther)) {

            continue;
        }

#ifndef RUN_LENGTH_COMPRESSION
        if (currentNode.getRanges().width() <= inTextSwitchPoint && !idxZero) {
            // crossing-over
            goToInTextVerificationEdit(currentNode, s, parts, occ, counters,
                                       pattern, idx, bpED, startMatch,
                                       descOther, initOther);
            continue; // continue with next node
        }
#endif // not  RUN_LENGTH_COMPRESSION
        extendFMPos(currentNode, stack, counters);
    }
}

bool IndexInterface::branchAndBound(
    IBitParallelED* bpED, Cluster& cluster, const FMPosExt& currentNode,
    const Search& s, const length_t& idx, const vector<Substring>& parts,
    Occurrences& occ, Counters& counters, const vector<uint16_t>& initOther,
    const vector<FMPosExt>& descOther, const vector<FMPosExt>& remainingDesc) {
    // compute, in a bit-parallel manner, a single row of the ED matrix
    const length_t row = currentNode.getDepth();
    bool validED = bpED->computeRow(row, currentNode.getCharacter());

    // check if we have reached the final column of the matrix
    if (bpED->inFinalColumn(row)) {
        // update the cluster
        length_t clusterIdx = cluster.size() + row - bpED->getNumberOfRows();
        cluster.setValue(clusterIdx, currentNode,
                         bpED->at(row, bpED->getNumberOfCols() - 1));

        if (!validED || bpED->onlyVerticalGapsLeft(row)) {
            // no need to further explore this branch for this part -> go to
            // next part

            goDeeper(cluster, idx + 1, s, parts, occ, counters, descOther,
                     initOther, remainingDesc);
            return true;
        }
    }

    return !validED;
}

void IndexInterface::goDeeper(Cluster& cluster, const length_t& nIdx,
                              const Search& s, const vector<Substring>& parts,
                              Occurrences& occ, Counters& counters,
                              const vector<FMPosExt>& descOtherD,
                              const vector<uint16_t>& initOtherD,
                              const vector<FMPosExt>& remDesc) {
    bool isEdge = s.isEdge(nIdx - 1);
    const auto& lowerBound = s.getLowerBound(nIdx - 1);

    if (isEdge) {
        // if this is final piece report highest minimum (to get shortest
        // match)
        if (nIdx == parts.size()) {
            auto matches = cluster.reportCentersAtEnd();
            for (auto& match : matches) {
                if (match.isValid() && match.getDistance() >= lowerBound) {
                    match.setStrand(strand);
                    match.setPairStatus(pairStatus);
                    occ.addFMOcc(match);
                }
            }
        } else {
            FMOcc match = cluster.reportDeepestMinimum(this->dir);
            if (match.isValid() && match.getDistance() >= lowerBound) {
                // go deeper in search
                Direction originalDir = this->dir;
                recApproxMatchEdit(s, match, occ, parts, counters, nIdx, {}, {},
                                   descOtherD, initOtherD);
                // set direction back again
                setDirection(originalDir,
                             s.isUnidirectionalBackwards(nIdx - 1));
            }
        }

        return;
    }

    // one of the later stages will return to this point, so keep track of
    // the descendants and eds at this branch
    vector<FMPosExt> descendants;
    vector<uint16_t> initEds;

    FMOcc newMatch = cluster.getClusterCentra(lowerBound, descendants, initEds);
    if (!newMatch.isValid()) {
        // no centre above lower bound found
        return;
    }

    // add the remaining descendants (copy)
    descendants.insert(descendants.end(), remDesc.begin(), remDesc.end());

    // reset the depth of all the descendants
    for (length_t i = 0; i < descendants.size(); i++) {
        descendants[i].setDepth(i + 1);
    }
    length_t maxEDNext = s.getUpperBound(nIdx);

    // remove trailing initEds that are higher than maxEDNext
    while (initEds.back() > maxEDNext) {
        initEds.pop_back();
    }

    // is the next direction equal to this direction?
    bool switchAfter = s.getDirectionSwitch(nIdx);

    if (switchAfter) {
        // switching direction as this is not the end of a search direction,
        // this means we'll get back here, thus range of newMatch should be
        // deepest point in branch
        if (!descendants.empty()) {
            newMatch.setRanges(descendants.back().getRanges());

            //  edit distance for search in other direction should be lowest
            //  value possible
            newMatch.setDistance(*min_element(initEds.begin(), initEds.end()));
        }

        Direction originalDir = this->dir;

        recApproxMatchEdit(s, newMatch, occ, parts, counters, nIdx, descendants,
                           initEds, descOtherD, initOtherD);

        // set direction back again
        // s.isUnidirectionalBackwards() is always false here
        setDirection(originalDir, s.isUnidirectionalBackwards(nIdx - 1));
    } else {
        // go deeper on next piece
        recApproxMatchEdit(s, newMatch, occ, parts, counters, nIdx, descendants,
                           initEds, descOtherD, initOtherD);
    }
}

// ----------------------------------------------------------------------------
// EXTEND ROUTINES
// ----------------------------------------------------------------------------

void IndexInterface::extendFMPos(const SARangePair& parentRanges,
                                 vector<FMPosExt>& stack, Counters& counters,
                                 length_t row) const {
    // iterate over the entire alphabet
    for (length_t i = 1; i < sigma.size(); ++i) {

        SARangePair pairForNewChar;

        // check if this character occurs in the specified range
        if ((this->*extraChar)(i, parentRanges, pairForNewChar)) {
            // push this range and character for the next iteration
            stack.emplace_back(sigma.i2c(i), pairForNewChar, row + 1);

            counters.inc(Counters::NODE_COUNTER);
        }
    }
}

void IndexInterface::extendFMPos(const FMPosExt& pos, vector<FMPosExt>& stack,
                                 Counters& counters) const {
    extendFMPos(pos.getRanges(), stack, counters, pos.getDepth());
}

// ----------------------------------------------------------------------------
// ROUTINES FOR ACCESSING THE DATA STRUCTURE
// ----------------------------------------------------------------------------

SeqNameFound IndexInterface::findSeqName(TextOcc& t, length_t& seqID,
                                         Counters& counters,
                                         length_t largestStratum,
                                         const DistanceMetric& metric,
                                         const string& pattern) const {
    // returns FOUND if successful without trimming
    // returns FOUND_WITH_TRIMMING if trimming was needed
    // return NOT_FOUND if no name could be found

    auto& range = t.getRange();
    const auto& begin = range.getBegin();
    const auto& end = range.getEnd();

    // find first start position greater than begin
    auto it = upper_bound(startPos.begin(), startPos.end(), begin);
    // as startPositions[0] == 0 we know that "it" points to at least the second
    // element we need to move one back to find the start position of the
    // sequence in which begin lies
    length_t index = distance(startPos.begin(), --it);

    // now we check if the end position lies in the same sequence
    if (end <= startPos[index + 1]) {
        // adapt the range to be relative to the start of the sequence
        t.setIndexBegin(begin);
        range = Range(begin - startPos[index], end - startPos[index]);
        seqID = index;
        return FOUND;
    }

    // begin and end of occurrence do not lie in the same sequence
    if (metric == HAMMING) {
        // do not allow trimming in case of hamming distance
        return NOT_FOUND;
    }

    // option 1: begin is just before the end of the sequence
    if ((startPos[index + 1] - begin) <= largestStratum) {
        // try the next sequence
        index++;
        // update the beginning of the range and clip end to end of
        // reference
        range = Range(startPos[index], min(end, startPos[index + 1]));
    }
    // option 2: end is just over the start of the next sequence
    else if ((end - startPos[index + 1]) <= largestStratum) {
        // try clipping end
        // update the end of the range
        range = Range(begin, startPos[index + 1]);
    } else {
        return NOT_FOUND;
    }

#ifndef RUN_LENGTH_COMPRESSION
    // do in-text verification with the changed range
    // This guarantees that the CIGAR string (if needed) and edit distance
    // are calculated correctly
    Occurrences occ;

    inTextVerificationOneString(range.getBegin(), range.getEnd(),
                                largestStratum, 0, occ, counters, pattern);

    // TODO update cigar for soft clipping
    if (occ.textOccSize() > 0) {
        // get minimal element using operator<
        const auto& textOccs = occ.getTextOccurrences();
        t = *std::min_element(textOccs.begin(), textOccs.end());
    } else {
        // no match found with clipping
        return NOT_FOUND;
    }
#else
    if (t.getRange().width() - range.width() + t.getDistance() >
        largestStratum) {
        return NOT_FOUND;
    } else {
        // as we cannot verify in text using bmove we assume that trimming X
        // characters increases the score by X characters. In some cases this
        // can be a too pessimistic assumption.
        t = TextOcc(range,
                    t.getDistance() + t.getRange().width() - range.width(),
                    t.getStrand(), t.getPairStatus());
    }
#endif

    seqID = index;
    length_t seqPosition = startPos[index];
    // adapt the range to be relative to the start of the sequence
    t.setIndexBegin(range.getBegin());
    range = Range(range.getBegin() - seqPosition, range.getEnd() - seqPosition);
    return FOUND_WITH_TRIMMING;
}

// ----------------------------------------------------------------------------
// ROUTINES FOR EXACT PATTERN MATCHING
// ----------------------------------------------------------------------------

#ifndef RUN_LENGTH_COMPRESSION

/**
 * @brief Verify the exact match in the text. Helper function for exact matching
 * in the FM Index
 * @param s the string to be matched
 * @param p the positions in the FM Index where the end of the string occurs
 * @param counters the performance counters to be updated
 * @param tOcc the occurrences to be updated
 * @param charsMatched the number of characters matched until now in-index
 * @param CIGAR the CIGAR string to be used
 * @param strand the strand of the read
 * @param pairStatus the pair status of the read
 * @param text the text to be searched
 */
void verifyInTextExact(const std::string& s, const std::vector<length_t>& p,
                       Counters& counters, std::vector<TextOcc>& tOcc,
                       length_t charsMatched, const std::string& CIGAR,
                       const Strand& strand, const PairStatus& pairStatus,
                       const std::string& text) {
    // for each of the positions check if the substring in the text
    // equals the string
    length_t remaining = s.size() - charsMatched; // the remaining characters
                                                  // to be matched
    const std::string prefix =
        s.substr(0, remaining); // the prefix that still needs to be matched

    for (const auto& pos : p) {
        // we need to match the remaining characters BEFORE the position to
        // the prefix as we are matching backwards
        length_t posInText = pos - remaining;

        if (pos >= remaining &&
            Substring(text, posInText, posInText + remaining).equals(prefix)) {
            tOcc.emplace_back(Range(posInText, posInText + s.size()), 0, CIGAR,
                              strand, pairStatus);
        }
    }
    counters.inc(Counters::IN_TEXT_STARTED, p.size());
    counters.inc(Counters::ABORTED_IN_TEXT_VERIF, p.size() - tOcc.size());
}

#endif // end non RLC

void IndexInterface::exactMatchesOutput(const string& s, Counters& counters,
                                        vector<TextOcc>& tOcc) const {
    if (s.size() == 0) {
        return;
    }

    SARangePair completeRanges = getCompleteRange();

#ifndef RUN_LENGTH_COMPRESSION
    SARangeBackwards range = completeRanges.getRangeSA();
#else
    SARangeBackwards range = SARangeBackwards(
        completeRanges.getRangeSA(), completeRanges.getToehold(),
        completeRanges.getToeholdRepresentsEnd(),
        completeRanges.getOriginalDepth());
#endif // end non RLC

    length_t i = s.size();

    for (; i-- > 0;) {
        const auto& c = s[i];
        auto pos = sigma.c2i(c);
        if (pos == -1 || !findRangeWithExtraCharBackward(pos, range, range)) {
            // c is not in alphabet or no updated range is empty, exact match
            // impossible
            return;
        }
        counters.inc(Counters::NODE_COUNTER);

#ifndef RUN_LENGTH_COMPRESSION
        // check for in-text verification
        if (range.width() <= getSwitchPoint()) {
            break;
        }
#endif // end non RLC
    }

    // get all positions from the range
    const auto& p = getBeginPositions(range, 0);

#ifndef RUN_LENGTH_COMPRESSION
    // exact match, so all M characters
    std::string CIGAR = getNoCIGAR() ? "*" : fmt::format("{}M", s.size());
    if (i != (length_t)-1) {
        verifyInTextExact(s, p, counters, tOcc, s.size() - i, CIGAR, strand,
                          pairStatus, getText());
    } else {
        // No verification needed, we can just add each occurrence
        for (const auto& pos : p) {
            tOcc.emplace_back(Range(pos, pos + s.size()), 0, CIGAR, strand,
                              pairStatus);
        }
    }
#else
    std::string CIGAR = "*"; // no CIGAR in RLC
    // directly add occurrences
    for (const auto& pos : p) {
        tOcc.emplace_back(Range(pos, pos + s.size()), 0, CIGAR, strand,
                          pairStatus);
    }
#endif // end RLC

    counters.inc(Counters::TOTAL_REPORTED_POSITIONS, tOcc.size());
}

SARangePair
IndexInterface::matchStringBidirectionally(const Substring& pattern,
                                           SARangePair rangesOfPrev,
                                           Counters& counters) const {
    assert(pattern.getDirection() == dir);

    for (length_t i = 0; i < pattern.size(); i++) {

        char c = pattern[i];
        if (!addChar(c, rangesOfPrev, counters)) {
            // rangesOfPrev was made empty
            break;
        }
    }

    return rangesOfPrev;
}

bool IndexInterface::addChar(const char& c, SARangePair& startRange,
                             Counters& counters) const {

    int posInAlphabet = sigma.c2i((unsigned char)c);
    if (posInAlphabet > -1) {

        if ((this->*extraChar)(posInAlphabet, startRange, startRange)) {
            // each character that we look at is a new node that is visited
            counters.inc(Counters::NODE_COUNTER);
            return true;
        }
    }
    // the range is now empty or c does not exist in alphabet
    startRange = SARangePair();
    return false;
}

// ----------------------------------------------------------------------------
// ROUTINES FOR APPROXIMATE MATCHING
// ----------------------------------------------------------------------------

void IndexInterface::approxMatchesNaive(const std::string& pattern,
                                        length_t maxED, Counters& counters,
                                        vector<TextOcc>& matches) {

    Occurrences occurrences;

    resetMatrices(1); // reset one matrix (both 64 and 128 bit)
    IBitParallelED* matrix;
    if (maxED > BitParallelED64::getMatrixMaxED()) {
        matrix = &matrices128.front();
    } else {
        matrix = &matrices.front();
    }
    matrix->setSequence(pattern);
    matrix->initializeMatrix(maxED);

    setDirection(FORWARD, false);

    vector<FMPosExt> stack;
    stack.reserve((pattern.size() + maxED + 1) * (sigma.size() - 1));

    extendFMPos(getCompleteRange(), stack, counters, 0);

    length_t lastCol = pattern.size();

    while (!stack.empty()) {
        const FMPosExt currentNode = stack.back();
        stack.pop_back();
        length_t row = currentNode.getDepth();

        if (row >= matrix->getNumberOfRows()) {
            continue;
        }

        bool valid = matrix->computeRow(row, currentNode.getCharacter());

        if (!valid) {
            continue; // backtrack
        }

        if (matrix->inFinalColumn(row)) {
            // full pattern was matched
            if (matrix->at(row, lastCol) <= maxED) {
                occurrences.addFMOcc(currentNode, matrix->at(row, lastCol),
                                     strand, pairStatus);
            }
        }

#ifndef RUN_LENGTH_COMPRESSION // check for in-text verification
        if (currentNode.getRanges().width() <= getSwitchPoint()) {
            // crossing-over to in-text verification
            // convert node to an occurrence in the text with the
            // correct edit distance
            const auto& startPos =
                getBeginPositions(currentNode.getRanges().getRangeSA(), 0, 0);

            inTextVerification(startPos, maxED, 0, occurrences, counters,
                               {pattern, FORWARD},
                               true); // naive search is forward only
            continue;
        }
#endif // end non RLC
        extendFMPos(currentNode, stack, counters);
    }

    const auto& bundle = ReadBundle::createBundle(pattern, strand);
    matches = getUniqueTextOccurrences(occurrences, maxED, counters, bundle);

#ifndef RUN_LENGTH_COMPRESSION // CIGAR generation
    generateCIGARS(matches, counters, bundle);
#endif // end non RLC
}

void IndexInterface::recApproxMatchHamming(const Search& s,
                                           const FMOcc& startMatch,
                                           Occurrences& occ,
                                           const vector<Substring>& parts,
                                           Counters& counters, const int& idx) {
    // shortcut variables
    const Substring& p = parts[s.getPart(idx)]; // the current part
    const length_t& pSize = p.size();           // the size of the current part
    const Direction& d = s.getDirection(idx);   // direction of current part
    const length_t& maxED = s.getUpperBound(idx); // upper bound of current part
    const length_t& minED = s.getLowerBound(idx); // lower bound of current part
    setDirection(d, s.isUnidirectionalBackwards(idx));

    // create vector for the scores
    vector<length_t> vec(p.size() + 1, 0);
    // set root element of vector to the distance of the start match
    vec[0] = startMatch.getDistance();
    // get stack for current part
    auto& stack = stacks[idx];

    extendFMPos(startMatch.getRanges(), stack, counters);

    while (!stack.empty()) {
        const FMPosExt node = stack.back();
        stack.pop_back();

#ifndef RUN_LENGTH_COMPRESSION // in-text verification
        if (node.getRanges().width() <= getSwitchPoint()) {
            inTextVerificationHamming(node, s, parts, idx, occ, counters);
            continue;
        }
#endif // end non RLC

        // update the vector
        length_t row = node.getRow();
        vec[row] = vec[row - 1] + (node.getCharacter() != p[row - 1]);

        if (vec[row] > maxED) {
            // backtrack
            continue;
        }

        if (row == pSize) {
            // end of part
            if (vec[row] >= minED) {
                // valid occurrence
                FMOcc match =
                    FMOcc(node.getRanges(), vec[row],
                          startMatch.getDepth() + pSize, strand, pairStatus);
                if (s.isEnd(idx)) {
                    // end of search
                    occ.addFMOcc(match);
                } else {
                    // continue search
                    recApproxMatchHamming(s, match, occ, parts, counters,
                                          idx + 1);
                    setDirection(s.getDirection(idx),
                                 s.isUnidirectionalBackwards(idx));
                }
            }
            continue;
        }
        extendFMPos(node, stack, counters);
    }
}

void IndexInterface::recApproxMatchEditEntry(
    const Search& search, const FMOcc& startMatch, Occurrences& occ,
    const vector<Substring>& parts, Counters& counters, const int& idx) {

#ifndef RUN_LENGTH_COMPRESSION
    if (startMatch.getRanges().width() > getSwitchPoint()) {
        counters.inc(Counters::SEARCH_STARTED);
        recApproxMatchEdit(search, startMatch, occ, parts, counters, idx);
        return;
    }
    // verify the partial match in text
    verifyExactPartialMatchInText(
        startMatch, parts[search.getLowestPartProcessedBefore(idx)].begin(),
        search.getMaxED(), occ, counters, search.getMinED(),
        {parts.back(), 0, parts.back().end(), FORWARD});
#else
    counters.inc(Counters::SEARCH_STARTED);
    recApproxMatchEdit(search, startMatch, occ, parts, counters, idx);
#endif
}

// ----------------------------------------------------------------------------
// POST-PROCESSING ROUTINES FOR APPROXIMATE PATTERN MATCHING
// ----------------------------------------------------------------------------

vector<TextOcc> IndexInterface::getTextOccHamming(Occurrences& occ,
                                                  Counters& counters) const {
    counters.inc(Counters::TOTAL_REPORTED_POSITIONS, occ.textOccSize());

    // erase in-index doubles
    occ.eraseDoublesFM();

    // Get the in-index occurrences
    const auto& fmOccs = occ.getFMOccurrences();

    // The size of the match
    length_t size = (fmOccs.empty()) ? 0 : fmOccs[0].getDepth();

    std::string CIGAR = "*"; // default CIGAR string

#ifndef RUN_LENGTH_COMPRESSION
    // Calculate the CIGAR string for hamming distance
    CIGAR = (getNoCIGAR()) ? "*" : fmt::format("{}M", size);
#endif

    for (const auto& fmOcc : fmOccs) {
        // Get the range
        const SARange& saRange = fmOcc.getRanges().getRangeSA();
        counters.inc(Counters::TOTAL_REPORTED_POSITIONS, saRange.width());

        auto fmStrand = fmOcc.getStrand();
        auto fmPairStatus = fmOcc.getPairStatus();

        // find the positions in the text for this index-occurrence
        std::vector<length_t> textPositions;
        getTextPositionsFromSARange(fmOcc.getRanges(), textPositions);

        // add the text occurrences
        for (length_t p : textPositions) {
            occ.addTextOcc(Range(p, p + size), fmOcc.getDistance(), CIGAR,
                           fmStrand, fmPairStatus);
        }
    }
    // remove doubles
    occ.eraseDoublesText();

    return occ.getTextOccurrences();
}

vector<TextOcc> IndexInterface::getUniqueTextOccurrences(
    Occurrences& occ, const length_t& maxED, Counters& counters,
    const ReadBundle& bundle) const {

    // TODO: add sorted requirement to doc
    // TODO: add assertion at the end that checks if it was sorted

    // increment report position counter
    counters.inc(Counters::TOTAL_REPORTED_POSITIONS, occ.textOccSize());

    // erase equal occurrences from the in-index occurrences
    occ.eraseDoublesFM();

    // convert the in-index occurrences to in-text occurrences
    const auto& fmOccs = occ.getFMOccurrences();
    for (const auto& f : fmOccs) {

        const SARange& saRange = f.getRanges().getRangeSA();

        // increment reported positions counter
        counters.inc(Counters::TOTAL_REPORTED_POSITIONS, saRange.width());

        // whether to calculate the CIGAR string. If the CIGAR_THRESHOLD is
        // lower than or equal to the in-text switch point, the CIGAR string is
        // calculated for all in-index occurrences. Otherwise, the CIGAR string
        // is calculated only for those in-index occurrences that occur more
        // than CIGAR_THRESHOLD times in the text.
        const bool calculateCIGAR =
#ifndef RUN_LENGTH_COMPRESSION
            (!getNoCIGAR()) && (saRange.width() > CIGAR_THRESHOLD);
#else
            false;
#endif

        std::string CIGARstring;
        auto depth = f.getDepth(), distance = f.getDistance();
        length_t shift = f.getShift();

        // find the positions in the text for this index-occurrence
        std::vector<length_t> textPositions;
        getTextPositionsFromSARange(f.getRanges(), textPositions);

        bool first = true;
        // add the text occurrences
        for (length_t p : textPositions) {
            length_t startPos = p + shift;
            assert(startPos < textLength);

            // convert to a text occurrence
            TextOcc tOcc = TextOcc(Range(startPos, startPos + depth), distance,
                                   f.getStrand(), f.getPairStatus());

            if (calculateCIGAR) {
                // this will always be false for RUN_LENGTH_COMPRESSION
                if (first) {
                    // calculate the CIGAR string for the first encountered
                    // in-text occurrence of this FM occurrence. The CIGAR
                    // string is equal for all text occurrences within this FM
                    // occurrence
                    const auto& seq = bundle.getSequence(tOcc.getStrand());
                    generateCIGAR(tOcc, counters, seq);
                    CIGARstring = tOcc.getCigar();
                    first = false;
                } else {
                    tOcc.setCigar(CIGARstring);
                }
            }
            occ.addTextOcc(tOcc);
        }
    }

    // erase equal occurrences from the in-text occurrences, note
    // that an in-text occurrence with calculated CIGAR string takes
    // preference over an equal one without CIGAR string
    occ.eraseDoublesText();

    // find the non-redundant occurrences
    vector<TextOcc> nonRedundantOcc;
    nonRedundantOcc.reserve(occ.textOccSize());

    length_t maxDiff = 2 * maxED;
    length_t prevBegin = numeric_limits<length_t>::max();
    length_t prevDepth = numeric_limits<length_t>::max();
    length_t prevED = maxED + 1;

    const auto& textOcc = occ.getTextOccurrences();
    for (const auto& o : textOcc) {
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
            if (o.getDistance() > prevED ||
                (o.getDistance() == prevED &&
                 o.getRange().width() >= prevDepth)) {
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

    // assert that nonRedundantOcc is sorted
    assert(is_sorted(nonRedundantOcc.begin(), nonRedundantOcc.end()));

    return nonRedundantOcc;
}
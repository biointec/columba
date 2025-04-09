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

#include "indexhelpers.h"
#ifndef RUN_LENGTH_COMPRESSION
#include "bitparallelmatrix.h" // for BitParallelED
#endif
#include "logger.h"     // for logger
#include <cstdint>      // for uint16_t, uint64_t
#include <fmt/core.h>   // for format
#include <fmt/format.h> // for to_string

using namespace std;

// ============================================================================
// CLASS RANGE
// ============================================================================

ostream& operator<<(ostream& os, const Range& r) {
    os << "[" << r.begin << ", " << r.end << ")";
    return os;
}

#ifdef RUN_LENGTH_COMPRESSION

// ============================================================================
// CLASS MOVERANGE
// ============================================================================

ostream& operator<<(ostream& os, const MoveRange& r) {
    os << "[" << r.begin << ", " << r.end << ") in runs [" << r.beginRun << ", "
       << r.endRun << "]";
    return os;
}

#endif

// ============================================================================
// CLASS TEXT OCC
// ============================================================================
void TextOcc::generateSAMSingleEnd(const string& seqID, const string& printSeq,
                                   const string& printQual, length_t nHits,
                                   length_t minScore, bool primaryAlignment,
                                   const vector<string>& seqNames) {

    // Estimate the length of the resulting string and reserve memory
    size_t estimated_length =
        seqID.length() + printSeq.length() + printQual.length() + 100;

    outputLine.reserve(estimated_length);

    // Precompute constant values to avoid repeated calculations
    uint16_t flags = getFlagsSE(primaryAlignment);
    int16_t mapQ = getMapQ(nHits, minScore);
    length_t pos = range.getBegin() + 1; // SAM is 1-based

    // Format the output line
    outputLine =
        fmt::format("{}\t{}\t{}\t{}\t{}\t{}\t*\t0\t0\t{}\t{}\tAS:i:{}\tNM:i:{}"
                    "\tPG:Z:Columba\n",
                    seqID,                        // read name
                    flags,                        // sam flags
                    seqNames[assignedSequenceID], // reference sequence name
                    pos,         // 1-based pos in reference sequence
                    mapQ,        // mapping quality
                    stringCIGAR, // CIGAR string
                    printSeq,    // sequence or *
                    printQual,   // quality or *
                    distance,    // AS:i: distance
                    distance     // NM:i: distance
        );
}

void TextOcc::generateSAMSingleEndXA(
    const std::string& seqID, const std::string& printSeq,
    const std::string& printQual, length_t nHits,
    std::vector<TextOcc>::const_iterator otherMatchesBegin,
    std::vector<TextOcc>::const_iterator otherMatchesEnd,
    const std::vector<std::string>& seqNames) {

    generateSAMSingleEnd(seqID, printSeq, printQual, nHits, distance, true,
                         seqNames);
    assert(outputLine.back() == '\n'); // generateSAMSingleEnd should add \n
    // remove \n from end of the line
    outputLine.pop_back();
    // add the X0 X1 and XA tag
    length_t x0 = nHits - 1; // number of co-optimal hits (nHits includes this)
    // number of suboptimal hits
    length_t x1 = std::distance(otherMatchesBegin, otherMatchesEnd) - x0;

    // Append the SAM line
    outputLine += fmt::format("\tX0:i:{}\tX1:i:{}\tXA:Z:", x0, x1);
    for (auto it = otherMatchesBegin; it != otherMatchesEnd; ++it) {
        outputLine += it->asXA(seqNames); // Use the iterator to access elements
    }
    outputLine += "\n";
}

void TextOcc::generateSAMPairedEnd(ReadBundle& bundle, uint32_t nPairs,
                                   uint32_t minScore, const TextOcc& mateOcc,
                                   const uint32_t fragSize, bool discordant,
                                   bool primaryAlignment,
                                   const vector<string>& seqNames) {
    // Precompute constant values to avoid repeated calculations
    const string& printSeq =
        (isRevCompl()) ? bundle.getRevComp() : bundle.getRead();
    string printQual =
        (isRevCompl()) ? bundle.getRevQuality() : bundle.getQual();
    if (printQual.empty()) {
        printQual = "*";
    }

    uint16_t flags = getFlagsPE(mateOcc, discordant, primaryAlignment);
    length_t pos = range.getBegin() + 1; // SAM is 1-based

    size_t estimated_length = bundle.getSeqID().length() + printSeq.length() +
                              printQual.length() + 150;
    outputLine.reserve(estimated_length);

    bool mateMapped = mateOcc.isValid();
    const auto& assignedSeq = seqNames[assignedSequenceID];
    const auto& mateAssignedSeq =
        (mateMapped) ? seqNames[mateOcc.assignedSequenceID] : "*";
    const auto& matePos = (mateMapped) ? mateOcc.range.getBegin() + 1 : 0;
    const auto& sign =
        (mateMapped && range.getBegin() > mateOcc.range.getBegin()) ? "-" : "";
    const auto insertSize = (mateMapped) ? fragSize : 0;

    outputLine = fmt::format(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}{}\t{}\t{}\tAS:i:{}\tNM:i:"
        "{}\tPG:Z:Columba\n",
        bundle.getSeqID(), // read name
        flags,             // SAM flags
        assignedSeq,       // reference sequence name
        pos,               // 1-based pos in reference sequence
        getMapQPairedEnd(nPairs, minScore,
                         mateOcc.getDistance()), // mapping quality
        stringCIGAR,                             // CIGAR string
        mateAssignedSeq,                         // mate ref seq name
        matePos,                                 // mate pos (1-based)
        sign, insertSize,                        // template length
        printSeq,                                // read sequence
        printQual,                               // read quality
        distance,                                // alignment score of this read
        distance // distance to reference of this read
    );
}

TextOcc TextOcc::createUnmappedSAMOccurrenceSE(const ReadBundle& bundle) {
    TextOcc t;
    // Estimate the final size and reserve space to avoid multiple
    // allocations
    size_t estimatedSize = bundle.getSeqID().length() +
                           bundle.getRead().length() +
                           bundle.getQual().length() +
                           50; // 50 is a rough estimate for the fixed parts

    t.outputLine.reserve(estimatedSize);

    // Construct the SAM line directly
    t.outputLine =
        fmt::format("{}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}\tPG:Z:Columba\n",
                    bundle.getSeqID(), // Read name
                    bundle.getRead(),  // Read sequence
                    bundle.getQual()   // Read quality
        );

    return t;
}

TextOcc TextOcc::createUnmappedSAMOccurrencePE(const ReadBundle& bundle,
                                               PairStatus pairStatus,
                                               bool mateMapped,
                                               Strand mateStrand) {
    auto t = TextOcc();
    t.setDistance(0); // an unmapped occurrence has distance 0
    t.setPairStatus(pairStatus);
    t.strand = FORWARD_STRAND;

    uint16_t flags = 1;                                 // paired flag
    flags |= 4;                                         // unmapped flag
    flags |= (mateMapped ? 0 : 8);                      // mate unmapped flag
    flags |= (mateStrand == REVERSE_C_STRAND ? 32 : 0); // mate rev comp flag
    flags |= (pairStatus == FIRST_IN_PAIR ? 64 : 128);  // first in pair flag

    t.outputLine.clear();
    t.outputLine.reserve(bundle.getSeqID().length() +
                         bundle.getRead().length() + bundle.getQual().length() +
                         100); // 100 is a rough estimate for the fixed parts
    t.outputLine =
        fmt::format("{}\t{}\t*\t0\t0\t*\t*\t0\t0\t{}\t{}\tPG:Z:Columba\n",
                    bundle.getSeqID(),     // read name
                    fmt::to_string(flags), // flags (unmapped)
                    bundle.getRead(),      // read sequence
                    bundle.getQual()       // read quality
        );

    return t;
}

void TextOcc::generateSAMUnpaired(ReadBundle& bundle, uint32_t nHits,
                                  uint32_t minScore, bool primaryAlignment,
                                  const vector<string>& seqNames) {
    // generate the SAM line for this occurrence
    // this occurrence was part of a pair, its mate is mapped, but no pair
    // could be found (if discordant not allowed) or there are too many
    // discordant pairs to handle

    // TODO: remove redundant and merge singleEnd and unpaired function
    const string& printSeq =
        (primaryAlignment)
            ? ((isRevCompl()) ? bundle.getRevComp() : bundle.getRead())
            : "*";
    string printQual =
        (primaryAlignment)
            ? ((isRevCompl()) ? bundle.getRevQuality() : bundle.getQual())
            : "*";

    if (printQual.empty()) {
        printQual = "*";
    }

    uint16_t flags = 1;                                // paired flag
    flags += (pairStatus == FIRST_IN_PAIR ? 64 : 128); // first in pair flag
    flags |= (!primaryAlignment << 8);                 // set secondary flag

    outputLine.clear();
    outputLine.reserve(bundle.getSeqID().length() + bundle.getRead().length() +
                       bundle.getQual().length() +
                       150); // 150 is a rough estimate for the fixed parts
    outputLine =
        fmt::format("{}\t{}\t{}\t{}\t{}\t{}\t*\t0\t0\t{}\t{}\tAS:i:{}\tNM:i:{}"
                    "\tPG:Z:Columba\n",
                    bundle.getSeqID(),            // read name
                    flags,                        // flags
                    seqNames[assignedSequenceID], // reference sequence name
                    range.getBegin() + 1,         // 1-based pos in ref seq
                    getMapQ(nHits, minScore),     // mapping quality
                    stringCIGAR,                  // CIGAR string
                    printSeq,                     // read sequence
                    printQual,                    // read quality
                    distance,                     // alignment score
                    distance                      // distance to reference
        );
}

// ============================================================================
// CLASS FM OCC
// ============================================================================

ostream& operator<<(ostream& o, const FMOcc& m) {
    return o << "SARange: " << m.getRanges().getRangeSA()
             << "\tEdit distance: " << m.getDistance()
             << "\tdepth: " << m.getDepth();
}

// ============================================================================
// CLASS CLUSTER
// ============================================================================

FMOcc Cluster::getClusterCentra(uint16_t lowerBound, vector<FMPosExt>& desc,
                                vector<uint16_t>& initEds) {
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
                if (initEds[k] < lowerBound && initEds[k] <= initEds[k - 1] &&
                    (k == initEds.size() - 1 || initEds[k] <= initEds[k + 1])) {
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
                                initEds[hC] =
                                    min(maxED + 1, initEds[hC - 1] + 1);
                                hC++;
                            } else {
                                initEds[lC] =
                                    min(maxED + 1, initEds[lC + 1] + 1);
                                lC--;
                            }
                            highest = !highest;
                        }
                        if (lC == hC) {
                            // change middle element of cluster
                            initEds[lC] =
                                min(initEds[lC + 1] + 1, initEds[lC - 1] + 1);
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

// ============================================================================
// CLASS COUNTERS
// ============================================================================

void Counters::reportStatistics(const SequencingMode& sMode) const {
    std::stringstream ss;

    bool zeroReads = (counters[NUMBER_OF_READS] == 0);
    if (!zeroReads) {
        ss << "Average no. nodes: "
           << counters[NODE_COUNTER] / (counters[NUMBER_OF_READS] * 1.0);
        logger.logVerbose(ss);
    }

    ss << "Total no. Nodes: " << counters[NODE_COUNTER];
    logger.logVerbose(ss);

    if (sMode == SINGLE_END) {
        if (!zeroReads) {
            ss << "Average no. unique matches per read: "
               << counters[TOTAL_UNIQUE_MATCHES] /
                      (counters[NUMBER_OF_READS] * 1.0);
            logger.logInfo(ss);
        }

        ss << "Total no. matches: " << counters[TOTAL_UNIQUE_MATCHES];
        logger.logInfo(ss);

        if (!zeroReads) {
            ss << "Average no. matches per read "
               << counters[TOTAL_REPORTED_POSITIONS] /
                      (counters[NUMBER_OF_READS] * 1.0);
            logger.logVerbose(ss);
        }

        ss << "Total no. reported matches: "
           << counters[TOTAL_REPORTED_POSITIONS];
        logger.logVerbose(ss);

        ss << "Mapped reads: " << counters[MAPPED_READS];
        logger.logInfo(ss);

        ss << "Number of reads: " << counters[NUMBER_OF_READS];
        logger.logInfo(ss);

        if (!zeroReads) {
            ss << "Percentage reads mapped: "
               << (counters[MAPPED_READS] * 100.0) / counters[NUMBER_OF_READS]
               << "%";
            logger.logInfo(ss);
        }
    } else {
        // Paired end
        if (!zeroReads) {
            ss << "Average no. matches per pair: "
               << counters[TOTAL_UNIQUE_PAIRS] /
                      (counters[NUMBER_OF_READS] / 2.0);
            logger.logInfo(ss);
        }

        ss << "Total no. matches : " << counters[TOTAL_UNIQUE_PAIRS];
        logger.logInfo(ss);

        ss << "Mapped pairs: " << counters[MAPPED_PAIRS];
        logger.logInfo(ss);

        if (!zeroReads) {
            ss << "Percentage of pairs mapped: "
               << (counters[MAPPED_PAIRS] * 100.0) /
                      (counters[NUMBER_OF_READS] / 2)
               << "%";
            logger.logInfo(ss);
        }

        ss << "Discordantly mapped pairs: "
           << counters[DISCORDANTLY_MAPPED_PAIRS];
        logger.logInfo(ss);

        if (!zeroReads) {
            ss << "Percentage of discordantly mapped pairs: "
               << (counters[DISCORDANTLY_MAPPED_PAIRS] * 100.0) /
                      (counters[NUMBER_OF_READS] / 2)
               << "%";
            logger.logInfo(ss);
        }

        ss << "No. unpaired reads that did match: "
           << counters[MAPPED_HALF_PAIRS];
        logger.logInfo(ss);

        if (!zeroReads) {
            ss << "Percentage pairs for which only 1 read matched: "
               << (counters[MAPPED_HALF_PAIRS] * 100.0) /
                      (counters[NUMBER_OF_READS] / 2)
               << "%";
            logger.logInfo(ss);
        }

        ss << "Total read pairs both mapped but unpaired: "
           << counters[UNPAIRED_BUT_MAPPED_PAIRS];
        logger.logInfo(ss);
    }

// reports on in-text verification, not needed in RLC
#ifndef RUN_LENGTH_COMPRESSION
    ss << "In text verification procedures " << counters[IN_TEXT_STARTED];
    logger.logVerbose(ss);

    ss << "Failed in-text verifications procedures: "
       << counters[ABORTED_IN_TEXT_VERIF];
    logger.logVerbose(ss);

    if (counters[IN_TEXT_STARTED] != 0) {
        ss << "Aborted in-text relative to started "
           << (counters[ABORTED_IN_TEXT_VERIF] * 1.0) /
                  counters[IN_TEXT_STARTED];
    } else {
        ss << "Aborted in-text relative to started: N/A (No in-text "
              "verifications started)";
    }
    logger.logVerbose(ss);

    ss << "Immediate switch after first part: " << counters[IMMEDIATE_SWITCH];
    logger.logVerbose(ss);

    ss << "Searches started (does not include immediate switches) : "
       << counters[SEARCH_STARTED];
    logger.logVerbose(ss);
#endif // end not RLC
}

#ifndef RUN_LENGTH_COMPRESSION // No in-text verification in RLC

template <typename WordType>
void InTextVerificationTask<WordType>::doTask(Counters& counters,
                                              Occurrences& occ) {
    assert(matrix->sequenceSet());
    counters.inc(Counters::IN_TEXT_STARTED, refs.size());

    for (const auto& ref : refs) {

        length_t refBegin = ref.begin();

        length_t i;
        const length_t size = ref.size();

        if (!matrix->inFinalColumn(size)) {
            // no possible match
            continue;
        }

        for (i = 0; i < size; ++i) {
            if (!matrix->computeRow(i + 1, ref.forwardAccessor(i))) {
                break;
            }
        }

        // did we break before a possible match?
        if (i <= size - matrix->getSizeOfFinalColumn()) {
            counters.inc(Counters::ABORTED_IN_TEXT_VERIF);
            continue;
        }

        vector<length_t> refEnds;
        matrix->findClusterCenters(i, refEnds, maxED, minED);

        if (refEnds.empty()) {
            counters.inc(Counters::ABORTED_IN_TEXT_VERIF);
            continue;
        }

        // for each valid end -> calculate CIGAR string and report
        for (const auto& refEnd : refEnds) {
            length_t bestScore = maxED + 1, bestBegin = 0;

            // we need to track back to find out from where (what position in
            // the text) this occurrence starts
            std::string cigar;
            matrix->traceBack(ref, refEnd, bestBegin, bestScore, cigar);
            counters.inc(Counters::CIGARS_IN_TEXT_VERIFICATION);

            if (noCIGAR) {
                cigar = "";
            }

            // make an occurrence
            occ.addTextOcc(Range(refBegin + bestBegin, refBegin + refEnd),
                           bestScore, std::move(cigar), strand, pairStatus);
        }
    }
}

template class InTextVerificationTask<uint64_t>;
template class InTextVerificationTask<UInt128>;

#endif // end not RLC

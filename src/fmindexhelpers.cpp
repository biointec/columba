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

#include "fmindexhelpers.h"
#include <limits>

using namespace std;

ostream& operator<<(ostream& os, const Range& r) {
    os << "[" << r.begin << ", " << r.end << ")";
    return os;
}

// ============================================================================
// CLASS BIFMOCC
// ============================================================================

ostream& operator<<(ostream& o, const FMOcc& m) {
    return o << "SARange: " << m.getRanges().getRangeSA()
             << "\tEdit distance: " << m.getDistance()
             << "\tdepth: " << m.getDepth();
}

// ============================================================================
// CLASS CLUSTER
// ============================================================================

FMOcc Cluster::getClusterCentra(uint lowerBound, std::vector<FMPosExt>& desc,
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
                                    std::min(maxED + 1, initEds[hC - 1] + 1);
                                hC++;
                            } else {
                                initEds[lC] =
                                    std::min(maxED + 1, initEds[lC + 1] + 1);
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

void IntextVerificationTask::doTask(Counters& counters, Occurrences& occ) {
    counters.inTextStarted += refs.size();
    for (const auto& ref : refs) {

        length_t i;
        length_t size = ref.size();

        for (i = 0; i < size; i++) {
            if (!matrix.computeRow(i + 1, ref[i])) {
                break;
            }
        }

        // did we break before a possible match?
        if (i <= size - matrix.getSizeOfFinalColumn()) {
            counters.abortedInTextVerificationCounter++;
            continue;
        }

        std::vector<length_t> refEnds;
        matrix.findClusterCenters(i, refEnds, maxED, minED);

        if (refEnds.empty()) {
            counters.abortedInTextVerificationCounter++;
            continue;
        }

        // for each valid end -> calculate CIGAR string and report
        for (const auto& refEnd : refEnds) {
            length_t bestScore = maxED + 1, bestBegin = 0;

            std::vector<std::pair<char, uint>> CIGAR;
            matrix.trackBack(ref, refEnd, bestBegin, bestScore, CIGAR);
            counters.cigarsInTextVerification++;

            // make an occurrence
            occ.addTextOcc(Range(ref.begin() + bestBegin, ref.begin() + refEnd),
                           bestScore, CIGAR);
        }
    }
}

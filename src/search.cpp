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

#include "search.h"

using namespace std;

ostream& operator<<(ostream& os, const Search& obj) {
    os << "{";
    length_t numParts = obj.getNumParts();
    for (length_t i = 0; i < numParts; i++) {
        os << obj.getPart(i) << ((i == numParts - 1) ? "}" : ",");
    }
    os << " {";
    for (length_t i = 0; i < numParts; i++) {
        os << obj.getLowerBound(i) << ((i == numParts - 1) ? "}" : ",");
    }
    os << " {";
    for (length_t i = 0; i < numParts; i++) {
        os << obj.getUpperBound(i) << ((i == numParts - 1) ? "}" : ",");
    }
    return os;
}
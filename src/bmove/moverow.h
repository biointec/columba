/******************************************************************************
 *  Columba: Approximate Pattern Matching using Search Schemes                *
 *  Copyright (C) 2020-2024 - Lore Depuydt <lore.depuydt@ugent.be> and        *
 *                            Luca Renders <luca.renders@ugent.be> and        *
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
#ifndef MOVEROW_H
#define MOVEROW_H

#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include "../definitions.h"

#include <sdsl/int_vector.hpp>

using namespace std;

struct MoveRow {

    // Run head for this move row. 'c' in the b-move paper. Encoded.
    uint8_t c;

    // The position index of the start of the input interval in the BWT. 'p' in
    // the b-move paper.
    length_t inputStartPos;

    // The position index of the mapping of the first element of the input
    // interval (first element of ouput interval). '\pi' in the b-move paper.
    length_t outputStartPos;

    // Index of the input interval containing the mapping. '\ksi' in the b-move
    // paper.
    length_t outputStartRun;

    MoveRow() {
    }

    MoveRow(uint8_t c, length_t inputStartPos, length_t outputStartPos,
            length_t outputStartRun)
        : c(c), inputStartPos(inputStartPos), outputStartPos(outputStartPos),
          outputStartRun(outputStartRun) {
    }

    MoveRow(ifstream& ifs) {
        ifs.read((char*)&c, sizeof(c));
        ifs.read((char*)&inputStartPos, sizeof(inputStartPos));
        ifs.read((char*)&outputStartPos, sizeof(outputStartPos));
        ifs.read((char*)&outputStartRun, sizeof(outputStartRun));
    }

    length_t serialize(ofstream& out) const {
        out.write((char*)&c, sizeof(c));
        out.write((char*)&inputStartPos, sizeof(inputStartPos));
        out.write((char*)&outputStartPos, sizeof(outputStartPos));
        out.write((char*)&outputStartRun, sizeof(outputStartRun));

        return sizeof(c) + sizeof(inputStartPos) + sizeof(outputStartPos) +
               sizeof(outputStartRun);
    }
};

#endif
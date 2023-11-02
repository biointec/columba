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

#ifndef WORDLENGTH_H
#define WORDLENGTH_H

#include <cstdint>

#if THIRTY_TWO
#define LENGTH_TYPE_NAME "32-bits"
typedef uint32_t length_t;
#else
#define LENGTH_TYPE_NAME "64-bits"
typedef uint64_t length_t;
#endif

#endif

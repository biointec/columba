/******************************************************************************
 *   Copyright (C) 2014 - 2021 Jan Fostier (jan.fostier@ugent.be)             *
 *   This file is part of Detox                                               *
 *                                                                            *
 *   This program is free software; you can redistribute it and/or modify     *
 *   it under the terms of the GNU Affero General Public License as published *
 *   by the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                      *
 *                                                                            *
 *   This program is distributed in the hope that it will be useful,          *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                             *
 *                                                                            *
 *   You should have received a copy of the GNU Affero General Public License *
 *   along with this program; if not, see <https://www.gnu.org/licenses/>.    *
 ******************************************************************************/

#include "nucleotide.h"
#include <iosfwd> // for std

using namespace std;

// ============================================================================
// NUCLEOTIDE CLASS
// ============================================================================

const NucleotideID Nucleotide::charToNucleotideLookup[4] = {0, 1, 3, 2};
const char Nucleotide::charMask = 3;
const char Nucleotide::nucleotideToCharLookup[4] = {65, 67, 71, 84};
const NucleotideID Nucleotide::nucleotideMask = 3;

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

#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <cstdint>
#include <string>


// ============================================================================
// DEFINITIONS OF CONSTANTS
// ============================================================================

#define VERSION_NUMBER_COLUMBA 2     // The version number of Columba
#define SUB_VERSION_NUMBER_COLUMBA 0 // The sub version number of Columba
#define SUB_SUB_VERSION_NUMBER_COLUMBA 2 // The sub version number of Columba
#define COLUMBA_BUILD_INDEX_TAG 21   // the tag for building the index
#ifndef RUN_LENGTH_COMPRESSION
#define COLUMBA_FLAVOUR "VANILLA" // the flavour of Columba (vanilla)
#define DEFAULT_SEED_LENGTH 0     // The default seed length for the FM-index
#else
#define COLUMBA_FLAVOUR "RLC"   // the flavour of Columba (rlc)
#define DEFAULT_SEED_LENGTH 100 // The default seed length for the FM-index
#endif
// increment these tags if you change buildhelpers.cpp, buildfmindex.cpp or
// buildbmove.cpps

#define DEFAULT_SPARSENESS                                                     \
    4 // The default sparseness factor fo the suffix array

#define MAX_MAPQ 60 // The maximal map quality score in SAM format
#define MAX_K 13 // the max number of allowed errors for Columba search scheme
#define MAX_K_EDIT 20 // the max number of allowed errors for edit  distance
#define BEST_CUTOFF_COLUMBA MAX_K // The cutoff for the best mode in Columba

#define PE_STD_DEV_CONSIDERED                                                  \
    6 // The number of standard deviations to consider for paired end reads

#define PE_NUMBER_PAIRS_FOR_INFERENCE 750 // The number of pairs to use for
#define PE_MAX_READS_FOR_INFERENCE 10000  // The maximum number of reads to use

#define CIGAR_THRESHOLD                                                        \
    10 // The threshold for the CIGAR string calculation when converting from FM
       // occurrence to in-text occurrence

// ============================================================================
// AT-COMPILE-TIME DEFINITIONS
// ============================================================================

// The type of the length_t depending on compiler option
#if THIRTY_TWO
#define LENGTH_TYPE_NAME "32-bits"
typedef uint32_t length_t;
#else
#define LENGTH_TYPE_NAME "64-bits"
typedef uint64_t length_t;
#endif // THIRTY_TWO

// Check if the compiler is GCC or Clang and if the architecture is 64-bit for
// the use of uint128_t
#if (defined(__GNUC__) || defined(__clang__)) &&                               \
    (defined(__x86_64__) || __SIZEOF_POINTER__ == 8)
#define HAS_UINT128_T 1
#else
#define HAS_UINT128_T 0
#endif // (__GNUC__ || __clang__) && (__x86_64__ || __SIZEOF_POINTER__ == 8)


inline std::string getVersionString() {
    std::string str = std::to_string(VERSION_NUMBER_COLUMBA) + "." +
                      std::to_string(SUB_VERSION_NUMBER_COLUMBA);
    if (SUB_SUB_VERSION_NUMBER_COLUMBA > 0) {
        str += "." + std::to_string(SUB_SUB_VERSION_NUMBER_COLUMBA);
    }
    return str;
}

// ============================================================================
// ENUMS
// ============================================================================

/**
 * An enum for the direction of the search/a substring
 */
enum Direction { FORWARD, BACKWARD };

/**
 * Helper enum for finding the sequence name corresponding to the match of a
 * read.
 */
enum SeqNameFound { FOUND, FOUND_WITH_TRIMMING, NOT_FOUND };

// An enum for the partition strategies
enum PartitionStrategy { UNIFORM, STATIC, DYNAMIC };
// An enum for which distance metric to use
enum DistanceMetric { HAMMING, EDIT };
// An enum for the mapping mode
enum MappingMode { BEST, ALL };

// An enum for single or paired end reads
enum SequencingMode { SINGLE_END, PAIRED_END };

// An enum for the orientation of the paired end reads
enum Orientation { FR, RF, FF };

// An enum for the strand of the read
enum Strand { FORWARD_STRAND = 0, REVERSE_C_STRAND = 1 };

// An enum for the status in the pair of the read
enum PairStatus { FIRST_IN_PAIR = 0, SECOND_IN_PAIR = 1 };

#endif // DEFINITIONS_H

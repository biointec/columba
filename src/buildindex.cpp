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

#include "alphabet.h"    // for Alphabet
#include "bitvec.h"      // for Bitvec, Bitref
#include "definitions.h" // for length_t, COLUMBA_BUILD_INDEX_TAG, COLUMBA_...
#include "libsais64.h"   // for libsais64
#include "logger.h"      // for Logger, logger

#include <algorithm>  // for max, transform, equal, max_element, min_ele...
#include <array>      // for array
#include <cctype>     // for toupper
#include <cstdint>    // for int64_t, uint8_t, uint32_t
#include <functional> // for function
#include <iostream>   // for ifstream, ofstream, operator<<, ios, basic_...
#include <limits>     // for numeric_limits
#include <random>     // for minstd_rand
#include <stdexcept>  // for runtime_error
#include <stdlib.h>   // for size_t
#include <string>     // for string
#include <vector>     // for vector

#ifndef RUN_LENGTH_COMPRESSION
#include "fmindex/bwtrepr.h"     // for BWTRepresentation
#include "fmindex/encodedtext.h" // for EncodedText
#include "fmindex/suffixArray.h" // for SparseSuffixArray
#else
#ifdef BIG_BWT_USABLE
#include "../external/Big-BWT/utils.h" // for SABYTES
#endif
#include "bmove/moverepr.h"     // for MoveRepr
#include "bmove/moverow.h"      // for MoveRow
#include "bmove/plcp.h"         // for PLCP
#include "bmove/sparsebitvec.h" // for SparseB...

#endif

using namespace std;

// ============================================================================
// PARSING
// ============================================================================

#include "parameters/buildparameters.h" // for BuildParameters

// ============================================================================
// FUNCTIONS FOR BOTH VANILLA AND RLC
// ============================================================================

/**
 * @brief Replace non-ACGT characters with a random ACGT character.
 * @param original The original character.
 * @param gen The random number generator.
 * @param seed The seed for the random number generator. (dummy variable)
 * @param seedIndex The index of the seed. (dummy variable)
 * @return The replaced character.
 */
char replaceNonACGT(char original, std::minstd_rand& gen,
                    const std::string& seed, size_t& seedIndex) {
    const std::string validChars = "ACGT";
    if (original != 'A' && original != 'C' && original != 'G' &&
        original != 'T') {
        std::uniform_int_distribution<size_t> distribution(
            0, validChars.length() - 1);
        return validChars[distribution(gen)];
    }
    return original;
}

/**
 * @brief Replace non-ACGT characters with a seeded ACGT character.
 * @param original The original character.
 * @param gen The random number generator.
 * @param seed The seed for the random number generator.
 * @param seedIndex The index of the seed.
 * @return The replaced character.
 */
char replaceNonACGTWithSeed(char original, std::minstd_rand& gen,
                            const std::string& seed, size_t& seedIndex) {
    if (original != 'A' && original != 'C' && original != 'G' &&
        original != 'T') {
        char replacement = seed[seedIndex];
        seedIndex = (seedIndex + 1) %
                    seed.length(); // Move to the next character in the seed
        return replacement;
    }
    seedIndex = 0; // Reset seed index
    return original;
}

/**
 * @brief Concatenates and transforms the sequences from a FASTA file.
 *
 * This function reads a FASTA file and concatenates the sequences into a single
 * std::string. It also replaces non-ACGT characters with a random ACGT
 * character. The start positions of each sequence and the sequence names are
 * stored in separate std::vectors. The resulting concatenated std::string is
 * terminated with a '$' character.
 *
 * @param fastaFile The path to the FASTA file.
 * @param concatenation The resulting concatenated std::string. (output)
 * @param positions The std::vector to store the start positions of each
 * sequence. (output)
 * @param seqNames The std::vector to store the sequence names. (output)
 * @param replaceFunc The function to use for replacing non-ACGT characters.
 * @param expectedNumber The expected number of sequences in the FASTA file.
 * Default is 12000.
 */
void concatenateAndTransform(const std::string& fastaFile,
                             std::string& concatenation,
                             std::vector<length_t>& positions,
                             std::vector<std::string>& seqNames,
                             std::vector<length_t>& firstSeqIDPerFile,
                             std::function<char(char, size_t&)> replaceFunc,
                             length_t expectedNumber = 12000) {

    logger.logInfo("Reading FASTA file " + fastaFile);

    // Open input file and get its size
    std::ifstream inputFile(fastaFile, std::ios::ate);

    if (!inputFile.is_open()) {
        throw std::runtime_error("Error opening file: " + fastaFile);
    }

    std::streamsize fileSize = inputFile.tellg();
    inputFile.seekg(0, std::ios::beg);

    // Reserve memory for the concatenation string based on file size
    concatenation.reserve(static_cast<size_t>(fileSize + concatenation.size()));

    // Temporary variables for processing
    std::string line;
    length_t startPosition = concatenation.size();

    // Index for seed
    size_t seedIndex = 0;

    // Reserve memory for positions and seqNames vectors
    positions.reserve(expectedNumber + positions.size());
    seqNames.reserve(expectedNumber + seqNames.size());

    // A hash table of sequence names for O(1) lookup
    std::unordered_set<string> seqNamesLookup;

    firstSeqIDPerFile.push_back(seqNames.size());

    std::string sequence;
    bool sequenceName = false;

    while (std::getline(inputFile, line)) {
        if (line.empty())
            continue; // Skip empty lines

        if (line[0] == '>') { // Sequence name line
            // Process the previous sequence
            if (!sequence.empty()) {
                // Process the previous sequence
                for (char& c : sequence) {
                    c = replaceFunc(std::toupper(c), seedIndex);
                    concatenation += c;
                }

                positions.emplace_back(
                    startPosition); // push back the start position of this
                                    // sequence
                startPosition +=
                    sequence
                        .length(); // create start position for next sequence
                sequence.clear();  // Clear the sequence for the next one
            }
            std::string description = line.substr(1);
            // get the sequence name (before the first space)
            std::string name = description.substr(0, description.find(' '));

            // Check if the sequence name is already present
            if (seqNamesLookup.count(name) > 0) {
                throw std::runtime_error("Error: Sequence name " + name +
                                         " in file " + fastaFile +
                                         " is not unique!");
            }

            seqNames.emplace_back(name); // Store sequence name
            seqNamesLookup.insert(name);
            sequenceName = true;

        } else {
            // Append to the current sequence
            sequence += line;
        }
    }

    // Process the last sequence in the file
    for (char& c : sequence) {
        c = replaceFunc(std::toupper(c), seedIndex);
        concatenation += c;
    }
    positions.emplace_back(startPosition);

    if (!sequenceName) {
        seqNames.emplace_back(fastaFile);
    }

    // Close input file
    inputFile.close();
}

/**
 * @brief Read the contents of a text file into a std::string buffer.
 *
 */
void readText(const string& filename, string& buf) {
    ifstream ifs(filename);
    if (!ifs)
        throw runtime_error("Cannot open file: " + filename);

    ifs.seekg(0, ios::end);
    buf.resize(ifs.tellg());
    ifs.seekg(0, ios::beg);
    ifs.read((char*)buf.data(), buf.size());
}

/**
 * @brief Perform a sanity check on the suffix array.
 * @param T The text.
 * @param sa The suffix array.
 */
void sanityCheck(const string& T, vector<length_t>& sa) {
    logger.logInfo("Performing sanity checks...");
    // check T for correctness
    if (T.back() == '\n')
        throw runtime_error("T should end with a \'$\' character, "
                            "not with a newline");

    if (T.back() != '$')
        throw runtime_error("T should end with a \'$\' character");

    if (sa.size() != T.size())
        throw runtime_error("Text and suffix array contain a "
                            "different number of elements");

    // briefly check the suffix array
    length_t min = *min_element(sa.begin(), sa.end());
    length_t max = *max_element(sa.begin(), sa.end());

    if (min == 1 && max == T.size()) { // rebase to [0..T.size()-1]
        for (auto& el : sa)
            el--;
        min--;
        max--;
    }

    if (min != 0 || max != T.size() - 1)
        throw runtime_error("Suffix array must contain numbers between "
                            "[0 and " +
                            to_string(T.size() - 1) + "]");

    // check if all numbers in the suffix array are present
    Bitvec bv(sa.size());
    for (length_t i : sa)
        bv[i] = true;

    for (size_t i = 0; i < bv.size(); i++)
        if (!bv[i])
            throw runtime_error("Suffix " + to_string(i) +
                                " seems "
                                "to be missing from suffix array");

    // extra check:
    //      we could check T to see if the SA correctly sorts suffixes of T

    logger.logInfo("\tSanity checks OK");
}

/**
 * @brief Create and write the header information for output SAM files.
 * @param baseFN The base filename.
 * @param seqNames The sequence names.
 * @param positions The start positions of each sequence.
 */
void createAndWriteHeaderInfo(const string& baseFN,
                              const vector<string>& seqNames,
                              const vector<length_t>& positions) {
    logger.logInfo("Writing SAM header info for reference text to " + baseFN +
                   ".headerSN.bin...");
    // create the header lines with these reference sequences
    std::ofstream headerStream(baseFN + ".headerSN.bin", ios::binary);
    for (length_t i = 0; i < seqNames.size(); i++) {
        headerStream << "@SQ\tSN:" << seqNames[i]
                     << "\tLN:" << positions[i + 1] - positions[i] << "\n";
    }
    headerStream.close();
}

/**
 * @brief Write the start positions and sequence names to disk.
 * @param baseFN The base filename.
 * @param positions The start positions of each sequence.
 * @param seqNames The sequence names.
 */
void writePositionsAndSequenceNames(const string& baseFN,
                                    const vector<length_t>& positions,
                                    const vector<string>& seqNames,
                                    const vector<length_t>& firstSeqIDPerFile) {
    logger.logInfo("Write positions and sequence names to " + baseFN +
                   ".pos and " + baseFN + ".sna...");
    // Write the positions to disk
    std::ofstream ofs2(baseFN + ".pos", ios::binary);
    ofs2.write((char*)positions.data(), positions.size() * sizeof(length_t));
    ofs2.close();

    // Write the sequence names to disk
    std::ofstream ofs3(baseFN + ".sna", std::ios::binary);
    for (const auto& str : seqNames) {
        size_t len = str.size();
        ofs3.write(reinterpret_cast<const char*>(&len), sizeof(len));
        ofs3.write(str.c_str(), len);
    }
    ofs3.close();

    // Write the first sequence ID per file to disk
    std::ofstream ofs4(baseFN + ".fsid", ios::binary);
    ofs4.write((char*)firstSeqIDPerFile.data(),
               firstSeqIDPerFile.size() * sizeof(length_t));
    ofs4.close();
}

/**
 * @brief Check if the text size is within the limits of the length_t type.
 * @param tSize The size of the text.
 * @throws runtime_error if the text size is too large.
 */
void checkTextSize(const size_t tSize) {
    if (tSize > std::numeric_limits<length_t>::max()) {
        throw std::runtime_error(
            "The size of the text ()" + std::to_string(tSize) +
            ") is too large. Maximum size for the current "
            "compiler option is: " +
            std::to_string(std::numeric_limits<length_t>::max()));
    }
    if ((sizeof(length_t) * 8 == 64) &&
        (tSize <= std::numeric_limits<uint32_t>::max())) {
        logger.logWarning(
            "Program was compiled with 64-bit words, but the text size (" +
            std::to_string(tSize) +
            ") fits in a 32-bit word. Consider recompiling the program to "
            "improve performance and save memory.");
    }
}

/**
 * @brief Count the frequency of each character in the text.
 * @param T The text.
 * @param charCounts The vector to store the character counts.
 */
void countChars(const string& T, vector<length_t>& charCounts) {
    // initialize the counts to 0 and 256 elements
    charCounts = vector<length_t>(256, 0);
    for (char c : T)
        charCounts[(unsigned char)c]++;
}

/**
 * @brief Write the character counts to disk.
 * @param baseFN The base filename.
 * @param charCounts The character counts.
 */
void writeCharCounts(const string& baseFN, const vector<length_t>& charCounts) {
    ofstream ofs(baseFN + ".cct", ios::binary);
    ofs.write((char*)charCounts.data(), charCounts.size() * sizeof(length_t));
    ofs.close();
    logger.logInfo("Wrote file " + baseFN + ".cct");
}

/**
 * @brief Create the alphabet based on the character counts.
 * @param T The text.
 * @param charCounts The character counts.
 * @param sigma The alphabet. (output)
 */
void createAlphabet(const string& T, const vector<length_t>& charCounts,
                    Alphabet<ALPHABET>& sigma) {
    // count the number of unique characters in T
    int nUniqueChar = 0;
    for (length_t count : charCounts)
        if (count > 0)
            nUniqueChar++;

    logger.logInfo("Text has length " + std::to_string(T.size()));
    logger.logInfo("Text has " + std::to_string(nUniqueChar) +
                   " unique characters");

    if (nUniqueChar > ALPHABET) {
        logger.logError("FATAL: the number of unique characters in the "
                        "text exceeds the alphabet size. Please recompile "
                        "Columba using a higher value for ALPHABET");
        exit(EXIT_FAILURE);
    }

    if (nUniqueChar < ALPHABET) {
        logger.logWarning("the number of unique characters in the "
                          "text is less than the ALPHABET size specified when "
                          "Columba was compiled. Performance may be affected.");
    }

    sigma = Alphabet<ALPHABET>(charCounts);
}

/**
 * @brief Create the suffix array using the libsais library.
 * @param T The text.
 * @param SA The suffix array. (output)
 */
void createSuffixArray(const string& T, vector<length_t>& SA) {
    logger.logInfo("Generating the suffix array using libsais...");

    // create the suffix array with libsais
    // Convert std::string to const uint8_t*
    const uint8_t* tPtr = reinterpret_cast<const uint8_t*>(T.c_str());

    // Length of the input string
    int64_t n = static_cast<int64_t>(T.size());

    // Suffix array, size should be n (or n+fs if you want extra space, but 0 is
    // enough for most cases)
    std::vector<int64_t> suffixArray64(n, 0);
    // Call the libsais64 function
    int64_t result = libsais64(tPtr, suffixArray64.data(), n, 0, nullptr);
    if (result != 0) {
        throw runtime_error("libsais64 failed with error code " +
                            to_string(result));
    }
    // cast int64_t to length_t
    SA = vector<length_t>(suffixArray64.begin(), suffixArray64.end());

    logger.logInfo("Suffix array generated successfully!");
}

/**
 * @brief Create the BWT of the reversed text.
 * @param baseFN The base filename.
 * @param T The text.
 * @param revSA The suffix array of the reversed text.
 * @param sigma The alphabet.
 * @param rBWT The BWT of the reversed text. (output)
 */
void createRevBWT(const string& baseFN, const string& T,
                  const vector<length_t>& revSA,
                  const Alphabet<ALPHABET>& sigma, string& rBWT) {

    // build the reverse BWT
    logger.logInfo("Generating BWT of reversed text...");
    rBWT.resize(T.size());
    for (size_t i = 0; i < revSA.size(); i++)
        rBWT[i] = (revSA[i] > 0) ? T[T.size() - revSA[i]] : T.front();
}

/**
 * @brief Write a std::string to a binary file
 * @param name The filename.
 * @param str The std::string to write.
 */
void writeStringToBinaryFile(const std::string& name, const std::string& str) {
    std::ofstream outFile(name, std::ios::binary);
    if (!outFile) {
        std::cerr << "Error opening file for writing: " << name << std::endl;
        return;
    }

    // First, write the size of the string
    length_t size = str.size();
    outFile.write(reinterpret_cast<const char*>(&size), sizeof(size));
    // Then, write the string data itself
    outFile.write(str.c_str(), size);
    outFile.close();
}

/**
 * @brief Preprocess the FASTA files.
 * @param fastaFiles The paths to the FASTA files.
 * @param baseFN The base filename.
 * @param T The text. (output)
 * @param seedLength The length of the seed for seeded replacement.
 * @param noWriting Whether to write the text to disk. Default is false.
 */
void preprocessFastaFiles(const std::vector<std::string>& fastaFiles,
                          const std::string& baseFN, string& T,
                          length_t seedLength, bool noWriting = false) {

    // Seeded random number generation
    std::minstd_rand gen(42);

    // Generate random seed for seeded replacement
    std::string seed;
    size_t seedIndex = 0;
    for (length_t i = 0; i < seedLength; i++) {
        seed += replaceNonACGT('N', gen, seed, seedIndex);
    }

    // Reset gen to original seed
    gen.seed(42);

    std::function<char(char, size_t&)> replaceFunc;

    if (seedLength == 0) {

        replaceFunc = [&gen, &seed](char c, size_t& seedIndex) -> char {
            return replaceNonACGT(c, gen, seed, seedIndex);
        };

        logger.logInfo(
            "Using random (non-seeded) replacement for non-ACGT characters...");

    } else {

        replaceFunc = [&gen, &seed](char c, size_t& seedIndex) -> char {
            return replaceNonACGTWithSeed(c, gen, seed, seedIndex);
        };

        logger.logInfo("Using seeded replacement for non-ACGT characters, with "
                       "a seed length of " +
                       std::to_string(seedLength));
    }

    std::vector<length_t> positions; // the start positions of each
    std::vector<string> seqNames;    // the name of each sequence
    std::vector<length_t> firstSeqIDPerFile;

    for (const auto& fastaFile : fastaFiles) {
        logger.logInfo("Preprocessing FASTA file " + fastaFile);
        concatenateAndTransform(fastaFile, T, positions, seqNames,
                                firstSeqIDPerFile, replaceFunc);
    }
    // add the end of the text
    positions.emplace_back(T.size());
    // Ensure the concatenation ends with a dollar sign
    if (T.back() != '$') {
        T += '$';
    }
    logger.logInfo("Read and concatenated " + std::to_string(seqNames.size()) +
                   " sequence(s) from " + std::to_string(fastaFiles.size()) +
                   " file(s)");
    checkTextSize(T.size());

    if (!noWriting) {
        logger.logInfo("Writing concatenated uppercase sequence to disk...");
        writeStringToBinaryFile(baseFN + ".txt.bin", T);
    }

    createAndWriteHeaderInfo(baseFN, seqNames, positions);
    writePositionsAndSequenceNames(baseFN, positions, seqNames,
                                   firstSeqIDPerFile);
}

/**
 * @brief Write the meta information to a file.
 * @param baseFN The base filename.
 */
void writeMetaInfo(const string& baseFN) {
    ofstream metaFile(baseFN + ".meta");

    // Write the build tag to a file
    metaFile << COLUMBA_BUILD_INDEX_TAG << endl;
    // Write the 64 or 32-bit compiled info
    metaFile << sizeof(length_t) << endl;
    // Write the flavour of Columba
    metaFile << COLUMBA_FLAVOUR << endl;
    metaFile.close();
}

/**
 * @brief Generate the BWT of the text.
 * @param T The text.
 * @param SA The suffix array.
 * @param BWT The BWT. (output)
 */
void generateBWT(const string& T, const vector<length_t>& SA, string& BWT) {
    // build the BWT
    logger.logInfo("Generating BWT...");
    BWT.resize(T.size());
    for (size_t i = 0; i < SA.size(); i++)
        BWT[i] = (SA[i] > 0) ? T[SA[i] - 1] : T.back();
}

/**
 * @brief Write the character counts and create the alphabet.
 * @param baseFN The base filename.
 * @param T The text.
 * @param sigma The alphabet. (output)
 * @param charCounts The character counts.
 */
void writeCharCountsAndCreateAlphabet(const string& baseFN, const string& T,
                                      Alphabet<ALPHABET>& sigma,
                                      vector<length_t>& charCounts) {
    // count the frequency of each characters in T
    countChars(T, charCounts);
    // write the character counts table
    writeCharCounts(baseFN, charCounts);
    // Create the alphabet
    createAlphabet(T, charCounts, sigma);
}

/**
 * @brief Create the suffix array with a sanity check.
 * @param SA The suffix array. (output)
 * @param T The text.
 */
void createSAWithSanityCheck(vector<length_t>& SA, const string& T) {
    createSuffixArray(T, SA);

    // perform a sanity check on the suffix array
    sanityCheck(T, SA);
}

/**
 * @brief Create the suffix array of the reversed text with a sanity check.
 * @param revSA The suffix array of the reversed text. (output)
 * @param T The text.
 */
void createRevSAWithSanityCheck(vector<length_t>& revSA, const string& T) {

    std::string revT = T;
    std::reverse(revT.begin(), revT.end());

    // create the suffix array of the reverse text
    createSuffixArray(revT, revSA);
    revT.clear();

    // perform a sanity check on the suffix array
    sanityCheck(T, revSA);
}

/**
 * @brief Create the suffix array and BWT
 * @param T The text.
 * @param SA The suffix array. (output)
 * @param BWT The BWT. (output)
 */
void createSuffixArrayAndBWT(const string& T, vector<length_t>& SA,
                             string& BWT) {
    createSAWithSanityCheck(SA, T);
    generateBWT(T, SA, BWT);
}

#ifdef RUN_LENGTH_COMPRESSION
// ============================================================================
// FUNCTIONALITY FOR RLC
// ============================================================================
/**
 * Get run index of a given position.
 * @param position The position for which to find the run index.
 * @param size The amount of input-output intervals.
 * @param dPair Vector with pairs of input-output interval start positions.
 * @returns The run index of the given position.
 */
length_t getRunIndex(const length_t position, const length_t size,
                     const vector<MoveRow>& rows) {
    length_t rightBoundary = size - 1;
    length_t leftBoundary = 0;
    // Iteratively make the possible range smaller by binary search, until only
    // 1 interval remains.
    while (rightBoundary - leftBoundary >= 1) {
        // Use the middle of the possible range as a test value.
        length_t testIndex = ((rightBoundary + leftBoundary) / 2) + 1;

        // Eliminate half of the possible range by comparing the value to the
        // test value.
        if (rows[testIndex].inputStartPos <= position) {
            leftBoundary = testIndex;
        } else {
            rightBoundary = testIndex - 1;
        }
    }

    assert(position >= rows[leftBoundary].inputStartPos);
    assert(leftBoundary == size - 1 ||
           position < rows[leftBoundary + 1].inputStartPos);

    return leftBoundary;
}

/**
 * Fill dIndex array.
 * @param tIn Balanced tree: stores pairs of I_out with increasing input
 * interval index.
 * @param rows The vector to fill, rows[i] = input-output interval i +
 * index j of input interval containing q_i.
 * @param size The amount of input-output intervals.
 * @param bwtSize The BWT size.
 */
void fillRows(const map<length_t, pair<uint8_t, length_t>>& tIn,
              vector<MoveRow>& rows, const length_t size,
              const length_t bwtSize) {
    rows.resize(size + 1);

    map<length_t, pair<uint8_t, length_t>>::const_iterator tInIt = tIn.begin();
    length_t i = 0;
    while (tInIt != tIn.end()) {
        length_t inputStart = tInIt->first;
        uint8_t c = tInIt->second.first;
        length_t outputStart = tInIt->second.second;

        rows[i] = MoveRow(c, inputStart, outputStart, 0);

        i++;
        tInIt++;
    }

    tInIt = tIn.begin();
    i = 0;
    while (tInIt != tIn.end()) {
        length_t outputRunIndex = getRunIndex(tInIt->second.second, size, rows);

        rows[i].outputStartRun = outputRunIndex;

        i++;
        tInIt++;
    }

    rows[size] = MoveRow(0, bwtSize, bwtSize, size);
}

/**
 * Creates the Move structure.
 * @param baseFN the base filename to write the move structure to.
 * @param BWT The burrows-wheeler transform in string format (uncompressed).
 * @param charCounts Accumulated number of characters in lex order.
 * @param sigma The alphabet.
 */
length_t createAndWriteMove(const string& baseFN, const string& BWT,
                            const vector<length_t>& charCounts,
                            const Alphabet<ALPHABET>& sigma) {

    // Get the alphabet size
    size_t S = sigma.size();

    // Set bwt size
    length_t bwtSize = BWT.size();

    // balanced tree: stores pairs of I_out with increasing input interval index
    map<length_t, pair<uint8_t, length_t>> tIn;

    // Create accumulated charCounts
    vector<length_t> charCountsAcc(S, 0);
    length_t total = 0;
    for (size_t i = 0; i < S; i++) {
        char c = sigma.i2c(i);
        charCountsAcc[i] = total;
        total += charCounts[c];
    }

    // fill tIn, tOut
    vector<length_t> charsVisited(S, 0);
    length_t prevC = S;
    length_t zeroCharPos;
    for (length_t i = 0; i < bwtSize; i++) {

        length_t c = sigma.c2i(BWT[i]);
        if (c == 0) {
            zeroCharPos = i;
        }
        if (prevC != c) {
            length_t lf = charCountsAcc[c] + charsVisited[c];
            tIn[i] = make_pair(c, lf);
        }

        charsVisited[c]++;
        prevC = c;
    }

    // Set dPair and dIndex size
    length_t arraySize = tIn.size();
    length_t size = arraySize;
    logger.logInfo("There are " + to_string(size) + " runs.");

    // Fill dPair and dIndex
    vector<MoveRow> rows;
    fillRows(tIn, rows, size, BWT.size());

    // Write Move structures to files
    string fileName = baseFN;
    fileName += ".move";
    ofstream ofs(fileName, ios::binary);

    // Write bwtSize, amount of input-output intervals and the alphabet size to
    // the output stream.
    ofs.write((char*)&bwtSize, sizeof(bwtSize));
    ofs.write((char*)&size, sizeof(size));
    ofs.write((char*)&zeroCharPos, sizeof(zeroCharPos));

    // Write rows to the output stream.
    for (length_t i = 0; i < size + 1; i++) {
        rows[i].serialize(ofs);
    }

    ofs.close();
    logger.logInfo("Wrote file " + fileName);

    return size;
}

void writeIntVectorBinary(const std::string& filename,
                          const std::vector<length_t>& array) {
    // convert to int_vector
    uint8_t width =
        (uint8_t)ceil(log2(*max_element(array.begin(), array.end())));
    sdsl::int_vector<> intVector(array.size(), 0, width);
    for (size_t i = 0; i < array.size(); i++) {
        intVector[i] = array[i];
    }
    std::ofstream ofs(filename, std::ios::binary);
    intVector.serialize(ofs);
    ofs.close();
}

/**
 * Builds samplesFirst and samplesLast.
 * @param samplesFirst Vector to fill. First sa samples of each input interval.
 * @param samplesLast Vector to fill. Last sa samples of each input interval.
 * @param SA The suffix array.
 * @param BWT The burrows-wheeler transform in string format (uncompressed).
 */
void buildSamples(vector<length_t>& samplesFirst, vector<length_t>& samplesLast,
                  const vector<length_t>& SA, const string& BWT) {

    samplesFirst.push_back(SA[0] > 0 ? SA[0] - 1 : BWT.size() - 1);
    for (size_t pos = 0; pos < BWT.size() - 1; pos++) {
        if (BWT[pos] != BWT[pos + 1]) {
            samplesLast.push_back(SA[pos] > 0 ? SA[pos] - 1 : BWT.size() - 1);
            samplesFirst.push_back(SA[pos + 1] > 0 ? SA[pos + 1] - 1
                                                   : BWT.size() - 1);
        }
    }
    samplesLast.push_back(SA[BWT.size() - 1] > 0 ? SA[BWT.size() - 1] - 1
                                                 : BWT.size() - 1);
}

/**
 * Generate predecessor structures
 * @param samplesFirst samplesFirst array
 * @param samplesLast samplesLast array
 * @param [out] predFirst predecessor bitvector of the samplesFirst array
 * @param [out] predLast predecessor bitvector of the samplesLast array
 * @param textLength length of the text
 * @param verbose print verbose output to stdout
 */
void generatePredecessors(const vector<length_t>& samplesFirst,
                          const vector<length_t>& samplesLast,
                          SparseBitvec& predFirst, SparseBitvec& predLast,
                          const length_t textLength) {

    vector<bool> predFirstBV(textLength, false);
    vector<bool> predLastBV(textLength, false);

    for (length_t SApos : samplesFirst) {
        predFirstBV[SApos] = true;
    }

    for (length_t SApos : samplesLast) {
        predLastBV[SApos] = true;
    }

    predFirst = SparseBitvec(predFirstBV);
    predLast = SparseBitvec(predLastBV);
}

/**
 * Generate predToRun arrays
 * @param samplesFirst samplesFirst array
 * @param samplesLast samplesLast array
 * @param [out] firstToRun mapping between rank of ones in predFirst bitvector
 * and run indices
 * @param [out] lastToRun mapping between rank of ones in predLast bitvector and
 * run indices
 */
void generatePredToRun(const vector<length_t>& samplesFirst,
                       const vector<length_t>& samplesLast,
                       vector<length_t>& firstToRun,
                       vector<length_t>& lastToRun) {

    firstToRun.resize(samplesFirst.size());
    iota(firstToRun.begin(), firstToRun.end(), 0);
    sort(firstToRun.begin(), firstToRun.end(),
         [&samplesFirst](length_t a, length_t b) {
             return samplesFirst[a] < samplesFirst[b];
         });

    lastToRun.resize(samplesLast.size());
    iota(lastToRun.begin(), lastToRun.end(), 0);
    sort(lastToRun.begin(), lastToRun.end(),
         [&samplesLast](length_t a, length_t b) {
             return samplesLast[a] < samplesLast[b];
         });
}

void createIndex(const string& T, const BuildParameters& params,
                 const Alphabet<ALPHABET>& sigma,
                 const vector<length_t>& charCounts) {
    // aliasing
    const auto& baseFN = params.baseFN;

    {
        // create SA and BWT
        vector<length_t> SA;
        string BWT;
        createSuffixArrayAndBWT(T, SA, BWT);
        { // generate and write PLCP array
            logger.logInfo(
                "Generating the permuted longest common prefix array...");
            PLCP plcp(T, SA);
            plcp.write(baseFN + ".plcp");
            logger.logInfo("Wrote file " + baseFN + ".plcp");
        }
        // Create samplesLast and samplesFirst
        logger.logInfo("Sampling the suffix array values at run boundaries...");
        vector<length_t> samplesFirst;
        vector<length_t> samplesLast;
        buildSamples(samplesFirst, samplesLast, SA, BWT);

        // Clear the suffix array
        SA.clear();

        // write samplesFirst and samplesLast to file
        writeIntVectorBinary(baseFN + ".smpf", samplesFirst);
        logger.logInfo("Wrote file " + baseFN + ".smpf");
        writeIntVectorBinary(baseFN + ".smpl", samplesLast);
        logger.logInfo("Wrote file " + baseFN + ".smpl");

        SparseBitvec predFirst;
        SparseBitvec predLast;
        vector<length_t> firstToRun;
        vector<length_t> lastToRun;

        // Generate the predecessor bitvectors
        logger.logInfo(
            "Generating the predecessor bitvectors for the samples...");
        generatePredecessors(samplesFirst, samplesLast, predFirst, predLast,
                             BWT.size());
        // Generate the predToRun arrays
        logger.logInfo(
            "Mapping the predecessor bits to their corresponding runs...");
        generatePredToRun(samplesFirst, samplesLast, firstToRun, lastToRun);

        samplesFirst.clear();
        samplesLast.clear();

        // Write predecessor structures
        predFirst.write(baseFN + ".prdf");
        logger.logInfo("Wrote file " + baseFN + ".prdf");
        predLast.write(baseFN + ".prdl");
        logger.logInfo("Wrote file " + baseFN + ".prdl");
        writeIntVectorBinary(baseFN + ".ftr", firstToRun);
        logger.logInfo("Wrote file " + baseFN + ".ftr");
        writeIntVectorBinary(baseFN + ".ltr", lastToRun);
        logger.logInfo("Wrote file " + baseFN + ".ltr");

        firstToRun.clear();
        lastToRun.clear();

        // Create the Move structure
        logger.logInfo("Creating the move table...");
        createAndWriteMove(baseFN, BWT, charCounts, sigma);
        BWT.clear();
    }

    { // read or create the reverse suffix array
        vector<length_t> revSA;
        createRevSAWithSanityCheck(revSA, T);

        // build the reverse BWT
        string revBWT;
        createRevBWT(baseFN, T, revSA, sigma, revBWT);

        // Create samplesLast and samplesFirst
        logger.logInfo(
            "Sampling the reverse suffix array values at run boundaries...");
        vector<length_t> revSamplesFirst;
        vector<length_t> revSamplesLast;
        buildSamples(revSamplesFirst, revSamplesLast, revSA, revBWT);

        // Clear the reverse suffix array
        revSA.clear();

        // write samplesFirst and samplesLast to file
        writeIntVectorBinary(baseFN + ".rev.smpf", revSamplesFirst);
        logger.logInfo("Wrote file " + baseFN + ".rev.smpf");
        writeIntVectorBinary(baseFN + ".rev.smpl", revSamplesLast);
        logger.logInfo("Wrote file " + baseFN + ".rev.smpl");

        revSamplesFirst.clear();
        revSamplesLast.clear();

        // Create the Move structure
        logger.logInfo("Creating the reverse move table...");
        createAndWriteMove(baseFN + ".rev", revBWT, charCounts, sigma);

        // Clear the reverse BWT
        revBWT.clear();
    }
}

#ifdef BIG_BWT_USABLE

void readSuffixArrayFile(const std::string& baseFN,
                         const std::string& extension,
                         vector<length_t>& samples, vector<length_t>& toRun,
                         length_t size, length_t nrOfRuns, bool reverse) {
    logger.logInfo("Reading suffix array samples from " + baseFN + extension +
                   "...");

    string fileName = baseFN + extension;

    // Open the file
    FILE* file = fopen(fileName.c_str(), "rb");

    samples = vector<length_t>(nrOfRuns, 0);
    toRun = vector<length_t>(nrOfRuns, 0);
    std::vector<std::pair<length_t, length_t>> pos_run_pairs(nrOfRuns);

    length_t* buf = new length_t[1];

    for (length_t i = 0; i < nrOfRuns; ++i) {
        // Read the first SABYTES bytes into buf from file
        if (fread(buf, SABYTES, 1, file) != 1)
            throw runtime_error((fileName + " read failed."));

        // Read the next SABYTES bytes into buf from file
        if (fread(buf, SABYTES, 1, file) != 1)
            throw runtime_error((fileName + " read failed."));

        // Calculate sa_val from buf: take the first byte, ensure it's within
        // range, and adjust as needed based on the data size
        length_t sa_val = buf[0] % (1UL << (8 * SABYTES));
        if (!reverse) {
            // Only for the regular suffix array, sa_val must be corrected.
            // For the reverse case, the value is already correct since Big-BWT
            // puts the sentinel character at the end of the reverse text as
            // well.
            sa_val = (sa_val > 0) ? (sa_val - 1) : (size - 1);
        }

        // Store sa_val in samples array at index i
        samples[i] = sa_val;

        // Store {sa_val, i} pair in pos_run_pairs array
        pos_run_pairs[i] = {sa_val, i};
    }

    delete[] buf;

    if (!reverse) {

        logger.logInfo(
            "Creating the mapping between the predecessor bits and the "
            "runs...");

        std::sort(pos_run_pairs.begin(), pos_run_pairs.end());

        std::vector<length_t> positions;
        for (length_t i = 0; i < nrOfRuns; ++i) {
            positions.push_back(pos_run_pairs[i].first);
            toRun[i] = pos_run_pairs[i].second;
        }
    }

    fclose(file);
}

#endif

void constructRunLengthEncodedPLCP(const SparseBitvec& first,
                                   const vector<length_t>& first_to_run,
                                   const std::vector<length_t>& charCounts,
                                   length_t size, length_t r, PLCP& plcp,
                                   const std::string& baseFN,
                                   const Alphabet<ALPHABET>& sigma) {
    logger.logInfo("Constructing run-length encoded PLCP...");

    // Construct cumulative char counts
    vector<length_t> charCountsCumulative(256, 0);
    for (size_t i = 1; i < 256; i++) {
        auto temp = charCountsCumulative[i - 1] + charCounts[i];
        charCountsCumulative[i] = temp;
    }

    logger.logDeveloper("Reading " + baseFN + ".move...");

    MoveRepr rows;
    rows.load(baseFN, false);

    logger.logDeveloper("Creating bitvectors to support select for each "
                        "character in the BWT...");

    sdsl::bit_vector char_bv[ALPHABET];
    sdsl::bit_vector::select_1_type select[ALPHABET];

    {

        for (length_t i = 0; i < ALPHABET; i++) {
            char_bv[i] = sdsl::bit_vector(size, 0);
        }

        for (length_t i = 0; i < r; i++) {
            size_t c = rows.getRunHead(i);
            for (length_t j = rows.getInputStartPos(i);
                 j < rows.getInputStartPos(i + 1); j++) {
                char_bv[c][j] = true;
            }
        }

        for (length_t i = 0; i < ALPHABET; i++) {
            select[i] = sdsl::bit_vector::select_1_type(&char_bv[i]);
        }
    }

    logger.logDeveloper("Creating PLCP...");

    std::vector<length_t> ones, zeros;

    char c, c0;
    length_t p, p0, pos, l, gap;
    length_t prev_pos = 0;
    length_t prev_l = 0;
    length_t acc0 = 0, acc1 = 0;

    try {

        for (length_t i = 0; i < r; ++i) {
            // Log a progress message that takes r into account (e.g 100 times
            // from 0 to r-1)
            if (i % (r / 100) == 0) {
                logger.logDeveloper("Progress: " + std::to_string(i) + " / " +
                                    std::to_string(r) + " (" +
                                    std::to_string((i * 100) / r) + "%)");
            }

            try {
                pos = first.select(i);
            } catch (const std::exception& e) {
                std::cerr << "Error: " << e.what() << '\n';
                std::cerr << "While trying to perform select on first, with i: "
                          << i << '\n';
                exit(1);
            }

            gap = i == 0 ? pos : (pos - prev_pos - 1);

            l = 0;

            length_t tempIdx = first_to_run[i];
            p = rows.getInputStartPos(tempIdx);
            c = sigma.i2c(rows.getRunHead(tempIdx));
            rows.findLF(p, tempIdx);
            if (p != 0) {
                p0 = p - 1;
                c0 = (std::upper_bound(charCountsCumulative.begin(),
                                       charCountsCumulative.end(), p0) -
                      charCountsCumulative.begin());
                while (c == c0) {
                    l++;
                    length_t i = p - charCountsCumulative[c - 1];
                    try {
                        p = select[sigma.c2i(c)](i + 1);
                    } catch (const std::exception& e) {
                        std::cerr << "Error: " << e.what() << '\n';
                        std::cerr << "While trying to perform select on "
                                     "sigma.c2i(c), with i: "
                                  << i << " and c: " << c << '\n';
                        exit(1);
                    }
                    c = (std::upper_bound(charCountsCumulative.begin(),
                                          charCountsCumulative.end(), p) -
                         charCountsCumulative.begin());
                    length_t i0 = p0 - charCountsCumulative[c0 - 1];
                    try {
                        p0 = select[sigma.c2i(c0)](i0 + 1);
                    } catch (const std::exception& e) {
                        std::cerr << "Error: " << e.what() << '\n';
                        std::cerr << "While trying to perform select on "
                                     "sigma.c2i(c0), with i0: "
                                  << i0 << " and c0: " << c0 << '\n';
                        exit(1);
                    }
                    c0 = (std::upper_bound(charCountsCumulative.begin(),
                                           charCountsCumulative.end(), p0) -
                          charCountsCumulative.begin());
                }
            }

            if (i == 0) {
                if (l + gap > 0) {
                    ones.push_back(acc1);
                    acc0 += l + gap - 1;
                    zeros.push_back(acc0);
                    acc1 += gap + 1;
                } else {
                    acc1 += gap + 1;
                }
            } else if (l + gap + 1 - prev_l) {
                ones.push_back(acc1);
                acc0 += l + gap + 1 - prev_l;
                zeros.push_back(acc0);
                acc1 += gap + 1;
            } else {
                acc1 += gap + 1;
            }

            prev_pos = pos;
            prev_l = l;
        }

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        std::cerr << "While trying to construct PLCP" << '\n';
        exit(1);
    }

    ones.push_back(acc1);

    logger.logDeveloper("Copying the PCLP into its correct object...");

    plcp = PLCP(size, ones, zeros);

    logger.logInfo("Constructed run-length encoded PLCP.");
}

#ifdef BIG_BWT_USABLE

void createIndexPFP(const string& baseFN) {

    // build the BWT
    string BWT;
    // Read the BWT from disk
    logger.logInfo("Reading " + baseFN + ".bwt...");
    readText(baseFN + ".bwt", BWT);

    // Replace '\0' with '$'
    std::replace(BWT.begin(), BWT.end(), '\0', '$');

    length_t bwtSize = BWT.size();

    // count the frequency of each characters in T
    vector<length_t> charCounts;
    countChars(BWT, charCounts);

    // write the character counts table
    writeCharCounts(baseFN, charCounts);

    // Create the alphabet
    Alphabet<ALPHABET> sigma;
    createAlphabet(BWT, charCounts, sigma);

    { // Create the Move structure
        logger.logInfo("Creating the move table...");
        length_t nrOfRuns = createAndWriteMove(baseFN, BWT, charCounts, sigma);

        BWT.clear();

        vector<length_t> samplesFirst;
        vector<length_t> firstToRun;
        readSuffixArrayFile(baseFN, ".ssa", samplesFirst, firstToRun, bwtSize,
                            nrOfRuns, false);

        vector<length_t> samplesLast;
        vector<length_t> lastToRun;
        readSuffixArrayFile(baseFN, ".esa", samplesLast, lastToRun, bwtSize,
                            nrOfRuns, false);

        // write samplesFirst and samplesLast to file
        writeIntVectorBinary(baseFN + ".smpf", samplesFirst);
        logger.logInfo("Wrote file " + baseFN + ".smpf");
        writeIntVectorBinary(baseFN + ".smpl", samplesLast);
        logger.logInfo("Wrote file " + baseFN + ".smpl");
        writeIntVectorBinary(baseFN + ".ftr", firstToRun);
        logger.logInfo("Wrote file " + baseFN + ".ftr");
        writeIntVectorBinary(baseFN + ".ltr", lastToRun);
        logger.logInfo("Wrote file " + baseFN + ".ltr");

        lastToRun.clear();

        SparseBitvec predFirst;
        SparseBitvec predLast;

        // Generate the predecessor bitvectors
        logger.logInfo(
            "Generating the predecessor bitvectors for the samples...");
        generatePredecessors(samplesFirst, samplesLast, predFirst, predLast,
                             bwtSize);

        samplesFirst.clear();
        samplesLast.clear();

        // Write predecessor structures
        predFirst.write(baseFN + ".prdf");
        logger.logInfo("Wrote file " + baseFN + ".prdf");
        predLast.write(baseFN + ".prdl");
        logger.logInfo("Wrote file " + baseFN + ".prdl");

        // generate and write PLCP array
        logger.logInfo(
            "Generating the permuted longest common prefix array...");
        PLCP plcp;
        constructRunLengthEncodedPLCP(predFirst, firstToRun, charCounts,
                                      bwtSize, nrOfRuns, plcp, baseFN, sigma);
        plcp.write(baseFN + ".plcp");
        logger.logInfo("Wrote file " + baseFN + ".plcp");

        firstToRun.clear();
    }

    logger.logInfo("Switching to reversed text...");

    {
        // build the BWT
        string revBWT;
        // Read the BWT from disk
        logger.logInfo("Reading " + baseFN + ".rev.bwt...");
        readText(baseFN + ".rev.bwt", revBWT);

        // Replace '\0' with '$'
        std::replace(revBWT.begin(), revBWT.end(), '\0', '$');

        // Create the Move structure
        logger.logInfo("Creating the move table...");
        length_t nrOfRuns =
            createAndWriteMove(baseFN + ".rev", revBWT, charCounts, sigma);

        // Clear the reverse BWT
        revBWT.clear();

        vector<length_t> revSamplesFirst;
        vector<length_t> revFirstToRun;
        readSuffixArrayFile(baseFN, ".rev.ssa", revSamplesFirst, revFirstToRun,
                            bwtSize, nrOfRuns, true);

        revFirstToRun.clear();

        vector<length_t> revSamplesLast;
        vector<length_t> revLastToRun;
        readSuffixArrayFile(baseFN, ".rev.esa", revSamplesLast, revLastToRun,
                            bwtSize, nrOfRuns, true);
        revLastToRun.clear();

        // write samplesFirst and samplesLast to file
        writeIntVectorBinary(baseFN + ".rev.smpf", revSamplesFirst);
        logger.logInfo("Wrote file " + baseFN + ".rev.smpf");
        writeIntVectorBinary(baseFN + ".rev.smpl", revSamplesLast);
        logger.logInfo("Wrote file " + baseFN + ".rev.smpl");

        revSamplesFirst.clear();
        revSamplesLast.clear();
    }

    writeMetaInfo(baseFN);
}

int indexConstructingAfterPFP(const BuildParameters& params) {
    logger.logInfo(
        "Starting index construction after prefix-free parsing step...");

    std::array<string, 6> requiredExtensionsPFP = {
        ".bwt", ".esa", ".ssa", ".rev.bwt", ".rev.esa", ".rev.ssa"};

    for (auto& ext : requiredExtensionsPFP) {
        string filename = params.baseFN + ext;

        // check if file with filename exists
        ifstream ifs(filename);
        if (!ifs) {
            logger.logError("Missing file: " + filename);
            logger.logError("Please run the prefix-free parsing step first.");
            return EXIT_FAILURE;
        }
    }

    createIndexPFP(params.baseFN);

    logger.logInfo("Index construction completed successfully!");
    logger.logInfo("Exiting... bye!");

    return EXIT_SUCCESS;
}

int preprocessingOnly(const BuildParameters& params) {
    logger.logInfo("Preprocessing only...");

    std::string T; // the (concatenated) text
    preprocessFastaFiles(params.fastaFiles, params.baseFN, T, params.seedLength,
                         true);

    // Remove the sentinel character
    T.pop_back();

    // Write the preprocessed text to disk
    logger.logInfo("Writing concatenated uppercase sequence to disk...");
    std::ofstream ofs(params.baseFN);
    ofs << T;
    ofs.close();

    // Reverse the text
    logger.logInfo("Reversing text...");
    std::reverse(T.begin(), T.end());

    // Write the reversed text to disk
    logger.logInfo("Writing reversed text to disk...");
    ofs = std::ofstream(params.baseFN + ".rev");
    ofs << T;
    ofs.close();

    logger.logInfo("Fasta input preprocessed successfully!");
    logger.logInfo("Exiting... bye!");

    return EXIT_SUCCESS;
}
#endif
#else

// ============================================================================
// FUNCTIONALITY FOR VANILLA
// ============================================================================

void writeBWT(const Alphabet<ALPHABET>& sigma, const string& baseFN,
              string& BWT) {

    // Encode bwt
    logger.logInfo("Encoding BWT...");
    EncodedText<ALPHABET> eBWT(sigma, BWT);

    eBWT.write(baseFN + ".bwt");

    logger.logInfo("Wrote file " + baseFN + ".bwt");
}

void writeSA(const string& baseFN, const vector<length_t>& SA, length_t saSF) {
    SparseSuffixArray sparseSA(SA, saSF);
    sparseSA.write(baseFN);
    logger.logInfo("Wrote sparse suffix array with factor " +
                   std::to_string(saSF));
}

void writeSA(const string& baseFN, const vector<length_t>& SA) {
    for (int saSF = 1; saSF <= 128; saSF *= 2) {
        writeSA(baseFN, SA, saSF);
    }
}

void writeBWTBitvectors(const string& baseFN, const string& BWT,
                        const Alphabet<ALPHABET>& sigma, bool reverse) {
    logger.logInfo("Creating BWT bitvectors...");
    BWTRepresentation<ALPHABET> bwt(sigma, BWT);
    string extension = (reverse) ? ".rev.brt" : ".brt";
    bwt.write(baseFN + extension);
    logger.logInfo("Wrote file " + baseFN + extension);
}

void createIndex(const string& T, const BuildParameters& params,
                 const Alphabet<ALPHABET>& sigma,
                 const vector<length_t>& charCounts) {

    // aliasing
    const auto& baseFN = params.baseFN;

    // create SA and BWT
    vector<length_t> SA;
    string BWT;
    createSuffixArrayAndBWT(T, SA, BWT);
    writeBWT(sigma, baseFN, BWT);

    // create sparse suffix arrays
    if (params.allSparsenessFactors) {
        writeSA(baseFN, SA);
    } else {
        writeSA(baseFN, SA, params.sparsenessFactor);
    }
    SA.clear();

    // create succinct BWT bitvector table
    writeBWTBitvectors(baseFN, BWT, sigma, false);
    BWT.clear();

    logger.logInfo("Switching to reversed text...");

    // read or create the reverse suffix array
    vector<length_t> revSA;
    createRevSAWithSanityCheck(revSA, T);

    // build the reverse BWT
    string rBWT;
    createRevBWT(baseFN, T, revSA, sigma, rBWT);
    revSA.clear();
    writeBWTBitvectors(baseFN, rBWT, sigma, true);
    rBWT.clear();
}
#endif

// ============================================================================
// INPUT FUNCTIONALITY FOR VANILLA AND RLC
// ============================================================================

/**
 * Process the input fasta files and create the index.
 * @param params the build parameters
 *
 */
void processFastaFiles(const BuildParameters& params) {
    std::string T; // the (concatenated) text
    preprocessFastaFiles(params.fastaFiles, params.baseFN, T,
                         params.seedLength);

    // create alphabet and charCounts
    Alphabet<ALPHABET> sigma;
    vector<length_t> charCounts;
    writeCharCountsAndCreateAlphabet(params.baseFN, T, sigma, charCounts);
    createIndex(T, params, sigma, charCounts);
    // write meta info
    writeMetaInfo(params.baseFN);
}

int main(int argc, char* argv[]) {

    logger.logInfo("Welcome to Columba's index construction!");

    BuildParameters params;
    if (!BuildParameters::parse(argc, argv, params)) {
        BuildParameters::showUsage();
        return EXIT_FAILURE;
    }

    logger.logInfo("Alphabet size is " + std::to_string(ALPHABET - 1) + " + 1");

#if defined(RUN_LENGTH_COMPRESSION) && defined(BIG_BWT_USABLE)
    if (params.preprocessOnly) {
        return preprocessingOnly(params);
    }

    if (params.pfp) {
        return indexConstructingAfterPFP(params);
    }
#elif !defined(RUN_LENGTH_COMPRESSION)
    if (params.allSparsenessFactors) {
        logger.logInfo("Using all sparseness factors 1 to 128");
    } else {
        logger.logInfo("Using sparseness factor " +
                       std::to_string(params.sparsenessFactor));
    }
#endif

    try {
        for (auto& fastaFile : params.fastaFiles) {
            if (!BuildParameters::validFastaExtension(fastaFile)) {
                throw runtime_error("Invalid fasta file extension for file " +
                                    fastaFile);
            }
        }
        processFastaFiles(params);
    } catch (const std::exception& e) {
        logger.logError("Fatal: " + string(e.what()));

        return EXIT_FAILURE;
    }

    logger.logInfo("Index construction completed successfully!");
    logger.logInfo("Exiting... bye!");

    return EXIT_SUCCESS;
}

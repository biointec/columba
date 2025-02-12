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
#include "logger.h"      // for Logger, logger
#include "util.h"        // for fileExists..

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

#ifdef THIRTY_TWO
#include "libsais.h" // for libsais
#include "libsais64.h"
#else
#include "divsufsort64.h"
#endif

#ifdef HAVE_ZLIB
#include "seqfile.h"
#include <zlib.h> // for gzFile, gzopen, gzclose, gzread, gzwrite
#endif

#ifndef RUN_LENGTH_COMPRESSION
#include "fmindex/bwtrepr.h"     // for BWTRepresentation
#include "fmindex/encodedtext.h" // for EncodedText
#include "fmindex/suffixArray.h" // for SparseSuffixArray
#else
#ifdef BIG_BWT_USABLE
#include "../external/Big-BWT/utils.h" // for SABYTES
#endif
#include "bmove/moverepr.h"     // for MoveLFRepr
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

bool isGzipped(FileType fileType, const string& filename) {
    if (fileType == FileType::FASTA_GZ) {
#ifndef HAVE_ZLIB
        // If the file is gzipped but zlib support is not available, throw an
        // error
        throw runtime_error("Error: " + fastaFile +
                            " is a gzipped FASTA file, but Columba was not "
                            "compiled with zlib support.");
#else
        // return true if the file is gzipped and zlib support is
        // available
        return true;
#endif
    }
    return false;
}

SeqFile openFastaFile(const string& fastaFile) {
    // Determine the file type of the provided FASTA file
    FileType fileType;
    tie(fileType, ignore) = getFileType(fastaFile);

    // Validate that the file is either a FASTA or gzipped FASTA file
    if (fileType != FileType::FASTA && fileType != FileType::FASTA_GZ) {
        throw runtime_error("Error: " + fastaFile +
                            " is not a valid FASTA file.");
    }

    // Validate that the FASTA file exists
    if (!Util::fileExists(fastaFile)) {
        throw runtime_error("Error: " + fastaFile + " does not exist.");
    }

    // Create a SeqFile object and open the FASTA file
    SeqFile file(isGzipped(fileType, fastaFile));
    file.open(fastaFile);

    return file;
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

    SeqFile file = openFastaFile(fastaFile);

    // Reserve memory for the concatenation string based on estimated file size
    concatenation.reserve(
        static_cast<size_t>(file.estimateASCIISize() + concatenation.size()));

    // Reserve memory for positions and seqNames vectors to avoid reallocations
    positions.reserve(expectedNumber + positions.size());
    seqNames.reserve(expectedNumber + seqNames.size());

    // Record the index of the first sequence ID in the file
    firstSeqIDPerFile.push_back(seqNames.size());

    // Temporary variables for processing sequence data
    std::string sequence;      // To hold the current sequence
    bool sequenceName = false; // Flag to indicate if a seq name was present
    std::string line;          // To hold each line read from the file
    length_t startPosition = concatenation.size(); // start position in concat.
    size_t seedIndex = 0; // Index for seed to replace non-ACGT characters

    while (file.good()) {
        file.getLine(line);

        if (line.empty() || (line.size() == 1 && line[0] == '\n'))
            continue; // Skip empty lines

        // pop back the new line character
        if (line.back() == '\n')
            line.pop_back();

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
            if (std::find(seqNames.begin(), seqNames.end(), name) !=
                seqNames.end()) {
                throw std::runtime_error("Error: Sequence name " + name +
                                         " in file " + fastaFile +
                                         " is not unique!");
            }

            seqNames.emplace_back(name); // Store sequence name
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
    file.close();
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
 * @param checkBothLimits Check both the lower and upper limits. Default is
 * true. If false, only the upper limit is checked.
 * @throws runtime_error if the text size is too large.
 */
void checkTextSize(const size_t tSize, bool checkBothLimits = true) {
    if (tSize > std::numeric_limits<length_t>::max()) {
        throw std::runtime_error(
            "The size of the reference (" + std::to_string(tSize) +
            ") is too large. Maximum size for the current "
            "compiler option is: " +
            std::to_string(std::numeric_limits<length_t>::max()) +
            "\nPlease recompile Columba with a larger word size.");
    }

    if (checkBothLimits && (sizeof(length_t) * 8 == 64) &&
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

    if (T.size() == 1) {
        logger.logError(
            "Reference text is empty except for sentinel character. Please "
            "provide a valid reference with the -f or -F flag.");
        exit(EXIT_FAILURE);
    }

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

#ifdef THIRTY_TWO
void createSuffixArrayLibsais(const string& T, vector<uint32_t>& SA) {
    logger.logInfo("Generating the suffix array using libsais...");

    // Convert std::string to const uint8_t*
    const uint8_t* tPtr = reinterpret_cast<const uint8_t*>(T.c_str());

    if (T.size() > INT32_MAX) {
        // We need libsais64 to handle this
        // create a temporary vector to store the 64-bit suffix array
        SA.clear();
        vector<int64_t> SA64(T.size());
        logger.logDeveloper("Calling libsais64");
        int64_t result = libsais64(tPtr, SA64.data(),
                                   static_cast<int64_t>(T.size()), 0, nullptr);
        if (result != 0) {
            throw runtime_error("libsais64 failed with error code " +
                                to_string(result));
        }
        SA.reserve(T.size());
        // convert the 64-bit suffix array to 32-bit
        for (size_t i = 0; i < T.size(); i++) {
            SA.emplace_back(static_cast<uint32_t>(SA64[i]));
        }
    } else {
        // Resize SA to the size of the input string to hold the suffix array
        SA.resize(T.size());
        logger.logDeveloper("Calling libsais");
        int32_t result = libsais(tPtr, reinterpret_cast<int32_t*>(SA.data()),
                                 static_cast<int32_t>(T.size()), 0, nullptr);

        if (result != 0) {
            throw runtime_error("libsais failed with error code " +
                                to_string(result));
        }
    }

    logger.logDeveloper("Libsais successful");
}
#else
void createSuffixArrayDivsufsort64(const string& T, vector<uint64_t>& SA) {
    logger.logInfo("Generating the suffix array using divsufsort...");

    // Resize SA to the size of the input string to hold the suffix array
    SA.resize(T.size());

    // Convert std::string to const uint8_t*
    const sauchar_t* tPtr = reinterpret_cast<const sauchar_t*>(T.c_str());

    logger.logDeveloper("Calling divsufsort64");
    saint_t result = divsufsort64(tPtr, reinterpret_cast<int64_t*>(SA.data()),
                                  static_cast<int64_t>(T.size()));

    if (result != 0) {
        throw std::runtime_error("divsufsort failed with error code " +
                                 std::to_string(result));
    }

    logger.logDeveloper("Divsufsort successful");
}
#endif

/**
 * @brief Create the suffix array using the libsais library.
 * @param T The text.
 * @param SA The suffix array. (output)
 */
void createSuffixArray(const string& T, vector<length_t>& SA) {
#ifdef THIRTY_TWO
    createSuffixArrayLibsais(T, SA);
#else
    createSuffixArrayDivsufsort64(T, SA);

#endif

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
        checkTextSize(T.size(), false);
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
void createRevSAWithSanityCheck(vector<length_t>& revSA, string& T) {

    if (T.size() < (2ull << 32)) {
        std::string revT = T;
        std::reverse(revT.begin(), revT.end());

        // create the suffix array of the reverse text
        createSuffixArray(revT, revSA);
        revT.clear();
    } else {
        std::reverse(T.begin(), T.end());
        createSuffixArray(T, revSA);
        std::reverse(T.begin(), T.end());
    }

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
                     const MoveLFReprBP& rows) {
    length_t rightBoundary = size - 1;
    length_t leftBoundary = 0;
    // Iteratively make the possible range smaller by binary search, until only
    // 1 interval remains.
    while (rightBoundary - leftBoundary >= 1) {
        // Use the middle of the possible range as a test value.
        length_t testIndex = ((rightBoundary + leftBoundary) / 2) + 1;

        // Eliminate half of the possible range by comparing the value to the
        // test value.
        if (rows.getInputStartPos(testIndex) <= position) {
            leftBoundary = testIndex;
        } else {
            rightBoundary = testIndex - 1;
        }
    }

    assert(position >= rows.getInputStartPos(leftBoundary));
    assert(leftBoundary == size - 1 ||
           position < rows.getInputStartPos(leftBoundary + 1));

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
              MoveLFReprBP& rows, const length_t size, const length_t bwtSize) {
    rows.initialize(tIn.size(), bwtSize);

    length_t i = 0;
    for (auto it = tIn.begin(); it != tIn.end(); ++it) {
        rows.setRowValues(i, it->second.first, it->first, it->second.second, 0);
        assert(rows.getRunHead(i) == it->second.first);
        assert(rows.getInputStartPos(i) == it->first);
        assert(rows.getOutputStartPos(i) == it->second.second);
        i++;
    }

    i = 0;
    for (auto it = tIn.begin(); it != tIn.end(); ++it) {
        length_t outputRunIndex = getRunIndex(it->second.second, size, rows);
        rows.setOutputStartRun(i, outputRunIndex);
        assert(rows.getOutputStartRun(i) == outputRunIndex);
        i++;
    }

    rows.setRowValues(size, 0, bwtSize, bwtSize, size);
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
    length_t zeroCharPos = bwtSize;
    for (length_t i = 0; i < bwtSize; i++) {

        char _c = BWT[i];
        length_t c = sigma.c2i(_c);
        assert(c < S);
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
    MoveLFReprBP rows;
    rows.setZeroCharPos(zeroCharPos);
    fillRows(tIn, rows, size, BWT.size());

    string moveFileName = baseFN + ".LFBP";
    rows.write(moveFileName);
    logger.logInfo("\tWrote file " + moveFileName);

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
 * Builds samplesFirst and samplesLast from the suffix array (SA) and BWT.
 * Fills samplesFirst and samplesLast at character changes in BWT for length_t.
 *
 * @param samplesFirst The samplesFirst array to fill.
 * @param samplesLast The samplesLast array to fill.
 * @param SA The suffix array.
 * @param BWT The Burrows-Wheeler transform.
 */
void buildSamples(vector<length_t>& samplesFirst, vector<length_t>& samplesLast,
                  const vector<length_t>& SA, const string& BWT) {

    samplesFirst.emplace_back(SA[0]);
    for (size_t pos = 0; pos < BWT.size() - 1; pos++) {
        if (BWT[pos] != BWT[pos + 1]) {
            samplesLast.emplace_back(SA[pos]);
            samplesFirst.emplace_back(SA[pos + 1]);
        }
    }
    samplesLast.emplace_back(SA[BWT.size() - 1]);
}

/**
 * Initiates the build for the MovePhi structure.
 *
 * @param tIn The tIn map.
 * @param tOut The tOut map.
 * @param SamplesFirst The samplesFirst array.
 * @param SamplesLast The samplesLast array.
 */
void startBuildForPhiMove(map<length_t, length_t>& tIn,
                          map<length_t, length_t>& tOut,
                          vector<length_t>& samplesFirst,
                          vector<length_t>& samplesLast) {

    tIn[samplesFirst.front()] = samplesLast.back();
    tOut[samplesLast.back()] = samplesFirst.front();
    // Remove the front of the samplesFirst array
    samplesFirst.erase(samplesFirst.begin());
    // Remove the back of the samplesLast array
    samplesLast.pop_back();
    while (!samplesFirst.empty()) {
        tIn[samplesFirst.back()] = samplesLast.back();
        tOut[samplesLast.back()] = samplesFirst.back();
        samplesFirst.pop_back();
        samplesLast.pop_back();
    }
}

/**
 * Generate predecessor structures
 * @param movePhi MovePhi array
 * @param movePhiInv MovePhiInv array
 * @param [out] predFirst predecessor bitvector of the samplesFirst array
 * @param [out] predLast predecessor bitvector of the samplesLast array
 * @param textLength length of the text
 */
void generatePredecessors(const vector<length_t>& movePhi,
                          const vector<length_t>& movePhiInv,
                          SparseBitvec& predFirst, SparseBitvec& predLast,
                          const length_t textLength) {

    vector<bool> predFirstBV(textLength + 1, false);
#ifndef PHI_MOVE
    vector<bool> predLastBV(textLength + 1, false);
#endif // PHI_MOVE

    for (const length_t& row : movePhi) {
        predFirstBV[row > 0 ? row - 1 : textLength - 1] = true;
    }

#ifndef PHI_MOVE
    for (const length_t& row : movePhiInv) {
        predLastBV[row > 0 ? row - 1 : textLength - 1] = true;
    }
#endif // PHI_MOVE

    predFirst = SparseBitvec(predFirstBV);
#ifndef PHI_MOVE
    predLast = SparseBitvec(predLastBV);
#endif // PHI_MOVE
}

/**
 * Generate predecessor structures
 * @param movePhi MovePhi array
 * @param [out] predFirst predecessor bitvector of the samplesFirst array
 * @param textLength length of the text
 */
void generatePredecessors(const MovePhiReprBP& move, SparseBitvec& pred,
                          const length_t textLength) {

    vector<bool> predFirstBV(textLength + 1, false);

    for (length_t i = 0; i < move.size(); i++) {
        predFirstBV[move.getInputStartPos(i) > 0 ? move.getInputStartPos(i) - 1
                                                 : textLength - 1] = true;
    }

    pred = SparseBitvec(predFirstBV);
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
                       vector<length_t>& lastToRun, length_t textLength) {

    firstToRun.resize(samplesFirst.size());
    iota(firstToRun.begin(), firstToRun.end(), 0);
    sort(firstToRun.begin(), firstToRun.end(),
         [&samplesFirst, &textLength](length_t a, length_t b) {
             return (samplesFirst[a] > 0 ? samplesFirst[a] - 1
                                         : textLength - 1) <
                    (samplesFirst[b] > 0 ? samplesFirst[b] - 1
                                         : textLength - 1);
         });

    lastToRun.resize(samplesLast.size());
    iota(lastToRun.begin(), lastToRun.end(), 0);
    sort(lastToRun.begin(), lastToRun.end(),
         [&samplesLast, &textLength](length_t a, length_t b) {
             return (samplesLast[a] > 0 ? samplesLast[a] - 1 : textLength - 1) <
                    (samplesLast[b] > 0 ? samplesLast[b] - 1 : textLength - 1);
         });
}

/**
 * @brief Fill tUnbalanced with unbalanced intervals
 *
 * @param tUnbalanced The map to fill with unbalanced intervals
 * @param tIn The input intervals
 * @param tOut The output intervals
 * @param maxOverlap The maximum overlap allowed
 */
void fillTUnbalanced(map<length_t, length_t>& tUnbalanced,
                     const map<length_t, length_t>& tIn,
                     const map<length_t, length_t>& tOut, length_t maxOverlap) {
    auto tInIt = tIn.begin();
    auto tOutIt = tOut.begin();

    while (tOutIt != tOut.end()) {
        length_t overlapCount = 0;
        length_t outStart = tOutIt->first;
        length_t outEnd = (std::next(tOutIt) != tOut.end())
                              ? std::next(tOutIt)->first
                              : numeric_limits<length_t>::max();

        // Skip input intervals that are before the current output interval
        while (tInIt != tIn.end() && tInIt->first < outStart) {
            ++tInIt;
        }

        // Determine if there is overlap or if the interval just started
        if (tInIt == tIn.end() || outStart < tInIt->first) {
            overlapCount = 1;
        }

        // Count overlapping input intervals
        while (tInIt != tIn.end() && tInIt->first < outEnd &&
               overlapCount < maxOverlap) {
            ++overlapCount;
            ++tInIt;
        }

        // Mark as unbalanced if overlap exceeds maxOverlap
        if (overlapCount >= maxOverlap) {
            tUnbalanced[tOutIt->second] = outStart;
        }

        ++tOutIt;
    }
}

/**
 * @brief Balance the intervals
 *
 * @param tUnbalanced The unbalanced intervals
 * @param rows The phi move array (output)
 * @param tIn The input intervals
 * @param tOut The output intervals
 * @param maxOverlap The maximum overlap allowed
 * @param textLength The length of the text
 */
void balanceIntervals(map<length_t, length_t>& tUnbalanced, MovePhiReprBP& rows,
                      map<length_t, length_t> tIn, map<length_t, length_t> tOut,
                      length_t maxOverlap, length_t textLength) {
    while (!tUnbalanced.empty()) {
        // Get the first unbalanced interval
        auto unbalancedIt = tUnbalanced.begin();
        length_t outStart = unbalancedIt->second;
        length_t inputStart = unbalancedIt->first;

        // Compute the interval difference for splitting
        auto nextInputIt = tIn.upper_bound(outStart);
        ++nextInputIt; // Todolore: this is only correct for maxOverlap = 4
        length_t intervalDiff = nextInputIt->first - outStart;

        // Insert the new intervals
        tIn.emplace(inputStart + intervalDiff, outStart + intervalDiff);
        tOut.emplace(outStart + intervalDiff, inputStart + intervalDiff);

        // Remove the balanced interval from unbalanced
        tUnbalanced.erase(inputStart);

        // Check overlapping intervals for rebalancing
        length_t outputsToCheck[] = {
            tOut.lower_bound(inputStart + intervalDiff) != tOut.begin()
                ? std::prev(tOut.lower_bound(inputStart + intervalDiff))->first
                : numeric_limits<length_t>::max(),
            outStart, outStart + intervalDiff};

        for (length_t outputStart : outputsToCheck) {
            auto endIt = tOut.upper_bound(outputStart);
            length_t endInterval = (endIt != tOut.end())
                                       ? endIt->first
                                       : numeric_limits<length_t>::max();

            // Count overlapping input intervals
            length_t overlapCount = 1;
            auto overlapInputIt = tIn.upper_bound(outputStart);
            while (overlapInputIt != tIn.end() &&
                   overlapInputIt->first < endInterval &&
                   overlapCount <= maxOverlap) {
                ++overlapCount;
                ++overlapInputIt;
            }

            // Rebalance if needed
            if (overlapCount > maxOverlap) {
                tUnbalanced[tOut[outputStart]] = outputStart;
            }
        }
    }

    rows.initialize(tIn.size(), textLength);
    length_t i = 0;
    for (auto it = tIn.begin(); it != tIn.end(); ++it) {
        rows.setRowValues(i, it->first, it->second, 0);
        assert(rows.getInputStartPos(i) == it->first);
        assert(rows.getOutputStartPos(i) == it->second);
        i++;
    }
    length_t size = tIn.size();
    rows.setRowValues(i, textLength, textLength, size);
}

void createUnbalancedMoveTable(MovePhiReprBP& rows,
                               const map<length_t, length_t>& tIn,
                               length_t textLength) {

    rows.initialize(tIn.size(), textLength);
    length_t i = 0;
    for (auto it = tIn.begin(); it != tIn.end(); ++it) {
        rows.setRowValues(i, it->first, it->second, 0);
        assert(rows.getInputStartPos(i) == it->first);
        assert(rows.getOutputStartPos(i) == it->second);
        i++;
    }
    length_t size = tIn.size();
    rows.setRowValues(i, textLength, textLength, size);
}

/**
 * Generate the mapping from start indices in both output intervals (phi and phi
 * inverse) to their corresponding run indices
 *
 * @param move Move structure
 * @param pred predecessor bitvector of the start samples
 * @param textLength length of the text
 */
void generatePhiRunMapping(MovePhiReprBP& move, const SparseBitvec& pred,
                           length_t textLength) {
    for (length_t i = 0; i < move.size(); i++) {
        length_t textPos = move.getOutputStartPos(i);
        length_t runIndexInText = pred.rank(textPos);
        move.setOutputStartRun(i, runIndexInText);
        assert(move.getOutputStartRun(i) == runIndexInText);
    }
}

#ifdef BIG_BWT_USABLE

void readSuffixArrayFile(const std::string& baseFN,
                         const std::string& extension,
                         vector<length_t>& samples, vector<length_t>& toRun,
                         length_t size, length_t nrOfRuns, bool reverse,
                         bool makeToRun) {
    logger.logInfo("Reading suffix array samples from " + baseFN + extension +
                   "...");

    string fileName = baseFN + extension;

    // Open the file
    FILE* file = fopen(fileName.c_str(), "rb");

    samples = vector<length_t>(nrOfRuns, 0);
    if (makeToRun) {
        toRun = vector<length_t>(nrOfRuns, 0);
    }
    std::vector<std::pair<length_t, length_t>> pos_run_pairs(nrOfRuns);

    uint64_t* buf = new uint64_t[1];

    for (length_t i = 0; i < nrOfRuns; ++i) {
        // Read the first SABYTES bytes into buf from file
        if (fread(buf, SABYTES, 1, file) != 1)
            throw runtime_error((fileName + " read failed."));

        // Read the next SABYTES bytes into buf from file
        if (fread(buf, SABYTES, 1, file) != 1)
            throw runtime_error((fileName + " read failed."));

        // Calculate sa_val from buf: take the first byte, ensure it's within
        // range, and adjust as needed based on the data size
        uint64_t sa_val = buf[0] % (1UL << (8 * SABYTES));
        if (reverse) {
            // Only for the reverse suffix array, sa_val must be corrected.
            // For the reverse case, the value is too small since Big-BWT
            // puts the sentinel character at the end of the reverse text as
            // well.
            sa_val = (sa_val < size - 1) ? (sa_val + 1) : 0;
        }

        // Store sa_val in samples array at index i
        assert(sa_val >= 0 && sa_val < size);
        samples[i] = (length_t)sa_val;

        // Store {sa_val, i} pair in pos_run_pairs array
        pos_run_pairs[i] = {sa_val > 0 ? (length_t)sa_val - 1 : size - 1, i};
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
            if (makeToRun) {
                toRun[i] = pos_run_pairs[i].second;
            }
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

    // Construct cumulative char counts more efficiently
    vector<length_t> charCountsCumulative(256, 0);
    for (size_t i = 1; i < 256; ++i) {
        charCountsCumulative[i] = charCountsCumulative[i - 1] + charCounts[i];
    }

    logger.logDeveloper("Reading " + baseFN + ".LFBP...");
    MoveLFReprBP rows;
    if (!rows.load(baseFN)) {
        logger.logError("Error: Could not load " + baseFN + ".LFBP");
        exit(1);
    }

    logger.logDeveloper("Creating bitvectors for select support...");

    // Create bitvectors and select structures
    sdsl::bit_vector char_bv[ALPHABET];
    sdsl::bit_vector::select_1_type select[ALPHABET];

    for (length_t i = 0; i < ALPHABET; ++i) {
        char_bv[i] = sdsl::bit_vector(size, 0); // Initialize bit vector
    }

    logger.logDeveloper("Initialization successful.");

    for (length_t i = 0; i < r; ++i) {
        size_t c = rows.getRunHead(i);
        for (length_t j = rows.getInputStartPos(i);
             j < rows.getInputStartPos(i + 1); ++j) {
            char_bv[c][j] = true; // Mark positions with true
        }
    }

    logger.logDeveloper("Marked positions in bitvectors.");

    for (length_t i = 0; i < ALPHABET; ++i) {
        select[i] = sdsl::bit_vector::select_1_type(
            &char_bv[i]); // Initialize select support
    }

    logger.logDeveloper("Creating PLCP...");

    std::vector<length_t> ones, zeros;
    ones.reserve(r + 1);
    zeros.reserve(r + 1);

    length_t pos = 0, gap = 0, l = 0, p = 0, p0 = 0, prev_pos = 0, prev_l = 0;
    length_t acc0 = 0, acc1 = 0;

    for (length_t i = 0; i < r; ++i) {
        if (r > 100 && i % (r / 100) == 0) {
            logger.logDeveloper("Progress: " + std::to_string(i) + " / " +
                                std::to_string(r) + " (" +
                                std::to_string((i * 100.0) / r) + "%)");
        }

        // Perform select on the sparse bit vector
        try {
            pos = first.select(i);
        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << " while selecting index " << i
                      << '\n';
            exit(1);
        }

        gap = (i == 0) ? pos : (pos - prev_pos - 1);

        l = 0;
        length_t tempIdx = first_to_run[i];
        p = rows.getInputStartPos(tempIdx);
        char c = sigma.i2c(rows.getRunHead(tempIdx));
        rows.findLF(p, tempIdx);

        if (p != 0) {
            p0 = p - 1;
            char c0 = std::upper_bound(charCountsCumulative.begin(),
                                       charCountsCumulative.end(), p0) -
                      charCountsCumulative.begin();

            // Loop until characters diverge
            while (c == c0) {
                l++;
                length_t i = p - charCountsCumulative[c - 1];

                // Perform select on the current character
                try {
                    p = select[sigma.c2i(c)](i + 1);
                } catch (const std::exception& e) {
                    std::cerr << "Error: " << e.what()
                              << " while selecting character " << c << '\n';
                    exit(1);
                }

                c = std::upper_bound(charCountsCumulative.begin(),
                                     charCountsCumulative.end(), p) -
                    charCountsCumulative.begin();
                length_t i0 = p0 - charCountsCumulative[c0 - 1];

                // Perform select on the previous character
                try {
                    p0 = select[sigma.c2i(c0)](i0 + 1);
                } catch (const std::exception& e) {
                    std::cerr << "Error: " << e.what()
                              << " while selecting character " << c0 << '\n';
                    exit(1);
                }

                c0 = std::upper_bound(charCountsCumulative.begin(),
                                      charCountsCumulative.end(), p0) -
                     charCountsCumulative.begin();
            }
        }

        // Update ones and zeros based on gap and l values
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

    // Final push to ensure last values are stored
    ones.push_back(acc1);

    logger.logDeveloper("Copying the PLCP into its correct object...");
    plcp = PLCP(size, ones, zeros);

    logger.logInfo("Constructed run-length encoded PLCP.");
}

void writePhiMoveStructures(MovePhiReprBP& phiMove, const string& baseFN,
                            const string& suffix, bool inverse,
                            length_t bwtSize) {
    // Generate and write predecessor bitvector
    logger.logInfo("\tGenerating the predecessor bitvector...");
    SparseBitvec pred;
    generatePredecessors(phiMove, pred, bwtSize);

    string predFileName =
        baseFN + ".move." + suffix + "." + (inverse ? "prdl" : "prdf");
    pred.write(predFileName);
    logger.logInfo("\tWrote file " + predFileName);

    // Generate and write phi run mapping
    logger.logInfo("\tGenerating the mapping between output interval start "
                   "positions and their run identifiers...");
    generatePhiRunMapping(phiMove, pred, bwtSize);

    string moveFileName = baseFN + ".phiBP." + suffix + (inverse ? ".inv" : "");
    phiMove.write(moveFileName);
    logger.logInfo("\tWrote file " + moveFileName);
}

void processPhiMoveStructures(const string& baseFN,
                              map<length_t, length_t>& inputIntervals,
                              map<length_t, length_t>& outputIntervals,
                              bool inverse, length_t bwtSize,
                              length_t maxOverlap = 4) {

    { // Balanced move table
        logger.logInfo("\tCollecting the unbalanced intervals...");
        map<length_t, length_t> unbalancedIntervals;
        fillTUnbalanced(unbalancedIntervals, inputIntervals, outputIntervals,
                        maxOverlap);

        logger.logInfo("\tBalancing the intervals...");
        MovePhiReprBP phiMove;
        balanceIntervals(unbalancedIntervals, phiMove, inputIntervals,
                         outputIntervals, maxOverlap, bwtSize);
        unbalancedIntervals.clear();

        writePhiMoveStructures(phiMove, baseFN, "balanced", inverse, bwtSize);
    }
}

void processPhiAndPhiInverse(const string& baseFN,
                             vector<length_t>& samplesFirst,
                             vector<length_t>& samplesLast, length_t bwtSize) {
    // Initialize move structures for phi
    map<length_t, length_t> inputIntervals;
    map<length_t, length_t> outputIntervals;

    logger.logInfo("Initiating the samples for the phi move structures...");
    startBuildForPhiMove(inputIntervals, outputIntervals, samplesFirst,
                         samplesLast);

    samplesFirst.clear();
    samplesLast.clear();

    length_t maxOverlap = 4;

    // Log creation of move table for phi
    logger.logInfo("Creating the move table for phi...");
    processPhiMoveStructures(baseFN, inputIntervals, outputIntervals, false,
                             bwtSize, maxOverlap);

    // Log creation of move table for phi inverse
    logger.logInfo("Creating the move table for phi inverse...");
    processPhiMoveStructures(baseFN, outputIntervals, inputIntervals, true,
                             bwtSize, maxOverlap);

    inputIntervals.clear();
    outputIntervals.clear();
}

#ifdef BIG_BWT_USABLE

// Optimized helper function to replace multiple characters and log their counts
void replaceSentinel(std::string& str) {
    logger.logInfo("Replacing special characters in the BWT...");

    size_t countNull = 0; // Count for '\0'
    size_t countOne = 0;  // Count for '\1'
    size_t countTwo = 0;  // Count for '\2'

    // Traverse the string once to count and replace
    for (auto& ch : str) {
        if (ch == '\0') {
            countNull++; // Count '\0'
            ch = '$';    // Replace with '$'
        } else if (ch == '\1') {
            countOne++; // Count '\1'
            ch = '$';   // Replace with '$'
        } else if (ch == '\2') {
            countTwo++; // Count '\2'
            ch = '$';   // Replace with '$'
        }
    }

    // Log results after the single pass
    logger.logDeveloper("Found " + std::to_string(countNull) +
                        " \\0 characters.");
    logger.logDeveloper("Found " + std::to_string(countOne) +
                        " \\1 characters.");
    logger.logDeveloper("Found " + std::to_string(countTwo) +
                        " \\2 characters.");

    size_t totalCount = countNull + countOne + countTwo;

    // Error check
    if (totalCount != 1) {
        logger.logError("Found " + std::to_string(totalCount) +
                        " special characters in the BWT. Expected 1.");
        exit(1);
    }
}

#endif // BIG_BWT_USABLE

void processSamplesAndPreds(const string& baseFN,
                            const vector<length_t>& samplesFirst,
                            const vector<length_t>& samplesLast,
                            const string& BWT) {
    vector<length_t> firstToRun;
    vector<length_t> lastToRun;

    // Generate the predToRun arrays
    logger.logInfo(
        "Mapping the predecessor bits to their corresponding runs...");
    generatePredToRun(samplesFirst, samplesLast, firstToRun, lastToRun,
                      BWT.size());

    SparseBitvec predFirst;
    SparseBitvec predLast;

    logger.logInfo("Generating the predecessor bitvectors for the samples...");
    generatePredecessors(samplesFirst, samplesLast, predFirst, predLast,
                         BWT.size());

    // Write predecessor structures
    predFirst.write(baseFN + ".og.prdf");
    logger.logInfo("Wrote file " + baseFN + ".og.prdf");
    predLast.write(baseFN + ".og.prdl");
    logger.logInfo("Wrote file " + baseFN + ".og.prdl");

    // Generate and write PLCP array
    writeIntVectorBinary(baseFN + ".ftr", firstToRun);
    logger.logInfo("Wrote file " + baseFN + ".ftr");
    writeIntVectorBinary(baseFN + ".ltr", lastToRun);
    logger.logInfo("Wrote file " + baseFN + ".ltr");

    firstToRun.clear();
    lastToRun.clear();
}

void createIndex(string& T, const BuildParameters& params,
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

#ifndef PHI_MOVE
        processSamplesAndPreds(baseFN, samplesFirst, samplesLast, BWT);
#else
        processPhiAndPhiInverse(baseFN, samplesFirst, samplesLast, BWT.size());
#endif // PHI_MOVE

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

void createIndexPFP(const string& baseFN) {

    // build the BWT
    string BWT;
    // Read the BWT from disk
    logger.logInfo("Reading " + baseFN + ".bwt...");
    readText(baseFN + ".bwt", BWT);

    replaceSentinel(BWT);

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

#ifdef PHI_MOVE
        bool makeToRun = false;
#else
        bool makeToRun = true;
#endif // PHI_MOVE

        vector<length_t> samplesFirst;
        vector<length_t> firstToRun;
        readSuffixArrayFile(baseFN, ".ssa", samplesFirst, firstToRun, bwtSize,
                            nrOfRuns, false, true);

        vector<length_t> samplesLast;
        vector<length_t> lastToRun;
        readSuffixArrayFile(baseFN, ".esa", samplesLast, lastToRun, bwtSize,
                            nrOfRuns, false, makeToRun);

        // write samplesFirst and samplesLast to file
        writeIntVectorBinary(baseFN + ".smpf", samplesFirst);
        logger.logInfo("Wrote file " + baseFN + ".smpf");
        writeIntVectorBinary(baseFN + ".smpl", samplesLast);
        logger.logInfo("Wrote file " + baseFN + ".smpl");

#ifndef PHI_MOVE
        // write firstToRun and lastToRun to file
        writeIntVectorBinary(baseFN + ".ftr", firstToRun);
        logger.logInfo("Wrote file " + baseFN + ".ftr");
        writeIntVectorBinary(baseFN + ".ltr", lastToRun);
        logger.logInfo("Wrote file " + baseFN + ".ltr");
#endif

        lastToRun.clear();

        { // Original phi operation + PLCP
            SparseBitvec predFirst;
            SparseBitvec predLast;

            // Generate the predecessor bitvectors
            logger.logInfo(
                "Generating the predecessor bitvectors for the samples...");
            generatePredecessors(samplesFirst, samplesLast, predFirst, predLast,
                                 bwtSize);

#ifndef PHI_MOVE
            // Write predecessor structures
            predFirst.write(baseFN + ".og.prdf");
            logger.logInfo("Wrote file " + baseFN + ".og.prdf");
            predLast.write(baseFN + ".og.prdl");
            logger.logInfo("Wrote file " + baseFN + ".og.prdl");
#endif

            // generate and write PLCP array
            logger.logInfo(
                "Generating the permuted longest common prefix array...");
            PLCP plcp;
            constructRunLengthEncodedPLCP(predFirst, firstToRun, charCounts,
                                          bwtSize, nrOfRuns, plcp, baseFN,
                                          sigma);
            plcp.write(baseFN + ".plcp");
            logger.logInfo("Wrote file " + baseFN + ".plcp");

            firstToRun.clear();
        }

#ifdef PHI_MOVE
        processPhiAndPhiInverse(baseFN, samplesFirst, samplesLast, bwtSize);
#endif
    }

    logger.logInfo("Switching to reversed text...");

    {
        // build the BWT
        string revBWT;
        // Read the BWT from disk
        logger.logInfo("Reading " + baseFN + ".rev.bwt...");
        readText(baseFN + ".rev.bwt", revBWT);

        replaceSentinel(revBWT);

        // Create the Move structure
        logger.logInfo("Creating the move table...");
        length_t nrOfRuns =
            createAndWriteMove(baseFN + ".rev", revBWT, charCounts, sigma);

        // Clear the reverse BWT
        revBWT.clear();

        vector<length_t> revSamplesFirst;
        vector<length_t> revFirstToRun;
        readSuffixArrayFile(baseFN, ".rev.ssa", revSamplesFirst, revFirstToRun,
                            bwtSize, nrOfRuns, true, false);

        revFirstToRun.clear();

        vector<length_t> revSamplesLast;
        vector<length_t> revLastToRun;
        readSuffixArrayFile(baseFN, ".rev.esa", revSamplesLast, revLastToRun,
                            bwtSize, nrOfRuns, true, false);
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

/**
 * @brief Write the sparse suffix array with a given sparseness factor.
 * @param baseFN The base filename.
 * @param SA The suffix array.
 * @param saSF The sparseness factor.
 */
void writeSA(const string& baseFN, const vector<length_t>& SA, length_t saSF) {
    SparseSuffixArray sparseSA(SA, saSF);
    sparseSA.write(baseFN);
    logger.logInfo("Wrote sparse suffix array with factor " +
                   std::to_string(saSF));
}

/**
 * @brief Write the sparse suffix array with all sparseness factors.
 * @param baseFN The base filename.
 * @param SA The suffix array.
 */
void writeSA(const string& baseFN, const vector<length_t>& SA) {
    for (int saSF = 1; saSF <= 128; saSF *= 2) {
        writeSA(baseFN, SA, saSF);
    }
}

/**
 * @brief Create the suffix array and BWT and write the sparse suffix array with
 * the sparsenessfactor defined in the BuildParameters.
 * @param T The text.
 * @param SA The suffix array. (output)
 * @param BWT The BWT. (output)
 * @param params The build parameters.
 */
void createSuffixArrayAndBWTAndWriteSparseSA(const string& T,
                                             vector<length_t>& SA, string& BWT,
                                             const BuildParameters& params) {
    createSAWithSanityCheck(SA, T);

    // create sparse suffix arrays
    if (params.allSparsenessFactors) {
        writeSA(params.baseFN, SA);
    } else {
        writeSA(params.baseFN, SA, params.sparsenessFactor);
    }

    generateBWT(T, SA, BWT);
}

void writeBWT(const Alphabet<ALPHABET>& sigma, const string& baseFN,
              string& BWT) {

    // Encode bwt
    logger.logInfo("Encoding BWT...");
    EncodedText<ALPHABET> eBWT(sigma, BWT);

    eBWT.write(baseFN + ".bwt");

    logger.logInfo("Wrote file " + baseFN + ".bwt");
}

void writeBWTBitvectors(const string& baseFN, const string& BWT,
                        const Alphabet<ALPHABET>& sigma, bool reverse) {
    logger.logInfo("Creating BWT bitvectors...");
    BWTRepresentation<ALPHABET> bwt(sigma, BWT);
    string extension = (reverse) ? ".rev.brt" : ".brt";
    bwt.write(baseFN + extension);
    logger.logInfo("Wrote file " + baseFN + extension);
}

void createIndex(string& T, const BuildParameters& params,
                 const Alphabet<ALPHABET>& sigma,
                 const vector<length_t>& charCounts) {

    // aliasing
    const auto& baseFN = params.baseFN;

    // create SA and BWT
    {
        vector<length_t> SA;
        string BWT;
        createSuffixArrayAndBWTAndWriteSparseSA(T, SA, BWT, params);
        SA.clear();
        vector<length_t>().swap(SA);
        writeBWT(sigma, baseFN, BWT);

        // create succinct BWT bitvector table
        writeBWTBitvectors(baseFN, BWT, sigma, false);

    } // string BWT goes out of scope here
    logger.logInfo("Switching to reversed text...");

    // read or create the reverse suffix array
    vector<length_t> revSA;
    createRevSAWithSanityCheck(revSA, T);

    // build the reverse BWT
    string rBWT;
    createRevBWT(baseFN, T, revSA, sigma, rBWT);
    revSA.clear();
    vector<length_t>().swap(revSA);
    writeBWTBitvectors(baseFN, rBWT, sigma, true);
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

#ifdef RUN_LENGTH_COMPRESSION
    bool noWriting = true;
#else
    bool noWriting = false;
#endif

    preprocessFastaFiles(params.fastaFiles, params.baseFN, T, params.seedLength,
                         noWriting);

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

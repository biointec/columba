/******************************************************************************
 *  Columba: Approximate Pattern Matching using Search Schemes                *
 *  Copyright (C) 2020-2022 - Luca Renders <luca.renders@ugent.be> and        *
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
#include "definitions.h"
#include "logger.h"
#include "parameters/alignparameters.h"
#include "reads.h"
#include "searchstrategy.h"

#ifdef RUN_LENGTH_COMPRESSION
#include "bmove/bmove.h" // for BMove
#else
#include "fmindex/fmindex.h" // for FMIndex
#endif

#include <algorithm>
#include <chrono>
#include <set>
#include <sstream> // used for splitting strings
#include <string.h>

using namespace std;

string getFileExt(const string& s) {

    size_t i = s.rfind('.', s.length());
    if (i != string::npos) {
        return (s.substr(i + 1, s.length() - i));
    }

    return ("");
}

bool containsNonACGT(const std::string& str) {
    for (char c : str) {
        // Check if the character is not one of A, C, G, or T
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
            return true; // Found a non-ACGT character
        }
    }
    return false; // All characters are A, C, G, or T
}

vector<ReadBundle> getReads(const string& file) {
    vector<ReadBundle> reads;
    reads.reserve(200000);

    const auto& extension = getFileExt(file);

    bool fasta =
        (extension == "FASTA") || (extension == "fasta") || (extension == "fa");
    bool fastq = (extension == "fq") || (extension == "fastq");

    ifstream iFile(file.c_str());
    if (!iFile) {
        throw runtime_error("Cannot open file " + file);
    }
    if (!fasta && !fastq) {
        // this is a not readable

        throw runtime_error("extension " + extension +
                            " is not a valid extension for the reads file");
    } else if (fasta) {
        // fasta file
        string read = "";
        string id = "";
        string qual = ""; // empty quality string for fasta
        string line;

        while (getline(iFile, line)) {
            if (line.empty()) {
                continue; // Skip empty lines
            }

            if (line[0] == '>' || line[0] == '@') {
                // This is an ID line
                if (!id.empty()) {
                    // If we already have data, process it and clear
                    reads.emplace_back(ReadBundle(Read(id, read, qual)));

                    id.clear();
                    read.clear();
                }
            } else {
                // This is a sequence line
                read += line;
            }
        }

        // Process the last entry if it exists
        if (!id.empty()) {
            reads.emplace_back(ReadBundle(Read(id, read, qual)));
        }
    } else {
        // fastQ
        string read = "";
        string id = "";
        string qual = "";
        string plusLine = ""; // Skip the '+' line
        string line;

        while (getline(iFile, id) && getline(iFile, read) &&
               getline(iFile, plusLine) && // Skip the '+' line
               getline(iFile, qual)) {
            if (!id.empty() && id[0] != '@') {
                throw runtime_error("File " + file +
                                    "doesn't appear to be in FastQ format");
            }

            if (id.back() == '\n') {
                id.pop_back();
            }

            if (!read.empty() && read.back() == '\n') {
                read.pop_back();
            }

            assert(id.size() > 1);

            reads.emplace_back(ReadBundle(Read(id, read, qual)));

            id.clear(), read.clear(), qual.clear();
        }
    }

    return reads;
}

void writeToOutput(const string& file, const vector<vector<TextOcc>>& mPerRead,
                   const std::string& baseFile,
                   const std::string& commandLineParameters) {
    stringstream ss;

    ss << "Writing to output file " << file << " ...";
    logger.logInfo(ss);
    ofstream f2;
    f2.open(file);

    f2 << "@HD\tVN:1.6"
       << "\t"
       << "SO:queryname"
       << "\n";

    f2 << "@PG\tID:Columba" << to_string(VERSION_NUMBER_COLUMBA) << "."
       << to_string(SUB_VERSION_NUMBER_COLUMBA)
       << "\tPN:Columba\tCL:" << commandLineParameters << "\n";

    ifstream ifs(baseFile + ".headerSN.bin", std::ios::binary);

    if (!ifs) {
        logger.logError("Unable to open the header file.");
    } else {
        std::string line;
        f2 << line << "\n";
    }

    ifs.close();

    for (unsigned int i = 0; i < mPerRead.size(); i++) {

        for (auto m : mPerRead[i]) {
            f2 << m.getOutputLine() << "\n";
        }
    }

    f2.close();
}

double findMedian(vector<length_t> a, int n) {

    // If size of the arr[] is even
    if (n % 2 == 0) {

        // Applying nth_element
        // on n/2th index
        nth_element(a.begin(), a.begin() + n / 2, a.end());

        // Applying nth_element
        // on (n-1)/2 th index
        nth_element(a.begin(), a.begin() + (n - 1) / 2, a.end());

        // Find the average of value at
        // index N/2 and (N-1)/2
        return (double)(a[(n - 1) / 2] + a[n / 2]) / 2.0;
    }

    // If size of the arr[] is odd
    else {

        // Applying nth_element
        // on n/2
        nth_element(a.begin(), a.begin() + n / 2, a.end());

        // Value at index (N/2)th
        // is the median
        return (double)a[n / 2];
    }
}

void doBench(vector<ReadBundle>& reads, IndexInterface& mapper,
             std::unique_ptr<SearchStrategy>& strategy,
             const string& outputFile, length_t distanceOrIdentity,
             const std::string& baseFile,
             const std::string& commandLineParameters) {

    stringstream ss;

    const auto& mMode = strategy->getMappingModeString();

    ss << "Benchmarking with " << strategy->getName() << " strategy\n";
    ss << "\tMapping mode: " << mMode << "\n";
    if (strategy->getMappingModeString() == "BEST") {
        ss << "\tMin identity: " << distanceOrIdentity << "\n";
    } else {
        ss << "\tMax distance: " << distanceOrIdentity << "\n";
    }
    ss << "\tPartitioning strategy: " << strategy->getPartitioningStrategy()
       << "\n";
    ss << "\tDistance metric: " << strategy->getDistanceMetric() << "\n";
#ifndef RUN_LENGTH_COMPRESSION
    ss << "\tIn-text verification switch point: " << strategy->getSwitchPoint();
#endif

    logger.logInfo(ss);

    ss.precision(2);

    vector<vector<TextOcc>> matchesPerRead = {};
    matchesPerRead.reserve(reads.size());

    std::vector<length_t> numberMatchesPerRead;
    numberMatchesPerRead.reserve(reads.size());

    Counters counters;

    auto start = chrono::high_resolution_clock::now();
    for (unsigned int i = 0; i < reads.size(); i += 1) {

        auto& bundle = reads[i];

        if (i % 1024 == 0) {
            // Progress printed to cout
            // is this necessary??
            cout << "Progress: " << i << "/" << reads.size() << "\r";
            cout.flush();
        }

        counters.inc(Counters::NUMBER_OF_READS);

        Counters recordCounters;
        std::vector<TextOcc> matches;

        strategy->matchApprox(bundle, distanceOrIdentity, recordCounters,
                              matches);

        bool mapped = !matches.empty() && matches.front().isValid();

        counters.inc(Counters::TOTAL_UNIQUE_MATCHES,
                     recordCounters.get(Counters::DROPPED_UNIQUE_MATCHES) +
                         ((mapped) ? (matches.size()) : 0));

        // keep track of the number of mapped reads
        counters.inc(Counters::MAPPED_READS, mapped);

        counters.addCounters(recordCounters);

        matchesPerRead.emplace_back(matches);
        numberMatchesPerRead.emplace_back(matches.size());
    }

    auto finish = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = finish - start;
    ss << "Progress: " << reads.size() << "/" << reads.size();
    logger.logInfo(ss);
    ss << "Results for " << strategy->getName();
    logger.logInfo(ss);

    ss << "Total duration: " << fixed << elapsed.count() << "s";
    logger.logInfo(ss);
    counters.reportStatistics(SINGLE_END);
    writeToOutput(outputFile, matchesPerRead, baseFile, commandLineParameters);
}

int main(int argc, char* argv[]) {

    if (argc == 2 && strcmp("help", argv[1]) == 0) {
        Parameters::printHelp();
        return EXIT_SUCCESS;
    }

    int requiredArguments = 4; // baseFile of files and file containing reads

    if (argc < requiredArguments) {
        logger.logError("Insufficient number of arguments\n");
        Parameters::printHelp();
        return EXIT_FAILURE;
    }

    Parameters params = Parameters::processOptionalArguments(argc, argv);

    if (params.logFile != "") {
        logger.setLogFile(params.logFile);
    }

    logger.logInfo("Welcome to Columba!\n");
    logger.logInfo("Using " + string(LENGTH_TYPE_NAME) + "\n");

#ifndef RUN_LENGTH_COMPRESSION
    FMIndex index(params.base, params.inTextVerificationPoint, false,
                  params.sparsenessFactor, true, params.kmerSize);
#else
    BMove index(params.base, true, params.kmerSize);
#endif
    logger.logInfo("Reading in reads from " + params.firstReadsFile + "\n");

    vector<ReadBundle> reads;
    try {
        reads = getReads(params.firstReadsFile);
    } catch (const exception& e) {
        string er = e.what();
        er += " Did you provide a valid reads file?";
        throw runtime_error(er);
    }

    auto strategy = params.createStrategy(index);

    doBench(reads, index, strategy, params.outputFile,
            params.getMaxOrIdentity(), params.base, params.command);

    logger.logInfo("Bye...!");
    return 0;
}
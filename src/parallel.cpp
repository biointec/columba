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

#include "definitions.h" // for length_t, Orientation, FF, FR, RF, ALL
#include "fastq.h"       // for OutputChunk, OutputRecord, Reader, Out...
#include "indexinterface.h"
#ifdef RUN_LENGTH_COMPRESSION
#include "bmove/bmove.h" // for BMove
#else
#include "fmindex/fmindex.h" // for FMIndex
#endif
#include "indexhelpers.h"               // for TextOcc, Counters
#include "logger.h"                     // for Logger, logger
#include "parameters/alignparameters.h" // for Parameters
#include "reads.h"                      // for ReadPair, Read, ReadBundle
#include "searchstrategy.h"             // for SearchStrategy

#include <algorithm>  // for max, copy, for_each, count, sort
#include <assert.h>   // for assert
#include <chrono>     // for duration, operator-, high_resolution_c...
#include <cmath>      // for sqrt, abs
#include <exception>  // for exception
#include <functional> // for ref, mem_fn, _Mem_fn
#include <iterator>   // for move_iterator, make_move_iterator
#include <memory>     // for allocator, allocator_traits<>::value_type
#include <mutex>      // for mutex, lock_guard, unique_lock
#include <numeric>    // for accumulate
#include <ostream>    // for operator<<, basic_ostream, stringstream
#include <stdint.h>   // for int32_t
#include <stdlib.h>   // for abs, EXIT_SUCCESS, EXIT_FAILURE
#include <string>     // for char_traits, operator<<, operator+
#include <thread>     // for thread
#include <utility>    // for move
#include <vector>     // for vector, _Bit_iterator, _Bit_reference

using namespace std;

//----------------------------------------------------------------------------
// Single-ended processing
//----------------------------------------------------------------------------

/**
 * Process a chunk of reads singled-ended using the given search strategy
 * @param input The input chunk of reads
 * @param output The output chunk (results will be added to this chunk)
 * @param strategy The search strategy to use
 * @param maxEDOrIdentity The maximum edit distance or minimum identity of the
 * matches
 */
void processChunk(vector<ReadBundle>& input, OutputChunk& output,
                  std::unique_ptr<SearchStrategy>& strategy,
                  const length_t& maxEDOrIdentity) {
    output.clear();
    output.reserve(input.size());
    for (auto& readBundle : input) {
        std::vector<TextOcc> matches;
        strategy->matchApprox(readBundle, maxEDOrIdentity,
                              output.getMutableCounters(), matches);
        output.addSingleEndRecord(std::move(matches));
    }
}

/**
 * Entry for a worker thread for processing reads single-ended.
 * @param myReader The reader to read the input from
 * @param myWriter The writer to write the output to
 * @param strategy The search strategy to use
 * @param maxEDOrIdentity The maximum edit distance or minimum identity of
 * the matches
 */
void threadEntrySingleEnd(Reader& myReader, OutputWriter& myWriter,
                          std::unique_ptr<SearchStrategy>& strategy,
                          const length_t& maxEDOrIdentity) {
    // local storage of reads
    vector<ReadBundle> input;
    OutputChunk output;

    size_t chunkID;
    while (myReader.getNextChunk(input, chunkID)) {

        auto start = chrono::high_resolution_clock::now();

        processChunk(input, output, strategy, maxEDOrIdentity);
        // mark output as ready
        output.set(chunkID);
        // add output to writer
        myWriter.commitChunk(output);

        auto end = chrono::high_resolution_clock::now();
        auto elapsed =
            std::chrono::duration_cast<std::chrono::microseconds>(end - start);

        // Record processing time
        myReader.addProcessingTime(elapsed);
    }

    myWriter.sendTermination(chunkID);
}

//----------------------------------------------------------------------------
// Paired-ended processing
//----------------------------------------------------------------------------

/**
 * Process a chunk of reads paired-ended using the given search strategy
 * @param input The input chunk of reads
 * @param output The output chunk (results will be added to this chunk)
 * @param strategy The search strategy to use
 * @param maxEDOrIdentity The maximum edit distance or minimum identity of
 the
 * matches
 * @param maxFragSize The maximum fragment size
 * @param minFragSize The minimum fragment size

 */
void processChunkPairedEnd(const vector<ReadBundle>& input, OutputChunk& output,
                           std::unique_ptr<SearchStrategy>& strategy,
                           const length_t& maxEDOrIdentity,
                           const length_t& maxFragSize,
                           const length_t& minFragSize) {
    output.clear();

    // In case of Paired-end processing we need to take two reads at a time
    // from the input vector. We can do this by iterating over the input
    // vector with a step size of 2.
    assert(input.size() % 2 == 0);
    output.reserve(input.size() / 2);

    for (length_t i = 0; i < input.size(); i += 2) {
        ReadPair readPair(input[i], input[i + 1]);

        std::vector<TextOcc> unpairedOccs;
        auto pairedMatches = strategy->matchApproxPE(
            readPair, maxEDOrIdentity, maxFragSize, minFragSize,
            output.getMutableCounters(), unpairedOccs);

        output.addPairedEndRecord(std::move(pairedMatches),
                                  std::move(unpairedOccs));
    }
}

/**
 * Entry for a worker thread for processing reads paired-ended.
 * @param myReader The reader to read the input from
 * @param myWriter The writer to write the output to
 * @param strategy The search strategy to use
 * @param maxEDOrIdentity The maximum edit distance or minimum identity of
 * the matches
 * @param maxFragSize The maximum fragment size
 * @param minFragSize The minimum fragment size
 */
void threadEntryPairedEnd(Reader& myReader, OutputWriter& myWriter,
                          std::unique_ptr<SearchStrategy>& strategy,
                          const length_t& maxEDOrIdentity,
                          const length_t& maxFragSize,
                          const length_t& minFragSize) {
    // local storage of reads
    vector<ReadBundle> input;
    OutputChunk output;

    size_t chunkID;
    while (myReader.getNextChunk(input, chunkID)) {

        auto start = chrono::high_resolution_clock::now();
        processChunkPairedEnd(input, output, strategy, maxEDOrIdentity,
                              maxFragSize, minFragSize);

        // mark output as ready
        output.set(chunkID);
        // add output to writer
        myWriter.commitChunk(output);
        output.clear();
        auto end = chrono::high_resolution_clock::now();
        auto elapsed =
            std::chrono::duration_cast<std::chrono::microseconds>(end - start);

        // Record processing time
        myReader.addProcessingTime(elapsed);
    }

    myWriter.sendTermination(chunkID);
}

//----------------------------------------------------------------------------
// Inferring paired-end parameters
//----------------------------------------------------------------------------

/**
 * A struct that bundles the read pairs and the single-ended output chunks
 * for those reads. Also keeps a vector with booleans that indicate whether
 * the second read was processed. The second read is not processed if the
 * processed first read does not have an unambiguous match.
 */
struct InferParametersProcessedChunk {
    vector<ReadPair> readPairs; // the read pairs processed in this chunk
    OutputChunk output; // the outputChunk with the single-ended occurrences
    vector<bool> boolsRead2Done; // whether the second read was processed
                                 // for each pair

    void clear() {
        readPairs.clear();
        output.clear();
        boolsRead2Done.clear();
    }

    InferParametersProcessedChunk() {
        clear();
    }
};

/**
 * Check if there is only one match in the first file for which the best
 * distance score is found. If this is the case, swap this match with the first
 * match in the vector.
 * @param matches The vector of matches
 * @param strategy The underlying search strategy to find if the match is in the
 * first file or not
 */
bool hasUnambiguousMatchInFirstFile(
    std::vector<TextOcc>& matches,
    const std::unique_ptr<SearchStrategy>& strategy) {

    // only check if front is valid
    if (matches.empty() || !matches.front().isValid()) {
        return false;
    }

    // check if only one match is found for which
    // strategy->isInFirstFile(match) is true
    length_t count = 0;
    size_t indexOfFirstMatch = 0;
    for (size_t i = 0; i < matches.size(); ++i) {
        const auto& match = matches[i];
        if (strategy->isInFirstFile(match)) {
            if (++count == 2) {
                return false;
            } else {
                indexOfFirstMatch = i;
            }
        }
    }
    // swap the element at index indexOfFirstMatch with the first element
    std::swap(matches[0], matches[indexOfFirstMatch]);
    return count == 1;
}

/**
 * Process a chunk of reads single-ended with the purpose of inferring
 * paired-end parameters. Processed chunks will afterwards be re-processed
 * to do the actual pairing.
 * @param input The input chunk of reads
 * @param output The output chunk (results will be added to this chunk)
 * @param strategy The search strategy to use
 * @param maxEDOrIdentity The maximum edit distance or minimum identity of
 * the
 * @param identity whether maxEDOrIdentity is identity or edit distance
 * matches
 */
void processChunkSingleEndForPairInferring(
    vector<ReadBundle>& input, InferParametersProcessedChunk& output,
    std::unique_ptr<SearchStrategy>& strategy, const length_t& maxEDOrIdentity,
    bool identity) {
    output.clear();

    bool first = true;
    bool firstHasUnambiguousMatch = false;

    ReadBundle firstBundle;
    length_t maxED = maxEDOrIdentity;

    for (auto& bundle : input) {

        if (identity) {
            maxED = strategy->getMaxED(maxEDOrIdentity, bundle.size());
        }

        std::vector<TextOcc> matches;
        if (first || firstHasUnambiguousMatch) {
            strategy->matchApprox(bundle, maxED,
                                  output.output.getMutableCounters(), matches);
        }

        // Keep the info of the paired reads together
        if (first) {
            firstBundle = std::move(bundle);
            firstHasUnambiguousMatch =
                hasUnambiguousMatchInFirstFile(matches, strategy);
        } else {
            output.readPairs.emplace_back(std::move(firstBundle),
                                          std::move(bundle));
            output.boolsRead2Done.push_back(firstHasUnambiguousMatch);
        }
        // add the matches to the output record
        output.output.addSingleEndRecord(std::move(matches));
        first = !first; // switch between first and second read
    }
}

/**
 * Helper function for inferPairedEndParameters. Adds the fragment size and
 * the orientation of the two matches to the vectors fragSizes and
 * orientations.
 * @param match1 The first match
 * @param match2 The second match
 * @param fragSizes The vector of fragment sizes (this will be updated
 * during the function)
 * @param orientations The vector of orientations (this will be updated
 * during the function)
 */

void addFragmentAndOrientation(const TextOcc& match1, const TextOcc& match2,
                               vector<length_t>& fragSizes,
                               vector<Orientation>& orientations) {

    // calculate the fragment size
    bool match1First = match1.getBegin() < match2.getBegin();
    length_t fragSize = 0;

    if (match1First) {
        fragSize = match2.getEnd() - match1.getBegin();
    } else {
        fragSize = match1.getEnd() - match2.getBegin();
    }

    fragSizes.push_back(fragSize);

    // calculate the orientation
    if (match1.isRevCompl() == match2.isRevCompl()) {
        // both reads same orientation
        orientations.push_back(FF);
    } else {
        // one read maps to the reverse complement of the reference
        // if this read is the first read, the orientation is RF,
        // else FR

        if (match1First == match1.isRevCompl()) {
            orientations.push_back(RF);
        } else {
            orientations.push_back(FR);
        }
    }
}

/**
 * Helper function for inferPairedEndParameters. Calculates the median of a
 * vector of length_t.
 */
length_t calcMedian(vector<length_t>& data) {
    sort(data.begin(), data.end());
    if (data.size() % 2 == 0) {
        return (data[data.size() / 2 - 1] + data[data.size() / 2]) / 2;
    } else {
        return data[data.size() / 2];
    }
}

/**
 * Helper function for inferPairedEndParameters. Calculates the average of a
 * vector of length_t.
 */
float average(std::vector<length_t>& v) {
    return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}

/**
 * Helper function for inferPairedEndParameters. Calculates the standard
 * deviation of a vector of length_t. With Bessel's correction.
 * @param data The vector of length_t
 * @param mean The mean of the data
 */
float stddev(std::vector<length_t>& data, float mean) {
    float accum = 0.0;
    for_each(data.begin(), data.end(),
             [&](const float d) { accum += (d - mean) * (d - mean); });
    return sqrt(accum / (data.size() - 1));
}

/**
 * Infers the paired-end parameters (orientation, max and min insert size)
 * by looking at the fragment sizes and orientations.
 * @param params The Parameters object
 * @param fragSizes The vector of fragment sizes
 * @param orientations The vector of orientations
 */
void inferPairedEndParameters(Parameters& params, vector<length_t>& fragSizes,
                              vector<Orientation>& orientations,
                              float& meanInsertSize, float& stddevInsertSize) {
    // calculate the mean and standard deviation of the fragment sizes
    // but remove outliers via median absolute deviation
    if (fragSizes.empty()) {
        logger.logWarning("No pairs mapped unambiguously. "
                          "Using default values!");

        return;
    }

    // calculate median fragSize
    int32_t medianFragSize = static_cast<int32_t>(calcMedian(fragSizes));

    // calculate median absolute deviation
    vector<length_t> absDeviations;
    absDeviations.reserve(fragSizes.size());
    std::transform(
        fragSizes.begin(), fragSizes.end(), std::back_inserter(absDeviations),
        [medianFragSize](const length_t fragSize) {
            return abs(static_cast<int32_t>(fragSize) - medianFragSize);
        });

    length_t medianAbsDev = calcMedian(absDeviations);

    // remove outliers
    vector<length_t> filteredFragSizes;
    filteredFragSizes.reserve(fragSizes.size());
    std::copy_if(fragSizes.begin(), fragSizes.end(),
                 std::back_inserter(filteredFragSizes),
                 [medianFragSize, medianAbsDev](const length_t fragSize) {
                     return abs(static_cast<int32_t>(fragSize) -
                                medianFragSize) <
                            PE_STD_DEV_CONSIDERED * medianAbsDev;
                 });

    // calculate the mean and standard deviation of the filtered fragment
    // sizes
    meanInsertSize = average(filteredFragSizes);
    stddevInsertSize = stddev(filteredFragSizes, meanInsertSize);

    if (meanInsertSize == 0 || stddevInsertSize == 0) {
        meanInsertSize = average(fragSizes);
        stddevInsertSize = stddev(fragSizes, meanInsertSize);
    }

    // set the max and min insert size
    length_t maxDev = PE_STD_DEV_CONSIDERED * stddevInsertSize;

    params.maxInsertSize = meanInsertSize + maxDev;
    params.minInsertSize =
        (meanInsertSize > maxDev) ? meanInsertSize - maxDev : 0;

    // count the occurrences of each orientation
    length_t rfCount = count(orientations.begin(), orientations.end(), RF);
    length_t frCount = count(orientations.begin(), orientations.end(), FR);
    length_t ffCount = count(orientations.begin(), orientations.end(), FF);

    // Determine the most common orientation
    params.orientation = (frCount >= rfCount && frCount >= ffCount) ? FR
                         : (rfCount >= ffCount)                     ? RF
                                                                    : FF;
}

/**
 * Class for inferring the paired-end parameters. This class is first used
 * to give out input chunks to threads for single-ended processing. The
 * threads then indicate if they have found an unambiguously mapped pair. If
 * enough such pairs have been found the queue stops giving out input
 * chunks. The worker threads then all dump their output from trying to find
 * the unambiguous pairs into this class. Then new workers can start
 * processing the single-end output in paired-end mode by getting that
 * output from this class.
 */
class InferParametersQueue {
  private:
    Reader& myReader; // The reader to read the input from

    mutex fragSizesAndOrientationsMutex;
    vector<length_t> fragSizes;       // vector of fragment sizes
    vector<Orientation> orientations; // vector of orientations

    mutex readsGivenMutex;
    length_t readsGiven = 0; // Number of reads given to the workers

    mutex unambiguousPairsMutex;
    length_t unambiguousPairs = 0; // Number of unambiguous pairs

    mutex stopMutex;
    bool stop = false; // Stop flag

    length_t logThresholdUnambiguousPairs = PE_NUMBER_PAIRS_FOR_INFERENCE / 2;
    length_t logThresholdReadsGiven = PE_MAX_READS_FOR_INFERENCE / 2;

    mutex numLogsMutex;
    length_t numLogsPairs = 0;
    length_t numLogsReads = 0;

    mutex processedChunksMutex;
    // The single-ended processed chunks
    vector<InferParametersProcessedChunk> processedChunks;
    length_t nextOutputToGive = 0; // the next outputChunk to give to a worker

    void logMessage() const {
        stringstream s;
        s << "Found " << unambiguousPairs
          << " unambiguous pairs while processing " << readsGiven
          << " reads for inferring paired-end parameters";
        logger.logInfo(s);
    }

    void log() {
        // acquire the numLogsMutex
        lock_guard<mutex> lock(numLogsMutex);
        if (unambiguousPairs >=
            logThresholdUnambiguousPairs * (numLogsPairs + 1)) {
            logMessage();
            numLogsPairs++;
        } else if (readsGiven >= logThresholdReadsGiven * (numLogsReads + 1)) {
            logMessage();
            numLogsReads++;
        }
    }

  public:
    InferParametersQueue(Reader& myReader) : myReader(myReader) {
    }
    /**
     * Thread-safe increment of the number of unambiguous pairs. If enough
     * such pairs have been found the stop flag is set to true.
     */
    void unambiguousPairsIncrement() {
        lock_guard<mutex> lock(unambiguousPairsMutex);
        unambiguousPairs++;
        log();
        if (unambiguousPairs >= PE_NUMBER_PAIRS_FOR_INFERENCE) {
            lock_guard<mutex> lockStop(stopMutex);
            stop = true;
        }
    }

    /**
     * Add processing time to the reader, such that the input chunk size can
     * change dynamically
     */
    void addProcessingTime(std::chrono::microseconds time) {
        myReader.addProcessingTime(time);
    }

    /**
     * Get the next input chunk. If enough unambiguous pairs have been found
     * or if too many reads have already been given for inference, false is
     * returned. If no more chunks are available, false is returned.
     * @param input The input chunk (output parameter)
     * @param chunkID The ID of the chunk (output parameter)
     * @return True if a chunk was given, false if no more chunks should be
     * processed of if no chunk exists
     */
    bool getNextChunk(vector<ReadBundle>& input, size_t& chunkID) {
        unique_lock<mutex> lock(stopMutex);
        if (stop) {
            return false;
        }
        lock.unlock(); // release stopMutex
        if (!myReader.getNextChunk(input, chunkID)) {
            return false;
        }
        lock_guard<mutex> lockReadsGiven(readsGivenMutex);
        log();
        readsGiven += input.size();

        if (readsGiven >= PE_MAX_READS_FOR_INFERENCE) {
            lock_guard<mutex> lockStop(stopMutex);
            stop = true;
        } // releases stopMutex
        return true;
    }

    /**
     * Thread-safely add processedChunks for inference to this
     * queue. Should be called after getNextChunk returns false.
     * @param processedChunks The vector of processedChunks. Will be made
     * invalid.
     */
    void addProcessedChunks(
        vector<InferParametersProcessedChunk>&& processedChunks) {
        lock_guard<mutex> lock(processedChunksMutex);
        // Use move semantics to avoid copying
        this->processedChunks.insert(
            this->processedChunks.end(),
            make_move_iterator(processedChunks.begin()),
            make_move_iterator(processedChunks.end()));
    }

    /**
     * Thread-safely add the fragment sizes and orientations to this queue.
     * Should be called after getNextChunk returns false.
     * @param fragSizes The vector of fragment sizes. Will be made invalid.
     * @param orientations The vector of orientations Will be made invalid.
     */
    void addFragmentsAndOrientations(vector<length_t>&& fragSizes,
                                     vector<Orientation>&& orientations) {
        assert(fragSizes.size() == orientations.size());
        lock_guard<mutex> lock(fragSizesAndOrientationsMutex);
        // Use move semantics to avoid copying
        this->fragSizes.insert(this->fragSizes.end(),
                               make_move_iterator(fragSizes.begin()),
                               make_move_iterator(fragSizes.end()));
        this->orientations.insert(this->orientations.end(),
                                  make_move_iterator(orientations.begin()),
                                  make_move_iterator(orientations.end()));
    }

    /**
     * Infer the paired-end parameters (orientation, max and min insert
     * size) and set it in the Parameters object. Should be called after all
     * worker threads have called addFragmentsAndOrientations.
     * @param params The Parameters object
     * @param meanInsertSize The mean insert size (output parameter)
     * @param stddevInsertSize The standard deviation of the insert size (output
     * parameter)
     *
     */
    void inferParameters(Parameters& params, float& meanInsertSize,
                         float& stddevInsertSize) {
        lock_guard<mutex> lock(fragSizesAndOrientationsMutex);
        inferPairedEndParameters(params, fragSizes, orientations,
                                 meanInsertSize, stddevInsertSize);
    }

    /**
     * Thread-safely gets the next InferParametersProcessedChunk to process
     * paired-end. If no more output chunks are available false is returned.
     * The very first call to this function ensures that the processedChunks
     * are sorted based on the chunkID of the encapsulated outputChunks.
     * Such that the chunks are processed in the same order.
     * @param output The output chunk to process in paired-end mode (output
     * parameter)
     */
    bool getNextOutput(InferParametersProcessedChunk& output) {

        lock_guard<mutex> lock(processedChunksMutex);

        if (nextOutputToGive >= processedChunks.size()) {
            return false;
        }
        output = std::move(processedChunks[nextOutputToGive]);
        nextOutputToGive++;
        return true;
    }
};

/**
 * Class for a worker thread that processes reads single-ended for inferring
 * the paired-end parameters. It keeps track of their own vectors of info
 * about reads, outputChunks from single-end processing, fragment sizes and
 * orientations.
 */
class InferParametersWorker {
  private:
    // the processed chunks for inferring the parameters
    vector<InferParametersProcessedChunk> processedChunks = {};
    vector<length_t> fragSizes = {};       // vector of fragment sizes
    vector<Orientation> orientations = {}; // vector of orientations

    std::unique_ptr<SearchStrategy>& strategy; // The search strategy to use
    length_t maxEDOrIdentity; // The maximum edit distance or minimum
                              // identity of the matches
    bool identity; // whether maxEDOrIdentity is identity or  distance

    std::thread thread; // The worker thread

    /**
     * The worker thread function. Processes the input chunks single-ended
     * and if a pair is mapped unambiguously, adds the fragment size and
     * orientation to the vectors as well as increments the number of
     * unambiguously mapped pairs in the queue. If no more input chunks need
     * to be processed. This thread will dump all its local info to the
     * queue such that parameters can be inferred.
     * @param queue The queue to get input chunks from and to dump the local
     * output to
     */
    void workerThread(InferParametersQueue& queue) {

        // local storage of reads
        vector<ReadBundle> input;
        size_t chunkID;

        while (queue.getNextChunk(input, chunkID)) {
            // Assert that these reads are paired end
            assert(input.size() % 2 == 0);

            // process the chunk singled ended
            auto startTime = chrono::high_resolution_clock::now();
            InferParametersProcessedChunk processedChunk;
            processChunkSingleEndForPairInferring(
                input, processedChunk, strategy, maxEDOrIdentity, identity);

            auto end = chrono::high_resolution_clock::now();
            auto elapsed =
                std::chrono::duration_cast<std::chrono::microseconds>(
                    end - startTime);

            // Record processing time
            queue.addProcessingTime(elapsed);

            // iterate over the output records and add the fragment sizes if
            // unique pair
            const auto& records = processedChunk.output.getRecords();
            for (length_t i = 0; i < records.size(); i += 2) {
                const auto& outRecord1 = records[i];
                const auto& outRecord2 = records[i + 1];

                // we consider the pairs for which both reads have an
                // unambiguous match in the first file
                if (processedChunk.boolsRead2Done[i / 2] &&
                    hasUnambiguousMatchInFirstFile(*outRecord2.outputOcc,
                                                   strategy)) {
                    const auto& match1 = outRecord1.outputOcc->front();
                    const auto& match2 = outRecord2.outputOcc->front();

                    addFragmentAndOrientation(match1, match2, fragSizes,
                                              orientations);
                    queue.unambiguousPairsIncrement();
                }
            }

            // mark the output as ready
            processedChunk.output.set(chunkID);
            processedChunks.emplace_back(std::move(processedChunk));
        }

        // termination signal received, add the processed chunks and
        // fragSizes and Orientations to the queue
        queue.addProcessedChunks(std::move(processedChunks));
        queue.addFragmentsAndOrientations(std::move(fragSizes),
                                          std::move(orientations));
    }

  public:
    InferParametersWorker(std::unique_ptr<SearchStrategy>& strategy,
                          length_t maxEDOrIdentity, bool identity)
        : strategy(strategy), maxEDOrIdentity(maxEDOrIdentity),
          identity(identity) {
    }

    void start(InferParametersQueue& queue) {
        thread =
            std::thread(&InferParametersWorker::workerThread, this, ref(queue));
    }

    void join() {
        thread.join();
    }
};

/**
 * Entry for a worker thread for processing reads paired-ended after
 * inferring the parameters. First process the chunks that were used for
 * inferring the parameters, then process the remaining chunks.
 * @param myReader The reader to read the input from
 * @param myWriter The writer to write the output to
 * @param strategy The search strategy to use
 * @param maxEDOrIdentity The maximum edit distance or minimum identity of
 * the matches
 * @param maxFragSize The maximum fragment size
 * @param minFragSize The minimum fragment size
 * @param q The queue that contains the read pairs and output chunks that
 * were used for inferring the paired end parameters. Must stay in scope
 * until all accessing threads have finished.
 */
void threadEntryPairedEndAfterInferring(
    Reader& myReader, OutputWriter& myWriter,
    std::unique_ptr<SearchStrategy>& strategy, const length_t& maxEDOrIdentity,
    const length_t& maxFragSize, const length_t& minFragSize,
    InferParametersQueue& q) {

    InferParametersProcessedChunk
        inferChunkCurrent; // The current chunk to process paired-end

    while (q.getNextOutput(inferChunkCurrent)) {

        OutputChunk output; // The output chunk to be written
        output.clear();
        output.set(inferChunkCurrent.output.getChunkID()); // copy the ID

        // combine the counters of all records
        output.getMutableCounters().addCounters(
            inferChunkCurrent.output.getMutableCounters());

        // Pair the records and generate SAM lines
        auto& records = inferChunkCurrent.output.getRecordsNonConst();
        for (length_t i = 0; i < records.size(); i += 2) {
            auto& matches1 = records[i].outputOcc;     // matches for read 1
            auto& matches2 = records[i + 1].outputOcc; // matches for read 2
            auto& readPair =
                inferChunkCurrent.readPairs[i / 2]; // The read pair with info
            // Whether read 2 has been processed
            auto read2done = inferChunkCurrent.boolsRead2Done[i / 2];

            std::vector<TextOcc> unpairedOccs;
            std::vector<PairedTextOccs> pairedMatches = strategy->pairSingleEnd(
                readPair, *matches1, *matches2, maxFragSize, minFragSize,
                maxEDOrIdentity, unpairedOccs, output.getMutableCounters(),
                read2done);

            // Add the paired matches and unpaired to the output chunk
            output.addPairedEndRecord(std::move(pairedMatches),
                                      std::move(unpairedOccs));
        }
        // send chunk to writer
        myWriter.commitChunk(output);
    }
    // After all chunks used for inferring have been processed, go to
    // threadEntryPairedEnd to process the remaining input chunks
    threadEntryPairedEnd(myReader, myWriter, strategy, maxEDOrIdentity,
                         maxFragSize, minFragSize);
}

/**
 * @brief Logs paired-end parameters.
 *
 * @param params Parameters object containing necessary parameters.

 */
void logPairedEndParameters(const Parameters& params) {
    std::stringstream ss;
    ss << "Command-line paired-end parameters:\n";
    ss << "\tMax fragment size: " << params.maxInsertSize << "\n";
    ss << "\tMin fragment size: " << params.minInsertSize << "\n";
    ss << "\tOrientation: "
       << ((params.orientation == FF)
               ? "FF"
               : ((params.orientation == FR) ? "FR" : "RF"));
    logger.logInfo(ss.str());
}

/**
 * @brief Logs paired-end parameters.
 *
 * @param params Parameters object containing necessary parameters.
 * @param meanInsertSize Mean detected insert size.
 * @param stddevInsertSize Standard deviation of mean detected insert size.
 */
void logInferredPairedEndParameters(const Parameters& params,
                                    float meanInsertSize,
                                    float stddevInsertSize) {
    std::stringstream ss;
    ss << "Inferred paired-end parameters:\n";
    ss << "\tDetected Mean fragment size: " << meanInsertSize << "\n";
    ss << "\tStandard deviation of detected fragment size: " << stddevInsertSize
       << "\n";
    ss << "\tAllowing " << PE_STD_DEV_CONSIDERED
       << " standard deviations from the mean\n";
    ss << "\tMax fragment size: " << params.maxInsertSize << "\n";
    ss << "\tMin fragment size: " << params.minInsertSize << "\n";
    ss << "\tOrientation: "
       << ((params.orientation == FF)
               ? "FF"
               : ((params.orientation == FR) ? "FR" : "RF"));
    logger.logInfo(ss.str());
}

/**
 * Infer the paired-end parameters (orientation, max and min insert size) by
 * processing a number of read pairs single-ended. The parameters are set in
 the
 * Parameters object. Afterwards the paired-end processing is started.

 * @param myReader The  Reader to read the input from
 * @param myWriter The writer to write the paired-end output to
 * @param strategy The search strategy to use
 * @param maxEDOrIdentity The maximum edit distance or minimum identity of
 the
 * matches
 * @param params The Parameters object
 */
void inferParametersAndStartWorkers(Reader& myReader, OutputWriter& myWriter,
                                    std::unique_ptr<SearchStrategy>& strategy,
                                    length_t maxEDOrIdentity,
                                    Parameters& params) {
    // Create the queue
    InferParametersQueue queue(myReader);

    // set the strategy to not generate SAM lines and in all mapping mode
    strategy->setGenerateSAMLines(false);
    strategy->setMappingMode(ALL);
    bool identity = params.mMode != ALL;

    // Create the workers and start them
    vector<InferParametersWorker> workers;
    for (size_t i = 0; i < params.nThreads; i++) {
        workers.emplace_back(strategy, maxEDOrIdentity, identity);
    }
    for (size_t i = 0; i < params.nThreads; i++) {
        workers[i].start(queue);
    }

    // Wait for all workers to finish
    for_each(workers.begin(), workers.end(),
             [](InferParametersWorker& w) { w.join(); });

    // Infer the parameters
    float meanInsertSize = 0;
    float stddevInsertSize = 0;
    queue.inferParameters(params, meanInsertSize, stddevInsertSize);

    // set strategy correct again and with inferred orientation
    strategy->setGenerateSAMLines(true);
    strategy->setMappingMode(params.mMode);
    strategy->setOrientation(params.orientation);

    // Start the workers for the paired-end processing after inferring
    vector<thread> afterInferringWorkers(params.nThreads);
    for (size_t i = 0; i < afterInferringWorkers.size(); i++) {
        afterInferringWorkers[i] = thread(
            threadEntryPairedEndAfterInferring, ref(myReader), ref(myWriter),
            ref(strategy), ref(maxEDOrIdentity), ref(params.maxInsertSize),
            ref(params.minInsertSize), ref(queue));
    }

    // report paired-end parameters
    logInferredPairedEndParameters(params, meanInsertSize, stddevInsertSize);

    // wait for worker threads to finish
    for_each(afterInferringWorkers.begin(), afterInferringWorkers.end(),
             mem_fn(&thread::join));
}

//----------------------------------------------------------------------------
// Logging
//----------------------------------------------------------------------------

/**
 * Log the benchmarking parameters
 * @param params The Parameters object
 * @param index The index
 * @param strategy The search strategy
 */
void logAlignmentParameters(const Parameters& params,
                            const IndexInterface& index,
                            const std::unique_ptr<SearchStrategy>& strategy) {
    stringstream ss;
    ss << "Aligning with " << strategy->getName() << " strategy\n";
    ss << "\tMapping mode: " << strategy->getMappingModeString() << "\n";
    if (strategy->getMappingModeString() == "BEST") {
        ss << "\tStrata after best: " << params.strataAfterBest << "\n";
        ss << "\tMin identity: " << params.minIdentity << "\n";
    } else {
        ss << "\tMax distance: " << params.maxDistance << "\n";
    }
    ss << "\tPartitioning strategy: " << strategy->getPartitioningStrategy()
       << "\n";
    ss << "\tDistance metric: " << strategy->getDistanceMetric() << "\n";
    ss << "\tNumber of threads: " << params.nThreads << "\n";
    ss << "\tSequencing mode: "
       << ((params.sMode == SINGLE_END) ? "Single-end" : "Paired-end") << "\n";
    ss << "\tReads file: " << params.firstReadsFile << "\n";
    ss << "\tOutput file: " << params.outputFile << "\n";
#ifndef RUN_LENGTH_COMPRESSION
    ss << "\tIn-text verification switch point: " << index.getSwitchPoint()
       << "\n";
    ss << "\tCIGAR string calculation: " << (params.noCIGAR ? "No" : "Yes")
       << "\n";
#endif
    ss << "\tReorder output: " << (params.reorder ? "Yes" : "No") << "\n";
    if (params.sMode == PAIRED_END) {
        ss << "\tSecond reads file: " << params.secondReadsFile << "\n";
        ss << "\tDiscordant allowed: "
           << (params.discordantAllowed ? "Yes" : "No") << "\n";
        if (params.discordantAllowed)
            ss << "\tMax number of discordant pairs: " << params.nDiscordant
               << "\n";
        ss << "\tInfer paired-end parameters: "
           << (params.inferPairedEndParameters ? "Yes" : "No");
    } else {
        // single end
        ss << "\tXA tag: " << (params.XATag ? "Yes" : "No") << "\n";
    }

    logger.logInfo(ss);
}

/**
 * @brief Processes command-line arguments.
 *
 * @param argc Argument count.
 * @param argv Argument vector.
 * @param params Reference to Parameters object to store parsed arguments.
 * @return true if arguments are processed successfully.
 * @return false if there is an error in processing arguments (error will be
 * logged).
 */
bool processArguments(int argc, char* argv[], Parameters& params) {
    try {
        params = Parameters::processOptionalArguments(argc, argv);
    } catch (const std::exception& e) {
        logger.logError("Problem with processing parameters: " +
                        std::string(e.what()));
        Parameters::printHelp();
        return false;
    }

    if (!params.logFile.empty()) {
        logger.setLogFile(params.logFile);
    }

    return true;
}

namespace std {
/**
 * @brief Provides a `make_unique` function for C++11.
 *
 * @tparam T The type of the object to create.
 * @tparam Args The types of the arguments to pass to the constructor.
 * @param args The arguments to pass to the constructor.
 * @return std::unique_ptr<T> A unique pointer to the created object.
 */
template <typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args) {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}
} // namespace std

#ifndef RUN_LENGTH_COMPRESSION
/**
 * @brief Creates an FMIndex object.
 *
 * @param params Parameters object containing necessary parameters.
 * @return std::unique_ptr<FMIndex> Pointer to the created FMIndex object, or
 * nullptr on failure (error will be logged).
 */
std::unique_ptr<FMIndex> createFMIndex(const Parameters& params) {
    try {
        return std::make_unique<FMIndex>(
            params.base, params.inTextVerificationPoint, params.noCIGAR,
            params.sparsenessFactor, true, params.kmerSize);
    } catch (const std::exception& e) {
        logger.logError("Problem with index loading: " + std::string(e.what()));
        return nullptr;
    }
}
#else

/**
 * @brief Creates a BMove object.
 *
 * @param params Parameters object containing necessary parameters.
 * @return std::unique_ptr<BMove> Pointer to the created BMove object, or
 * nullptr on failure (error will be logged).
 */
std::unique_ptr<BMove> createBMove(const Parameters& params) {
    try {
        return std::make_unique<BMove>(params.base, true, params.kmerSize);
    } catch (const std::exception& e) {
        logger.logError("Problem with index loading: " + std::string(e.what()));
        return nullptr;
    }
}
#endif // not RUN_LENGTH_COMPRESSION

/**
 * @brief Creates a SearchStrategy object.
 *
 * @param params Parameters object containing necessary parameters.
 * @param index Reference to IndexInterface object.
 * @return std::unique_ptr<SearchStrategy> Pointer to the created SearchStrategy
 * object, or nullptr on failure (error will be logged).
 */
std::unique_ptr<SearchStrategy> createSearchStrategy(const Parameters& params,
                                                     IndexInterface& index) {
    try {
        return params.createStrategy(index);
    } catch (const std::exception& e) {
        logger.logError("Problem with search scheme loading: " +
                        std::string(e.what()));
        return nullptr;
    }
}

/**
 * @brief Creates a Reader object.
 *
 * @param params Parameters object containing necessary parameters.
 * @return std::unique_ptr<Reader> Pointer to the created Reader object, or
 * nullptr on failure (error will be logged).
 */
std::unique_ptr<Reader> createReader(const Parameters& params) {
    try {
        auto reader = std::make_unique<Reader>(params.firstReadsFile,
                                               params.secondReadsFile);
        size_t targetChunk = 64 * 1;
        reader->startReaderThread(targetChunk, params.nThreads);
        return reader;
    } catch (const std::exception& e) {
        logger.logError("Problem with reading the reads file(s): " +
                        std::string(e.what()));
        return nullptr;
    }
}

/**
 * @brief Creates an OutputWriter object.
 *
 * @param params Parameters object containing necessary parameters.
 * @return std::unique_ptr<OutputWriter> Pointer to the created OutputWriter
 * object, or nullptr on failure (error will be logged).
 */
std::unique_ptr<OutputWriter> createOutputWriter(const Parameters& params) {
    try {
        auto writer = std::make_unique<OutputWriter>(
            params.outputFile, params.base + ".headerSN.bin", params.command,
            params.sMode, params.reorder);
        writer->start(max<int>(5, (params.nThreads * 3) / 2), params.nThreads);
        return writer;
    } catch (const std::exception& e) {
        logger.logError("Problem with writing to the output file: " +
                        std::string(e.what()));
        return nullptr;
    }
}

/**
 * @brief Starts worker threads for single-end mode and waits on them
 *
 * @param workers Vector of worker threads.
 * @param reader Reference to Reader object.
 * @param writer Reference to OutputWriter object.
 * @param strategy Pointer to SearchStrategy object.
 * @param maxOrIdentity Reference to length_t for max or identity parameter.
 */
void startSingleEndWorkers(std::vector<std::thread>& workers, Reader& reader,
                           OutputWriter& writer,
                           std::unique_ptr<SearchStrategy>& strategy,
                           length_t& maxOrIdentity) {
    for (size_t i = 0; i < workers.size(); i++) {
        workers[i] = std::thread(threadEntrySingleEnd, std::ref(reader),
                                 std::ref(writer), std::ref(strategy),
                                 std::ref(maxOrIdentity));
    }
    for (auto& worker : workers) {
        worker.join();
    }
}

/**
 * @brief Starts worker threads for paired-end mode and waits for them to
 * finish.
 *
 * @param params Parameters object containing necessary parameters.
 * @param workers Vector of worker threads.
 * @param reader Reference to Reader object.
 * @param writer Reference to OutputWriter object.
 * @param strategy Pointer to SearchStrategy object.
 * @param maxOrIdentity Reference to length_t for max or identity parameter.
 */
void startPairedEndWorkers(Parameters& params,
                           std::vector<std::thread>& workers, Reader& reader,
                           OutputWriter& writer,
                           std::unique_ptr<SearchStrategy>& strategy,
                           length_t& maxOrIdentity) {
    if (params.inferPairedEndParameters) {
        inferParametersAndStartWorkers(reader, writer, strategy, maxOrIdentity,
                                       params);
    } else {
        strategy->setOrientation(params.orientation);
        for (size_t i = 0; i < workers.size(); i++) {
            workers[i] = std::thread(
                threadEntryPairedEnd, std::ref(reader), std::ref(writer),
                std::ref(strategy), std::ref(maxOrIdentity),
                std::ref(params.maxInsertSize), std::ref(params.minInsertSize));
        }
        logPairedEndParameters(params);
        for (auto& worker : workers) {
            worker.join();
        }
    }
}

/**
 * @brief Starts worker threads for processing and waits on them as well as the
 * reader and writer to finish.
 *
 * @param params Parameters object containing necessary parameters.
 * @param reader Reference to Reader object.
 * @param writer Reference to OutputWriter object.
 * @param strategy Pointer to SearchStrategy object.
 * @return true if worker threads start and process successfully.
 * @return false if there is an error during processing (error will be logged).
 */
bool startWorkerThreads(Parameters& params, Reader& reader,
                        OutputWriter& writer,
                        std::unique_ptr<SearchStrategy>& strategy) {
    try {
        auto start = std::chrono::high_resolution_clock::now();

        std::vector<std::thread> workers(params.nThreads);
        length_t maxOrIdentity = params.getMaxOrIdentity();

        if (params.sMode == SINGLE_END) {
            startSingleEndWorkers(workers, reader, writer, strategy,
                                  maxOrIdentity);
        } else {
            startPairedEndWorkers(params, workers, reader, writer, strategy,
                                  maxOrIdentity);
        }

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;

        logger.logInfo(std::to_string(params.nThreads) +
                       " mapping threads finished mapping in " +
                       std::to_string(elapsed.count()) + "s");

        reader.joinReaderThread();
        writer.joinWriterThread();
    } catch (const std::exception& e) {
        logger.logError("Problem during mapping: " + std::string(e.what()));
        return false;
    }
    return true;
}

/**
 * @brief Cleans up resources and exits the program.
 *
 * @param indexPtr Pointer to the index interface object.
 * @param readerPtr Pointer to Reader object.
 * @param writerPtr Pointer to OutputWriter object.
 */
void cleanupAndExit(std::unique_ptr<IndexInterface>& indexPtr,
                    std::unique_ptr<Reader>& readerPtr,
                    std::unique_ptr<OutputWriter>& writerPtr) {
    indexPtr.reset();
    readerPtr.reset();
    writerPtr.reset();
    logger.logInfo("Bye...\n");
}

/**
 * @brief Get a time as a string.
 *
 * @param time The time to convert to a string.
 * @return std::string The time as a string.
 */
std::string getTimeAsString(std::chrono::system_clock::time_point time) {
    std::time_t now_time = std::chrono::system_clock::to_time_t(time);
    std::stringstream ss;
    ss << std::ctime(&now_time);
    std::string timeStr = ss.str();
    timeStr.pop_back(); // remove the newline character
    return timeStr;
}

/**
 * @brief Get the current time as a string.
 *
 * @return std::string The current time.
 */
std::string getCurrentTime() {
    auto now = std::chrono::system_clock::now();
    return getTimeAsString(now);
}

/**
 * @brief Entry point of the program.
 *
 * @param argc Argument count.
 * @param argv Argument vector.
 * @return int Exit status.
 */
int main(int argc, char* argv[]) {
    logger.logDeveloper("Developer mode is on");
    Parameters params;
    if (!processArguments(argc, argv, params)) {
        return EXIT_FAILURE;
    }

    // set logger verbosity
    logger.setVerbose(params.verbose);

    std::unique_ptr<OutputWriter> writerPtr = createOutputWriter(params);
    if (!writerPtr) {
        return EXIT_FAILURE;
    }

    std::unique_ptr<Reader> readerPtr = createReader(params);
    if (!readerPtr) {
        writerPtr.reset();
        return EXIT_FAILURE;
    }

    std::unique_ptr<IndexInterface> indexPtr;

    // Conditional compilation
#ifdef RUN_LENGTH_COMPRESSION
    indexPtr = createBMove(params);
#else
    indexPtr = createFMIndex(params);
#endif

    if (!indexPtr) {
        return EXIT_FAILURE;
    }

    std::unique_ptr<SearchStrategy> strategy =
        createSearchStrategy(params, *indexPtr);
    if (!strategy) {
        return EXIT_FAILURE;
    }

    logAlignmentParameters(params, *indexPtr, strategy);

    auto startTime = std::chrono::high_resolution_clock::now();
    logger.logInfo("Start alignment at " + getTimeAsString(startTime));
    writerPtr->setStartTime(startTime);

    if (!startWorkerThreads(params, *readerPtr, *writerPtr, strategy)) {
        return EXIT_FAILURE;
    }

    cleanupAndExit(indexPtr, readerPtr, writerPtr);

    return EXIT_SUCCESS;
}

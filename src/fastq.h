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

#ifndef FASTQ_H
#define FASTQ_H

#include "definitions.h"  // for SequencingMode, length_t
#include "indexhelpers.h" // for TextOcc, PairedTextOccs, Counters
#include "logger.h"       // for Logger
#include "reads.h"        // for Read
#include "seqfile.h"      // for FileType

#include <algorithm> // for max, transform
#include <atomic>
#include <chrono> // for minutes, seconds, high_resolution_clock, time_p...
#include <condition_variable> // for condition_variable
#include <ctype.h>            // for toupper
#include <deque>              // for deque
#include <functional>         // for bind, function, _1
#include <memory>             // for allocator_traits<>::value_type
#include <mutex>              // for mutex, lock_guard
#include <stddef.h>           // for size_t
#include <string>             // for string, basic_string
#include <thread>             // for thread
#include <utility>            // for move
#include <vector>             // for vector

#include <parallel_hashmap/btree.h>

// ============================================================================
// SEQUENCE RECORD
// ============================================================================

/**
 * A sequence record object contains a single read that is read from file.
 */
class SequenceRecord : public ReadBundle {
  private:
    /**
     * Read record from read file
     * @param Opened read file
     * @return True upon successful read, false otherwise
     */
    bool readFromFileFASTA(SeqFile& input);
    /**
     * Read record from read file
     * @param Opened read file
     * @return True upon successful read, false otherwise
     */
    bool readFromFileFASTQ(SeqFile& input);

    void writeToFileFASTA(SeqFile& input) const;
    void writeToFileFASTQ(SeqFile& input) const;

    bool doTrim = false;
    length_t trimStart = 0;
    length_t trimEnd = 0;

    void trimRead() {
        // Store original values before modifying
        std::string originalRead = read;
        std::string originalQual = qual;

        try {
            assert(read.size() == qual.size());
            read = read.substr(trimStart, trimEnd - trimStart);
            qual = qual.substr(trimStart, trimEnd - trimStart);
     
        } catch (const std::out_of_range& e) {
            // Reset to original values if an error occurs
            read = originalRead;
            qual = originalQual;

            logger.logWarning("Trimming failed for read " + seqID +
                              ". The read is too short for the specified "
                              "trimming positions. The read will not be "
                              "trimmed.");
        }
    }

    void dummyTrim() {
    }

    std::function<void()> trimFunc;
    std::function<bool(SeqFile&)> readFunc;
    std::function<void(SeqFile&)> writeFunc;

  public:
    bool readFromFile(SeqFile& input) {
        if (readFunc(input)) {
            trimFunc(); // trim read if necessary
            makeReverseComplement();
            return true;
        }
        return false;
    }

    void writeToFile(SeqFile& input) const {
        writeFunc(input);
    }

    SequenceRecord(bool fastq, bool trim, length_t trimStart, length_t trimEnd)
        : ReadBundle(), doTrim(trim), trimStart(trimStart), trimEnd(trimEnd) {
        readFunc = fastq ? std::bind(&SequenceRecord::readFromFileFASTQ, this,
                                     std::placeholders::_1)
                         : std::bind(&SequenceRecord::readFromFileFASTA, this,
                                     std::placeholders::_1);
        writeFunc = fastq ? std::bind(&SequenceRecord::writeToFileFASTQ, this,
                                      std::placeholders::_1)
                          : std::bind(&SequenceRecord::writeToFileFASTA, this,
                                      std::placeholders::_1);

        trimFunc = doTrim ? std::bind(&SequenceRecord::trimRead, this)
                          : std::bind(&SequenceRecord::dummyTrim, this);
    }
};

// ============================================================================
// READ BLOCK CLASS
// ============================================================================

/**
 * A read block object contains a number of single-end or paired-end reads
 * that are read from file(s). In case of paired-end reads, the reads are
 * stored in an interleaved manner. In general, a readBlock consists of
 * multiple read chunks.
 */
class ReadBlock : public std::vector<SequenceRecord> {

  private:
    size_t nextChunkOffset; // offset of the next record to be read
    size_t targetChunkSize; // the used targetChunkSize for this block

  public:
    /**
     * Construct read block from a vector of FastQ records (deep copy)
     * @param buffer A vector of FastQ records
     */
    ReadBlock(const std::vector<SequenceRecord>& buffer)
        : std::vector<SequenceRecord>(buffer), nextChunkOffset(0) {
    }

    /**
     * Default, deleted copy and default move constructor
     */
    ReadBlock() : nextChunkOffset(0) {
    }
    ReadBlock(const ReadBlock&) = delete;
    ReadBlock(ReadBlock&&) = default;

    /**
     * Deleted copy and default move assignment operator
     */
    ReadBlock& operator=(const ReadBlock&) = delete;
    ReadBlock& operator=(ReadBlock&& rhs) = default;

    /**
     * Return true of at least one chunk is available for reading
     * @return true or false
     */
    bool chunkAvailable() const {
        return nextChunkOffset < this->size();
    }

    /**
     * Get next available input record chunk
     * @param buffer Record buffer to write to (contents will be appended)
     * @param pairedEnd If true, buffer.size() will always be even
     */
    void getNextChunk(std::vector<ReadBundle>& buffer, bool pairedEnd);

    /**
     * Read a block from file (paired-end reads)
     * @param file1 First fastq file
     * @param file2 First fastq file
     * @param targetBlockSize Desired number of nucleotides in this block
     * @param fastq1 True if the first file is a fastq file, false if it is a
     * fasta file
     * @param fastq2 True if the second file is a fastq file, false if it is a
     * fasta file
     * @param doTrim True if the reads should be trimmed, false otherwise
     * @param trimStart Start position for trimming
     * @param trimEnd End position for trimming
     */
    void readFromFile(SeqFile& file1, SeqFile& file2, size_t targetBlockSize,
                      bool fastq1, bool fastq2, bool doTrim, length_t trimStart,
                      length_t trimEnd);

    /**
     * Read a block from file (paired-end reads)
     * @param file1 Fastq file
     * @param targetBlockSize Desired number of nucleotides in this block
     * @param fastq True if the file is a fastq file, false if it is a fasta
     * file
     * @param doTrim True if the reads should be trimmed, false otherwise
     * @param trimStart Start position for trimming
     * @param trimEnd End position for trimming
     */
    void readFromFile(SeqFile& file1, size_t targetBlockSize, bool fastq,
                      bool doTrim, length_t trimStart, length_t trimEnd);

    /**
     * Write a block to file (paired-end reads)
     * @param file1 First fastq file
     * @param file2 First fastq file
     */
    void writeToFile(SeqFile& file1, SeqFile& file2);

    /**
     * Write a block to file (paired-end reads)
     * @param file Fastq file
     */
    void writeToFile(SeqFile& file);

    /**
     * Return the target chunk size that was used in this block
     * @return The target chunk size
     */
    size_t getTargetChunkSize() const {
        return targetChunkSize;
    }

    /**
     * Set the target chunk size that was used in this block
     * @param targetChunkSize The target chunk size
     */
    void setTargetChunkSize(size_t targetChunkSize) {
        this->targetChunkSize = targetChunkSize;
    }
};

// ============================================================================
// FASTQ READER
// ============================================================================

/**
 * @class Reader
 * @brief Class for reading FASTQ or FASTA files.
 *
 * The Reader class provides functionality for reading FASTQ/FASTA files,
 * both single-end and paired-end reads. It supports multithreaded reading,
 * allowing for efficient processing of large FASTQ/FASTA datasets.
 */
class Reader {
  private:
    std::string filename1;     // name of the input file /1
    FileType fileType1;        // file type (FASTQ, FASTQ.GZ, FASTA, FASTA.GZ)
    std::string baseFilename1; // file name w/o extension /1

    std::string filename2;     // name of the input file /2
    FileType fileType2;        // file type (FASTQ, FASTQ.GZ, FASTA, FASTA.GZ)
    std::string baseFilename2; // file name w/o extension /2

    bool fastq1 = true; // true if first file is fastq, false otherwise
    bool fastq2 = true; // true if second file is fastq, false otherwise

    bool doTrim = false; // true if reads should be trimmed, false otherwise
    length_t trimStart;  // start position for trimming
    length_t trimEnd;    // end position for trimming

    bool pairedEnd; // true for paired-end reads, false
                    // for single-end reads (ignore file 2)

    const int numReadBlocks = 2; // number of read blocks

    size_t targetChunkSize; // target size for a single chunk
    size_t targetBlockSize; // target size for a single block

    size_t numberOfWorkerThreads = 0; // number of worker threads

    std::thread iThread; // input thread

    // input thread variables (protect by inputMutex)
    std::mutex inputMutex;              // input thread mutex
    std::condition_variable inputReady; // input blocks ready condition
    std::vector<ReadBlock> inputBlocks; // input blocks that are empty

    // worker thread variables (protect by workMutex)
    std::mutex workMutex;              // worker thread mutex
    std::condition_variable workReady; // work ready condition
    std::deque<ReadBlock> workBlocks;  // blocks being processed
    size_t chunkID;                    // unique ID per chunk

    std::mutex processMutex; // mutex for processing time

    std::chrono::microseconds lowEndProcessingTime{5000};
    std::chrono::microseconds highEndProcessingTime{10000};

    std::chrono::microseconds midProcessingTime =
        (highEndProcessingTime + lowEndProcessingTime) / 2;
    std::vector<std::chrono::microseconds> processingTimes;

    std::exception_ptr exceptionPtr = nullptr; // Store exception

    /**
     * Entry routine for the input thread.
     * This function is responsible for reading the input files and
     * populating the input blocks.
     */
    void readerThread();

  public:
    /**
     * Default constructor.
     * @param filename1 Name of the input file /1.
     * @param filename2 Name of the input file /2.
     */
    Reader(const std::string& filename1, const std::string& filename2);

    /**
     * Set trimming positions for the reads.
     * All reads will be trimmed to the specified positions.
     * @param trimStart Start position for trimming.
     * @param trimEnd End position for trimming.
     */
    void setTrimming(length_t trimStart, length_t trimEnd) {
        doTrim = true;
        this->trimStart = trimStart;
        this->trimEnd = trimEnd;
    }

    /**
     * Return the filename without the extension /1.
     * @return The filename without the extension.
     */
    std::string getBaseFilename1() const {
        return baseFilename1;
    }

    /**
     * Return the filename without the extension /2.
     * @return The filename without the extension.
     */
    std::string getBaseFilename2() const {
        return baseFilename2;
    }

    /**
     * Start the reader thread.
     * @param targetChunkSize Target size for a single chunk.
     * @param numWorkerThreads Number of worker threads.
     */
    void startReaderThread(size_t targetChunkSize, size_t numWorkerThreads);

    /**
     * Get the next read chunk from the input.
     * @param chunk Buffer in which to store the records (output).
     * @param chunkID Unique chunk identifier (output).
     * @return False if no more reads are available, true otherwise.
     */
    bool getNextChunk(std::vector<ReadBundle>& chunk, size_t& chunkID);

    /**
     * Join the reader thread.
     * This function blocks until the reader thread has finished.
     */
    void joinReaderThread();

    /**
     * Add processing time of a single chunk
     * @param time Processing time of a single chunk
     */
    void addProcessingTime(std::chrono::microseconds time) {
        std::lock_guard<std::mutex> lock(processMutex);
        processingTimes.push_back(time);
    }
};

// ============================================================================
// OUTPUT RECORD
// ============================================================================

/**
 * An output record object contains the non-redundant occurrences for a single
 * read or read-pair. In case of paired-end reads, both paired and unpaired
 * occurrences are present.
 */
class OutputRecord {
  public:
    std::shared_ptr<std::vector<TextOcc>>
        outputOcc; // The non-redundant occurrences for this record if SE
    std::shared_ptr<std::vector<TextOcc>>
        unpairedOcc; // unpaired occurrences if PE
    std::shared_ptr<std::vector<PairedTextOccs>>
        pairOcc; // The non-redundant occurrences for this record if PE

    /**
     * Constructor with single-end occurrences
     */
    OutputRecord(std::vector<TextOcc>&& occurrences)
        : outputOcc(
              std::make_shared<std::vector<TextOcc>>(std::move(occurrences))) {
    }

    /**
     * Constructor with paired-end occurrences
     */
    OutputRecord(std::vector<PairedTextOccs>&& pairedOccurrences,
                 std::vector<TextOcc>&& unpairedOcc)
        : unpairedOcc(
              std::make_shared<std::vector<TextOcc>>(std::move(unpairedOcc))),
          pairOcc(std::make_shared<std::vector<PairedTextOccs>>(
              std::move(pairedOccurrences))){};

    // Move constructor
    OutputRecord(OutputRecord&& other) noexcept
        : outputOcc(std::move(other.outputOcc)),
          unpairedOcc(std::move(other.unpairedOcc)),
          pairOcc(std::move(other.pairOcc)){};

    // Move assignment operator
    OutputRecord& operator=(OutputRecord&& other) noexcept {
        if (this != &other) {
            outputOcc = std::move(other.outputOcc);
            unpairedOcc = std::move(other.unpairedOcc);
            pairOcc = std::move(other.pairOcc);
        }
        return *this;
    }

    // Destructor (if needed)
    ~OutputRecord() = default;
};

// ============================================================================
// OUTPUT CHUNK
// ============================================================================

/**
 * An output chunk object contains a number of output records that need to be
 * written to disk. As well as the counters that were gathered during the
 * processing of the records. A special end output chunk is used to signal
 * termination.
 */
class OutputChunk {
  private:
    bool empty = true; // If the chunk is filled in
    size_t chunkID;    // The unique id of this chunk
    bool end = false;  // is this an end signal chunk
    std::vector<OutputRecord>
        records; // the records for the input of this chunk

    Counters counters; // the performance counters for this chunk

    /**
     * Set this as an end outputChunk and mark the chunk as non-empty
     */
    void setEnd(size_t id) {
        end = true;
        set(id);
    }

  public:
    /**
     * @returns, true if the chunk has not been filled in yet (= worker has
     * not finished yet), false if the chunk is filled in
     */
    bool isEmpty() const {
        return empty;
    }

    /**
     * Get the id
     * @returns the chunk id
     */
    size_t getChunkID() const {
        return chunkID;
    }
    /**
     * Default constructor. Creates an empty output chunk with no records
     */
    OutputChunk() : empty(true), chunkID(-1), end(false), records{} {
    }

    /**
     * Clear the output chunk to reset it to a default constructed output
     * chunk
     */
    void clear() {
        empty = true;
        chunkID = -1;
        end = false;
        records.clear();
        counters.resetCounters();
    }

    /**
     * Set the chunkID and mark this OutputChunk as filled in
     * @param id the required id for this chunk
     */
    void set(size_t id) {
        chunkID = id;
        empty = false;
    }

    /**
     * Static function to create an end output with the provided id. An end
     * output is *never* empty.
     * @param id the id for this chunk
     * @returns an end output chunk with the provided id
     */
    static OutputChunk makeEndOutput(size_t id) {
        OutputChunk out = OutputChunk();
        out.setEnd(id);
        return out;
    }

    void addSingleEndRecord(std::vector<TextOcc>&& occ) {
        records.emplace_back(std::move(occ));
    }

    void addPairedEndRecord(std::vector<PairedTextOccs>&& pairedOcc,
                            std::vector<TextOcc>&& unpairedOcc) {
        records.emplace_back(std::move(pairedOcc), std::move(unpairedOcc));
    }
    /**
     * @returns true if this is an end OutputChunk
     */
    bool isEnd() const {
        return end;
    }

    /**
     * @returns the records for this chunk
     */
    const std::vector<OutputRecord>& getRecords() const {
        return records;
    }

    /**
     * @returns the records for this chunk
     */
    std::vector<OutputRecord>& getRecordsNonConst() {
        return records;
    }

    /**
     * Reserve space for the records
     */
    void reserve(size_t s) {
        records.reserve(s);
    }

    // Move constructor
    OutputChunk(OutputChunk&& other) noexcept
        : empty(other.empty), chunkID(other.chunkID), end(other.end),
          records(std::move(other.records)),
          counters(std::move(other.counters)) {
        other.clear();
    }

    // Move assignment operator
    OutputChunk& operator=(OutputChunk&& other) noexcept {
        if (this != &other) {
            empty = other.empty;
            chunkID = other.chunkID;
            end = other.end;
            records = std::move(other.records);
            counters = std::move(other.counters);

            other.clear();
        }
        return *this;
    }

    Counters& getMutableCounters() {
        return counters;
    }

    // Delete copy constructor
    OutputChunk(const OutputChunk&) = delete;

    ~OutputChunk() {
        clear();
    }
};

// ============================================================================
// OUTPUT WRITER
// ============================================================================

/**
 * Class that handles the writing to disk of occurrences.
 */
class OutputWriter {
  private:
    std::string writeFile; // The file to write to
    FileType fileType;     // The type of the writeFile
    bool reorder = false;  // whether we keep the original order
    size_t maxNumChunks; // the maximum no. of unprocessed chunks stored in the
                         // writer
    std::thread wThread; // thread to write to disk
    std::string headerFile; // The file to get the info about the index for
                            // the header from
    std::string commandLineParameters; // The command line parameters with
                                       // which the program was started
    const SequencingMode sMode; // The sequencing mode (paired- or single-end)

    // output thread variables (protect by outputMutex)
    std::mutex outputMapMutex;           // output thread mutex
    std::condition_variable outputReady; // output blocks ready condition
    std::condition_variable outputSpace; // output blocks ready condition
    phmap::btree_map<size_t, OutputChunk>
        outputChunks;   // output blocks to be written
    size_t nextChunkID; // next chunk to be written

    length_t terminationMessagesSent = 0;
    length_t terminationMessagesReceived = 0;
    length_t threads = 0;

    // variables for logging
    length_t logIntervalRecords =
        1024 * 8; // after how many processed records do we log
    std::chrono::time_point<std::chrono::high_resolution_clock> lastLogTime;
    const std::chrono::seconds logIntervalSeconds =
        std::chrono::seconds(60); // after how many seconds must we log
    const std::chrono::seconds minimumLogSeconds =
        std::chrono::seconds(5); // minimum time between two logs

    std::chrono::time_point<std::chrono::high_resolution_clock> startTime;

    // variables for stopping
    std::atomic<bool> stopRequested{false}; // Add stop flag

    /**
     * Writes chunks to file
     */
    void writerThread();

    /**
     * Log the processing of the records
     * @param processedRecords the number of records that have been processed
     * @param alwaysLog whether to always log, independent of how much time has
     * passed/records have been processed. (Default = false)
     */
    void logProcess(length_t processedRecords, bool alwaysLog = false);

  public:
    /**
     * Constructor
     * @param writeFile the file to write to
     * @param headerFile the file to get the info about the index for the
     * header from
     * @param commandLineParameters the command line parameters with which the
     * program was started
     * @param sMode the sequencing mode (paired- or single-end)
     * @param reorder whether to keep original order
     */
    OutputWriter(const std::string& writeFile, const std::string& headerFile,
                 const std::string& commandLineParameters,
                 const SequencingMode& sMode, bool reorder);

    ~OutputWriter() {
        stop();
        joinWriterThread();
    }

    /**
     * Start writer thread.
     *
     * @param maxNumChunks the maximum number of chunks to store concurrently in
     * the writer. Must be at least the number of threads.
     * @param threads the number of threads
     */
    void start(const size_t maxNumChunks, const size_t threads);

    /**
     * Stop the writer thread
     */
    void stop() {
        stopRequested = true;
        outputReady.notify_all(); // Wake up the thread if it's waiting
    }
    /**
     * Join the writer thread
     */
    void joinWriterThread();

    /**
     * Commit chunk of data to be written
     * @param chunk Buffer in which to store the records (output)
     */
    void commitChunk(OutputChunk& chunk);

    /**
     * Adds a termination message at id chunkID if no termination has been
     * set yet
     */
    void sendTermination(size_t chunkID);

    void setStartTime(
        const std::chrono::time_point<std::chrono::high_resolution_clock>&
            startTime) {
        this->startTime = startTime;
        lastLogTime = startTime;
    }
};
#endif

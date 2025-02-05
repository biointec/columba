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

#include "fastq.h"
#include "definitions.h" // for SINGLE_END, length_t, SUB_VERSION_NUMB...
#include "logger.h"      // for Logger, logger
#include "seqfile.h"     // for SeqFile, (anonymous), FileType, UNKNOW...
#include "util.h"        // for Util

#include <algorithm> // for transform, max, copy
#include <cassert>   // for assert
#include <fstream>   // for operator<<, basic_ostream, stringstream
#include <iomanip>   // for operator<<, setprecision
#include <numeric>   // for accumulate
#include <parallel_hashmap/btree.h> // for btree_map, btree_iterator, btree...
#include <ratio>                    // for ratio
#include <stdexcept>                // for runtime_error
#include <tuple>                    // for tie, tuple

using namespace std;

// ============================================================================
// FASTQ RECORD
// ============================================================================

bool SequenceRecord::readFromFileFASTQ(SeqFile& inputFile) {
    // read sequence identifier
    char c = inputFile.peekCharacter();

    // end of file might be reached
    if (!inputFile.good())
        return false;
    if (c != '@')
        throw ios::failure("File " + inputFile.getFileName() +
                           " doesn't appear to be in FastQ format");
    inputFile.getLine(seqID);
    if (!seqID.empty() && seqID.back() == '\n')
        seqID.pop_back();

    // read the actual read
    inputFile.getLine(read);
    if (!read.empty() && read.back() == '\n')
        read.pop_back();

    // read the + line
    string dummy;
    inputFile.getLine(dummy);
    // read the quality scores
    inputFile.getLine(qual);
    if (!qual.empty() && qual.back() == '\n')
        qual.pop_back();

    cleanUpRecord();
    ReadBundle::makeReverseComplement();

    return !read.empty();
}

void SequenceRecord::writeToFileFASTQ(SeqFile& inputFile) const {
    inputFile.writeLine(seqID);
    inputFile.writeChar('\n');
    inputFile.writeLine(read);
    inputFile.writeLine("\n+\n");
    inputFile.writeLine(qual);
    inputFile.writeChar('\n');
}

bool SequenceRecord::readFromFileFASTA(SeqFile& inputFile) {
    clear();

    // end of file might be reached
    if (!inputFile.good())
        return false;
    // read sequence identifier
    char c = inputFile.peekCharacter();

    if (c != '>')
        throw ios::failure("File doesn't appear to be in Fasta format");
    inputFile.getLine(seqID);
    if (!seqID.empty() && seqID.back() == '\n')
        seqID.pop_back();

    // read the actual read line by line until a new sequence identifier is
    // found or the end of the file is reached
    string line;
    while (inputFile.good()) {
        c = inputFile.peekCharacter();
        if (c == '>') {
            break;
        }
        inputFile.getLine(line);
        if (!line.empty() && line.back() == '\n')
            line.pop_back();
        read += line;
    }

    this->qual = "*";
    cleanUpRecord();
    ReadBundle::makeReverseComplement();

    return !read.empty();
}

void SequenceRecord::writeToFileFASTA(SeqFile& inputFile) const {
    inputFile.writeLine(seqID);
    inputFile.writeChar('\n');
    // write read in 80 character lines
    const size_t lineLength = 80;
    for (size_t i = 0; i < read.size(); i += lineLength) {
        inputFile.writeLine(read.substr(i, lineLength));
    }

    inputFile.writeChar('\n');
}

// ============================================================================
// READ BLOCK CLASS
// ============================================================================

void ReadBlock::getNextChunk(vector<ReadBundle>& buffer, bool pairedEnd) {

    size_t remainingElements = this->size() - nextChunkOffset;
    size_t bufferSize =
        std::min(remainingElements,
                 targetChunkSize + (pairedEnd && (targetChunkSize % 2 == 1)));

    buffer.resize(bufferSize);
    auto sourceStart = this->begin() + nextChunkOffset;
    auto sourceEnd = sourceStart + bufferSize;
    std::move(sourceStart, sourceEnd, buffer.begin());

    nextChunkOffset += bufferSize;
}

void ReadBlock::readFromFile(SeqFile& file1, SeqFile& file2,
                             size_t targetBlockSize, bool fastq1, bool fastq2) {
    // clear the block
    this->clear();
    nextChunkOffset = 0;

    // read in new contents
    SequenceRecord recordA(fastq1), recordB(fastq2);
    size_t thisBlockSize = 0;

    while (thisBlockSize < targetBlockSize) {
        bool flag1 = recordA.readFromFile(file1);
        bool flag2 = recordB.readFromFile(file2);
        if (!flag1 || !flag2)
            break;

        this->push_back(recordA);
        this->push_back(recordB);
        thisBlockSize += 1;
        thisBlockSize += 1;
    }

    if (file1.good() != file2.good())
        logger.logWarning(
            "paired-end FastQ files contain different  number of reads\n");
}

void ReadBlock::readFromFile(SeqFile& file, size_t targetBlockSize,
                             bool fastq) {
    // clear the block
    this->clear();
    nextChunkOffset = 0;

    // read in new contents
    SequenceRecord record(fastq);
    size_t thisBlockSize = 0;

    while (thisBlockSize < targetBlockSize) {
        if (!record.readFromFile(file))
            break;
        this->push_back(record);
        thisBlockSize += 1;
    }
}

void ReadBlock::writeToFile(SeqFile& file1, SeqFile& file2) {
    // make sure we have an even number of read records
    assert(this->size() % 2 == 0);

    for (size_t i = 0; i < this->size(); i += 2) {
        (*this)[i].writeToFile(file1);
        (*this)[i + 1].writeToFile(file2);
    }
}

void ReadBlock::writeToFile(SeqFile& file) {
    for (const auto& r : *this)
        r.writeToFile(file);
}

// ============================================================================
// FASTQ READER
// ============================================================================

Reader::Reader(const string& filename1, const string& filename2)
    : filename1(filename1), fileType1(UNKNOWN_FT), filename2(filename2),
      fileType2(UNKNOWN_FT), pairedEnd(!filename2.empty()) {
    // try to figure out the file format based on the ext
    string ext;

    tie(fileType1, baseFilename1) = getFileType(filename1);
    if (fileType1 == UNKNOWN_FT) {
        string msg = "don't know how to open file: \"" + filename1 +
                     "\"\nExpected one of the following extensions: "
                     ".fastq, .fq, .fasta or .fa  (or .gz variants thereof)\n";
        throw runtime_error(msg);
    }

    if (!Util::fileExists(filename1))
        throw runtime_error("cannot open file " + filename1);

    fastq1 = (fileType1 == FASTQ || fileType1 == FASTQ_GZ);

    // no second file specified, assume single-ended reads
    if (!pairedEnd)
        return;

    tie(fileType2, baseFilename2) = getFileType(filename2);
    if (fileType2 == UNKNOWN_FT) {
        string msg = "don't know how to open file: \"" + filename2 +
                     "\"\nExpected one of the following extensions: "
                     ".fastq or .fq (or .gz variants thereof)\n";
        throw runtime_error(msg);
    }

    if (!Util::fileExists(filename2))
        throw runtime_error("cannot open file " + filename2);

    fastq2 = (fileType2 == FASTQ || fileType2 == FASTQ_GZ);
}

void Reader::readerThread() {
    // open the read file(s)
    SeqFile file1(fileType1 == FASTQ_GZ || fileType1 == FASTA_GZ);
    file1.open(filename1);

    SeqFile file2(fileType2 == FASTQ_GZ || fileType2 == FASTA_GZ);
    if (pairedEnd)
        file2.open(filename2);

    // recheck the chunksize every 8 blocks to start
    size_t checkPoint = 8;
    size_t blockCounter = 0;

    while (file1.good()) {
        // A) wait until an input block is available
        unique_lock<mutex> inputLock(inputMutex);
        inputReady.wait(inputLock, [this] { return !inputBlocks.empty(); });

        // get the input block
        ReadBlock block = std::move(inputBlocks.back());
        inputBlocks.pop_back();
        inputLock.unlock();

        block.setTargetChunkSize(targetChunkSize);

        // B) fill up the block (only this thread has access)
        // no mutexes are held by this thread at this point
        if (pairedEnd)
            block.readFromFile(file1, file2, targetBlockSize, fastq1, fastq2);
        else
            block.readFromFile(file1, targetBlockSize, fastq1);

        // empty block: end-of-file reached
        if (block.empty())
            break;

        // C) push the record block onto the worker queue
        unique_lock<mutex> workLock(workMutex);
        workBlocks.push_back(std::move(block));
        workLock.unlock();
        workReady.notify_all();

        // D) check if recalculation of target chunk size is necessary
        blockCounter++;
        if (blockCounter >= checkPoint) {
            // Acquire the processing time mutex
            unique_lock<mutex> processLock(processMutex);

            if (!processingTimes.empty()) {
                std::chrono::microseconds averageProcessingTime =
                    accumulate(processingTimes.begin(), processingTimes.end(),
                               std::chrono::microseconds(0)) /
                    processingTimes.size();

                // Determine if we need to adjust based on processing times
                bool change = false;
                auto prevSize = targetChunkSize;

                if (averageProcessingTime < lowEndProcessingTime ||
                    averageProcessingTime > highEndProcessingTime) {
                    int multFactor =
                        (midProcessingTime / averageProcessingTime);
                    targetChunkSize *= multFactor;
                    targetChunkSize = (targetChunkSize / 64) * 64;
                    targetChunkSize =
                        max(targetChunkSize, 1ul); // Ensure it's at least 1
                    targetBlockSize = targetChunkSize * numberOfWorkerThreads;
                    change = targetChunkSize != prevSize;
                }

                // Adjust checkpoint size based on processing times
                if (change) {
                    checkPoint = 8; // Reset checkpoint if there was a change
                    logger.logDeveloper(
                        "New target chunk size based on processing time: " +
                        to_string(targetChunkSize));
                } else if (checkPoint < 1024) {
                    checkPoint *= 2; // Increase checkpoint size
                    logger.logDeveloper("Check point increased to: " +
                                        to_string(checkPoint));
                }
            }

            processingTimes.clear(); // Clear processing times
            blockCounter = 0;        // Reset block counter after processing
        }
    }

    file1.close();
    if (pairedEnd)
        file2.close();

    // send a termination message to the workers
    unique_lock<mutex> workLock(workMutex);
    workBlocks.push_back(ReadBlock()); // empty block == termination
    workReady.notify_all();
    workLock.unlock();

    logger.logDeveloper("Reader thread finished");
}

bool Reader::getNextChunk(vector<ReadBundle>& chunk, size_t& chunkID) {
    // wait until work becomes available
    unique_lock<mutex> workLock(workMutex);
    workReady.wait(workLock, [this] { return !workBlocks.empty(); });

    // check if it is a termination message (empty block)
    if (workBlocks.front().empty()) {
        chunkID = this->chunkID;
        return false;
    }

    // get a chunk of work (contents of chunk will be overwritten)
    workBlocks.front().getNextChunk(chunk, pairedEnd);
    chunkID = this->chunkID++;

    // last chunk: move block from work queue back to the input queue
    if (!workBlocks.front().chunkAvailable()) {
        ReadBlock tmp = std::move(workBlocks.front());
        workBlocks.pop_front();
        workLock.unlock();

        unique_lock<mutex> inputLock(inputMutex);
        inputBlocks.push_back(std::move(tmp));
        inputReady.notify_one();
        inputLock.unlock();
    }

    assert(!chunk.empty());
    return true;
}

void Reader::startReaderThread(size_t targetChunkSize, size_t numWorkThreads) {
    this->targetChunkSize = targetChunkSize;
    this->numberOfWorkerThreads = numWorkThreads;
    this->targetBlockSize = targetChunkSize * numWorkThreads;

    chunkID = 0;

    // initialize (empty) blocks on the input stack
    inputBlocks.resize(numReadBlocks * 2);

    // start reader thread
    iThread = thread(&Reader::readerThread, this);
}

void Reader::joinReaderThread() {
    iThread.join();

    inputBlocks.clear();
    workBlocks.clear();
}

// ============================================================================
// OUTPUT WRITER
// ============================================================================

OutputWriter::OutputWriter(const string& writeFile, const string& headerFile,
                           const string& commandLineParameters,
                           const SequencingMode& sequencingMode, bool reorder)
    : writeFile(writeFile), fileType(UNKNOWN_FT), reorder(reorder),
      headerFile(headerFile), commandLineParameters(commandLineParameters),
      sMode(sequencingMode), lastLogTime(chrono::steady_clock::now()),
      startTime(chrono::steady_clock::now()) {

    // try to figure out the file format based on the ext
    string ext;
    fileType = UNKNOWN_FT;

    if (writeFile.length() >= 4)
        ext = writeFile.substr(writeFile.length() - 4);
    transform(ext.begin(), ext.end(), ext.begin(), ::toupper);
    if (ext == ".SAM") {
        fileType = SAM;
    } else if (ext == ".RHS") {
        fileType = RHS;
    }

    if (writeFile.length() >= 7)
        ext = writeFile.substr(writeFile.length() - 7);
    transform(ext.begin(), ext.end(), ext.begin(), ::toupper);
    if (ext == ".SAM.GZ") {
        fileType = SAM_GZ;
    } else if (ext == ".RHS.GZ") {
        fileType = RHS_GZ;
    }

    if (fileType == UNKNOWN_FT) {
        string msg = "don't know how to open file: \"" + writeFile +
                     "\"\nExpected one of the following extensions: "
                     ".sam or .sam.gz\n";
        throw runtime_error(msg);
    }
}

void OutputWriter::start(const size_t size, const size_t threads) {
    maxNumChunks = size;
    terminationMessagesReceived = 0, terminationMessagesSent = 0;
    this->threads = threads;
    nextChunkID = 0; // waiting on zeroth chunk
    // start writer thread
    wThread = thread(&OutputWriter::writerThread, this);
}

void OutputWriter::joinWriterThread() {
    wThread.join();
    outputChunks.clear();
}

void OutputWriter::commitChunk(OutputChunk& chunk) {
    // wait until space becomes available to store the read block
    // there is *always* space to store the next-to-be-written block
    size_t chunkID = chunk.getChunkID();
    // acquire the mutex
    unique_lock<mutex> outputLock(outputMapMutex);
    outputSpace.wait(outputLock, [this, chunkID] {
        return (reorder) ? ((chunkID == nextChunkID) ||
                            (outputChunks.size() < maxNumChunks - 1))
                         : outputChunks.size() < maxNumChunks;
    });
    // move chunk into map
    outputChunks[chunkID] = std::move(chunk);

    bool wasEmpty = outputChunks.size() == 1;
    bool wasNext = chunkID == nextChunkID;

    outputLock.unlock();

    // signal the output thread if necessary
    if ((!reorder && wasEmpty) || wasNext) {
        outputReady.notify_one();
    }
}

void OutputWriter::logProcess(length_t processedRecords, bool alwaysLog) {
    auto now = chrono::steady_clock::now();
    auto timeElapsedSinceLastLog =
        chrono::duration_cast<chrono::seconds>(now - lastLogTime);

    string readsOrPairs = (sMode == SINGLE_END) ? "reads" : "read pairs";
    string recordsOrPairs =
        (sMode == SINGLE_END) ? "read records" : "read pairs";

    if (alwaysLog || (timeElapsedSinceLastLog >= minimumLogSeconds &&
                      (processedRecords % logIntervalRecords == 0 ||
                       timeElapsedSinceLastLog >= logIntervalSeconds))) {

        auto timeElapsed =
            chrono::duration_cast<chrono::milliseconds>(now - startTime)
                .count();

        // timeElapsed cannot be zero as we exceed the minimum logTime
        double rate = (processedRecords / static_cast<double>(timeElapsed));

        stringstream ss;
        ss << "Processed " << processedRecords << " " << recordsOrPairs
           << ". Rate: " << fixed << setprecision(0) << rate * 1000 << " "
           << readsOrPairs << "/second";

        logger.logInfo(ss);
        lastLogTime = now;
    }
}

void OutputWriter::writerThread() {
    // open the read file(s)
    SeqFile writeSeq(fileType == SAM_GZ || fileType == RHS_GZ);
    writeSeq.open(writeFile, WRITE);

    bool isSAM = (fileType == SAM || fileType == SAM_GZ);

    Counters counters;
    counters.resetCounters();

    // write the  headers
    if (isSAM) {
        writeSeq.writeLine("@HD\tVN:1.6\tSO:queryname\n");
        writeSeq.writeLine("@PG\tID:Columba" +
                           to_string(VERSION_NUMBER_COLUMBA) + "." +
                           to_string(SUB_VERSION_NUMBER_COLUMBA) +
                           "\tPN:Columba\tCL:" + commandLineParameters + "\n");

        ifstream ifs(headerFile, ios::binary);

        if (!ifs) {
            logger.logError("Unable to open the header file\n");
        } else {
            string line;
            while (getline(ifs, line)) {
                writeSeq.writeLine(line + "\n");
            }
        }
    } else {
        if (sMode == SINGLE_END) {
            writeSeq.writeLine("ReadID\tHits\n");
        } else {
            throw std::runtime_error(
                "Output format not supported for paired-end alignment\n");
        }
    }

    length_t processedRecords = 0;

    while (true) {
        // A) wait until the next chunk is available
        OutputChunk chunk;
        {
            unique_lock<mutex> outputLock(outputMapMutex);
            outputReady.wait(outputLock, [this] {
                return !outputChunks.empty() &&
                       (!reorder ||
                        (outputChunks.begin()->first == nextChunkID));
            });

            // get the output chunk
            chunk = std::move(outputChunks.begin()->second);
            outputChunks.erase(outputChunks.begin());
            nextChunkID++; // not necessary if !reorder but does not matter
            outputSpace.notify_all();
            outputLock.unlock();
        }

        // check for termination message
        if (chunk.isEnd()) {
            terminationMessagesReceived++;
            if (terminationMessagesReceived == threads) {
                break;
            }
            continue;
        }

        // add the counters
        counters.addCounters(chunk.getMutableCounters());

        // B) write the chunk (only this thread has access)
        // no thread is held at this point
        for (const OutputRecord& record : chunk.getRecords()) {

            if (sMode == SINGLE_END) {
                counters.inc(Counters::NUMBER_OF_READS);

                bool mapped = !record.outputOcc->empty() &&
                              record.outputOcc->front().isValid();
                counters.inc(Counters::MAPPED_READS, mapped);

                counters.inc(Counters::TOTAL_UNIQUE_MATCHES,
                             (mapped) ? record.outputOcc->size() : 0);

                // write the occurrences  to file
                for (const auto& match : *record.outputOcc) {
                    writeSeq.writeLine(match.getOutputLine());
                }

            } else {
                counters.inc(Counters::NUMBER_OF_READS, 2); // two read per pair

                const auto& pairedOccs = record.pairOcc;
                const auto& unpairedOccs = record.unpairedOcc;

                bool mapped = !pairedOccs->empty() &&
                              pairedOccs->front().getUpStream().isValid() &&
                              pairedOccs->front().getDownStream().isValid();
                bool discordantlyMapped =
                    mapped && pairedOccs->front().isDiscordant();
                bool mappedHalf =
                    !(mapped) &&
                    (!pairedOccs->empty() &&
                     (pairedOccs->front().getUpStream().isValid() ||
                      pairedOccs->front().getDownStream().isValid()));
                bool unpairedButMapped = !unpairedOccs->empty();

                counters.inc(Counters::TOTAL_UNIQUE_PAIRS,
                             (mapped) ? pairedOccs->size() : 0);

                counters.inc(Counters::MAPPED_PAIRS, mapped);
                counters.inc(Counters::DISCORDANTLY_MAPPED_PAIRS,
                             discordantlyMapped);
                counters.inc(Counters::MAPPED_HALF_PAIRS, mappedHalf);
                counters.inc(Counters::UNPAIRED_BUT_MAPPED_PAIRS,
                             unpairedButMapped);

                // write the pairs to file
                bool firstWrite = true;
                for (const auto& pair : *pairedOccs) {
                    writeSeq.writeLine(pair.getUpStream().getOutputLine());
                    if (!mappedHalf || firstWrite) {
                        writeSeq.writeLine(
                            pair.getDownStream().getOutputLine());
                    }
                    firstWrite = false;
                }
                // then write unpaired occurrences to file
                for (const auto& match : *unpairedOccs) {
                    writeSeq.writeLine(match.getOutputLine());
                }
            }

            logProcess(++processedRecords);
        }
    }

    writeSeq.close();

    logProcess(processedRecords, true);
    logger.logInfo("Finished writing to output file: " + writeFile);

    // add dropped unique matches
    counters.inc(Counters::TOTAL_UNIQUE_MATCHES,
                 counters.get(Counters::DROPPED_UNIQUE_MATCHES));
    counters.reportStatistics(sMode);
}

void OutputWriter::sendTermination(size_t chunkID) {
    // send a termination message (empty block)
    unique_lock<mutex> outputLock(outputMapMutex);

    outputChunks[chunkID + terminationMessagesSent] =
        OutputChunk::makeEndOutput(chunkID + terminationMessagesSent);
    terminationMessagesSent++;

    // signal the output thread if necessary
    if (!reorder || chunkID == nextChunkID) {

        outputReady.notify_one();
    }

} // unique_lock goes out of scope

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

#ifndef READS_H
#define READS_H

#include "definitions.h"
#include "nucleotide.h"

#include <string>

/**
 * A class to represent a read with a sequence identifier, the read itself and a
 * quality string
 */
class Read {
  protected:
    std::string seqID; /// sequence identifier
    std::string read;  /// read itself
    std::string qual;  /// quality string

    /**
     * Clean up a read record by removing the @ from the sequence identifier
     * and replacing all non-ACGT characters with N.
     */
    void cleanUpRecord() {
        // Remove everything after the first space in seqID, if space exists
        size_t spacePos = seqID.find(' ');
        if (spacePos != std::string::npos) {
            seqID.erase(spacePos);
        }

        // Remove the first character (@ or >)
        seqID = seqID.substr(1);

        // convert read to upper case
        std::transform(read.begin(), read.end(), read.begin(), [](unsigned char c) { return std::toupper(c); });

        // Replace all non-ACGT characters with N
        replaceNonACTGWithN();
    }

  public:
    /**
     * Clear the record
     */
    void clear() {
        seqID.clear();
        read.clear();
        qual.clear();
    }

    /**
     * Get a const-reference to the sequence identifier
     * @return a const-reference to the sequence identifier
     */
    const std::string& getSeqID() const {
        return seqID;
    }

    /**
     * Get a const-reference to the read
     * @return a const-reference to the read
     */
    const std::string& getRead() const {
        return read;
    }

    /**
     * Get a const-reference to the quality string
     * @return a const-reference to the quality string
     */
    const std::string& getQual() const {
        return qual;
    }

    /**
     * Replaces all non ACTG characters of the read with N
     */
    void replaceNonACTGWithN() {
        std::replace_if(
            read.begin(), read.end(),
            [](char c) { return !Nucleotide::isACGT(c); }, 'N');
    }

    /**
     * Constructor
     */
    Read(std::string& seqID, const std::string& read, const std::string& qual)
        : seqID(seqID), read(Nucleotide::uppercase(read)), qual(qual) {
        cleanUpRecord();
    }

    /**
     * Copy constructor
     */
    Read(const Read& other)
        : seqID(other.seqID), read(other.read), qual(other.qual) {
    }

    /**
     * Default constructor
     */
    Read() : qual("*") {
    }
};

/**
 * A class to bundle a read with its reverse complement and reversed quality
 */
class ReadBundle : public Read {
  private:
    std::string reverseComplement; // the reverse complement of the read
    std::string revQuality;        // the reversed quality of the read

    ReadBundle(const std::string& seq, bool isRevComp) {
        if (isRevComp) {
            reverseComplement = seq;
        } else {
            read = seq;
        }
    }

  protected:
    void makeReverseComplement() {
        reverseComplement = Nucleotide::getRevComplWithN(read);
    }

  public:
    ReadBundle() : Read(), reverseComplement(""), revQuality("") {
    }
    ReadBundle(const Read& read)
        : Read(read),
          reverseComplement(Nucleotide::getRevComplWithN(read.getRead())),
          revQuality("") {
        // leave reverse quality untouched for now, only consider it when needed
    }

    const std::string& getRevComp() const {
        return reverseComplement;
    }

    const std::string& getRevQuality() {
        if (revQuality.empty()) {
            revQuality = getQual();
            std::reverse(revQuality.begin(), revQuality.end());
        }
        return revQuality;
    }

    const size_t size() const {
        return getRead().size();
    }

    /**
     * Get the sequence of the read along the correct strand
     * @param strand the strand to get the sequence for
     * @return the sequence of the read along the correct strand
     */
    const std::string& getSequence(Strand strand) const {
        return (strand == Strand::FORWARD_STRAND) ? getRead() : getRevComp();
    }

    /**
     * Create a ReadBundle from a sequence and a strand. WARNING: This will only
     * have the sequence along this strand!
     * @param seq the sequence of the read
     * @param strand the strand of the read
     */
    static ReadBundle createBundle(const std::string& seq,
                                   const Strand strand) {
        return ReadBundle(seq, strand == Strand::REVERSE_C_STRAND);
    }
};

/**
 * A pair of reads for paired-end alignment.
 */
class ReadPair {
  private:
    ReadBundle read1;
    ReadBundle read2;

  public:
    ReadBundle& getBundle1() {
        return read1;
    }

    ReadBundle& getBundle2() {
        return read2;
    }

    const ReadBundle& getBundle1() const {
        return read1;
    }

    const ReadBundle& getBundle2() const {
        return read2;
    }

    const std::string& getRead1() const {
        return read1.getRead();
    }

    const std::string& getRead2() const {
        return read2.getRead();
    }

    const std::string& getRevComp1() const {
        return read1.getRevComp();
    }

    const std::string& getRevComp2() const {
        return read2.getRevComp();
    }

    ReadPair(const ReadBundle& read1, const ReadBundle& read2)
        : read1(read1), read2(read2) {
    }

    /**
     * Get the sequence of the read along the correct strand and with the
     * correct pair status
     * @param strand the strand to get the sequence for
     * @param status the pair status (first or second in pair) to get the
     * sequence for
     * @return the sequence of the read along the correct strand and with the
     * correct pair status
     */
    const std::string& getSequence(Strand strand, PairStatus status) const {
        return (status == FIRST_IN_PAIR) ? read1.getSequence(strand)
                                         : read2.getSequence(strand);
    }
};

#endif
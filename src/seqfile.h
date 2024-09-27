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
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *   GNU General Public License for more details.                             *
 *                                                                            *
 *   You should have received a copy of the GNU Affero General Public License *
 *   along with this program; if not, see <https://www.gnu.org/licenses/>.    *
 ******************************************************************************/

#ifndef READFILE_H
#define READFILE_H

#include <cstdio>  // for fgetc, feof, fgets, fputc, fputs, ungetc, NULL, FILE
#include <cstdlib> // for free, malloc
#include <fstream> // for ostream
#include <string>  // for string

#ifdef HAVE_ZLIB
#include "zlib.h" // for gzeof, gzgets, gzputc, gzungetc, gzwrite, gzgetc
#endif

// ============================================================================
// ENUM TYPES
// ============================================================================

typedef enum { READ, WRITE } ReadFileMode;
typedef enum {
    FASTQ,
    FASTQ_GZ,
    FASTA,
    FASTA_GZ,
    SAM,
    SAM_GZ,
    RHS, // read hit summary
    RHS_GZ,
    UNKNOWN_FT
} FileType;

std::ostream& operator<<(std::ostream& out, const FileType& fileType);

// ============================================================================
// READFILE HANDLER
// ============================================================================

/**
 * Abstract class for reading files
 */
class ReadFileHandler {

  protected:
    ReadFileMode mode;

  public:
    /**
     * Virtual destructor for your convenience
     */
    virtual ~ReadFileHandler() {};

    /**
     * Open a file with a given filename
     * @param filename File to open
     * @param mode Read (default) or write
     */
    virtual void open(const std::string& filename,
                      ReadFileMode mode = READ) = 0;

    /**
     * Check if the input is still valid
     * @return True of false
     */
    virtual bool good() = 0;

    /**
     * Read a line from file
     * @return Pointer to an internal buffer containing the line
     */
    virtual void getLine(std::string& line) = 0;

    /**
     * Write a line to file
     * @param line String to write
     */
    virtual void writeLine(const std::string& line) = 0;

    /**
     * Write a character to file
     * @param c Character to write
     */
    virtual void writeChar(char c) = 0;

    /**
     * Peek at the next character in the stream
     * @return The character
     */
    virtual char peekCharacter() = 0;

    /**
     * Get one character from the stream
     * @return The character
     */
    virtual char getCharacter() = 0;

    /**
     * Close the file
     */
    virtual void close() = 0;

    /**
     * Reset the input file
     */
    virtual void reset() = 0;
};

// ============================================================================
// REGULAR READFILE HANDLER
// ============================================================================

/**
 * Read file handler for regular files
 */
class RegularReadFileHandler : public ReadFileHandler {

  protected:
    FILE* fh;                 // file handler
    const int bufSize = 1024; // buffer size
    char* buffer;             // buffer

  public:
    RegularReadFileHandler() {
        buffer = (char*)malloc(bufSize * sizeof(char));
    }

    /**
     * Virtual destructor for you convenience
     */
    ~RegularReadFileHandler() {
        free(buffer);
    }

    /**
     * Open a file with a given filename
     * @param filename File to open
     * @param mode Read (default) or write
     */
    void open(const std::string& filename, ReadFileMode mode = READ);

    /**
     * Check if the input is still valid
     * @return True of false
     */
    bool good() {
        return !feof(fh);
    }

    /**
     * Read a line from file
     */
    void getLine(std::string& line) {
        line.clear();

        while (fgets(buffer, bufSize, fh) != NULL) {
            line.append(buffer);
            if (line.back() == '\n')
                break;
        }
    }

    /**
     * Write a line to file
     * @param line String to write
     */
    void writeLine(const std::string& line) {
        fputs(line.c_str(), fh);
    }

    /**
     * Write a character to file
     * @param c Character to write
     */
    void writeChar(char c) {
        fputc(c, fh);
    }

    /**
     * Peek at the next character in the stream
     * @return The character
     */
    char peekCharacter() {
        char c = fgetc(fh);
        ungetc(c, fh);
        return c;
    }

    /**
     * Get one character from the stream
     * @return The character
     */
    char getCharacter() {
        char c = fgetc(fh);
        return c;
    }

    /**
     * Close the file
     */
    void close();

    /**
     * Reset the input file
     */
    void reset();
};

// ============================================================================
// GZIPPED READFILE HANDLER
// ============================================================================

#ifdef HAVE_ZLIB

/**
 * Read file handler for gzipped files
 */
class GZipReadFileHandler : public ReadFileHandler {

  protected:
    gzFile ifs;               // gzipped input file stream
    const int bufSize = 1024; // buffer size
    char* buffer;             // buffer

  public:
    /**
     * Default constructor
     */
    GZipReadFileHandler() {
        buffer = (char*)malloc(bufSize * sizeof(char));
    }

    /**
     * Destructor
     */
    ~GZipReadFileHandler() {
        free(buffer);
    }

    /**
     * Open a file with a given filename
     * @param filename File to open
     * @param mode Read (default) of write
     * @return True upon success, false otherwise
     */
    void open(const std::string& filename, ReadFileMode mode);

    /**
     * Check if the input is still valid
     * @return True of false
     */
    bool good() {
        return !gzeof(ifs);
    }

    void getLine(std::string& line) {
        line.clear();

        while (gzgets(ifs, buffer, bufSize) != NULL) {
            line.append(buffer);
            if (line.back() == '\n')
                break;
        }
    }

    /**
     * Write a line to file
     * @param line String to write
     */
    void writeLine(const std::string& line) {
        gzwrite(ifs, line.c_str(), line.length());
    }

    /**
     * Write a character to file
     * @param c Character to write
     */
    void writeChar(char c) {
        gzputc(ifs, c);
    }

    /**
     * Peek at the next character in the stream
     * @return The character
     */
    char peekCharacter() {
        char c;
        c = gzgetc(ifs);
        gzungetc(c, ifs);
        return c;
    }

    /**
     * Get one character from the stream
     * @return The character
     */
    char getCharacter() {
        return gzgetc(ifs);
    }

    /**
     * Close the file
     */
    void close();

    /**
     * Reset the input file
     */
    void reset();
};

#endif

// ============================================================================
// SEQUENCE FILE
// ============================================================================

/**
 * Class that represents a sequence file
 */
class SeqFile {

  protected:
    ReadFileHandler* rfHandler; // read file handler
    std::string filename;       // filename

  public:
    /**
     * Default constructor
     * @param gzipped True if the file is gzipped
     */
    SeqFile(bool gzipped);

    /**
     * Destructor
     */
    virtual ~SeqFile() {
        delete rfHandler;
    }

    /**
     * Open a file
     * @param filename File to open
     * @param mode READ or WRITE
     */
    void open(const std::string& filename, ReadFileMode mode = READ) {
        rfHandler->open(filename, mode);
        this->filename = filename;
    }

    /**
     * Read a line from file
     * @param line Line that was read (empty if failed)
     */
    void getLine(std::string& line) {
        rfHandler->getLine(line);
    }

    /**
     * Write a line to file
     * @param line String to write
     */
    void writeLine(const std::string& line) {
        rfHandler->writeLine(line);
    }

    /**
     * Write a character to file
     * @param c Character to write
     */
    void writeChar(char c) {
        rfHandler->writeChar(c);
    }

    /**
     * Peek at the next character in the stream
     * @return The character
     */
    char peekCharacter() {
        return rfHandler->peekCharacter();
    }

    /**
     * Get one character from the stream
     * @return The character
     */
    char getCharacter() {
        return rfHandler->getCharacter();
    }

    /**
     * Check if the input is still valid
     * @return True of false
     */
    bool good() {
        return rfHandler->good();
    }

    /**
     * Reset file to starting position
     */
    void reset() {
        rfHandler->reset();
    }

    /**
     * Close a read file
     */
    void close() {
        rfHandler->close();
    }

    std::string getFileName() const {
        return filename;
    }
};

#endif

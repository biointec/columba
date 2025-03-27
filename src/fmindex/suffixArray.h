#ifndef SUFFIXARRAY_H
#define SUFFIXARRAY_H

#include "../bitvec.h"
#include "../definitions.h"
#include "../logger.h"

#include <cstdint>
#include <fstream>
#include <iostream> // used for printing
#include <string>
#include <vector>

#ifdef _WIN32
#include <windows.h>
#else
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>
#endif

/**
 * @class SparseSuffixArray
 * @brief Represents a sparse suffix array with support for memory-mapped file
 * reading.
 */
class SparseSuffixArray {
  private:
    length_t sparseNessFactor; /**< Sparsity factor of the suffix array */
    Bitvec bitvector; /**< Bitvector to mark positions in the suffix array */
    std::vector<length_t>
        sparseSA; /**< Sparse suffix array data (only used in index building)*/

    // Memory-mapped file related members
#ifdef _WIN32
    HANDLE hFile = INVALID_HANDLE_VALUE; /**< Handle to the file (Windows) */
    HANDLE hMapFile = NULL; /**< Handle to the file mapping object (Windows) */
#else
    int fd = -1; /**< File descriptor (Unix) */
#endif
    length_t* mappedData = nullptr; /**< Pointer to the mapped data */
    size_t mappedSize = 0;          /**< Size of the mapped data */

    // Function to pre-fault memory-mapped file
    static void preFaultMappedData(length_t* mappedData, size_t size) {
        logger.logDeveloper("Pre-faulting memory-mapped data");
        for (size_t i = 0; i < size; ++i) {
            // Accessing each element to load it into memory
            volatile length_t temp = mappedData[i];
            (void)temp; // Prevents optimization that could skip this access
        }
    }

    /**
     * @brief Reads the suffix array from a file using memory mapping.
     * @param name The name of the file to read from.
     */
    void readArrayWithMemoryMapping(const std::string& name) {
#ifdef _WIN32
        hFile = CreateFile(name.c_str(), GENERIC_READ, 0, NULL, OPEN_EXISTING,
                           FILE_ATTRIBUTE_NORMAL, NULL);
        if (hFile == INVALID_HANDLE_VALUE) {
            throw std::runtime_error("Cannot open file: " + name);
        }

        hMapFile = CreateFileMapping(hFile, NULL, PAGE_READONLY, 0, 0, NULL);
        if (hMapFile == NULL) {
            CloseHandle(hFile);
            throw std::runtime_error("Cannot create file mapping: " + name);
        }

        mappedData = (length_t*)MapViewOfFile(hMapFile, FILE_MAP_READ, 0, 0, 0);
        if (mappedData == NULL) {
            CloseHandle(hMapFile);
            CloseHandle(hFile);
            throw std::runtime_error("Cannot map view of file: " + name);
        }

        LARGE_INTEGER fileSize;
        if (!GetFileSizeEx(hFile, &fileSize)) {
            UnmapViewOfFile(mappedData);
            CloseHandle(hMapFile);
            CloseHandle(hFile);
            throw std::runtime_error("Cannot get file size: " + name);
        }

        mappedSize = fileSize.QuadPart / sizeof(length_t);
#else
        fd = open(name.c_str(), O_RDONLY);
        if (fd == -1) {
#ifdef DEVELOPER_MODE
            throw std::runtime_error("Cannot open file: " + name);
#else
            throw std::runtime_error("Problem reading file: " + name);
#endif
        }

        off_t fileSize = lseek(fd, 0, SEEK_END);
        if (fileSize == -1) {
            close(fd);
#ifdef DEVELOPER_MODE
            throw std::runtime_error("Cannot get file size: " + name);
#else
            throw std::runtime_error("Problem reading file: " + name);
#endif
        }

        mappedData =
            (length_t*)mmap(NULL, fileSize, PROT_READ, MAP_PRIVATE, fd, 0);
        if (mappedData == MAP_FAILED) {
            close(fd);
#ifdef DEVELOPER_MODE
            throw std::runtime_error("Cannot map file: " + name);
#else
            throw std::runtime_error("Problem reading file: " + name);
#endif
        }

        mappedSize = fileSize / sizeof(length_t);
#endif
        // Pre-faulting all elements
        preFaultMappedData(mappedData, mappedSize);
    }

  public:
    /**
     * @brief Overloaded subscript operator to access the bitvector.
     * @param i Index to access.
     * @return True if the bitvector at index i is set, false otherwise.
     */
    bool inline operator[](const length_t i) const {
        return bitvector[i];
    }

    /**
     * @brief Gets the value from the suffix array using the appropriate
     * implementation.
     * @param i Index to retrieve.
     * @return The value from the suffix array.
     */
    length_t inline get(const length_t i) const {
        assert(bitvector[i]);
#ifdef BUILD_INDEX
        return sparseSA[bitvector.rank(i)];
#else
        return mappedData[bitvector.rank(i)];
#endif
    }

    /**
     * @brief Constructor for creating a sparse suffix array from an
     * existing suffix array.
     * @param sa The original suffix array.
     * @param sparseNessFactor The sparsity factor to use.
     */
    SparseSuffixArray(const std::vector<length_t>& sa,
                      const length_t sparseNessFactor)
        : sparseNessFactor(sparseNessFactor), bitvector(sa.size()) {

        sparseSA.reserve(sa.size() / sparseNessFactor);
        for (length_t i = 0; i < sa.size(); i++) {
            const auto& el = sa[i];
            if (el % sparseNessFactor == 0) {
                sparseSA.emplace_back(el);
                bitvector[i] = true;
            }
        }
        bitvector.index();
    }

    /**
     * @brief Constructor for reading a sparse suffix array from file
     * @param basename The base name of the file to read.
     * @param sparseNess The sparsity factor to use.
     */
    SparseSuffixArray(const std::string& basename, const length_t sparseNess)
        : sparseNessFactor(sparseNess) {

        using namespace std;
        stringstream ss;
        {
            ss << "Reading "
               << basename + ".sa.bv." + std::to_string(sparseNessFactor)
               << "...";
            logger.logInfo(ss);
            auto name = basename + ".sa.bv." + std::to_string(sparseNessFactor);
            std::ifstream ifs(name, std::ios::binary);
            if (!ifs) {
                throw std::runtime_error(
                    "Cannot open file: " + name +
                    ". Did you set an incorrect suffix array sparseness "
                    "factor using the -s flag or move your index files?");
            }
            bitvector.read(ifs);
        }
        {
            ss << "Reading "
               << basename + ".sa." + std::to_string(sparseNessFactor) << "..";
            logger.logInfo(ss);
            auto name = basename + ".sa." + std::to_string(sparseNessFactor);

            readArrayWithMemoryMapping(name);
        }
    }

    /**
     * @brief Destructor to clean up memory-mapped resources.
     */
    ~SparseSuffixArray() {
#ifdef _WIN32
        if (mappedData)
            UnmapViewOfFile(mappedData);
        if (hMapFile)
            CloseHandle(hMapFile);
        if (hFile != INVALID_HANDLE_VALUE)
            CloseHandle(hFile);
#else
        if (mappedData)
            munmap(mappedData, mappedSize * sizeof(length_t));
        if (fd != -1)
            close(fd);
#endif
    }

    /**
     * @brief Writes the sparse suffix array to files.
     * @param basename The base name of the files to write.
     */
    void write(const std::string& basename) const {
        {
            std::ofstream ofs(basename + ".sa.bv." +
                                  std::to_string(sparseNessFactor),
                              std::ios::binary);
            bitvector.write(ofs);
        }
        {
            std::ofstream ofs(basename + ".sa." +
                                  std::to_string(sparseNessFactor),
                              std::ios::binary);
            ofs.write((char*)sparseSA.data(),
                      sparseSA.size() * sizeof(length_t));
        }
    }

    length_t inline getFactor() const {
        return sparseNessFactor;
    }
};

#endif

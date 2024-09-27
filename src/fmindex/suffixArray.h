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
    std::vector<length_t> sparseSA; /**< Sparse suffix array data */

    /**
     * @brief Function pointer to access the suffix array, based on if the data
     * is memory mapped or not.
     */
    length_t (SparseSuffixArray::*getImpl)(const length_t) const;

    void readArray(const std::string& name) {

        std::ifstream ifs(name, std::ios::binary);
        if (!ifs) {
            throw std::runtime_error("Cannot open file: " + name);
        }
        ifs.seekg(0, std::ios::end);
        sparseSA.resize(ifs.tellg() / sizeof(length_t));
        ifs.seekg(0, std::ios::beg);
        ifs.read((char*)&sparseSA[0], sparseSA.size() * sizeof(length_t));
    }

  public:
    /**
     * @brief Overloaded subscript operator to access the bitvector.
     * @param i Index to access.
     * @return True if the bitvector at index i is set, false otherwise.
     */
    bool operator[](const length_t i) const {
        return bitvector[i];
    }

    /**
     * @brief Gets the value from the suffix array using the appropriate
     * implementation.
     * @param i Index to retrieve.
     * @return The value from the suffix array.
     */
    length_t get(const length_t i) const {
        assert(bitvector[i]);
        return sparseSA[bitvector.rank(i)];
    }

    /**
     * @brief Constructor for creating a sparse suffix array from an existing
     * suffix array. This wil set the getImpl function to use the vector with
     * existing data.
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
     * @brief Constructor for reading a sparse suffix array from file using
     * memory mapping. This will set the getImpl function to use the memory
     * mapped data.
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
            readArray(name);
        }
    }

    /**
     * @brief Destructor to clean up memory-mapped resources.
     */
    ~SparseSuffixArray() {
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
};

#endif

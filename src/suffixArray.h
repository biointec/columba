#ifndef SUFFIXARRAY_H
#define SUFFIXARRAY_H

#include "bitvec.h"
#include <fstream>
#include <iostream> // used for printing
#include <stdint.h>
#include <string>
#include <vector>

typedef uint32_t length_t;

class SparseSuffixArray {
  private:
    length_t sparseNessFactor;
    Bitvec bitvector;
    std::vector<length_t> sparseSA;

  public:
    bool operator[](const length_t i) const {
        return bitvector[i];
    }

    length_t get(const length_t i) const {
        assert(bitvector[i]);
        return sparseSA[bitvector.rank(i)];
    }

    SparseSuffixArray(const std::vector<length_t>& sa,
                      const length_t sparseNessFactor)
        : sparseNessFactor(sparseNessFactor) {
        bitvector = Bitvec(sa.size());
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

    SparseSuffixArray(const std::string& basename, const length_t sparseNess)
        : sparseNessFactor(sparseNess) {

        using namespace std;
        {

            cout << "Reading "
                 << basename + ".sa.bv." + std::to_string(sparseNessFactor)
                 << "...";
            cout.flush();
            auto name = basename + ".sa.bv." + std::to_string(sparseNessFactor);
            std::ifstream ifs(name);
            if (!ifs) {
                throw std::runtime_error("Cannot open file: " + name);
            }
            bitvector.read(ifs);
            cout << "done " << endl;
        }
        {
            cout << "Reading "
                 << basename + ".sa." + std::to_string(sparseNessFactor)
                 << "..";
            cout.flush();
            auto name = basename + ".sa." + std::to_string(sparseNessFactor);
            std::ifstream ifs(name);
            if (!ifs) {
                throw std::runtime_error("Cannot open file: " + name);
            }
            ifs.seekg(0, std::ios::end);
            sparseSA.resize(ifs.tellg() / sizeof(length_t));
            ifs.seekg(0, std::ios::beg);
            ifs.read((char*)&sparseSA[0], sparseSA.size() * sizeof(length_t));
            cout << "done" << endl;
        }
    }

    void write(const std::string& basename) const {
        {
            std::ofstream ofs(basename + ".sa.bv." +
                              std::to_string(sparseNessFactor));
            bitvector.write(ofs);
        }
        {
            std::ofstream ofs(basename + ".sa." +
                              std::to_string(sparseNessFactor));
            ofs.write((char*)sparseSA.data(),
                      sparseSA.size() * sizeof(length_t));
        }
    }
};

#endif
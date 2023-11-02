/******************************************************************************
 *  Columba 1.2: Approximate Pattern Matching using Search Schemes            *
 *  Copyright (C) 2020-2023 - Luca Renders <luca.renders@ugent.be> and        *
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

#ifndef SUFFIXARRAY_H
#define SUFFIXARRAY_H

#include "bitvec.h"
#include <fstream>
#include <iostream> // used for printing
#include <stdint.h>
#include <string>
#include <vector>

#include "wordlength.h"

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
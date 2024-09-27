/******************************************************************************
 *  This file is a slightly modified version of the sparse_sd_vector class    *
 *  from the r-index project.                                                 *
 *  Copyright (c) 2017 - Nicola Prezza <nicola.prezza@gmail.com>              *
 *                                                                            *
 *  MIT License                                                               *
 *                                                                            *
 *  Permission is hereby granted, free of charge, to any person obtaining a   *
 *  copy of this software and associated documentation files (the "Software"),*
 *  to deal in the Software without restriction, including without limitation *
 *  the rights to use, copy, modify, merge, publish, distribute, sublicense,  *
 *  and/or sell copies of the Software, and to permit persons to whom the     *
 *  Software is furnished to do so, subject to the following conditions:      *
 *                                                                            *
 *  The above copyright notice and this permission notice shall be included   *
 *  in all copies or substantial portions of the Software.                    *
 *                                                                            *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS   *
 *  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF                *
 *  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN *
 *  NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,  *
 *  DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR     *
 *  OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE *
 *  USE OR OTHER DEALINGS IN THE SOFTWARE.                                    *
 ******************************************************************************/

#ifndef SPARSEBITVEC_HPP
#define SPARSEBITVEC_HPP

#include <sdsl/sd_vector.hpp>

class SparseBitvec {
  private:
    size_t N; // bitvector length

    sdsl::sd_vector<> sdv;
    sdsl::sd_vector<>::rank_1_type rank1;
    sdsl::sd_vector<>::select_1_type select1;

  public:
    /*
     * empty constructor. Initialize bitvector with length 0.
     */
    SparseBitvec() : N(0) {
    }

    /*
     * constructor. build bitvector given a vector of bools
     */
    SparseBitvec(std::vector<bool>& b) {

        N = b.size();

        if (b.size() == 0) {
            return;
        }

        sdsl::bit_vector bv(b.size());

        for (uint64_t i = 0; i < b.size(); ++i) {
            bv[i] = b[i];
        }

        sdv = sdsl::sd_vector<>(bv);
        rank1 = sdsl::sd_vector<>::rank_1_type(&sdv);
        select1 = sdsl::sd_vector<>::select_1_type(&sdv);
    }

    /*
     * constructor. build bitvector given a bit_vector
     */
    SparseBitvec(sdsl::bit_vector& bv) {

        sdv = sdsl::sd_vector<>(bv);
        rank1 = sdsl::sd_vector<>::rank_1_type(&sdv);
        select1 = sdsl::sd_vector<>::select_1_type(&sdv);
    }

    SparseBitvec& operator=(const SparseBitvec& other) {

        N = other.sdv.size();
        sdv = sdsl::sd_vector<>(other.sdv);
        rank1 = sdsl::sd_vector<>::rank_1_type(&sdv);
        select1 = sdsl::sd_vector<>::select_1_type(&sdv);

        return *this;
    }

    /**
     * @param position i in the bitvector
     * @return bit in position i
     * only access! the bitvector is static.
     */
    bool operator[](size_t i) const {
        assert(i < N);
        return sdv[i];
    }

    /**
     * @param i position in the bitvector
     * @return number of bits equal to 1 before position i excluded
     */
    size_t rank(size_t i) const {
        assert(i <= N);
        return rank1(i);
    }

    /**
     * @param i position in the bitvector
     * @return the predecessor of i (position i excluded) in
     * rank space (apply select to get bitvector space)
     */
    size_t predecessorRank(size_t i) const {
        // i must have a predecessor
        assert(rank(i) > 0);
        return rank(i) - 1;
    }

    /**
     * @param i position in the bitvector
     * @return predecessor of i (i excluded) in
     * bitvector space
     */
    size_t predecessor(size_t i) const {
        // i must have a predecessor
        assert(rank(i) > 0);
        return select(rank(i) - 1);
    }

    /**
     * @param position in the bitvector
     * @return the rank of predecessor of i (i excluded) in
     * bitvector space. If i does not have a predecessor,
     * return rank of the last bit set in the bitvector
     */
    size_t predecessorRankCircular(size_t i) const {
        return rank(i) == 0 ? rank(N) - 1 : rank(i) - 1;
    }

    /**
     * retrieve length of the i-th gap (i>=0). gap length includes the leading 1
     * @param i i < number_of_1()
     *
     */
    size_t gapAt(size_t i) const {
        assert(i < rank(N));

        if (i == 0) {
            return select(0) + 1;
        }
        return select(i) - select(i - 1);
    }

    /**
     * @param i i>0
     * @return the position of the i-th one in the bitvector. i starts from 0!
     */
    size_t select(size_t i) const {
        assert(i < rank(N));
        return select1(i + 1); // in sd_vector, i starts from 1
    }

    /**
     * @return the size of the bitvector
     */
    size_t size() const {
        return N;
    }

    /**
     * @return the number of bits set to 1 in the bitvector
     */
    size_t numberOf1() const {
        return rank(N);
    }

    /**
     * serialize the structure to the ostream
     * @param out the output stream
     * @return number of bytes written to the output stream
     */
    size_t serialize(std::ostream& out) const {

        size_t w_bytes = 0;

        out.write((char*)&N, sizeof(N));

        w_bytes += sizeof(N);

        if (N == 0)
            return w_bytes;

        w_bytes += sdv.serialize(out);

        return w_bytes;
    }

    /**
     * load the structure from the istream
     * @param in the input stream
     */
    void load(std::istream& in) {

        in.read((char*)&N, sizeof(N));

        if (N == 0)
            return;

        sdv.load(in);
        rank1 = sdsl::sd_vector<>::rank_1_type(&sdv);
        select1 = sdsl::sd_vector<>::select_1_type(&sdv);
    }

    /**
     * write the structure to a file
     * @param filename name of the output file
     * @return number of bytes written to the file
     */
    size_t write(const std::string& filename) const {
        std::ofstream ofs(filename);
        if (!ofs) {
            throw std::runtime_error("Cannot open file: " + filename);
        }
        size_t w_bytes = serialize(ofs);
        ofs.close();
        return w_bytes;
    }

    /**
     * read the structure from a file
     * @param filename name of the input file
     * @return true if succeeded, false otherwise
     */
    bool read(const std::string& filename) {
        std::ifstream ifs(filename);
        if (!ifs) {
            return false;
        }
        load(ifs);
        ifs.close();
        return true;
    }
};

#endif /* SPARSEBITVEC_HPP */

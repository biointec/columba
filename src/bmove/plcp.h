/******************************************************************************
 *  This file includes code inspired by the full-br-index repository          *
 *  Original source: https://github.com/U-Ar/full-br-index                    *
 *  Author: Yuma Arakawa                                                      *
 *                                                                            *
 *  This program is free software: you can redistribute it and/or modify      *
 *  it under the terms of the GNU Affero General Public License as            *
 *  published by the Free Software Foundation, either version 3 of the        *
 *  License, or (at your option) any later version.                           *
 *                                                                            *
 *  This program is distributed in the hope that it will be useful,           *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the              *
 *  GNU Affero General Public License for more details.                       *
 *                                                                            *
 *  You should have received a copy of the GNU Affero General Public License  *
 *  along with this program. If not, see <https://www.gnu.org/licenses/>.     *
 ******************************************************************************/


#ifndef PLCP_H
#define PLCP_H

#include "../definitions.h"
#include "../logger.h"
#include "sparsebitvec.h"

#include <fstream>  // used for reading in files
#include <iostream> // used for printing
#include <cstdint>

/**
 * class for the PLCP (permuted longest common prefix) array
*/
class PLCP {
    private:
        // length of PLCP
        size_t n = 0;

        // run length encoded PLCP representation
        SparseBitvec ones;
        SparseBitvec zeros;

    public:

        /**
         * Default constructor
        */
        PLCP() {}

        /**
         * Constructor
         * Constructs the PLCP array using Kasai's algorithm and 
         * compresses it to a run-length encoded representation
         * @param T the text
         * @param SA the suffix array
         * @param verbose Print info to cout, default false
        */
        PLCP(const std::string& T, const std::vector<length_t>& SA, bool verbose = false) {
            n = SA.size();

            // construct inverse suffix array
            std::vector<length_t> ISA(n, 0);
            for (length_t i = 0; i < n; i++) {
                ISA[SA[i]] = i;
            }

            // compute S = unary encoding of j + PLCP[j] in single loop without constructing the PLCP array explicitly
            // this is based on Kasai's algorithm 
            std::vector<bool> S(2*n+1, false);
            S[0] = true;
            length_t plcp_prev = 1;
            length_t pos = 0;

            length_t k = 0;
            if (verbose) {
                std::cout << "\tPLCP:";
            }
            for (length_t i = 0; i < n - 1; i++) {
                length_t j = SA[ISA[i]-1];
                while (T[i+k] == T[j+k]) {
                    k++;
                }
                // PLCP[i] = k
                pos += k + 2 - plcp_prev;
                S[pos] = true;
                plcp_prev = k;
                if (verbose) {
                    std::cout << " " << k;
                }
                if (k > 0) {
                    k--;
                }
            }
            // PLCP[n-1] = 0
            pos += 2 - plcp_prev;
            S[pos] = true;
            if (verbose) {
                std::cout << " 0\n";
            }

            // create run length encoding of PLCP by using sparse bitvectors
            // -> generate bitvectors 'ones' & 'zeros' from S
            size_t u = pos + 1;
            n = n + 2;
            std::vector<bool> ones_bv(n, false);
            std::vector<bool> zeros_bv(n, false);
            bool bit_1 = true; // tell if we are counting 1 bits
            length_t cont_1 = 0, cont_0 = 0; // continuous number of 1 bits & 0 bits

            for (length_t i = 0; i < u; i++) {
                if (S[i]) {
                    if (bit_1) {
                        cont_1++;
                    } else {
                        bit_1 = true;
                        cont_1++;
                        zeros_bv[cont_0-1] = true;
                    }
                } else {
                    if (bit_1) {
                        bit_1 = false;
                        cont_0++;
                        ones_bv[cont_1-1] = true;
                    } else {
                        cont_0++;
                    }
                }
            }
            if (bit_1) {
                ones_bv[cont_1-1] = true;
            } else {
                zeros_bv[cont_0-1] = true;
            }

            ones = SparseBitvec(ones_bv);
            zeros = SparseBitvec(zeros_bv);
        }

        /*
         * constructor receiving positional vectors of One-Zero encoding
         * directly
         */
        PLCP(length_t _n, std::vector<length_t> const& O,
             std::vector<length_t> const& Z) {
            this->n = _n + 2;

            std::vector<bool> ones_bv(n, false);
            std::vector<bool> zeros_bv(n, false);
            for (length_t i = 0; i < O.size(); i++) {
                assert(O[i] <= n);
                ones_bv[O[i]] = true;
            }
            for (length_t i = 0; i < Z.size(); i++) {
                assert(O[i] <= n);
                zeros_bv[Z[i]] = true;
            }
            ones = SparseBitvec(ones_bv);
            zeros = SparseBitvec(zeros_bv);
        }

        /**
         * Access operator
         * @param i index
         * @return PLCP value at index i
        */
        length_t operator[](length_t i) const {
            assert(i < n);
            length_t rank_0 = ones.rank(i + 1);
            if (rank_0 > 0) {
                return zeros.select(rank_0 - 1) - i + 1;
            }
            return 0;
        }

        /**
         * Write PLCP to output stream
         * @param os Output stream
         */
        void write(std::ostream& os) {
            os.write((char*)&n, sizeof(n));
            ones.serialize(os);
            zeros.serialize(os);
        }

        /**
         * Load PLCP from input stream
         * @param is Input stream
         * @return true if succeeded, false otherwise
         */
        bool read(std::istream& is) {
            is.read((char*)&n, sizeof(n));
            ones.load(is);
            zeros.load(is);
            return true;
        }

        /**
         * Write PLCP to disk
         * @param filename File name
         */
        void write(const std::string& filename) {
            std::ofstream ofs(filename);
            if (!ofs) {
                throw std::runtime_error("Cannot open file: " + filename);
            }
            write(ofs);
        }

        /**
         * Load PLCP from disk
         * @param filename File name
         * @return true if succeeded, false otherwise
         */
        bool read(const std::string& filename) {
            std::ifstream ifs(filename);
            if (!ifs) {
                throw false;
            }
            return read(ifs);
        }

        length_t printMemSize() {
            // n
            size_t totalSize = sizeof(size_t);

            std::ofstream devnull("/dev/null");
            totalSize += ones.serialize(devnull);
            totalSize += zeros.serialize(devnull);
            std::cout << "total PLCP size: " << totalSize << "\n";
            return totalSize;
        }
};


#endif // PLCP_H
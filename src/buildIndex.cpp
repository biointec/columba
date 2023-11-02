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

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

#include "bwtrepr.h"
#include "encodedtext.h"
#include "suffixArray.h"

using namespace std;

#include "wordlength.h"

void showUsage() {
    cout << "Usage: ./fmidx-build <base filename>\n\n";
    cout << "Following files are required:\n";
    cout << "\t<base filename>.txt: input text T\n";
    cout << "\t<base filename>.sa: suffix array of T\n";
    cout << "\t<base filename>.rev.sa: suffix array of reverse of T\n\n";

    cout << "Report bugs to jan.fostier@ugent.be" << endl;
}

bool parseArguments(int argc, char* argv[], string& baseFN) {
    if (argc != 2)
        return false;

    baseFN = argv[1];
    return true;
}

void readText(const string& filename, string& buf) {
    ifstream ifs(filename);
    if (!ifs)
        throw runtime_error("Cannot open file: " + filename);

    ifs.seekg(0, ios::end);
    buf.resize(ifs.tellg());
    ifs.seekg(0, ios::beg);
    ifs.read((char*)buf.data(), buf.size());
}

void readSATextMode(const string& filename, vector<length_t>& sa,
                    size_t saSizeHint) {
    ifstream ifs(filename);
    if (!ifs)
        throw runtime_error("Cannot open file: " + filename);

    sa.reserve(saSizeHint);
    length_t el;
    while (ifs >> el)
        sa.push_back(el);
}

void readSA(const string& filename, vector<length_t>& sa, size_t saSizeHint) {
    ifstream ifs(filename, ios::binary);
    if (!ifs)
        throw runtime_error("Cannot open file: " + filename);

    ifs.seekg(0, ios::end);
    size_t numElements = ifs.tellg() / sizeof(length_t);

    if (numElements == saSizeHint) { // file is likely binary
        sa.resize(ifs.tellg() / sizeof(length_t));
        ifs.seekg(0, ios::beg);
        ifs.read((char*)sa.data(), sa.size() * sizeof(length_t));
    } else { // try to read SA in text mode
        readSATextMode(filename, sa, saSizeHint);
    }
}

void sanityCheck(const string& T, vector<length_t>& sa) {
    // check T for correctness
    if (T.back() == '\n')
        throw runtime_error("T should end with a \'$\' character, "
                            "not with a newline");

    if (T.back() != '$')
        throw runtime_error("T should end with a \'$\' character");

    if (sa.size() != T.size())
        throw runtime_error("Text and suffix array contain a "
                            "different number of elements");

    // briefly check the suffix array
    length_t min = *min_element(sa.begin(), sa.end());
    length_t max = *max_element(sa.begin(), sa.end());

    if (min == 1 && max == T.size()) { // rebase to [0..T.size()-1]
        for (auto& el : sa)
            el--;
        min--;
        max--;
    }

    if (min != 0 || max != T.size() - 1)
        throw runtime_error("Suffix array must contain numbers between "
                            "[0 and " +
                            to_string(T.size() - 1) + "]");

    // check if all numbers in the suffix array are present
    Bitvec bv(sa.size());
    for (length_t i : sa)
        bv[i] = true;

    for (size_t i = 0; i < bv.size(); i++)
        if (!bv[i])
            throw runtime_error("Suffix " + to_string(i) +
                                " seems "
                                "to be missing from suffix array");

    // extra check:
    //      we could check T to see if the SA correctly sorts suffixes of T
}

void createFMIndex(const string& baseFN) {
    // read the text file from disk
    cout << "Reading " << baseFN << ".txt..." << endl;
    string T;
    readText(baseFN + ".txt", T);

    // count the frequency of each characters in T
    vector<length_t> charCounts(256, 0);
    for (char c : T)
        charCounts[(unsigned char)c]++;

    // count the number of unique characters in T
    int nUniqueChar = 0;
    for (length_t count : charCounts)
        if (count > 0)
            nUniqueChar++;

    cout << "\tText has length " << T.size() << "\n";
    cout << "\tText has " << nUniqueChar << " unique characters\n";

    if (nUniqueChar > ALPHABET) {
        cerr << "FATAL ERROR: the number of unique characters in the "
             << "text exceeds the alphabet size. Please recompile"
             << "Columba using a higher value for ALPHABET " << endl;
        exit(EXIT_FAILURE);
    }

    if (nUniqueChar < ALPHABET) {
        cout << "WARNING: the number of unique characters in the "
             << "text is less than the ALPHABET size specified when "
             << "Columba was compiled. Performance may be affected\n";
    }

    Alphabet<ALPHABET> sigma(charCounts);

    // read the suffix array
    cout << "Reading " << baseFN << ".sa..." << endl;
    vector<length_t> SA;
    readSA(baseFN + ".sa", SA, T.size());

    // perform a sanity check on the suffix array
    cout << "\tPerforming sanity checks..." << endl;
    sanityCheck(T, SA);
    cout << "\tSanity checks OK" << endl;

    // build the BWT
    cout << "Generating BWT..." << endl;
    string BWT(T.size(), '\0');
    for (size_t i = 0; i < SA.size(); i++)
        if (SA[i] > 0)
            BWT[i] = T[SA[i] - 1];
        else
            BWT[i] = T.back();

    // Encode bwt
    cout << "Encoding BWT.." << endl;
    EncodedText<ALPHABET> eBWT(sigma, BWT);

    eBWT.write(baseFN + ".bwt");

    cout << "Wrote file " << baseFN << ".bwt\n";

    // write the character counts table
    {
        ofstream ofs(baseFN + ".cct", ios::binary);
        ofs.write((char*)charCounts.data(),
                  charCounts.size() * sizeof(length_t));
        ofs.close();
    }

    cout << "Wrote file " << baseFN << ".cct\n";

    // create sparse suffix arrays
    for (int saSF = 1; saSF <= 128; saSF *= 2) {
        SparseSuffixArray sparseSA(SA, saSF);
        sparseSA.write(baseFN);
        cout << "Wrote sparse suffix array with factor " << saSF << endl;
    }
    SA.clear();

    // create succint BWT bitvector table
    BWTRepr<ALPHABET> fwdBWT(sigma, BWT);
    fwdBWT.write(baseFN + ".brt");
    cout << "Wrote file: " << baseFN << ".brt" << endl;

    BWT.clear();

    // read the reverse suffix array
    cout << "Reading " << baseFN << ".rev.sa..." << endl;
    vector<length_t> revSA;
    readSA(baseFN + ".rev.sa", revSA, T.size());

    // perform a sanity check on the suffix array
    cout << "\tPerforming sanity checks..." << endl;
    sanityCheck(T, revSA);
    cout << "\tSanity checks OK" << endl;

    // build the reverse BWT
    string rBWT(T.size(), '\0');
    rBWT.resize(T.size());
    for (size_t i = 0; i < revSA.size(); i++)
        if (revSA[i] > 0)
            rBWT[i] = T[T.size() - revSA[i]];
        else
            rBWT[i] = T.front();
    revSA.clear();

    // create succint reverse BWT bitvector table
    BWTRepr<ALPHABET> revBWT(sigma, rBWT);
    revBWT.write(baseFN + ".rev.brt");
    cout << "Wrote file: " << baseFN << ".rev.brt" << endl;
}

int main(int argc, char* argv[]) {
    string baseFN;

    if (!parseArguments(argc, argv, baseFN)) {
        showUsage();
        return EXIT_FAILURE;
    }

    cout << "Welcome to Columba's index construction!\n";
    cout << "Alphabet size is " << ALPHABET - 1 << " + 1\n";

    try {
        createFMIndex(baseFN);
    } catch (const std::exception& e) {
        cerr << "Fatal error: " << e.what() << endl;
        return EXIT_FAILURE;
    }

    cout << "Exiting... bye!" << endl;
    return EXIT_SUCCESS;
}

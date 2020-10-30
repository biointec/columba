/******************************************************************************
 *  Columba: Approximate Pattern Matching using Search Schemes                *
 *  Copyright (C) 2020 - Luca Renders <luca.renders@ugent.be> and             *
 *                       Jan Fostier <jan.fostier@ugent.be>                   *
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
#include <cstdlib>
#include <fstream>
#include <iostream>

#include "customtypedefs.h"

void printHelp() {

    std::cout << "Usage: ./columba-build <base filename>\n\n";
    std::cout << "Following files are required:\n";
    std::cout << "\t<base filename>.txt: input text T\n";
    std::cout << "\t<base filename>.sa: suffix array of T\n";
    std::cout << "\t<base filename>.rev.sa: suffix array of reverse of T\n\n";

    std::cout << "Report bugs to luca.renders@ugent.be" << std::endl;
}

void createBWT(const std::string& text, const std::vector<length_t>& sa,
               std::string& bwt, std::vector<length_t>& counts,
               std::vector<unsigned char>& indexToChar,
               std::vector<int>& charToIndex) {

    std::vector<length_t> charToCount(256, 0);
    bwt.reserve(sa.size());

    for (auto saIt = sa.cbegin(); saIt != sa.cend(); saIt++) {
        // get the previous character
        char c;
        if ((*saIt) != 0) {
            c = text[*saIt - 1];
        } else {
            // if the entry in the suffix array is zero, then the final
            // character is the special end character
            c = '$';
        }

        // append this char to the bwt
        bwt += c;

        // increase count of c
        charToCount[c]++;
    }
    // create the alphabet
    std::vector<int> allChars(256, -1);
    length_t charNumber = 0;
    for (int c = 0; c < 256; c++) {
        if (charToCount[(unsigned char)c] > 0) {
            allChars[(unsigned char)c] = charNumber++;
            indexToChar.push_back(c);
        }
    }

    charToIndex.swap(allChars);
    counts.swap(charToCount);
}

void createIndex(const std::string& name, const std::string& output) {

    std::cout << "Read in the textfile " << name << ".txt" << std::endl;
    std::string text = readString(name + ".txt");
    if (text[text.size() - 1] == '\n') {
        text.erase(text.end() - 1);
    }
    std::cout << "Textsize: " << text.size() << std::endl;
    std::vector<length_t> sa;
    bool dollarPresent = text[text.size() - 1] == '$';
    if (!dollarPresent) {
        throw std::runtime_error("No sentinel ($) at end of file");
    }

    std::cout << "Read in suffix array" << std::endl;
    try {
        readArray(name + ".sa", text.size(), sa);
    } catch (const std::exception& e) {
        // retry with .1 extensions
        readArray(name + ".sa.1", text.size(), sa);
    }

    std::cout << "Size of SA: " << sa.size() << std::endl;
    if (sa.size() != text.size()) {
        throw std::runtime_error(
            "The suffix array and text do not have the same size!");
    }

    std::cout << "Creating alphabet and BWT" << std::endl;
    std::string bwt;

    std::vector<size_t> charToCount(256, 0);
    bwt.reserve(sa.size());

    for (auto saIt = sa.cbegin(); saIt != sa.cend(); saIt++) {
        // get the previous character
        char c;
        if ((*saIt) != 0) {
            c = text[*saIt - 1];
        } else {
            // if the entry in the suffix array is zero, then the final
            // character is the special end character
            c = '$';
        }

        // append this char to the bwt
        bwt += c;

        // increase count of c
        charToCount[c]++;
    }

    // initialize the alphabet
    Alphabet<ALPHABET> sigma(charToCount);

    std::cout << "Alphabet contains " << sigma.size() << " characters"
              << std::endl;

    std::cout << "Writing counts to disk" << std::endl;
    std::ofstream ofs(output + ".cct", std::ios::binary);

    ofs.write((char*)&charToCount[0], charToCount.size() * sizeof(size_t));
    ofs.close();

    std::cout << "Writing BWT to disk" << std::endl;
    ofs = std::ofstream(output + ".bwt");
    ofs.write((char*)&bwt[0], bwt.size());
    ofs.close();

    // create sparse suffix arryas
    std::cout << "Create sparse versions of the suffix array" << std::endl;
    for (int saSF = 1; saSF <= 256; saSF *= 2) {
        std::vector<length_t> spSA((sa.size() + saSF - 1) / saSF);

        for (size_t i = 0; i < spSA.size(); i++) {
            spSA[i] = sa[i * saSF];
        }

        ofs = std::ofstream(output + ".sa." + std::to_string(saSF));
        ofs.write((char*)&spSA[0], spSA.size() * sizeof(length_t));

        ofs.close();
        std::cout << "Wrote sparse file for factor: " << saSF << std::endl;
    }

    // no need for SA anymore
    sa.clear();

    // create occ en cumocc tab
    std::cout << "Writing occ en cumOcc tables" << std::endl;
    BWTRepr<5> fwdRepr(sigma, bwt);

    fwdRepr.write(output + ".brt");

    std::cout << "Read reverse sa" << std::endl;
    std::vector<length_t> revSA;
    readArray(name + ".rev.sa", text.size(), revSA);
    if (revSA.size() != text.size()) {
        throw std::runtime_error(
            "The reverse suffix array and text do no have the same size!");
    }
    std::cout << "Build reverse BWT" << std::endl;
    std::string revBWT;
    revBWT.reserve(text.size());

    for (auto saIt = revSA.begin(); saIt != revSA.end(); saIt++) {
        // get the previous character
        char c;
        if ((*saIt) != 0) {
            c = text[text.size() - *saIt];
        } else {
            // if the entry in the suffix array is zero, take the first
            // character of the text (= final of rev text)
            c = text[0];
        }

        // append this char to the reversed bwt
        revBWT += c;
    }

    std::cout << "Writing reversed occurrences" << std::endl;
    BWTRepr<ALPHABET> revRepr(sigma, revBWT);

    revRepr.write(output + ".rev.brt");
}

int main(int argc, char* argv[]) {

    if (argc != 2 && argc != 3) {
        std::cout << argc;
        printHelp();
        return EXIT_FAILURE;
    }

    std::string name = argv[1];
    std::string output = name;
    if (argc == 3) {
        output = argv[2];
    }

    try {
        createIndex(name, output);
    } catch (const std::exception& e) {
        std::cerr << "Error while creating index: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    std::cout << "successs" << std::endl;
    return EXIT_SUCCESS;
}
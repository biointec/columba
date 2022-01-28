/******************************************************************************
 *  Columba: Approximate Pattern Matching using Search Schemes                *
 *  Copyright (C) 2020-2022 - Luca Renders <luca.renders@ugent.be> and        *
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
#include "searchstrategy.h"
#include <algorithm>
#include <chrono>
#include <set>
#include <string.h>

using namespace std;
vector<string> schemes = {"kuch1",  "kuch2", "kianfar", "manbest",
                          "pigeon", "01*0",  "custom",  "naive"};

int editDistDP(string P, string O, int maxED) {

    int n = (int)P.length();
    int m = (int)O.length();
    string* horizontal = &P;
    string* vertical = &O;
    if (n > m) {
        horizontal = &O;
        vertical = &P;
        int temp = n;
        n = m;
        m = temp;
    }

    // check the dimensions of s1 and s2
    if ((max(m, n) - min(m, n)) > maxED)
        return numeric_limits<int>::max();

    BandMatrix mat(m + 1 + maxED, maxED);

    // fill in the rest of the matrix
    for (int i = 1; i <= m; i++) {
        for (int j = mat.getFirstColumn(i); j <= mat.getLastColumn(i) && j <= n;
             j++) {
            mat.updateMatrix(vertical->at(i - 1) != horizontal->at(j - 1), i,
                             j);
        }
    }
    return mat(m, n);
}

string getFileExt(const string& s) {

    size_t i = s.rfind('.', s.length());
    if (i != string::npos) {
        return (s.substr(i + 1, s.length() - i));
    }

    return ("");
}

vector<pair<string, string>> getReads(const string& file) {
    vector<pair<string, string>> reads;
    reads.reserve(200000);

    const auto& extension = getFileExt(file);

    bool fasta =
        (extension == "FASTA") || (extension == "fasta") || (extension == "fa");
    bool fastq = (extension == "fq") || (extension == "fastq");

    ifstream ifile(file.c_str());
    if (!ifile) {
        throw runtime_error("Cannot open file " + file);
    }
    if (!fasta && !fastq) {
        // this is a csv file
        if (extension != "csv") {
            throw runtime_error("extension " + extension +
                                " is not a valid extension for the readsfile");
        }
        string line;
        // get the first line we do not need this
        getline(ifile, line);

        while (getline(ifile, line)) {
            istringstream iss{line};

            vector<string> tokens;
            string token;

            while (getline(iss, token, ',')) {
                tokens.push_back(token);
            }
            string position = tokens[1];
            string read =
                tokens[2]; // ED + 2 column contains a read with ED compared
                           // to the read at position position with length
            string p = position;
            reads.push_back(make_pair(p, read));
        }
    } else if (fasta) {
        // fasta file
        string read = "";
        string p = "";
        string line;
        while (getline(ifile, line)) {
            if (!line.empty() && line[0] == '>') {

                if (!read.empty()) {

                    reads.push_back(make_pair(p, read));
                    reads.push_back(
                        make_pair(p, Nucleotide::getRevCompl(read)));
                    read.clear();
                }

                p = (line.substr(1));

            } else {
                read += line;
            }
        }
        if (!read.empty()) {

            reads.push_back(make_pair(p, read));
            reads.push_back(make_pair(p, Nucleotide::getRevCompl(read)));
            read.clear();
        }
    } else {
        // fastQ
        string read = "";
        string id = "";
        string line;
        bool readLine = false;
        while (getline(ifile, line)) {
            if (!line.empty() && line[0] == '@') {
                if (!read.empty()) {

                    reads.push_back(make_pair(id, read));
                    reads.push_back(
                        make_pair(id, Nucleotide::getRevCompl(read)));
                    read.clear();
                }
                id = (line.substr(1));
                readLine = true;
            } else if (readLine) {
                read = line;
                readLine = false;
            }
        }
        if (!read.empty()) {

            reads.push_back(make_pair(id, read));
            reads.push_back(make_pair(id, Nucleotide::getRevCompl(read)));
            read.clear();
        }
    }

    return reads;
}

void writeToOutput(const string& file, const vector<vector<TextOcc>>& mPerRead,
                   const vector<pair<string, string>>& reads) {

    cout << "Writing to output file " << file << " ..." << endl;
    ofstream f2;
    f2.open(file);

    f2 << "identifier\tposition\tlength\tED\tCIGAR\treverseComplement\n";
    for (unsigned int i = 0; i < reads.size(); i += 2) {
        auto id = reads[i].first;

        for (auto m : mPerRead[i]) {
            f2 << id << "\t" << m.getOutput() << "\t0\n";
        }

        for (auto m : mPerRead[i + 1]) {
            f2 << id << "\t" << m.getOutput() << "\t1\n";
        }
    }

    f2.close();
}

double findMedian(vector<length_t> a, int n) {

    // If size of the arr[] is even
    if (n % 2 == 0) {

        // Applying nth_element
        // on n/2th index
        nth_element(a.begin(), a.begin() + n / 2, a.end());

        // Applying nth_element
        // on (n-1)/2 th index
        nth_element(a.begin(), a.begin() + (n - 1) / 2, a.end());

        // Find the average of value at
        // index N/2 and (N-1)/2
        return (double)(a[(n - 1) / 2] + a[n / 2]) / 2.0;
    }

    // If size of the arr[] is odd
    else {

        // Applying nth_element
        // on n/2
        nth_element(a.begin(), a.begin() + n / 2, a.end());

        // Value at index (N/2)th
        // is the median
        return (double)a[n / 2];
    }
}

void doBench(vector<pair<string, string>>& reads, FMIndex& mapper,
             SearchStrategy* strategy, string readsFile, length_t ED) {

    size_t totalUniqueMatches = 0, sizes = 0, mappedReads = 0;

    cout << "Benchmarking with " << strategy->getName()
         << " strategy for max distance " << ED << " with "
         << strategy->getPartitioningStrategy() << " partitioning and using "
         << strategy->getDistanceMetric() << " distance " << endl;
    cout << "Switching to in text verification at "
         << strategy->getSwitchPoint() << endl;
    cout.precision(2);

    vector<vector<TextOcc>> matchesPerRead = {};
    matchesPerRead.reserve(reads.size());

    std::vector<length_t> numberMatchesPerRead;
    numberMatchesPerRead.reserve(reads.size());

    Counters counters;

    auto start = chrono::high_resolution_clock::now();
    for (unsigned int i = 0; i < reads.size(); i += 2) {

        const auto& p = reads[i];

        auto originalPos = p.first;
        string read = p.second;
        string revCompl = reads[i + 1].second;

        if (((i >> 1) - 1) % (8192 / (1 << ED)) == 0) {
            cout << "Progress: " << i / 2 << "/" << reads.size() / 2 << "\r";
            cout.flush();
        }

        sizes += read.size();

        auto matches = strategy->matchApprox(read, ED, counters);
        totalUniqueMatches += matches.size();

        // do the same for the reverse complement
        vector<TextOcc> matchesRevCompl =
            strategy->matchApprox(revCompl, ED, counters);
        totalUniqueMatches += matchesRevCompl.size();

        // keep track of the number of mapped reads
        mappedReads += !(matchesRevCompl.empty() && matches.empty());

        matchesPerRead.emplace_back(matches);
        matchesPerRead.emplace_back(matchesRevCompl);

        numberMatchesPerRead.emplace_back(matches.size() +
                                          matchesRevCompl.size());

        // correctness check, comment this out if you want to check
        // For each reported match the reported edit distance is checked and
        // compared to a recalculated value using a single banded matrix this
        // is slow WARNING: this checks the EDIT DISTANCE, for it might be that
        // the hamming distance is higher
        /*for (auto match : matches) {

            string O = text.substr(match.getRange().getBegin(),
                                   match.getRange().getEnd() -
                                       match.getRange().getBegin());

            int trueED = editDistDP(read, O, ED);
            int foundED = match.getDistance();
            if (foundED != trueED) {
                cout << i << "\n";
                cout << "Wrong ED!!"
                     << "\n";
                cout << "P: " << read << "\n";
                cout << "O: " << O << "\n";
                cout << "true ED " << trueED << ", found ED " << foundED <<
        "\n"
                     << match.getRange().getBegin() << "\n";
            }
        }*/

        // this block checks if at least one occurrence is found and if the
        // identifier is a number and then checks if this position is found as
        // a match (for checking correctness) if you want to check if the
        // position is found as a match make sure that the identifier of the
        // read is the position. Out comment this block for the check
        /*  bool originalFound = true;
          try {
              length_t pos = stoull(originalPos);
              originalFound = false;

              for (auto match : matches) {

                  if (match.getRange().getBegin() >= pos - (ED + 2) &&
                      match.getRange().getBegin() <= pos + (ED + 2)) {
                      originalFound = true;
                      break;
                  }
              }
          } catch (const std::exception& e) {
              // nothing to do, identifier is  not the orignal position
          }

          // check if at least one occurrence was found (for reads that were
          // sampled from actual reference) Out-cooment this block if you want
          // to do this.
          if (matches.size() == 0 || (!originalFound)) {
              cout << "Could not find occurrence for " << originalPos << endl;
          }*/
    }

    auto finish = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = finish - start;
    cout << "Progress: " << reads.size() << "/" << reads.size() << "\n";
    cout << "Results for " << strategy->getName() << endl;

    cout << "Total duration: " << fixed << elapsed.count() << "s\n";
    cout << "Average no. nodes: " << counters.nodeCounter / (reads.size() / 2.0)
         << endl;
    cout << "Total no. Nodes: " << counters.nodeCounter << "\n";

    cout << "Average no. unique matches: "
         << totalUniqueMatches / (reads.size() / 2.0) << endl;
    cout << "Total no. unique matches: " << totalUniqueMatches << "\n";
    cout << "Average no. reported matches "
         << counters.totalReportedPositions / (reads.size() / 2.0) << endl;
    cout << "Total no. reported matches: " << counters.totalReportedPositions
         << "\n";
    cout << "Mapped reads: " << mappedReads << endl;
    cout << "Median number of occurrences per read "
         << findMedian(numberMatchesPerRead, numberMatchesPerRead.size())
         << endl;
    cout << "Reported matches via in-text verification: "
         << counters.cigarsInTextVerification << endl;
    cout << "Unique matches via (partial) in-text verification "
         << counters.usefulCigarsInText << endl;
    cout << "Unique matches via pure in-index matching "
         << counters.cigarsInIndex << endl;
    cout << "In text verification procedures " << counters.inTextStarted
         << endl;
    cout << "Failed in-text verifications procedures: "
         << counters.abortedInTextVerificationCounter << endl;

    cout << "Aborted in-text relative to started "
         << (counters.abortedInTextVerificationCounter * 1.0) /
                counters.inTextStarted
         << endl;
    cout << "Immediate switch after first part: " << counters.immediateSwitch
         << endl;
    cout << "Searches started (does not include immediate switches) : "
         << counters.approximateSearchStarted << endl;

    cout << "Average size of reads: " << sizes / (reads.size() / 2.0) << endl;

    writeToOutput(readsFile + "_output.txt", matchesPerRead, reads);
}

void showUsage() {
    cout << "Usage: ./columba [options] basefilename readfile.[ext]\n\n";
    cout << " [options]\n";
    cout << "  -e  --max-ed\t\tmaximum edit distance [default = 0]\n";
    cout << "  -s  --sa-sparseness\tsuffix array sparseness factor "
            "[default = "
            "1]\n";
    cout << "  -p  --partitioning \tAdd flag to do uniform/static/dynamic "
            "partitioning [default = "
            "dynamic]\n";
    cout << "  -m   --metric\tAdd flag to set distance metric "
            "(editnaive/editopt/hamming) [default = "
            "editopt]\n";
    cout << "  -i  --in-text\tThe tipping point for in-text verification "
            "[default = 5]\n";
    cout << "  -ss --search-scheme\tChoose the search scheme\n  options:\n\t"
         << "kuch1\tKucherov k + 1\n\t"
         << "kuch2\tKucherov k + 2\n\t"
         << "kianfar\t Optimal Kianfar scheme\n\t"
         << "manbest\t Manual best improvement for kianfar scheme (only for ed "
            "= 4)\n\t"
         << "pigeon\t Pigeon hole scheme\n\t"
         << "01*0\t01*0 search scheme\n\t"
         << "custom\tcustom search scheme, the next parameter should be a path "
            "to the folder containing this search scheme\n\n";

    cout << "[ext]\n"
         << "\tone of the following: fq, fastq, FASTA, fasta, fa\n";

    cout << "Following input files are required:\n";
    cout << "\t<base filename>.txt: input text T\n";
    cout << "\t<base filename>.cct: character counts table\n";
    cout << "\t<base filename>.sa.[saSF]: suffix array sample every [saSF] "
            "elements\n";
    cout << "\t<base filename>.bwt: BWT of T\n";
    cout << "\t<base filename>.brt: Prefix occurrence table of T\n";
    cout << "\t<base filename>.rev.brt: Prefix occurrence table of the "
            "reverse "
            "of T\n";
}

int main(int argc, char* argv[]) {

    int requiredArguments = 2; // baseFile of files and file containing reads

    if (argc < requiredArguments) {
        cerr << "Insufficient number of arguments" << endl;
        showUsage();
        return EXIT_FAILURE;
    }
    if (argc == 2 && strcmp("help", argv[1]) == 0) {
        showUsage();
        return EXIT_SUCCESS;
    }

    cout << "Welcome to Columba!\n";

    string saSparse = "1";
    string maxED = "0";
    string searchscheme = "kuch1";
    string customFile = "";
    string inTextPoint = "5";

    PartitionStrategy pStrat = DYNAMIC;
    DistanceMetric metric = EDITOPTIMIZED;

    // process optional arguments
    for (int i = 1; i < argc - requiredArguments; i++) {
        const string& arg = argv[i];

        if (arg == "-p" || arg == "--partitioning") {
            if (i + 1 < argc) {
                string s = argv[++i];
                if (s == "uniform") {
                    pStrat = UNIFORM;
                } else if (s == "dynamic") {
                    pStrat = DYNAMIC;
                } else if (s == "static") {
                    pStrat = STATIC;
                } else {
                    throw runtime_error(
                        s + " is not a partitioning option\nOptions are: "
                            "uniform, static, dynamic");
                }

            } else {
                throw runtime_error(arg + " takes 1 argument as input");
            }
        } else if (arg == "-s" || arg == "--sa-sparseness") {
            if (i + 1 < argc) {
                saSparse = argv[++i];

            } else {
                throw runtime_error(arg + " takes 1 argument as input");
            }
        } else if (arg == "-e" || arg == "--max-ed") {
            if (i + 1 < argc) {
                maxED = argv[++i];
            }
        } else if (arg == "-ss" || arg == "--search-scheme") {
            if (i + 1 < argc) {
                searchscheme = argv[++i];
                if (find(schemes.begin(), schemes.end(), searchscheme) ==
                    schemes.end()) {
                    throw runtime_error(searchscheme +
                                        " is not on option as search scheme");
                }
                if (searchscheme == "custom") {
                    if (i + 1 < argc) {
                        customFile = argv[++i];
                    } else {
                        throw runtime_error(
                            "custom search scheme takes a folder as argument");
                    }
                }

            } else {
                throw runtime_error(arg + " takes 1 argument as input");
            }

        } else if (arg == "-m" || arg == "-metric") {
            if (i + 1 < argc) {
                string s = argv[++i];
                if (s == "editopt") {
                    metric = EDITOPTIMIZED;
                } else if (s == "editnaive") {
                    metric = EDITNAIVE;
                } else if (s == "hamming") {
                    metric = HAMMING;
                } else {
                    throw runtime_error(s +
                                        " is not a metric option\nOptions are: "
                                        "editopt, editnaive, hamming");
                }

            } else {
                throw runtime_error(arg + " takes 1 argument as input");
            }
        } else if (arg == "-i" || arg == "--in-text") {
            if (i + 1 < argc) {
                inTextPoint = argv[++i];
            } else {
                throw runtime_error(arg + " takes 1 argument as input");
            }
        }

        else {
            cerr << "Unknown argument: " << arg << " is not an option" << endl;
            return false;
        }
    }

    length_t ed = stoi(maxED);
    if (ed < 0 || ed > 6) {
        cerr << ed << " is not allowed as maxED should be in [0, 4]" << endl;

        return EXIT_FAILURE;
    }
    length_t saSF = stoi(saSparse);
    if (saSF == 0 || saSF > 256 || (saSF & (saSF - 1)) != 0) {
        cerr << saSF
             << " is not allowed as sparse factor, should be in 2^[0, 8]"
             << endl;
    }

    length_t inTextSwitchPoint = stoi(inTextPoint);

    if (ed != 4 && searchscheme == "manbest") {
        throw runtime_error("manbest only supports 4 allowed errors");
    }

    string baseFile = argv[argc - 2];
    string readsFile = argv[argc - 1];

    cout << "Reading in reads from " << readsFile << endl;
    vector<pair<string, string>> reads;
    try {
        reads = getReads(readsFile);
    } catch (const exception& e) {
        string er = e.what();
        er += " Did you provide a valid reads file?";
        throw runtime_error(er);
    }

    FMIndex bwt = FMIndex(baseFile, inTextSwitchPoint, saSF);

    SearchStrategy* strategy;
    if (searchscheme == "kuch1") {
        strategy = new KucherovKplus1(bwt, pStrat, metric);
    } else if (searchscheme == "kuch2") {
        strategy = new KucherovKplus2(bwt, pStrat, metric);
    } else if (searchscheme == "kianfar") {
        strategy = new OptimalKianfar(bwt, pStrat, metric);
    } else if (searchscheme == "manbest") {
        strategy = new ManBestStrategy(bwt, pStrat, metric);
    } else if (searchscheme == "01*0") {
        strategy = new O1StarSearchStrategy(bwt, pStrat, metric);
    } else if (searchscheme == "pigeon") {
        strategy = new PigeonHoleSearchStrategy(bwt, pStrat, metric);
    } else if (searchscheme == "custom") {
        strategy = new CustomSearchStrategy(bwt, customFile, pStrat, metric);
    } else if (searchscheme == "naive") {
        strategy = new NaiveBackTrackingStrategy(bwt, pStrat, metric);
    } else {
        // should not get here
        throw runtime_error(searchscheme +
                            " is not on option as search scheme");
    }

    doBench(reads, bwt, strategy, readsFile, ed);
    delete strategy;
    cout << "Bye...\n";
}

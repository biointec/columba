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

#include "searchstrategy.h"
#include "wordlength.h"
#include <algorithm>
#include <chrono>
#include <set>
#include <sstream> // used for splitting strings
#include <string.h>

using namespace std;
vector<string> schemes = {"kuch1", "kuch2",  "kianfar", "manbest", "pigeon",
                          "01*0",  "custom", "naive",   "multiple"};

int editDistDP(string P, string O, int maxED) {

    length_t n = P.length();
    length_t m = O.length();
    string* horizontal = &P;
    string* vertical = &O;
    if (n > m) {
        horizontal = &O;
        vertical = &P;
        swap(m, n);
    }

    // check the dimensions of s1 and s2
    if ((max(m, n) - min(m, n)) > (length_t)maxED)
        return numeric_limits<int>::max();

    BitParallelED mat;
    mat.setSequence(horizontal);
    mat.initializeMatrix(maxED);

    // fill in matrix
    for (length_t i = 1; i <= m; i++) {
        bool valid = mat.computeRow(i, vertical->at(i - 1));
        if (!valid) {
            return numeric_limits<int>::max();
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

struct ReadRecord {
    string id;
    string read;

    string qual;

    /**
     * Deep copies the strings
     */
    ReadRecord(string id, string read, string qual)
        : id(id), read(read), qual(qual) {
    }
};

vector<ReadRecord> getReads(const string& file) {
    vector<ReadRecord> reads;
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
        // this is a not readable

        throw runtime_error("extension " + extension +
                            " is not a valid extension for the readsfile");
    } else if (fasta) {
        // fasta file
        string read = "";
        string id = "";
        string qual = ""; // empty quality string for fasta
        string line;

        bool idline = false;

        while (getline(ifile, line)) {
            if (line.empty()) {
                continue; // Skip empty lines
            }

            if (line[0] == '>' || line[0] == '@') {
                // This is an ID line

                if (idline) {
                    // If we already have data, process it and clear
                    reads.emplace_back(id, read, qual);
                    reads.emplace_back(id, Nucleotide::getRevCompl(read), qual);
                    id.clear();
                    read.clear();
                }
                id = line.substr(1); // Extract ID (skip '>')
                idline = true;
            } else {
                // This is a sequence line
                read += line;
            }
        }

        // Process the last entry if it exists
        if (idline) {
            reads.emplace_back(id, read, qual);
            reads.emplace_back(id, Nucleotide::getRevCompl(read), qual);
        }
    } else {
        // fastQ
        string read = "";
        string id = "";
        string qual = "";
        string plusLine = ""; // Skip the '+' line
        string line;

        while (getline(ifile, id) && getline(ifile, read) &&
               getline(ifile, plusLine) && // Skip the '+' line
               getline(ifile, qual)) {
            if (!id.empty() && id[0] != '@') {
                throw runtime_error("File " + file +
                                    "doesn't appear to be in FastQ format");
            }

            if (id.back() == '\n') {
                id.pop_back();
            }
            if (!read.empty() && read.back() == '\n') {
                read.pop_back();
            }

            assert(id.size() > 1);
            id = (id.substr(1));
            reads.emplace_back(id, read, qual);
            reverse(qual.begin(), qual.end());
            reads.emplace_back(id, Nucleotide::getRevCompl(read), qual);
            id.clear(), read.clear(), qual.clear();
        }
    }

    return reads;
}

void writeToOutput(const string& file, const vector<vector<TextOcc>>& mPerRead,
                   const vector<ReadRecord>& reads, length_t textLength,
                   const std::string& basefile) {

    cout << "Writing to output file " << file << " ..." << endl;
    ofstream f2;
    f2.open(file);

    f2 << "@HD"
       << "\t"
       << "VN:1.6"
       << "\t"
       << "SO:queryname"
       << "\n";
    f2 << "@SQ"
       << "\t"
       << "SN:" << basefile << "\t"
       << "LN:" << textLength << "\n";

    for (unsigned int i = 0; i < reads.size(); i += 2) {
        auto id = reads[i].id;

        for (auto m : mPerRead[i]) {
            f2 << m.getSamLine() << "\n";
        }

        for (auto m : mPerRead[i + 1]) {
            f2 << m.getSamLine() << "\n";
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

void doBench(const vector<ReadRecord>& reads, FMIndex& mapper,
             SearchStrategy* strategy, const string& outputfile, length_t ED,
             const std::string& basefile) {

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
    counters.resetCounters();

    auto start = chrono::high_resolution_clock::now();
    for (unsigned int i = 0; i < reads.size(); i += 2) {

        const auto& readrecord = reads[i];
        const auto& revReadrecord = reads[i + 1];

        const string& id = readrecord.id;
        const string& read = readrecord.read;
        const string& revCompl = revReadrecord.read;
        const string& qual = readrecord.qual;
        const string& revQual = revReadrecord.qual;

        if (((i >> 1) - 1) % (8192 / (1 << ED)) == 0) {
            cout << "Progress: " << i / 2 << "/" << reads.size() / 2 << "\r";
            cout.flush();
        }

        sizes += read.size();

        auto matches =
            strategy->matchApprox(read, ED, counters, id, qual, false);
        totalUniqueMatches += matches.size();

        // do the same for the reverse complement
        vector<TextOcc> matchesRevCompl =
            strategy->matchApprox(revCompl, ED, counters, id, revQual, true);
        totalUniqueMatches += matchesRevCompl.size();

        // keep track of the number of mapped reads
        mappedReads += !(matchesRevCompl.empty() && matches.empty());

        matchesPerRead.emplace_back(matches);
        matchesPerRead.emplace_back(matchesRevCompl);

        numberMatchesPerRead.emplace_back(matches.size() +
                                          matchesRevCompl.size());
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

    writeToOutput(outputfile, matchesPerRead, reads, strategy->getText().size(),
                  basefile);
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
    cout << "  -m   --metric\t\tAdd flag to set distance metric "
            "(editnaive/editopt/hamming) [default = "
            "editopt]\n";
    cout << "  -i  --in-text\t\tThe tipping point for in-text verification "
            "[default = 5]\n";
    cout << "  -ks --kmer-size\tThe size of the seeds for dynamic partitioning "
            "[default = 10]\n";
    cout
        << "  -o  --output\t\tThe name of the outputfile. This file must be in "
           ".sam format. [default = ColumbaOutput.sam]\n";
    cout << "  -ss --search-scheme\tChoose the search scheme\n  options:\n\t"
         << "kuch1\t\tKucherov k + 1\n\t"
         << "kuch2\t\tKucherov k + 2\n\t"
         << "kianfar\t\tOptimal Kianfar scheme\n\t"
         << "manbest\t\tManual best improvement for kianfar scheme (only for "
            "ed "
            "= 4)\n\t"
         << "pigeon\t\tPigeon hole scheme\n\t"
         << "01*0\t\t01*0 search scheme\n\t"
         << "custom\t\tcustom search scheme, the next parameter should be a "
            "path "
            "to the folder containing this search scheme\n\t"
         << "multiple\tmultiple search scheme, the next parameter should be a "
            "path "
            "to the folder containing the different search schemes to choose "
            "from with dynamic selection.\n\n";

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

    std::cout << "Using " << LENGTH_TYPE_NAME << std::endl;

    int requiredArguments = 3; // baseFile of files and file containing reads

    if (argc == 2 &&
        (strcmp("help", argv[1]) == 0 || strcmp("--help", argv[1]) == 0 ||
         strcmp("-h", argv[1]) == 0)) {
        showUsage();
        return EXIT_SUCCESS;
    }

    if (argc < requiredArguments) {
        cerr << "Insufficient number of arguments" << endl;
        showUsage();
        return EXIT_FAILURE;
    }

    cout << "Welcome to Columba!\n";

    string saSparse = "1";
    string maxED = "0";
    string searchscheme = "kuch1";
    string customFile = "";
    string inTextPoint = "5";

    string kmerSize = "10";
    string outputfile = "ColumbaOutput.sam";

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
                if (searchscheme == "custom" || searchscheme == "multiple") {
                    if (i + 1 < argc) {
                        customFile = argv[++i];
                    } else {
                        throw runtime_error("custom/multiple search scheme "
                                            "takes a folder as argument");
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
        } else if (arg == "-ks" || arg == "--kmer-size") {
            if (i + 1 < argc) {
                kmerSize = argv[++i];
            } else {
                throw runtime_error(arg + " takes 1 argument as input");
            }
        } else if (arg == "-o" || arg == "--output") {
            if (i + 1 < argc) {
                outputfile = argv[++i];
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
    if (ed < 0 || ed > MAX_K) {
        cerr << ed << " is not allowed as maxED should be in [0, " << MAX_K
             << "]" << endl;

        return EXIT_FAILURE;
    }

    length_t saSF = stoi(saSparse);
    if (saSF == 0 || saSF > 256 || (saSF & (saSF - 1)) != 0) {
        cerr << saSF
             << " is not allowed as sparse factor, should be in 2^[0, 8]"
             << endl;
    }

    length_t inTextSwitchPoint = stoi(inTextPoint);
    length_t kMerSize = stoi(kmerSize);

    if (ed > 6 && inTextSwitchPoint > 0 && metric != HAMMING) {
        cout
            << "Warning in-text verification currently only supports up to k = "
            << 6
            << " for the edit distance. Switching of in-text verification..."
            << endl;
        inTextSwitchPoint = 0;
    }

    if (ed > 4 && searchscheme != "custom" && searchscheme != "multiple" &&
        searchscheme != "naive") {
        throw runtime_error(
            "Hard-coded search schemes are only available for "
            "up to 4 errors. Use a custom search scheme instead.");
    }

    if (ed != 4 && searchscheme == "manbest") {
        throw runtime_error("manbest only supports 4 allowed errors");
    }

    string baseFile = argv[argc - 2];
    string readsFile = argv[argc - 1];

    FMIndex index = FMIndex(baseFile, inTextSwitchPoint, saSF, true, kMerSize);

    cout << "Reading in reads from " << readsFile << endl;
    vector<ReadRecord> reads;
    try {
        reads = getReads(readsFile);
    } catch (const exception& e) {
        string er = e.what();
        er += " Did you provide a valid reads file?";
        throw runtime_error(er);
    }

    SearchStrategy* strategy;
    if (searchscheme == "kuch1") {
        strategy = new KucherovKplus1(index, pStrat, metric);
    } else if (searchscheme == "kuch2") {
        strategy = new KucherovKplus2(index, pStrat, metric);
    } else if (searchscheme == "kianfar") {
        strategy = new OptimalKianfar(index, pStrat, metric);
    } else if (searchscheme == "manbest") {
        strategy = new ManBestStrategy(index, pStrat, metric);
    } else if (searchscheme == "01*0") {
        strategy = new O1StarSearchStrategy(index, pStrat, metric);
    } else if (searchscheme == "pigeon") {
        strategy = new PigeonHoleSearchStrategy(index, pStrat, metric);
    } else if (searchscheme == "custom") {
        strategy = new CustomSearchStrategy(index, customFile, pStrat, metric);
    } else if (searchscheme == "multiple") {
        strategy =
            new MultipleSchemesStrategy(index, customFile, pStrat, metric);
    } else if (searchscheme == "naive") {
        strategy = new NaiveBackTrackingStrategy(index, pStrat, metric);
    } else {
        // should not get here
        throw runtime_error(searchscheme +
                            " is not on option as search scheme");
    }

    doBench(reads, index, strategy, outputfile, ed, baseFile);
    delete strategy;
    cout << "Bye...\n";
}

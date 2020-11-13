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
#include "searchstrategy.h"
#include <algorithm>
#include <chrono>
#include <set>
#include <string.h>

using namespace std;
vector<string> schemes = {"kuch1",   "kuch2",  "kianfar",
                          "manbest", "pigeon", "01*0"};
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

void printMatches(vector<AppMatch> matches, string text, bool printLine,
                  string duration, FMIndex* mapper, string name) {

    cout << endl;

    tuple<length_t, length_t, length_t> counters = mapper->getCounters();

    cout << name << ":\tduration: " << duration
         << "µs\t nodes visited: " << get<0>(counters)
         << "\t matrix elements written: " << get<1>(counters)
         << "\t startpositions reported: " << get<2>(counters)
         << " #matches: " << matches.size() << endl;

    for (auto match : matches) {
        cout << "Found match at position " << match.range.begin << " with ED "
             << match.editDist << endl;

        cout << "\tCorresponding substring:\t"
             << text.substr(match.range.begin,
                            match.range.end - match.range.begin)
             << endl;
    }
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
    reads.reserve(100000);

    const auto& extension = getFileExt(file);

    bool fasta =
        (extension == "FASTA") || (extension == "fasta") || (extension == "fa");
    bool fastq = (extension == "fq") || (extension == "fastq");

    ifstream ifile = getStream(file);
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
                    read.clear();
                }

                p = (line.substr(1));

            } else {
                read += line;
            }
        }
        if (!read.empty()) {

            reads.push_back(make_pair(p, read));
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
            read.clear();
        }
    }

    return reads;
}

double
avgVec(vector<length_t> const& v) // note: the average must not be an integer
{
    return v.empty() ? 0.0 : accumulate(v.begin(), v.end(), 0.0) / v.size();
    ;
}

length_t sum(vector<length_t> const& v) {
    return accumulate(v.begin(), v.end(), 0.0);
}
void writeToOutput(const string& file, const vector<vector<AppMatch>>& mPerRead,
                   const vector<pair<string, string>>& reads) {

    cout << "Writing to output file " << file << " ..." << endl;
    ofstream f2;
    f2.open(file);

    f2 << "identifier\tposition\tlength\tED\n";
    for (unsigned int i = 0; i < reads.size(); i++) {
        auto id = reads[i].first;

        for (auto m : mPerRead[i]) {
            f2 << id << "\t" << m.range.begin << "\t" << m.range.width() << "\t"
               << m.editDist << "\n";
        }
    }
    f2.close();
}

void doBench(vector<pair<string, string>>& reads, BidirecFMIndex* mapper,
             SearchStrategy* strategy, string readsFile, length_t ED) {
    string text = mapper->getText();
    vector<length_t> durations;
    vector<length_t> nodes;
    vector<length_t> matrixElements;
    vector<length_t> uniqueMatches;
    vector<length_t> totalreportedmatches;

    vector<AppMatch> matches;

    cout << "Benchmarking with " << strategy->getName()
         << " strategy for max distance " << ED << " with "
         << strategy->getPartitioningStrategy() << " partitioning and using "
         << strategy->getEditOrHamming() << " distance " << endl;
    cout.precision(2);

    vector<vector<AppMatch>> matchesPerRead = {};
    matchesPerRead.reserve(reads.size());

    auto start = chrono::high_resolution_clock::now();
    for (unsigned int i = 0; i < reads.size(); i++) {

        const auto& p = reads[i];

        auto originalPos = p.first;
        string read = p.second;

        if ((i - 1) % (8192 / (1 << ED)) == 0) {
            cout << "Progress: " << i << "/" << reads.size() << "\r";
            cout.flush();
        }

        matches = strategy->matchApprox(read, ED);
        auto counters = mapper->getCounters();

        nodes.push_back(get<0>(counters));
        matrixElements.push_back(get<1>(counters));
        totalreportedmatches.push_back(get<2>(counters));

        uniqueMatches.push_back(matches.size());
        matchesPerRead.push_back(matches);

        // correctness check, comment this out if you want to check
        // For each reported match the reported edit distance is checked and
        // compared to a recalculated value using a single banded matrix this is
        // slow
        // WARNING: this checks the EDIT DISTANCE, for it might be that the
        // hamming distance is higher
        /*for (auto match : matches) {

            string O = text.substr(match.range.begin,
                                   match.range.end - match.range.begin);

            int trueED = editDistDP(read, O, ED);
            int foundED = match.editDist;
            if (foundED != trueED) {
                cout << i << "\n";
                cout << "Wrong ED!!"
                     << "\n";
                cout << "P: " << read << "\n";
                cout << "O: " << O << "\n";
                cout << "true ED " << trueED << ", found ED " << foundED << "\n"
                     << match.range.begin << "\n";
            }
        }*/

        // this block checks if at least one occurrence is found and if the
        // identifier is a number and then checks if this position is found as a
        // match (for checking correctness) if you want to check if the position
        // is found as a match make sure that the identifier of the read is the
        // position. Out comment this block for the check
        /* bool originalFound = true;
         try {
              length_t pos = stoull(originalPos);
              originalFound = false;

              for (auto match : matches) {

                  if (match.range.begin >= pos - (ED + 2) &&
                      match.range.begin <= pos + (ED + 2)) {
                      originalFound = true;
                      break;
                  }
              }
          } catch (const std::exception& e) {
              // nothing to do, identifier is  not the orignal position
          }

         // check if at least one occurrence was found (for reads that were
         // sampled from actual reference) Out-cooment this block if you want to
         // do this.
         if (matches.size() == 0 || (!originalFound)) {
              cout << "Could not find occurrence for " << originalPos << endl;
          }*/
    }
    auto finish = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = finish - start;
    cout << "Progress: " << reads.size() << "/" << reads.size() << "\n";
    cout << "Results for " << strategy->getName() << endl;

    cout << "Total duration: " << fixed << elapsed.count() << "s\n";

    cout << "Average no. nodes: " << avgVec(nodes) << endl;
    cout << "Total no. Nodes: " << sum(nodes) << "\n";
    cout << "Average no. matrix elements written: " << avgVec(matrixElements)
         << endl;
    cout << "Total no. matrix elements: " << sum(matrixElements) << "\n";
    cout << "Average no. unique matches: " << avgVec(uniqueMatches) << endl;
    cout << "Total no. unique matches: " << sum(uniqueMatches) << "\n";
    cout << "Average no. reported matches " << avgVec(totalreportedmatches)
         << endl;
    cout << "Total no. reported matches: " << sum(totalreportedmatches) << "\n";

    size_t lastindex = readsFile.find_last_of(".");
    string rawname = readsFile.substr(0, lastindex);

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
    cout << "  -h   --hamming\tAdd flag to use hamming distance [default = "
            "false]\n";
    cout << "  -ss --search-scheme\tChoose the search scheme\n  options:\n\t"
         << "kuch1\tKucherov k + 1\n\t"
         << "kuch2\tKucherov k + 2\n\t"
         << "kianfar\t Optimal Kianfar scheme\n\t"
         << "manbest\t Manual best improvement for kianfar scheme (only for ed "
            "= 4)\n\t"
         << "pigeon\t Pigeon hole scheme\n\t"
         << "01*0\t01*0 search scheme\n\n";

    cout << "[ext]\n"
         << "\tone of the following: fq, fastq, FASTA, fasta, fa";

    cout << "Following input files are required:\n";
    cout << "\t<base filename>.txt: input text T\n";
    cout << "\t<base filename>.cct: charachter counts table\n";
    cout << "\t<base filename>.sa.[saSF]: suffix array sample every [saSF] "
            "elements\n";
    cout << "\t<base filename>.bwt: BWT of T\n";
    cout << "\t<base filename>.brt: Prefix occurrence table of T\n";
    cout << "\t<base filename>.rev.brt: Prefix occurrence table of the "
            "reverse "
            "of T\n";
}

int main(int argc, char* argv[]) {

    int requiredArguments = 2; // prefix of files and file containing reads

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

    PartitionStrategy pStrat = DYNAMIC;
    bool edit = true;

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

            } else {
                throw runtime_error(arg + " takes 1 argument as input");
            }

        } else if (arg == "-h" || arg == "-hamming") {
            edit = false;
        }

        else {
            cerr << "Unknown argument: " << arg << " is not an option" << endl;
            return false;
        }
    }

    length_t ed = stoi(maxED);
    if (ed < 0 || ed > 4) {
        cerr << ed << " is not allowed as maxED should be in [0, 4]" << endl;

        return EXIT_FAILURE;
    }
    length_t saSF = stoi(saSparse);
    if (saSF == 0 || saSF > 256 || (saSF & (saSF - 1)) != 0) {
        cerr << saSF
             << " is not allowed as sparse factor, should be in 2^[0, 8]"
             << endl;
    }

    if (ed != 4 && searchscheme == "manbest") {
        throw runtime_error("manbest only supports 4 allowed errors");
    }

    string prefix = argv[argc - 2];
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
    cout << "Start creation of BWT approximate matcher" << endl;

    BidirecFMIndex bwt = BidirecFMIndex(prefix, saSF);

    SearchStrategy* strategy;
    if (searchscheme == "kuch1") {
        strategy = new KucherovKplus1(&bwt, pStrat, edit);
    } else if (searchscheme == "kuch2") {
        strategy = new KucherovKplus2(&bwt, pStrat, edit);
    } else if (searchscheme == "kianfar") {
        strategy = new OptimalKianfar(&bwt, pStrat, edit);
    } else if (searchscheme == "manbest") {
        strategy = new ManBestStrategy(&bwt, pStrat, edit);
    } else if (searchscheme == "01*0") {
        strategy = new O1StarSearchStrategy(&bwt, pStrat, edit);
    } else if (searchscheme == "pigeon") {
        strategy = new PigeonHoleSearchStrategy(&bwt, pStrat, edit);
    } else {
        // should not get here
        throw runtime_error(searchscheme +
                            " is not on option as search scheme");
    }
    doBench(reads, &bwt, strategy, readsFile, ed);
    delete strategy;
    cout << "Bye...\n";
}

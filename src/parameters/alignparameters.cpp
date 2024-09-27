#include "alignparameters.h"
#include "../indexinterface.h"
#include "../logger.h"
#include "../searchstrategy.h"
#include <algorithm>
#include <thread> // for thread

#include <vector>

const std::vector<std::string> schemes = {
    "kuch1",  "kuch2", "kianfar",  "pigeon", "01*0",
    "custom", "naive", "multiple", "minU",   "columba"};

/**
 * Option for the command line arguments considering the partitioning strategy.
 */
class PartitioningOption : public ParameterOption {
  public:
    PartitioningOption()
        : ParameterOption("p", "partitioning", true, STRING, ADVANCED) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        if (arg == "uniform") {
            params.pStrategy = UNIFORM;
        } else if (arg == "dynamic") {
            params.pStrategy = DYNAMIC;
        } else if (arg == "static") {
            params.pStrategy = STATIC;
        } else {
            logger.logWarning(arg +
                              " is not a partitioning option\nOptions "
                              "are: uniform, static, dynamic" +
                              ignoreMessage());
        }
    }

    std::string getDescription() const override {
        return "Partitioning strategy to use. Options are: uniform, "
               "static, "
               "dynamic. Default is dynamic.";
    }
};

#ifndef RUN_LENGTH_COMPRESSION
/**
 * Option for the command line arguments considering the sparseness factor of
 * the suffix array.
 */
class SparsenessOption : public ParameterOption {
  public:
    SparsenessOption()
        : ParameterOption("s", "sa-sparseness", true, INTEGER, ADVANCED) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        length_t sFactor = params.sparsenessFactor;
        try {
            sFactor = std::stoi(arg);
        } catch (...) {
            logger.logWarning("Sparseness factor should be an integer." +
                              ignoreMessage());
        }
        if (sFactor == 0 || ((sFactor & (sFactor - 1)) != 0)) {
            logger.logWarning(std::to_string(sFactor) +
                              " is not allowed as sparse factor, should be "
                              "in a power of 2." +
                              ignoreMessage());
            sFactor = params.sparsenessFactor;
        }
        params.sparsenessFactor = sFactor;
    }

    std::string getDescription() const override {
        return "Sparseness factor to use. Should be an integer in 2^[0, "
               "8]. Default is 4. The suffix array with this sparseness factor "
               "must have been constructed.";
    }
};

/**
 * Option for the command line arguments considering the in-text verification
 * switch point.
 */
class InTextOption : public ParameterOption {
  public:
    InTextOption() : ParameterOption("i", "in-text", true, INTEGER, ADVANCED) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        length_t iFactor = params.inTextVerificationPoint;
        try {
            iFactor = std::stoi(arg);
        } catch (...) {
            logger.logWarning("In-text verification switch point should be a "
                              "positive integer." +
                              ignoreMessage());
        }
        if (iFactor < 0) {
            logger.logWarning(std::to_string(iFactor) +
                              " is not allowed as the in-text verification "
                              "switch point, should be in [0, inf]." +
                              ignoreMessage());
            iFactor = params.inTextVerificationPoint;
        }
        params.inTextVerificationPoint = iFactor;
    }

    std::string getDescription() const override {
        return "In-text verification switch point. Should be a positive "
               "integer. Default is 5.";
    }
};

/**
 * Option for the command line arguments to turn of CIGAR string generation.
 */
class NoCigarOption : public ParameterOption {
  public:
    NoCigarOption() : ParameterOption("nC", "no-CIGAR", false, NONE, OUTPUT) {
    }
    void process(const std::string& arg, Parameters& params) const override {
        params.noCIGAR = true;
    }

    std::string getDescription() const override {
        return "Do not output CIGAR strings for SAM format.";
    }
};
#endif

/**
 * Option for the command line arguments considering the search scheme to be
 * used.
 */
class SearchSchemeOption : public ParameterOption {
  public:
    SearchSchemeOption()
        : ParameterOption("S", "search-scheme", true, STRING, ADVANCED) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        params.searchScheme = arg;
        if (std::find(schemes.begin(), schemes.end(), params.searchScheme) ==
            schemes.end()) {
            logger.logWarning(params.searchScheme +
                              " is not an option as a search scheme." +
                              ignoreMessage());
            params.searchScheme = "columba";
        }
    }

    std::string getDescription() const override {
        std::string help = "Search scheme to use. Options are: ";
        help += schemes.front();
        for (size_t i = 1; i < schemes.size(); ++i) {
            help += ", " + schemes[i];
        }
        help += ". Default is columba.";
        return help;
    }
};

/**
 * Option for the command line arguments considering the custom search scheme to
 * be used.
 */
class CustomOption : public ParameterOption {
  public:
    CustomOption() : ParameterOption("c", "custom", true, STRING, ADVANCED) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        params.custom = arg;
    }

    std::string getDescription() const override {
        return "Path to custom search scheme (overrides default "
               "search scheme).";
    }
};

/**
 * Option for the command line arguments considering the dynamic selection
 * custom collection of search schemes option to be used.
 */
class DynamicSelectionOption : public ParameterOption {
  public:
    DynamicSelectionOption()
        : ParameterOption("d", "dynamic-selection", true, STRING, ADVANCED) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        params.dynamicSelectionPath = arg;
    }

    std::string getDescription() const override {
        return "Path to custom search scheme with dynamic selection "
               "(overrides "
               "default search scheme).";
    }
};

/**
 * Option for the command line arguments considering the distance metric to be
 * used.
 */
class MetricOption : public ParameterOption {
  public:
    MetricOption() : ParameterOption("m", "metric", true, STRING, ALIGNMENT) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        if (arg == "edit") {
            params.metric = EDIT;
        } else if (arg == "hamming") {
            params.metric = HAMMING;
        } else {
            logger.logWarning(arg +
                              " is not a distance metric option\nOptions are: "
                              "edit, hamming" +
                              ignoreMessage());
        }
    }

    std::string getDescription() const override {
        return "Distance metric to use. Options are: edit, hamming. Default is "
               "edit.";
    }
};

/**
 * Option for the command line arguments considering the log file to be used.
 */
class LogFileOption : public ParameterOption {
  public:
    LogFileOption() : ParameterOption("l", "log-file", true, STRING, OUTPUT) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        params.logFile = arg;
    }

    std::string getDescription() const override {
        return "Path to the log file. Default is stdout.";
    }
};

/**
 * Option for the command line arguments considering the first reads file to be
 * used.
 */
class FirstReadsOption : public ParameterOption {
  public:
    FirstReadsOption()
        : ParameterOption("f", "first-reads-file", true, STRING, REQUIRED) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        params.firstReadsFile = arg;
    }

    std::string getDescription() const override {
        return "Path to the (first) reads file.";
    }
};

/**
 * Option for the command line arguments considering the second reads file to be
 * used.
 */
class SecondReadsOption : public ParameterOption {
  public:
    SecondReadsOption()
        : ParameterOption("F", "second-reads-file", true, STRING,
                          PAIRED_END_OPTION) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        params.secondReadsFile = arg;
        params.sMode = PAIRED_END;
    }

    std::string getDescription() const override {
        return "Path to the second reads file (optional). If this is set "
               "Columba will use paired-end alignment.";
    }
};

/**
 * Option for the command line arguments considering the reference to be
 * used.
 */
class ReferenceOption : public ParameterOption {
  public:
    ReferenceOption()
        : ParameterOption("r", "reference", true, STRING, REQUIRED) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        params.base = arg;
    }

    std::string getDescription() const override {
        return "Path to the basename of the index.";
    }
};

/**
 * Option for the command line arguments considering the output file to be
 * used.
 */
class OutputFileOption : public ParameterOption {
  public:
    OutputFileOption()
        : ParameterOption("o", "output-file", true, STRING, OUTPUT) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        params.outputFile = arg;
        if (params.outputFile.size() > 3 &&
            params.outputFile.substr(params.outputFile.size() - 4) != ".sam") {
            params.outputIsSAM = false;
        }
    }

    std::string getDescription() const override {
        return "Path to the output file. Should be .sam or .rhs. Default is "
               "ColumbaOutput.sam.";
    }
};

/**
 * Option for the command line arguments considering the maximal distance to be
 * used.
 */
class MaxDistanceOption : public ParameterOption {
  public:
    MaxDistanceOption()
        : ParameterOption("e", "max-distance", true, STRING, ALIGNMENT) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        try {
            params.maxDistance = std::stoi(arg);
        } catch (...) {
            logger.logWarning("Max Distance should be a positive integer" +
                              ignoreMessage());
        }
        if (params.maxDistance < 0 || params.maxDistance > MAX_K) {
            logger.logWarning("Max Distance should be in [0, " +
                              std::to_string(MAX_K) + "]" + ignoreMessage());
        }
    }

    std::string getDescription() const override {
        return "The maximum allowed distance (for ALL mode). Default is 0.";
    }
};

/**
 * Option for the command line arguments considering the k-mer size to be
 * used.
 */
class KmerSizeOption : public ParameterOption {
  public:
    KmerSizeOption()
        : ParameterOption("K", "kmer-size", true, INTEGER, ADVANCED) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        length_t defaultKmerSize = params.kmerSize;
        try {
            params.kmerSize = std::stoi(arg);
        } catch (...) {
            logger.logWarning("kmerSkipSize should be an integer." +
                              ignoreMessage());
        }
        if (params.kmerSize < 0 || params.kmerSize > 15) {
            logger.logWarning("kmerSkipSize should be in [0, 15]." +
                              ignoreMessage());
            params.kmerSize = defaultKmerSize;
        }
    }

    std::string getDescription() const override {
        return "The size of k-mers in the hash table (used as seeds during "
               "partitioning). Default is 10.";
    }
};

/**
 * Option for the command line arguments considering the alignment mode to be
 * used.
 */
class AlignModeOption : public ParameterOption {
  public:
    AlignModeOption()
        : ParameterOption("a", "align-mode", true, STRING, ALIGNMENT) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        if (arg == "all") {
            params.mMode = ALL;
        } else if (arg == "best") {
            params.mMode = BEST;
        } else {
            logger.logWarning(
                arg + " is not a valid mode option\nOptions are: all, best" +
                ignoreMessage());
        }
    }

    std::string getDescription() const override {
        return "Alignment mode to use. Options are: all, best. Default is "
               "best.";
    }
};

/**
 * Option for the command line arguments considering the minimal identity to be
 * used.
 */
class MinIdentityOption : public ParameterOption {
  public:
    MinIdentityOption()
        : ParameterOption("I", "minIdentity", true, INTEGER, ALIGNMENT) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        length_t defaultMinIdentity = params.minIdentity;
        try {
            params.minIdentity = std::stoi(arg);
        } catch (...) {
            logger.logWarning("minIdentity should be an integer." +
                              ignoreMessage());
        }
        if (params.minIdentity < 50 || params.minIdentity > 100) {
            logger.logWarning("minIdentity should be in [50, 100]." +
                              ignoreMessage());
            params.minIdentity = defaultMinIdentity;
        }
    }

    std::string getDescription() const override {
        return "The minimum identity for alignments in BEST mode. Default is "
               "95.";
    }
};

/**
 * Option for the command line arguments considering the number of threads  to
 * be used.
 */
class ThreadsOption : public ParameterOption {
  public:
    ThreadsOption()
        : ParameterOption("t", "threads", true, INTEGER, PARALLELIZATION) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        try {
            params.nThreads = std::stoi(arg);
            if (params.nThreads <= 0) {
                logger.logWarning("Only positive values are allowed for the "
                                  "number of threads. Using 1 thread instead.");
                params.nThreads = 1;
            }

        } catch (...) {
            logger.logWarning(
                "Number of threads should be an integer. Set to 1 by default.");
            params.nThreads = 1;
        }
    }

    std::string getDescription() const override {
        return "The number of threads to be used. Default is 1.";
    }
};

/**
 * Option for the command line arguments considering the orientation for
 * paired-end reads file to be used.
 */
class OrientationOption : public ParameterOption {
  public:
    OrientationOption()
        : ParameterOption("O", "orientation", true, STRING, PAIRED_END_OPTION) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        if (arg == "fr") {
            params.orientation = FR;
        } else if (arg == "rf") {
            params.orientation = RF;
        } else if (arg == "ff") {
            params.orientation = FF;
        } else {
            logger.logWarning(arg +
                              " is not a valid orientation option\nOptions "
                              "are: fr, rf, ff, rr" +
                              ignoreMessage());
        }
    }

    std::string getDescription() const override {
        return "Orientation of the paired end reads. Options are: fr, rf, ff. "
               "Default is fr.";
    }
};

/**
 * Option for the command line arguments considering the max insert size to be
 * used.
 */
class MaxInsertSizeOption : public ParameterOption {
  public:
    MaxInsertSizeOption()
        : ParameterOption("X", "max-insert-size", true, INTEGER,
                          PAIRED_END_OPTION) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        length_t defaultMaxInsertSize = params.maxInsertSize;
        try {
            params.maxInsertSize = std::stoi(arg);
        } catch (...) {
            logger.logWarning("maxInsertSize should be an integer." +
                              ignoreMessage());
        }
        if (params.maxInsertSize < 0) {
            logger.logWarning("maxInsertSize should be a positive integer." +
                              ignoreMessage());
            params.maxInsertSize = defaultMaxInsertSize;
        }
    }

    std::string getDescription() const override {
        return "The maximum insert size for paired end reads. Default is 500.";
    }
};

/**
 * Option for the command line arguments considering the minimal insert size to
 * be used.
 */
class MinInsertSizeOption : public ParameterOption {
  public:
    MinInsertSizeOption()
        : ParameterOption("N", "min-insert-size", true, INTEGER,
                          PAIRED_END_OPTION) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        length_t defaultMinInsertSize = params.minInsertSize;
        try {
            params.minInsertSize = std::stoi(arg);
        } catch (...) {
            logger.logWarning("minInsertSize should be an integer." +
                              ignoreMessage());
        }
        if (params.minInsertSize < 0) {
            logger.logWarning("minInsertSize should be a positive integer." +
                              ignoreMessage());
            params.minInsertSize = defaultMinInsertSize;
        }
    }

    std::string getDescription() const override {
        return "The minimum insert size for paired end reads. Default is 0.";
    }
};

/**
 * Option for the command line arguments to turn of inference of paired-end
 * parameters.
 */
class NotInferPairedEndParametersOption : public ParameterOption {
  public:
    NotInferPairedEndParametersOption()
        : ParameterOption("nI", "no-inferring", false, NONE,
                          PAIRED_END_OPTION) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        params.inferPairedEndParameters = false;
    }

    std::string getDescription() const override {
        return "Do not infer paired-end parameters. Do not infer the paired "
               "end parameters. By default the parameters are inferred. If "
               "this option is set the values provided by, -O (default FR), -X "
               "(default 500) and "
               "-N (default 0) are used.";
    }
};

/**
 * Option for the command line arguments to allow discordant pairings
 */
class DiscordantOption : public ParameterOption {
  public:
    DiscordantOption()
        : ParameterOption("D", "discordant", false, INTEGER,
                          PAIRED_END_OPTION) {
    }
    void process(const std::string& arg, Parameters& params) const override {
        params.discordantAllowed = true;
        if (!arg.empty()) {
            try {
                params.nDiscordant = std::stoi(arg);
            } catch (...) {
                logger.logWarning(
                    "Number of allowed discordant pairs should be an integer." +
                    ignoreMessage());
                params.discordantAllowed = false;
            }
        }
    }
    std::string getDescription() const override {
        return "Allow discordant alignments. Optionally you can provide "
               "the maximal number of discordant alignments per pair to "
               "allow. Default is 100000.";
    }

    bool argumentRequired() const override {
        return false;
    }
};

/**
 * Option for the command line arguments to turn of output of unmapped reads.
 */
class NoUnmappedOption : public ParameterOption {
  public:
    NoUnmappedOption()
        : ParameterOption("nU", "no-unmapped", false, NONE, OUTPUT) {
    }
    void process(const std::string& arg, Parameters& params) const override {
        params.noUnmappedRecord = true;
    }

    std::string getDescription() const override {
        return "Do not output unmapped reads.";
    }
};

/**
 * Option for the command line arguments to turn on XA tags for multiple
 * alignments.
 */
class XATagOption : public ParameterOption {
  public:
    XATagOption() : ParameterOption("XA", "XA-tag", false, NONE, OUTPUT) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        params.XATag = true;
    }

    std::string getDescription() const override {
        return "Output secondary alignments in XA tag for SAM format.";
    }
};

/**
 * Option for the command line arguments to trigger the help.
 */
class HelpOption : public ParameterOption {
  public:
    HelpOption() : ParameterOption("h", "help", false, NONE, HELP) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        ParametersInterface::printHelp();
        exit(0);
    }

    std::string getDescription() const override {
        return "Print this help message.";
    }
};

/**
 * Option for the command line arguments to ensure that the output is in the
 * same order as the input.
 */
class ReorderOption : public ParameterOption {
  public:
    ReorderOption() : ParameterOption("R", "reorder", false, NONE, OUTPUT) {
    }

    void process(const std::string& arg, Parameters& params) const override {
        params.reorder = true;
    }

    std::string getDescription() const override {
        return "Guarantees that output SAM or RHS records are printed in the "
               "order corresponding to the order of reads in the original "
               "file. Setting this will cause Columba to be somewhat slower "
               "and use somewhat more memory.";
    }
};

const std::vector<std::shared_ptr<Option>> ParametersInterface::options = {
    std::make_shared<PartitioningOption>(),
    std::make_shared<MetricOption>(),
    std::make_shared<LogFileOption>(),
    std::make_shared<FirstReadsOption>(),
    std::make_shared<SecondReadsOption>(),
    std::make_shared<ReferenceOption>(),
    std::make_shared<OutputFileOption>(),
    std::make_shared<KmerSizeOption>(),
    std::make_shared<AlignModeOption>(),
    std::make_shared<MinIdentityOption>(),
    std::make_shared<MaxDistanceOption>(),
    std::make_shared<ThreadsOption>(),
    std::make_shared<NotInferPairedEndParametersOption>(),
    std::make_shared<OrientationOption>(),
    std::make_shared<MaxInsertSizeOption>(),
    std::make_shared<MinInsertSizeOption>(),
    std::make_shared<DiscordantOption>(),
    std::make_shared<NoUnmappedOption>(),

    std::make_shared<XATagOption>(),
    std::make_shared<SearchSchemeOption>(),
    std::make_shared<CustomOption>(),
    std::make_shared<DynamicSelectionOption>(),
    std::make_shared<ReorderOption>(),
    std::make_shared<HelpOption>(),
#ifndef RUN_LENGTH_COMPRESSION // options that are not available with run-length
                               // compression
    std::make_shared<InTextOption>(),
    std::make_shared<SparsenessOption>(),
    std::make_shared<NoCigarOption>()
#endif
};

void ParametersInterface::printHeader() {
    std::cout << "Usage: columba [OPTIONS]\n\n";
}

Parameters Parameters::processOptionalArguments(int argc, char** argv) {

    Parameters params;
    params.command = getCommand(argc, argv);

    params.parse(argc, argv);

    if (params.maxDistance > 4 && params.custom.empty() &&
        params.dynamicSelectionPath.empty() && params.searchScheme != "naive") {
        if (params.searchScheme != "minU" && params.searchScheme != "columba") {

            throw std::runtime_error(
                "Hard-coded search schemes are only available for "
                "up to 4 errors. Use a custom search scheme instead.");
        } else if (params.searchScheme == "minU" && params.maxDistance > 7) {
            throw std::runtime_error(
                "MinU search scheme is only available for up to 7 errors.");
        }
    }

    if (params.metric == EDIT && params.maxDistance > MAX_K_EDIT) {
        throw std::runtime_error("Edit distance only supports up to " +
                                 std::to_string(MAX_K_EDIT) + " errors");
    }

    if (params.maxDistance != 4 && params.searchScheme == "manbest") {
        throw std::runtime_error("manbest only supports 4 allowed errors");
    }

    if (params.maxDistance > MAX_K) {
        throw std::runtime_error("Max distance cannot be higher than " +
                                 std::to_string(MAX_K));
    }

    if (!params.inferPairedEndParameters &&
        params.maxInsertSize < params.minInsertSize) {
        throw std::runtime_error(
            "Max insert size cannot be smaller than min insert size");
    }

    // check if threads does not exceed hardware options
    if (params.nThreads > std::thread::hardware_concurrency()) {
        std::stringstream ss;
        ss << "The entered number of threads: " << params.nThreads
           << " is higher than the number of threads "
              "available: "
           << std::thread::hardware_concurrency()
           << ". Setting number of threads to maximum";
        logger.logWarning(ss);
        params.nThreads = std::thread::hardware_concurrency();
    }

    if (params.secondReadsFile.empty()) {
        if (params.discordantAllowed) {
            logger.logWarning("Discordant pairs does not apply to single-end "
                              "alignment. " +
                              DiscordantOption().ignoreMessage());
            params.discordantAllowed = false;
        }
        if (!params.inferPairedEndParameters) {
            logger.logWarning(
                "Paired-end parameters are never inferred for "
                "single-end alignment. " +
                NotInferPairedEndParametersOption().ignoreMessage());
        }
    } else {
        // paired-end reads
        if (params.XATag) {
            logger.logWarning("XA tag is not supported for paired-end reads. " +
                              XATagOption().ignoreMessage());
            params.XATag = false;
        }
        if (!params.outputIsSAM) {
            params.outputFile =
                params.outputFile.substr(0, params.outputFile.size() - 4) +
                ".sam";
            logger.logWarning("RHS output is not supported for paired-end "
                              "reads. Renaming output file to " +
                              params.outputFile);
            params.outputIsSAM = true;
        }

        if (!params.outputIsSAM) {
            if (params.XATag) {
                logger.logWarning("XA tag is not supported for RHS output. " +
                                  XATagOption().ignoreMessage());
                params.XATag = false;
            }
#ifndef RUN_LENGTH_COMPRESSION
            if (params.noCIGAR) {
                logger.logWarning(
                    "CIGAR string is not supported for RHS output. " +
                    NoCigarOption().ignoreMessage());
                params.noCIGAR = false;
            }
#endif
        }
    }

    return params;
}

std::unique_ptr<SearchStrategy>
Parameters::createStrategy(IndexInterface& index) const {
    std::unique_ptr<SearchStrategy> strategy;

    if (!dynamicSelectionPath.empty()) {
        strategy.reset(new MultipleSchemesStrategy(
            index, dynamicSelectionPath, pStrategy, metric, mMode, sMode));
    } else if (!custom.empty()) {
        DynamicCustomStrategy* customStrategy =
            DynamicCustomStrategy::createDynamicCustomSearchStrategy(
                index, custom, pStrategy, metric, mMode, sMode);

        strategy.reset(customStrategy);

        if (mMode == ALL) {
            if (!strategy->supportsDistanceScore(maxDistance)) {
                throw std::runtime_error(
                    "Custom search scheme with path " + custom +
                    " does not support the given distance score: " +
                    std::to_string(maxDistance));
            }

        } else if (mMode == BEST) {
            int maxForBest = 0;
            if (!strategy->supportsBestMapping(maxForBest)) {
                throw std::runtime_error("Custom search scheme with path " +
                                         custom +
                                         " does not support best mapping");
            }
            if (maxForBest < MAX_K) {
                logger.logWarning(
                    "Custom search scheme with path " + custom +
                    " does not support best mapping for more than " +
                    std::to_string(maxForBest) + " errors");
            }
        }

    } else {
        // Handle other search schemes
        if (searchScheme == "kuch1") {
            strategy.reset(
                new KucherovKPlus1(index, pStrategy, metric, mMode, sMode));
        } else if (searchScheme == "kuch2") {
            strategy.reset(
                new KucherovKPlus2(index, pStrategy, metric, mMode, sMode));
        } else if (searchScheme == "kianfar") {
            strategy.reset(
                new OptimalKianfar(index, pStrategy, metric, mMode, sMode));
        } else if (searchScheme == "01*0") {
            strategy.reset(new O1StarSearchStrategy(index, pStrategy, metric,
                                                    mMode, sMode));
        } else if (searchScheme == "pigeon") {
            strategy.reset(new PigeonHoleSearchStrategy(index, pStrategy,
                                                        metric, mMode, sMode));
        } else if (searchScheme == "naive") {
            strategy.reset(new NaiveBackTrackingStrategy(index, pStrategy,
                                                         metric, mMode, sMode));
        } else if (searchScheme == "minU") {
            strategy.reset(
                new MinUSearchStrategy(index, pStrategy, metric, mMode, sMode));
        } else if (searchScheme == "columba") {
            DynamicColumbaStrategy* columbaStrategy =
                DynamicColumbaStrategy::createDynamicColumbaStrategy(
                    index, pStrategy, metric, mMode, sMode);

            // Assuming DynamicCustomStrategy is a subclass of SearchStrategy
            strategy.reset(columbaStrategy);
        } else {
            // should not get here
            throw std::runtime_error(searchScheme +
                                     " is not an option as search scheme");
        }
    }

    // Set common parameters for all strategies
    strategy->setDiscordantAllowed(discordantAllowed);
    strategy->setNumDiscordantAllowed(nDiscordant);
    strategy->setUnmappedSam(!noUnmappedRecord);
    strategy->setXATag(XATag);
    strategy->setSamOutput(outputIsSAM);

    return strategy;
}

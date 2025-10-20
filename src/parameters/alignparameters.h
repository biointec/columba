#ifndef ALIGN_PARAMETERS_H
#define ALIGN_PARAMETERS_H

#include "../definitions.h"
#include "parameters.h"

// forward declaration
class SearchStrategy;
class IndexInterface;

/**
 * Struct Parameters contains all the parameters that can be set by the user
 * Provides a method to process the optional arguments.
 * Provides a method to create a search strategy based on the parameters.
 * Provides a method to get the maximum distance or identity based on the mode.
 */
struct Parameters : public ParametersInterface {

    // Parameters for the search strategy mode
    PartitionStrategy pStrategy = DYNAMIC; // default dynamic partitioning
    DistanceMetric metric = EDIT;          // default edit distance (optimized)
    MappingMode mMode = BEST;              // Default mode is best-map
    SequencingMode sMode = SINGLE_END;     // Default mode is single end reads

    // Paired-end parameters
    Orientation orientation = FR;         // default orientation is FR
    length_t maxInsertSize = 500;         // default maximum insert size
    length_t minInsertSize = 0;           // default minimum insert size
    bool inferPairedEndParameters = true; // default infer the parameters
    bool discordantAllowed = false;       // default no discordant pairs allowed
    length_t nDiscordant = 100 * 1000;    // default number of discordant pairs

    // Parameters for the search scheme to use
    std::string searchScheme = "columba"; // default search scheme
    std::string custom =
        ""; // path to custom search scheme (overrides default search scheme)
    std::string dynamicSelectionPath =
        ""; // path to custom search scheme with dynamic selection (overrides
            // default search scheme)
    bool customDynamicSelection = true; // use dynamic selection with custom
                                        // search scheme

    // Input output parameters
    std::string firstReadsFile = "";  // path to first reads file
    std::string secondReadsFile = ""; // path to second reads file (optional)
    std::string base = "";            // path to the basename of the index
    std::string outputFile = "ColumbaOutput.sam"; // path to the output file
    bool outputIsSAM = true;                      // output in SAM format
    std::string logFile = "";                     // Path to the log file
    std::string command = ""; // The command that was used to run the program
    bool noUnmappedRecord = false; // Do not output the unmapped reads
    bool XATag = false;            // Add secondary alignments via XA tag
    bool reorder = false;          // Ensure order is same as original file

    bool doTrim = false;    // Trim the reads before aligning
    length_t trimStart = 0; // The start position of the bases to keep
    length_t trimEnd = 0;   // The end position of the bases to keep

    bool verbose =
#ifdef DEVELOPER_MODE
        true; // Print verbose stats
#else
        false; // Print verbose stats
#endif

    // Numerical parameters
    length_t maxDistance = 0; // the maximum allowed distance (for ALL mode)
    length_t nThreads = 1;    // the number of threads to be used
    length_t kmerSize = 10;   // The size of k-mers in the hash table (used as
                              // seeds during partitioning)
    length_t minIdentity =
        95; // The minimum identity for alignments in BEST mode
    length_t strataAfterBest = 0; // the number of strata above the best stratum
    // to explore in BEST mode

#ifndef RUN_LENGTH_COMPRESSION
    bool noCIGAR = false; // By default calculate the CIGAR string
#else
    bool noCIGAR = true; // By default do not calculate the CIGAR string
#endif // RUN_LENGTH_COMPRESSION

    bool cigarBehaviourChanged = false; // If the CIGAR behaviour has changed

    // Non-bmove parameters
#ifndef RUN_LENGTH_COMPRESSION
    // Input/Output

    // Numerical
    length_t sparsenessFactor =
        DEFAULT_SPARSENESS; // the sparseness factor for the suffix array
    length_t inTextVerificationPoint =
        4; // the tipping point for in-text verification
#endif
    /**
     * Get the command that was used to run the program
     */
    static std::string getCommand(int argc, char** argv) {
        std::string command;
        for (int i = 0; i < argc; ++i) {
            command += argv[i];
            if (i < argc - 1) {
                command += " ";
            }
        }
        return command;
    }

    /**
     * Process the optional arguments given by the user
     * @param argc The number of arguments
     * @param argv The arguments
     * @return The parameters struct
     */
    static Parameters processOptionalArguments(int argc, char** argv);

    /**
     * Create a search strategy based on the parameters
     * @param index Pointer to he index to search in
     * @return The search strategy
     */
    std::unique_ptr<SearchStrategy> createStrategy(IndexInterface& index) const;

    /**
     * Get the maximum distance or identity based on the mode
     * @return The maximum distance or minimum identity
     */
    length_t getMaxOrIdentity() const {
        return (mMode == BEST) ? minIdentity : maxDistance;
    }
};

class ParameterOption : public Option {
  public:
    ParameterOption(const std::string& shortOpt, const std::string& longOpt,
                    bool hasArgument, ArgumentType argumentType,
                    OptionType optionType)
        : Option(shortOpt, longOpt, hasArgument, argumentType, optionType) {
    }

    void process(const std::string& arg,
                 ParametersInterface& params) const override {
        Parameters& p = dynamic_cast<Parameters&>(params);
        process(arg, p);
    };

    virtual void process(const std::string& arg, Parameters& params) const = 0;
};
#endif
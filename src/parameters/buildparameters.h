#ifndef BUILD_PARAMETERS_H
#define BUILD_PARAMETERS_H

#include "../definitions.h"
#include "parameters.h"

struct BuildParameters : public ParametersInterface {

    std::string baseFN;
    std::vector<std::string> fastaFiles;
    int seedLength = DEFAULT_SEED_LENGTH;
#ifdef RUN_LENGTH_COMPRESSION
    bool pfp = false;
    bool preprocessOnly = false;
#else
    int sparsenessFactor = 4;
    bool allSparsenessFactors = false;

#endif
    bool valid = false;

    static bool parse(int argc, char** argv, BuildParameters& params) {
        params = BuildParameters::parse(argc, argv);
        return params.valid;
    }

    static BuildParameters parse(int argc, char** argv);

    static bool validFastaExtension(std::string& fastaFile);
    static void removeFastaExtension(std::string& baseFN);

    static std::array<std::string, 6> allowedExtensionsFasta;

    static void showUsage() {
        printHelp();
    }
};

class BuildParameterOption : public Option {
  public:
    BuildParameterOption(const std::string& shortOpt,
                         const std::string& longOpt, bool hasArgument,
                         ArgumentType argumentType, OptionType optionType)
        : Option(shortOpt, longOpt, hasArgument, argumentType, optionType) {
    }

    void process(const std::string& arg,
                 ParametersInterface& params) const override {
        BuildParameters& p = dynamic_cast<BuildParameters&>(params);
        process(arg, p);
    };

    virtual void process(const std::string& arg,
                         BuildParameters& params) const = 0;
};

#endif
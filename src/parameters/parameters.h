#ifndef OPTION_H
#define OPTION_H

#include "../logger.h"
#include <algorithm> // for copy_if
#include <cassert>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric> // for accumulate
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

/**
 * Enum to represent the different types of arguments
 */
enum ArgumentType { STRING, INTEGER, NONE };

/**
 * Enum to represent the different types of options
 */
enum OptionType {
    REQUIRED,
    PARALLELIZATION,
    ALIGNMENT,
    PAIRED_END_OPTION,
    OUTPUT,
    BUILD_INPUT_FILE,
    BUILD_OUTPUT_FILE,
    ADVANCED,
    HELP,

    NUM_OPTION_TYPES // helper to get the number of option types
};

/**
 * Function to convert an option type to a string for the help function
 */
inline std::string optionTypeToString(OptionType type) {
    switch (type) {
    case REQUIRED:
        return "Required";
    case PARALLELIZATION:
        return "Parallelization";
    case ALIGNMENT:
        return "Alignment";
    case PAIRED_END_OPTION:
        return "Paired-end";
    case OUTPUT:
        return "Output";
    case ADVANCED:
        return "Advanced";
    case HELP:
        return "Help";
    case BUILD_INPUT_FILE:
        return "Input file";
    case BUILD_OUTPUT_FILE:
        return "Index name";
    default:
        return "Unknown type";
    }
}

// Forward declaration
struct ParametersInterface;

/**
 * A class to represent a command line option.
 */
class Option {
  protected:
    std::string shortOpt; ///<  The short flag for the option
    std::string longOpt;  ///< The long flag for the option

    bool hasArgument;          ///< If the option has an argument
    ArgumentType argumentType; ///< The type of the argument (NONE if no)

    OptionType optionType; ///< The type of the option

    bool noPrint = false; ///< If the option should not be printed in the help

    /**
     * @brief Get the description of the option
     */
    virtual std::string getDescription() const = 0;

    /**
     * @brief Convert an argument type to a string
     */
    static std::string argumentTypeToString(ArgumentType argumentType) {
        switch (argumentType) {
        case STRING:
            return "STR";
        case INTEGER:
            return "INT";
        case NONE:
            return "   ";
        default:
            return "   ";
        }
    }

  public:
    /**
     * @brief Process the argument of the option
     *
     * @param arg the argument
     * @param params the parameters struct, if the argument is valid the
     * relevant parameters are updated
     */
    virtual void process(const std::string& arg,
                         ParametersInterface& params) const = 0;

    /**
     * @brief Construct a new Option object
     *
     * @param shortOpt the short flag for the option
     * @param longOpt the long flag for the option
     * @param hasArgument if the option has an argument
     * @param argumentType the type of the argument (NONE if no argument)
     * @param optionType the type of the option
     */
    Option(const std::string& shortOpt, const std::string& longOpt,
           bool hasArgument, ArgumentType argumentType, OptionType optionType)
        : shortOpt(shortOpt), longOpt(longOpt), hasArgument(hasArgument),
          argumentType(argumentType), optionType(optionType) {
    }

    virtual ~Option() {
    } // Virtual destructor

    /**
     * @brief Check if the option has an argument
     * @return true if the option has an argument
     */
    bool hasArg() const {
        return hasArgument;
    }

    /**
     * @brief Check if the option requires an argument
     * @return true if the option requires an argument
     */
    virtual bool argumentRequired() const {
        return hasArgument;
    }

    /**
     * @brief return the message to print when the option is ignored due to
     * invalid configuration.
     * @return the message to print when the option is ignored
     */
    std::string ignoreMessage(const std::string& arg = "") const {
        std::string extra = "";
        if (arg.size() > 0)
            extra = " with argument \"" + arg + "\"";
        return "Ignoring -" + shortOpt + "/--" + longOpt + " option" + extra;
    }

    /**
     * @brief Get the short flag for the option
     * @return the short flag for the option
     */
    std::string getShortOpt() const {
        return shortOpt;
    }

    /**
     * @brief Get the long flag for the option
     * @return the long flag for the option
     */
    std::string getLongOpt() const {
        return longOpt;
    }

    /**
     * @brief Check if the option is required
     * @return true if the option is required
     */
    bool isRequired() const {
        return optionType == REQUIRED;
    }

    /**
     * @brief get the size of the option tag (for pretty printing)
     */
    size_t sizeOptionTag() const {
        return 1 + shortOpt.size() + 4 + longOpt.size();
    }

    /**
     * @brief get the size of the option argument type (for pretty printing)
     */
    size_t sizeOptionArgumentType() const {
        return argumentTypeToString(argumentType).size();
    }

    /**
     * @brief get the size of the option description (for pretty printing)
     */
    size_t sizeOptionDescription() const {
        return getDescription().size();
    }

    /**
     * @brief get option tag (for pretty printing)
     */
    std::string getOptionTag() const {
        return "  -" + shortOpt + ", --" + longOpt;
    }

    /**
     * @brief get option argument type (for pretty printing)
     */
    std::string getOptionArgumentType() const {
        return argumentTypeToString(argumentType);
    }

    /**
     * @brief get option description (for pretty printing)
     */
    std::string getOptionDescription() const {
        return getDescription();
    }

    /**
     * @brief get the option type
     */
    OptionType getOptionType() const {
        return optionType;
    }

    void printOption(size_t tagWidth, size_t typeWidth,
                     size_t descWidth) const {
        std::cout << std::left << std::setw(tagWidth) << getOptionTag()
                  << std::setw(typeWidth) << getOptionArgumentType();

        // Handle description to wrap within maxWidth
        std::string description = getOptionDescription();

        // remove leading spaces
        while (description.size() > 0 && description[0] == ' ') {
            description = description.substr(1);
        }

        while (description.size() > 0) {
            if (description.size() <= descWidth) {
                std::cout << description << std::endl;
                return;
            }

            // find last space before end
            size_t end = std::min(description.size(), descWidth);
            size_t lastSpace = description.find_last_of(' ', end);
            if (lastSpace == std::string::npos) {
                lastSpace = end;
            }

            // split into line and remaining
            std::string line = description.substr(0, lastSpace);
            description = description.substr(lastSpace);

            std::cout << line << std::endl;

            // remove leading spaces
            while (description.size() > 0 && description[0] == ' ') {
                description = description.substr(1);
            }
            if (description.size() > 0) {
                std::cout << std::setw(tagWidth + typeWidth) << "";
            }
        }
    }

    bool toPrint() const {
        return !noPrint;
    }
};

/**
 * Function to add an option to the options map.
 */
inline void addOption(Option* option,
                      std::unordered_map<std::string, Option*>& options) {
    // assert not in map yet
    assert(options.find(option->getShortOpt()) == options.end());
    assert(options.find(option->getLongOpt()) == options.end());
    options[option->getShortOpt()] = option;
    options[option->getLongOpt()] = option;
}

/**
 * ParametersInterface structure to hold program parameters and provide help
 * functionality.
 */
struct ParametersInterface {

    virtual ~ParametersInterface() = default;
    static void printHeader(); ///< Print the header of the program help
    static void printFooter() {
        std::cout << "Report bugs to "
                  <<
#ifdef RUN_LENGTH_COMPRESSION
            "lore.depuydt@ugent.be and "
                  <<
#endif
            "luca.renders@ugent.be" << std::endl;
    }
    static void printHelp() {
        printHeader();
        const size_t maxWidth = 100;

        // get the maximum value of sizeCol1 of all options
        size_t tagWidth = 0, typeWidth = 0;
        for (const auto& option : options) {
            tagWidth = std::max(tagWidth, option->sizeOptionTag());
            typeWidth = std::max(typeWidth, option->sizeOptionArgumentType());
        }

        tagWidth += 3;  // padding
        typeWidth += 3; // padding

        // Determine the maximum substring length that fits within
        // remaining space
        size_t descWidth = maxWidth - tagWidth - typeWidth;

        // print each option by option type
        for (uint16_t i = 0; i < OptionType::NUM_OPTION_TYPES; i++) {
            std::vector<std::shared_ptr<Option>> toPrint;
            toPrint.reserve(options.size());

            std::copy_if(options.begin(), options.end(),
                         std::back_inserter(toPrint),
                         [i](const std::shared_ptr<Option>& opt) {
                             return opt->getOptionType() == (OptionType)i &&
                                    opt->toPrint();
                         });

            if (!toPrint.empty()) {
                std::cout << optionTypeToString((OptionType)i) << " options:\n";
                for (const auto& opt : toPrint) {
                    opt->printOption(tagWidth, typeWidth, descWidth);
                }
                std::cout << "\n";
            }
        }
    }

    const static std::vector<std::shared_ptr<Option>> options;

    void parse(int argc, char** argv) {
        std::unordered_map<std::string, Option*> optionMap;
        std::unordered_set<Option*> requiredOptions;

        std::unordered_set<std::string> usedOptions;

        for (const auto& option : ParametersInterface::options) {
            addOption(option.get(), optionMap);
            if (option->isRequired()) {
                requiredOptions.insert(option.get());
            }
        }

        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];
            if (arg.size() < 2 || arg[0] != '-') {
                logger.logWarning("Invalid option: " + arg);
                continue;
            }

            std::string optionName =
                arg[1] == '-' ? arg.substr(2) : arg.substr(1);
            auto optionIt = optionMap.find(optionName);
            if (optionIt == optionMap.end()) {
                logger.logWarning("Invalid option: " + optionName);
                continue;
            }

            Option* option = optionIt->second;

            // Construct the argument
            std::string argument = "";
            while (i < argc - 1 && argv[i + 1][0] != '-') {
                argument += argv[++i];
                argument += " ";
            }

            // remove trailing space
            if (!argument.empty() && argument.back() == ' ') {
                argument.pop_back();
            }

            if (option->argumentRequired() && argument.empty()) {
                logger.logWarning("Missing argument for option: " + arg + ". " +
                                  option->ignoreMessage());
                continue;
            }

            if (!option->hasArg() && !argument.empty()) {
                logger.logWarning("Ignoring argument \"" + argument +
                                  "\" for option: " + arg);
                argument = "";
            }

            // check if the option was already used
            if (usedOptions.find(option->getShortOpt()) != usedOptions.end()) {
                logger.logWarning(
                    "Option " + option->getOptionTag() +
                    " is used more than once. Ignoring all but the first.");
                continue;
            }

            option->process(argument, *this);
            usedOptions.insert(option->getShortOpt());

            if (option->isRequired() &&
                requiredOptions.find(option) != requiredOptions.end()) {
                requiredOptions.erase(option);
            }
        }

        if (!requiredOptions.empty()) {
            std::string missingOptions = std::accumulate(
                requiredOptions.begin(), requiredOptions.end(), std::string(),
                [](const std::string& acc, const Option* o) {
                    return acc + "\n\t" + o->getOptionTag() + " " +
                           o->getOptionDescription();
                });
            throw std::runtime_error("Missing required options:" +
                                     missingOptions);
        }
    }
};

#endif

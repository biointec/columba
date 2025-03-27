/******************************************************************************
 *  Columba: Approximate Pattern Matching using Search Schemes                *
 *  Copyright (C) 2020-2024 - Luca Renders <luca.renders@ugent.be> and        *
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

#ifndef LOGGER_H
#define LOGGER_H

#include <chrono>  // for steady_clock, time_point
#include <cstdio>  // for stderr
#include <fstream> // for ofstream, ostream
#include <mutex>   // for mutex
#include <sstream> // for stringstream
#include <string>  // for string

/**
 * @brief The Logger class
 *
 * This class is used to log messages to a file. It is thread-safe.
 */
class Logger {
  public:
    /**
     * @brief Construct a new Logger object.
     *
     * @param filename The name of the file to log to.
     */
    Logger(const std::string& filename);

    /**
     * Default constructor. Sets the output stream to std::cout.
     */
    Logger();

    /**
     * @brief Set the log file to a new file.
     *
     * @param filename The name of the file to log to.
     */
    void setLogFile(const std::string& filename);

    /**
     * @brief Set the verbosity of the logger.
     * @param verbose Whether to log verbose messages.
     */
    void setVerbose(bool verbose) {
        this->verbose = verbose;
    }

    /**
     * @brief Log an info message to the log file.
     *
     * @param message The message to log.
     */
    void logInfo(const std::string& message) {
        logWithTimestamp("[INFO]", message);
    }

    /**
     * @brief Log an info message to the log file.
     *
     * The string stream is cleared after logging.
     * @param stream The stream to log.
     */
    void logInfo(std::stringstream& stream) {
        logInfo(stream.str());
        // clear the stream
        stream.str(std::string());
    }

    /**
     * @brief Log a verbose message to the log file.
     *
     * @param message The message to log.
     */
    void logVerbose(const std::string& message) {
        if (verbose)
            logWithTimestamp("[VERBOSE]", message);
    }

    /**
     * @brief Log a verbose message to the log file.
     *
     * The string stream is cleared after logging.
     * @param stream The stream to log.
     */
    void logVerbose(std::stringstream& stream) {
        logVerbose(stream.str());
        // clear the stream
        stream.str(std::string());
    }

    /**
     * @brief Log a developer message to the log file.
     *
     * @param message The developer message to log.
     */
    void logDeveloper(const std::string& message) {
#ifdef DEVELOPER_MODE
        logWithTimestamp("[DEVELOPER]", message);
#endif
    }

    /**
     * @brief Log a developer message to the log file.
     *
     * The string stream is cleared after logging.
     *
     * @param stream The stream to log.
     */
    void logDeveloper(std::stringstream& stream) {
        logDeveloper(stream.str());

        // clear the stream
        stream.str(std::string());
    }

    /**
     * @brief Log an error message to the log file.
     *
     * @param message The error message to log.
     */
    void logError(const std::string& message) {
        logWithTimestamp("[ERROR]", message, true);
    }

    /**
     * @brief Log an error message to the log file.
     *
     * The string stream is cleared after logging.
     *
     * @param stream The stream to log.
     */
    void logError(std::stringstream& stream) {
        logError(stream.str());

        // clear the stream
        stream.str(std::string());
    }

    /**
     * @brief Log a warning message to the log file.
     *
     * @param message The warning message to log.
     */
    void logWarning(const std::string& message) {
        logWithTimestamp("[WARNING]", message);
    }

    /**
     * @brief Log a warning message to the log file.
     *
     * The string stream is cleared after logging.
     *
     * @param stream The stream to log.
     */
    void logWarning(std::stringstream& stream) {
        logWarning(stream.str());

        // clear the stream
        stream.str(std::string());
    }

    ~Logger() {
        // wait for the mutex to be released
        std::lock_guard<std::mutex> lock(mutex);

        // close the file stream if it is open
        if (fileStream.is_open()) {
            fileStream.close();
        }
    }

  private:
    std::ostream* outputStream; // The output stream
    std::ofstream fileStream;   // The file stream
    std::mutex mutex;           // The mutex to ensure thread-safety
    std::chrono::time_point<std::chrono::steady_clock>
        startTime; // Start time of the logger

    bool verbose = false; // whether to log verbose messages

    /**
     * @brief Log a message with a timestamp.
     *
     * @param level The log level (e.g., [INFO], [ERROR], [WARNING]).
     * @param message The message to log.
     * @param stderrBool Whether to log to stderr. If this is false logging will
     * happen to the log file if it is set or to stdout.  (default: false)
     */
    void logWithTimestamp(const std::string& level, const std::string& message,
                          bool stderrBool = false);
};

// Declare a global instance of the Logger
extern Logger logger;

#endif // LOGGER_H

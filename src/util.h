/******************************************************************************
 *   Copyright (C) 2014 - 2021 Jan Fostier (jan.fostier@ugent.be)             *
 *   This file is part of Detox                                               *
 *                                                                            *
 *   This program is free software; you can redistribute it and/or modify     *
 *   it under the terms of the GNU Affero General Public License as published *
 *   by the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                      *
 *                                                                            *
 *   This program is distributed in the hope that it will be useful,          *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *   GNU General Public License for more details.                             *
 *                                                                            *
 *   You should have received a copy of the GNU Affero General Public License *
 *   along with this program; if not, see <https://www.gnu.org/licenses/>.    *
 ******************************************************************************/

#ifndef UTIL_H
#define UTIL_H

#include <algorithm> // for max, abs
#include <chrono>    // for system_clock, time_point
#include <ctime>     // for size_t
#include <fstream>   // for ifstream, ios
#include <mutex>     // for mutex
#include <string>    // for string
// ============================================================================
// DEFINITIONS
// ============================================================================

#define MAX_TIMERS 16

// ============================================================================
// WORK LOAD BALANCER
// ============================================================================

class WorkLoadBalancer {
  private:
    size_t idxBegin;  // begin index to process
    size_t idxEnd;    // end index to process
    size_t chunkSize; // target chunk size
    size_t idxNext;   // next node offset

    std::mutex mutex;    // mutex
    std::string message; // message to print during progress

  public:
    /**
     * Get a chunk of nodes (thread-safe)
     * @param chunkBegin Begin of chunk to handle
     * @param chunkEnd End of chunk to handle
     * @return False if there is no work left
     */
    bool getChunk(size_t& chunkBegin, size_t& chunkEnd);

    /**
     * Default constructor
     * @param idxBegin Begin index to process
     * @param idxEnd End index to process
     * @param chunkSize Maximum size per chunk
     * @param message Message to print during progress
     */
    WorkLoadBalancer(size_t idxBegin, size_t idxEnd, size_t chunkSize,
                     const std::string& message = "Processing")
        : idxBegin(idxBegin), idxEnd(idxEnd), chunkSize(chunkSize),
          idxNext(idxBegin), message(message) {
    }
};

// ============================================================================
// UTILITY CLASS WITH DIVERSE AUXILIARY ROUTINES
// ============================================================================

class Util {
  private:
    static double prevProgress;
    static int currentTimer;
    static std::chrono::time_point<std::chrono::system_clock>
        startTime[MAX_TIMERS];

  public:
    /**
     * Create a string with a human readable version of a time period
     * @param time Time period (expressed in s)
     * @return String with a human readable version of a time period
     */
    static std::string humanTime(double time);

    /**
     * Create a string with a human readable version of a (file) size
     * @param size Size (expressed in bytes)
     * @return String with a human readable version of a (file) size
     */
    static std::string humanSize(size_t size);

    /**
     * Write a progress indicator to stdout
     * @param str String to print before percentage
     * @param curr Current progress thus far
     * @param max Maximum progress
     */
    static void progress(const std::string& str, double curr, double max);

    /**
     * Write a progress indicator to stdout followed by "(100.0%)"
     * @param str String to print before percentage
     * @param elapsed Elapsed amount of time (optional)
     */
    static void progressEnd(const std::string& str, double elapsed = -1.0);

    /**
     * Start a chronometer
     */
    static void startChrono();

    /**
     * Stop the chronometer
     * @return The time in
     */
    static double stopChrono();

    /**
     * Stop the chronometer and return a human readable string
     * @return A human readable string containing the elapsed time
     */
    static std::string stopChronoStr() {
        return humanTime(stopChrono());
    }

    /**
     * Get a string containing the date and time
     * @return string containing date and time
     */
    static std::string getDateTime();

    /**
     * Check whether a file exists
     * @param filename
     * @return True of false
     */
    static bool fileExists(const std::string& filename) {
        std::ifstream file(filename.c_str(), std::ios::in);
        bool OK = file.good();
        file.close();
        return OK;
    }

    /**
     * Compute the relative distance between two numbers
     * @param a number one
     * @param b number two
     * @return abs((a-b)/min(a,b)) OR 0 when both a=0 and b=0
     */
    static double relDiff(double a, double b) {
        if (a == 0.0 && b == 0.0)
            return 0.0;
        return std::max(std::abs((a - b) / a), std::abs((a - b) / b));
    }

    /**
     * Compute the percentage of two size_t numbers
     * @param nom Nominator
     * @param den Denominator
     * @return The percentage
     */
    static double toPercentage(size_t nom, size_t den) {
        if (den == 0)
            return 0;
        return 100.0 * double(nom) / double(den);
    }
};

#endif

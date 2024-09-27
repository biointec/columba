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

#include "logger.h"
#include <iomanip>  // for operator<<, setfill, setw
#include <iostream> // for cout

using namespace std;

Logger::Logger(const string& filename)
    : startTime(chrono::steady_clock::now()) {
    fileStream.open(filename);
    outputStream = &fileStream;
}

Logger::Logger() : startTime(chrono::steady_clock::now()) {
    outputStream = &cout;
}

void Logger::setLogFile(const string& filename) {
    lock_guard<std::mutex> lock(mutex);
    if (fileStream.is_open()) {
        fileStream.close();
    }
    fileStream.open(filename);
    outputStream = &fileStream;
}

void Logger::logWithTimestamp(const string& level, const string& message,
                              bool stderrBool) {
    lock_guard<std::mutex> lock(mutex);
    auto now = chrono::steady_clock::now();
    auto elapsed = now - startTime;

    auto hours = chrono::duration_cast<chrono::hours>(elapsed).count();
    elapsed -= chrono::hours(hours);
    auto minutes = chrono::duration_cast<chrono::minutes>(elapsed).count();
    elapsed -= chrono::minutes(minutes);
    auto seconds = chrono::duration_cast<chrono::seconds>(elapsed).count();
    elapsed -= chrono::seconds(seconds);
    auto milliseconds =
        chrono::duration_cast<chrono::milliseconds>(elapsed).count();

    ostringstream timestamp;

    if (hours > 0 || minutes > 0) {
        if (hours > 0) {
            timestamp << setfill('0') << setw(2) << hours << ":";
        }
        timestamp << setfill('0') << setw(2) << minutes << ":";
    }

    timestamp << setfill('0') << setw(2) << seconds << "." << setfill('0')
              << setw(3) << milliseconds;

    (stderrBool ? std::cerr : *outputStream)
        << "[" << timestamp.str() << "] " << level << "\t" << message
        << std::endl;
}

// Define the global instance of the Logger
Logger logger; // Use default constructor with cout

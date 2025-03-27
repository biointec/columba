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

#include "util.h"

#include <cassert>  // for assert
#include <cstdint>  // for uint64_t
#include <ctime>    // for size_t, ctime, time_t
#include <iomanip>  // for operator<<, setprecision
#include <iostream> // for cout
#include <sstream>  // for ostringstream

using namespace std;
using namespace std::chrono;

// ============================================================================
// PARALLEL GRAPH CLASSES
// ============================================================================

bool WorkLoadBalancer::getChunk(size_t& chunkBegin, size_t& chunkEnd) {
    // lock the mutex (RAII)
    std::unique_lock<std::mutex> lock(mutex);

    // no more work: send termination signal
    if (idxNext >= idxEnd)
        return false;

    // assign a workload
    chunkBegin = idxNext;
    idxNext += min(chunkSize, idxEnd - idxNext);
    chunkEnd = idxNext;

    size_t total = idxEnd - idxBegin;
    size_t done = idxNext - idxBegin;
    double perc = 100.0 * (double)done / (double)total;

    // update message (only if non-empty)
    if (!message.empty()) {
        cout << fixed << setprecision(1) << message << " (" << perc << "%)";
        if (idxNext >= idxEnd) {
            cout << endl;
        } else {
            cout << "\r";
            cout.flush();
        }
    }

    return true;
}

// ============================================================================
// UTILITY CLASS WITH DIVERSE AUXILIARY ROUTINES
// ============================================================================

time_point<system_clock> Util::startTime[MAX_TIMERS];
int Util::currentTimer = 0;
double Util::prevProgress = 0.0;

string Util::humanTime(double time) {
    uint64_t timeInt = uint64_t(time);

    uint64_t days = timeInt / 86400;
    timeInt -= days * 86400;
    uint64_t hours = timeInt / 3600;
    timeInt -= hours * 3600;
    uint64_t min = timeInt / 60;
    timeInt -= min * 60;
    uint64_t sec = timeInt;
    uint64_t ms = uint64_t((time - timeInt) * 1000);

    ostringstream iss;
    if (days > 0) {
        iss << days << "d " << hours << "h " << min << "min";
    } else if (hours > 0) {
        iss << hours << "h " << min << "min " << sec << "s";
    } else if (min > 0) {
        iss << min << "min " << sec << "s";
    } else if (sec > 0) {
        iss << sec << "s " << ms << "ms";
    } else {
        iss << ms << "ms";
    }

    return iss.str();
}

string Util::humanSize(size_t size) {
    ostringstream iss;
    iss.precision(4);

    if (size > (1ull << 40))
        iss << (double)size / (1ull << 40) << " TB";
    else if (size > (1ull << 30))
        iss << (double)size / (1ull << 30) << " GB";
    else if (size > (1ull << 20))
        iss << (double)size / (1ull << 20) << " MB";
    else if (size > (1ull << 10))
        iss << (double)size / (1ull << 10) << " kB";
    else
        iss << size << " B" << endl;

    return iss.str();
}

void Util::progress(const std::string& str, double curr, double max) {
    double currProgress = Util::toPercentage(curr, max);
    cout.precision(2);

    // only write to stdout when enough progress has been made
    if (currProgress == 0.0 || ((currProgress - prevProgress) >= 1)) {
        cout << str << " (" << Util::toPercentage(curr, max) << "%)  \r";
        cout.flush();
        prevProgress = currProgress;
    }
}

void Util::progressEnd(const std::string& str, double elapsed) {
    if (elapsed > 0.0)
        cout << str << " (" << Util::humanTime(elapsed) << ")  \n";
    else
        cout << str << " (100%)  \n";

    prevProgress = 0.0;
}

void Util::startChrono() {
    // make sure we don't use too many timers
    assert(currentTimer < MAX_TIMERS);
    startTime[currentTimer] = system_clock::now();
    currentTimer++;
}

double Util::stopChrono() {
    // make sure stopChrono isn't called too often
    currentTimer--;
    assert(currentTimer >= 0);

    duration<double> elapsed = system_clock::now() - startTime[currentTimer];
    return (elapsed.count());
}

string Util::getDateTime() {
    time_t time = system_clock::to_time_t(system_clock::now());
    return string(ctime(&time));
}

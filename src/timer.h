/******************************************************************************
!@ This file is part of the OpenRBC code for protein-resolution red blood cells
!@ Copyright (c) 2016 Yu-Hang Tang, Lu Lu
!@ Licensed under the Apache License, Version 2.0 (the "License")
!@ you may not use this file except in compliance with the License.
!@ You may obtain a copy of the License at
!@             http://www.apache.org/licenses/LICENSE-2.0
!@ Unless required by applicable law or agreed to in writing, software
!@ distributed under the License is distributed on an "AS IS" BASIS,
!@ WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!@ See the License for the specific language governing permissions and
!@ limitations under the License.
******************************************************************************/
#ifndef TIMER_H_
#define TIMER_H_

#include <iostream>
#include <map>
#include <unordered_map>
#include <omp.h>
#include "display.h"

namespace openrbc {

// A single timer
struct Timer {
    inline void start() {
        _start_time = omp_get_wtime();
    }
    inline double get_start_time() {
        return _start_time;
    }
    inline double stop() {
        if ( _start_time )
            return _cumu_time += omp_get_wtime() - _start_time;
        else
            return 0.0;
    }
    inline double read() {
        return _cumu_time;
    }
    inline void reset() {
        _cumu_time = 0.;
    }
    Timer() : _start_time( 0. ), _cumu_time( 0. ) {}
private:
    double _start_time;
    double _cumu_time;
};

// Group of timers indexed by label
struct Timers {
    Timer    &    operator [] ( std::string const & label )       { return timers[label]; }
    Timer const & operator [] ( std::string const & label ) const { return timers[label]; }
    ~Timers() { report( true ); }
    void report( bool destructive = false ) {
        std::map<std::string, double> readings; // use map to do sort timers by label
        for ( auto & item : timers ) {
            readings[item.first] = item.second.read();
        }
        for ( auto & read : readings ) {
            display_timing( std::cout, read.second, read.first.c_str(), false );
        }
        if ( destructive ) timers.clear();
    }

protected:
    mutable std::unordered_map<std::string, Timer> timers;
};

}

#endif /* TIMER_H_ */

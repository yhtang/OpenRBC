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
#ifndef OPENRBC_DISPLAY_H_
#define OPENRBC_DISPLAY_H_

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>

namespace openrbc {

inline void display_header( std::ostream & out ) {
    static bool header = false;
    if ( !header ) {
        std::ios::fmtflags f( out.flags() );
        out << std::fixed << std::setprecision( 7 ) << std::setw( 10 ) << "Time" << '\t'
            << std::fixed << std::setw( 10 ) << "Temperature" << '\t'
            << std::fixed << std::setw( 10 ) << "Wall Time" << '\t'
            << std::fixed << std::setw( 10 ) << "Lost lipid" << '\t'
            << std::endl;
        out.flags( f );
        header = true;
    }
}

inline void display( std::ostream & out, double t, double temp, double wall_time ) {
    display_header( out );

    std::ios::fmtflags f( out.flags() );

    out << std::fixed << std::setprecision( 7 ) << std::setw( 10 ) << t << '\t'
        << std::fixed << std::setw( 10 ) << temp << '\t'
        << std::fixed << std::setw( 10 ) << wall_time << '\t'
        << std::fixed << std::setw( 10 ) << Service<Variable<int,0> >::call().value << '\t'
        << std::endl;

    out.flags( f );
}

inline void display_timing( std::ostream & out, double reading, const char msg[], bool human = true ) {
    std::ios::fmtflags f( out.flags() );

    std::string unit = "  s ";
    if ( human ) {
        if ( reading < 9e-4 ) {
            reading *= 1e6;
            unit = " us ";
        } else if ( reading < 9e-1 ) {
            reading *= 1e3;
            unit = " ms ";
        }
    }

    out << std::fixed << std::setprecision( 7 ) << std::setw( 12 )
        << reading << unit << "on " << msg << "\n";

    out.flags( f );
}

}

#endif

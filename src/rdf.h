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
#ifndef OPENRBC_RDF_H_
#define OPENRBC_RDF_H_

#include <iostream>
#include <fstream>

#include "aligned_array.h"
#include "math_vector.h"

namespace openrbc {

template<typename SCALAR>
void rdf( const AlignedArray<vector<SCALAR, 3>, true> & x ) {
    std::cout << "calculating rdf ..." << std::flush;

    AlignedArray<int> rdf;
    rdf.assign( 10000, 0 );
    const SCALAR dr = 0.1;

    for ( int i = 0; i < x.size(); ++i )
        for ( int j = i + 1; j < x.size(); ++j ) {
            const SCALAR dx = norm( x[i] - x[j] );
            rdf[static_cast<int>( dx / dr )] += 1;
        }

    std::ofstream fout( "rdf.txt" );
    for ( int i = 0; i < rdf.size(); ++i ) {
        const SCALAR r = dr * i;
        fout << r << ' ' << rdf[i] / ( 4.0 * 3.14 * r * r * dr ) << '\n';
    }
    fout.close();

    std::cout << "done." << std::endl;
}

template<typename SCALAR>
void rdf_1peak_r( const ProteContainer & c ) {
    const SCALAR dr = 0.5;
    const SCALAR r_cut = 25;

    AlignedArray<int> rdf;
    rdf.assign( static_cast<int>(r_cut / dr) + 1, 0 );

    for ( int i = 0; i < c.size(); ++i )
        for ( int j = i + 1; j < c.size(); ++j )
            if (c.type[i] == 4 && c.type[j] == 4) {
                const SCALAR dx = norm( c.x[i] - c.x[j] );
                if (dx < r_cut)
                    rdf[static_cast<int>( dx / dr )] += 1;
            }

    SCALAR max_r = 0;
    SCALAR max_d = 0;
    std::ofstream fout( "rdf.txt" );
    for ( int i = 1; i < rdf.size(); ++i ) {
        const SCALAR r = dr * i;
        const SCALAR d = rdf[i] / ( 4.0 * 3.14 * r * r * dr );
        fout << r << ' ' << d << '\n';

        if (d > max_d) {
            max_r = r;
            max_d = d;
        }
    }
    fout.close();

    std::cout << "actin: r(1st peak) = " << max_r << std::endl;
}

}

#endif

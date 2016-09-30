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
#ifndef RDF_H_
#define RDF_H_

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

}

#endif

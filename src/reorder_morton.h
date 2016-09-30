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
#ifndef REORDER_MORTON_H_
#define REORDER_MORTON_H_

#include <algorithm>
#include <vector>
#include "aligned_array.h"
#include "math_vector.h"
#include "forcefield.h"

namespace openrbc {

inline unsigned int bit_space3( unsigned int x ) {
    x = ( x | ( x << 12 ) ) & 0X00FC003FU;
    x = ( x | ( x <<  6 ) ) & 0X381C0E07U;
    x = ( x | ( x <<  4 ) ) & 0X190C8643U;
    x = ( x | ( x <<  2 ) ) & 0X49249249U;
    return x;
}

inline unsigned int interleave3( unsigned int i, unsigned int j, unsigned int k ) {
    return bit_space3( i ) | ( bit_space3( j ) << 1 ) | ( bit_space3( k ) << 2 ) ;
}

template<typename SCALAR>
inline unsigned int morton_encode( const SCALAR x, const SCALAR y, const SCALAR z, RTParameter const & param ) {
    return interleave3( static_cast<unsigned int>( 2 * x + param.bsize[0] ),
                        static_cast<unsigned int>( 2 * y + param.bsize[1] ),
                        static_cast<unsigned int>( 2 * z + param.bsize[2] ) );
}

template<typename VECTOR, std::size_t ALIGN>
void reorder_morton( AlignedArray<VECTOR, true, ALIGN> & pts, RTParameter const & param, LoadBalancer load = LoadBalancer() ) {
    Service<Timers>::call()["|reorder_morton"].start();

    static AlignedArray<int,    true> indices;
    static AlignedArray<unsigned int> values;

    const auto n = pts.size();
    indices.resize( n );
    values.resize( n );
    if ( !load.set() ) load.set_range( n );

    static std::vector<std::pair<int, int> > chunk_now, chunk_new;
    chunk_now.resize( load.ntd() );
    chunk_new.resize( ( chunk_now.size() + 1 ) / 2 );
    for ( int i = 0; i < load.ntd(); i++ ) chunk_now[i] = std::make_pair( load.beg( i ), load.end( i ) );

    auto const & comp = []( int i, int j ) { return values[i] < values[j]; };

    #pragma omp parallel
    {
        // build key
        auto beg = load.beg();
        auto end = load.end();

        for ( int i = beg; i < end; ++i ) {
            indices[i] = i;
            values[i] = morton_encode( pts[i][0], pts[i][1], pts[i][2], param );
        }

        // local std::sort
        std::sort( indices.data() + beg, indices.data() + end, comp );

        #pragma omp barrier

        // parallel merge sort
        int tid = omp_get_thread_num();

        while ( chunk_now.size() > 1 ) {
            if ( tid % 2 == 0 && tid < chunk_now.size() - 1 ) {
                auto beg = chunk_now[tid].first;
                int mid = chunk_now[tid].second;
                auto end = chunk_now[tid + 1].second;
                assert( mid == chunk_now[tid + 1].first );

                int i = beg, j = mid, k = 0;
                for ( int k = 0; k < end - beg; ++k ) {
                    if ( ( j >= end ) || ( i < mid && comp( indices[i], indices[j] ) ) ) {
                        indices.shadow( beg + k ) = indices[i++];
                    } else {
                        indices.shadow( beg + k ) = indices[j++];
                    }
                }
                chunk_new[ tid / 2 ] = std::make_pair( beg, end );
            } else if ( chunk_now.size() % 2 && tid == chunk_now.size() - 1 ) {
                auto beg = chunk_now.back().first;
                auto end = chunk_now.back().second;
                for ( int i = beg; i < end; ++i ) indices.shadow( i ) = indices[i];
                chunk_new.back() = std::make_pair( beg, end );
            }

            #pragma omp barrier
            #pragma omp master
            {
                indices.swap();
                chunk_now.swap( chunk_new );
                chunk_new.resize( ( chunk_now.size() + 1 ) / 2 );
            }
            #pragma omp barrier
        }

        // permute data according to sorting result
        for ( int i = beg; i < end; ++i ) pts.shadow( i ) = pts[ indices[i] ];
    }

    pts.swap();

    Service<Timers>::call()["|reorder_morton"].stop();
}

}

#endif

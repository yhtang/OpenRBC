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
#ifndef STRAY_CLEANUP_H_
#define STRAY_CLEANUP_H_

#include "aligned_array.h"
#include "container.h"
#include "math_vector.h"
#include "service.h"
#include "timer.h"
#include "util_numa.h"

namespace openrbc {

using namespace config;

// clean up lipid particles that fell off the membrane
void delete_lipid( LipidContainer & lipid, VoronoiDiagram const & voronoi, VCellList & cell_lipid, RTParameter const & param ) {
    Service<Timers>::call()["cleanup_stray"].start();

    static AlignedArray<int> keep;
    static AlignedArray<int> new_position;

    const auto n = lipid.size();
    keep.resize( n );
    new_position.resize( n + 1 );

    #pragma omp parallel
    {
        const int beg = voronoi.load_balancer().beg();
        const int end = voronoi.load_balancer().end();
        std::vector<real> dr2;
        for ( int i = beg; i < end; ++i ) {
            const auto c = voronoi.centroids[i];
            dr2.clear();
            // calculate cell median displacement
            for ( int j = cell_lipid.cell_start[i]; j < cell_lipid.cell_start[i + 1]; ++j ) {
                const vector<real, 3> pt { lipid.x[j][0], lipid.x[j][1], lipid.x[j][2] };
                dr2.push_back( normsq( pt - c ) );
            }
            std::sort( dr2.begin(), dr2.end() );
            // mark lipids too far away as to delete
            real threshold = dr2[ dr2.size() / 2 ] * param.stray_tolerance * param.stray_tolerance;
            for ( int j = cell_lipid.cell_start[i]; j < cell_lipid.cell_start[i + 1]; ++j ) {
                const vector<real, 3> pt { lipid.x[j][0], lipid.x[j][1], lipid.x[j][2] };
                keep[j] = normsq( pt - c ) < threshold ? 1 : 0;
            }
        }
    }

    new_position[0] = 0;
    std::partial_sum( keep.data(), keep.data() + n, new_position.data() + 1 );
    int size_new = new_position[ n ];

    if ( size_new < n ) {
        #pragma omp parallel
        {
            auto beg = Service<BalancerMap>::call()[ lipid.id() ].beg();
            auto end = Service<BalancerMap>::call()[ lipid.id() ].end();
            for ( auto i = beg; i < end; ++i ) {
                if ( keep[i] ) {
                    const int j = new_position[i];
                    lipid.x.shadow( j ) = lipid.x[i];
                    lipid.v.shadow( j ) = lipid.v[i];
                    lipid.n.shadow( j ) = lipid.n[i];
                    lipid.o.shadow( j ) = lipid.o[i];
                }
            }
        }

        lipid.resize( size_new );
        lipid.swap();

        cell_lipid.update( lipid, voronoi, param );

        Service<Variable<int, 0> >::call().value += n - size_new;
    }

    Service<Timers>::call()["cleanup_stray"].stop();
}

}

#endif

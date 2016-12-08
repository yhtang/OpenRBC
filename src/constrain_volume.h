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
#ifndef OPENRBC_CONSTRAIN_VOLUME_H_
#define OPENRBC_CONSTRAIN_VOLUME_H_

#include "integrate_nh.h"
#include "voronoi.h"
#include "runtime_parameter.h"
#include "math_vector.h"
#include "aligned_array.h"

namespace openrbc {

//-------------------------- Kernel --------------------------
void constrain_volume( LipidContainer & lipid, ProteContainer & prote, VoronoiDiagram & voronoi, VCellList & cell_lipid, VCellList & cell_prote, RTParameter & param, config::real target_volume, config::real strength ) {
    Service<Timers>::call()["volume_constraint"].start();

    using namespace config;

    vect global_center = 0;
    real volume = 0;

    static AlignedArray<vect> cell_center, cell_normal;
    cell_center.resize( voronoi.centroids.size() );
    cell_normal.resize( voronoi.centroids.size() );

    #pragma omp parallel
    {
        const int cell_beg = voronoi.load_balancer().beg();
        const int cell_end = voronoi.load_balancer().end();

        vect local_center = 0;
        for ( int i = cell_beg; i < cell_end; ++i ) {
            for ( int d = 0; d < 3; ++d ) cell_center[i][d] = voronoi.centroids[i][d];
            local_center += cell_center[i];
        }
        #pragma omp critical
        global_center += local_center;
        #pragma omp barrier
        auto center = global_center / voronoi.n_cells;

        real local_volume = 0;
        for ( int i = cell_beg; i < cell_end; ++i ) {
            for ( int j = cell_lipid.cell_start[i]; j < cell_lipid.cell_start[i + 1]; ++j ) cell_normal[i] += lipid.n[j];
            cell_normal[i] = normalize( cell_normal[i] );
            const vect distance = cell_center[i] - center;
            if ( dot( cell_normal[i], distance ) < 0 ) cell_normal[i] = -cell_normal[i];

            const real height = dot( distance, cell_normal[i] );
            local_volume += height * ( cell_lipid.cell_start[i + 1] - cell_lipid.cell_start[i] ) * 3.1415926 * 1.26 / 4.0 / 3.0 * 1e-6;
        }
        #pragma omp atomic
        volume += local_volume;
        #pragma omp barrier

        real f = strength * ( target_volume - volume ) / target_volume;
        for ( int i = cell_beg; i < cell_end; ++i ) {
            for ( int j = cell_lipid.cell_start[i]; j < cell_lipid.cell_start[i + 1]; ++j ) {
                lipid.f[j] += f * cell_normal[i] * ForceField::mass[ lipid.type[i] ];
            }
            for ( int j = cell_prote.cell_start[i]; j < cell_prote.cell_start[i + 1]; ++j ) {
                prote.f[j] += f * cell_normal[i] * ForceField::mass[ prote.type[i] ];
            }
        }
    }

    if ( param.nstep % 100 == 0 ) {
        fprintf( stdout, "target_volume=%lf, global_volume=%f volume_force=%f\n", target_volume, volume, strength * ( target_volume - volume ) / target_volume );
    }

    Service<Timers>::call()["volume_constraint"].stop();
}

}

#endif

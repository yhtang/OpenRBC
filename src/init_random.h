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
#ifndef OPENRBC_INIT_RANDOM_H_
#define OPENRBC_INIT_RANDOM_H_

#include <cmath>
#include "forcefield.h"
#include "rng.h"
#include "container.h"
#include "reorder_morton.h"
#include "util_numa.h"

namespace openrbc {

// randomly place particles on the surface of a shpere with num density rho
inline void init_random_sphere( LipidContainer & model, RTParameter const & param, const config::real r ) {
    Service<Timers>::call()["init_random_sphere"].start();

    constexpr config::real PI = 3.14159265358;
    const int num_particles = static_cast<int>( param.rho * 4 * PI * r * r );

    model.resize( num_particles );
    #pragma omp parallel
    {
        auto beg = Service<BalancerMap>::call()[ model.id() ].beg();
        auto end = Service<BalancerMap>::call()[ model.id() ].end();
        MT19937 rng( omp_get_thread_num() + beg * end );

        for ( int i = beg; i < end; ++i ) {
            const double theta = 2 * PI * rng.u01();
            const double u = 2 * rng.u01() - 1;
            model.x[i]    = 0;
            model.x[i][0] = r * std::sqrt( 1 - u * u ) * std::cos( theta );
            model.x[i][1] = r * std::sqrt( 1 - u * u ) * std::sin( theta );
            model.x[i][2] = r * u;
        }
    }

    reorder_morton( model.x, param );

    #pragma omp parallel
    {
        auto beg = Service<BalancerMap>::call()[ model.id() ].beg();
        auto end = Service<BalancerMap>::call()[ model.id() ].end();

        for ( int i = beg; i < end; ++i ) model.n[i] = normalize( model.x[i] );
    }

    model.v.assign( num_particles, real( 0 ), Service<BalancerMap>::call()[ model.id() ] );
    model.o.assign( num_particles, real( 0 ), Service<BalancerMap>::call()[ model.id() ] );

    Service<Timers>::call()["init_random_sphere"].stop();

    model.tag2idx.build_map( model );
}

}

#endif

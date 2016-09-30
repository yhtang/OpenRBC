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
#ifndef COMPUTE_BONDED_H_
#define COMPUTE_BONDED_H_

#include <cmath>
#include "container.h"
#include "runtime_parameter.h"
#include "omp.h"
#include "forcefield.h"

namespace openrbc {

using namespace config;

#ifdef EXPLICIT_SIMD

template<class CONTAINER>
void compute_bonded_simd( CONTAINER & model ) {
    Service<Timers>::call()["compute_bonded"].start();

    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        const int nthreads = omp_get_num_threads();
        const int chunk_size = model.size() / nthreads;
        const int pbegin = chunk_size * tid;
        const int pend = ( tid == nthreads - 1 ) ? model.size() : ( pbegin + chunk_size );

        for ( std::size_t l = 0; l < model.bonds.size(); ++l ) {
            const int p1 = model.tag2idx[model.bonds[l].i];
            const int p2 = model.tag2idx[model.bonds[l].j];

            const bool p1_in = p1 >= pbegin && p1 < pend;
            const bool p2_in = p2 >= pbegin && p2 < pend;
            if ( ( !p1_in ) && ( !p2_in ) ) continue;

            const auto type  = model.bonds[l].type;
            const vect r     = model.x[p2] - model.x[p1];
            const auto cur_r = vect( ForceField::K[type] ) * ( vect( 1.0f ) - vect( ForceField::r0[type] ) * _rsqrt( _normsq( r ) ) );
            const vect force = cur_r * r;

            if ( p1_in ) model.f[p1] += force;
            if ( p2_in ) model.f[p2] -= force;
        }
    }

    Service<Timers>::call()["compute_bonded"].stop();
}

#else

template<class CONTAINER>
void compute_bonded( CONTAINER & model ) {
    Service<Timers>::call()["compute_bonded"].start();

    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        const int nthreads = omp_get_num_threads();
        const int chunk_size = model.size() / nthreads;
        const int pbegin = chunk_size * tid;
        const int pend = ( tid == nthreads - 1 ) ? model.size() : ( pbegin + chunk_size );

        for ( std::size_t l = 0; l < model.bonds.size(); ++l ) {
            const int p1 = model.tag2idx[model.bonds[l].i];
            const int p2 = model.tag2idx[model.bonds[l].j];

            const bool p1_in = p1 >= pbegin && p1 < pend;
            const bool p2_in = p2 >= pbegin && p2 < pend;
            if ( ( !p1_in ) && ( !p2_in ) ) continue;

            const auto type = model.bonds[l].type;
            const auto dx = model.x[p2] - model.x[p1];
            const auto rinv = _rsqrt( normsq( dx ) );
            const auto cur_r = ForceField::K[type] * ( 1 - ForceField::r0[type] * rinv );
            const auto force = cur_r * dx;

            if ( p1_in ) model.f[p1] += force;
            if ( p2_in ) model.f[p2] -= force;
        }
    }

    Service<Timers>::call()["compute_bonded"].stop();
}

#endif

}

#endif /* COMPUTE_BONDED_H_ */

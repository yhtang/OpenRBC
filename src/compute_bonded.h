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
#ifndef OPENRBC_COMPUTE_BONDED_H_
#define OPENRBC_COMPUTE_BONDED_H_

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
    rt_assert( Service<BalancerMap>::call()[ model.id() + "-bonds" ].set(), "The OpenMP load scheduler for bonds was not set up before compute_bonded" );

    #pragma omp parallel
    {
        const int beg = Service<BalancerMap>::call()[ model.id() ].beg();
        const int end = Service<BalancerMap>::call()[ model.id() ].end();
        for ( std::size_t i = beg; i < end; ++i ) model.f.shadow( i ) = constant::zero;
    }

    #pragma omp parallel
    {
        const int bond_beg = Service<BalancerMap>::call()[ model.id() + "-bonds" ].beg();
        const int bond_end = Service<BalancerMap>::call()[ model.id() + "-bonds" ].end();
        const int particle_beg = Service<BalancerMap>::call()[ model.id() ].beg();
        const int particle_end = Service<BalancerMap>::call()[ model.id() ].end();

        for ( std::size_t l = bond_beg; l < bond_end; ++l ) {
            const int p1 = model.tag2idx[model.bonds[l].i];
            const int p2 = model.tag2idx[model.bonds[l].j];
            const auto type  = model.bonds[l].type;
            const vect r     = model.x[p2] - model.x[p1];
            const auto cur_r = vect( ForceField::K[type] ) * ( vect( 1.0f ) - vect( ForceField::r0[type] ) * _rsqrt( _normsq( r ) ) );
            const vect force = cur_r * r;

            if ( p1 >= particle_beg && p1 < particle_end ) model.f[p1] += force;
            else {
                #pragma omp atomic
                model.f.shadow( p1 )[0] += force[0];
                #pragma omp atomic
                model.f.shadow( p1 )[1] += force[1];
                #pragma omp atomic
                model.f.shadow( p1 )[2] += force[2];
            }
            if ( p2 >= particle_beg && p2 < particle_end ) model.f[p2] -= force;
            else {
                #pragma omp atomic
                model.f.shadow( p2 )[0] -= force[0];
                #pragma omp atomic
                model.f.shadow( p2 )[1] -= force[1];
                #pragma omp atomic
                model.f.shadow( p2 )[2] -= force[2];
            }
        }
    }

    #pragma omp parallel
    {
        const int beg = Service<BalancerMap>::call()[ model.id() ].beg();
        const int end = Service<BalancerMap>::call()[ model.id() ].end();
        for ( std::size_t i = beg; i < end; ++i ) model.f[i] += model.f.shadow( i );
    }

    Service<Timers>::call()["compute_bonded"].stop();
}

#else

template<class CONTAINER>
void compute_bonded( CONTAINER & model ) {
    Service<Timers>::call()["compute_bonded"].start();
    rt_assert( Service<BalancerMap>::call()[ model.id() + "-bonds" ].set(), "The OpenMP load scheduler for bonds was not set up before compute_bonded" );

    #pragma omp parallel
    {
        const int beg = Service<BalancerMap>::call()[ model.id() ].beg();
        const int end = Service<BalancerMap>::call()[ model.id() ].end();
        for ( std::size_t i = beg; i < end; ++i ) model.f.shadow( i ) = constant::zero;
    }

    #pragma omp parallel
    {
        const int bond_beg = Service<BalancerMap>::call()[ model.id() + "-bonds" ].beg();
        const int bond_end = Service<BalancerMap>::call()[ model.id() + "-bonds" ].end();
        const int particle_beg = Service<BalancerMap>::call()[ model.id() ].beg();
        const int particle_end = Service<BalancerMap>::call()[ model.id() ].end();

        for ( std::size_t l = bond_beg; l < bond_end; ++l ) {
            const int p1 = model.tag2idx[model.bonds[l].i];
            const int p2 = model.tag2idx[model.bonds[l].j];
            const auto type = model.bonds[l].type;
            const auto dx = model.x[p2] - model.x[p1];
            const auto rinv = norm_inv( dx );
            const auto cur_r = ForceField::K[type] * ( 1 - ForceField::r0[type] * rinv );
            const auto force = cur_r * dx;

            if ( p1 >= particle_beg && p1 < particle_end ) model.f[p1] += force;
            else {
                #pragma omp atomic
                model.f.shadow( p1 )[0] += force[0];
                #pragma omp atomic
                model.f.shadow( p1 )[1] += force[1];
                #pragma omp atomic
                model.f.shadow( p1 )[2] += force[2];
            }
            if ( p2 >= particle_beg && p2 < particle_end ) model.f[p2] -= force;
            else {
                #pragma omp atomic
                model.f.shadow( p2 )[0] -= force[0];
                #pragma omp atomic
                model.f.shadow( p2 )[1] -= force[1];
                #pragma omp atomic
                model.f.shadow( p2 )[2] -= force[2];
            }
        }
    }

    #pragma omp parallel
    {
        const int beg = Service<BalancerMap>::call()[ model.id() ].beg();
        const int end = Service<BalancerMap>::call()[ model.id() ].end();
        for ( std::size_t i = beg; i < end; ++i ) model.f[i] += model.f.shadow( i );
    }

    Service<Timers>::call()["compute_bonded"].stop();
}

#endif

}

#endif

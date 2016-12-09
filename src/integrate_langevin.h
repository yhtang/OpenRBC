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
#ifndef OPENRBC_INTEGRATE_LANGEVIN_H_
#define OPENRBC_INTEGRATE_LANGEVIN_H_

#include <cmath>
#include <iostream>
#include <omp.h>

#include "integrate_nh.h"

namespace openrbc {

using namespace config;

struct verlet_langevin_initial {
    static std::string name() { return "verlet_langevin_initial"; }

    verlet_langevin_initial( RTParameter & p ) : parameter( p ) {}
    ~verlet_langevin_initial() {}

    template<class CONTAINER> void operator () ( CONTAINER & container, int beg, int end ) {
        const real dt        = parameter.dt;
        const real inertia   = constant::one; // 2/5 * m * R^2
        const real dt_over_m = dt / inertia;

        for ( int i = beg; i < end; ++i ) {
            // Position update
            // v = (v + 0.5 * f / m * dt) / (1 + 0.5 * dt * zeta)
            container.v[i] += constant::half * dt * container.f[i] / ForceField::mass[container.type[i]];
            // x = x + v * dt
            container.x[i] += container.v[i] * dt;

            // Rotation update
            // omega = omega + 0.5 * T / I * dt
            container.o[i] += constant::half / inertia * dt * container.t[i];

            // mu = mu + (omega X mu) * dt
            container.n[i] = normalize( container.n[i] + cross( container.o[i], container.n[i] ) * dt );

            // Clear force
            container.f[i] = constant::zero;
            container.t[i] = constant::zero;

        }
    }

    RTParameter & parameter;
};

struct verlet_langevin_final {
    static std::string name() { return "verlet_langevin_final"; }

    verlet_langevin_final( RTParameter & p ) : parameter( p ) {}
    ~verlet_langevin_final() {}

    template<class CONTAINER> void operator () ( CONTAINER & container, int beg, int end ) {
        const real dt        = parameter.dt;
        const real inertia   = constant::one; // 2/5 * m * R^2
        const real dt_over_m = dt / inertia;

        real gamma[ ForceField::n_type ], sigma[ ForceField::n_type ];
        for ( int i = 0; i < ForceField::n_type; ++i ) {
            gamma[i] = 6.0 * M_PI * parameter.eta * ForceField::radius[ i ];
            sigma[i] = std::sqrt( 2 * parameter.kBT * gamma[i] ) * std::sqrt( 3.0 / parameter.dt );
        }

        MT19937 & rng = parameter.prng[ omp_get_thread_num() ];

        for ( int i = beg; i < end; ++i ) {

            container.t[i] = cross( container.n[i], container.t[i] );

            // v = v + 0.5 * dt * f / m;
            auto type = container.type[i];
            vect r { rng.u11(), rng.u11(), rng.u11(), 0 };
            container.f[i] -= gamma[type] * container.v[i] + sigma[type] * r;
            container.v[i] += constant::half * dt * container.f[i] / ForceField::mass[container.type[i]];

            // omega = omega + 0.5 * T / I * dt
            container.o[i] += constant::half * dt_over_m * container.t[i];
        }
    }

    RTParameter & parameter;
};

struct verlet_langevin {
    static std::string name() { return "verlet_langevin"; }

    verlet_langevin( RTParameter & p ) : parameter( p ) {}
    ~verlet_langevin() {}

    template<class CONTAINER> void operator () ( CONTAINER & container, int beg, int end ) {
        const real dt        = parameter.dt;
        const real inertia   = constant::one; // 2/5 * m * R^2
        const real dt_over_m = dt / inertia;

        real gamma[ ForceField::n_type ], sigma[ ForceField::n_type ];
        for ( int i = 0; i < ForceField::n_type; ++i ) {
            gamma[i] = 6.0 * M_PI * parameter.eta * ForceField::radius[ i ];
            sigma[i] = std::sqrt( 2 * parameter.kBT * gamma[i] ) * std::sqrt( 3.0 / parameter.dt );
        }

        MT19937 & rng = parameter.prng[ omp_get_thread_num() ];
        vector<uint, vect::d>
        x { rng.uint(), rng.uint(), rng.uint(), 0 },
        y { rng.uint(), rng.uint(), rng.uint(), 0 },
        z { rng.uint(), rng.uint(), rng.uint(), 0 },
        w { rng.uint(), rng.uint(), rng.uint(), 0 };

        for ( int i = beg; i < end; ++i ) {
            // ROTATION
            container.t[i]  = cross( container.n[i], container.t[i] );
            container.o[i] += dt_over_m * container.t[i];
            container.n[i]  = normalize( container.n[i] + cross( container.o[i], container.n[i] ) * dt );
            container.t[i]  = constant::zero;

            // xorshift128
            auto t = x;
            t ^= t << 11;
            t ^= t >> 8;
            x = y; y = z; z = w;
            w ^= w >> 19;
            w ^= t;
            vect r = uint2u11( w );

            // TRANSLATION
            auto type = container.type[i];
            container.f[i] -= gamma[type] * container.v[i] + sigma[type] * r;
            container.v[i] += container.f[i] * ( dt / ForceField::mass[container.type[i]] );
            container.x[i] += container.v[i] * dt;
            container.f[i]  = constant::zero;
        }
    }

    RTParameter & parameter;
};

}

#endif

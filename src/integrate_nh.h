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
#ifndef OPENRBC_INTEGRATE_NH_H_
#define OPENRBC_INTEGRATE_NH_H_

#include <cmath>
#include <iostream>
#include <omp.h>

#include "container.h"
#include "forcefield.h"

namespace openrbc {

using namespace config;

//-------------------------- Driver --------------------------
template<class KERNEL, class ... CONTAINERS>
inline void integrate( KERNEL && kernel, CONTAINERS & ... containers ) {
    Service<Timers>::call()[ kernel.name() ].start();

    #pragma omp parallel
    for_each_container( kernel, containers... );

    Service<Timers>::call()[ kernel.name() ].stop();
}

template<class KERNEL, class CONTAINER, class ...MORE_CONTAINERS>
inline void for_each_container( KERNEL && kernel, CONTAINER & container, MORE_CONTAINERS & ... more_containers ) {
    action( kernel, container );
    for_each_container( kernel, more_containers... );
}

template<class KERNEL, class CONTAINER>
inline void for_each_container( KERNEL && kernel, CONTAINER & container ) {
    action( kernel, container );
}

template<class KERNEL, class CONTAINER>
inline void action( KERNEL && kernel, CONTAINER & container ) {
    auto beg = Service<BalancerMap>::call()[ container.id() ].beg();
    auto end = Service<BalancerMap>::call()[ container.id() ].end();
    kernel( container, beg, end );
}

//-------------------------- Kernel --------------------------
struct clear_force {
    static std::string name() { return "clear_force"; }

    template<class CONTAINER> void operator () ( CONTAINER & container, int beg, int end ) const {
        for ( int i = beg; i < end; ++i ) {
            container.f[i] = constant::zero;
            container.t[i] = constant::zero;
        }
    }
};

struct verlet_nh_update {
    static std::string name() { return "verlet_nh_update"; }

    verlet_nh_update( RTParameter & p ) : parameter( p ) {}
    // constructor serve as a finalization hook
    ~verlet_nh_update() {
        if ( !parameter.Q ) parameter.set_Q( n );
        parameter.zeta += constant::half * parameter.dt / parameter.Q * ( ke - constant::onehalf * n * parameter.kBT );
    }

    template<class CONTAINER> void operator () ( CONTAINER & container, int beg, int end ) {
        // zeta = zeta + 0.5 * dt / Q * (Sum(0.5*m*v^2) - 0.5*(3N+1)*kBT)
        double ke_local = 0;
        for ( int i = beg; i < end; ++i ) {
            ke_local += constant::half * ForceField::mass[container.type[i]] * normsq( container.v[i] );
        }
        #pragma omp atomic
        ke += ke_local;
        #pragma omp atomic
        n += end - beg;
    }

    RTParameter & parameter;
    double ke = 0;
    int n = 0;
};

struct verlet_nh_initial {
    static std::string name() { return "verlet_nh_initial"; }

    verlet_nh_initial( RTParameter & p ) : parameter( p ) {}

    template<class CONTAINER> void operator () ( CONTAINER & container, int beg, int end ) const {
        const real inertia  = constant::one; // 2/5 * m * R^2
        const real dt       = parameter.dt;

        for ( int i = beg; i < end; ++i ) {
            // v = (v + 0.5 * f / m * dt) / (1 + 0.5 * dt * zeta)
            container.v[i] = ( container.v[i] + constant::half * container.f[i] / ForceField::mass[container.type[i]] * dt ) / ( constant::one + constant::half * dt * parameter.zeta );
            // x = x + v * dt
            container.x[i] += container.v[i] * dt;

            // omega = omega + 0.5 * T / I * dt
            container.o[i] += constant::half * container.t[i] / inertia * dt;

            // mu = mu + (omega X mu) * dt
            auto g = container.n[i] + cross( container.o[i], container.n[i] ) * dt;
            auto scale = rsqrt( normsq( g ) );
            container.n[i] = g * scale;
        }
    }

    RTParameter & parameter;
};

struct bounce_back {
    static std::string name() { return "bounce_back"; }

    bounce_back( RTParameter & p ) : parameter( p ) {}

    template<class CONTAINER> void operator () ( CONTAINER & container, int beg, int end ) const {
        for ( int i = beg; i < end; ++i ) {
            for ( int d = 0; d < 3; ++d ) {
                if ( container.x[i][d] < parameter.box[d][0] ) {
                    container.x[i][d] = parameter.box[d][0] + ( parameter.box[d][0] - container.x[i][d] );
                    container.v[i][d] = -container.v[i][d];
                } else if ( container.x[i][d] > parameter.box[d][1] ) {
                    container.x[i][d] = parameter.box[d][1] - ( container.x[i][d] - parameter.box[d][1] );
                    container.v[i][d] = -container.v[i][d];
                }
            }
        }
    }

    RTParameter & parameter;
};

struct post_torque {
    static std::string name() { return "post_torque"; }

    template<class CONTAINER> void operator () ( CONTAINER & container, int beg, int end ) const {
        for ( int i = beg; i < end; ++i ) {
            container.t[i] = cross( container.n[i], container.t[i] );
        }
    }
};

struct verlet_nh_final {
    static std::string name() { return "verlet_nh_final"; }

    verlet_nh_final( RTParameter & p ) : parameter( p ) {}

    template<class CONTAINER> void operator () ( CONTAINER & container, int beg, int end ) const {
        const size_t number_atoms = container.size();
        const real dt       = parameter.dt;
        const real inertia  = constant::one; // 2/5 * m * R^2

        for ( int i = beg; i < end; ++i ) {
            // v = v + 0.5 * dt * (f / m - zeta * v)
            container.v[i] += constant::half * ( container.f[i] / ForceField::mass[container.type[i]] - parameter.zeta * container.v[i] ) * dt;

            // omega = omega + 0.5 * T / I * dt
            container.o[i] += constant::half * container.t[i] / inertia * dt;
        }
    }

    RTParameter & parameter;
};

struct verlet_initial_bounce_clearforce_update {
    static std::string name() { return "verlet_initial_bounce_clearforce"; }

    verlet_initial_bounce_clearforce_update( RTParameter & p ) : parameter( p ) {}
    ~verlet_initial_bounce_clearforce_update() {
        if ( !parameter.Q ) parameter.set_Q( n );
        parameter.zeta += 0.5 * parameter.dt / parameter.Q * ( ke - 0.5 * 3.0 * n * parameter.kBT );
    }

    template<class CONTAINER> void operator () ( CONTAINER & container, int beg, int end ) {
        const real inertia = constant::one; // 2/5 * m * R^2
        const real dt      = parameter.dt;
        const real gamma   = constant::one / ( constant::one + constant::half * dt * parameter.zeta );
        double ke_local    = 0;

        for ( int i = beg; i < end; ++i ) {
            // Position update
            // v = (v + 0.5 * f / m * dt) / (1 + 0.5 * dt * zeta)
            container.v[i] = ( container.v[i] + real( 0.5 ) / ForceField::mass[container.type[i]] * dt * container.f[i] ) * gamma;
            // x = x + v * dt
            container.x[i] += container.v[i] * dt;

            // Bounce back
            for ( int d = 0; d < 3; ++d ) {
                if ( container.x[i][d] < parameter.box[d][0] ) {
                    container.x[i][d] = parameter.box[d][0] + ( parameter.box[d][0] - container.x[i][d] );
                    container.v[i][d] = -container.v[i][d];
                } else if ( container.x[i][d] > parameter.box[d][1] ) {
                    container.x[i][d] = parameter.box[d][1] - ( container.x[i][d] - parameter.box[d][1] );
                    container.v[i][d] = -container.v[i][d];
                }
            }

            // compute temperature
            ke_local += constant::half * ForceField::mass[container.type[i]] * normsq( container.v[i] );

            // Rotation update
            // omega = omega + 0.5 * T / I * dt
            container.o[i] += constant::half / inertia * dt * container.t[i];

            // mu = mu + (omega X mu) * dt
            container.n[i] = normalize( container.n[i] + cross( container.o[i], container.n[i] ) * dt );

            // Clear force
            container.f[i] = constant::zero;
            container.t[i] = constant::zero;
        }

        #pragma omp atomic
        ke += ke_local;
        #pragma omp atomic
        n += end - beg;
    }

    RTParameter & parameter;
    double ke = 0;
    int n = 0;
};

struct post_toque_final_update {
    static std::string name() { return "post_toque_final_update"; }

    post_toque_final_update( RTParameter & p ) : parameter( p ) {}
    ~post_toque_final_update() {
        if ( !parameter.Q ) parameter.set_Q( n );
        parameter.zeta += 0.5 * parameter.dt / parameter.Q * ( ke - 0.5 * 3.0 * n * parameter.kBT );
    }

    template<class CONTAINER> void operator () ( CONTAINER & container, int beg, int end ) {
        const real dt        = parameter.dt;
        const real inertia   = constant::one; // 2/5 * m * R^2
        const real dt_over_m = dt / inertia;
        double ke_local      = 0;

        for ( int i = beg; i < end; ++i ) {
            container.t[i] = cross( container.n[i], container.t[i] );

            // v = v + 0.5 * dt * (f / m - zeta * v)
            container.v[i] += constant::half * dt * ( container.f[i] / ForceField::mass[container.type[i]] - parameter.zeta * container.v[i] );

            // omega = omega + 0.5 * T / I * dt
            container.o[i] += constant::half * dt_over_m * container.t[i];

            ke_local += constant::half * ForceField::mass[container.type[i]] * normsq( container.v[i] );
        }

        #pragma omp atomic
        ke += ke_local;
        #pragma omp atomic
        n += end - beg;
    }

    RTParameter & parameter;
    double ke = 0;
    int n = 0;
};

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

#endif /* INTEGRATE_NH_H_ */

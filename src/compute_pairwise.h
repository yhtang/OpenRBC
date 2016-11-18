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
#ifndef OPENRBC_PAIRWISE_LL_H
#define OPENRBC_PAIRWISE_LL_H

#include <cmath>
#include <iostream>
#include <omp.h>

#include "aligned_array.h"
#include "container.h"
#include "math_vector.h"
#include "forcefield.h"
#include "pairwise_kernel.h"
#include "voronoi.h"

namespace openrbc {

// lipid-lipid interactions
// particles must be ordered.

#if 1 // atomics-free

void compute_pairwise( LipidContainer & c, const VCellList & lst, const VoronoiDiagram & voronoi ) {
    Service<Timers>::call()["compute_pairwise_omp"].start();

    vect const * const __restrict x  = c.x.data();
    vect const * const __restrict mu = c.n.data();

    static std::vector<AlignedArray<int, true> > stencils( omp_get_max_threads() );

    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        const int nthreads = omp_get_num_threads();
        const int chunk_size = lst.n_cells / nthreads;
        const int cell_begin = chunk_size * tid;
        const int cell_end = ( tid == nthreads - 1 ) ? lst.n_cells : ( cell_begin + chunk_size );
        auto & stencil = stencils[ omp_get_thread_num() ];

        for ( int cell1 = cell_begin; cell1 < cell_end; ++cell1 ) {

            // self interaction

            for ( int p1 = lst.cell_start[cell1]; p1 < lst.cell_start[cell1 + 1]; ++p1 ) {
                const vect x1  = x [p1];
                const vect mu1 = mu[p1];

                for ( int p2 = p1 + 1; p2 < lst.cell_start[cell1 + 1]; ++p2 ) {
                    const vect x2  = x [p2];
                    const vect mu2 = mu[p2];

                    const vect dx = x1 - x2; // 3 FLOPS
                    const real r_sq = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]; // 5 FLOPS

                    if ( r_sq < ForceField::cutsqll && r_sq > 1e-5 )
                        lipid_lipid_omp( dx, r_sq, mu1, mu2, c, p1, p2, newton_on );
                }
            }

            // interaction between cells

            const int n_neighbor = voronoi.get_stencil_whole( cell1, stencil, 6 );

            for ( int i = 0; i < n_neighbor; ++i ) {
                const int cell2 = stencil[i];

                bool same_chunk;
                if ( cell2 >= cell_begin && cell2 < cell_end ) {
                    if ( cell2 <= cell1 ) continue;

                    for ( int p1 = lst.cell_start[cell1]; p1 < lst.cell_start[cell1 + 1]; ++p1 ) {
                        const vect x1  = x [p1];
                        const vect mu1 = mu[p1];

                        for ( int p2 = lst.cell_start[cell2]; p2 < lst.cell_start[cell2 + 1]; ++p2 ) {
                            const vect x2  = x [p2];
                            const vect mu2 = mu[p2];

                            const vect dx = x1 - x2; // 3 FLOPS
                            const real r_sq = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]; // 5 FLOPS

                            if ( r_sq < ForceField::cutsqll && r_sq > 1e-5 )
                                lipid_lipid_omp( dx, r_sq, mu1, mu2, c, p1, p2, newton_on );
                        }
                    }
                } else {
                    for ( int p1 = lst.cell_start[cell1]; p1 < lst.cell_start[cell1 + 1]; ++p1 ) {
                        const vect x1  = x [p1];
                        const vect mu1 = mu[p1];

                        for ( int p2 = lst.cell_start[cell2]; p2 < lst.cell_start[cell2 + 1]; ++p2 ) {
                            const vect x2  = x [p2];
                            const vect mu2 = mu[p2];

                            const vect dx = x1 - x2; // 3 FLOPS
                            const real r_sq = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]; // 5 FLOPS

                            if ( r_sq < ForceField::cutsqll && r_sq > 1e-5 )
                                lipid_lipid_omp( dx, r_sq, mu1, mu2, c, p1, p2, newton_off );
                        }
                    }
                }
            }
        }
    }

    Service<Timers>::call()["compute_pairwise_omp"].stop();
}

#else // atomics

inline void lipid_lipid_omp_atomics(
    const vect dx, const real r_sq,
    const vect mu1, const vect mu2,
    Container & c, RTParameter const & param,
    const int p1, const int p2 ) {
    // 98 FLOPS/pair
    const real r = std::sqrt( r_sq );
    const vect unitdel = dx * ( real( 1 ) / r );

    const real nidotnj = dot( mu1, mu2 );
    const real nidotunitdel = dot( mu1, unitdel );
    const real njdotunitdel = dot( mu2, unitdel );
    const real a = nidotnj - nidotunitdel * njdotunitdel;
    const real A = real( 1.0 ) + param.alphall * ( a - real( 1.0 ) );

    const real rcminusr = param.cutll - r;
    const real rcminusr3 = rcminusr * rcminusr * rcminusr;
    const real rcminusr4 = rcminusr * rcminusr3;
    const real rcminusr7 = rcminusr3 * rcminusr4;

    const vect pni = mu1 - nidotunitdel * unitdel;
    const vect pnj = mu2 - njdotunitdel * unitdel;

    const real ua = param.attll * rcminusr4;
    const real alphaua = param.alphall * ua;
    const real alphauar = alphaua / r;

    const real fra = real( 8.0 ) * param.repll * rcminusr7 + A * real( 4.0 ) * param.attll * rcminusr3;
    const vect f = fra * unitdel + alphauar * ( njdotunitdel * pni + nidotunitdel * pnj );

    for ( int d = 0; d < 3; ++d ) {
        #pragma omp atomic update
        c.f[p1][d] += f[d];             // 3 FLOPS, force
        #pragma omp atomic update
        c.t[p1][d] -= ( alphaua * pnj )[d]; // 6 FLOPS, t
    }

    for ( int d = 0; d < 3; ++d ) {
        #pragma omp atomic update
        c.f[p2][d] -= f[d];             // 3 FLOPS, force
        #pragma omp atomic update
        c.t[p2][d] -= ( alphaua * pni )[d]; // 6 FLOPS, t
    }
}

void compute_pairwise( Container & c, RTParameter const & param, const Voronoi<T4> & lst ) {
    Service<Timers>::call()["compute_pairwise_omp"].start();

    vect const * const __restrict x  = c.x.data();
    vect const * const __restrict mu = c.n.data();

    static std::vector<AlignedArray<int, true> > stencils( omp_get_max_threads() );

    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        const int nthreads = omp_get_num_threads();
        const int chunk_size = lst.n_cells / nthreads;
        const int cell_begin = chunk_size * tid;
        const int cell_end = ( tid == nthreads - 1 ) ? lst.n_cells : ( cell_begin + chunk_size );
        auto & stencil = stencils[ omp_get_thread_num() ];

        for ( int b = cell_begin; b < cell_end; ++b ) {

            // self interaction

            for ( int p1 = lst.cell_start[b]; p1 < lst.cell_start[b + 1]; ++p1 ) {
                const vect x1  = x [p1];
                const vect mu1 = mu[p1];

                for ( int p2 = p1 + 1; p2 < lst.cell_start[b + 1]; ++p2 ) {
                    const vect x2  = x [p2];
                    const vect mu2 = mu[p2];

                    const vect dx = x1 - x2; // 3 FLOPS
                    const real r_sq = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]; // 5 FLOPS

                    if ( r_sq < ForceField::cutsqll && r_sq > 1e-5 ) lipid_lipid_omp_atomics( dx, r_sq, mu1, mu2, c, param, p1, p2 );
                }
            }

            // interaction between cells

            const int ncell = lst.get_stencil( b, stencil );

            for ( int i = 0; i < ncell; ++i ) {
                const int cell2 = stencil[i];

                for ( int p1 = lst.cell_start[b]; p1 < lst.cell_start[b + 1]; ++p1 ) {
                    const vect x1  = x [p1];
                    const vect mu1 = mu[p1];

                    for ( int p2 = lst.cell_start[cell2]; p2 < lst.cell_start[cell2 + 1]; ++p2 ) {
                        const vect x2  = x [p2];
                        const vect mu2 = mu[p2];

                        const vect dx = x1 - x2; // 3 FLOPS
                        const real r_sq = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]; // 5 FLOPS

                        if ( r_sq < ForceField::cutsqll && r_sq > 1e-5 )
                            lipid_lipid_omp_atomics( dx, r_sq, mu1, mu2, c, param, p1, p2 );
                    }
                }
            }
        }
    }

    Service<Timers>::call()["compute_pairwise_omp"].stop();
}

#endif

}

#endif

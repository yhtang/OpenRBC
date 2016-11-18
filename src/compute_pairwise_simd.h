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
#ifndef OPENRBC_PAIRWISE_SIMD_H
#define OPENRBC_PAIRWISE_SIMD_H

#include <cmath>
#include <iostream>
#include <omp.h>

#include "aligned_array.h"
#include "container.h"
#include "math_vector.h"
#include "forcefield.h"
#include "pairwise_kernel_simd.h"
#include "voronoi.h"

namespace openrbc {

#ifdef EXPLICIT_SIMD

// pairwise interactions
// particles must be ordered.
struct lipid_lipid_simd {
    constexpr static config::real rmax = 6.0;

    static inline void intra_cell( LipidContainer & lipid, const int p1_beg, const int p1_end ) {
        vect const * const __restrict x = lipid.x.data();
        vect const * const __restrict n = lipid.n.data();

        for ( int i = p1_beg; i < p1_end; ++i ) {
            const vect x1 = x[i];
            const vect n1 = n[i];
            vect f1 {0.f}, t1 {0.f};

            for ( int j = i + 1; j < p1_end; ++j ) {
                const vect x2 = x[j];

                const vect dx = x1 - x2;
                const vect r_sq = _normsq( dx );

                if ( r_sq[0] < ForceField::cutsqll && r_sq[0] > 1e-5 )
                    lipid_lipid_omp( dx, r_sq, n1, n[j], lipid, i, j, f1, t1, newton_on );
            }

            lipid.f[i] += f1;
            lipid.t[i] += t1;
        }
    }

    template<class NEWTON>
    static inline void inter_cell( LipidContainer & lipid,
                                   const int p1_beg, const int p1_end,
                                   const int p2_beg, const int p2_end,
                                   NEWTON const newton ) {
        vect const * const __restrict x = lipid.x.data();
        vect const * const __restrict n = lipid.n.data();

        for ( int i = p1_beg; i < p1_end; ++i ) {
            const vect x1 = x[i];
            const vect n1 = n[i];
            vect f1 {0.f}, t1 {0.f};

            for ( int j = p2_beg; j < p2_end; ++j ) {
                const vect x2 = x[j];

                const vect dx = x1 - x2;
                const vect r_sq = _normsq( dx );

                if ( r_sq[0] < ForceField::cutsqll && r_sq[0] > 1e-5 )
                    lipid_lipid_omp( dx, r_sq, n1, n[j], lipid, i, j, f1, t1, newton );
            }

            lipid.f[i] += f1;
            lipid.t[i] += t1;
        }
    }
};

struct prote_lipid_simd {
    constexpr static config::real rmax = 8.0;

    template<typename ...T>
    static inline void intra_cell( T const && ...whatever ) {}

    template<class NEWTON1, class NEWTON2>
    static inline void inter_cell( ProteContainer & protein,
                                   LipidContainer & lipid,
                                   const int p1_beg, const int p1_end,
                                   const int p2_beg, const int p2_end,
                                   NEWTON1 const newton1, NEWTON2 const newton2 ) {
        vect const * const __restrict px = protein.x.data();
        vect const * const __restrict pn = protein.n.data();
        vect const * const __restrict lx = lipid.x.data();
        vect const * const __restrict ln = lipid.n.data();

        for ( int i = p1_beg; i < p1_end; ++i ) {
            const vect x1 = px[i];
            const vect n1 = pn[i];
            const int type = protein.type[i];

            for ( int j = p2_beg; j < p2_end; ++j ) {
                const vect x2 = lx[j];
                const vect n2 = ln[j];

                const vect dx = x1 - x2;
                const vect r_sq = _normsq( dx );

                if ( r_sq[0] < ForceField::cutsqlp[type] && r_sq[0] > 1e-5 )
                    protein_lipid_omp( dx, r_sq, n1, n2, protein, lipid, i, j, type, newton1, newton2 );
                else if ( r_sq[0] < ForceField::lj_cutsq[type] && r_sq[0] > 1e-5 )
                    lennard_jones_omp( dx, r_sq, protein, lipid, i, j, type, newton1, newton2 );
            }
        }
    }
};

struct prote_prote_simd {
    constexpr static config::real rmax = 9.0;

    static inline void intra_cell( ProteContainer & protein, const int p1_beg, const int p1_end ) {

        vect const * const __restrict x    = protein.x .data();
        int   const * const __restrict type = protein.type.data();

        for ( int i = p1_beg; i < p1_end; ++i ) {
            const vect x1    = x   [i];
            const int   type1 = type[i];

            for ( int j = i + 1; j < p1_end; ++j ) {
                const vect x2    = x[j];
                const int   type2 = type[j];

                const vect dx = x1 - x2;
                const vect r_sq = _normsq( dx );

                const int type12 = type1 + type2 * ForceField::n_type;
                if ( r_sq[0] < ForceField::cutsqpp[type12] && r_sq[0] > 1e-5 )
                    protein_protein_omp( dx, r_sq, protein, i, j, type12, newton_on );
                else if ( r_sq[0] < ForceField::lj_cutsq[type12] && r_sq[0] > 1e-5 )
                    lennard_jones_omp( dx, r_sq, protein, protein, i, j, type12, newton_on, newton_on );
            }
        }
    }

    template<class NEWTON>
    static inline void inter_cell( ProteContainer & protein,
                                   const int p1_beg, const int p1_end,
                                   const int p2_beg, const int p2_end,
                                   NEWTON const newton ) {

        vect const * const __restrict x    = protein.x   .data();
        int   const * const __restrict type = protein.type.data();

        for ( int i = p1_beg; i < p1_end; ++i ) {
            const vect x1    = x   [i];
            const int   type1 = type[i];

            for ( int j = p2_beg; j < p2_end; ++j ) {
                const vect x2    = x[j];
                const int   type2 = type[j];

                const vect dx = x1 - x2;
                const vect r_sq = _normsq( dx );

                const int type12 = type1 + type2 * ForceField::n_type;
                if ( r_sq[0] < ForceField::cutsqpp[type12] && r_sq[0] > 1e-5 )
                    protein_protein_omp( dx, r_sq, protein, i, j, type12, newton );
                else if ( r_sq[0] < ForceField::lj_cutsq[type12] && r_sq[0] > 1e-5 )
                    lennard_jones_omp( dx, r_sq, protein, protein, i, j, type12, newton_on, newton );
            }
        }
    }
};

void compute_pairwise_simd( VoronoiDiagram const & voronoi,
                            LipidContainer & lipid,
                            ProteContainer & protein,
                            VCellList const & cell_lipid,
                            VCellList const & cell_protein ) {
    Service<Timers>::call()["compute_pairwise_simd"].start();

    static std::vector<AlignedArray<int, true> > stencils( omp_get_max_threads() );

    #pragma omp parallel
    {
        int beg = voronoi.load_balancer().beg();
        int end = voronoi.load_balancer().end();
        auto & stencil = stencils[ omp_get_thread_num() ];

        for ( int cell1 = beg; cell1 < end; ++cell1 ) {

            // self interaction
            lipid_lipid_simd::intra_cell( lipid,   cell_lipid  .cell_start[cell1], cell_lipid  .cell_start[cell1 + 1] );
            prote_prote_simd::intra_cell( protein, cell_protein.cell_start[cell1], cell_protein.cell_start[cell1 + 1] );

            // interaction between cells
            int n_neighbor = voronoi.get_stencil_whole( cell1, stencil, prote_prote_simd::rmax );

            for ( int i = 0; i < n_neighbor; ++i ) {
                const int cell2 = stencil[i];
                if ( cell2 >= beg && cell2 < end ) {
                    if ( cell2 > cell1 )
                        prote_prote_simd::inter_cell( protein,
                                                      cell_protein.cell_start[cell1], cell_protein.cell_start[cell1 + 1],
                                                      cell_protein.cell_start[cell2], cell_protein.cell_start[cell2 + 1],
                                                      newton_on );
                } else {
                    prote_prote_simd::inter_cell( protein,
                                                  cell_protein.cell_start[cell1], cell_protein.cell_start[cell1 + 1],
                                                  cell_protein.cell_start[cell2], cell_protein.cell_start[cell2 + 1],
                                                  newton_off );
                }
            }

            n_neighbor = voronoi.refine_stencil( cell1, n_neighbor, stencil, prote_lipid_simd::rmax );

            for ( int i = 0; i < n_neighbor; ++i ) {
                const int cell2 = stencil[i];
                if ( cell2 >= beg && cell2 < end ) {
                    prote_lipid_simd::inter_cell( protein, lipid,
                                                  cell_protein.cell_start[cell1], cell_protein.cell_start[cell1 + 1],
                                                  cell_lipid  .cell_start[cell2], cell_lipid  .cell_start[cell2 + 1],
                                                  newton_on, newton_on );
                } else {
                    prote_lipid_simd::inter_cell( protein, lipid,
                                                  cell_protein.cell_start[cell1], cell_protein.cell_start[cell1 + 1],
                                                  cell_lipid  .cell_start[cell2], cell_lipid  .cell_start[cell2 + 1],
                                                  newton_on, newton_off );
                    prote_lipid_simd::inter_cell( protein, lipid,
                                                  cell_protein.cell_start[cell2], cell_protein.cell_start[cell2 + 1],
                                                  cell_lipid  .cell_start[cell1], cell_lipid  .cell_start[cell1 + 1],
                                                  newton_off, newton_on );
                }
            }

            n_neighbor = voronoi.refine_stencil( cell1, n_neighbor, stencil, lipid_lipid_simd::rmax );

            for ( int i = 0; i < n_neighbor; ++i ) {
                const int cell2 = stencil[i];
                if ( cell2 >= beg && cell2 < end ) {
                    if ( cell2 > cell1 )
                        lipid_lipid_simd::inter_cell( lipid,
                                                      cell_lipid.cell_start[cell1], cell_lipid.cell_start[cell1 + 1],
                                                      cell_lipid.cell_start[cell2], cell_lipid.cell_start[cell2 + 1],
                                                      newton_on );
                } else {
                    lipid_lipid_simd::inter_cell( lipid,
                                                  cell_lipid.cell_start[cell1], cell_lipid.cell_start[cell1 + 1],
                                                  cell_lipid.cell_start[cell2], cell_lipid.cell_start[cell2 + 1],
                                                  newton_off );
                }
            }
        }
    }

    Service<Timers>::call()["compute_pairwise_simd"].stop();
}

#endif

}

#endif

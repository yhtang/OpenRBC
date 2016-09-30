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
#ifndef PAIRWISE_FUSED_H
#define PAIRWISE_FUSED_H

#include <cmath>
#include <iostream>
#include <omp.h>

#include "aligned_array.h"
#include "container.h"
#include "math_vector.h"
#include "forcefield.h"
#include "pairwise_kernel_fused.h"
#include "voronoi.h"

namespace openrbc {

// pairwise interactions
// particles must be ordered.

#ifdef IMPLICIT_SIMD

struct lipid_lipid {
    constexpr static config::real rmax = 6.0;
    constexpr static int simd_width = 8;

    static inline void intra_cell( LipidContainer & c, const int p1_beg, const int p1_end ) {
        vect const * const __restrict x = c.x.data();
        vect const * const __restrict n = c.n.data();

        int queue = 0;
        int p1[simd_width], p2[simd_width];
        real drx[simd_width], dry[simd_width], drz[simd_width];
        real n1x[simd_width], n1y[simd_width], n1z[simd_width];
        real n2x[simd_width], n2y[simd_width], n2z[simd_width];
        real rsq[simd_width];

        for ( int i = p1_beg; i < p1_end; ++i ) {
            const vect x1 = x[i];
            const vect n1 = n[i];

            for ( int j = i + 1; j < p1_end; ++j ) {
                const vect x2  = x[j];

                const vect dx = x1 - x2;
                const real r_sq = normsq( dx );

                if ( r_sq < ForceField::cutsqll && r_sq > 1e-5f ) {
                    p1[queue] = i;
                    p2[queue] = j;
                    drx[queue] = dx[0],   dry[queue] = dx[1],   drz[queue] = dx[2];
                    n1x[queue] = n1[0],   n1y[queue] = n1[1],   n1z[queue] = n1[2];
                    n2x[queue] = n[j][0], n2y[queue] = n[j][1], n2z[queue] = n[j][2];
                    rsq[queue] = r_sq;
                    if ( ++queue == simd_width ) {
                        real fx[simd_width], fy[simd_width], fz[simd_width];
                        real t1x[simd_width], t1y[simd_width], t1z[simd_width];
                        real t2x[simd_width], t2y[simd_width], t2z[simd_width];
                        for ( int k = 0; k < simd_width; ++k ) {
                            real r         = std::sqrt( rsq[k] );
                            real ex        = drx[k] / r;
                            real ey        = dry[k] / r;
                            real ez        = drz[k] / r;
                            real nidotnj   = n1x[k] * n2x[k] + n1y[k] * n2y[k] + n1z[k] * n2z[k];
                            real nidote    = n1x[k] *  ex + n1y[k] *  ey + n1z[k] *  ez;
                            real njdote    = n2x[k] *  ex + n2y[k] *  ey + n2z[k] *  ez;
                            real a         = nidotnj - nidote * njdote;
                            real A         = 1.f + ForceField::alphall * ( a - 1.f );
                            real rcminusr  = ForceField::cutll - r;
                            real rcminusr3 = rcminusr * rcminusr * rcminusr;
                            real rcminusr4 = rcminusr * rcminusr3;
                            real rcminusr7 = rcminusr3 * rcminusr4;
                            real pnix      = n1x[k] - nidote * ex;
                            real pniy      = n1y[k] - nidote * ey;
                            real pniz      = n1z[k] - nidote * ez;
                            real pnjx      = n2x[k] - njdote * ex;
                            real pnjy      = n2y[k] - njdote * ey;
                            real pnjz      = n2z[k] - njdote * ez;
                            real ua        = ForceField::attll * rcminusr4;
                            real alphaua   = ForceField::alphall * ua;
                            real alphauar  = alphaua / r;
                            real fra       = 8.0f * ForceField::repll * rcminusr7 + A * 4.f * ForceField::attll * rcminusr3;
                            fx [k] = fra * ex + alphauar * ( njdote * pnix + nidote * pnjx );
                            fy [k] = fra * ey + alphauar * ( njdote * pniy + nidote * pnjy );
                            fz [k] = fra * ez + alphauar * ( njdote * pniz + nidote * pnjz );
                            t1x[k] = alphaua * pnjx;
                            t1y[k] = alphaua * pnjy;
                            t1z[k] = alphaua * pnjz;
                            t2x[k] = alphaua * pnix;
                            t2y[k] = alphaua * pniy;
                            t2z[k] = alphaua * pniz;
                        }
                        for ( int k = 0; k < simd_width; ++k ) {
                            c.f[ p1[k] ] += vect( fx [k], fy [k], fz [k] );
                            c.t[ p1[k] ] -= vect( t1x[k], t1y[k], t1z[k] );
                            c.f[ p2[k] ] -= vect( fx [k], fy [k], fz [k] );
                            c.t[ p2[k] ] -= vect( t2x[k], t2y[k], t2z[k] );
                        }
                        queue = 0;
                    }
                }
            }
        }
        if ( queue ) {
            real fx[simd_width], fy[simd_width], fz[simd_width];
            real t1x[simd_width], t1y[simd_width], t1z[simd_width];
            real t2x[simd_width], t2y[simd_width], t2z[simd_width];
            for ( int k = 0; k < queue; ++k ) {
                real r         = std::sqrt( rsq[k] );
                real ex        = drx[k] / r;
                real ey        = dry[k] / r;
                real ez        = drz[k] / r;
                real nidotnj   = n1x[k] * n2x[k] + n1y[k] * n2y[k] + n1z[k] * n2z[k];
                real nidote    = n1x[k] *  ex + n1y[k] *  ey + n1z[k] *  ez;
                real njdote    = n2x[k] *  ex + n2y[k] *  ey + n2z[k] *  ez;
                real a         = nidotnj - nidote * njdote;
                real A         = 1.f + ForceField::alphall * ( a - 1.f );
                real rcminusr  = ForceField::cutll - r;
                real rcminusr3 = rcminusr * rcminusr * rcminusr;
                real rcminusr4 = rcminusr * rcminusr3;
                real rcminusr7 = rcminusr3 * rcminusr4;
                real pnix      = n1x[k] - nidote * ex;
                real pniy      = n1y[k] - nidote * ey;
                real pniz      = n1z[k] - nidote * ez;
                real pnjx      = n2x[k] - njdote * ex;
                real pnjy      = n2y[k] - njdote * ey;
                real pnjz      = n2z[k] - njdote * ez;
                real ua        = ForceField::attll * rcminusr4;
                real alphaua   = ForceField::alphall * ua;
                real alphauar  = alphaua / r;
                real fra       = 8.0f * ForceField::repll * rcminusr7 + A * 4.f * ForceField::attll * rcminusr3;
                fx [k] = fra * ex + alphauar * ( njdote * pnix + nidote * pnjx );
                fy [k] = fra * ey + alphauar * ( njdote * pniy + nidote * pnjy );
                fz [k] = fra * ez + alphauar * ( njdote * pniz + nidote * pnjz );
                t1x[k] = alphaua * pnjx;
                t1y[k] = alphaua * pnjy;
                t1z[k] = alphaua * pnjz;
                t2x[k] = alphaua * pnix;
                t2y[k] = alphaua * pniy;
                t2z[k] = alphaua * pniz;
            }
            for ( int k = 0; k < queue; ++k ) {
                c.f[ p1[k] ] += vect( fx [k], fy [k], fz [k] );
                c.t[ p1[k] ] -= vect( t1x[k], t1y[k], t1z[k] );
                c.f[ p2[k] ] -= vect( fx [k], fy [k], fz [k] );
                c.t[ p2[k] ] -= vect( t2x[k], t2y[k], t2z[k] );
            }
            queue = 0;
        }
    }

    template<bool NEWTON>
    static inline void inter_cell( LipidContainer & c,
                                   const int p1_beg, const int p1_end,
                                   const int p2_beg, const int p2_end,
                                   boolean<NEWTON> const newton ) {
        vect const * const __restrict x = c.x.data();
        vect const * const __restrict n = c.n.data();

        int queue = 0;
        int p1[simd_width], p2[simd_width];
        real drx[simd_width], dry[simd_width], drz[simd_width];
        real n1x[simd_width], n1y[simd_width], n1z[simd_width];
        real n2x[simd_width], n2y[simd_width], n2z[simd_width];
        real rsq[simd_width];

        for ( int i = p1_beg; i < p1_end; ++i ) {
            const vect x1 = x[i];
            const vect n1 = n[i];

            for ( int j = p2_beg; j < p2_end; ++j ) {
                const vect x2  = x[j];

                const vect dx = x1 - x2;
                const real r_sq = normsq( dx );

                if ( r_sq < ForceField::cutsqll && r_sq > 1e-5f ) {
                    p1[queue] = i;
                    p2[queue] = j;
                    drx[queue] = dx[0],   dry[queue] = dx[1],   drz[queue] = dx[2];
                    n1x[queue] = n1[0],   n1y[queue] = n1[1],   n1z[queue] = n1[2];
                    n2x[queue] = n[j][0], n2y[queue] = n[j][1], n2z[queue] = n[j][2];
                    rsq[queue] = r_sq;
                    if ( ++queue == simd_width ) {
                        real fx[simd_width], fy[simd_width], fz[simd_width];
                        real t1x[simd_width], t1y[simd_width], t1z[simd_width];
                        real t2x[simd_width], t2y[simd_width], t2z[simd_width];
                        for ( int k = 0; k < simd_width; ++k ) {
                            real r         = std::sqrt( rsq[k] );
                            real ex        = drx[k] / r;
                            real ey        = dry[k] / r;
                            real ez        = drz[k] / r;
                            real nidotnj   = n1x[k] * n2x[k] + n1y[k] * n2y[k] + n1z[k] * n2z[k];
                            real nidote    = n1x[k] *  ex + n1y[k] *  ey + n1z[k] *  ez;
                            real njdote    = n2x[k] *  ex + n2y[k] *  ey + n2z[k] *  ez;
                            real a         = nidotnj - nidote * njdote;
                            real A         = 1.f + ForceField::alphall * ( a - 1.f );
                            real rcminusr  = ForceField::cutll - r;
                            real rcminusr3 = rcminusr * rcminusr * rcminusr;
                            real rcminusr4 = rcminusr * rcminusr3;
                            real rcminusr7 = rcminusr3 * rcminusr4;
                            real pnix      = n1x[k] - nidote * ex;
                            real pniy      = n1y[k] - nidote * ey;
                            real pniz      = n1z[k] - nidote * ez;
                            real pnjx      = n2x[k] - njdote * ex;
                            real pnjy      = n2y[k] - njdote * ey;
                            real pnjz      = n2z[k] - njdote * ez;
                            real ua        = ForceField::attll * rcminusr4;
                            real alphaua   = ForceField::alphall * ua;
                            real alphauar  = alphaua / r;
                            real fra       = 8.0f * ForceField::repll * rcminusr7 + A * 4.f * ForceField::attll * rcminusr3;
                            fx [k] = fra * ex + alphauar * ( njdote * pnix + nidote * pnjx );
                            fy [k] = fra * ey + alphauar * ( njdote * pniy + nidote * pnjy );
                            fz [k] = fra * ez + alphauar * ( njdote * pniz + nidote * pnjz );
                            t1x[k] = alphaua * pnjx;
                            t1y[k] = alphaua * pnjy;
                            t1z[k] = alphaua * pnjz;
                            if ( NEWTON ) t2x[k] = alphaua * pnix;
                            if ( NEWTON ) t2y[k] = alphaua * pniy;
                            if ( NEWTON ) t2z[k] = alphaua * pniz;
                        }
                        for ( int k = 0; k < simd_width; ++k ) {
                            c.f[ p1[k] ] += vect( fx [k], fy [k], fz [k] );
                            c.t[ p1[k] ] -= vect( t1x[k], t1y[k], t1z[k] );
                            if ( NEWTON ) c.f[ p2[k] ] -= vect( fx [k], fy [k], fz [k] );
                            if ( NEWTON ) c.t[ p2[k] ] -= vect( t2x[k], t2y[k], t2z[k] );
                        }
                        queue = 0;
                    }
                }
            }
        }

        if ( queue ) {
            real fx[simd_width], fy[simd_width], fz[simd_width];
            real t1x[simd_width], t1y[simd_width], t1z[simd_width];
            real t2x[simd_width], t2y[simd_width], t2z[simd_width];
            for ( int k = 0; k < queue; ++k ) {
                real r         = std::sqrt( rsq[k] );
                real ex        = drx[k] / r;
                real ey        = dry[k] / r;
                real ez        = drz[k] / r;
                real nidotnj   = n1x[k] * n2x[k] + n1y[k] * n2y[k] + n1z[k] * n2z[k];
                real nidote    = n1x[k] *  ex + n1y[k] *  ey + n1z[k] *  ez;
                real njdote    = n2x[k] *  ex + n2y[k] *  ey + n2z[k] *  ez;
                real a         = nidotnj - nidote * njdote;
                real A         = 1.f + ForceField::alphall * ( a - 1.f );
                real rcminusr  = ForceField::cutll - r;
                real rcminusr3 = rcminusr * rcminusr * rcminusr;
                real rcminusr4 = rcminusr * rcminusr3;
                real rcminusr7 = rcminusr3 * rcminusr4;
                real pnix      = n1x[k] - nidote * ex;
                real pniy      = n1y[k] - nidote * ey;
                real pniz      = n1z[k] - nidote * ez;
                real pnjx      = n2x[k] - njdote * ex;
                real pnjy      = n2y[k] - njdote * ey;
                real pnjz      = n2z[k] - njdote * ez;
                real ua        = ForceField::attll * rcminusr4;
                real alphaua   = ForceField::alphall * ua;
                real alphauar  = alphaua / r;
                real fra       = 8.0f * ForceField::repll * rcminusr7 + A * 4.f * ForceField::attll * rcminusr3;
                fx [k] = fra * ex + alphauar * ( njdote * pnix + nidote * pnjx );
                fy [k] = fra * ey + alphauar * ( njdote * pniy + nidote * pnjy );
                fz [k] = fra * ez + alphauar * ( njdote * pniz + nidote * pnjz );
                t1x[k] = alphaua * pnjx;
                t1y[k] = alphaua * pnjy;
                t1z[k] = alphaua * pnjz;
                if ( NEWTON ) t2x[k] = alphaua * pnix;
                if ( NEWTON ) t2y[k] = alphaua * pniy;
                if ( NEWTON ) t2z[k] = alphaua * pniz;
            }
            for ( int k = 0; k < queue; ++k ) {
                c.f[ p1[k] ] += vect( fx [k], fy [k], fz [k] );
                c.t[ p1[k] ] -= vect( t1x[k], t1y[k], t1z[k] );
                if ( NEWTON ) c.f[ p2[k] ] -= vect( fx [k], fy [k], fz [k] );
                if ( NEWTON ) c.t[ p2[k] ] -= vect( t2x[k], t2y[k], t2z[k] );
            }
        }
    }
};

#else

struct lipid_lipid {
    constexpr static config::real rmax = 6.0;

    static inline void intra_cell( LipidContainer & c, const int p1_beg, const int p1_end ) {
        vect const * const __restrict x = c.x.data();
        vect const * const __restrict n = c.n.data();

        for ( int i = p1_beg; i < p1_end; ++i ) {
            const vect x1  = x[i];
            const vect mu1 = n[i];

            for ( int j = i + 1; j < p1_end; ++j ) {
                const vect x2  = x[j];
                const vect mu2 = n[j];

                const vect dx = x1 - x2;
                const real r_sq = normsq( dx );

                if ( r_sq < ForceField::cutsqll && r_sq > 1e-5 )
                    lipid_lipid_omp( dx, r_sq, mu1, mu2, c, i, j, newton_on );
            }
        }
    }

    template<class NEWTON>
    static inline void inter_cell( LipidContainer & c,
                                   const int p1_beg, const int p1_end,
                                   const int p2_beg, const int p2_end,
                                   NEWTON const newton ) {
        vect const * const __restrict x = c.x.data();
        vect const * const __restrict n = c.n.data();

        for ( int i = p1_beg; i < p1_end; ++i ) {
            const vect x1  = x[i];
            const vect mu1 = n[i];

            for ( int j = p2_beg; j < p2_end; ++j ) {
                const vect x2  = x[j];
                const vect mu2 = n[j];

                const vect dx = x1 - x2;
                const real r_sq = normsq( dx );

                if ( r_sq < ForceField::cutsqll && r_sq > 1e-5 )
                    lipid_lipid_omp( dx, r_sq, mu1, mu2, c, i, j, newton );
            }
        }
    }
};

#endif

struct prote_lipid {
    constexpr static config::real rmax = 9.0;

    template<typename ...T>
    static inline void intra_cell( T const && ...whatever ) {}

    template<class NEWTON1, class NEWTON2>
    static inline void inter_cell( ProteContainer & protein,
                                   LipidContainer & lipid,
                                   const int p1_beg, const int p1_end,
                                   const int p2_beg, const int p2_end,
                                   NEWTON1 const newton1, NEWTON2 const newton2 ) {
        vect const * const __restrict lx = lipid.x.data();
        vect const * const __restrict ln = lipid.n.data();
        vect const * const __restrict px = protein.x.data();
        vect const * const __restrict pn = protein.n.data();

        for ( int i = p1_beg; i < p1_end; ++i ) {
            const vect x1  = px[i];
            const vect mu1 = pn[i];
            const int type = protein.type[i];

            for ( int j = p2_beg; j < p2_end; ++j ) {
                const vect x2  = lx[j];
                const vect mu2 = ln[j];

                const vect dx = x1 - x2;
                const real r_sq = normsq( dx );

                if ( r_sq < ForceField::cutsqlp[type] && r_sq > 1e-5 )
                    protein_lipid_omp( dx, r_sq, mu1, mu2, protein, lipid, i, j, type, newton1, newton2 );
                else if ( r_sq < ForceField::lj_cutsq[type] && r_sq > 1e-5 )
                    lennard_jones_omp( dx, r_sq, protein, lipid, i, j, type, newton1, newton2 );
            }
        }
    }
};

struct prote_prote {
    constexpr static config::real rmax = 9.0;

    static inline void intra_cell( ProteContainer & protein, const int p1_beg, const int p1_end ) {
        vect const * const __restrict x    = protein.x   .data();
        int  const * const __restrict type = protein.type.data();

        for ( int i = p1_beg; i < p1_end; ++i ) {
            const vect x1    = x   [i];
            const int          type1 = type[i];

            for ( int j = i + 1; j < p1_end; ++j ) {
                const vect x2    = x[j];
                const int          type2 = type[j];

                const vect dx = x1 - x2;
                const real r_sq = normsq( dx );

                const int type12 = type1 + type2 * ForceField::n_type;
                if ( r_sq < ForceField::cutsqpp[type12] && r_sq > 1e-5 )
                    protein_protein_omp( dx, r_sq, protein, i, j, type12, newton_on );
                else if ( r_sq < ForceField::lj_cutsq[type12] && r_sq > 1e-5 )
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
        int  const * const __restrict type = protein.type.data();

        for ( int i = p1_beg; i < p1_end; ++i ) {
            const vect x1    = x   [i];
            const int          type1 = type[i];

            for ( int j = p2_beg; j < p2_end; ++j ) {
                const vect x2    = x[j];
                const int          type2 = type[j];

                const vect dx = x1 - x2;
                const real r_sq = normsq( dx );

                const int type12 = type1 + type2 * ForceField::n_type;
                if ( r_sq < ForceField::cutsqpp[type12] && r_sq > 1e-5 )
                    protein_protein_omp( dx, r_sq, protein, i, j, type12, newton );
                else if ( r_sq < ForceField::lj_cutsq[type12] && r_sq > 1e-5 )
                    lennard_jones_omp( dx, r_sq, protein, protein, i, j, type12, newton_on, newton );
            }
        }
    }
};

void compute_pairwise_fused( VoronoiDiagram const & voronoi,
                             LipidContainer & lipid,
                             ProteContainer & protein,
                             VCellList const & cell_lipid,
                             VCellList const & cell_protein ) {
    Service<Timers>::call()["compute_pairwise_fused"].start();

    static std::vector<AlignedArray<int, true> > stencils( omp_get_max_threads() );

    #pragma omp parallel
    {
        int beg = voronoi.load_balancer().beg();
        int end = voronoi.load_balancer().end();
        auto & stencil = stencils[ omp_get_thread_num() ];

        for ( int cell1 = beg; cell1 < end; ++cell1 ) {

            // self interaction
            lipid_lipid::intra_cell( lipid,   cell_lipid  .cell_start[cell1], cell_lipid  .cell_start[cell1 + 1] );
            prote_prote::intra_cell( protein, cell_protein.cell_start[cell1], cell_protein.cell_start[cell1 + 1] );

            // interaction between cells
            int n_neighbor = voronoi.get_stencil_whole( cell1, stencil, prote_prote::rmax );

            for ( int i = 0; i < n_neighbor; ++i ) {
                const int cell2 = stencil[i];
                if ( cell2 >= beg && cell2 < end ) {
                    if ( cell2 > cell1 )
                        prote_prote::inter_cell( protein,
                                                 cell_protein.cell_start[cell1], cell_protein.cell_start[cell1 + 1],
                                                 cell_protein.cell_start[cell2], cell_protein.cell_start[cell2 + 1],
                                                 newton_on );
                    prote_lipid::inter_cell( protein, lipid,
                                             cell_protein.cell_start[cell1], cell_protein.cell_start[cell1 + 1],
                                             cell_lipid  .cell_start[cell2], cell_lipid  .cell_start[cell2 + 1],
                                             newton_on, newton_on );
                } else {
                    prote_prote::inter_cell( protein,
                                             cell_protein.cell_start[cell1], cell_protein.cell_start[cell1 + 1],
                                             cell_protein.cell_start[cell2], cell_protein.cell_start[cell2 + 1],
                                             newton_off );
                    prote_lipid::inter_cell( protein, lipid,
                                             cell_protein.cell_start[cell1], cell_protein.cell_start[cell1 + 1],
                                             cell_lipid  .cell_start[cell2], cell_lipid  .cell_start[cell2 + 1],
                                             newton_on, newton_off );
                    prote_lipid::inter_cell( protein, lipid,
                                             cell_protein.cell_start[cell2], cell_protein.cell_start[cell2 + 1],
                                             cell_lipid  .cell_start[cell1], cell_lipid  .cell_start[cell1 + 1],
                                             newton_off, newton_on );
                }
            }

            n_neighbor = voronoi.refine_stencil( cell1, n_neighbor, stencil, lipid_lipid::rmax );

            for ( int i = 0; i < n_neighbor; ++i ) {
                const int cell2 = stencil[i];
                if ( cell2 >= beg && cell2 < end ) {
                    if ( cell2 > cell1 )
                        lipid_lipid::inter_cell( lipid,
                                                 cell_lipid.cell_start[cell1], cell_lipid.cell_start[cell1 + 1],
                                                 cell_lipid.cell_start[cell2], cell_lipid.cell_start[cell2 + 1],
                                                 newton_on );
                } else {
                    lipid_lipid::inter_cell( lipid,
                                             cell_lipid.cell_start[cell1], cell_lipid.cell_start[cell1 + 1],
                                             cell_lipid.cell_start[cell2], cell_lipid.cell_start[cell2 + 1],
                                             newton_off );
                }
            }
        }
    }

    Service<Timers>::call()["compute_pairwise_fused"].stop();
}

}

#endif /* PAIRWISE_FUSED_H */

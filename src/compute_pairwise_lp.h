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
#ifndef OPENRBC_PAIRWISE_LP_H
#define OPENRBC_PAIRWISE_LP_H

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

// lipid-protein interactions
void compute_pairwise_lp( LipidContainer & lipid, ProteContainer & protein, const VoronoiDiagram & voronoi, const VCellList & cl_lipid, const VCellList & cl_protein ) {
    if ( !protein.size() ) return;

    Service<Timers>::call()["compute_pairwise_lp"].start();

    vect const * const __restrict lx  = lipid.x.data();
    vect const * const __restrict ln = lipid.n.data();
    vect const * const __restrict px  = protein.x.data();
    vect const * const __restrict pn = protein.n.data();

    static std::vector<AlignedArray<int, true> > stencils( omp_get_max_threads() );

    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        const int nthreads = omp_get_num_threads();
        const int chunk_size = voronoi.n_cells / nthreads;
        const int cell_begin = chunk_size * tid;
        const int cell_end = ( tid == nthreads - 1 ) ? voronoi.n_cells : ( cell_begin + chunk_size );
        auto & stencil = stencils[ omp_get_thread_num() ];

        for ( int b = cell_begin; b < cell_end; ++b ) {

            if ( cl_protein.cell_start[b] == cl_protein.cell_start[b + 1] ) continue;

            // interaction between cells

            const int ncell = voronoi.get_stencil_whole( b, stencil, 9 );

            for ( int i = 0; i < ncell; ++i ) {
                const int cell2 = stencil[i];

                for ( int p1 = cl_protein.cell_start[b]; p1 < cl_protein.cell_start[b + 1]; ++p1 ) {
                    const vect x1 = px [p1];
                    const vect n1 = pn[p1];
                    const int type = protein.type[p1];

                    for ( int p2 = cl_lipid.cell_start[cell2]; p2 < cl_lipid.cell_start[cell2 + 1]; ++p2 ) {
                        const vect x2 = lx [p2];
                        const vect n2 = ln[p2];

                        const vect dx = x2 - x1; // 3 FLOPS
                        const real r_sq = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]; // 5 FLOPS

                        if ( r_sq < ForceField::cutsqlp[type] && r_sq > 1e-5 )
                            lipid_protein_omp( dx, r_sq, n2, n1, lipid, protein, p2, p1, type );
                        else if ( r_sq < ForceField::lj_cutsq[type] && r_sq > 1e-5 )
                            lennard_jones_omp( dx, r_sq, lipid, protein, p2, p1, type );
                    }
                }
            }
        }
    }

    Service<Timers>::call()["compute_pairwise_lp"].stop();
}

// protein-protein interactions
void compute_pairwise_pp( ProteContainer & protein, const VoronoiDiagram & voronoi, const VCellList & cl_protein ) {
    if ( !protein.size() ) return;

    Service<Timers>::call()["compute_pairwise_pp"].start();

    vect const * const __restrict x  = protein.x.data();
    int const * const __restrict type = protein.type.data();

    static std::vector<AlignedArray<int, true> > stencils( omp_get_max_threads() );

    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        const int nthreads = omp_get_num_threads();
        const int chunk_size = voronoi.n_cells / nthreads;
        const int cell_begin = chunk_size * tid;
        const int cell_end = ( tid == nthreads - 1 ) ? voronoi.n_cells : ( cell_begin + chunk_size );
        auto & stencil = stencils[ omp_get_thread_num() ];

        for ( int b = cell_begin; b < cell_end; ++b ) {

            if ( cl_protein.cell_start[b] == cl_protein.cell_start[b + 1] ) continue;

            // interaction between cells

            const int ncell = voronoi.get_stencil_whole( b, stencil, 9 );

            for ( int i = 0; i < ncell; ++i ) {
                const int cell2 = stencil[i];

                for ( int p1 = cl_protein.cell_start[b]; p1 < cl_protein.cell_start[b + 1]; ++p1 ) {
                    const vect x1 = x [p1];
                    const int type1 = type[p1];

                    for ( int p2 = cl_protein.cell_start[cell2]; p2 < cl_protein.cell_start[cell2 + 1]; ++p2 ) {
                        if ( p2 <= p1 ) continue;

                        const vect x2 = x[p2];
                        const int type2 = type[p2];

                        const vect dx = x1 - x2; // 3 FLOPS
                        const real r_sq = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]; // 5 FLOPS

                        const int type12 = type1 + type2 * ForceField::n_type;
                        if ( r_sq < ForceField::cutsqpp[type12] && r_sq > 1e-5 )
                            protein_protein_omp( dx, r_sq, protein, p1, p2, type12 );
                        else if ( r_sq < ForceField::lj_cutsq[type12] && r_sq > 1e-5 )
                            lennard_jones_omp( dx, r_sq, protein, protein, p1, p2, type12 );
                    }
                }
            }
        }
    }

    Service<Timers>::call()["compute_pairwise_pp"].stop();
}

}

#endif

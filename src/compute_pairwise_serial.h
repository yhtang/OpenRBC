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
#ifndef FORCEKERNEL_H
#define FORCEKERNEL_H

#include <cmath>
#include <iostream>
#include <omp.h>

#include "aligned_array.h"
#include "cell_list.h"
#include "container.h"
#include "math_vector.h"
#include "neigh_kdtree.h"
#include "forcefield.h"
#include "stencil.h"
#include "voronoi.h"
#include "pairwise_kernel.h"

//void compute_pairwise( Container& c, const Parameters& param, const CellList& lst, const Stencil& stencil )
//{
//  Service<Timers>::call()["compute_pairwise_1"].start();
//
//  real const * const __restrict x_x = c.x[0].data();
//    real const * const __restrict x_y = c.x[1].data();
//    real const * const __restrict x_z = c.x[2].data();
//    real const * const __restrict mu_x = c.n[0].data();
//    real const * const __restrict mu_y = c.n[1].data();
//    real const * const __restrict mu_z = c.n[2].data();
//
//    for( int b = 0; b < lst.get_total_bins(); ++b ) {
//        const int bin_size = lst.cell_start[b + 1] - lst.cell_start[b];
//        const int stencil_size = stencil.get_bin_size( b );
//        int const * current_bin = &( lst.cells[ lst.cell_start[b] ] );
//        for( int n = 0; n < bin_size; ++n ) {
//            const int p1 = *current_bin++;
//            const real x1 = x_x[p1];
//            const real y1 = x_y[p1];
//            const real z1 = x_z[p1];
//            const real mux1 = mu_x[p1];
//            const real muy1 = mu_y[p1];
//            const real muz1 = mu_z[p1];
//
//            int const *stcl = &( stencil._bins[0] ) + stencil._bin_starts[b] + n + 1;
//            for( int m = n + 1; m < stencil_size; ++m ) {
//                const int p2 = *stcl++;
//                const real x2 = x_x[p2];
//                const real y2 = x_y[p2];
//                const real z2 = x_z[p2];
//                const real mux2 = mu_x[p2];
//                const real muy2 = mu_y[p2];
//                const real muz2 = mu_z[p2];
//
//                const real dx = x1 - x2; // 1 FLOPS
//                const real dy = y1 - y2; // 1 FLOPS
//                const real dz = z1 - z2; // 1 FLOPS
//                const real r_sq = dx * dx + dy * dy + dz * dz; // 5 FLOPS
//
//                if( r_sq >= param.cutsqll || r_sq < 1e-5 ) continue;
//
//                lipid_lipid_serial(dx, dy, dz, r_sq, mux1, muy1, muz1,
//                                        mux2, muy2, muz2, c, param, p1, p2);
//            }
//        }
//    }
//
//    Service<Timers>::call()["compute_pairwise_1"].stop();
//}
//
//// particles must be ordered.
//void compute_pairwise( Container& c, const Parameters& param, const Voronoi<T4>& lst)
//{
//    Service<Timers>::call()["compute_pairwise_2"].start();
//
//    real const * const __restrict x_x = c.x[0].data();
//    real const * const __restrict x_y = c.x[1].data();
//    real const * const __restrict x_z = c.x[2].data();
//    real const * const __restrict mu_x = c.n[0].data();
//    real const * const __restrict mu_y = c.n[1].data();
//    real const * const __restrict mu_z = c.n[2].data();
//
//    AlignedArray<int, 16> stencil;
//
//    for( int b = 0; b < lst.n_cells; ++b ) {
//
//        // self interaction
//
//        for (int p1 = lst.cell_start[b]; p1 < lst.cell_start[b + 1]; ++p1) {
//            const real x1 = x_x[p1];
//            const real y1 = x_y[p1];
//            const real z1 = x_z[p1];
//            const real mux1 = mu_x[p1];
//            const real muy1 = mu_y[p1];
//            const real muz1 = mu_z[p1];
//
//            for (int p2 = p1 + 1; p2 < lst.cell_start[b + 1]; ++p2) {
//                const real x2 = x_x[p2];
//                const real y2 = x_y[p2];
//                const real z2 = x_z[p2];
//                const real mux2 = mu_x[p2];
//                const real muy2 = mu_y[p2];
//                const real muz2 = mu_z[p2];
//
//                const real dx = x1 - x2; // 1 FLOPS
//                const real dy = y1 - y2; // 1 FLOPS
//                const real dz = z1 - z2; // 1 FLOPS
//                const real r_sq = dx * dx + dy * dy + dz * dz; // 5 FLOPS
//
//                if( r_sq >= Parameters::cutsqll || r_sq < 1e-5 ) continue;
//
//                lipid_lipid_serial(dx, dy, dz, r_sq, mux1, muy1, muz1, mux2, muy2, muz2,
//                                        c, param, p1, p2);
//            }
//        }
//
//        // interaction between cells
//
//        const int ncell = lst.get_stencil(b, stencil);
//
//        for (int p1 = lst.cell_start[b]; p1 < lst.cell_start[b + 1]; ++p1) {
//            const real x1 = x_x[p1];
//            const real y1 = x_y[p1];
//            const real z1 = x_z[p1];
//            const real mux1 = mu_x[p1];
//            const real muy1 = mu_y[p1];
//            const real muz1 = mu_z[p1];
//
//            /* const real rx1 = x1 - lst.centroids[b][0];
//            const real ry1 = y1 - lst.centroids[b][1];
//            const real rz1 = z1 - lst.centroids[b][2]; */
//
//            for (int i = 0; i < ncell; ++i) {
//                const int cell2 = stencil[i];
//
//                // vect dcx = lst.centroids[b] - lst.centroids[cell2];
//
//                for (int p2 = lst.cell_start[cell2]; p2 < lst.cell_start[cell2 + 1]; ++p2) {
//                    const real x2 = x_x[p2];
//                    const real y2 = x_y[p2];
//                    const real z2 = x_z[p2];
//                    const real mux2 = mu_x[p2];
//                    const real muy2 = mu_y[p2];
//                    const real muz2 = mu_z[p2];
//
//                    const real dx = x1 - x2; // 1 FLOPS
//                    const real dy = y1 - y2; // 1 FLOPS
//                    const real dz = z1 - z2; // 1 FLOPS
//                    /* const real rx2 = x2 - lst.centroids[cell2][0];
//                    const real ry2 = y2 - lst.centroids[cell2][1];
//                    const real rz2 = z2 - lst.centroids[cell2][2];
//                    const real dx = rx1 - rx2 + dcx[0];
//                    const real dy = ry1 - ry2 + dcx[1];
//                    const real dz = rz1 - rz2 + dcx[2]; */
//
//                    const real r_sq = dx * dx + dy * dy + dz * dz; // 5 FLOPS
//
//                    if( r_sq >= param.cutsqll || r_sq < 1e-5 ) continue;
//
//                    lipid_lipid_serial(dx, dy, dz, r_sq, mux1, muy1, muz1, mux2, muy2, muz2,
//                                            c, param, p1, p2);
//                }
//            }
//        }
//    }
//
//    Service<Timers>::call()["compute_pairwise_2"].stop();
//}
//
//template<typename real, typename real, typename real, typename T4>
//void compute_pairwise( Container& c, const Parameters& param, NeighKdtree<T4>& neigh )
//{
//  Service<Timers>::call()["compute_pairwise_3"].start();
//
//  constexpr int nmax = 64;
//
//    real const * const __restrict x_x = c.x[0].data();
//    real const * const __restrict x_y = c.x[1].data();
//    real const * const __restrict x_z = c.x[2].data();
//    real const * const __restrict mu_x = c.n[0].data();
//    real const * const __restrict mu_y = c.n[1].data();
//    real const * const __restrict mu_z = c.n[2].data();
//
//    for (int p1 = 0; p1 < c.size(); ++p1) {
//        const real x1 = x_x[p1];
//        const real y1 = x_y[p1];
//        const real z1 = x_z[p1];
//        const real mux1 = mu_x[p1];
//        const real muy1 = mu_y[p1];
//        const real muz1 = mu_z[p1];
//
//        int neigh_pts[nmax];
//        const int npts = neigh.find_neigh(c, p1, param.cutll, neigh_pts, nmax);
//
//        for (int j = 0; j < npts; ++j) {
//            const int p2 = neigh_pts[j];
//            if (p2 <= p1) continue;
//
//            const real x2 = x_x[p2];
//            const real y2 = x_y[p2];
//            const real z2 = x_z[p2];
//            const real mux2 = mu_x[p2];
//            const real muy2 = mu_y[p2];
//            const real muz2 = mu_z[p2];
//
//            const real dx = x1 - x2; // 1 FLOPS
//            const real dy = y1 - y2; // 1 FLOPS
//            const real dz = z1 - z2; // 1 FLOPS
//            const real r_sq = dx * dx + dy * dy + dz * dz; // 5 FLOPS
//
//            if( r_sq >= param.cutsqll || r_sq < 1e-5 ) continue;
//
//            lipid_lipid_serial(dx, dy, dz, r_sq, mux1, muy1, muz1, mux2, muy2, muz2,
//                                    c, param, p1, p2);
//        }
//    }
//
//    Service<Timers>::call()["compute_pairwise_3"].stop();
//}

#endif /* FORCEKERNEL_H */

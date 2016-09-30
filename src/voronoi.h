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
#ifndef VORONOI_H_
#define VORONOI_H_

#include "aligned_array.h"
#include "container.h"
#include "kdtree.h"
#include "reorder.h"
#include "reorder_morton.h"
#include "util_numa.h"
#include "service.h"
#include <iostream>
#include <omp.h>
#include <vector>
#include <numeric>

namespace openrbc {

using namespace config;

struct VoronoiDiagram {

    const int n_cells; // number of voronoi cells

    KDTree<real, 3, int> tree;
    AlignedArray<vector<real, 3>, true> centroids;

    VoronoiDiagram( int n_cells_ ) : n_cells( n_cells_ ) {
        load_balancer( true  ).set_range( n_cells + 1 );
        load_balancer( false ).set_range( n_cells     );
    }
    ~VoronoiDiagram() {}
    LoadBalancer & load_balancer( bool extra_cell = false ) const {
        if ( extra_cell )
            return Service<BalancerMap>::call()[ "voronoi-e" ];
        else
            return Service<BalancerMap>::call()[ "voronoi" ];
    }

    // Initialize centroids position and set particle affiliation
    template<class CONTAINER, class CELL>
    void init( CONTAINER & cont, CELL & cell_list, RTParameter const & param, const int n_iterate ) {
        Service<Timers>::call()["VoronoiDiagram::init"].start();

        // Randomly picking initial guess
        const int delta = cont.size() / n_cells;
        centroids.resize( n_cells );
        for ( int i = 0, idxcentroid = delta / 2; i < n_cells; ++i, idxcentroid = ( idxcentroid + delta ) % cont.size() ) {
            for ( int d = 0; d < 3; ++d ) centroids[i][d] = cont.x[idxcentroid][d];
        }
        // Iterative clustering
        for ( int k = 0; k < n_iterate; ++k ) {
            std::cout << '.' << std::flush;
            reorder_morton( centroids, param );
            tree.build( centroids );
            cell_list.partition( cont, *this );
            update_centroid( cont, cell_list, [&cell_list]( int i ) {return cell_list.cells[i];} );
            if ( ( k & ( ~k + 1 ) ) == k ) reorder( cont, param, cell_list ); // magic: k is-power-of-2
        }
        Service<Timers>::call()["VoronoiDiagram::init"].stop();

        reorder( cont, param, cell_list );
    }

    template<class CELL>
    void update( Container & cont, CELL & cell_list, RTParameter const & param ) {
        Service<Timers>::call()["VoronoiDiagram::update"].start();

        update_centroid( cont, cell_list, []( int i ) {return i;} );
        if ( param.nstep % param.freq_sort_ctrd == 0 ) reorder_morton( centroids, param );
        tree.build( centroids );

        Service<Timers>::call()["VoronoiDiagram::update"].stop();
    }


    int get_stencil_half( const int cell, AlignedArray<int, true> & stencil, const real rmax ) const {
        constexpr int nmax = 32;

        stencil.resize( nmax );
        int n = 0;

        int cells[nmax];
        const int n_neigh = tree.find_within( centroids[cell], rmax, cells, nmax );
        for ( int i = 0; i < n_neigh; ++i ) {
            const int cell2 = cells[i];
            if ( cell2 > cell )
                stencil[n++] = cell2;
        }
        return n;
    }

    int get_stencil_whole( const int cell, AlignedArray<int, true> & stencil, const real rmax ) const {
        return tree.find_within( centroids[cell], rmax, stencil );
    }

    int refine_stencil( const int cell, const int n, AlignedArray<int, true> & stencil, const real rmax ) const {
        int j = 0;
        for ( int i = 0; i < n ; ++i )
            if ( normsq( centroids[ stencil[i] ] - centroids[cell] ) < rmax * rmax )
                stencil.shadow( j++ ) = stencil[i];
        stencil.swap();

        return j;
    }

protected:

    // iuput: points and a mapping from points to centroids
    // side effect: updates centroid location
    template<class MAP, class CELL>
    void update_centroid( Container const & cont, CELL & cell_list, MAP const & map ) {
        Service<Timers>::call()["|Voronoi::update_centroid"].start();

        auto & cell_start = cell_list.cell_start;
        #pragma omp parallel for
        for ( int i = 0; i < n_cells; ++i ) {
            vector<real, 3> center( 0, 0, 0 );
            for ( int j = cell_start[i]; j < cell_start[i + 1]; ++j ) {
                auto k = map( j );
                center += { cont.x[k][0], cont.x[k][1], cont.x[k][2] };
            }
            center *= real( 1.0 ) / ( cell_start[i + 1] - cell_start[i] );
            centroids[i] = center;
        }

        Service<Timers>::call()["|Voronoi::update_centroid"].stop();
    }
};

struct VCellList {

    int n_cells = 0; // number of voronoi cells
    AlignedArray<int> cells;       // store indices of particles for each cell in a compact layout
    AlignedArray<int> affiliation; // affiliation[i]: cell id of atom i
    AlignedArray<int> cell_start;  // the starting index of each bin in cells
    AlignedArray<int> local_index; // local_index[i] will list the ith particle's position within its cell

    //std::vector<std::vector<int> > stencil;

    template<class CONTAINER>
    void update( CONTAINER & cont, VoronoiDiagram const & voronoi, RTParameter const & param ) {
        Service<Timers>::call()["VCellList::update"].start();

        partition( cont, voronoi );

        Service<Timers>::call()["VCellList::update"].stop();

        reorder( cont, param, *this );
        if ( param.nstep % param.freq_sort_bond == 0 ) reorder_bond( cont, *this );
    }

    // update affiliation after reordering particles
    void update_particle_affiliation( const Container & cont ) {
        Service<Timers>::call()["VCellList::update_particle_affiliation"].start();

        #pragma omp parallel for
        for ( int i = 0; i < n_cells; ++i )
            for ( int j = cell_start[i]; j < cell_start[i + 1]; ++j )
                affiliation[j] = i;

        Service<Timers>::call()["VCellList::update_particle_affiliation"].stop();
    }

    // input: centroids, tree
    // side effect: affiliation, cells, cell_start
    void partition( Container const & cont, VoronoiDiagram const & voronoi ) {
        Service<Timers>::call()["|VCellList::tessellation"].start();

        Service<Timers>::call()["|Tessellation-1"].start();

        const size_t n = cont.size();
        n_cells = voronoi.n_cells;
        cell_start .assign( n_cells + 1, 0, voronoi.load_balancer( true ) );
        cells      .resize( n );
        affiliation.resize( n );
        local_index.resize( n );

        Service<Timers>::call()["|Tessellation-1"].stop();

        Service<Timers>::call()["|Tessellation-2"].start();

        // update affiliation
        #pragma omp parallel
        {
            int last_nearest = 0;
            real r;

            for ( auto itr = Service<BalancerMap>::call()[ cont.id() ].range(); !itr.exhausted ; ++itr ) {
                auto i = *itr;
                vector<real, 3> pt { cont.x[i][0], cont.x[i][1], cont.x[i][2] };
                int  guess_i = last_nearest;
                real guess_r = norm( pt - voronoi.centroids[last_nearest] );
                auto update_guess = [&]( int g ) {
                    real r = norm( pt - voronoi.centroids[ g ] );
                    if ( r < guess_r ) guess_r = r, guess_i = g;
                };
                // perform addition initial guess, huge speedup at large particle count
                if ( affiliation[i] >= 0 && affiliation[i] < n_cells ) update_guess( affiliation[i] );
                if ( last_nearest + 1 < n_cells ) update_guess( last_nearest + 1 );
                if ( last_nearest - 1 >= 0 ) update_guess( last_nearest - 1 );
                std::tie( r, affiliation[i] ) = voronoi.tree.find_nearest( pt, std::make_tuple( guess_r, guess_i ) );
                #pragma omp atomic capture
                local_index[i] = cell_start[ affiliation[i] + 1 ]++;
                last_nearest = affiliation[i];
            }
        }

        Service<Timers>::call()["|Tessellation-2"].stop();

        Service<Timers>::call()["|Tessellation-3"].start();

        // make cell_start cumulative (number of particles in cells before i)
        std::partial_sum( cell_start.data() + 1, cell_start.data() + n_cells + 1, cell_start.data() + 1 );

        Service<Timers>::call()["|Tessellation-3"].stop();

        Service<Timers>::call()["|Tessellation-4"].start();
        // update cells array to put atoms in order
        #pragma omp parallel for
        for ( std::size_t i = 0; i < n; ++i )
            cells[local_index[i] + cell_start[affiliation[i]]] = i;
        Service<Timers>::call()["|Tessellation-4"].stop();

        Service<Timers>::call()["|VCellList::tessellation"].stop();
    }
};

}

#endif

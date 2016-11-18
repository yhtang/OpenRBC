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
#ifndef OPENRBC_REORDER_H_
#define OPENRBC_REORDER_H_

#include "container.h"
#include "forcefield.h"
#include "runtime_parameter.h"
#include "timer.h"
#include "util_numa.h"
#include "voronoi.h"
#include <numeric>

namespace openrbc {

// reorder the bonds of particles within the container by the position
// need to move the bonds data
// out-of-place reorder
template<class C>
inline void reorder_bond( LipidContainer & model, C & cell_list ) {}

template<class C>
inline void reorder_bond( ProteContainer & model, C & cell_list ) {
    Service<Timers>::call()["|reorder_bond"].start();

    static AlignedArray<int> values;
    static AlignedArray<int> cell_start;
    static AlignedArray<int> local_index;

    const auto n = model.bonds.size();
    const auto n_cells = cell_list.n_cells;

    values.resize( n );
    cell_start.assign( n_cells + 1, 0 );
    local_index.resize( n );

    cell_list.update_particle_affiliation( model );

    #pragma omp parallel for
    for ( int i = 0; i < n; ++i ) {
        values[i] = cell_list.affiliation[model.tag2idx[model.bonds[i].i]];
        #pragma omp atomic capture
        local_index[i] = cell_start[ values[i] + 1 ]++;
    }

    std::partial_sum( cell_start.data() + 1, cell_start.data() + n_cells + 1, cell_start.data() + 1 );

    #pragma omp parallel for
    for ( int i = 0; i < n; ++i )
        model.bonds.shadow( local_index[i] + cell_start[values[i]] ) = model.bonds[i];

    model.bonds.swap();

    Service<BalancerMap>::call()[ model.id() + "-bonds" ].set_range( model.bonds.size() );

    Service<Timers>::call()["|reorder_bond"].stop();
}

// reorder the position of particles within the container by their bin ID
// need to move the x, v, mu, o, type and tag data
// out-of-place reorder
template<class C>
inline void reorder( ProteContainer & model, RTParameter const & param, const C & cell_list ) {
    Service<Timers>::call()["reorder"].start();

    model.reserve( model.size() + ( 29 - ( model.size() % 16 ) ) );

    Service<BalancerMap>::call()["reorder"].set_range( cell_list.cell_start.size() - 1 );

    std::size_t q = 0;

    #pragma omp parallel reduction( max : q )
    {
        auto beg = Service<BalancerMap>::call()["reorder"].beg();
        auto end = Service<BalancerMap>::call()["reorder"].end();
        q = 0;
        for ( auto i = beg; i < end; i++ ) {
            for ( int j = cell_list.cell_start[i]; j < cell_list.cell_start[i + 1]; j++ ) {
                model.x.shadow( j ) = model.x[cell_list.cells[j]];
                model.v.shadow( j ) = model.v[cell_list.cells[j]];
                model.n.shadow( j ) = model.n[cell_list.cells[j]];
                model.o.shadow( j ) = model.o[cell_list.cells[j]];
                model.tag .shadow( j ) = model.tag  [cell_list.cells[j]];
                model.type.shadow( j ) = model.type [cell_list.cells[j]];
            }
            q = cell_list.cell_start[i + 1];
        }
    }

    if ( q != model.size() ) {
        fprintf( stderr, "<OpenRBC> Reorder failed, expecting %ld protein, got %ld\n", model.size(), q );
        raise( SIGSEGV );
    }

    model.swap();

    Service<Timers>::call()["reorder"].stop();

    model.tag2idx.build_map( model );
}

template<class C>
inline void reorder( LipidContainer & model, RTParameter const & param, const C & cell_list ) {
    Service<Timers>::call()["reorder"].start();

    model.reserve( model.size() + ( 29 - ( model.size() % 16 ) ) );

    Service<BalancerMap>::call()["reorder"].set_range( cell_list.cell_start.size() - 1 );

    std::size_t q = 0;

    #pragma omp parallel reduction( max : q )
    {
        auto beg = Service<BalancerMap>::call()["reorder"].beg();
        auto end = Service<BalancerMap>::call()["reorder"].end();
        q = 0;
        for ( auto i = beg; i < end; i++ ) {
            for ( int j = cell_list.cell_start[i]; j < cell_list.cell_start[i + 1]; j++ ) {
                model.x.shadow( j ) = model.x[cell_list.cells[j]];
                model.v.shadow( j ) = model.v[cell_list.cells[j]];
                model.n.shadow( j ) = model.n[cell_list.cells[j]];
                model.o.shadow( j ) = model.o[cell_list.cells[j]];
            }
            q = cell_list.cell_start[i + 1];
        }
    }

    if ( q != model.size() ) {
        fprintf( stderr, "<OpenRBC> Reorder failed, expecting %ld lipid, got %ld\n", model.size(), q );
        raise( SIGSEGV );
    }

    model.swap();

    Service<Timers>::call()["reorder"].stop();

    model.tag2idx.build_map( model );
}

}

#endif

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
#ifndef OPENRBC_CONTAINER_H_
#define OPENRBC_CONTAINER_H_

#include <unordered_map>
#include "aligned_array.h"
#include "service.h"
#include "timer.h"
#include "math_vector.h"
#include "math_vector.h"
#include "util_numa.h"
#include "config_static.h"

namespace openrbc {

using namespace config;

struct Bond {
    int type, i, j;
    Bond() = default;
    Bond( int type_, int i_, int j_ ) : type( type_ ), i( i_ ), j( j_ ) {}
};

struct TagIndexMap {
    // construct the map
    template<class CONTAINER>
    void build_map( CONTAINER const & c ) {
        Service<Timers>::call()["TagIndexMap::build_map"].start();

        #pragma omp parallel
        {
            auto beg = Service<BalancerMap>::call()[ c.id() ].beg();
            auto end = Service<BalancerMap>::call()[ c.id() ].end();

            for ( std::size_t i = beg; i < end; ++i )
                if ( c.tag[i] >= map.size() )
                    #pragma omp critical (TagIndexMap_build_map)
                    if ( c.tag[i] >= map.size() ) map.resize( std::max<std::size_t>( c.tag[i] + 1, map.size() * 2 ) );

            #pragma omp barrier

            for ( int i = beg; i < end; ++i ) map[ c.tag[i] ] = i;
        }

        Service<Timers>::call()["TagIndexMap::build_map"].stop();
    }
    // return index given tag
    int operator []( int tag ) const { return map[tag]; }
private:
    AlignedArray<int> map;
};

struct Container {

    // double-buffer for efficient particle reordering
    AlignedArray<vect, true> x; // position
    AlignedArray<vect, true> v; // velocity
    AlignedArray<vect, true> f; // force
    AlignedArray<vect, true> n; // angular position
    AlignedArray<vect, true> o; // angular velocity
    AlignedArray<vect, true> t; // t

    Container( std::string id ) : id_( id ) {
        Service<BalancerMap>::call()[ id_ ].set_range( n_ );
    }

    std::size_t size() const { return n_; }
    std::string id() const { return id_; }
    std::size_t capacity() const { return x.capacity(); }

    void resize( std::size_t sz ) {
        x.resize( sz );
        v.resize( sz );
        f.resize( sz );
        n.resize( sz );
        o.resize( sz );
        t.resize( sz );
        n_ = sz;
        Service<BalancerMap>::call()[ id_ ].set_range( n_ );
    }
    void reserve( std::size_t sz ) {
        x.reserve( sz );
        v.reserve( sz );
        f.reserve( sz );
        n.reserve( sz );
        o.reserve( sz );
        t.reserve( sz );
    }
    void swap( bool xforce = false, bool xtorque = false ) {
        x.swap();
        v.swap();
        n.swap();
        o.swap();
        if ( xforce )  f.swap();
        if ( xtorque ) t.swap();
    }

    TagIndexMap tag2idx;

protected:
    std::size_t n_ = 0;
    std::string id_;
};

struct LipidContainer : Container {
    using Base = Container;

    LipidContainer( std::string id ) : Base( id ) {}

    struct {
        int operator [] ( int i ) const { return base + i; }
        void set_base( int b ) { base = b; }
        int base = 1;
    } tag;

    struct {
        int operator [] ( int i ) const { return 0; }
    } const type {};
};

struct ProteContainer : Container {
    using Base = Container;

    AlignedArray<int, true> type;
    AlignedArray<int, true> tag;
    AlignedArray<Bond, true> bonds;

    ProteContainer( std::string id ) : Base( id ) {}

    void resize( std::size_t n ) {
        type.resize( n );
        tag .resize( n );
        Base::resize( n );
    }
    void reserve( std::size_t n ) {
        type.reserve( n );
        tag .reserve( n );
        Base::reserve( n );
    }
    void swap( bool xforce = false, bool xtorque = false ) {
        type.swap();
        tag .swap();
        Base::swap( xforce, xtorque );
    }
};

}

#endif

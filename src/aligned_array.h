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
#ifndef ALIGNED_ARRAY
#define ALIGNED_ARRAY

#include "util_numa.h"
#include "util_misc.h"

#include <cstdlib>
#include <cstring>
#include <cassert>
#include <algorithm>
#include <csignal>

namespace openrbc {

template<typename TYPE, bool SHADOW = false, std::size_t ALIGN = 64>
struct alignas(128) AlignedArray  {
protected:
    //static_assert( std::is_pod<TYPE>::value, "AlignedArray can only be used for POD types" );
    std::size_t size_ = 0, capacity_ = 0;
    alignas( ALIGN ) TYPE * __restrict data_ = nullptr, * __restrict shadow_ = nullptr;

public:
    AlignedArray( std::size_t s = 0 ) {
        resize( s );
    }
    ~AlignedArray() {
        if ( data_ ) free( data_ );
        data_ = shadow_ = nullptr;
        capacity_ = size_ = 0;
    }

    AlignedArray( const AlignedArray & a ) = delete; // forbid automatic copy-initialize
    AlignedArray & operator= ( const AlignedArray & a ) = delete; // forbid automatic copy-assign

    std::size_t size    () const { return size_;     }
    std::size_t capacity() const { return capacity_; }
    TYPE * data  () const { return data_; }
    TYPE * shadow() const { if ( !SHADOW ) raise( SIGSEGV ); return shadow_; }
    TYPE    &    operator []( std::size_t i )       { return data_[i]; }
    TYPE const & operator []( std::size_t i ) const { return data_[i]; }
    TYPE    &    shadow( std::size_t i )       { if ( !SHADOW ) raise( SIGSEGV ); return shadow_[i]; }
    TYPE const & shadow( std::size_t i ) const { if ( !SHADOW ) raise( SIGSEGV ); return shadow_[i]; }
    template<typename ...T> void emplace_back( T && ... initializer ) {
        std::size_t i = size();
        resize( size() + 1 );
        new( data_ + i ) TYPE( initializer... );
    }

    void resize( std::size_t sz, bool copy = true ) {
        if ( capacity_ >= sz ) {
            size_ = sz;
        } else {
            std::size_t cap_ = std::max<decltype( sz )>( sz, capacity_ * 1.4 );
            extend( cap_, copy );
            capacity_ = cap_;
            size_ = sz;
        }
    }
    void reserve( std::size_t cp, bool copy = true ) {
        if ( capacity_ < cp ) {
            extend( cp, copy );
            capacity_ = cp;
        }
    }
    void assign( std::size_t sz, TYPE val, LoadBalancer const & balancer ) {
        resize( sz );
        #pragma omp parallel
        for ( std::size_t i = balancer.beg(); i < balancer.end(); ++i ) data_[i] = val;
    }
    void assign( std::size_t sz, TYPE val ) {
        resize( sz );
        #pragma omp parallel for
        for ( std::size_t i = 0; i < sz; i++ ) data_[i] = val;
    }
    void swap( AlignedArray & other ) {
        std::swap( size_, other.size_ );
        std::swap( capacity_, other.capacity_ );
        std::swap( data_, other.data_ );
        if ( SHADOW ) std::swap( shadow_, other.shadow_ );
    }
    void swap() {
        if ( SHADOW ) std::swap( data_, shadow_ );
    }

protected:
    void extend( std::size_t sz, bool copy ) {
        alignas( ALIGN ) TYPE * __restrict old_data = data_;
        alignas( ALIGN ) void * new_data = aligned_alloc( ALIGN, sz * sizeof( TYPE ) );
        data_ = ( TYPE * __restrict )( new_data );
        if ( copy && old_data ) {
            memcpy( data_, old_data, ( ( size_ > sz ) ? sz : size_ ) * sizeof( TYPE ) );
        }
        if ( old_data ) free ( old_data );
        if ( SHADOW ) {
            if ( shadow_ ) free ( shadow_ );
            shadow_ = ( TYPE * __restrict ) aligned_alloc( ALIGN, sz * sizeof( TYPE ) );
        }
    }
};

}

#endif /* ALIGNED_ARRAY_H_ */

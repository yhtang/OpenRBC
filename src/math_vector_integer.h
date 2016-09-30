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
#ifndef MATH_VECTOR_INTEGER_H_
#define MATH_VECTOR_INTEGER_H_

#include "math_vector_base.h"

namespace openrbc {

template<uint D> struct vector<uint, D> {
protected:
    uint x[D];
public:
    using TYPE = uint;
    static const uint d = D;

    inline vector() {}
    inline vector( uint const & v ) {
        for ( uint i = 0; i < D; ++i ) x[i] = v;
    }
    template<typename ...T>
    inline vector( uint const & first, uint const & second, T const & ... tail ) {
        static_assert( sizeof...( T ) + 2 >= D, "Too few initializers for openrbc::vector" );
        std::array < TYPE, sizeof...( T ) + 2 > s { first, second, static_cast<uint>( tail )... };
        for ( uint i = 0 ; i < D ; ++i ) x[i] = s[i];
    }
    inline vector( vector const & other ) = default;
    inline vector & operator = ( vector const & other ) = default;

    // point must be assignable, while other expressions may not
    inline uint    &    operator [] ( uint i )       { return x[i]; }
    inline uint const & operator [] ( uint i ) const { return x[i]; }

    /*---------------------------------------------------------------------------
                             Bitwise operators
    ---------------------------------------------------------------------------*/
    inline vector & operator ^= ( vector const & u ) {
        for ( uint i = 0 ; i < D ; ++i ) x[i] ^= u[i];
        return *this;
    }
    friend inline vector operator >> ( vector const & u, int disp ) {
        vector w;
        for ( uint i = 0; i < D; ++i ) w[i] = u[i] >> disp;
        return w;
    }
    friend inline vector operator << ( vector const & u, int disp ) {
        vector w;
        for ( uint i = 0; i < D; ++i ) w[i] = u[i] << disp;
        return w;
    }
    friend inline vector<float, D> uint2u11( vector const & u ) {
        vector<float, D> w;
        for ( uint i = 0; i < D; ++i ) w[i] = u[i] * 4.6566129e-10f - 1.f;
        return w;
    }
};

}

#endif

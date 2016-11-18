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
#ifndef OPENRBC_MATH_VECTOR_BASE_H_
#define OPENRBC_MATH_VECTOR_BASE_H_

#include <array>
#include <cmath>
#include <cassert>
#include <type_traits>
#include <algorithm>
#include <cstdint>

namespace openrbc {

using uint = std::uint32_t;

template<typename SCALAR, uint D>
struct vector {
    /*---------------------------------------------------------------------------
                                     Container
    ---------------------------------------------------------------------------*/
protected:
    SCALAR x[D];
public:
    using TYPE = SCALAR;
    static const uint d = D;

    inline vector() {}
    inline vector( SCALAR const & v ) {
        for ( uint i = 0; i < D; ++i ) x[i] = v;
    }
    inline vector( vector const & other ) {
        for ( uint i = 0; i < D; ++i ) x[i] = other.x[i];
    }
    template<typename ...T>
    inline vector( SCALAR const & first, SCALAR const & second, T const & ... tail ) {
        static_assert( sizeof...( T ) + 2 >= D, "Too few initializers for openrbc::vector" );
        std::array < TYPE, sizeof...( T ) + 2 > s { first, second, static_cast<SCALAR>( tail )... };
        for ( uint i = 0 ; i < D ; ++i ) x[i] = s[i];
    }

    // copy-assign
    inline vector & operator = ( vector const & other ) {
        for ( uint i = 0 ; i < D ; ++i ) x[i] = other[i];
        return *this;
    }
    template<typename T> inline vector( vector<T, D> const & other ) {
        for ( uint i = 0; i < D; ++i ) x[i] = static_cast<T>( other[i] );
    }

    // point must be assignable, while other expressions may not
    inline SCALAR    &   operator [] ( uint i )       { return x[i]; }
    inline SCALAR const & operator [] ( uint i ) const { return x[i]; }

    static inline vector zero() {
        static vector zero_ {0.0f};
        return zero_;
    }

    /*---------------------------------------------------------------------------
                             Operator Overloads
    ---------------------------------------------------------------------------*/
    // OP-Assign operators
    inline vector & operator += ( vector const & u ) {
        for ( uint i = 0 ; i < D ; ++i ) x[i] += u[i];
        return *this;
    }
    inline vector & operator -= ( vector const & u ) {
        for ( uint i = 0 ; i < D ; ++i ) x[i] -= u[i];
        return *this;
    }
    inline vector & operator *= ( vector const & u ) {
        for ( uint i = 0 ; i < D ; ++i ) x[i] *= u[i];
        return *this;
    }
    inline vector & operator /= ( vector const & u ) {
        for ( uint i = 0 ; i < D ; ++i ) x[i] /= u[i];
        return *this;
    }
    // Vector-Scalar operators
    inline vector & operator += ( SCALAR const u ) {
        for ( uint i = 0 ; i < D ; ++i ) x[i] += u;
        return *this;
    }
    inline vector & operator -= ( SCALAR const u ) {
        for ( uint i = 0 ; i < D ; ++i ) x[i] -= u;
        return *this;
    }
    inline vector & operator *= ( SCALAR const u ) {
        for ( uint i = 0 ; i < D ; ++i ) x[i] *= u;
        return *this;
    }
    inline vector & operator /= ( SCALAR const u ) {
        for ( uint i = 0 ; i < D ; ++i ) x[i] /= u;
        return *this;
    }

    friend inline vector operator + ( vector const & u, vector const & v ) {
        vector w;
        for ( uint i = 0; i < D; ++i ) w[i] = u[i] + v[i];
        return w;
    }
    friend inline vector operator - ( vector const & u, vector const & v ) {
        vector w;
        for ( uint i = 0; i < D; ++i ) w[i] = u[i] - v[i];
        return w;
    }
    friend inline vector operator - ( vector const & u ) {
        vector w;
        for ( uint i = 0; i < D; ++i ) w[i] = -u[i];
        return w;
    }
    friend inline vector operator + ( vector const & u, SCALAR a ) {
        vector w;
        for ( uint i = 0; i < D; ++i ) w[i] = u[i] + a;
        return w;
    }
    friend inline vector operator - ( vector const & u, SCALAR a ) {
        vector w;
        for ( uint i = 0; i < D; ++i ) w[i] = u[i] - a;
        return w;
    }
    friend inline vector operator * ( vector const & u, vector const & v ) {
        vector w;
        for ( uint i = 0; i < D; ++i ) w[i] = u[i] * v[i];
        return w;
    }
    friend inline vector operator * ( vector const & u, SCALAR a ) {
        vector w;
        for ( uint i = 0; i < D; ++i ) w[i] = u[i] * a;
        return w;
    }
    friend inline vector operator * ( SCALAR a, vector const & u ) {
        return u * a;
    }
    friend inline vector operator / ( vector const & u, vector const & v ) {
        vector w;
        for ( uint i = 0; i < D; ++i ) w[i] = u[i] / v[i];
        return w;
    }
    friend inline vector operator / ( vector const & u, SCALAR a ) {
        vector w;
        for ( uint i = 0; i < D; ++i ) w[i] = u[i] / a;
        return w;
    }
    friend inline vector operator / ( SCALAR a, vector const & u ) {
        vector w;
        for ( uint i = 0; i < D; ++i ) w[i] = a / u[i];
        return w;
    }
    template<typename OP> friend inline vector apply( vector const & u, OP const & o ) {
        vector w;
        for ( uint i = 0; i < D; ++i ) w[i] = o( u[i] );
        return w;
    }
    template<typename OP> friend inline vector apply( vector const & u, vector const & v, OP const & o ) {
        vector w;
        for ( uint i = 0; i < D; ++i ) w[i] = o( u[i], v[i] );
        return w;
    }

    /*---------------------------------------------------------------------------
                             Math functions
    ---------------------------------------------------------------------------*/
    // generic reduction template
    template<class OP> friend inline SCALAR reduce( vector const & u, OP const & op ) {
        SCALAR core( u[0] );
        for ( uint i = 1 ; i < D ; i++ ) core = op( core, u[i] );
        return core;
    }
    // biggest element within a vector
    friend inline SCALAR max( vector const & u ) {
        return reduce( u, []( SCALAR a, SCALAR b ) {return a > b ? a : b;} );
    }
    // index of biggest element within a vector
    friend inline uint argmax( vector const & u ) {
        SCALAR max = std::numeric_limits<SCALAR>::min();
        uint m;
        for ( uint i = 0; i < D; ++i ) if ( u[i] > max ) { max = u[i]; m = i; }
        return m;
    }
    // smallest element within a vector
    friend inline SCALAR min( vector const & u ) {
        return reduce( u, []( SCALAR a, SCALAR b ) {return a < b ? a : b;} );
    }
    // sum of elements
    friend inline SCALAR sum( vector const & u ) {
        return reduce( u, []( SCALAR a, SCALAR b ) {return a + b;} );
    }
    // mean of elements
    friend inline SCALAR mean( vector const & u ) {
        return sum( u ) / double( D );
    }
    // inner product
    friend inline SCALAR dot( vector const & u, vector const & v ) {
        SCALAR s = 0;
        for ( uint i = 0; i < D; ++i ) s += u[i] * v[i];
        return s;
    }
    // square of L2 norm
    friend inline SCALAR normsq( vector const & u ) {
        SCALAR s = 0;
        for ( uint i = 0; i < D; ++i ) s += u[i] * u[i];
        return s;
    }
    // L2 norm
    friend inline SCALAR norm( vector const & u ) {
        SCALAR s = 0;
        for ( uint i = 0; i < D; ++i ) s += u[i] * u[i];
        return std::sqrt( s );
    }
    // inverse of L2 norm
    friend inline SCALAR norm_inv( vector const & u ) {
        SCALAR s = 0;
        for ( uint i = 0; i < D; ++i ) s += u[i] * u[i];
        return SCALAR( 1 ) / std::sqrt( s );
    }
    // normalize
    friend inline vector normalize( vector const & u ) {
        return u * norm_inv( u );
    }
    // functions that should never be used, but necessary for compiler syntax parsing
    friend inline vector _dot   ( vector const &, vector const & ) { return 0xDEAD; }
    friend inline vector _rsqrt ( vector const & ) { return 0xDEAD; }
    friend inline vector _normsq( vector const & ) { return 0xDEAD; }
};

template<typename SCALAR, uint D> inline
vector<SCALAR, D> cross( vector<SCALAR, D> const & u, vector<SCALAR, D> const & v ) {
    static_assert( D >=3, "Cross product only defined between 3D vectors" );
	vector<SCALAR, D> x;
    for ( uint i = 0; i < 3; ++i ) x[i] = u[( i + 1U ) % 3U] * v[( i + 2U ) % 3U] - u[( i + 2U ) % 3U] * v[( i + 1U ) % 3U];
    for ( uint i = 3; i < D; ++i ) x[i] = 0;
    return x;
}

namespace functors {
template<typename T> struct min { inline T operator () ( T const & u, T const & v ) const { return std::min<T>( u, v ); } };
template<typename T> struct max { inline T operator () ( T const & u, T const & v ) const { return std::max<T>( u, v ); } };
}

}

// force std to treat vector as a POD type
namespace std {
template<typename SCALAR, uint D> struct is_pod<openrbc::vector<SCALAR, D> > { const static bool value = true; };
}

#endif

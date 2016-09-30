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
#ifndef MATH_VECTOR_SIMD_H_
#define MATH_VECTOR_SIMD_H_

#include "config_static.h"
#include "math_vector_base.h"
#include <array>
#include <cmath>
#include <cassert>
#include <type_traits>
#include <algorithm>
#include <cstdint>

namespace openrbc {

#ifdef EXPLICIT_SIMD

#if defined( __x86_64__ ) || defined( _M_X64 )

#include <x86intrin.h>

template<> struct vector<uint, 4> {
protected:
    union { __v4su m128; uint x[4]; };
    friend vector<float, 4>;
public:
    using TYPE = uint;
    static const uint d = 4;

    inline vector() {}
    inline vector( __v4su && m ) : m128( m ) {}
    inline vector( uint const & v ) {
        for ( uint i = 0; i < 4; ++i ) x[i] = v;
    }
    inline vector( uint x, uint y, uint z, uint w ) {
        this->x[0] = x, this->x[1] = y, this->x[2] = z, this->x[3] = w;
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
        m128 ^= u.m128;
        return *this;
    }
    friend inline vector operator >> ( vector const & u, int disp ) {
        return u.m128 >> disp;
    }
    friend inline vector operator << ( vector const & u, int disp ) {
        return u.m128 << disp;
    }
    friend inline vector<float, 4> uint2u11( vector const & u );
};

template<> struct vector<float, 4> {
protected:
    union { __m128 m128; float x[4]; };
public:
    static const uint d = 4;

    // default constructor
    inline vector() {}
    inline vector( __m128 && m ) : m128( m ) {
    }
    inline vector( float u ) {
        m128 = _mm_set1_ps( u );
    }
    inline vector( float x, float y, float z, float w ) {
        m128 = _mm_set_ps( w, z, y, x );
    }
    inline vector( vector const & other ) {
        m128 = other.m128;
    }
    inline vector & operator = ( vector const & other ) {
        m128 = other.m128;
        return *this;
    }
    inline vector & operator = ( float u ) {
        m128 = _mm_set1_ps( u );
        return *this;
    }

    inline float    &    operator [] ( uint i )       { return x[i]; }
    inline float const & operator [] ( uint i ) const { return x[i]; }

    static inline vector zero() {
        static vector zero_ {0.0f};
        return zero_;
    }

    /*---------------------------------------------------------------------------
                             Operator Overloads
    ---------------------------------------------------------------------------*/
    // OP-Assign operators
    inline vector & operator += ( vector const & u ) {
        m128 = _mm_add_ps( m128, u.m128 );
        return *this;
    }
    inline vector & operator -= ( vector const & u ) {
        m128 = _mm_sub_ps( m128, u.m128 );
        return *this;
    }
    inline vector & operator *= ( vector const & u ) {
        m128 = _mm_mul_ps( m128, u.m128 );
        return *this;
    }
    inline vector & operator /= ( vector const & u ) {
        m128 = _mm_div_ps( m128, u.m128 );
        return *this;
    }
    friend inline vector operator + ( vector const & u, vector const & v ) {
        return _mm_add_ps( u.m128, v.m128 );
    }
    friend inline vector operator - ( vector const & u, vector const & v ) {
        return _mm_sub_ps( u.m128, v.m128 );
    }
    friend inline vector operator * ( vector const & u, vector const & v ) {
        return _mm_mul_ps( u.m128, v.m128 );
    }
    friend inline vector operator / ( vector const & u, vector const & v ) {
        return _mm_div_ps( u.m128, v.m128 );
    }

    /*---------------------------------------------------------------------------
                             Math functions
    ---------------------------------------------------------------------------*/
    friend inline vector _rsqrt( vector const & u ) {
        return _mm_rsqrt_ps( u.m128 );
    }
    friend inline vector _normsq( vector const & u ) {
        return _mm_dp_ps( u.m128, u.m128, 0xFF );
    }
    friend inline vector _dot( vector const & u, vector const & v ) {
        return _mm_dp_ps( u.m128, v.m128, 0xFF );
    }
    friend inline vector _cross( vector const & u, vector const & v ) {
        __m128 a_yzx = _mm_shuffle_ps( u.m128, u.m128, _MM_SHUFFLE( 3, 0, 2, 1 ) );
        __m128 b_yzx = _mm_shuffle_ps( v.m128, v.m128, _MM_SHUFFLE( 3, 0, 2, 1 ) );
        __m128 c = _mm_sub_ps( _mm_mul_ps( u.m128, b_yzx ), _mm_mul_ps( a_yzx, v.m128 ) );
        return _mm_shuffle_ps( c, c, _MM_SHUFFLE( 3, 0, 2, 1 ) );
    }
    /*---------------------------------------------------------------------------
                Math functions NOT TO BE USED IN FORCE KERNEL
    ---------------------------------------------------------------------------*/
    friend inline float normsq( vector const & u ) {
        return _mm_dp_ps( u.m128, u.m128, 0xF1 )[0];
    }
    friend inline float dot( vector const & u, vector const & v ) {
        return _mm_dp_ps( u.m128, v.m128, 0xF1 )[0];
    }
    friend inline vector cross( vector const & a, vector const & b ) {
        __m128 a_yzx = _mm_shuffle_ps( a.m128, a.m128, _MM_SHUFFLE( 3, 0, 2, 1 ) );
        __m128 b_yzx = _mm_shuffle_ps( b.m128, b.m128, _MM_SHUFFLE( 3, 0, 2, 1 ) );
        __m128 c = _mm_sub_ps( _mm_mul_ps( a.m128, b_yzx ), _mm_mul_ps( a_yzx, b.m128 ) );
        return _mm_shuffle_ps( c, c, _MM_SHUFFLE( 3, 0, 2, 1 ) );
    }
    friend inline float norm( vector const & u ) {
        __m128 nsq = _mm_dp_ps( u.m128, u.m128, 0xF1 );
        return ( nsq * _mm_rsqrt_ss( nsq ) )[0];
    }
    friend inline vector normalize( vector const & u ) {
        return u * _rsqrt( _normsq( u ) );
    }
    template<typename OP> friend inline vector apply( vector const & u, OP const & o ) {
        vector w;
        for ( uint i = 0; i < 4; ++i ) w[i] = o( u[i] );
        return w;
    }
    template<typename OP> friend inline vector apply( vector const & u, vector const & v, OP const & o ) {
        vector w;
        for ( uint i = 0; i < 4; ++i ) w[i] = o( u[i], v[i] );
        return w;
    }
    // index of biggest element within a vector
    friend inline uint argmax( vector const & u ) {
        float max = -std::numeric_limits<float>::max();
        uint m;
        for ( uint i = 0; i < 4; ++i ) if ( u[i] > max ) { max = u[i]; m = i; }
        return m;
    }
};

inline vector<float, 4> uint2u11( vector<uint, 4> const & u ) {
    return _mm_cvtepi32_ps( u.m128 ) * 4.6566129e-10f;
}

#endif

#if defined(__powerpc64__) || defined(__ppc64__) || defined(__PPC64__)

#include <altivec.h>
#undef vector // use __vector instead to avoid conflict with std::vector
#undef bool

template<> struct vector<float, 4> {
protected:
    union { __vector float m128; float x[4]; };
public:
    static const uint d = 4;

    // default constructor
    inline vector() {}
    inline vector( __vector float m ) : m128( m ) {
    }
    inline vector( float u ) {
        m128 = vec_splats( u );
    }
    inline vector( float x, float y, float z, float w ) {
        m128[0] = x;
        m128[1] = y;
        m128[2] = z;
        m128[3] = w;
    }
    explicit inline vector( float * x ) {
        m128 = vec_ld( 0, x );
    }
    inline vector( vector const & other ) {
        m128 = other.m128;
    }
    inline vector & operator = ( vector const & other ) {
        m128 = other.m128;
        return *this;
    }
    inline vector & operator = ( float u ) {
        m128 = vec_splats( u );
        return *this;
    }

    inline float    &    operator [] ( uint i )       { return x[i]; }
    inline float const & operator [] ( uint i ) const { return x[i]; }

    static inline vector zero() {
        static vector zero_ { -0.0f};
        return zero_;
    }

    /*---------------------------------------------------------------------------
                             Operator Overloads
    ---------------------------------------------------------------------------*/
    // OP-Assign operators
    inline vector & operator += ( vector const & u ) {
        m128 = vec_add( m128, u.m128 );
        return *this;
    }
    inline vector & operator -= ( vector const & u ) {
        m128 = vec_sub( m128, u.m128 );
        return *this;
    }
    inline vector & operator *= ( vector const & u ) {
        m128 = vec_madd( m128, u.m128, zero().m128 );
        return *this;
    }
//    inline vector & operator /= ( vector const & u ) {
//        m128 = _mm_div_ps( m128, u.m128 );
//        return *this;
//    }
    friend inline vector operator + ( vector const & u, vector const & v ) {
        return vec_add( u.m128, v.m128 );
    }
    friend inline vector operator - ( vector const & u, vector const & v ) {
        return vec_sub( u.m128, v.m128 );
    }
    friend inline vector operator * ( vector const & u, vector const & v ) {
        return vec_madd( u.m128, v.m128, zero().m128 );
    }
//    friend inline vector operator / ( vector const & u, vector const & v ) {
//        return _mm_div_ps( u.m128, v.m128 );
//    }

    /*---------------------------------------------------------------------------
                             Math functions
    ---------------------------------------------------------------------------*/
    friend inline vector _rsqrt( vector const & u ) {
        return vec_rsqrt( u.m128 );
    }
    friend inline vector _normsq( vector const & u ) {
        __vector float sqr = vec_madd( u.m128, u.m128, zero().m128 );
        sqr = vec_add( sqr, vec_sld( sqr, sqr, 8 ) );
        sqr = vec_add( sqr, vec_sld( sqr, sqr, 4 ) );
        return sqr;
    }
    friend inline vector _dot( vector const & u, vector const & v ) {
        __vector float sqr = vec_madd( u.m128, v.m128, zero().m128 );
        sqr = vec_add( sqr, vec_sld( sqr, sqr, 8 ) );
        sqr = vec_add( sqr, vec_sld( sqr, sqr, 4 ) );
        return sqr;
    }
    friend inline vector _cross( vector const & u, vector const & v ) {
        vector x;
        for ( uint i = 0; i < 3; i++ ) x[i] = u[( i + 1U ) % 3U] * v[( i + 2U ) % 3U] - u[( i + 2U ) % 3U] * v[( i + 1U ) % 3U];
        x[3] = 0;
        return x;
    }
    /*---------------------------------------------------------------------------
                Math functions NOT TO BE USED IN FORCE KERNEL
    ---------------------------------------------------------------------------*/
    friend inline vector operator / ( vector const & u, float v ) {
        return vec_mul( u.m128, vec_re( vec_splats( v ) ) );
    }
    friend inline vector operator * ( vector const & u, float v ) {
        return vec_mul( u.m128, vec_splats( v ) );
    }
    friend inline vector operator * ( float v, vector const & u ) {
        return vec_mul( u.m128, vec_splats( v ) );
    }
    friend inline float dot( vector const & u, vector const & v ) {
        return _dot( u, v )[0];
    }
    friend inline float normsq( vector const & u ) {
        return _normsq( u )[0];
    }
    friend inline float norm( vector const & u ) {
        vector nsq { _normsq( u ) };
        return ( nsq * _rsqrt( nsq ) )[0];
    }
    friend inline vector cross( vector const & a, vector const & b ) {
        return _cross( a, b );
    }
    friend inline vector normalize( vector const & u ) {
        return u * _rsqrt( _normsq( u ) );
    }
};

#endif

#endif

}

#endif

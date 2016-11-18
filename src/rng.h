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
#ifndef OPENRBC_RNG_H_
#define OPENRBC_RNG_H_

#include <cstdlib>

namespace openrbc {

using namespace config;

const static real SQRT3 = 1.7320508075688773;
const static real TWOSQRT3 = 3.4641016151377545871;
const static real SQRT2 = 1.41421356237309504880;
const static real TWOSQRT2 = 2.8284271247461900976;

//alignas(config::cacheline)
struct  MT19937 {
public:
    MT19937() {
        //srand( time( NULL ) );
        srand( 0 );
        init( rand() );
    }

    MT19937( unsigned int seed ) {
        init( seed );
    }

    inline void init( unsigned int seed ) {
        state[0] = seed;
        for ( int i = 1; i < 624; ++i ) {
            state[i] = ( 1812433253 * ( state[i - 1] ^ ( state[i - 1] >> 30 ) ) + i );
        }

        make_random( 0U );
        make_random( real( 0 ) );
        saved = false;
    }

    // uniform 32bit unsigned integer
    inline unsigned int uint() {
        if ( pidata == piend ) make_random( 0U );
        return *pidata++;
    }

    // uniform random number between 0 and 1
    inline real u01() {
        if ( prdata == prend ) make_random( real( 0 ) );
        return *prdata++;
    }

    // uniform random number between -1 and 1
    inline real u11() {
        if ( prdata == prend ) make_random( real( 0 ) );
        return *prdata++ * constant::two - constant::one;
    }

    // Gaussian random number with zero mean and unit variance
    inline real gaussian() {
        if ( saved ) {
            saved = !saved;
            return gsave;
        } else {
            real x, y, s;
            do {
                x = u11();
                y = u11();
                s = x * x + y * y;
            } while ( s >= 1.0 );
            real r = sqrt( -2.0 * log( s ) / s );
            gsave = y * r;
            saved = !saved;
            return x * r;
        }
    }

private:
    unsigned int * pidata, *piend;
    unsigned int idata[624];
    real * prdata, *prend;
    real rdata[624];
    unsigned int state[624];

    // for the Polar method of Gaussian
    bool saved;
    real gsave;

    inline void make_random( unsigned int hint ) {
        for ( int i = 0; i < 624; i++ ) {
            int y = ( state[i] & 0x80000000 ) + ( state[( i + 1 ) % 624] & 0x7fffffff );
            state[i] = state[( i + 397 ) % 624] ^ ( y >> 1 );
            if ( y % 2 ) {
                state[i] ^= 0x9908b0df;
            }

            unsigned int z = state[i];
            z ^= ( z >> 11 );
            z ^= ( ( z << 7 ) & 0x9d2c5680 );
            z ^= ( ( z << 15 ) & 0xefc60000 );
            z ^= ( z >> 18 );
            idata[i] = z;
        }

        pidata = idata;
        piend = idata + 624;
    }

    inline void make_random( real hint ) {
        for ( int i = 0; i < 624; i++ ) {
            int y = ( state[i] & 0x80000000 ) + ( state[( i + 1 ) % 624] & 0x7fffffff );
            state[i] = state[( i + 397 ) % 624] ^ ( y >> 1 );
            if ( y % 2 ) {
                state[i] ^= 0x9908b0df;
            }

            unsigned int z = state[i];
            z ^= ( z >> 11 );
            z ^= ( ( z << 7 ) & 0x9d2c5680 );
            z ^= ( ( z << 15 ) & 0xefc60000 );
            z ^= ( z >> 18 );
            rdata[i] = z / real( 0xFFFFFFFF );
        }

        prdata = rdata;
        prend = rdata + 624;
    }
};

}

#endif

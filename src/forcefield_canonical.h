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
#ifndef OPENRBC_FORCE_FIELD_H_
#define OPENRBC_FORCE_FIELD_H_

#include <algorithm>
#include <cmath>
#include "config_static.h"
#include "rng.h"

namespace openrbc {

template<typename base_t>
constexpr base_t POW( base_t base, int expo ) {return ( expo != 0 ) ? base * POW( base, expo - 1 ) : 1;}
constexpr config::real cutsq( const config::real cut ) {return cut * cut;}
constexpr config::real rep( const config::real cut, const config::real req, const config::real epsilon ) {return epsilon / POW( cut - req, 8 );}
constexpr config::real att( const config::real cut, const config::real req, const config::real epsilon ) {return -2.0 * epsilon / POW( cut - req, 4 );}

struct ForceField {
    // 0 -- lipid
    // 1 -- mobile band III
    // 2 -- immobile band III
    // 3 -- glycophorin
    // 4 -- actin
    // 5 -- spectrin
    static constexpr int n_type = 6;
    static constexpr int spectrin_resolution = 19;
    // atom mass
    static constexpr config::real mass  [n_type] = {1, 4, 4, 1, 10, 10};
    static constexpr config::real radius[n_type] = {0.56125, 1.12375, 1.12375, 0.56125, 1.5, 0.5};

    // pairwise interaction coefficients
    // poly 4-8
    // lipid - lipid
    static constexpr config::real cutll = 2.6;
    static constexpr config::real reqll = 1.1225;
    static constexpr config::real epsilonll = 1.2;
    static constexpr config::real alphall = 1.55;
    // lipid - protein
    static constexpr config::real cutlp[n_type] = {cutll, 2.6, 2.6, 2.6};
    static constexpr config::real reqlp[n_type] = {reqll, 1.685, 1.685, 1.1225};
    static constexpr config::real epsilonlp[n_type] = {epsilonll, 1.4, 2.8, 2.8};
    static constexpr config::real alphalp[n_type] = {alphall, 5, 5, 5};
    // protein - protein
    static constexpr config::real cutpp[n_type * n_type] = {
        cutlp[0],     cutlp[1],                             cutlp[2],                                 cutlp[3],                                 0, 0,
        cutlp[1],     2.6,                                  2.6,                                      2.6,                                      0, 0,
        cutlp[2],     2.6,                                  2.6,                                      2.6,                                      0, 0,
        cutlp[3],     2.6,                                  2.6,                                      2.6
    };
    static constexpr config::real reqpp[n_type * n_type] = {
        reqlp[0],     reqlp[1],                             reqlp[2],                                 reqlp[3],                                 0, 0,
        reqlp[1],     2.245,                                2.245,                                    1.685,                                    0, 0,
        reqlp[2],     2.245,                                2.245,                                    1.685,                                    0, 0,
        reqlp[3],     1.685,                                1.685,                                    1.1225
    };
    static constexpr config::real epsilonpp[n_type * n_type] = {
        epsilonlp[0], epsilonlp[1],                         epsilonlp[2],                             epsilonlp[3],                             0, 0,
        epsilonlp[1], 1,                                    1,                                        1,                                        0, 0,
        epsilonlp[2], 1,                                    1,                                        1,                                        0, 0,
        epsilonlp[3], 1,                                    1,                                        1
    };

    // lj
    static constexpr config::real lj_epsilon[n_type * n_type] = {
        0, 0, 0, 0, 1, 1,
        0, 0, 0, 0, 1, 1,
        0, 0, 0, 0, 1, 1,
        0, 0, 0, 0, 0, 0,
        1, 1, 1, 0, 0, 0,
        1, 1, 1, 0, 0, 0
    };
    static constexpr config::real lj_sigma[n_type * n_type] = {
        0, 0,   0,   0, 1,   1,
        0, 0,   0,   0, 3.4, 3.4,
        0, 0,   0,   0, 3.4, 1,
        0, 0,   0,   0, 0,   0,
        1, 3.4, 3.4, 0, 3,   1.8,
        1, 3.4, 1,   0, 1.8, 1
    };
    static constexpr config::real lj_cut[n_type * n_type] = {
        0,      0,          0,          0, 1.1225,         1.1225,
        0,      0,          0,          0, 3.4 * 1.1225,     3.4 * 1.1225,
        0,      0,          0,          0, 3.4 * 1.1225,     1.1225,
        0,      0,          0,          0, 0,              0,
        1.1225, 3.4 * 1.1225, 3.4 * 1.1225, 0, 0,              0,
        1.1225, 3.4 * 1.1225, 1.1225,     0, 0,              0
    };

    // bond coefficients
    // 0: actin - glycophorin
    // 1: spectrin - immobile band 3
    // 2: actin - spectrin
    // 3: spectrin - spectrin
    static constexpr int n_bondtype = 4;
    //static constexpr config::real r0[n_bondtype] = {2.25, 1.1225, 2.25, 1.12};
    static constexpr config::real r0[n_bondtype] = {2.25, 1.1225, 2.25, 2.24};
    static constexpr config::real K[n_bondtype] = {57, 57, 57, 57};

    // auxiliary interaction coefficients
    // poly
    static constexpr config::real cutsqll = cutsq( cutll );
    static constexpr config::real repll = rep( cutll, reqll, epsilonll );
    static constexpr config::real attll = att( cutll, reqll, epsilonll );
    static constexpr config::real cutsqlp[n_type] = {cutsqll, cutsq( cutlp[1] ), cutsq( cutlp[2] ), cutsq( cutlp[3] )};
    static constexpr config::real replp[n_type] = {
        repll, rep( cutlp[1], reqlp[1], epsilonlp[1] ), rep( cutlp[2], reqlp[2], epsilonlp[2] ), rep( cutlp[3], reqlp[3], epsilonlp[3] )
    };
    static constexpr config::real attlp[n_type] = {
        attll, att( cutlp[1], reqlp[1], epsilonlp[1] ), att( cutlp[2], reqlp[2], epsilonlp[2] ), att( cutlp[3], reqlp[3], epsilonlp[3] )
    };
    static constexpr config::real cutsqpp[n_type * n_type] = {
        cutsqlp[0],   cutsqlp[1],                           cutsqlp[2],                               cutsqlp[3],                               0, 0,
        cutsqlp[1],   cutsq( cutpp[7] ),                      cutsq( cutpp[8] ),                          cutsq( cutpp[9] ),                          0, 0,
        cutsqlp[2],   cutsq( cutpp[13] ),                     cutsq( cutpp[14] ),                         cutsq( cutpp[15] ),                         0, 0,
        cutsqlp[3],   cutsq( cutpp[19] ),                     cutsq( cutpp[20] ),                         cutsq( cutpp[21] )
    };
    static constexpr config::real reppp[n_type * n_type] = {
        replp[0], replp[1],                                 replp[2],                                 replp[3],                                 0, 0,
        replp[1], rep( cutpp[7], reqpp[7], epsilonpp[7] ),    rep( cutpp[8], reqpp[8], epsilonpp[8] ),    rep( cutpp[9], reqpp[9], epsilonpp[9] ),    0, 0,
        replp[2], rep( cutpp[13], reqpp[13], epsilonpp[13] ), rep( cutpp[14], reqpp[14], epsilonpp[14] ), rep( cutpp[15], reqpp[15], epsilonpp[15] ), 0, 0,
        replp[3], rep( cutpp[19], reqpp[19], epsilonpp[19] ), rep( cutpp[20], reqpp[20], epsilonpp[20] ), rep( cutpp[21], reqpp[21], epsilonpp[21] )
    };
    // lj
    static config::real lj_cutsq[n_type * n_type];
    static config::real lj_lj1[n_type * n_type];
    static config::real lj_lj2[n_type * n_type];
    // max cut for tree::find_within
    static config::real cut_max[n_type];
    static bool initialized;

    static bool init() {
        for ( int i = 0; i < n_type * n_type; ++i ) {
            lj_cutsq[i] = lj_cut[i] * lj_cut[i];
            lj_lj1[i] = 48.0 * lj_epsilon[i] * std::pow( lj_sigma[i], 12 );
            lj_lj2[i] = 24.0 * lj_epsilon[i] * std::pow( lj_sigma[i], 6 );
        }
        for ( int i = 0; i < n_type; ++i ) {
            auto max1 = *std::max_element( cutpp + i * n_type, cutpp + i * n_type + n_type );
            auto max2 = *std::max_element( lj_cut + i * n_type, lj_cut + i * n_type + n_type );
            cut_max[i] = ( max1 > max2 ) ? max1 : max2;
        }
        return true;
    }
};

constexpr config::real ForceField::mass[ForceField::n_type];
constexpr config::real ForceField::radius[ForceField::n_type];
constexpr config::real ForceField::cutlp[ForceField::n_type];
constexpr config::real ForceField::reqlp[ForceField::n_type];
constexpr config::real ForceField::epsilonlp[ForceField::n_type];
constexpr config::real ForceField::alphalp[ForceField::n_type];
constexpr config::real ForceField::cutsqlp[ForceField::n_type];
constexpr config::real ForceField::replp[ForceField::n_type];
constexpr config::real ForceField::attlp[ForceField::n_type];
constexpr config::real ForceField::cutpp[ForceField::n_type * ForceField::n_type];
constexpr config::real ForceField::reqpp[ForceField::n_type * ForceField::n_type];
constexpr config::real ForceField::epsilonpp[ForceField::n_type * ForceField::n_type];
constexpr config::real ForceField::cutsqpp[ForceField::n_type * ForceField::n_type];
constexpr config::real ForceField::reppp[ForceField::n_type * ForceField::n_type];
constexpr config::real ForceField::lj_epsilon[ForceField::n_type * ForceField::n_type];
constexpr config::real ForceField::lj_sigma[ForceField::n_type * ForceField::n_type];
constexpr config::real ForceField::lj_cut[ForceField::n_type * ForceField::n_type];
constexpr config::real ForceField::r0[ForceField::n_bondtype];
constexpr config::real ForceField::K[ForceField::n_bondtype];
config::real ForceField::lj_cutsq[ForceField::n_type * ForceField::n_type];
config::real ForceField::lj_lj1[ForceField::n_type * ForceField::n_type];
config::real ForceField::lj_lj2[ForceField::n_type * ForceField::n_type];
config::real ForceField::cut_max[ForceField::n_type];
bool ForceField::initialized = ForceField::init();

}

#endif

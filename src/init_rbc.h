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
#ifndef OPENRBC_INIT_RBC_H_
#define OPENRBC_INIT_RBC_H_

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <set>
#include "container.h"
#include "forcefield.h"
#include "rng.h"
#include "kdtree.h"

namespace openrbc {

// Reference: Li & Lykotrafitis, Biophysical Journal, 2014
template<class C1, class C2>
double init_rbc( C1 & lipid, C2 & protein, RTParameter & param, const int vertex, const int step = 200, const int seed = 0 ) {
    using namespace config;
    using victor = vector<real, 3>;

    Service<Timers>::call()["init_rbc"].start();

    // Generate triangular mesh with python script
    if ( param.init == "vesicle" ) {
		char cmd[512];
		param.mesh = format( "vesicle-%08x", param.rng.uint() );
		sprintf( cmd, "python ../util/rbc_mesh.py -s %d -v %d -a %s -e %u", step, vertex, param.mesh.c_str(), param.rng.uint() );
		int error = system( cmd );
    }

    // Load the generated mesh
    printf("Loading mesh data from %s.{vert|bond|face}.txt\n", param.mesh.c_str() );
    std::ifstream fvert( param.mesh + ".vert.txt" );
    std::vector<victor> vert;
    while ( !fvert.eof() ) {
        double x, y, z;
        fvert >> x >> y >> z;
        if ( !fvert.eof() ) vert.emplace_back( x, y, z );
    }

    std::ifstream fbond( param.mesh + ".bond.txt" );
    std::vector<vector<int, 2> > bond;
    while ( !fbond.eof() ) {
        int i, j;
        fbond >> i >> j;
        if ( !fbond.eof() ) bond.emplace_back( i, j );
    }

    std::ifstream fface( param.mesh + ".face.txt" );
    std::vector<vector<int, 3> > face;
    while ( !fface.eof() ) {
        int i, j, k;
        fface >> i >> j >> k;
        if ( !fface.eof() ) face.emplace_back( i, j, k );
    }

    // Determine scaling/mapping by assuming that
    // the average bond length should be 80 nm
    // 1 length unit in simulation = 5 / 2^1/6 = 4.4545 nm
    double bond_length_avg = 0;
    for ( auto & b : bond ) bond_length_avg += norm( vert[b[0]] - vert[b[1]] );
    bond_length_avg /= bond.size();
    double bond_length_nm = 80.0;
    double length_unit_nm = 4.4545;
    double bond_length_reduced = bond_length_nm / length_unit_nm;
    double unit2reduced = bond_length_reduced / bond_length_avg;
    double layer_dist = ForceField::r0[0]; // actin-glycophorin bond length

    // Scale vertices coordinate from unit sphere to simulation units
    for ( auto & v : vert ) v *= unit2reduced;

    // Compute vertex- and edge- averaged normal
    std::vector<victor > normal_face( face.size(), 0 ),
        normal_vert( vert.size(), 0 ),
        normal_bond( bond.size(), 0 );
    std::vector<std::set<int> >    vert2face( vert.size() ),
        bond2face( bond.size() );
    // Calculate face normal
    for ( std::size_t i = 0; i < face.size(); i++ ) {
        normal_face[i] = cross( vert[face[i][0]] - vert[face[i][1]], vert[face[i][2]] - vert[face[i][1]] );
        normal_face[i] *= norm_inv( normal_face[i] );
        //if ( dot( normal_face[i], 0.333 * ( vert[face[i][0]] + vert[face[i][1]] + vert[face[i][2]] ) ) < 0 ) normal_face[i] *= -1.0; // this 'smart' flipping caused more trouble than good
        for ( int d = 0; d < 3; d++ ) vert2face[ face[i][d] ].insert( i );
    }
    // Find neighboring faces for each bond
    for ( std::size_t i = 0; i < bond.size(); i++ ) {
        auto & f1 = vert2face[bond[i][0]];
        auto & f2 = vert2face[bond[i][1]];
        std::set_intersection( f1.begin(), f1.end(), f2.begin(), f2.end(), std::inserter( bond2face[i], bond2face[i].begin() ) );
        assert( bond2face[i].size() == 2 );
    }
    // Compute averaged vertex normal from all neighboring faces
    for ( std::size_t i = 0; i < vert.size(); i++ ) {
        for ( auto & f : vert2face[i] ) {
            normal_vert[i] += normal_face[ f ];
        }
        normal_vert[i] *= norm_inv( normal_vert[i] );
    }
    // Compute averaged bond normal from the two neighboring faces
    for ( std::size_t i = 0; i < bond.size(); i++ ) {
        for ( auto & f : bond2face[i] ) {
            normal_bond[i] += normal_face[ f ];
        }
        normal_bond[i] *= norm_inv( normal_bond[i] );
    }

    int base;

    printf( "Generating glycophorin\n" );
    std::map<int, int> vert2gly;
    protein.resize( 0 );
    protein.resize( vert.size() );
    for ( std::size_t i = 0; i < vert.size(); i++ ) {
        protein.x   [i][0] = vert[i][0];
        protein.x   [i][1] = vert[i][1];
        protein.x   [i][2] = vert[i][2];
        protein.n   [i][0] = normal_vert[i][0];
        protein.n   [i][1] = normal_vert[i][1];
        protein.n   [i][2] = normal_vert[i][2];
        protein.v   [i] = 0;
        protein.o   [i] = 0;
        protein.type[i] = 3;
        protein.tag [i] = i + 1;
        vert2gly     [i] = i;
    }

    printf( "Generating actin\n" );
    base = protein.size();
    std::map<int, int> vert2actin;
    protein.resize( protein.size() + vert.size() );
    for ( std::size_t i = 0; i < vert.size(); i++ ) {
        protein.x   [base + i][0] = vert[i][0] - layer_dist * normal_vert[i][0];
        protein.x   [base + i][1] = vert[i][1] - layer_dist * normal_vert[i][1];
        protein.x   [base + i][2] = vert[i][2] - layer_dist * normal_vert[i][2];
        protein.n   [base + i][0] = normal_vert[i][0];
        protein.n   [base + i][1] = normal_vert[i][1];
        protein.n   [base + i][2] = normal_vert[i][2];
        protein.v   [base + i] = 0;
        protein.o   [base + i] = 0;
        protein.type[base + i] = 4;
        protein.tag [base + i] = base + i + 1;
        vert2actin   [       i] = base + i;
        protein.bonds.emplace_back( 0, protein.tag[ vert2gly[i] ], protein.tag [base + i] );
    }

    printf( "Generating immobile band-III\n" );
    base = protein.size();
    std::map<int, int> bond2ib3;
    protein.resize( protein.size() + bond.size() );
    for ( std::size_t i = 0; i < bond.size(); i++ ) {
        auto mid = 0.5 * ( vert[ bond[i][0] ] + vert[ bond[i][1] ] );
        protein.x   [base + i][0] = mid[0];
        protein.x   [base + i][1] = mid[1];
        protein.x   [base + i][2] = mid[2];
        protein.n   [base + i][0] = normal_bond[i][0];
        protein.n   [base + i][1] = normal_bond[i][1];
        protein.n   [base + i][2] = normal_bond[i][2];
        protein.v   [base + i] = 0;
        protein.o   [base + i] = 0;
        protein.type[base + i] = 2;
        protein.tag [base + i] = base + i + 1;
        bond2ib3    [       i] = base + i;
    }

    /*
    printf("Generating mobile band-III\n");
    base = protein.size();
    protein.resize( protein.size() + face.size() * 3 );
    for ( std::size_t i = 0; i < face.size(); i++ ) {
        for ( int j = 0; j < 3; j++ ) {
            // Uniformly sampling triangle
            // Ref: Osada, et al. ACM Transactions on Graphics 2002
            double r1 = param.rng.u01();
            double r2 = param.rng.u01();
            auto rt = vert[ face[i][0] ] * ( 1.0 - std::sqrt( r1 ) ) +
                      vert[ face[i][1] ] * ( std::sqrt( r1 ) * ( 1.0 - r2 ) ) +
                      vert[ face[i][2] ] * ( std::sqrt( r1 ) * r2 );
            protein.x   [base + i * 3 + j] = rt;
            protein.v   [base + i * 3 + j] = 0;
            protein.n   [base + i * 3 + j] = normal_face[i];
            protein.type[base + i * 3 + j] = 1;
            protein.tag [base + i * 3 + j] = base + i * 3 + j + 1;
        }
    }
    */

    printf( "Generating spectrin network\n" );
    base = protein.size();
    const static int spc_res = ForceField::spectrin_resolution;
    protein.resize( protein.size() + bond.size() * spc_res );
    for ( std::size_t i = 0; i < bond.size(); i++ ) {
        auto beg = vert[bond[i][0]] - layer_dist * normal_vert[bond[i][0]];
        auto end = vert[bond[i][1]] - layer_dist * normal_vert[bond[i][1]];
        double span = norm( beg - end );
        auto inc = 1.f / ( spc_res + 1 ) * ( end - beg );
        beg += inc / norm( inc ) * 1.0;
        end -= inc / norm( inc ) * 1.0;
        span = norm( beg - end );
        inc = 1.f / ( spc_res + 1 ) * ( end - beg );
        for ( int j = 0; j < spc_res; j++ ) {
            auto x = beg + inc * ( j + 1 );
            victor noise( param.rng.u11(), param.rng.u11(), param.rng.u11() );
            protein.x   [base + i * spc_res + j][0] = x[0] + 0.05 * noise[0];
            protein.x   [base + i * spc_res + j][1] = x[1] + 0.05 * noise[1];
            protein.x   [base + i * spc_res + j][2] = x[2] + 0.05 * noise[2];
            protein.n   [base + i * spc_res + j][0] = x[0];
            protein.n   [base + i * spc_res + j][1] = x[1];
            protein.n   [base + i * spc_res + j][2] = x[2];
            protein.v   [base + i * spc_res + j] = 0;
            protein.o   [base + i * spc_res + j] = 0;
            protein.type[base + i * spc_res + j] = 5;
            protein.tag [base + i * spc_res + j] = base + i * spc_res + j + 1;
            if ( j > 0 ) protein.bonds.emplace_back( 3, protein.tag[base + i * spc_res + j - 1], protein.tag[base + i * spc_res + j] );
            switch ( j ) {
            case         0: protein.bonds.emplace_back( 2, protein.tag[ vert2actin[bond[i][0]] ], protein.tag[base + i * spc_res + j] ); break;
            case spc_res/2: protein.bonds.emplace_back( 1, protein.tag[ bond2ib3[i]            ], protein.tag[base + i * spc_res + j] ); break;
            case spc_res-1: protein.bonds.emplace_back( 2, protein.tag[ vert2actin[bond[i][1]] ], protein.tag[base + i * spc_res + j] ); break;
            }
        }
    }

    printf( "Generating lipids\n" );
    KDTree<real, 3> tree;
    AlignedArray<vector<real, 3> > points( protein.size() );
    for ( int i = 0; i < protein.size(); ++i ) points[i] = { protein.x[i][0], protein.x[i][1], protein.x[i][2] };
    tree.build( points );
    lipid.resize( 0 );
    //double area_per_lipid = 2770.0 / 140 / 4.4545 / 4.4545;
    //printf("lipid number density %lf\n", 1.0 / area_per_lipid);
    const static double min_dist = 0.55;
    const static double min_dist_lp = ForceField::reqlp[1];

    double total_area = 0;

    for ( std::size_t i = 0; i < face.size(); i++ ) {
        if ( ( i + 1 ) % ( face.size() / 10 ) == 0 ) {
            printf( "%lu...", i );
            fflush( stdout );
        }
        auto v1 = vert[ face[i][0] ];
        auto v2 = vert[ face[i][1] ];
        auto v3 = vert[ face[i][2] ];
        auto area = 0.5 * norm( cross( v1 - v2, v3 - v2 ) );
        total_area += area;
        int n = area * param.rho;
        int base = lipid.size();
        lipid.resize( lipid.size() + n );
        std::vector<victor > pts;
        for ( int j = 0; j < n; j++ ) {
            victor r;
            bool clutter;
            do {
                double r1 = param.rng.u01();
                double r2 = param.rng.u01();
                r = v1 * ( 1.0 - std::sqrt( r1 ) ) +
                    v2 * ( std::sqrt( r1 ) * ( 1.0 - r2 ) ) +
                    v3 * ( std::sqrt( r1 ) * r2 );
                clutter = false;
                // Check clutter with protein
                double dp;
                std::tie( dp, std::ignore ) = tree.find_nearest( r, std::make_tuple( norm( r - v1 ), face[i][0] )  );
                if ( dp < min_dist_lp ) {
                    clutter = true;
                    continue;
                }
                // Check clutter with existing lipids
                for ( std::size_t k = 0; k < pts.size(); k++ ) {
                    if ( normsq( pts[k] - r ) < min_dist * min_dist ) {
                        clutter = true;
                        break;
                    }
                }
            } while ( clutter );
            pts.push_back( r );
            lipid.x[base + j][0] = r[0];
            lipid.x[base + j][1] = r[1];
            lipid.x[base + j][2] = r[2];
            lipid.n[base + j][0] = normal_face[i][0];
            lipid.n[base + j][1] = normal_face[i][1];
            lipid.n[base + j][2] = normal_face[i][2];
            lipid.v[base + j] = 0;
            lipid.o[base + j] = 0;
        }
    }
    lipid.tag.set_base( protein.size() + 1 );
    printf( "Lipid sphere total area: %lf\n", total_area );
    printf( "Done\n" );

    for ( int d = 3; d < vect::d; ++d ) {
        #pragma omp parallel for
        for ( int i = 0; i < lipid.size(); ++i ) {
            lipid.x[i][d] = 0;
            lipid.v[i][d] = 0;
            lipid.f[i][d] = 0;
            lipid.n[i][d] = 0;
            lipid.o[i][d] = 0;
            lipid.t[i][d] = 0;
        }
        #pragma omp parallel for
        for ( int i = 0; i < protein.size(); ++i ) {
            protein.x[i][d] = 0;
            protein.v[i][d] = 0;
            protein.f[i][d] = 0;
            protein.n[i][d] = 0;
            protein.o[i][d] = 0;
            protein.t[i][d] = 0;
        }
    }

    lipid  .tag2idx.build_map( lipid   );
    protein.tag2idx.build_map( protein );

    Service<Timers>::call()["init_rbc"].stop();

    return total_area;
}

}

#endif

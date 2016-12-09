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
#ifndef OPENRBC_CONFIG_RT_H_
#define OPENRBC_CONFIG_RT_H_

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <omp.h>
#include <cmath>
#include <cstdlib>
#include <vector>
#include "config_static.h"
#include "rng.h"
#include "util_misc.h"

namespace openrbc {

namespace DumpField {
const static unsigned int position = 1;
const static unsigned int rotation = 2;
const static unsigned int voronoi  = 4;
const static unsigned int velocity = 8;
const static unsigned int force    = 16;
};

struct RTParameter {
    using real = config::real;

    std::string init = "vesicle"; // initialization mode
    std::string mesh;             // mesh file prefix

    // system configuration
    double box[3][2];
    double bsize[3];
    int voronoi_cell_size = 14;
    real rho = 1.05;
    real kBT = 0.22;
    real eta = 0.01;
    real stray_tolerance = 2.5;

    // time marching
    double t_total      = 10;   // target total simulation time
    double dt           = 1E-2; // time step size
    int opt_nstep       = 1000; // number of optimization time steps
    double dr_opt       = 5E-2; // maximal displacement in enery minimization
    double dn_opt       = 5E-2; // maximal angular movement in enery minimization
    int freq_voronoi    = 2;
    int freq_sort_ctrd  = 24;
    int freq_sort_bond  = 120;
    int freq_cleanup    = 60;
    int nstep           = 0;  // current time step
    int rseed           = 0xBAD5EED;

    // I/O
    int freq_display    = 100;
    int freq_dump       = 10000;
    int dump_field      = DumpField::position + DumpField::rotation + DumpField::voronoi;
    std::string file_topo = "cell.data";
    std::string file_traj = "cell.orbc";

    // Nose-Hoover thermostat
    real Q; // ~ number of particles / 100
    real zeta = 0;
    void set_Q( std::size_t n ) { Q = n * 0.01; }

    // RNG
    mutable MT19937 rng;
    std::vector<MT19937> prng;

    RTParameter( int argc, char ** argv ) {
        auto atos = []( const char f[] ) { return std::string( f ); };
        parse( init,              "i",  "init",                atos, argc, argv, "Model initialization mode, must be 'lipid' or 'vesicle' or 'trimesh'" );
        parse( mesh,              "m",  "mesh",                atos, argc, argv, "Prefix for the 3 mesh files, e.g. [mesh].bond.txt, [mesh].face.txt, [mesh].vert.txt" );
        parse( t_total,           "t",  "total-time",          atof, argc, argv, "Total time in simulation units" );
        parse( dt,                "k",  "dt",                  atof, argc, argv, "Simulation time step size" );
        parse( eta,               "e",  "eta",                 atof, argc, argv, "Solvent viscosity for Langevin integrator" );
        parse( stray_tolerance,   "",   "stray-tolerance",     atof, argc, argv, "A lipid can stray at most stray-tolerance times than the median particles before being removed" );
        parse( opt_nstep,         "E",  "opt-nstep",           atoi, argc, argv, "Number of optimization steps" );
        parse( dr_opt,            "",   "opt-max-move",        atof, argc, argv, "Maximum displacement of a particle during each optimization step" );
        parse( dn_opt,            "",   "opt-max-rotate",      atof, argc, argv, "Maximum rotation of a particle during each optimization step" );
        parse( freq_voronoi,      "V",  "voronoi-update-freq", atoi, argc, argv, "Frequency (#steps) of Voronoi cell update" );
        parse( freq_sort_ctrd,    "",   "centroid-reorder-freq", atoi, argc, argv, "Frequency (#steps) of Voronoi cell centroid reorder" );
        parse( freq_sort_bond,    "",   "bond-reorder-freq",   atoi, argc, argv, "Frequency (#steps) of bond locality maintenance" );
        parse( freq_display,      "D",  "display-frequency",   atoi, argc, argv, "Frequency (#steps) of screen info display" );
        parse( freq_dump,         "d",  "dump-frequency",      atoi, argc, argv, "Frequency (#steps) of trajectory file output" );
        parse( freq_sort_bond,    "",   "bond-sort-frequency", atoi, argc, argv, "Frequency (#steps) to reorder bonds to improve data locality" );
        parse( freq_cleanup,      "",   "cleanup-frequency",   atoi, argc, argv, "Frequency (#steps) to remove stray lipid particles" );
        parse( dump_field,        "",   "dump-field",          atoi, argc, argv, "Fields to dump, 1: position, 2: rotation, 4: voronoi affiliation, 8: velocity, 16: force" );
        parse( Q,                 "Q",  "nose-hoover-Q",       atof, argc, argv, "Nose-Hoover relaxation parameter" );
        parse( rho,               "r",  "rho",                 atof, argc, argv, "Lipid number density" );
        parse( kBT,               "T",  "kBT",                 atof, argc, argv, "Thermal noise level" );
        parse( voronoi_cell_size, "",   "voronoi-cell-size",   atoi, argc, argv, "Average size (#lipids) of Voronoi cells" );
        parse( file_topo,         "o",  "topology-file",       atos, argc, argv, "Filename for saving topology" );
        parse( file_traj,         "x",  "trajectory-file",     atos, argc, argv, "Filename for saving trajectory" );
        parse( rseed,             "r" , "random-seed",         atoi, argc, argv, "Seed for random number generation" );
        rng.init( rseed );
        prng.resize( omp_get_max_threads() );
        for ( auto & r : prng ) r.init( rng.uint() );
        for ( int d = 0; d < 3; ++d ) {
            box[d][0] = -1000;
            box[d][1] =  1000;
            bsize[d] = box[d][1] - box[d][0];
        }
        rt_assert( init == "lipid" || init == "vesicle" || init == "trimesh", "Initialization mode must be 'lipid', 'vesicle', or 'trimesh'" );
        if ( init == "trimesh" ) rt_assert( mesh.size(), "Initialization mode 'trimesh' specified, but mesh files are not given" );
        rt_assert( freq_sort_ctrd % freq_voronoi == 0, "Voronoi centroid reorder frequency must be a multiple of voronoi update frequency" );
        rt_assert( freq_sort_bond % freq_voronoi == 0, "bond reorder frequency must be a multiple of voronoi update frequency" );
        rt_assert( freq_cleanup   % freq_voronoi == 0, "stray lipid cleanup frequency must be a multiple of voronoi update frequency" );
    }

    template<typename T, class CONVERT>
    static void parse( T & var, std::string short_name, std::string long_name, CONVERT const & cvt, int argc, char ** argv, const char explanation[] ) {
        bool found = false;
        for ( int i = 1; i < argc; ++i ) {
            if ( ( "-" + short_name ) == argv[i] && i < argc - 1 ) {
                var = cvt( argv[++i] );
                std::cout << "From command line: " << std::left << std::setw( 25 ) << short_name << " = " << std::setw( 16 ) << var << ' ' << explanation << std::endl;
                found = true;
            } else if ( ( "--" + long_name ) == argv[i] && i < argc - 1 ) {
                var = cvt( argv[++i] );
                std::cout << "From command line: " << std::left << std::setw( 25 ) << long_name << " = " << std::setw( 16 ) << var << ' ' << explanation << std::endl;
                found = true;
            } else if ( !strcmp( argv[i], "-h" ) || !strcmp( argv[i], "--help" ) ) {
                std::cout << std::left << std::setw( 6 ) << ( short_name.size() ? ( "-" + short_name ) : "" )
                          << std::left << std::setw( 25 ) << "--" + long_name
                          << "   " << explanation << std::endl;
                return;
            }
        }
        if ( !found ) std::cout << "Using default arg: " << std::left << std::setw( 25 ) << long_name << " = " << std::setw( 16 ) << var << ' ' << explanation << std::endl;
    }
};

}

#endif

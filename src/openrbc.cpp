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
#include <iostream>
#include "compute_bonded.h"
#include "compute_pairwise.h"
#include "compute_pairwise_simd.h"
#include "compute_pairwise_fused.h"
#include "compute_temperature.h"
#include "config_static.h"
#include "cleanup.h"
#include "runtime_parameter.h"
#include "container.h"
#include "display.h"
#include "init_random.h"
#include "init_rbc.h"
#include "citation.h"
#include "integrate_nh.h"
#include "integrate_langevin.h"
#include "forcefield.h"
#include "rdf.h"
#include "remove_bonds.h"
#include "rng.h"
#include "timer.h"
#include "topology.h"
#include "trajectory.h"
#include "voronoi.h"
#include "constrain_volume.h"
#include "zero_bulk_velocity.h"
#include "assign_temperature.h"

int main( int argc, char ** argv ) {
    using namespace openrbc;
    using namespace openrbc::config;

    omp_set_nested( 1 );

    RTParameter param( argc, argv );

    // Initialize system
    std::cout << "Initializing system ..." << std::flush;
    LipidContainer lipid( "lipid" );
    ProteContainer protein( "protein" );
    if ( param.init == "lipid" ) {
    	init_random_sphere( lipid, param, 100 );
    	save_topology( std::ofstream( param.file_topo ), param, protein, lipid );
    } else if ( param.init == "vesicle" ) {
    	init_rbc( lipid, protein, param, 500 );
        //remove_bonds(protein, UnaryPredicate());
    	save_topology( std::ofstream( param.file_topo ), param, protein, lipid );
    } else if ( param.init == "trimesh" ) {
    	init_rbc( lipid, protein, param, 0 );
    	save_topology( std::ofstream( param.file_topo ), param, protein, lipid );
    }
    std::cout << "Done." << std::endl;

    std::cout << "Initializing Voronoi cells " << std::flush;
    VoronoiDiagram voronoi( std::max<std::size_t>( 1, lipid.size() / param.voronoi_cell_size ) );
    VCellList cell_lipid, cell_protein;
    voronoi.init( lipid, cell_lipid, param, 64 );
    cell_lipid.update_particle_affiliation( lipid );
    cell_protein.update( protein, voronoi, param );
    cell_protein.update_particle_affiliation( protein );

    std::cout << "Done." << std::endl;
    std::ofstream traj( param.file_traj );
    save_frame( traj, lipid, protein, cell_lipid, cell_protein, param );

    std::cout << "Initialization complete." << std::endl;
    Service<Timers>::call().report( true );

    std::cout << "Opt..." << std::endl;
    Service<Timers>::call()["+optimization"].start();

    // energy minimization
    //while ( param.dt_opt * nopt++ < param.t_opt ) {
    for ( int nopt = 0; nopt < param.opt_nstep; ++nopt ) {
        voronoi.update( lipid, cell_lipid, param );
        cell_lipid.update( lipid, voronoi, param );
        cell_protein.update( protein, voronoi, param );

        // force evaluation
        integrate( clear_force(), lipid, protein );
        #ifdef EXPLICIT_SIMD
        compute_pairwise_simd( voronoi, lipid, protein, cell_lipid, cell_protein );
        compute_bonded_simd( protein );
        #else
        #ifdef FUSED_PAIRWISE
        compute_pairwise_fused( voronoi, lipid, protein, cell_lipid, cell_protein );
        #else
        compute_pairwise( lipid, cell_lipid, voronoi );
        compute_pairwise_lp( lipid, protein, voronoi, cell_lipid, cell_protein );
        compute_pairwise_pp( protein, voronoi, cell_protein );
        #endif
        compute_bonded( protein );
        #endif

        // time integration again
        integrate( post_torque(), lipid, protein );

        Service<Timers>::call()["OptIntegration"].start();

        #pragma omp parallel for
        for ( std::size_t i = 0; i < lipid.size(); ++i ) {
            auto dx = lipid.f[i] / ForceField::mass[lipid.type[i]];
            auto dn = cross( lipid.t[i], lipid.n[i] );
            auto dt = ( norm( dx ) > param.dr_opt || norm( dn ) > param.dn_opt ) ? std::min( param.dr_opt / norm( dx ), param.dn_opt / norm( dn ) ) : param.dt;
            lipid.x[i] += lipid.f[i] * ( dt / ForceField::mass[lipid.type[i]] );
            lipid.n[i] += cross( lipid.t[i], lipid.n[i] ) * dt;
            lipid.n[i]  = normalize( lipid.n[i] );
        }
        #pragma omp parallel for
        for ( std::size_t i = 0; i < protein.size(); ++i ) {
            auto dx = protein.f[i] / ForceField::mass[protein.type[i]];
            auto dn = cross( protein.t[i], protein.n[i] );
            auto dt = ( norm( dx ) > param.dr_opt || norm( dn ) > param.dn_opt ) ? std::min( param.dr_opt / norm( dx ), param.dn_opt / norm( dn ) ) : param.dt;
            protein.x[i] += protein.f[i] * ( dt / ForceField::mass[protein.type[i]] );
            protein.n[i] += cross( protein.t[i], protein.n[i] ) * dt;
            protein.n[i]  = normalize( protein.n[i] );
        }
        Service<Timers>::call()["OptIntegration"].stop();

        integrate( bounce_back( param ), lipid, protein );

        // I/O
        if ( ( nopt + 1 ) % param.freq_dump == 0 ) {
            cell_lipid.update_particle_affiliation( lipid );
            cell_protein.update_particle_affiliation( protein );
            save_frame( traj, lipid, protein, cell_lipid, cell_protein, param );
        }
        if ( ( nopt + 1 ) % param.freq_display == 0 ) {
            display( std::cout, nopt + 1, compute_temperature( lipid, protein, param ),
                     omp_get_wtime() - Service<Timers>::call()["+optimization"].get_start_time() );
        }
    }

    Service<Timers>::call()["+optimization"].stop();
    Service<Timers>::call().report( true );

    // run
    std::cout << "Run ... " << std::endl;
    Service<Timers>::call()["+main-loop"].start();

    voronoi.update( lipid, cell_lipid, param );
    cell_lipid.update( lipid, voronoi, param );
    cell_protein.update( protein, voronoi, param );

    integrate( clear_force(), lipid, protein );
    integrate( assign_temperature(param), lipid, protein );

    #ifdef EXPLICIT_SIMD
    compute_pairwise_simd( voronoi, lipid, protein, cell_lipid, cell_protein );
    compute_bonded_simd( protein );
    #else
    #ifdef FUSED_PAIRWISE
    compute_pairwise_fused( voronoi, lipid, protein, cell_lipid, cell_protein );
    #else
    compute_pairwise( lipid, cell_lipid, voronoi );
    compute_pairwise_lp( lipid, protein, voronoi, cell_lipid, cell_protein );
    compute_pairwise_pp( protein, voronoi, cell_protein );
    #endif
    compute_bonded( protein );
    #endif

    #ifdef LANGEVIN
    integrate( verlet_langevin( param ), lipid, protein );
    #else
    #ifdef FUSED_INTEGRATOR
    integrate( post_toque_final_update( param ), lipid, protein );
    #else
    integrate( post_torque(), lipid, protein );
    integrate( verlet_nh_final( param ), lipid, protein );
    integrate( verlet_nh_update( param ), lipid, protein );
    #endif
    #endif

    param.nstep = 0.0;
    while ( param.nstep * param.dt < param.t_total ) {
        // time integration
        #ifndef LANGEVIN
        #ifdef FUSED_INTEGRATOR
        integrate( verlet_initial_bounce_clearforce_update( param ), lipid, protein );
        #else
        integrate( verlet_nh_initial( param ), lipid, protein );
        integrate( bounce_back( param ), lipid, protein );
        #endif
        #endif

        if ( param.nstep % param.freq_voronoi == 0 ) {
            if ( param.nstep % param.freq_cleanup == 0) delete_lipid(lipid, voronoi, cell_lipid, param);
            voronoi.update( lipid, cell_lipid, param );
            cell_lipid.update( lipid, voronoi, param );
            cell_protein.update( protein, voronoi, param );
        }

        // force evaluation
        #ifndef LANGEVIN
        #ifndef FUSED_INTEGRATOR
        integrate( clear_force(), lipid, protein );
        #endif
        #endif

        #ifdef EXPLICIT_SIMD
        compute_pairwise_simd( voronoi, lipid, protein, cell_lipid, cell_protein );
        compute_bonded_simd( protein );
        #else
        #ifdef FUSED_PAIRWISE
        compute_pairwise_fused( voronoi, lipid, protein, cell_lipid, cell_protein );
        #else
        compute_pairwise( lipid, cell_lipid, voronoi );
        compute_pairwise_lp( lipid, protein, voronoi, cell_lipid, cell_protein );
        compute_pairwise_pp( protein, voronoi, cell_protein );
        #endif
        compute_bonded( protein );
        #endif

        //integrate( zero_bulk_velocity( param ), lipid );
        //constrain_volume( lipid, protein, voronoi, cell_lipid, cell_protein, param, 3.15, 0.05 );

        // time integration again
        #ifdef LANGEVIN
        integrate( verlet_langevin( param ), lipid, protein );
        #else
        #ifdef FUSED_INTEGRATOR
        integrate( post_toque_final_update( param ), lipid, protein );
        #else
        integrate( post_torque(), lipid, protein );
        integrate( verlet_nh_final( param ), lipid, protein );
        integrate( verlet_nh_update( param ), lipid, protein );
        #endif
        #endif

        ++param.nstep;

        // I/O
        if ( param.nstep % param.freq_dump == 0 ) {
            cell_lipid.update_particle_affiliation( lipid );
            cell_protein.update_particle_affiliation( protein );
            save_frame( traj, lipid, protein, cell_lipid, cell_protein, param );
        }
        if ( param.nstep % param.freq_display == 0 ) {
            display( std::cout, param.nstep * param.dt, compute_temperature( lipid, protein, param ),
                     omp_get_wtime() - Service<Timers>::call()["+main-loop"].get_start_time() );
        }
    }

    // finalization
    std::stringstream msg;
    msg << param.nstep << " steps * " << lipid.size() + protein.size() << " particles on " << omp_get_max_threads() << " threads.";
    display_timing( std::cout, Service<Timers>::call()["+main-loop"].stop(), msg.str().c_str() );

    return 0;
}

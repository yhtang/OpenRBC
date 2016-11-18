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
#ifndef OPENRBC_INIT_FILE_H_
#define OPENRBC_INIT_FILE_H_

#include <fstream>
#include <iostream>
#include "config.h"
#include "container.h"
#include "forcefield.h"

namespace openrbc {

// place particles from input file
void init_file( Container & model, RTParameter const & param, const char * const filename ) {
    Service<Timers>::call()["init_file"].start();

    std::ifstream fin( filename );

    int num_particles;
    fin >> num_particles;
    double scale[3];
    fin >> scale[0] >> scale[1] >> scale[2];

    model.resize( num_particles );

    for ( int i = 0; i < num_particles; ++i ) {
        fin >> model.tag[i];
        fin >> model.type[i];
        double tmp;
        for ( int d = 0; d < 3; ++d ) {
            fin >> tmp;
            model.x[i][d] = tmp * scale[d];
        }

        #if CASE == 1
        for ( int d = 0; d < 3; ++d ) fin >> tmp;
        #endif
        for ( int d = 0; d < 3; ++d ) fin >> model.n[i][d];
        #if CASE == 1
        for ( int d = 0; d < 3; ++d ) fin >> tmp;
        #endif
    }
    fin.close();

    model.v    .assign( num_particles, real( 0 ) );
    model.omega.assign( num_particles, real( 0 ) );

    Service<Timers>::call()["init_file"].stop();

    model.tag2idx.build_map( model );
}

void init_file_bond( Container & model, RTParameter const & param, const char * const filename ) {
    Service<Timers>::call()["init_file_bond"].start();

    std::ifstream fin( filename );

    int n_bonds;
    fin >> n_bonds;

    model.bonds.resize( n_bonds );
    for ( int i = 0; i < n_bonds; ++i )
        fin >> model.bonds[i].type >> model.bonds[i].i >> model.bonds[i].j;
    fin.close();

    Service<Timers>::call()["init_file_bond"].stop();
}

}

#endif

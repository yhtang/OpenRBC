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
#ifndef TOPOLOGY_H_
#define TOPOLOGY_H_

#include "container.h"
#include "forcefield.h"

namespace openrbc {

// write both lipid and protein
inline void save_topology( std::ostream && file, RTParameter const & param, const ProteContainer & prote, const LipidContainer & lipid ) {
    file << "WholeCell" << std::endl << std::endl;
    file << prote.size() + lipid.size() << " " <<  "atoms" << std::endl;
    file << prote.bonds.size() << " " << "bonds" << std::endl; // what is num_bonds?
    file << std::endl << std::endl;

    // dim
    file << param.box[0][0] << " " << param.box[0][1] << " " << "xlo xhi" << std::endl;
    file << param.box[1][0] << " " << param.box[1][1] << " " << "ylo yhi" << std::endl;
    file << param.box[2][0] << " " << param.box[2][1] << " " << "zlo zhi" << std::endl;

    file << std::endl          << std::endl;
    file << ForceField::n_type << " " << "atom types" << std::endl;
    file << "1 "               << " " << "bond types" << std::endl;

    file << std::endl   << std::endl;
    file << "Atoms"     << std::endl;

    for ( std::size_t i = 0; i < prote.size(); i++ ) {
        // (tag) (mol id) (type) pos(x y z)
        file << prote.tag[i] << " " << "1 " << prote.type[i] + 1 << " ";
        for ( int d = 0; d < 3 ; d++ ) file << prote.x[i][d] << " ";
        file << std::endl;
    }
    for ( std::size_t i = 0; i < lipid.size(); i++ ) {
        file << lipid.tag[i] << " " << "1 " << lipid.type[i] + 1 << " ";
        for ( int d = 0; d < 3 ; d++ ) file << lipid.x[i][d] << " ";
        file << std::endl;
    }

    file << std::endl << std::endl;
    file << "Bonds" << std::endl << std::endl;

    for ( std::size_t i = 0; i < prote.bonds.size(); i++ ) {
        // (serial #) (type) (id first bond) (id second bond)
        file << ( i + 1 ) << " " << prote.bonds[i].type  << " ";
        file << prote.bonds[i].i << " " << prote.bonds[i].j << std::endl;
    }
}

// write only lipid
inline void save_topology( const LipidContainer & model, RTParameter const & param, std::ostream && file ) {
    file << "Lipid" << std::endl << std::endl;
    file << model.size() << " " <<  "atoms" << std::endl;
    file << std::endl << std::endl;

    // dim
    file << param.box[0][0] << " " << param.box[0][1] << " " << "xlo xhi" << std::endl;
    file << param.box[1][0] << " " << param.box[1][1] << " " << "ylo yhi" << std::endl;
    file << param.box[2][0] << " " << param.box[2][1] << " " << "zlo zhi" << std::endl;

    file << std::endl          << std::endl;
    file << ForceField::n_type << " " << "atom types" << std::endl;
    file << "1 "               << " " << "bond types" << std::endl;

    file << std::endl   << std::endl;
    file << "Atoms"     << std::endl;

    for ( int i = 0; i < model.size(); i++ ) {
        // (tag) (mol id) (type) pos(x y z)
        file << model.tag[i] << " " << "1 " << model.type[i] + 1 << " ";
        for ( int d = 0; d < 3 ; d++ ) file << model.x[i][d] << " ";
        file << std::endl;
    }
}

// write only protein
inline void save_topology( const ProteContainer & model, RTParameter const & param, std::ostream && file ) {
    file << "Protein" << std::endl << std::endl;
    file << model.size() << " " <<  "atoms" << std::endl;
    file << model.bonds.size() << " " << "bonds" << std::endl; // what is num_bonds?
    file << std::endl << std::endl;

    // dim
    file << param.box[0][0] << " " << param.box[0][1] << " " << "xlo xhi" << std::endl;
    file << param.box[1][0] << " " << param.box[1][1] << " " << "ylo yhi" << std::endl;
    file << param.box[2][0] << " " << param.box[2][1] << " " << "zlo zhi" << std::endl;

    file << std::endl          << std::endl;
    file << ForceField::n_type << " " << "atom types" << std::endl;
    file << "1 "               << " " << "bond types" << std::endl;

    file << std::endl   << std::endl;
    file << "Atoms"     << std::endl;

    for ( int i = 0; i < model.size(); i++ ) {
        // (tag) (mol id) (type) pos(x y z)
        file << model.tag[i] << " " << "1 " << model.type[i] + 1 << " ";
        for ( int d = 0; d < 3 ; d++ ) file << model.x[i][d] << " ";
        file << std::endl;
    }

    file << std::endl << std::endl;
    file << "Bonds" << std::endl << std::endl;

    for ( int i = 0; i < model.bonds.size(); i++ ) {
        // (serial #) (type) (id first bond) (id second bond)
        file << ( i + 1 ) << " " << model.bonds[i].type  << " ";
        file << model.bonds[i].i << " " << model.bonds[i].j << std::endl;
    }
}

}

#endif /* TOPOLOGY_H_ */

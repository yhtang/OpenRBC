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
#ifndef OPENRBC_COMPUTE_TEMPERATURE_H_
#define OPENRBC_COMPUTE_TEMPERATURE_H_

#include "container.h"
#include "forcefield.h"

namespace openrbc {

// compute the temperature of the system
template<class CONTAINER1, class CONTAINER2>
double compute_temperature( CONTAINER1 const & model1, CONTAINER2 const & model2, RTParameter const & param ) {
    double ek = 0.0;
    for ( std::size_t i = 0; i < model1.size(); ++i ) ek += ForceField::mass[model1.type[i]] * normsq( model1.v[i] );
    for ( std::size_t i = 0; i < model2.size(); ++i ) ek += ForceField::mass[model2.type[i]] * normsq( model2.v[i] );
    return ek / ( 3.0 * ( model1.size() + model2.size() ) );
}

}

#endif

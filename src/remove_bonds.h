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
#ifndef REMOVE_BONDS_H_
#define REMOVE_BONDS_H_

#include <algorithm>
#include <iostream>
#include "container.h"

namespace openrbc {

template<class UnaryPredicate>
inline void remove_bonds( ProteContainer & protein, UnaryPredicate const & p ) {
    std::cout << "Number of bond " << protein.bonds.size();
    Bond * end = std::remove_if( protein.bonds.data(), protein.bonds.data() + protein.bonds.size(), p );
    protein.bonds.resize( end - protein.bonds.data() );
    std::cout << " --> " << protein.bonds.size() << std::endl;
}

struct UnaryPredicate {
    // return true if you want to remove this bond, otherwise false
    bool operator () ( Bond const & b ) const {
        return false;
    }
};

}

#endif

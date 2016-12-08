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
#ifndef ZERO_BULK_VELOCITY_H_
#define ZERO_BULK_VELOCITY_H_


#include "integrate_nh.h"
#include "voronoi.h"
#include "runtime_parameter.h"
#include "math_vector.h"
#include "aligned_array.h"

namespace openrbc {

using namespace config;

//-------------------------- Kernel --------------------------
struct zero_bulk_velocity {

    zero_bulk_velocity( RTParameter & p ) : parameter( p ) {}

    static std::string name() { return "zero_bulk_velocity"; }

    template<class CONTAINER>
    void operator () ( CONTAINER & container, int beg, int end ) {

    	config::vect local_velocity = 0;
        for ( int i = beg; i < end; ++i ) local_velocity += container.v[i];

        #pragma omp critical
		global_velocity += local_velocity;
        #pragma omp barrier

        auto com_velocity = global_velocity / container.size();
        for ( int i = beg; i < end; ++i ) container.v[i] -= com_velocity;

		if ( parameter.nstep % 1000 == 0 ) {
			#pragma omp master
			fprintf( stdout, "global velocity x=%lf y=%f z=%f\n", com_velocity[0], com_velocity[1], com_velocity[2] );
		}
    }

    config::vect global_velocity = 0;
    RTParameter & parameter;
};

}

#endif /* INTEGRATE_NH_H_ */

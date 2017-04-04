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
#ifndef OPENRBC_CONFIG_STATIC_H_
#define OPENRBC_CONFIG_STATIC_H_

#include "math_vector.h"

namespace openrbc {

namespace config {

#if 0
using real = float; // x, mu
using real = float; // v, omega
using real = float; // f, t
using T4 = float; // k-d tree
#endif

#define FUSED_INTEGRATOR
#define FUSED_PAIRWISE
#define LANGEVIN

using real = float;

#if   defined( EXPLICIT_SIMD )
	using vect = vector<real, 4>;
#elif defined( IMPLICIT_SIMD )
	using vect = vector<real, 3>;
#elif defined(_VEC4)
	using vect = vector<real, 4>;
#else
	using vect = vector<real, 3>;
#endif

const static std::size_t cacheline = 64;
}

namespace constant {

const static config::real zero    = 0;
const static config::real half    = 0.5;
const static config::real one     = 1;
const static config::real onehalf = 1.5;
const static config::real two     = 2;
const static config::real four    = 4;

}

}

#endif

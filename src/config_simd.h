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
#ifndef OPENRBC_CONFIG_SIMD_H_
#define OPENRBC_CONFIG_SIMD_H_

#if   0 || defined(_ESIMD)   // EXPLICIT SIMD, 4-COMPONENT VECTOR
	#define EXPLICIT_SIMD
	#if defined(_M_X64) || defined(_M_AMD64) || defined(__amd64__) || defined(__amd64) || defined(__x86_64__) || defined(__x86_64)
		#pragma message( "Explicit SIMD vectorization enabled for x86")
	#else
		#pragma message( "Explicit SIMD vectorization enabled, to get maximum performance on Power8 processors try 'make SIMD=3'")
	#endif
#elif 0 || defined(_ISIMD)   // IMPLICIT SIMD, 3-COMPONENT VECTOR
	#undef  EXPLICIT_SIMD
	#define IMPLICIT_SIMD
	#if defined(__powerpc64__) || defined(__ppc64__) || defined(_M_PPC) || defined(_ARCH_PPC64) || defined(__ppc)
		#pragma message( "Implicit SIMD vectorization enabled for Power8")
	#else
		#pragma message( "Implicit SIMD vectorization enabled, to get maximum performance on x86 processors try 'make SIMD=2'")
	#endif
#elif 0 || defined(_VEC4)    // SIMD OFF, 4-COMPONENT VECTOR
	#undef EXPLICIT_SIMD
	#undef IMPLICIT_SIMD
	#pragma message( "SIMD vectorization is not enabled, to get maximum performance try 'make SIMD=1'")
#else                        // SIMD OFF, 3-COMPONENT VECTOR
	#undef EXPLICIT_SIMD
	#undef IMPLICIT_SIMD
	#pragma message( "SIMD vectorization is not enabled, to get maximum performance try 'make SIMD=1'")
#endif

#endif

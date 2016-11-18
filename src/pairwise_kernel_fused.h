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
#ifndef OPENRBC_FORCE_KERNEL_FUSED_H_
#define OPENRBC_FORCE_KERNEL_FUSED_H_

#include "math_vector.h"
#include "pairwise_kernel.h"

namespace openrbc {

template<bool NEWTON1, bool NEWTON2>
inline void protein_lipid_omp(
    const vect dx, const real r_sq,
    const vect mu1, const vect mu2,
    ProteContainer & protein, LipidContainer & lipid,
    const int p1, const int p2, const int type, boolean<NEWTON1> newton1, boolean<NEWTON2> newton2 ) {

    const real r = std::sqrt( r_sq );
    const vect unitdel = dx * ( 1 / r );

    const real nidotnj = dot( mu1, mu2 );
    const real nidotunitdel = dot( mu1, unitdel );
    const real njdotunitdel = dot( mu2, unitdel );
    const real a = nidotnj - nidotunitdel * njdotunitdel;
    const real A = real( 1.0 ) + ForceField::alphalp[type] * ( a - real( 1.0 ) );

    const real rcminusr = ForceField::cutlp[type] - r;
    const real rcminusr3 = rcminusr * rcminusr * rcminusr;
    const real rcminusr4 = rcminusr * rcminusr3;
    const real rcminusr7 = rcminusr3 * rcminusr4;

    const vect pni = mu1 - nidotunitdel * unitdel;
    const vect pnj = mu2 - njdotunitdel * unitdel;

    const real ua = ForceField::attlp[type] * rcminusr4;
    const real alphaua = ForceField::alphalp[type] * ua;
    const real alphauar = alphaua / r;

    const real fra = real( 8.0 ) * ForceField::replp[type] * rcminusr7 + A * real( 4.0 ) * ForceField::attlp[type] * rcminusr3;
    const vect f = fra * unitdel + alphauar * ( njdotunitdel * pni + nidotunitdel * pnj );

    if ( NEWTON1 ) {
        protein.f[p1] += f;
        protein.t[p1] -= alphaua * pnj;
    }
    if ( NEWTON2 ) {
        lipid.f[p2] -= f;
        lipid.t[p2] -= alphaua * pni;
    }
}

template<bool NEWTON1, bool NEWTON2>
inline void lennard_jones_omp(
    const vect dx, const real rsq,
    Container & c1, Container & c2,
    const int p1, const int p2, const int type, boolean<NEWTON1> newton1, boolean<NEWTON2> newton2 ) {
    const real r2inv = real( 1 ) / rsq;
    const real r6inv = r2inv * r2inv * r2inv;
    const real forcelj = r6inv * ( ForceField::lj_lj1[type] * r6inv - ForceField::lj_lj2[type] );
    const real fpair = forcelj * r2inv;

    const vect f = dx * fpair;

    if ( NEWTON1 ) c1.f[p1] += f;
    if ( NEWTON2 ) c2.f[p2] -= f;
}

template<bool NEWTON>
inline void protein_protein_omp(
    const vect dx, const real r_sq,
    ProteContainer & protein,
    const int p1, const int p2, const int type, boolean<NEWTON> newton ) {
    const real r = std::sqrt( r_sq );
    const vect unitdel = dx / r;

    const real rcminusr = ForceField::cutpp[type] - r;
    const real rcminusr3 = rcminusr * rcminusr * rcminusr;
    const real rcminusr4 = rcminusr * rcminusr3;
    const real rcminusr7 = rcminusr3 * rcminusr4;

    const real fr = real( 8.0 ) * ForceField::reppp[type] * rcminusr7;
    const vect f = fr * unitdel;

    protein.f[p1] += f;
    if ( NEWTON ) protein.f[p2] -= f;
}

}

#endif

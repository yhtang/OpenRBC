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
#ifndef OPENRBC_FORCE_KERNEL_SIMD_H_
#define OPENRBC_FORCE_KERNEL_SIMD_H_

#include "math_vector.h"
#include "math_vector.h"
#include "pairwise_kernel.h"

namespace openrbc {

template<bool NEWTON>
inline void lipid_lipid_omp(
    const vect dx, const vect r_sq,
    const vect n1, const vect n2,
    LipidContainer & lipid,
    const int p1, const int p2,
    vect & f1, vect & t1,
    boolean<NEWTON> newton ) {

    const vect rinv    = _rsqrt( r_sq );
    const vect r       = r_sq * rinv;
    const vect unitdel = dx * rinv;

    const vect nidotnj = _dot( n1, n2 );
    const vect nidotunitdel = _dot( n1, unitdel );
    const vect njdotunitdel = _dot( n2, unitdel );
    const vect a = nidotnj - nidotunitdel * njdotunitdel;
    const vect A = vect( 1.0f ) + vect( ForceField::alphall ) * ( a - vect( 1.0f ) );

    const vect rcminusr = vect( ForceField::cutll ) - r;
    const vect rcminusr3 = rcminusr * rcminusr * rcminusr;
    const vect rcminusr4 = rcminusr * rcminusr3;
    const vect rcminusr7 = rcminusr3 * rcminusr4;

    const vect pni = n1 - nidotunitdel * unitdel;
    const vect pnj = n2 - njdotunitdel * unitdel;

    const vect ua = vect( ForceField::attll ) * rcminusr4;
    const vect alphaua = vect( ForceField::alphall ) * ua;
    const vect alphauar = alphaua * rinv;

    const vect fra = vect( 8.0f * ForceField::repll ) * rcminusr7 + A * vect( 4.0f * ForceField::attll ) * rcminusr3;
    const vect f = fra * unitdel + alphauar * ( njdotunitdel * pni + nidotunitdel * pnj );

    f1 += f;
    t1 -= alphaua * pnj;

    if ( NEWTON ) {
        lipid.f[p2] -= f;
        lipid.t[p2] -= alphaua * pni;
    }
}

template<bool NEWTON1, bool NEWTON2>
inline void protein_lipid_omp(
    const vect dx, const vect r_sq,
    const vect n1, const vect n2,
    ProteContainer & protein, LipidContainer & lipid,
    const int p1, const int p2, const int type, boolean<NEWTON1> newton1, boolean<NEWTON2> newton2 ) {

    const vect rinv = _rsqrt( r_sq );
    const vect r = r_sq * rinv;
    const vect unitdel = dx * rinv;

    const vect nidotnj = _dot( n1, n2 );
    const vect nidotunitdel = _dot( n1, unitdel );
    const vect njdotunitdel = _dot( n2, unitdel );
    const vect a = nidotnj - nidotunitdel * njdotunitdel;
    const vect A = vect( 1.0f ) + vect( ForceField::alphalp[type] ) * ( a - vect( 1.0f ) );

    const vect rcminusr = vect( ForceField::cutlp[type] ) - r;
    const vect rcminusr3 = rcminusr * rcminusr * rcminusr;
    const vect rcminusr4 = rcminusr * rcminusr3;
    const vect rcminusr7 = rcminusr3 * rcminusr4;

    const vect pni = n1 - nidotunitdel * unitdel;
    const vect pnj = n2 - njdotunitdel * unitdel;

    const vect ua = vect( ForceField::attlp[type] ) * rcminusr4;
    const vect alphaua = vect( ForceField::alphalp[type] ) * ua;
    const vect alphauar = alphaua * rinv;

    const vect fra = vect( 8.0f * ForceField::replp[type] ) * rcminusr7 + A * vect( 4.0f * ForceField::attlp[type] ) * rcminusr3;
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
    const vect dx, const vect rsq,
    Container & c1, Container & c2,
    const int p1, const int p2, const int type, boolean<NEWTON1> newton1, boolean<NEWTON2> newton2 ) {

    const vect rinv = _rsqrt( rsq );
    const vect r2inv = rinv * rinv;
    const vect r6inv = r2inv * r2inv * r2inv;
    const vect forcelj = r6inv * ( vect( ForceField::lj_lj1[type] ) * r6inv - vect( ForceField::lj_lj2[type] ) );
    //const vect forcelj_truncated = _min( forcelj, 5.0f );
    const vect fpair = forcelj * r2inv;

    const vect f = dx * fpair;

    if ( NEWTON1 ) c1.f[p1] += f;
    if ( NEWTON2 ) c2.f[p2] -= f;
}

template<bool NEWTON>
inline void protein_protein_omp(
    const vect dx, const vect r_sq,
    ProteContainer & protein,
    const int p1, const int p2, const int type, boolean<NEWTON> newton ) {

    const vect rinv = _rsqrt( r_sq );
    const vect r = r_sq * rinv;
    const vect unitdel = dx * rinv;

    const vect rcminusr = vect( ForceField::cutpp[type] ) - r;
    const vect rcminusr3 = rcminusr * rcminusr * rcminusr;
    const vect rcminusr4 = rcminusr * rcminusr3;
    const vect rcminusr7 = rcminusr3 * rcminusr4;

    const vect fr = vect( 8.0f * ForceField::reppp[type] ) * rcminusr7;
    const vect f = fr * unitdel;

    protein.f[p1] += f;
    if ( NEWTON ) protein.f[p2] -= f;
}

}

#endif

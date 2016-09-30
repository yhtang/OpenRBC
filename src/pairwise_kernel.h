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
#ifndef FORCE_KERNEL_H_
#define FORCE_KERNEL_H_

#include "math_vector.h"
#include "math_vector.h"

namespace openrbc {

// helper to achieve static dispatching
template<bool TRUTH> struct boolean {};
const static boolean<true>  newton_on = boolean<true>();
const static boolean<false> newton_off = boolean<false>();

/********************* OpenMP version ***************************/

using namespace config;

template<bool NEWTON>
inline void lipid_lipid_omp(
    const vect dx, const real r_sq,
    const vect mu1, const vect mu2,
    LipidContainer & c,
    const int p1, const int p2, boolean<NEWTON> newton ) {

    const real r = std::sqrt( r_sq );
    const vect unitdel = dx * ( real( 1 ) / r );

    const real nidotnj = dot( mu1, mu2 );
    const real nidotunitdel = dot( mu1, unitdel );
    const real njdotunitdel = dot( mu2, unitdel );
    const real a = nidotnj - nidotunitdel * njdotunitdel;
    const real A = real( 1.0 ) + ForceField::alphall * ( a - real( 1.0 ) );

    const real rcminusr = ForceField::cutll - r;
    const real rcminusr3 = rcminusr * rcminusr * rcminusr;
    const real rcminusr4 = rcminusr * rcminusr3;
    const real rcminusr7 = rcminusr3 * rcminusr4;

    const vect pni = mu1 - nidotunitdel * unitdel;
    const vect pnj = mu2 - njdotunitdel * unitdel;

    const real ua = ForceField::attll * rcminusr4;
    const real alphaua = ForceField::alphall * ua;
    const real alphauar = alphaua / r;

    const real fra = real( 8.0 ) * ForceField::repll * rcminusr7 + A * real( 4.0 ) * ForceField::attll * rcminusr3;
    const vect f = fra * unitdel + alphauar * ( njdotunitdel * pni + nidotunitdel * pnj );

    c.f[p1] += f;
    c.t[p1] -= alphaua * pnj;

    if ( NEWTON ) {
        c.f[p2] -= f;
        c.t[p2] -= alphaua * pni;
    }
}

inline void lipid_protein_omp(
    const vect dx, const real r_sq,
    const vect mu1, const vect mu2,
    Container & lipid, Container & protein,
    const int p1, const int p2, const int type ) {

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

    #pragma omp atomic
    lipid.f[p1][0] += f[0];
    #pragma omp atomic
    lipid.f[p1][1] += f[1];
    #pragma omp atomic
    lipid.f[p1][2] += f[2];

    #pragma omp atomic
    protein.f[p2][0] -= f[0];
    #pragma omp atomic
    protein.f[p2][1] -= f[1];
    #pragma omp atomic
    protein.f[p2][2] -= f[2];

    // t
    #pragma omp atomic
    lipid.t[p1][0] -= alphaua * pnj[0];
    #pragma omp atomic
    lipid.t[p1][1] -= alphaua * pnj[1];
    #pragma omp atomic
    lipid.t[p1][2] -= alphaua * pnj[2];

    #pragma omp atomic
    protein.t[p2][0] -= alphaua * pni[0];
    #pragma omp atomic
    protein.t[p2][1] -= alphaua * pni[1];
    #pragma omp atomic
    protein.t[p2][2] -= alphaua * pni[2];
}

inline void lennard_jones_omp(
    const vect dx, const real rsq,
    Container & c1, Container & c2,
    const int p1, const int p2, const int type ) {
    const real r2inv = real( 1 ) / rsq;
    const real r6inv = r2inv * r2inv * r2inv;
    const real forcelj = r6inv * ( ForceField::lj_lj1[type] * r6inv - ForceField::lj_lj2[type] );
    const real fpair = forcelj * r2inv;

    const vect f = dx * fpair;

    #pragma omp atomic
    c1.f[p1][0] += f[0];
    #pragma omp atomic
    c1.f[p1][1] += f[1];
    #pragma omp atomic
    c1.f[p1][2] += f[2];

    #pragma omp atomic
    c2.f[p2][0] -= f[0];
    #pragma omp atomic
    c2.f[p2][1] -= f[1];
    #pragma omp atomic
    c2.f[p2][2] -= f[2];
}

inline void protein_protein_omp(
    const vect dx, const real r_sq,
    Container & protein,
    const int p1, const int p2, const int type ) {
    const real r = std::sqrt( r_sq );
    const vect unitdel = dx / r;

    const real rcminusr = ForceField::cutpp[type] - r;
    const real rcminusr7 = std::pow( rcminusr, 7 );

    const real fr = real( 8.0 ) * ForceField::reppp[type] * rcminusr7;
    const vect f = fr * unitdel;

    #pragma omp atomic
    protein.f[p1][0] += f[0];
    #pragma omp atomic
    protein.f[p1][1] += f[1];
    #pragma omp atomic
    protein.f[p1][2] += f[2];
    #pragma omp atomic
    protein.f[p2][0] -= f[0];
    #pragma omp atomic
    protein.f[p2][1] -= f[1];
    #pragma omp atomic
    protein.f[p2][2] -= f[2];
}

/********************* Serial version ***************************/
inline void lipid_lipid_serial(
    const real dx, const real dy, const real dz, const real r_sq,
    const real mux1, const real muy1, const real muz1,
    const real mux2, const real muy2, const real muz2,
    Container & c,
    const int p1, const int p2 ) {

    const real r = std::sqrt( r_sq );
    const real unitdelx = dx / r;
    const real unitdely = dy / r;
    const real unitdelz = dz / r;

    const real nidotnj = mux1 * mux2 + muy1 * muy2 + muz1 * muz2;
    const real nidotunitdel = mux1 * unitdelx + muy1 * unitdely + muz1 * unitdelz;
    const real njdotunitdel = mux2 * unitdelx + muy2 * unitdely + muz2 * unitdelz;
    const real a = nidotnj - nidotunitdel * njdotunitdel;
    const real A = real( 1.0 ) + ForceField::alphall * ( a - real( 1.0 ) );

    const real rcminusr = ForceField::cutll - r;
    const real rcminusr3 = rcminusr * rcminusr * rcminusr;
    const real rcminusr4 = rcminusr * rcminusr3;
    const real rcminusr7 = rcminusr3 * rcminusr4;

    const real pnix = mux1 - nidotunitdel * unitdelx;
    const real pniy = muy1 - nidotunitdel * unitdely;
    const real pniz = muz1 - nidotunitdel * unitdelz;
    const real pnjx = mux2 - njdotunitdel * unitdelx;
    const real pnjy = muy2 - njdotunitdel * unitdely;
    const real pnjz = muz2 - njdotunitdel * unitdelz;

    const real ua = ForceField::attll * rcminusr4;
    const real alphaua = ForceField::alphall * ua;
    const real alphauar = alphaua / r;

    const real fra = real( 8.0 ) * ForceField::repll * rcminusr7 + A * real( 4.0 ) * ForceField::attll * rcminusr3;
    const real fx = fra * unitdelx + alphauar * ( njdotunitdel * pnix + nidotunitdel * pnjx );
    const real fy = fra * unitdely + alphauar * ( njdotunitdel * pniy + nidotunitdel * pnjy );
    const real fz = fra * unitdelz + alphauar * ( njdotunitdel * pniz + nidotunitdel * pnjz );

    c.f[p1][0] += fx;
    c.f[p1][1] += fy;
    c.f[p1][2] += fz;
    c.f[p2][0] -= fx;
    c.f[p2][1] -= fy;
    c.f[p2][2] -= fz;

    c.t[p1][0] -= alphaua * pnjx; // t
    c.t[p1][1] -= alphaua * pnjy;
    c.t[p1][2] -= alphaua * pnjz;
    c.t[p2][0] -= alphaua * pnix;
    c.t[p2][1] -= alphaua * pniy;
    c.t[p2][2] -= alphaua * pniz;
}

inline void lipid_protein_serial(
    const vect dx, const real r_sq,
    const vect mu1, const vect mu2,
    Container & lipid, Container & protein,
    const int p1, const int p2, const int type ) {

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

    lipid  .f[p1] += f;
    protein.f[p2] -= f;

    // t
    lipid  .t[p1] -= alphaua * pnj;
    protein.t[p2] -= alphaua * pni;
}

inline void lennard_jones_serial(
    const vect dx, const real rsq,
    Container & c1, Container & c2,
    const int p1, const int p2, const int type, RTParameter const & param ) {
    const real r2inv = real( 1 ) / rsq;
    const real r6inv = r2inv * r2inv * r2inv;
    const real forcelj = r6inv * ( ForceField::lj_lj1[type] * r6inv - ForceField::lj_lj2[type] );
    const real fpair = forcelj * r2inv;

    const vect f = dx * fpair;

    c1.f[p1] += f;
    c2.f[p2] -= f;
}

inline void protein_protein_serial(
    const vect dx, const real r_sq,
    Container & protein,
    const int p1, const int p2, const int type ) {
    const real r = std::sqrt( r_sq );
    const vect unitdel = dx / r;

    const real rcminusr = ForceField::cutpp[type] - r;
    const real rcminusr7 = std::pow( rcminusr, 7 );

    const real fr = real( 8.0 ) * ForceField::reppp[type] * rcminusr7;
    const vect f = fr * unitdel;

    protein.f[p1] += f;
    protein.f[p2] -= f;
}

}

#endif /* FORCE_KERNEL_H_ */

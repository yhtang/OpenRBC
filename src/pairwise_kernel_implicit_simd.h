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
#ifndef OPENRBC_FORCE_KERNEL_IMPLICIT_SIMD_H_
#define OPENRBC_FORCE_KERNEL_IMPLICIT_SIMD_H_

#include "math_vector.h"
#include "pairwise_kernel.h"

namespace openrbc {

const static boolean<true>  flush_all {};
const static boolean<false> flush_any {};

template<bool NEWTON>
struct lipid_lipid_implicit_simd {
    constexpr static int simd_width = 8;

    int queue = 0;
    int p1[simd_width], p2[simd_width];
    real drx[simd_width], dry[simd_width], drz[simd_width];
    real n1x[simd_width], n1y[simd_width], n1z[simd_width];
    real n2x[simd_width], n2y[simd_width], n2z[simd_width];
    real rsq[simd_width];
    LipidContainer & lipid_;

    lipid_lipid_implicit_simd( LipidContainer & lipid ) : lipid_( lipid ) {}
    ~lipid_lipid_implicit_simd() {
        flush( flush_any );
    }

    inline void enqueue( const int i, const int j, const vect dx, const real r_sq, const vect n1, const vect n2 ) {
        p1[queue] = i;
        p2[queue] = j;
        drx[queue] = dx[0], dry[queue] = dx[1], drz[queue] = dx[2];
        n1x[queue] = n1[0], n1y[queue] = n1[1], n1z[queue] = n1[2];
        n2x[queue] = n2[0], n2y[queue] = n2[1], n2z[queue] = n2[2];
        rsq[queue] = r_sq;
        if ( ++queue == simd_width ) flush( flush_all );
    }

    template<bool COMPLETE>
    inline void flush( boolean<COMPLETE> foo ) {
    	//printf("L-L COMPLETE %d\n", COMPLETE);
        real fx[simd_width], fy[simd_width], fz[simd_width];
        real t1x[simd_width], t1y[simd_width], t1z[simd_width];
        real t2x[simd_width], t2y[simd_width], t2z[simd_width];
        for ( int k = 0; k < ( COMPLETE ? simd_width : queue ); ++k ) {
            real r         = std::sqrt( rsq[k] );
            real ex        = drx[k] / r;
            real ey        = dry[k] / r;
            real ez        = drz[k] / r;
            real nidotnj   = n1x[k] * n2x[k] + n1y[k] * n2y[k] + n1z[k] * n2z[k];
            real nidote    = n1x[k] *  ex + n1y[k] *  ey + n1z[k] *  ez;
            real njdote    = n2x[k] *  ex + n2y[k] *  ey + n2z[k] *  ez;
            real a         = nidotnj - nidote * njdote;
            real A         = 1.f + ForceField::alphall * ( a - 1.f );
            real rcminusr  = ForceField::cutll - r;
            real rcminusr3 = rcminusr * rcminusr * rcminusr;
            real rcminusr4 = rcminusr * rcminusr3;
            real rcminusr7 = rcminusr3 * rcminusr4;
            real pnix      = n1x[k] - nidote * ex;
            real pniy      = n1y[k] - nidote * ey;
            real pniz      = n1z[k] - nidote * ez;
            real pnjx      = n2x[k] - njdote * ex;
            real pnjy      = n2y[k] - njdote * ey;
            real pnjz      = n2z[k] - njdote * ez;
            real ua        = ForceField::attll * rcminusr4;
            real alphaua   = ForceField::alphall * ua;
            real alphauar  = alphaua / r;
            real fra       = 8.0f * ForceField::repll * rcminusr7 + A * 4.f * ForceField::attll * rcminusr3;
            fx [k] = fra * ex + alphauar * ( njdote * pnix + nidote * pnjx );
            fy [k] = fra * ey + alphauar * ( njdote * pniy + nidote * pnjy );
            fz [k] = fra * ez + alphauar * ( njdote * pniz + nidote * pnjz );
            t1x[k] = alphaua * pnjx;
            t1y[k] = alphaua * pnjy;
            t1z[k] = alphaua * pnjz;
            if ( NEWTON ) t2x[k] = alphaua * pnix;
            if ( NEWTON ) t2y[k] = alphaua * pniy;
            if ( NEWTON ) t2z[k] = alphaua * pniz;
        }
        for ( int k = 0; k < ( COMPLETE ? simd_width : queue ); ++k ) {
            lipid_.f[ p1[k] ] += vect( fx [k], fy [k], fz [k] );
            lipid_.t[ p1[k] ] -= vect( t1x[k], t1y[k], t1z[k] );
            if ( NEWTON ) lipid_.f[ p2[k] ] -= vect( fx [k], fy [k], fz [k] );
            if ( NEWTON ) lipid_.t[ p2[k] ] -= vect( t2x[k], t2y[k], t2z[k] );
        }
        queue = 0;
    }
};

template<bool NEWTON1, bool NEWTON2>
struct protein_lipid_implicit_simd {
    constexpr static int simd_width = 4;

    int queue = 0;
    int p1[simd_width], p2[simd_width];
    real drx[simd_width], dry[simd_width], drz[simd_width];
    real n1x[simd_width], n1y[simd_width], n1z[simd_width];
    real n2x[simd_width], n2y[simd_width], n2z[simd_width];
    real rsq[simd_width];
    real alphalp[simd_width], cutlp[simd_width], attlp[simd_width], replp[simd_width];
    ProteContainer & prote_;
    LipidContainer & lipid_;

    protein_lipid_implicit_simd( ProteContainer & prote, LipidContainer & lipid ) : prote_( prote ), lipid_( lipid ) {}
    ~protein_lipid_implicit_simd() {
        flush( flush_any );
    }

    inline void enqueue( const int i, const int j, const int t, const vect dx, const real r_sq, const vect n1, const vect n2 ) {
        p1[queue] = i;
        p2[queue] = j;
        drx[queue] = dx[0], dry[queue] = dx[1], drz[queue] = dx[2];
        n1x[queue] = n1[0], n1y[queue] = n1[1], n1z[queue] = n1[2];
        n2x[queue] = n2[0], n2y[queue] = n2[1], n2z[queue] = n2[2];
        rsq[queue] = r_sq;
        alphalp[queue] = ForceField::alphalp[t];
        cutlp  [queue] = ForceField::cutlp  [t];
        attlp  [queue] = ForceField::attlp  [t];
        replp  [queue] = ForceField::replp  [t];
        if ( ++queue == simd_width ) flush( flush_all );
    }

    template<bool COMPLETE>
    inline void flush( boolean<COMPLETE> foo ) {
        real fx[simd_width], fy[simd_width], fz[simd_width];
        real t1x[simd_width], t1y[simd_width], t1z[simd_width];
        real t2x[simd_width], t2y[simd_width], t2z[simd_width];
        for ( int k = 0; k < ( COMPLETE ? simd_width : queue ); ++k ) {
            real r         = std::sqrt( rsq[k] );
            real ex        = drx[k] / r;
            real ey        = dry[k] / r;
            real ez        = drz[k] / r;
            real nidotnj   = n1x[k] * n2x[k] + n1y[k] * n2y[k] + n1z[k] * n2z[k];
            real nidote    = n1x[k] *  ex + n1y[k] *  ey + n1z[k] *  ez;
            real njdote    = n2x[k] *  ex + n2y[k] *  ey + n2z[k] *  ez;
            real a         = nidotnj - nidote * njdote;
            real A         = 1.f + alphalp[k] * ( a - 1.f );
            real rcminusr  = cutlp[k] - r;
            real rcminusr3 = rcminusr * rcminusr * rcminusr;
            real rcminusr4 = rcminusr * rcminusr3;
            real rcminusr7 = rcminusr3 * rcminusr4;
            real pnix      = n1x[k] - nidote * ex;
            real pniy      = n1y[k] - nidote * ey;
            real pniz      = n1z[k] - nidote * ez;
            real pnjx      = n2x[k] - njdote * ex;
            real pnjy      = n2y[k] - njdote * ey;
            real pnjz      = n2z[k] - njdote * ez;
            real ua        = attlp[k] * rcminusr4;
            real alphaua   = alphalp[k] * ua;
            real alphauar  = alphaua / r;
            real fra       = 8.0f * replp[k] * rcminusr7 + A * 4.f * attlp[k] * rcminusr3;
            fx [k] = fra * ex + alphauar * ( njdote * pnix + nidote * pnjx );
            fy [k] = fra * ey + alphauar * ( njdote * pniy + nidote * pnjy );
            fz [k] = fra * ez + alphauar * ( njdote * pniz + nidote * pnjz );
            if ( NEWTON1 ) t1x[k] = alphaua * pnjx;
            if ( NEWTON1 ) t1y[k] = alphaua * pnjy;
            if ( NEWTON1 ) t1z[k] = alphaua * pnjz;
            if ( NEWTON2 ) t2x[k] = alphaua * pnix;
            if ( NEWTON2 ) t2y[k] = alphaua * pniy;
            if ( NEWTON2 ) t2z[k] = alphaua * pniz;
        }
        for ( int k = 0; k < ( COMPLETE ? simd_width : queue ); ++k ) {
            if ( NEWTON1 ) prote_.f[ p1[k] ] += vect( fx [k], fy [k], fz [k] );
            if ( NEWTON1 ) prote_.t[ p1[k] ] -= vect( t1x[k], t1y[k], t1z[k] );
            if ( NEWTON2 ) lipid_.f[ p2[k] ] -= vect( fx [k], fy [k], fz [k] );
            if ( NEWTON2 ) lipid_.t[ p2[k] ] -= vect( t2x[k], t2y[k], t2z[k] );
        }
        queue = 0;
    }
};

template<bool NEWTON1, bool NEWTON2>
struct lennard_jones_implicit_simd {
    constexpr static int simd_width = 4;

    int queue = 0;
    int p1[simd_width], p2[simd_width];
    real drx[simd_width], dry[simd_width], drz[simd_width];
    real rsq[simd_width];
    real lj1[simd_width], lj2[simd_width];
    Container & c1_;
    Container & c2_;

    lennard_jones_implicit_simd( Container & c1, Container & c2 ) : c1_( c1 ), c2_( c2 ) {}
    ~lennard_jones_implicit_simd() {
        flush( flush_any );
    }

    inline void enqueue( const int i, const int j, const int t, const vect dx, const real r_sq ) {
        p1[queue] = i;
        p2[queue] = j;
        drx[queue] = dx[0], dry[queue] = dx[1], drz[queue] = dx[2];
        rsq[queue] = r_sq;
        lj1[queue] = ForceField::lj_lj1[t];
        lj2[queue] = ForceField::lj_lj2[t];
        if ( ++queue == simd_width ) flush( flush_all );
    }

    template<bool COMPLETE>
    inline void flush( boolean<COMPLETE> foo ) {
        real fx[simd_width], fy[simd_width], fz[simd_width];
        for ( int k = 0; k < ( COMPLETE ? simd_width : queue ); ++k ) {
            real r2inv     = 1.f / rsq[k];
            real r6inv     = r2inv * r2inv * r2inv;
            real forcelj   = r6inv * ( lj1[k] * r6inv - lj2[k] );
            real fpair     = forcelj * r2inv;
            fx[k] = fpair * drx[k];
            fy[k] = fpair * dry[k];
            fz[k] = fpair * drz[k];
        }
        for ( int k = 0; k < ( COMPLETE ? simd_width : queue ); ++k ) {
            if ( NEWTON1 ) c1_.f[ p1[k] ] += vect( fx[k], fy[k], fz[k] );
            if ( NEWTON2 ) c2_.f[ p2[k] ] -= vect( fx[k], fy[k], fz[k] );
        }
        queue = 0;
    }
};

template<bool NEWTON>
struct protein_protein_implicit_simd {
    constexpr static int simd_width = 8;

    int queue = 0;
    int p1[simd_width], p2[simd_width];
    real drx[simd_width], dry[simd_width], drz[simd_width];
    real rsq[simd_width];
    real cutpp[simd_width], reppp[simd_width];
    ProteContainer & protein_;

    protein_protein_implicit_simd( ProteContainer & protein ) : protein_( protein ) {}
    ~protein_protein_implicit_simd() {
        flush( flush_any );
    }

    inline void enqueue( const int i, const int j, const int t, const vect dx, const real r_sq ) {
        p1[queue] = i;
        p2[queue] = j;
        drx[queue] = dx[0], dry[queue] = dx[1], drz[queue] = dx[2];
        rsq[queue] = r_sq;
        cutpp[queue] = ForceField::cutpp[t];
        reppp[queue] = ForceField::reppp[t];
        if ( ++queue == simd_width ) flush( flush_all );
    }

    template<bool COMPLETE>
    inline void flush( boolean<COMPLETE> foo ) {
        real fx[simd_width], fy[simd_width], fz[simd_width];
        for ( int k = 0; k < ( COMPLETE ? simd_width : queue ); ++k ) {
            real r         = std::sqrt( rsq[k] );
            real ex        = drx[k] / r;
            real ey        = dry[k] / r;
            real ez        = drz[k] / r;
            real rcminusr  = cutpp[k] - r;
            real rcminusr3 = rcminusr * rcminusr * rcminusr;
            real rcminusr4 = rcminusr * rcminusr3;
            real rcminusr7 = rcminusr3 * rcminusr4;
            real fr        = 8.0f * reppp[k] * rcminusr7;
            fx[k] = fr * ex;
            fy[k] = fr * ey;
            fz[k] = fr * ez;
        }
        for ( int k = 0; k < ( COMPLETE ? simd_width : queue ); ++k ) {
            protein_.f[ p1[k] ] += vect( fx[k], fy[k], fz[k] );
            if ( NEWTON ) protein_.f[ p2[k] ] -= vect( fx[k], fy[k], fz[k] );
        }
        queue = 0;
    }
};

}

#endif

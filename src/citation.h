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
#ifndef OPENRBC_CITATION_H_
#define OPENRBC_CITATION_H_

namespace openrbc {

static struct Citation {
    Citation() {
        print_citation();
        atexit( print_citation );
    }
    static void print_citation() {
        printf( "|*******************************************************************************\n"
                "|* OpenRBC http://openrbc.io/                                                  \n"
                "|* Yu-Hang Tang, Lu Lu 2016 All Rights Reserved                                \n"
                "|* This code is released under the Apache 2.0 license                          \n"
                "|* In addition, any published document using OpenRBC must cite:                \n"
                "|* Yu-Hang Tang, Lu Lu, Leopold Grinberg, He Li, Vipin Sachdeva, Constantinos  \n"
                "|* Evangelinos, and George Em Karniadakis, Biophysical Journal, manuscript submitted.\n"
                "|*******************************************************************************\n"
              );
    }
} citation;

}

#endif

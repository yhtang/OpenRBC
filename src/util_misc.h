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
#ifndef OPENRBC_UTIL_MISC_H_
#define OPENRBC_UTIL_MISC_H_

#include<cstring>
#include<string>
#include<cstdlib>

namespace openrbc {

#ifdef _USE_CUSTOM_ALIGNED_ALLOC
void * aligned_alloc( size_t __alignment, size_t __size ) {
    return malloc( __size );
}
#endif

// aligned allocation
void * aligned_allocate( size_t alignment, size_t size ) {
    if ( size % alignment ) size += alignment - size % alignment;
    return aligned_alloc( alignment, size );
}

void * aligned_realloc ( void * ptr, std::size_t alignment, std::size_t old_size, std::size_t new_size ) {
    void * nptr = aligned_allocate ( alignment, new_size );
    std::memcpy ( nptr, ptr, std::min ( new_size, old_size ) );
    std::free ( ptr );
    return nptr;
}

inline double get_time_posix() {
    struct timespec time;
    clock_gettime ( CLOCK_REALTIME, &time );
    return ( double ) time.tv_sec + ( double ) time.tv_nsec * 1.0e-9 ;
}

float _rsqrt( float n ) {
    std::int32_t i = 0x5f3759df - ( ( * ( long * ) &n ) >> 1 );
    float y        = * ( float * ) &i;
    return y * ( 1.5f - ( 0.5f * n * y * y ) );
}

#define rt_assert_impl(predicate,message,file,line) \
        if ( !(predicate) ) {\
            printf("<OpenRBC> %s %d RUNTIME ASSERTION FAILED\n<OpenRBC>     CONDITION NOT SATISFIED: %s\n<OpenRBC>     REASON: %s\n",file,line,#predicate, message);\
            exit(0);\
        }
#define rt_assert(predicate,message) rt_assert_impl(predicate,message,__FILE__,__LINE__)

}

#endif

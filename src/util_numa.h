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
#ifndef UTIL_NUMA_H_
#define UTIL_NUMA_H_

#include <map>
#include <vector>
#include <cstdint>
#include <omp.h>

namespace openrbc {

struct LoadBalancer {
    LoadBalancer() = default;
    LoadBalancer( std::size_t range ) {
        set_range( range );
    }

    void set_range( std::size_t range ) {
        if ( omp_in_parallel() ) {
            fprintf( stderr, "<OpenRBC> LoadBalancer must be initialized outside of any parallel region\n" );
        }
        if ( range == range_ && omp_get_max_threads() == ntd_ ) return;

        range_ = range;
        ntd_ = omp_get_max_threads();
        beg_.resize( ntd_ );
        end_.resize( ntd_ );
        for ( std::size_t tid = 0; tid < ntd_; ++tid ) {
            beg_[tid] = tid * range_ / ntd_;
            end_[tid] = ( tid == ntd_ - 1 ) ? range_ : ( tid + 1 ) * range_ / ntd_;
        }
        set_ = true;
    }

    bool        set()                                const { return set_; }
    int         ntd()                                const { return ntd_; }
    std::size_t beg( int tid = omp_get_thread_num() ) const { return beg_[tid]; }
    std::size_t end( int tid = omp_get_thread_num() ) const { return end_[tid]; }

    struct iterator {
        const std::size_t n_, block_size_;
        const int tid_, ntd_;
        std::size_t block_beg_;
        int i, end;
        bool exhausted = false;
        iterator( std::size_t range, std::size_t block_size, int tid, int ntd ) : n_( range ), block_size_( block_size ), tid_( tid ), ntd_( ntd ), block_beg_( 0 ) {
            if ( block_beg_ < range ) compute_block_range();
            else exhausted = true;
        }
        void operator ++ () {
            if ( ++i >= end ) {
                block_beg_ += block_size_;
                if ( block_beg_ < n_ ) compute_block_range();
                else exhausted = true;
            }
        }
        int  operator * () const { return i; }
    protected:
        void compute_block_range() {
            int actual_size = std::min( n_ - block_beg_, block_size_ );
            i   = block_beg_ + actual_size * tid_ / ntd_;
            end = ( tid_ == ntd_ - 1 ) ? ( block_beg_ + actual_size ) : ( block_beg_ + actual_size * ( tid_ + 1 ) / ntd_ );
        }
    };

    iterator range( int tid = omp_get_thread_num(), std::size_t share = 2048 ) const {
        return iterator { range_, share * ntd_, tid, ntd_ };
    }

protected:
    bool        set_   = false;
    int         ntd_   = 0;
    std::size_t range_ = 0;
    std::vector<std::size_t> beg_, end_;
};

using BalancerMap = std::map<std::string, LoadBalancer>;

}

#endif

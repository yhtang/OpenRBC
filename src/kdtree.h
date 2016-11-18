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
#ifndef OPENRBC_KDTREE_H_
#define OPENRBC_KDTREE_H_

#include "math_vector.h"
#include "service.h"
#include "timer.h"
#include "util_misc.h"
#include <vector>
#include <string>
#include <functional>
#include <deque>

namespace openrbc {

template<typename SCALAR, uint DIM, typename SIZET = std::int64_t>
class KDTree {
protected:
    using point = vector<SCALAR, DIM>;
    AlignedArray<point> pts_;
    AlignedArray<SIZET, true> idx_;
    SIZET npts_ = 0;
    const static SIZET leaf_size = 16;

public:
    struct node {
        uint d = 0;
        SCALAR lo, hi;
        union { SCALAR lim[2]; SIZET plr[2]; };
        node * parent = nullptr, *child = nullptr, *next = nullptr;

        node() {
            lim[0] = -std::numeric_limits<SCALAR>::max();
            lim[1] =  std::numeric_limits<SCALAR>::max();
        }
    } root;

public:
    KDTree() : allocator( omp_get_max_threads() ) {}
    KDTree( KDTree const & ) = delete;
    KDTree & operator = ( KDTree const & ) = delete;
    ~KDTree() { clear(); }

    template<bool SHADOW, std::size_t ALIGN>
    void build_serial( std::tuple<node *, SIZET, SIZET> seed, AlignedArray<point, SHADOW, ALIGN> const & points ) {
        std::vector<std::tuple<node *, SIZET, SIZET> > stack( 1, seed );
        while ( stack.size() ) {
            node * n;
            SIZET l, r;
            std::tie( n, l, r ) = stack.back();
            stack.pop_back();

            // adaptively find bisecting dimension
            point lo =  std::numeric_limits<SCALAR>::max();
            point hi = -std::numeric_limits<SCALAR>::max();
            for ( int i = l; i < r; i++ ) {
                lo = apply( lo, points[ idx_[i] ], functors::min<SCALAR>() );
                hi = apply( hi, points[ idx_[i] ], functors::max<SCALAR>() );
            }
            uint d = argmax( hi - lo );
            n->d   = d;
            n->lo  = lo[d];
            n->hi  = hi[d];

            if ( r - l > leaf_size ) {
                SIZET m, qsel_left = l, qsel_right = r;
                do {
                    SIZET g = qsel_left + rand() % ( qsel_right - qsel_left );
                    SCALAR pivot = points[ idx_[ g ] ][d];

                    auto i = idx_.data() + qsel_left;
                    auto j = idx_.data() + qsel_right - 1;

                    do {
                        while ( i - idx_.data() < r && points[*i][d] <  pivot ) ++i;
                        while ( j - idx_.data() > l && points[*j][d] >= pivot ) --j;
                        if ( j - i > 0 ) std::swap( *i, *j );
                    } while ( j - i > 1 );

                    m = std::min( i - idx_.data(), j - idx_.data() ) + 1;
                    if ( m == qsel_left || m == qsel_right ) {
                        bool same = true;
                        for ( int i = qsel_left; same && i < qsel_right; i++ ) same &= ( points[ idx_[qsel_left] ][d] == points[ idx_[i] ][d] );
                        if ( same ) break;
                    } else if ( m == ( l + r ) / 2 ) break;
                    else if ( m < ( l + r ) / 2 ) qsel_left = m;
                    else qsel_right = m;

                } while ( true );

                for ( int k = l; k < m; k++ ) n->lim[0] = std::max( n->lim[0], points[idx_[k]][d] );
                for ( int k = m; k < r; k++ ) n->lim[1] = std::min( n->lim[1], points[idx_[k]][d] );
                n->child = allocator[ omp_get_thread_num() ].allocate( 2 );
                n->child[0].parent = n->child[1].parent = n;
                stack.emplace_back( n->child + 0, l, m );
                stack.emplace_back( n->child + 1, m, r );
            } else {
                n->plr[0] = l;
                n->plr[1] = r;
            }
        }
    }

    template<bool SHADOW, std::size_t ALIGN>
    void build( AlignedArray<point, SHADOW, ALIGN> const & points ) {
        Service<Timers>::call()["|KDTree::build"].start();

        npts_ = points.size();
        pts_.resize( points.size() );
        idx_.resize( points.size() );
        #pragma omp parallel for
        for ( SIZET i = 0; i < idx_.size(); i++ ) idx_[i] = i;

        root.child = root.parent = root.next = nullptr;
        for ( auto & pool : allocator ) pool.refresh();

        SIZET ntd = omp_get_max_threads();

        if ( points.size() < ntd * leaf_size ) {
            build_serial( std::make_tuple( &root, 0, points.size() ), points );
        } else {
            std::vector<std::tuple<node *, SIZET, SIZET> > stack( 1, std::make_tuple( &root, 0, points.size() ) ), stack_new( 2 );

            while ( stack.size() < ntd ) {

                #pragma omp parallel
                {
                    SIZET stride = omp_get_num_threads() / stack.size();
                    SIZET gid = omp_get_thread_num() / stride;
                    SIZET tid = omp_get_thread_num() % stride;

                    if ( tid == 0 && gid < stack.size() ) {
                        node * n;
                        SIZET l, r;
                        std::tie( n, l, r ) = stack[ gid ];

                        // adaptively find bisecting dimension
                        point lo =  std::numeric_limits<SCALAR>::max();
                        point hi = -std::numeric_limits<SCALAR>::max();
                        for ( SIZET i = l; i < r; i++ ) {
                            lo = apply( lo, points[ idx_[i] ], functors::min<SCALAR>() );
                            hi = apply( hi, points[ idx_[i] ], functors::max<SCALAR>() );
                        }
                        uint d = argmax( hi - lo );
                        n->d   = d;
                        n->lo  = lo[d];
                        n->hi  = hi[d];

                        // bisection
                        SCALAR me = 0.5 * ( n->lo + n->hi );
                        SIZET below = l, above = r;
                        for ( SIZET i = l; i < r; i++ ) {
                            if ( points[ idx_[i] ][d] < me )
                                idx_.shadow( below++ ) = idx_[i];
                            else
                                idx_.shadow( --above ) = idx_[i];
                        }
                        SIZET m = below;

                        // determine bounce and spawn child nodes
                        for ( SIZET k = l; k < m; k++ ) n->lim[0] = std::max( n->lim[0], points[idx_.shadow( k )][d] );
                        for ( SIZET k = m; k < r; k++ ) n->lim[1] = std::min( n->lim[1], points[idx_.shadow( k )][d] );

                        n->child = allocator[ omp_get_thread_num() ].allocate( 2 );
                        n->child[0].parent = n->child[1].parent = n;
                        stack_new[ gid * 2 + 0 ] = std::make_tuple( n->child + 0, l, m );
                        stack_new[ gid * 2 + 1 ] = std::make_tuple( n->child + 1, m, r );
                    }
                }

                idx_.swap();
                stack.swap( stack_new );
                stack_new.resize( stack.size() * 2 );
            }

            #pragma omp parallel for
            for ( SIZET i = 0; i < stack.size(); i++ ) build_serial( stack[i], points );
        }

        // build continuous storage
        #pragma omp parallel for
        for ( SIZET i = 0; i < npts_; i++ ) pts_[i] = points[ idx_[i] ];

        // build shortcut to the next branch
        visit( &root, [&]( node * n ) {
            while ( n->parent && n == n->parent->child + 1 ) n = n->parent; // go up until being a left child
            if ( n->parent ) return n->parent->child + 1;
            else return ( node * )nullptr;
        } );

        Service<Timers>::call()["|KDTree::build"].stop();
    }

    std::tuple<SCALAR, SIZET> find_nearest( point const & q, std::tuple<SCALAR, SIZET> ig = std::make_tuple( std::numeric_limits<SCALAR>::max(), -1 ) ) const {
        SCALAR min_r, r;
        SIZET min_i, i;

        std::tie( min_r, min_i ) = ig;
        if ( min_r == std::numeric_limits<SCALAR>::max() ) {
            SIZET x = rand() % npts_;
            min_i = idx_[ x ];
            min_r = norm( pts_[ x ] - q );
        }
        node const * p = &root;
        do {
            while ( p->child && in_box( q, p, min_r ) ) { // descend to the nearest leave
                if ( q[p->d] - min_r < p->lim[0] ) p = p->child + 0; // trying left tree
                else if ( q[p->d] + min_r > p->lim[1] ) p = p->child + 1; // tryign right tree
                else break; // both end no worth trying
            }
            if ( !p->child && in_box( q, p, min_r ) ) {
                std::tie( r, i ) = nearest( p->plr[0], p->plr[1], q );
                if ( r < min_r ) {
                    min_r = r;
                    min_i = i;
                }
            }
            do { // march to the next unsearched branch
                p = p->next;
            } while ( p && q[p->parent->d] + min_r <= p->parent->lim[1] );
        } while ( p ); // exit search when back to root

        return std::make_tuple( min_r, min_i );
    }

    // return first nmax points within max_r
    SIZET find_within( point const & q, SCALAR max_r, SIZET * __restrict result, SIZET nmax ) const {
        node const * p = &root;
        SIZET nfound = 0;
        do {
            while ( p->child && in_box( q, p, max_r ) ) { // descend to the nearest leave
                if ( q[p->d] - max_r < p->lim[0] ) p = p->child + 0; // trying left tree
                else if ( q[p->d] + max_r > p->lim[1] ) p = p->child + 1; // tryign right tree
                else break; // both end no worth trying
            }
            if ( !p->child && in_box( q, p, max_r ) ) {
                for ( SIZET i = p->plr[0]; i < p->plr[1]; i++ ) {
                    if ( normsq( pts_[ i ] - q ) < max_r * max_r ) {
                        result[nfound++] = idx_[i];
                        if ( nfound == nmax ) return nfound;
                    }
                }
            }
            do { // march to the next unsearched branch
                p = p->next;
            } while ( p && q[p->parent->d] + max_r <= p->parent->lim[1] );
        } while ( p ); // exit search when back to root
        return nfound;
    }

    SIZET find_within( point const & q, SCALAR max_r, AlignedArray<SIZET, true> & result ) const {
        node const * p = &root;
        //SIZET nfound = 0;
        result.resize( 0 );
        do {
            while ( p->child && in_box( q, p, max_r ) ) { // descend to the nearest leave
                if ( q[p->d] - max_r < p->lim[0] ) p = p->child + 0; // trying left tree
                else if ( q[p->d] + max_r > p->lim[1] ) p = p->child + 1; // tryign right tree
                else break; // both end no worth trying
            }
            if ( !p->child && in_box( q, p, max_r ) ) {
                for ( SIZET i = p->plr[0]; i < p->plr[1]; i++ ) {
                    if ( normsq( pts_[ i ] - q ) < max_r * max_r ) {
                        result.emplace_back( idx_[i] );
                    }
                }
            }
            do { // march to the next unsearched branch
                p = p->next;
            } while ( p && q[p->parent->d] + max_r <= p->parent->lim[1] );
        } while ( p ); // exit search when back to root
        return result.size();
    }

    void clear() {
        npts_ = 0;
        root.child = root.parent = root.next = nullptr;
    }

protected:
    struct alignas( config::cacheline ) NodeAllocator {
        std::vector<node *> ptrs_free, ptrs_used;
        node * lead = nullptr;
        SIZET stock = 0;
        const static SIZET chunk_size = 2048;

        ~NodeAllocator() {
            for ( auto & p : ptrs_used ) delete [] p;
            for ( auto & p : ptrs_free ) delete [] p;
        }

        node * allocate( SIZET n ) {
            if ( stock < n ) {
                if ( lead ) ptrs_used.push_back( lead );
                if ( !ptrs_free.size() ) {
                    lead = ( node * ) aligned_allocate( config::cacheline, sizeof( node ) * chunk_size );
                } else {
                    lead = ptrs_free.back();
                    ptrs_free.pop_back();
                }
                stock = chunk_size;
            }
            node * ret = new( lead + chunk_size - stock ) node[ n ];
            stock -= n;
            return ret;
        }

        void refresh() {
            for ( auto & p : ptrs_used ) ptrs_free.push_back( p );
            ptrs_used.clear();
            if ( lead ) ptrs_free.push_back( lead );
            lead = nullptr;
            stock = 0;
        }
    };

    std::vector<NodeAllocator> allocator;

    inline bool in_box( point const & q, node const * p, SCALAR min_r ) const {
        return ( q[p->d] - min_r < p->hi ) & ( q[p->d] + min_r > p->lo ); // non-short-circuit version faster by saving branches
    }

    // small-range brute-force search
    std::tuple<SCALAR, SIZET> nearest( SIZET l, SIZET r, point q ) const {
        SIZET p;
        auto r2min = std::numeric_limits<SCALAR>::max();
        for ( SIZET i = l; i < r; i++ ) {
            auto r2 = normsq( pts_[ i ] - q );
            if ( r2 < r2min ) {
                p = i;
                r2min = r2;
            }
        }
        return std::make_tuple( std::sqrt( r2min ), idx_[p] );
    }

    template<class VISITOR> void visit( node * n, VISITOR const & visitor ) {
        if ( n->child ) {
            visit( n->child + 0, visitor );
            visit( n->child + 1, visitor );
        }
        n->next = visitor( n );
    }

    static inline uint rand() {
        static uint seed = 0;
        return ( seed = seed * 1664525U + 1013904223U );
    }
};

}

#endif

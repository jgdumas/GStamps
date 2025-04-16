// ==========================================================================
// GStamps: C++ routines for the Global Postage Stamp Problem
// Authors: J-G. Dumas, A. Galan, B. Grenet, A. Maignan, D. S. Roche
// ==========================================================================

/****************************************************************
 * GStamps Library
 ****************************************************************/

#ifndef _GSTAMPS_LIBRARY_H_
#define _GSTAMPS_LIBRARY_H_

#include <iostream>
#include <stdlib.h>
#include <cstdint>
#include <cassert>
#include <vector>
#include <set>
#include <map>
#include <numeric>
#include <boost/dynamic_bitset.hpp>


// // Using Givaro::Integer arbitrary precision, or int64_t on 64 bits

#ifdef __GSTAMPS_EXTENDED_PRECISION

#  include <givaro/givtimer.h>
#  include <givaro/givinteger.h>
using StTimer=Givaro::Timer;
using bint=Givaro::Integer;
#  define __St_Zero bint::zero
#  define __St_One bint::one

#else

#  include <gstamps_chrono.h>
using StTimer=Chrono;
using bint=int64_t;
#  define __St_Zero 0l
#  define __St_One 1l

#endif

// Maximal RAM usage
#ifndef __GSTAMPS_MAXCOVER
#define __GSTAMPS_MAXCOVER 137438953472
#endif


// ============================================
// Printing Tools

// log[2](v)
uint64_t MSB(uint64_t v) { uint64_t r(0); for( ; v>>=1; ++r); return r; }

template<typename Cont>
std::ostream& rangeprint(std::ostream& out, const Cont& v);

template<typename Cont>
typename Cont::value_type firstrange(const Cont& v);

template<typename Cont>
std::ostream& firstrangeprint(std::ostream& out,
                              const Cont& v, const bint& max);


// ============================================
// Cover: Local Postage Stamp Problem

        // Low memory
template<typename List>
bint ICover(const List& points, const size_t s, const int verbose = 0);

        // Fast
template<typename Iterator>
bint _Cover(const Iterator& start, const Iterator& end, const size_t s);

        // Fast with checks
template<typename List>
bint Cover(const List& points, const size_t s, const int verbose);

        // Stores only reached values, not whole array
template<typename List>
bint Reach(List& points, const size_t h, const int verbose = 0);


// ============================================
// Exhaustive search

template<typename List>
bint FixedPoints(List& pmax,
                 const List& points, const size_t s, const size_t i=0);

template<typename List>
bint BruteForce(List& points, size_t k, size_t s, const int verbose=0);


// ============================================
// Special cases

        // Every two Fibonnaci numbers
bint Fibonacci(std::vector<bint>& points, const size_t k,
               const int verbose = 0);


        // S. Mossige Algorithms for Computing the $h$-Range
        // of the Postage Stamp Problem, Math. Comp 1981
        // Appendix by Torleiv Kltve and Svein Mossige
        // Error: n2 is not (4*k*x-14*x*x+12*x-4)
        // it is (4*k*x-14*x*x+12*x-6)
template<typename List>
bint KloveMossige(List& points, const size_t k, const int verbose = 0);


        // Alter, R., & Barnett, J. A. (1977).
        // Remarks on the postage stamp problem with applications to computers.
        // University of Kentucky. Department of Computer Science.
        // https://notatt.com/lsu-stamp.pdf
bint AlterBernett(std::vector<bint>& points, const size_t k, const size_t s,
                  const int verbose = 0);


        // AlterBarnett improvement with balanced ranges
bint BalGreedy(std::vector<bint>& points, const size_t k, const size_t s,
               const int verbose = 0);


// ============================================
// Known extremal cases
// Challis, Michael F, et John P Robinson. 2010.
// Some Extremal Postage Stamp Bases.
// https://cs.uwaterloo.ca/journals/JIS/VOL13/Challis/challis6.html

bint kThree(std::vector<bint>& points, const size_t s, const int verbose = 0);

bint kFour(std::vector<bint>& points, const size_t s, const int verbose = 0);

bint kFive(std::vector<bint>& points, const size_t s, const int verbose = 0);

bint kSix(std::vector<bint>& points, const size_t s, const int verbose = 0);

bint kSeven(std::vector<bint>& points, const size_t s, const int verbose = 0);

bint kEight(std::vector<bint>& points, const size_t s, const int verbose = 0);

bint sTwo(std::vector<bint>& points, const size_t k, const int verbose = 0);

bint sThree(std::vector<bint>& points, const size_t k, const int verbose = 0);

bint sFour(std::vector<bint>& points, const size_t k, const int verbose = 0);

bint sFive(std::vector<bint>& points, const size_t k, const int verbose = 0);

bint sSix(std::vector<bint>& points, const size_t k, const int verbose = 0);


// ============================================
// Precomputed base cases, but not extremal

#include <gstamps_basis.h>



// ============================================
// Hybrid methods

// Generic algorithm, calling all methods, with memoization of solutions
bint FSelect(std::vector<bint>& points, const size_t k, const size_t s,
             const int rlevel = 0, const int verbose = 0);

// Mrose Divide & Conquer
template<typename List>
bint CutSelect(List& points, const size_t k, const size_t kotwo,
               const size_t s, const size_t sotwo,
               const int rlevel = 0, const int verbose = 0);

// Recursive (rlevel) Quadratic exploration, or midpoints only
template<typename List>
bint RecSelect(List& points, const size_t k, const size_t s,
               const int rlevel = 0, const int verbose = 0);

// Switching among different solutions, known extremal first
bint DSelect(std::vector<bint>& points, const size_t k, const size_t s,
             const int rlevel = 0, const int verbose = 0);


// ============================================
// Exhausting from the end

// Exhaust additional denominations between prescribed and (k-prescribed.size())
bint complement(std::vector<bint>& prescribed,
                const size_t k, const size_t s, const int verbose = 0);

// Exhaust additional denominations, parallel version (first level)
bint par_complement(std::vector<bint>& prescribed,
                    const size_t k, const size_t s, const int verbose = 0);


// ============================================
// Inline implementations

#include <gstamps.inl>


// ============================================

#endif

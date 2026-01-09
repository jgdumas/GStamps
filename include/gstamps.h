// ==========================================================================
// GStamps: C++ routines for the Local & the Global Postage Stamp Problem
// Authors:
//   L. Colisson, J-G. Dumas, A. Galan, B. Grenet, A. Maignan, D. S. Roche
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

// ============================================
// By default, early terminate LPSP
#ifndef __GSTAMPS_SELMER_LEMMA
#define __GSTAMPS_SELMER_LEMMA 1
#endif

// ============================================
// Choosing precision:
//   using Givaro::Integer arbitrary precision, or int64_t on 64 bits
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

// ============================================
// Printing Tools

template<typename Cont>
std::ostream& ScopePrint(std::ostream& out, const Cont& v);

// ============================================
// Range: Local Postage Stamp Problem

        // Fast small s
template<typename Iterator, typename stype_t>
bint _SRange(const Iterator& start, const Iterator& end, const stype_t s);

        // Fast large s
template<typename Iterator, typename stype_t>
bint _KRange(const Iterator& start, const Iterator& end, const stype_t s);

        // Fast with checks
template<typename List, typename stype_t>
bint Range(const List& points, const stype_t s, const int verbose);

// ============================================
// Exhaustive search

template<typename List, typename stype_t>
bint FixedPoints(List& pmax,
                 const List& points, const stype_t s, const size_t i=0);

template<typename List, typename stype_t>
bint BruteForce(List& points, size_t k, stype_t s, const int verbose=0);


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
template<typename stype_t>
bint AlterBarnett(std::vector<bint>& points, const size_t k, const stype_t s,
                  const int verbose = 0);


        // AlterBarnett improvement with balanced ranges
template<typename stype_t>
bint BalGreedy(std::vector<bint>& points, const size_t k, const stype_t s,
               const int verbose = 0);


// ============================================
// Known extremal cases
// Challis, Michael F, et John P Robinson. 2010.
// Some Extremal Postage Stamp Bases.
// https://cs.uwaterloo.ca/journals/JIS/VOL13/Challis/challis6.html

bint kThree(std::vector<bint>& points, const size_t s, const int verbose = 0);

//  approx:	if false the range is computed (via an exponential algorithm)
//			otherwise only a lower bound on the range is returned
bint kFour(std::vector<bint>& points, const size_t s, const bool approx,
           const int verbose = 0);

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
// Hybrid methods: produce a k-basis in points for s-stamps
//	rlevel:	D&C splits (k,s) into (k1,s1) and (k-k1,s-s1) then
//			for (rlevel) recursive levels: (k-1)*(s-1) possibilites explored
//			then; midpoint (k1=k/2,s1=s/2) is chosen.
//  approx:	if false the range is computed (via an exponential algorithm)
//			otherwise only a lower bound on the range is returned
//  verbose:controls the level of verbosity


// Tries memoization of solutions or calls DSelect
template<typename stype_t>
bint FSelect(std::vector<bint>& points, const size_t k, const stype_t s,
             const int rlevel = 0, const bool approx=false,
             const int verbose = 0);

// Generic algorithm, calling all methods, known extremal first
template<typename stype_t>
bint DSelect(std::vector<bint>& points, const size_t k, const stype_t s,
             const int rlevel = 0, const bool approx=false,
             const int verbose = 0);

// Mrose Recursive Divide & Conquer
//   (rlevel) Quadratic exploration of CutSelect, or midpoints only
template<typename List, typename stype_t>
bint RecSelect(List& points, const size_t k, const stype_t s,
               const int rlevel = 0, const bool approx=false,
               const int verbose = 0);

// Divide & Conquer splitting (k,s) into (k1,s1) and (k-k1,s-s1)
template<typename List, typename stype_t>
bint CutSelect(List& points, const size_t k, const size_t k1,
               const stype_t s, const stype_t s1,
               const int rlevel = 0, const bool approx=false,
               const int verbose = 0);


// ============================================
// Exhausting from the end

// Exhaust additional denominations between prescribed and (k-prescribed.size())
template<typename stype_t>
bint complement(std::vector<bint>& prescribed,
                const size_t k, const stype_t s, const int verbose = 0);

// Exhaust additional denominations, parallel version (first level)
template<typename stype_t>
bint par_complement(std::vector<bint>& prescribed,
                    const size_t k, const stype_t s, const int verbose = 0);


// ============================================
// Inline implementations

#include <gstamps_print.inl>
#include <gstamps_range.inl>
#include <gstamps_extreme.inl>
#include <gstamps_special.inl>
#include <gstamps_exhaust.inl>
#include <gstamps_recursive.inl>

// ============================================

#endif

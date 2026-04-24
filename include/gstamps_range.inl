// ==========================================================================
// GStamps: C++ routines for the Local & the Global Postage Stamp Problem
// Authors:
//   L. Colisson, J-G. Dumas, A. Galan, B. Grenet, A. Maignan, D. S. Roche
// ==========================================================================

/****************************************************************
 * GStamps Library, LPSP inline implementations
 ****************************************************************/

// ============================================
// Masking Tools
// upmask: Round up to the next (highest power of 2, minus 1) of (input+1)

#ifdef __GSTAMPS_EXTENDED_PRECISION
Givaro::Integer upmask(const Givaro::Integer& w) {
    Givaro::Integer v(w);
    uint32_t exp(1);
    for(Givaro::Integer shi(1); shi>0; exp <<=1) {
        shi = v >> exp;
        v |= shi;
    }
    return std::move(v);
}

#else
// See: https://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
uint64_t upmask(const uint64_t& w) {
    bint v(w);
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v |= v >> 32;
    return std::move(v);
}
#endif



// ============================================
// Range: local postage stamp problem

// Loop from 1 to s, with a binary vector
template<typename Iterator, typename stype_t>
inline bint _SRange(const Iterator& start, const Iterator& end,
                    const stype_t s, const int verbose) {
	// Binary range
    const bint& back(*std::prev(end));				// k>=1
    if (back == __St_One) return s;
    bint vs(back);
    const bint upper(s*vs+1);
    boost::dynamic_bitset<> reached(upper,false);
    reached[0]=true;                                // 0 reached
    for(auto it=start; it!=end; ++it) reached[*it]=true;   // points reached

#if __GSTAMPS_SELMER_LEMMA
    const bint& penult(*std::prev(std::prev(end))); // back>1 => k>=2
    bint lb(penult-back);
#endif

    bint notin(1);
    for(stype_t d(1); d<s; ++d, vs += back) {
#if __GSTAMPS_SELMER_LEMMA
        notin = vs+1;
#endif
        for(size_t i=vs; i>=d; --i) {
            if(reached[i]) {
                for(auto right=start; right!=end; ++right) {
                    reached[i+(*right)]=true;
                }
            }
#if __GSTAMPS_SELMER_LEMMA
            else {
                notin=i;
            }
#endif
        }

#if __GSTAMPS_SELMER_LEMMA
            // Test Selmer's Lemma for early termination
        lb += penult;
        if ( (notin > back) && (notin > lb) ) {
                    // Range will now surely attain c+(s-d)ak
            if (verbose>0) std::clog << "#[ET(" << (size_t)d << '|' << lb
                                     << ")]: " << notin << " --> "
                                     << (notin-1+(s-d)*back) << std::endl;
            return --notin += (s-d)*back;
        }
#endif
    }

    size_t max(reached.size()-1);
    for(size_t jr(notin); jr<reached.size(); ++jr){
        if (! reached[jr]) {
            max = jr-1;
            break;
        }
    }

    return max;
}


// Loop from 1 to n
template<typename List, typename stype_t>
inline bint _KRange(const List& points, const size_t k, const stype_t s,
                    const int verbose) {
    assert( (k>=1) && (k<=points.size()) );
    const auto& back(points.back());				// k>=1
    if (back == __St_One) return s;

    const size_t window(upmask(back));	// highest 1-full mask gt
    const stype_t spu(s+1);				// s+1 is unreachable
    std::vector<std::pair<stype_t,stype_t>> reached(window+1u,
                                                    std::make_pair(spu,0u));

    for(stype_t i=0; i<k; ++i)
        reached[points[i]]=std::make_pair(1u,i);

#if __GSTAMPS_SELMER_LEMMA > 1
    const auto& penult(points[k-2]); // back>1 => k>=2
    std::vector<size_t> selmer(s+1);
    std::iota(selmer.begin(), selmer.end(), 1);
    stype_t maxs(0), mins(0);
    for(auto& is: selmer) {
        is *= penult;
        if (is <= back) ++mins;
        is -= back;
    }
#endif

    size_t index(1);
    for(; reached[index & window].first<=s; ++index) {
        auto& slocal(reached[index & window]);
        const stype_t slfirst(slocal.first);
        const stype_t vlocal(slfirst+1u);
        for(auto right=slocal.second; right<k; ++right) {
            auto& starget(reached[ (index+points[right]) & window]);
            if (starget.first > vlocal) {
                starget.first = vlocal;
                starget.second = right;
            }
        }

#if __GSTAMPS_SELMER_LEMMA > 1
        maxs = (slfirst>maxs?slfirst:maxs);
        if ( (index & 1048575u) == 1048575u) { // Reduce overhead
                // Selmer's lemma
            if ((maxs>mins) && (maxs<s) && (index > selmer[maxs])) {
                    // Find maxs_range
                for(size_t i=1; i<=back; ++i) {
                        // Complete s_range
                    if (reached[(index+i) & window].first>maxs) {
                        if (verbose>0) std::clog << "#[ET(" << (size_t)maxs
                                                 << '|' << selmer[maxs] << ")]: "
                                                 << i << " -> "
                                                 << (index+i-1+(s-maxs)*back)
                                                 << std::endl;
                        index += i;
                        return --index += (s-maxs)*back;
                    }
                }
            }
        }
#endif
        slocal.first=spu;		// clean up sliding window
    }

    return --index;

}

template<typename List, typename stype_t>
inline bint _KRange(const List& points, const stype_t s,
                    const int verbose) {
    return _KRange(points,points.size(),s,verbose);
}


template<typename List, typename stype_t>
inline bint Range(const List& points, const stype_t s, const int verbose) {
    if (verbose>1)
        ScopePrint(std::clog << "#[Range] Basis: ", points) << std::endl;

    StTimer chrono; chrono.start();
    const bint max( (s<6u) ?
                    _SRange(points.begin(), points.end(), s, verbose) :
                    _KRange(points, s, verbose)
                    );
    chrono.stop();

    if (verbose>1) {
        std::clog << "#[Range(" << size_t(s) << ")]: 1.." << max
                  << " ..." << std::endl;
    }

    if (verbose>0) std::clog << "#[Range(" << size_t(s) << ")]: " << max
                             << ' ' << chrono <<std::endl;
    return max;
}




// ==========================================================================
// Reads a basis from std::cin
// Computes the range and decompositions that basis with argv[1] stamps
// Modification of the _KRange algorithm without the sliding window
template<typename List, typename stype_t>
std::vector<PsVs<stype_t>> _Decompose(const List& points, const stype_t s,
                                      const int verbose) {
    const auto& back(points.back());				// k>=1
    const stype_t spu(s+1);				// s+1 is unreachable

    if (back == __St_One) {
        std::vector<PsVs<stype_t>> reached(spu);
        for(size_t i(1); i<spu; ++i) {
            reached[i]={0,std::vector<stype_t>(i,0)};
        }
        return reached;
    }

    const PsVs<stype_t> spustart{spu, {0u}};
        // Maximal valid index is s*back
        //   thus maximal tested in loop is s*back+1
        //   and maximal starget is at s*back+back=spu*back>=s*back+1
        //   with 0 indexing this gives a table of size: spu*back+1
    std::vector<PsVs<stype_t>> reached(spu*back+1,spustart);

    for(stype_t i=0; i<points.size(); ++i)
        reached[points[i]]= PsVs<stype_t>{1u,{i}};

    size_t index(1);
    for(; reached[index].first<=s; ++index) {
        const auto& slocal(reached[index]);
        const stype_t slfirst(slocal.first);
        const stype_t vlocal(slfirst+1u);
        for(auto right=slocal.second.back(); right<points.size(); ++right) {
            auto& starget(reached[index+points[right]]);
            if (starget.first>vlocal) {
                starget.first = vlocal;
                starget.second.resize(0);
                starget.second.assign(slocal.second.begin(),
                                      slocal.second.end());
                starget.second.push_back(right);
            }
        }

    }

    reached.resize(index); // range is --index

    return reached;
}
// ==========================================================================

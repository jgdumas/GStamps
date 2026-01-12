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
template<typename Iterator, typename stype_t>
inline bint _KRange(const Iterator& start, const Iterator& end,
                    const stype_t s, const int verbose) {
    const bint& back(*std::prev(end));	// k>=1
    if (back == __St_One) return s;

    const bint window(upmask(back));	// highest 1-full mask gt
    const stype_t spu(s+1);				// s+1 is unreachable
    std::vector<stype_t> reached(window+1u, spu);

    for(auto it=start; it!=end; ++it) reached[*it]=1u;

#if __GSTAMPS_SELMER_LEMMA > 1
    const bint& penult(*std::prev(std::prev(end))); // back>1 => k>=2
    std::vector<bint> selmer(s+1);
    std::iota(selmer.begin(), selmer.end(), 1);
    stype_t maxs(0), mins(0);
    for(auto& is: selmer) {
        is *= penult;
        is -= back;
        if (is <= back) {
            ++mins;
        }
    }
#endif

    bint index(1);
    for(; reached[index & window]<=s; ++index) {
        stype_t& slocal(reached[index & window]);
        const stype_t vlocal(slocal+1u);
        for(auto it=start; it!=end; ++it) {
            stype_t& starget(reached[ (index+(*it)) & window]);
            if (starget>vlocal) {
                starget = vlocal;
            }
        }

#if __GSTAMPS_SELMER_LEMMA > 1
            // Selmer's lemma
        maxs = (slocal>maxs?slocal:maxs);
        if ((maxs>mins) && (maxs<s) && (index > selmer[maxs])) {
                // Find maxs_range
            for(bint i=1; i<=back; ++i) {
                // Complete s_range
                if (reached[(index+i) & window]>maxs) {
                    if (verbose>0) std::clog << "#[ET(" << (size_t)maxs << '|'
                                             << selmer[maxs] << ")]: " << i
                                             << " -> "
                                             << (index+i-1+(s-maxs)*back)
                                             << std::endl;
                    index += i;
                    return --index += (s-maxs)*back;
                }
            }
        }
#endif
        slocal=spu;		// clean up sliding window
    }

    return --index;

}

template<typename List, typename stype_t>
inline bint Range(const List& points, const stype_t s, const int verbose) {
    if (verbose>1)
        ScopePrint(std::clog << "#[Range] Basis: ", points) << std::endl;

    StTimer chrono; chrono.start();
    const bint max( (s<6u) ?
                    _SRange(points.begin(), points.end(), s, verbose) :
                    _KRange(points.begin(), points.end(), s, verbose)
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

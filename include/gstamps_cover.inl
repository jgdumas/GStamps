// ==========================================================================
// GStamps: C++ routines for the Global Postage Stamp Problem
// Authors: 
//   L. Colisson, J-G. Dumas, A. Galan, B. Grenet, A. Maignan, D. S. Roche
// ==========================================================================

/****************************************************************
 * GStamps Library, LPSP inline implementations
 ****************************************************************/

// ============================================
// Cover: local postage stamp problem

// Loop from 1 to s, with a binary vector
template<typename Iterator, typename stype_t>
inline bint _SCover(const Iterator& start, const Iterator& end, const stype_t s) {
	// Binary cover
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
                    // Cover will now surely attain c+(s-d)ak
//                 std::clog << "#[Selmer(" << (size_t)d << '|' << lb << ")]: "
//                           << notin << " --> " << (notin-1+(s-d)*back)
//                           << std::endl;
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
inline bint _KCover(const Iterator& start, const Iterator& end,
                    const stype_t s) {
    const bint& back(*std::prev(end));	// k>=1
    if (back == __St_One) return s;

    const bint window(upmask(back));	// highest 1-full mask gt
    std::vector<stype_t> reached(window+1u, 0u);

    for(auto it=start; it!=end; ++it) reached[*it]=1u;

#if __GSTAMPS_SELMER_LEMMA
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
            const stype_t vtarget(starget);
            if ( (vtarget == 0u) || (vtarget>vlocal) ) {
                starget = vlocal;
            }
        }

#if __GSTAMPS_SELMER_LEMMA
            // Selmer's lemma
        maxs = (slocal>maxs?slocal:maxs);
        slocal=0u;		// clean up sliding window
        if ((maxs>mins) && (maxs<s) && (index > selmer[maxs])) {
                // Find maxs_range
            stype_t vlocal(0u);
            do {
                ++index;
                vlocal = reached[index & window];
            } while( (vlocal <= maxs) && (vlocal > 0) );
                // Complete s_range
//             std::clog << "#[Selmer(" << (size_t)maxs << '|' << selmer[maxs]
//                       << ")]: " << index << " -> " << (index-1+(s-maxs)*back)
//                       << std::endl;
            return --index += (s-maxs)*back;
        }
#else
        slocal=0u;		// clean up sliding window
#endif
    }

    return --index;

}

template<typename List, typename stype_t>
inline bint Cover(const List& points, const stype_t s, const int verbose) {
    if (verbose>1)
        rangeprint(std::clog << "#[Cover] Basis: ", points) << std::endl;

    StTimer chrono; chrono.start();
    const bint max( (s<6u) ? _SCover(points.begin(), points.end(), s)
                    : _KCover(points.begin(), points.end(), s)
                    );
    chrono.stop();

    if (verbose>1) {
        std::clog << "#[Cover(" << size_t(s) << ")]: 1.." << max
                  << " ..." << std::endl;
    }

    if (verbose>0) std::clog << "#[Cover(" << size_t(s) << ")]: " << max
                             << ' ' << chrono <<std::endl;
    return max;
}

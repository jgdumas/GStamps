// ==========================================================================
// GStamps: C++ routines for the Global Postage Stamp Problem
// Authors:
//   L. Colisson, J-G. Dumas, A. Galan, B. Grenet, A. Maignan, D. S. Roche
// ==========================================================================

/****************************************************************
 * GStamps Library, inline implementations
 ****************************************************************/


// ============================================
// Exhaustive search

template<typename List, typename stype_t>
inline bint FixedPoints(List& pmax,
			const List& points, const stype_t s, const size_t i) {
    pmax.resize(i); pmax.reserve(points.size());
    for(auto it=points.begin()+i; it!=points.end(); ++it) {
	pmax.push_back( _KCover(points.begin(), std::next(it), s) );
    }
    return pmax.back();
}



template<typename List, typename stype_t>
inline bint BruteForce(List& points, size_t k, stype_t s, const int verbose) {
    assert(k>0);
    if (k == 1) {
        points = {1};
        return s;
    }
    points.resize(k); std::iota(points.begin(), points.end(), 1);
    List pointsmax(points);

    List covsmax;
    bint max( FixedPoints(covsmax, points, s) );

        // Next set of points
    while (true){
            // Find where next increment is needed
        size_t i=points.size()-1;
        while ((i > 0) && (points[i] == (covsmax[i-1]+1))) { --i; }

            // Nothing left
        if (i == 0) {
            points.assign(pointsmax.begin(), pointsmax.end());
            return max;
        }

            // Increment is needed points[i] and reset all higher denominations
        std::iota(points.begin()+i,points.end(),points[i]+1);

            // Compute new maxima for the denominations after i
        const bint max2 = FixedPoints(covsmax, points, s, i-1);

            // Save new points if better
        if (max2 > max) {
            if (verbose>0) {
                std::clog << "#[Brute] " << max2 << ": ";
                for(const auto& it: points) std::clog << it << ' ';
                std::clog << std::endl;
            }

            max = max2;
            pointsmax.assign(points.begin(), points.end());
        }
    }
}




// ============================================
// Exhaust additional denominations ...

#ifndef __GSTAMPS_Half_search__
#  define __GSTAMPS_AMX(a,p) (a+__St_One)
#else
	// Search only within the largest half
#  define __GSTAMPS_AMX(a,p) ( ((p+a)>>1) + __St_One )
#endif

template<typename stype_t>
inline bint complement(std::vector<bint>& prescribed,
                       const size_t k, const stype_t s, const int verbose) {
    assert(prescribed.size()>=1);
    assert(prescribed.size()<=k);

    bint bmax(_KCover(prescribed.begin(),prescribed.end(),s));
    const bint pc(bmax+__St_One);
    const bint amx( __GSTAMPS_AMX(prescribed.back(), bmax) );
    if (verbose>0) {
        std::clog << "#[Cpmt(" << prescribed.size() << ")] amx: " << amx
                  << " bm: " << bmax << " pc: " << pc << " with prescribed: ";
        for(const auto& it: prescribed) std::clog << it << ' ';
        std::clog << std::endl;
    }

    if ( (prescribed.size()>=k) || (amx>pc) ) return bmax;

    std::vector<bint> points; points.reserve(k);
    points.assign(prescribed.begin(), prescribed.end());

    std::vector<bint> bfound; bfound.reserve(k);

    for(bint u(pc); u>=amx; --u) {
        points.push_back(u);
        const bint bu = complement(points, k, s, verbose-1);
        if (bu > bmax) {
            bfound.resize(0);
            bfound.assign(points.begin(), points.end());
            bmax = bu;
            if (verbose>0) {
                std::clog << "#[Cpmt(" << prescribed.size() << ")] max: "
                          << bmax << " with basis: ";
                for(const auto& it: points) std::clog << it << ' ';
                std::clog << std::endl;
            }
        }
        points.resize(prescribed.size());
        if (verbose>0) std::clog << "#[Cpmt(" << prescribed.size() << ")] "
                                 << amx << " <= " << u << " <= " << pc
                                 << " : " << bu << " <= " << bmax << std::endl;
    }

    prescribed.resize(0);
    prescribed.assign(bfound.begin(),bfound.end());

    return bmax;
}


// Exhaust additional denominations ...
template<typename stype_t>
inline bint par_complement(std::vector<bint>& prescribed,
                           const size_t k, const stype_t s, const int verbose) {

    assert(prescribed.size()>=1);
    assert(prescribed.size()<k);

    bint bmax(_KCover(prescribed.begin(),prescribed.end(),s));
    const bint pc(bmax + __St_One);
    const bint amx( __GSTAMPS_AMX(prescribed.back(), pc) );

    std::clog << "#[PCt(" << prescribed.size() << ")]"
              << " amx: " << amx << " pc: " << pc << std::endl;

    std::vector<bint> bfound; bfound.reserve(k);

    const int64_t maxu(pc-amx); // Should not loop more than 2^63 anyway ...

    if (maxu>=0) {

#pragma omp parallel for shared(prescribed,bfound,bmax,amx,pc,k,s,verbose) schedule(dynamic)
        for(int64_t iu=maxu; iu>=0; --iu) {
            std::vector<bint> points; points.reserve(k);
            points.assign(prescribed.begin(), prescribed.end());
            points.push_back(amx+bint(iu));
            const bint bu = complement(points, k, s, verbose-1);
#pragma omp critical
            {
                if (bu > bmax) {
                    bfound.resize(0);
                    bfound.assign(points.begin(), points.end());
                    bmax = bu;
                    if (verbose>0) {
                        std::clog << "#[PCt(" << prescribed.size()
                                  << ")] max: " << bmax << " with basis: ";
                        for(const auto& it: points) std::clog << it << ' ';
                        std::clog << std::endl;
                    }
                }
                if (verbose>0)
                    std::clog << "#[PCt(" << prescribed.size() << ")] " << amx
                              << " <= " << points[prescribed.size()] << " <= "
                              << pc << " : " << bu << " <= " << bmax
                              << std::endl;
            }
        }

        prescribed.resize(0);
        prescribed.assign(bfound.begin(),bfound.end());
    }

    return bmax;
}

// ==========================================================================
// GStamps: C++ routines for the Local & the Global Postage Stamp Problem
// Authors:
//   L. Colisson, J-G. Dumas, A. Galan, B. Grenet, A. Maignan, D. S. Roche
// ==========================================================================

/**********************************************************************
 * GStamps: Local Postage Stamp Problem
 **********************************************************************/

#include <gstamps.h>

// Reads a basis from std::cin
// Computes the range of that basis with argv[1] stamps
// Using _KRange algorithm

// Loop from 1 to n
template<typename Iterator, typename stype_t>
inline bint _ET_KRange(const Iterator& start, const Iterator& end,
                    const stype_t s, const int verbose) {
    const bint& back(*std::prev(end));	// k>=1
    if (back == __St_One) return s;

    const bint window(upmask(back));	// highest 1-full mask gt
    const stype_t spu(s+1);				// s+1 is unreachable
    const auto spustart(std::make_pair(spu,start));
    std::vector<std::pair<stype_t,Iterator>> reached(window+1u, spustart);

    for(auto it=start; it!=end; ++it) reached[*it]=std::make_pair(1u,it);

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
    for(; reached[index & window].first<=s; ++index) {
        auto& slocal(reached[index & window]);
        const stype_t slfirst(slocal.first);
        const stype_t vlocal(slfirst+1u);
        for(auto right=slocal.second; right!=end; ++right) {
            auto& starget(reached[ (index+(*right)) & window]);
            if (starget.first>vlocal) {
                starget.first = vlocal;
                starget.second = right;
            }
        }

#if __GSTAMPS_SELMER_LEMMA > 1
            // Selmer's lemma
        maxs = (slfirst>maxs?slfirst:maxs);
        if ((maxs>mins) && (maxs<s) && (index > selmer[maxs])) {
                // Find maxs_range
            for(bint i=1; i<=back; ++i) {
                // Complete s_range
                if (reached[(index+i) & window].first>maxs) {
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

        slocal=spustart;		// clean up sliding window
    }

    return --index;

}

template<typename stype_t>
int tmain(int argc, char **argv, stype_t s) {
    const int verbose(argc>2?atoi(argv[2]):0);

    std::vector<bint> points;
    bint tmp;
    std::clog << "# Enter a basis (end by non-digit) ... " << std::endl;
    while(std::cin >> tmp) {
        points.push_back(tmp);
    };

    ScopePrint(std::clog << "# Basis: ", points) << std::endl;

    StTimer chrono; chrono.start();
    const bint max( _ET_KRange(points.begin(), points.end(), s, verbose) );
    chrono.stop();

    if (verbose>0) std::clog << "#[WRange(" << size_t(s) << ")]: " << max
                             << ' ' << chrono <<std::endl;

    std::cout << max << std::endl;

    return 0;
}

int main(int argc, char **argv) {
    if (argc<=1) {
        std::cerr << "usage: " << argv[0] << " #s(stamps) [#](verbosity).\n";
        exit(1);
    }

    const size_t s(atoi(argv[1]));

    if (std::numeric_limits<uint8_t>::max()>s)
        return tmain(argc,argv,uint8_t(s));
    if (std::numeric_limits<uint16_t>::max()>s)
        return tmain(argc,argv,uint16_t(s));
    if (std::numeric_limits<uint32_t>::max()>s)
        return tmain(argc,argv,uint32_t(s));

    return tmain(argc,argv,s);
}

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
// Using _KRange algorithm without the sliding window

template<typename Iterator, typename stype_t>
inline bint _LReach(const Iterator& start, const Iterator& end,
                    const stype_t s, const int verbose) {
    const bint& back(*std::prev(end));				// k>=1
    if (back == __St_One) return s;

    const stype_t spu(s+1);				// s+1 is unreachable
        // Maximal valid index is s*back
        //   thus maximal tested in loop is s*back+1
        //   and maximal starget is at s*back+back=spu*back>=s*back+1
        //   with 0 indexing this gives a table of size: spu*back+1
    std::vector<stype_t> reached(spu*back+1,spu);

    for(auto it=start; it!=end; ++it) reached[*it]=1u;

#if __GSTAMPS_SELMER_LEMMA
    bint vs(back+2);
    const bint& penult(*std::prev(std::prev(end))); // back>1 => k>=2
    std::vector<bint> selmer(s+1);
    std::iota(selmer.begin(), selmer.end(), 1);
    stype_t maxs(0), mins(0);
    for(auto& is: selmer) {
        is *= penult;
        if (is <= back) {
            ++mins;
        }
        is -= back;
    }
#endif

    bint index(1);
    for(; reached[index]<=s; ++index) {
        stype_t& slocal(reached[index]);
        const stype_t vlocal(slocal+1u);
        for(auto right=start; right!=end; ++right) {
            stype_t& starget(reached[index+(*right)]);
            if (starget>vlocal) {
                starget = vlocal;
            }
        }

#if __GSTAMPS_SELMER_LEMMA
            // Selmer's lemma
        maxs = (slocal>maxs?slocal:maxs);
        if ((maxs>mins) && (maxs<s) && (index > selmer[maxs])) {
                // Find maxs-range
            size_t i=index+1;
            for(;(i<vs)&&(reached[i]<=maxs); ++i);
            if (i < vs) {
                // Complete s-range
                if (verbose>0) std::clog << "#[ET(" << (size_t)maxs << '|'
                                         << selmer[maxs] << ")]: " << i
                                         << " -> " << (i-1+(s-maxs)*back)
                                         << std::endl;
                return --i += (s-maxs)*back;
            }
        }
        ++vs;
#endif

    }

    return --index;
}

template<typename List, typename stype_t>
inline bint LReach(const List& points, const stype_t s, const int verbose) {
    StTimer chrono; chrono.start();
    const bint max( _LReach(points.begin(), points.end(), s, verbose) );
    chrono.stop();

    if (verbose>0) std::clog << "#[RRange(" << size_t(s) << ")]: " << max
                             << ' ' << chrono <<std::endl;

    return max;
}


// Reads a basis from std::cin
// Computes the range of that basis with argv[1] stamps
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

    std::cout << LReach(points, s, verbose) << std::endl;

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

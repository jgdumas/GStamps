// ==========================================================================
// GStamps: C++ routines for the Global Postage Stamp Problem
// Authors:
//   L. Colisson, J-G. Dumas, A. Galan, B. Grenet, A. Maignan, D. S. Roche
// ==========================================================================

/**********************************************************************
 * GStamps: Local Postage Stamp Problem
 **********************************************************************/

#include <gstamps.h>

// Loop from 1 to s, with a binary vector
template<typename Iterator, typename stype_t>
inline std::vector<size_t>& Depth_Cover(std::vector<size_t>& numreached,
                                        const Iterator& start,
                                        const Iterator& end,
                                        const stype_t s, const bint maxc,
                                        const int verbose) {
	// Binary cover
    const bint& back(*std::prev(end));						// k>1
    if (back == __St_One) {
        numreached.resize(0); numreached.push_back(1u);
        for(size_t i(1); i<s; i<<=1) {
            numreached.push_back(i);
        }
        if ( (s & (s-1)) != 0u) {
            numreached.back() = s - numreached.back();
        }
        return numreached;
    }

    bint vs(back);
    const bint upper(s*vs+1);
    boost::dynamic_bitset<> reached(upper,false);
    reached[0]=true;										// 0 reached
    for(auto it=start; it!=end; ++it) reached[*it]=true;	// points reached

    numreached.resize(0);
    size_t allreached(end-start);
    size_t prevreached(allreached);
    numreached.push_back(allreached);

    if (verbose>0)
        std::clog << "#[depth 1]: " << numreached.back() << std::endl;

    bint notin(1);
    for(stype_t d(1); d<s; ++d, vs += back) {
        for(size_t i=vs; i>=d; --i) {
            if(reached[i]) {
                for(auto right=start; right!=end; ++right) {
                    const auto& iright(i+(*right));
                    if ( (iright <= maxc) && (! reached[iright]) )
                        ++allreached;
                    reached[i+(*right)]=true;
                }
            }
        }
        if ( (d & (d+1)) == 0u) {
            numreached.push_back(allreached-prevreached);
            prevreached = allreached;
            if (verbose>0)
                std::clog << "#[depth " << (d+1) << "]: " << numreached.back()
                          << ' ' << allreached << std::endl;
            if (((d+1)<<1)>=s) {
                allreached = maxc;
                break;
            }
        }
    }

    numreached.push_back(allreached-prevreached);
    if (verbose>0)
        std::clog << "#[until " << (size_t)s << "]: " << numreached.back()
                  << ' ' << allreached << std::endl;

    return numreached;
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

    rangeprint(std::clog << "# Basis: ", points) << std::endl;

    const bint cmax( Cover(points, s, verbose) );

    std::vector<size_t> numreached;
    Depth_Cover(numreached, points.begin(), points.end(), s, cmax, verbose);

    for(const auto& it: numreached) std::cout << it << ' ';
    std::cout << std::endl;

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

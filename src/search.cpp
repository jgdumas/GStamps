// ==========================================================================
// GStamps: C++ routines for the Global Postage Stamp Problem
// Authors: J-G. Dumas, A. Galan, B. Grenet, A. Maignan, D. S. Roche
// ==========================================================================

/**********************************************************************
 * GStamps: Searches the smallest basis reaching N with s positions,
            Exhaust additional denominations one value at a time
 **********************************************************************/

#include <gstamps.h>

template<typename stype_t>
int tmain(int argc, char **argv, stype_t s) {
    const size_t N(atoi(argv[2]));
    const int verbose(argc>3?atoi(argv[3]):0);
    const int r(argc>4?atoi(argv[4]):1);
    const bool approx(argc>5?atoi(argv[5]):false);

        // Find first k s.t. numbcomb(k+s-1,s)>=N
        // If the converse is true, k can never reach N
    size_t startk(0), n(s);
    for(bint NCb(1); NCb<N; ++startk) {
            // NCb(n+1,s) = NCb(n,s)*(n+1)/(n+1-s) \in\N
        NCb *= ++n;
        NCb /= (n-s);
    }

    if(verbose>0) {
        std::clog << "#[Search] lower bound: " << startk << std::endl;
    }

    std::vector<bint> points;
    bint alpha(0);

        // Search smallest k s.t. FSelect(k,s)>=N
    for(size_t k(startk); alpha < N; ++k) {
        StTimer chrono; chrono.start();
        points.resize(0);
        alpha = FSelect(points, k, s, r, approx, verbose-1);
        chrono.stop();
        if (verbose>0)
            std::clog << "#[Search(" << k << ")] " << (approx?" >=\t":":\t")
                      << alpha << '\t' << chrono << std::endl;
    }

    for(const auto& it: points) std::cout << it << ' ';
    std::cout << std::endl;
    return 0;
}

int main(int argc, char **argv) {
    if (argc<=2) {
        std::cerr << "usage: " << argv[0]
                  << " #s(stamps) #N(range) [#](verbosity) [#r](level) [0/1](approximate).\n";
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

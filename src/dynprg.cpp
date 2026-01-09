// ==========================================================================
// GStamps: C++ routines for the Local & the Global Postage Stamp Problem
// Authors:
//   L. Colisson, J-G. Dumas, A. Galan, B. Grenet, A. Maignan, D. S. Roche
// ==========================================================================

/**********************************************************************
 * GStamps: Dynamic programming search of best Divide & Conquer cut
 **********************************************************************/

#include <gstamps.h>

// For k elements and s positions,
// computes a basis and its associated range
// r governs the number of initial recursive levels searched
template<typename stype_t>
int tmain(int argc, char **argv, stype_t s) {
    const size_t k(atoi(argv[1]));
    const int verbose(argc>3?atoi(argv[3]):0);
    const int r(argc>4?atoi(argv[4]):1);
    const bool approx(argc>5?atoi(argv[5]):false);


    std::vector<bint> points;
    StTimer chrono; chrono.start();
    const bint max = FSelect(points, k, s, r, approx, verbose);
    chrono.stop();
    std::clog << "#[DynProg(" << r <<")] nmax" << (approx?" >= ":": ")
              << max << ' ' << chrono << std::endl;
    for(const auto& it: points) std::cout << it << ' '; std::cout << std::endl;

    return 0;
}

int main(int argc, char **argv) {
    if (argc<=2) {
        std::cerr << "usage: " << argv[0]
                  << " #k(dim.) #s(stamps) [#](verbosity) [#r](level) [0/1](approximate).\n";
        exit(1);
    }

    const size_t s(atoi(argv[2]));

    if (std::numeric_limits<uint8_t>::max()>s)
        return tmain(argc,argv,uint8_t(s));
    if (std::numeric_limits<uint16_t>::max()>s)
        return tmain(argc,argv,uint16_t(s));
    if (std::numeric_limits<uint32_t>::max()>s)
        return tmain(argc,argv,uint32_t(s));

    return tmain(argc,argv,s);
}

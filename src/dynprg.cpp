// ==========================================================================
// GStamps: C++ routines for the Global Postage Stamp Problem
// Authors: J-G. Dumas, A. Galan, B. Grenet, A. Maignan, D. S. Roche
// ==========================================================================

/**********************************************************************
 * GStamps: Dynamic programming search of best Divide & Conquer cut
 **********************************************************************/

#include <gstamps.h>

// For k elements and s positions,
// computes a basis and its associated cover
// r governs the number of initial recursive levels searched
int main(int argc, char **argv) {
    if (argc<=2) {
        std::cerr << "usage: " << argv[0]
                  << " #k(dim.) #s(stamps) [#](verbosity) [#r](level).\n";
        exit(1);
    }

    const size_t k(atoi(argv[1]));
    const size_t s(atoi(argv[2]));
    const int verbose(argc>3?atoi(argv[3]):0);
    const int r(argc>4?atoi(argv[4]):1);


    std::vector<bint> points;
    StTimer chrono; chrono.start();
    const bint max = FSelect(points, k, s, r, verbose);
    chrono.stop();
    std::clog << "#[DynProg(" << r <<")] nmax: " << max
              << ' ' << chrono << std::endl;
    for(const auto& it: points) std::cout << it << ' '; std::cout << std::endl;

    return 0;
}

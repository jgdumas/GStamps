// ==========================================================================
// GStamps: C++ routines for the Global Postage Stamp Problem
// Authors: J-G. Dumas, A. Galan, B. Grenet, A. Maignan, D. S. Roche
// ==========================================================================

/**********************************************************************
 * GStamps: computes the Fibonacci basis and its associated cover
 **********************************************************************/

#include <gstamps.h>

int main(int argc, char **argv) {
    if (argc<=2) {
        std::cerr << "usage: " << argv[0]
                  << " #k(dim.) #s(stamps) [#](verbosity).\n";
        exit(1);
    }

    const size_t k(atoi(argv[1]));
    const size_t s(atoi(argv[2]));
    const int verbose(argc>3?atoi(argv[3]):0);

    std::vector<bint> points;
    const bint max = Fibonacci(points, k);

    if (verbose>0) std::clog << "#[Fibonacci(" << k << ',' << k << ")] nmax: "
                             << max << std::endl;

    std::clog << "#[Fibonacci(" << k << ',' << s << ")] nmax: "
              << Cover(points,s, verbose-1) << std::endl;

    for(const auto& it: points) std::cout << it << ' '; std::cout << std::endl;

    return 0;
}

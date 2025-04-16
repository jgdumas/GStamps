// ==========================================================================
// GStamps: C++ routines for the Global Postage Stamp Problem
// Authors: J-G. Dumas, A. Galan, B. Grenet, A. Maignan, D. S. Roche
// ==========================================================================

/**********************************************************************
 * GStamps: Local Postage Stamp Problem
 **********************************************************************/

#include <gstamps.h>

// Reads a basis from std::cin
// Computes the cover of that basis with argv[1] stamps
int main(int argc, char **argv) {
    if (argc<=1) {
        std::cerr << "usage: " << argv[0] << " #s(stamps) [#](verbosity).\n";
        exit(1);
    }

    const size_t s(atoi(argv[1]));
    const int verbose(argc>2?atoi(argv[2]):0);

    std::vector<bint> points;
    bint tmp;
    std::clog << "# Enter a basis (end by non-digit) ... " << std::endl;
    while(std::cin >> tmp) {
        points.push_back(tmp);
    };

    rangeprint(std::clog << "# Basis: ", points) << std::endl;

    std::cout << Cover(points, s, verbose) << std::endl;

    return 0;
}

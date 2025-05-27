// ==========================================================================
// GStamps: C++ routines for the Global Postage Stamp Problem
// Authors:
//   L. Colisson, J-G. Dumas, A. Galan, B. Grenet, A. Maignan, D. S. Roche
// ==========================================================================

/****************************************************************
 * GStamps Library
 * Alter & Barnett basis
 ****************************************************************/

// Alter, R., & Barnett, J. A. (1977).
// Remarks on the postage stamp problem with applications to computers.
// University of Kentucky. Department of Computer Science.
// https://notatt.com/lsu-stamp.pdf

#include <gstamps.h>

// For k elements and s positions,
// computes Alter & Barnett basis and its associated cover
template<typename stype_t>
int tmain(int argc, char **argv, stype_t s) {
    const size_t k(atoi(argv[1]));
    const int verbose(argc>3?atoi(argv[3]):0);

    std::vector<bint> points;
    const bint max = AlterBernett(points, k,s, verbose);
    std::clog << "#[AlterBernett] nmax: " << max << std::endl;

    for(const auto& it: points) std::cout << it << ' '; std::cout << std::endl;

    return 0;
}

int main(int argc, char **argv) {
    if (argc<=2) {
        std::cerr << "usage: " << argv[0]
                  << " #k(dim.) #s(stamps) [#](verbosity).\n";
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

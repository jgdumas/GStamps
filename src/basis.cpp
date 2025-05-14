// ==========================================================================
// GStamps: C++ routines for the Global Postage Stamp Problem
// Authors: J-G. Dumas, A. Galan, B. Grenet, A. Maignan, D. S. Roche
// ==========================================================================

/****************************************************************
 * GStamps basis: Generic algorithm calling all methods,
 *                with memoization of solutions
 ****************************************************************/

#include <gstamps.h>

// For k elements and s positions,
// computes a basis and its associated cover
template<typename stype_t>
int tmain(int argc, char **argv, stype_t s) {
    const size_t k(atoi(argv[1]));
    const int verbose(argc>3?atoi(argv[3]):0);
    const bool approx(argc>4?atoi(argv[4]):true);

    std::vector<bint> points;
    const bint max = FSelect(points, k, s, 0, approx, verbose);
    std::clog << "#[FSelect] nmax" << (approx?">= ":": ") << max << std::endl;

    for(const auto& it: points) std::cout << it << ' '; std::cout << std::endl;

    return 0;
}

int main(int argc, char **argv) {
    if (argc<=2) {
        std::cerr << "usage: " << argv[0]
                  << " #k(dim.) #s(stamps) [#](verbosity) [0/1](approximate).\n";
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

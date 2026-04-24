// ==========================================================================
// GStamps: C++ routines for the Local & the Global Postage Stamp Problem
// Authors:
//   L. Colisson, J-G. Dumas, A. Galan, B. Grenet, A. Maignan, D. S. Roche
// ==========================================================================

/**********************************************************************
 * GStamps: Local Postage Stamp Problem
 **********************************************************************/
#include <gstamps.h>

// ==========================================================================
// Computes the range & decompositions of a basis with a stamps
template<typename List, typename stype_t>
inline bint LDecompose(const List& points, const stype_t s, const int verbose) {
    StTimer chrono; chrono.start();
    const auto sums( _Decompose(points, s, verbose) );
    chrono.stop();

    const size_t range(sums.size()-1);

    if (verbose>0) std::clog << "#[RRange(" << size_t(s) << ")]: "
                             << range << ' ' << chrono <<std::endl;

    std::clog << std::string(30,'#') << std::endl;
    for(size_t i=1; i<sums.size(); ++i) {
        std::cout << "# " << i << '=';
        const auto& Ui(sums[i].second);
        const size_t endj(Ui.size()-1);
        for(size_t j=0; j<Ui.size(); ++j) {
            std::cout << points[Ui[j]];
            if (j<endj) std::cout << '+';
        }
        std::cout << std::endl;
    }

    return range;
}
// ==========================================================================



// ==========================================================================
// Reads a basis from std::cin
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

    LDecompose(points, s, verbose);

    return 0;
}
// ==========================================================================



// ==========================================================================
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
// ==========================================================================

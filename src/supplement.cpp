// ==========================================================================
// GStamps: C++ routines for the Global Postage Stamp Problem
// Authors: J-G. Dumas, A. Galan, B. Grenet, A. Maignan, D. S. Roche
// ==========================================================================

/****************************************************************
 * GStamps: Exhaust additional denominations one value at a time
 ****************************************************************/

#include <gstamps.h>
template<typename stype_t>
int tmain(int argc, char **argv, stype_t s) {
    const size_t k(atoi(argv[1]));
    const int verbose(argc>3?atoi(argv[3]):0);

    std::vector<bint> prescb;
    int g(argc>4?atoi(argv[4]):1); if (g<=0) g=2;

    if (argc <= 5) {
        bint tmp;
        std::clog << "# Enter prescribed denominations (end by non-digit) ... "
                  << std::endl;
        while(std::cin >> tmp) {
            prescb.push_back(tmp);
        };
    } else {
        const int r(argc>5?atoi(argv[5]):1);
        const size_t redk(k>g?k-g:1);
        FSelect(prescb, redk, s, r, false, verbose);
    }

    if (prescb.size()>k) {
		prescb.resize(k);
	}
    if (prescb.size()<1) {
        prescb.push_back(1);
    }

    StTimer chrono; chrono.start();
    bint max(_KCover(prescb.begin(),prescb.end(),s));

    std::vector<size_t> suppvec;
    for(int j = k; j > prescb.size(); j -= g) {
        suppvec.push_back(j);
    }

    for (auto it = suppvec.rbegin(); it != suppvec.rend(); ++it) {
        max = par_complement(prescb, *it, s, verbose-1);
        if (verbose>0)
            std::clog << "#[Spmt(" << *it << ")] max: " << max << std::endl;
    }
    chrono.stop();

    std::clog << "#[Spmt] nmax: " << max
              << ' ' << chrono << std::endl;
    for(const auto& it: prescb) std::cout << it << ' '; std::cout << std::endl;

    return 0;
}

int main(int argc, char **argv) {
    if (argc<=2) {
        std::cerr << "usage: " << argv[0]
                  << " #k(dim.) #s(stamps) [#](verbosity)"
                  << " [#g](gap) [#r](level).\n";
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

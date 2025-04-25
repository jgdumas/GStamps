// ==========================================================================
// GStamps: C++ routines for the Global Postage Stamp Problem
// Authors: J-G. Dumas, A. Galan, B. Grenet, A. Maignan, D. S. Roche
// ==========================================================================

/****************************************************************
 * GStamps: Exhaust all additional denominations in parallel
 ****************************************************************/

#include <omp.h>
#include <gstamps.h>

int main(int argc, char **argv) {
    if (argc<=2) {
        std::cerr << "usage: " << argv[0]
                  << " #k(dim.) #s(stamps) [#](verbosity)"
                  << " [#g](gap) [#r](level).\n";
        exit(1);
    }

    const size_t k(atoi(argv[1]));
    const size_t s(atoi(argv[2]));
    const int verbose(argc>3?atoi(argv[3]):0);

    std::vector<bint> prescb;
    if (argc <= 4) {
        bint tmp;
        std::clog << "# Enter prescribed denominations (end by non-digit) ... "
                  << std::endl;
        while(std::cin >> tmp) {
            prescb.push_back(tmp);
        };

    } else {
        int g(atoi(argv[4])); if (g<=0) g=3;
        const int r(argc>5?atoi(argv[5]):1);
        const size_t redk(k>g?k-g:1);
        FSelect(prescb, redk, s, r, false, verbose);
    }

    if (prescb.size()<1) {
        prescb.push_back(1);
    }



    int64_t numthreads(1);
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp single
        {
            numthreads = omp_get_num_threads();
        }
#endif


    bint max;

    StTimer chrono; chrono.start();
    if (prescb.size()>=k) {
        max = Cover(prescb, s, verbose);
	} else {
        max = par_complement(prescb, k, s, verbose);
    }
    chrono.stop();


    std::clog << "#[Cpmt " << numthreads << "t] nmax: " << max
              << ' ' << chrono << std::endl;
    for(const auto& it: prescb) std::cout << it << ' '; std::cout << std::endl;

    return 0;
}

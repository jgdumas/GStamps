// ==========================================================================
// GStamps: C++ routines for the Global Postage Stamp Problem
// Authors: J-G. Dumas, A. Galan, B. Grenet, A. Maignan, D. S. Roche
// ==========================================================================

/**********************************************************************
 * GStamps: basis in geometric progression
 **********************************************************************/

#include <gstamps.h>

template<typename stype_t>
int tmain(int argc, char **argv, stype_t s) {
    const size_t k(atoi(argv[1]));
    const int verbose(argc>3?atoi(argv[3]):0);
    assert(k>0); assert(s>0);

    const bint q1( (k>(s+1)?2: (k+s+1)/k) ), q2(q1+1);

    bint t1(q1), t2(q2);
    std::vector<bint> p1{1,q1}, p2{1,q2};
    p1.resize(k); p2.resize(k);

    for(size_t i(2); i<k; ++i) {
        p1[i] = (t1 *= q1) ;
        p2[i] = (t2 *= q2);
    }
    const bint i1(s/(q1-1)), i2(s/(q2-1));

    if (i1 < (k-1)) {
        t1 = p1[i1]; t1 *= (s-i1*(q1-1)+2); t1 -= 2;
    } else {
        t1 = p1[k-1]; t1 *= (s-(k-1)*(q1-1)+2); t1 -= 2;
    }
    if (i2 < (k-1)) {
        t2 = p2[i2]; t2 *= (s-i2*(q2-1)+2); t2 -= 2;
    } else {
        t2 = p2[k-1]; t2 *= (s-(k-1)*(q2-1)+2); t2 -= 2;
    }

    if (t2 > t1) {
        std::clog << "#[Geom(" << q1 << ")] less: " << t1 << std::endl;
        std::clog << "#[Geom(" << q2 << ")] nmax: " << t2 << std::endl;
        for(const auto& it: p2) std::cout << it << ' '; std::cout << std::endl;
    } else {
        std::clog << "#[Geom(" << q2 << ")] less: " << t2 << std::endl;
        std::clog << "#[Geom(" << q1 << ")] nmax: " << t1 << std::endl;
        for(const auto& it: p1) std::cout << it << ' '; std::cout << std::endl;
    }

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

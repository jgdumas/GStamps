// ==========================================================================
// GStamps: C++ routines for the Local & the Global Postage Stamp Problem
// Authors:
//   L. Colisson, J-G. Dumas, A. Galan, B. Grenet, A. Maignan, D. S. Roche
// ==========================================================================

/**********************************************************************
 * GStamps: Local Postage Stamp Problem
 **********************************************************************/
#include <gstamps.h>

template<typename stype_t> using PsVs=std::pair<stype_t, std::vector<stype_t> >;

// ==========================================================================
// Reads a basis from std::cin
// Computes the range and decompositions that basis with argv[1] stamps
// Modification of the _KRange algorithm without the sliding window
template<typename List, typename stype_t>
std::vector<PsVs<stype_t>> _Decompose(const List& points, const stype_t s,
                                      const int verbose) {
    const auto& back(points.back());				// k>=1
    const stype_t spu(s+1);				// s+1 is unreachable

    if (back == __St_One) {
        std::vector<PsVs<stype_t>> reached(spu);
        for(size_t i(1); i<spu; ++i) {
            reached[i]={0,std::vector<stype_t>(i,0)};
        }
        return reached;
    }

    const PsVs<stype_t> spustart{spu, {0u}};
        // Maximal valid index is s*back
        //   thus maximal tested in loop is s*back+1
        //   and maximal starget is at s*back+back=spu*back>=s*back+1
        //   with 0 indexing this gives a table of size: spu*back+1
    std::vector<PsVs<stype_t>> reached(spu*back+1,spustart);

    for(stype_t i=0; i<points.size(); ++i)
        reached[points[i]]= PsVs<stype_t>{1u,{i}};

    size_t index(1);
    for(; reached[index].first<=s; ++index) {
        const auto& slocal(reached[index]);
        const stype_t slfirst(slocal.first);
        const stype_t vlocal(slfirst+1u);
        for(auto right=slocal.second.back(); right<points.size(); ++right) {
            auto& starget(reached[index+points[right]]);
            if (starget.first>vlocal) {
                starget.first = vlocal;
                starget.second.resize(0);
                starget.second.assign(slocal.second.begin(),
                                      slocal.second.end());
                starget.second.push_back(right);
            }
        }

    }

    reached.resize(index); // range is --index

    return reached;
}
// ==========================================================================




// ==========================================================================
// Computes the range & decompositions of a basis with a stamps
template<typename List, typename stype_t>
inline bint LDecompose(const List& points, const stype_t s, const int verbose) {
    StTimer chrono; chrono.start();
    const auto sums( _Decompose(points, s, verbose) );
    chrono.stop();
    
    std::clog << "# sums:\t";
    for(const auto& it: sums) {
        std::clog << it.first << ':' << '[';
        for(const auto& e: it.second) {
            std::clog << (size_t)e << ' ';
        }
        std::clog << ']' << '\t';
    }
    std::clog << std::endl;


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

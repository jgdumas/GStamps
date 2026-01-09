// ==========================================================================
// GStamps: C++ routines for the Local & the Global Postage Stamp Problem
// Authors:
//   L. Colisson, J-G. Dumas, A. Galan, B. Grenet, A. Maignan, D. S. Roche
// ==========================================================================

/****************************************************************
 * GStamps Library, inline implementations
 ****************************************************************/


// ============================================
// Special cases

inline bint Fibonacci(std::vector<bint>& points, const size_t k,
	       const int verbose) {
    points.resize(0); points.reserve(k);
    bint f1(1), f2(1); points.push_back(f2);
    for(size_t i=1; i<k; ++i) {
	f1 += f2;
	f2 += f1;
	points.push_back(f2);
    }
    f1 += f2; --f1;
    if (verbose>0) rangeprint(std::clog << "#[Fibonacci] max: " << f1
			      << ", points: ", points) << std::endl;
    return f1;
}


template<typename List>
inline bint KloveMossige(List& points, const size_t k, const int verbose) {
    const size_t x( (k+6)/7 ); const size_t xmu(x-1);
    const size_t y(k-(x<<2)+3);
    const bint xy(bint(x)*bint(y));

    if (verbose>0) std::clog << "#[KM] Precomp. : " << k
			     << ", x: " << x << ", y: " << y << std::endl;

    points.push_back(__St_One);

    for(size_t i=2; i<x; ++i) points.emplace_back(i);

    for(size_t i=(x>1?1:2); i<y; ++i) points.push_back(bint(x)*i+xmu);

    for(size_t i=0; i<xmu; ++i) points.push_back(xy+i);

    bint u = Range(points,2,verbose);
    for(size_t i=0; i<x; ++i) points.push_back(u+xmu*i);

    bint v = Range(points,2,verbose);
    for(size_t i=0; i<x; ++i) points.push_back(v+xmu*i);

    bint c(k<<1); c-=(7*x); c<<=1; c+=12; c*=x; c-=6;
    return c;
}



// AlterBernett sub-range computation sub-routine
template<typename stype_t>
inline std::vector<bint>& SubRange(std::vector<bint>& B,
                                   const size_t& q, const stype_t& s,
                                   const size_t& r) {
    assert(B.size()>=2);
    auto B1 = B.back();
    auto B0 = B.end()[-2];
    for(stype_t i=2;i<=s; ++i) {
        B.push_back((q+2)*B1-B0+q);
        B1 = B.back();
        B0 = B.end()[-2];
    }
    if (r != 0) B.push_back(B1+r*(B0+1));
    return B;
}


// log[2](v)
uint64_t MSB(uint64_t v) { uint64_t r(0); for( ; v>>=1; ++r) {}; return r; }

template<typename stype_t>
inline bint AlterBernett(std::vector<bint>& points,
			 const size_t k, const stype_t s, const int verbose) {
    points.resize(0); points.reserve(k);
    if (k<=s) {
        const bint m1 ( Fibonacci(points,k,verbose-1) );
        bint max(s-k); max *= points.back(); max += m1;
        if (verbose>0) rangeprint(std::clog << "#[AB] max: " << max
                                  << ", points: ", points) << std::endl;
        return max;
    }

    const size_t r(k%s), q( (k-r)/s );

    if (verbose>0) std::clog << "#[AlterBernett] Precomp. : " << k
                             << ", Prof. : " << MSB(s)
                             << ", Stamps: " << size_t(s)
                             << ", q: " << q << ", r: " << r << std::endl;

    std::vector<bint> B(2); B[0]=__St_Zero; B[1]=q;
    SubRange(B,q,s,r) ;

    if (verbose>0) {
        std::clog << "#[AlterBernett] B: ";
        for(const auto& it:B) std::clog << it << ' ';
        std::clog << std::endl;
    }

    for(size_t j=1; j<=q; ++j) points.push_back(j);

    bint Fi(1), Di(1);
    for(size_t i=2; i<=s; ++i) {
        Fi = (B[i-1]<<1)-B[i-2]+1;
        Di = B[i-1]+1;
        for(size_t j=0; j<q; ++j) points.push_back(Fi+Di*j);
        if (verbose>0)
            std::clog << "#[AB] F[" << i << "]: " << Fi << '\t'
                      << "D[" << i << "]: " << Di << '\t'
                      << "B[-1]: " << B[i-1] << '\t'
                      << "B[-2]: " << B[i-2] << '\t'
                      << std::endl;
    }

    const bint Fs = Fi+Di*q;
    if (verbose>0)
        std::clog << "#[AB] F[" << size_t(s) << "]: " << Fs  << '\t'
                  << "D[" << size_t(s) << "]: " << Di << std::endl;

    for(size_t j=0; j<r; ++j) points.push_back(Fs+Di*j);

    return B.back();
}

template<typename stype_t>
inline bint BalGreedy(std::vector<bint>& points,
                      const size_t k, const stype_t s, const int verbose) {
    assert(k>=s);
    points.resize(0); points.reserve(k);
    const stype_t s1(k%s), s0(s-s1);
    const size_t q0( (k-s1)/s ), q1(q0+1);

    if (s1 == 0) return AlterBernett(points,k,s,verbose-1);

    if (verbose>0) std::clog << "#[BalGreedys] Precomp. : " << k
                             << ", Prof. : " << MSB(s)
                             << ", Stamps: " << size_t(s)
                             << ", q0: " << q0 << ", s0: " << s0
                             << ", q1: " << q1 << ", s1: " << s1
                             << std::endl;

    std::vector<bint> B(2); B[0]=__St_Zero; B[1]=q0;
    SubRange(B,q0,s0,0) ;
//     std::clog << "B: "; for(const auto& it:B) std::clog << it << ' ';
//     std::clog << std::endl;
    for(size_t j=1; j<=q0; ++j) points.push_back(j);

    bint Fi(1), Di(1);
    for(stype_t i=2; i<=s0; ++i) {
        Fi = (B[i-1]<<1)-B[i-2]+1;
        Di = B[i-1]+1;
        for(size_t j=0; j<q0; ++j) points.push_back(Fi+Di*j);
// std::clog << "F[" << i << "]: " << Fi << '\t'
//           << "D[" << i << "]: " << Di << '\t'
//           << ", B[-1]: " << B[i-1] << '\t'
//           << ", B[-2]: " << B[i-2] << '\t'
//           << std::endl;
    }

//     for(const auto& it: points) std::clog <<it<< ' '; std::clog << std::endl;
//     std::clog << "Range B: "; Range(points,s);

    SubRange(B,q1,s1+1,0);
//     std::clog << "B: "; for(const auto& it:B) std::clog << it << ' ';
//     std::clog << std::endl;

    for(stype_t i=0; i<s1; ++i) {
        Fi = (B[s0+i]<<1)-B[s0+i-1]+1;
        Di = B[s0+i]+1;
        for(size_t j=0; j<q1; ++j) points.push_back(Fi+Di*j);
// std::clog << "F[" << i << "]: " << Fi << '\t'
//           << "D[" << i << "]: " << Di  << '\t'
//           << ", B[" << (s0+i) << "]: " << B[s0+i] << '\t'
//           << ", B[" << (s0+i-1) << "]: " << B[s0+i-1] << '\t'
//           << std::endl;
    }

//     const bint Fs = Fi+Di*q0;
// // std::clog << "L[" << size_t(s) << "]: " << Fs << std::endl;
// // std::clog << "D[" << size_t(s) << "]: " << Di << std::endl;
//     for(size_t j=0; j<r; ++j) points.push_back(Fs+Di*j);

    return B.back();
}

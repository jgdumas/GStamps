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
// ============================================

// ============================================
// Fibonacci sequence
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
    if (verbose>0) ScopePrint(std::clog << "#[Fibonacci] max: " << f1
			      << ", points: ", points) << std::endl;
    return f1;
}
// ============================================


// ============================================
// Klove & Mossige, s=2
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
// ============================================



// ============================================
// Alter & Barnett

// AlterBarnett sub-range computation sub-routine
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

// Alter & Barnett core algorithm
template<typename stype_t>
inline bint AlterBarnett(std::vector<bint>& points,
			 const size_t k, const stype_t s, const int verbose) {
    points.resize(0); points.reserve(k);
    if (k<=s) {
        const bint m1 ( Fibonacci(points,k,verbose-1) );
        bint max(s-k); max *= points.back(); max += m1;
        if (verbose>0) ScopePrint(std::clog << "#[AB] max: " << max
                                  << ", points: ", points) << std::endl;
        return max;
    }

    const size_t r(k%s), q( (k-r)/s );

    if (verbose>0) std::clog << "#[AlterBarnett] Precomp. : " << k
                             << ", Prof. : " << MSB(s)
                             << ", Stamps: " << size_t(s)
                             << ", q: " << q << ", r: " << r << std::endl;

    std::vector<bint> B(2); B[0]=__St_Zero; B[1]=q;
    SubRange(B,q,s,r) ;

    if (verbose>0) {
        std::clog << "#[AlterBarnett] B: ";
        for(const auto& it:B) std::clog << it << ' ';
        std::clog << std::endl;
    }

    for(size_t j=1; j<=q; ++j) points.push_back(j);

    bint Fi(1), Di(1);
    for(size_t i=2; i<=s; ++i) {
        Fi = (B[i-1]<<1)-B[i-2]+1;
        Di = B[i-1]+1;
        for(size_t j=0; j<q; ++j) points.push_back(Fi+Di*j);
        if (verbose>1)
            std::clog << "#[AB] F[" << i << "]: " << Fi << '\t'
                      << "D[" << i << "]: " << Di << '\t'
                      << "B[-1]: " << B[i-1] << '\t'
                      << "B[-2]: " << B[i-2] << '\t'
                      << std::endl;
    }

    const bint Fs = Fi+Di*q;
    if (verbose>1)
        std::clog << "#[AB] F[" << size_t(s) << "]: " << Fs  << '\t'
                  << "D[" << size_t(s) << "]: " << Di << std::endl;

    for(size_t j=0; j<r; ++j) points.push_back(Fs+Di*j);

    return B.back();
}
// ============================================


// ============================================
// Improvement over Alter & Barnett, with balanced basis
template<typename stype_t>
inline bint Balanced(std::vector<bint>& points,
                      const size_t k, const stype_t s, const int verbose) {
    assert(k>=s);
    points.resize(0); points.reserve(k);
    if (k<=s) {
        const bint m1 ( Fibonacci(points,k,verbose-1) );
        bint max(s-k); max *= points.back(); max += m1;
        if (verbose>0) ScopePrint(std::clog << "#[BG] max: " << max
                                  << ", points: ", points) << std::endl;
        return max;
    }
    const stype_t s0(k%s), s1(s-s0);
    const size_t q1( (k-s0)/s ), q0(q1+1);

    if (s0 == 0) return AlterBarnett(points,k,s,verbose-1);

    if (verbose>0) std::clog << "#[Balanced] Precomp. : " << k
                             << ", Prof. : " << MSB(s)
                             << ", Stamps: " << size_t(s)
                             << ", q0: " << q0 << ", s0: " << (size_t)s0
                             << ", q1: " << q1 << ", s1: " << (size_t)s1
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
        if (verbose>1)
            std::clog << "#[BG " << q0 << "] F[" << i << "]: " << Fi << '\t'
                      << ", D[" << i << "]: " << Di << '\t'
                      << ", B[-1]: " << B[i-1] << '\t'
                      << ", B[-2]: " << B[i-2] << std::endl;
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
        if (verbose>1)
            std::clog << "#[BG " << q1 << "] F[" << i << "]: " << Fi << '\t'
                      << ", D[" << i << "]: " << Di  << '\t'
                      << ", B[" << (size_t)(s0+i) << "]: " << B[s0+i] << '\t'
                      << ", B[" << (size_t)(s0+i-1) << "]: " << B[s0+i-1]
                      << std::endl;
    }


    return B.back();
}
// ============================================



// ============================================
// J. Robinson, Some Postage Stamp 2-Bases, J. of integer sequences (12), 09.1.1
#define __GSTAMPS_PAs_MAX_ 11
inline bint PreambleAmble(std::vector<bint>& points,
                          const size_t k, const int verbose) {
    const static std::vector<bint> PAs[__GSTAMPS_PAs_MAX_]={
{1,3,4,7,8,9,16,17,21,24,35},
{1,2,5,7,10,11,19,21,22,25,29,30,43},
{1,2,5,6,8,9,13,19,22,27,29,33,40,41,56},
{1,2,3,7,8,9,12,15,22,26,30,36,37,43,45,61},
{1,2,5,6,7,12,13,16,26,28,31,37,38,42,44,49,66},
{1,2,3,6,9,11,12,15,16,27,32,37,45,48,52,55,61,62,80},
{1,2,4,5,11,13,14,19,29,35,37,43,46,47,50,52,56,58,68,88},
{1,2,3,6,10,14,17,19,26,29,36,41,49,51,54,55,58,60,67,74,95},
{1,3,5,7,8,12,14,18,26,32,33,42,43,50,60,63,68,79,81,83,97,105},
{1,3,5,6,13,15,16,18,22,38,41,44,47,52,55,58,59,60,74,80,81,91,93,117},
{1,3,4,6,7,14,16,19,20,28,36,38,39,48,49,60,61,70,76,77,89,93,95,99,109,135}
};
    const size_t ko2(k>>1);
    bint max(0);
    for(size_t i=0; i<__GSTAMPS_PAs_MAX_; ++i) {
        std::vector<bint> cand;
        cand.assign(PAs[i].begin(),PAs[i].end());
        const size_t preas(cand.size());
        if (preas>ko2) break;
        const size_t dblps(preas<<1);
        const size_t ambls(k-dblps);
        for(size_t j=0; j<ambls; ++j) cand.push_back(cand.back()+preas);
        const bint ak(cand[ko2-1]+cand[k-ko2-1]);
        const bint nk(ak<<1);
        if (verbose>0) std::clog << "#[PA] pa[" << i << "], n: " << nk << std::endl;
        if (nk>max) {
            for(size_t j=cand.size()+1;j<k;++j) {
                cand.push_back(ak-cand[k-j-1]);	// Symmetric base
            }
            cand.push_back(ak);
            points.assign(cand.begin(), cand.end());
            max=nk;
        }
    }
    return max;
}
// ============================================

// ==========================================================================
// GStamps: C++ routines for the Global Postage Stamp Problem
// Authors: J-G. Dumas, A. Galan, B. Grenet, A. Maignan, D. S. Roche
// ==========================================================================

/****************************************************************
 * GStamps Library, inline implementations
 ****************************************************************/

#include <gstamps.h>



// ============================================
// Masking Tools
// upmask: Round up to the next (highest power of 2, minus 1) of (input+1)

// log[2](v)
uint64_t MSB(uint64_t v) { uint64_t r(0); for( ; v>>=1; ++r) {}; return r; }


#ifdef __GSTAMPS_EXTENDED_PRECISION
Givaro::Integer upmask(const Givaro::Integer& w) {
    Givaro::Integer v(w);
    uint32_t exp(1);
    for(Givaro::Integer shi(1); shi>0; exp <<=1) {
        shi = v >> exp;
        v |= shi;
    }
    return std::move(v);
}

#else
// See: https://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
uint64_t upmask(const uint64_t& w) {
    bint v(w);
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v |= v >> 32;
    return std::move(v);
}
#endif



// ============================================
// Printing Tools

template<typename Cont>
inline std::ostream& rangeprint(std::ostream& out, const Cont& v) {
    typename Cont::value_type prev(0), curr(0);
    bool ft(true); size_t pt(0);
    for(const auto& it: v) {
        if (ft) {
            prev = it;
            out << prev;
            ft = false;
        } else {
            curr = it;
            if (prev+1 == curr) {
                ++pt;
                if (pt == 2) {
                    out << "..";
                }
            } else {
                if (pt >= 1) {
                    if (pt == 1) out << ' ';
                    pt = 0;
                    out << prev;
                }
                out << ' ' << curr;
            }
            prev = curr;
        }
    }
    if (pt >= 1) {
        if (pt == 1) out << ' ';
        return out << curr;
    }
    return out;
}


template<typename Cont>
inline typename Cont::value_type firstrange(const Cont& v) {
    typename Cont::value_type max(0);
    for(const auto& it: v)
	if (++max != it) break;
    if (v.count(max) == 0) --max;
    return max;
}

template<typename Cont>
inline std::ostream& firstrangeprint(std::ostream& out,
			      const Cont& v, const bint& max) {
    if (max > 1) out << "1..";
    out << max;
    if (v.size() > max) out << " ...";
    return out;
}


// ============================================
// Cover: local postage start problem

template<typename Iterator, typename stype_t>
inline bint _SCover(const Iterator& start, const Iterator& end, const stype_t s) {
	// Binary cover
    const bint& back(*std::prev(end));				// k>=1
    if (back == __St_One) return s;
    const bint& penult(*std::prev(std::prev(end))); // back>1 => k>=2
    bint vs(back), lb(penult);
    const bint upper(s*vs+1);
    boost::dynamic_bitset<> reached(upper,false);
    reached[0]=true;                                // 0 reached
    for(auto it=start; it!=end; ++it) reached[*it]=true;   // points reached

    for(stype_t d(1); d<s; ++d) {
        bint notin(vs+1);
        for(size_t i=vs; i>0; --i) {
            if(reached[i]) {
                for(auto right=start; right!=end; ++right) {
                    reached[i+(*right)]=true;
                }
            } else {
                notin=i;
            }
        }

        vs += back;

            // Test Selmer's Lemma for early termination
        lb += penult;
        if (notin>back) {
                // c is larger than ak, try early termination
            if (notin >= (lb-back) ) {
                    // Cover will now surely attain c+(s-d)ak
//                 std::clog << "#[_C(" << size_t(s) << ")|" << d << "] "
//                           << (notin-1+(s-d)*back) << std::endl;
                return --notin += (s-d)*back;
            }
        }
    }

    size_t max(reached.size()-1);
    for(size_t jr(1); jr<reached.size(); ++jr){
        if (! reached[jr]) {
            max = jr-1;
            break;
        }
    }

    return max;
}


template<typename Iterator, typename stype_t>
inline bint _Cover(const Iterator& start, const Iterator& end, const stype_t s) {
    const bint& back(*std::prev(end));	// k>=1
    if (back == __St_One) return s;

    const bint window(upmask(back));	// highest 1-full mask gt
    std::vector<stype_t> reached(window+1u, 0u);

    for(auto it=start; it!=end; ++it) reached[*it]=1u;

    bint index(1);
    for(; reached[index & window]<=s; ++index) {
        const stype_t slocal(reached[index & window]+1u);
        for(auto it=start; it!=end; ++it) {
            stype_t& starget(reached[ (index+(*it)) & window]);
            const stype_t vtarget(starget);
            if ( (vtarget == 0u) || (vtarget>slocal) ) {
                starget = slocal;
            }
        }
        reached[index & window]=0u;		// clean up sliding window
    }

    return --index;

}

template<typename List, typename stype_t>
inline bint Cover(const List& points, const stype_t s, const int verbose) {
    if (verbose>1)
        rangeprint(std::clog << "#[Cover] Basis: ", points) << std::endl;

    StTimer chrono; chrono.start();
    const bint max( _Cover(points.begin(), points.end(), s) );
    chrono.stop();

    if (verbose>1) {
        std::clog << "#[Cover(" << size_t(s) << ")]: 1.." << max
                  << " ..." << std::endl;
    }

    if (verbose>0) std::clog << "#[Cover(" << size_t(s) << ")]: " << max
                             << ' ' << chrono <<std::endl;
    return max;
}

// ============================================
// Exhaustive search

template<typename List, typename stype_t>
inline bint FixedPoints(List& pmax,
			const List& points, const stype_t s, const size_t i) {
    pmax.resize(i); pmax.reserve(points.size());
    for(auto it=points.begin()+i; it!=points.end(); ++it) {
	pmax.push_back( _Cover(points.begin(), std::next(it), s) );
    }
    return pmax.back();
}



template<typename List, typename stype_t>
inline bint BruteForce(List& points, size_t k, stype_t s, const int verbose) {
    assert(k>0);
    if (k == 1) {
        points = {1};
        return s;
    }
    points.resize(k); std::iota(points.begin(), points.end(), 1);
    List pointsmax(points);

    List covsmax;
    bint max( FixedPoints(covsmax, points, s) );

        // Next set of points
    while (true){
            // Find where next increment is needed
        size_t i=points.size()-1;
        while ((i > 0) && (points[i] == (covsmax[i-1]+1))) { --i; }

            // Nothing left
        if (i == 0) {
            points.assign(pointsmax.begin(), pointsmax.end());
            return max;
        }

            // Increment is needed points[i] and reset all higher denominations
        std::iota(points.begin()+i,points.end(),points[i]+1);

            // Compute new maxima for the denominations after i
        const bint max2 = FixedPoints(covsmax, points, s, i-1);

            // Save new points if better
        if (max2 > max) {
            if (verbose>0) {
                std::clog << "#[Brute] " << max2 << ": ";
                for(const auto& it: points) std::clog << it << ' ';
                std::clog << std::endl;
            }

            max = max2;
            pointsmax.assign(points.begin(), points.end());
        }
    }
}


// Stores only reached values, not whole array
template<typename List>
inline bint Reach(List& points, const size_t h, const int verbose) {
    StTimer chrono; chrono.start();

    std::set<bint> reached;
    for(const auto& it: points)
	reached.insert(it);
    if (verbose>2)
	rangeprint(std::clog << "#[Reach] Basis: ", reached) << std::endl;

    for(size_t d(0); d<h; ++d) {
	std::vector<bint> v(reached.begin(), reached.end());
	if (verbose>2) rangeprint(std::clog << "#[Reach] h=" << d
				  << ": ", reached) << std::endl;
	else if (verbose>1)
	    firstrangeprint(std::clog << "#[Reach] h=" << d << ": ",
			    reached, firstrange(reached)) << std::endl;
	for(const auto& left: v) {
	    for(const auto& right: v) {
		reached.insert(left+right);
	    }
	}
    }
    chrono.stop();

    const bint max(firstrange(reached));
    if (verbose>1) firstrangeprint(std::clog << "#[Reach] h=" << h
				   << ": ", reached, max) << std::endl;
    if (verbose>0) std::clog << "#[Reach<" << h << ">]: " << chrono <<std::endl;
    return max;
}

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

    bint u = Cover(points,2,verbose);
    for(size_t i=0; i<x; ++i) points.push_back(u+xmu*i);

    bint v = Cover(points,2,verbose);
    for(size_t i=0; i<x; ++i) points.push_back(v+xmu*i);

    bint c(k<<1); c-=(7*x); c<<=1; c+=12; c*=x; c-=6;
    return c;
}



// AlterBernett range computation sub-routine
template<typename stype_t>
inline std::vector<bint>& Range(std::vector<bint>& B,
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
    Range(B,q,s,r) ;

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
    Range(B,q0,s0,0) ;
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
//     std::clog << "Cover B: "; Cover(points,s);

    Range(B,q1,s1+1,0);
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


// ============================================
// Known extremal cases

inline bint kThree(std::vector<bint>& points, const size_t s,
		   const int verbose) {
    assert(s>0);
    points.resize(0); points.reserve(3);
    bint max(0);
    if (s<23) {
	const static size_t kThreeLow[22][4]={
	    {3,1,2,3},
	    {8,1,3,4},
	    {15,1,4,5},
	    {26,1,5,8},
	    {35,1,6,7},
	    {52,1,7,12},
	    {69,1,8,13},
	    {89,1,9,14},
	    {112,1,9,20},
	    {146,1,10,26},
	    {172,1,10,26},
	    {212,1,11,37},
	    {259,1,13,34},
	    {302,1,12,52},
	    {354,1,12,52},
	    {418,1,15,54},
	    {476,1,14,61},
	    {548,1,15,80},
	    {633,1,18,65},
	    {714,1,17,91},
	    {805,1,17,91},
	    {902,1,20,92}};
	points.push_back(kThreeLow[s-1][1]);
	points.push_back(kThreeLow[s-1][2]);
	points.push_back(kThreeLow[s-1][3]);
	max=kThreeLow[s-1][0];
    } else {
	const static size_t kThreeTab[9][6]={
	    {3,1,1,0,0,0},
	    {3,1,1,0,0,1},
	    {5,2,1,1,0,1},
	    {5,2,1,1,0,2},
	    {7,3,1,2,0,2},
	    {6,2,2,2,1,2},
	    {8,3,2,3,1,2},
	    {8,3,2,3,1,3},
	    {10,4,2,4,1,3}};

	const size_t t(s/9);
	const size_t r(s-9*t);
	bint a2(t); a2 *= 6; a2+=kThreeTab[r][0];
	bint a3(t); a3 <<=1; a3+=kThreeTab[r][1];
	bint a32(t); a32 <<=1; a32+=kThreeTab[r][2];
	a32 *= a2; a3 += a32;
	max = t; max <<=2; max+=kThreeTab[r][3];
	bint a42(t); a42 <<=1; a42+=kThreeTab[r][4];
	a42 *= a2; max += a42;
	bint a43(t); a43 *= 3; a43+=kThreeTab[r][5];
	a43 *= a3; max += a43;
	points.push_back(__St_One);
	points.push_back(std::move(a2));
	points.push_back(std::move(a3));
    }
    return max;
}

#define __GSTAMPS_kFour_sKNOWN_ 54
inline bint kFour(std::vector<bint>& points, const size_t s,
                  const bool approx, const int verbose) {
    assert(s>0);
    points.resize(0); points.reserve(4);
    const static size_t kFourLow[54][5]={
	{4,1,2,3,4},
	{12,1,3,5,6},
	{24,1,4,7,8},
	{44,1,3,11,18},
	{71,1,4,12,21},
	{114,1,4,19,33},
	{165,1,5,24,37},
	{234,1,6,25,65},
	{326,1,5,34,60},
	{427,1,6,41,67},
	{547,1,7,48,85},
	{708,1,7,48,126},
	{873,1,9,56,155},
	{1094,1,8,61,164},
	{1383,1,12,65,240},
	{1650,1,11,78,216},
	{1935,1,11,90,252},
	{2304,1,16,73,338},
	{2782,1,10,99,360},
	{3324,1,16,103,488},
	{3812,1,16,103,488},
	{4368,1,12,121,561},
	{5130,1,14,142,659},
	{5892,1,16,163,757},
	{6745,1,20,149,860},
	{7880,1,16,194,734},
	{8913,1,21,177,1006},
	{9919,1,21,177,1006},
	{11081,1,19,230,870},
	{12376,1,18,254,969},
	{13932,1,25,211,1410},
	{15657,1,25,236,1585},
	{17242,1,25,236,1585},
	{18892,1,24,225,1734},
	{21061,1,28,264,1773},
	{23445,1,22,355,1700},
	{25553,1,29,303,2346},
	{27978,1,22,355,2361},
	{31347,1,30,343,2634},
	{33981,1,30,343,2634},
	{36806,1,31,353,3092},
	{39914,1,27,465,2692},
	{43592,1,34,389,3376},
	{47536,1,34,423,3682},
	{51218,1,34,423,3682},
	{54900,1,28,564,3261},
	{59702,1,37,460,4004},
	{63891,1,38,473,4590},
	{69362,1,38,509,4986},
	{74348,1,38,509,4986},
	{81303,1,39,563,5448},
	{86751,1,39,563,5448},
	{92199,1,39,563,5448},
	{97836,1,41,630,6147}};

    if (s<=54) {
	points.push_back(kFourLow[s-1][1]);
	points.push_back(kFourLow[s-1][2]);
	points.push_back(kFourLow[s-1][3]);
	points.push_back(kFourLow[s-1][4]);
	return kFourLow[s-1][0];
    }

	// ON EXTREMAL h-BASES A4
	// S Mossige
	// Mathematica Scandinavica, 1987
    points.push_back(__St_One);
    bint a(s); a/=12;
    bint a2(a); a2*=9;
    bint a3(a); a3*=27; a3+=22; a3*=a;               // 27a^2 + 22a
    bint a4(a); a4*=54; a4+=62; a4*=a; a4+=7; a4*=a; // 54a^3 + 62a^2 + 7a
	// 162a^4+186a^3+30a^2+8a-2
    bint N(a); N*=162; N+=186; N*=a; N+=30; N*=a; N+=8; N*=a; N-=2;
    bint max(s%12); max*=a4; max+= N;
    points.push_back(std::move(a2));
    points.push_back(std::move(a3));
    points.push_back(std::move(a4));
    if (verbose>0) {
        std::clog << "#[k4] est.: " << max << ", points: ";
        for(const auto& it: points) std::clog << it << ' ';
        std::clog << std::endl;
    }
    return (approx?max:Cover(points, s, verbose));
}


#define __GSTAMPS_kFive_sMAX_ 90
inline bint kFive(std::vector<bint>& points, const size_t s,
		  const int verbose) {
    assert(s>0); points.resize(0); points.reserve(5);
    assert(s<=__GSTAMPS_kFive_sMAX_);
    const static size_t kFiveLow[__GSTAMPS_kFive_sMAX_][6]={
{5,1,2,3,4,5},
{16,1,3,5,7,8},
{36,1,4,6,14,15},
{70,1,3,11,15,32},
{126,1,4,9,31,51},
{216,1,7,12,43,52},
{345,1,8,11,64,102},
{512,1,9,15,78,115},
{797,1,9,23,108,181},
{1055,1,8,27,119,194},
{1475,1,10,34,165,270},
{2047,1,10,26,195,320},
{2659,1,13,34,242,409},
{3403,1,11,48,278,720},
{4422,1,14,50,325,782},
{5629,1,14,61,381,984},
{6865,1,13,67,326,1191},
{8669,1,14,75,500,1306},
{10835,1,14,89,523,1892},
{12903,1,14,102,589,1912},
{15785,1,14,88,727,2060},
{18801,1,18,97,858,2156},
{22456,1,20,91,894,3330},
{26469,1,16,148,843,3894},
{31108,1,16,148,975,4554},
{36949,1,22,136,1168,4227},
{42744,1,22,162,1372,4889},
{49436,1,25,139,1510,5657},
{57033,1,23,170,1610,5811},
{66771,1,24,201,1718,7596},
{75558,1,23,192,1976,7018},
{86303,1,25,180,1916,8793},
{96852,1,28,202,2150,9867},
{110253,1,29,209,2434,11256},
{123954,1,27,231,2495,11464},
{140688,1,30,227,2839,12993},
{158389,1,31,234,2926,13391},
{178811,1,30,275,2947,16472},
{197293,1,29,300,3671,16677},
{223580,1,29,266,3382,18856},
{247194,1,32,294,3739,20847},
{273443,1,34,325,4133,23063},
{300747,1,33,342,4560,25414},
{331461,1,32,393,4562,25751},
{368894,1,32,457,5137,28671},
{401350,1,37,421,5602,31205},
{443231,1,33,336,5224,34431},
{490325,1,36,515,6304,35400},
{536399,1,34,422,6065,38741},
{586322,1,42,444,6906,45542},
{634430,1,36,482,7132,40052},
{699698,1,35,570,6602,50446},
{754166,1,40,462,7666,50721},
{823136,1,39,474,7840,51893},
{892139,1,42,511,8494,56238},
{968914,1,39,489,8580,65052},
{1052562,1,43,617,10091,66380},
{1150377,1,46,606,9531,72397},
{1236682,1,44,552,10237,77846},
{1325927,1,41,631,11205,74232},
{1420882,1,44,623,10432,89278},
{1547688,1,49,646,12050,91649},
{1678695,1,49,664,12338,93848},
{1782370,1,52,705,13100,99644},
{1888725,1,48,698,12988,111755},
{2036874,1,48,746,13252,113747},
{2165553,1,51,793,14087,120914},
{2330896,1,47,1000,16255,123331},
{2496702,1,53,752,16196,139500},
{2653201,1,52,789,15540,133254},
{2846834,1,53,804,16640,143028},
{3047485,1,57,866,17811,153105},
{3250580,1,55,910,17943,171807},
{3429203,1,58,960,18929,181248},
{3629795,1,59,1013,18977,182285},
{3864527,1,63,980,20223,194059},
{4103963,1,56,906,20468,196209},
{4416370,1,59,1013,21001,201513},
{4643287,1,62,1065,22080,211867},
{4975426,1,65,1141,23658,227046},
{5223883,1,63,1145,23756,227980},
{5519971,1,66,1200,24898,263837},
{5796515,1,62,1150,25011,264322},
{6139689,1,60,1090,26793,256691},
{6513282,1,60,1090,26853,284084},
{6912409,1,69,1255,28480,301525},
{7258582,1,63,1168,29884,316651},
{7677138,1,72,1336,31584,334950},
{8029729,1,67,1285,31617,366107},
{8525267,1,70,1412,32085,371775},
};
    points.assign(std::begin(kFiveLow[s-1])+1,std::end(kFiveLow[s-1]));
    return kFiveLow[s-1][0];
}

#define __GSTAMPS_kSix_sMAX_ 26
inline bint kSix(std::vector<bint>& points, const size_t s,
		 const int verbose) {
    assert(s>0); points.resize(0); points.reserve(6);
    assert(s<=__GSTAMPS_kSix_sMAX_);
    const static size_t kSixLow[__GSTAMPS_kSix_sMAX_][7]={
{6,1,2,3,4,5,6},
{20,1,2,5,8,9,10},
{52,1,3,7,9,19,24},
{108,1,4,9,16,38,49},
{211,1,4,13,24,56,61},
{388,1,7,11,48,83,115},
{664,1,7,12,64,113,193},
{1045,1,9,14,65,170,297},
{1617,1,6,31,48,256,373},
{2510,1,9,31,96,366,411},
{3607,1,7,41,105,490,815},
{5118,1,6,47,120,565,946},
{7066,1,10,35,133,759,1304},
{9748,1,11,49,188,810,2109},
{12793,1,8,71,192,1215,1993},
{17061,1,15,49,285,1292,3043},
{22342,1,13,82,387,1723,4789},
{28874,1,13,94,354,1968,5062},
{36560,1,16,87,408,2351,6452},
{45754,1,17,93,436,2898,6897},
{57814,1,14,129,469,3585,8757},
{72997,1,17,109,624,3998,9618},
{87555,1,12,117,541,4487,11496},
{106888,1,19,138,782,5346,13991},
{129783,1,19,157,896,5656,19313},
{156744,1,19,177,816,6708,18060},
};
    points.assign(std::begin(kSixLow[s-1])+1,std::end(kSixLow[s-1]));
    return kSixLow[s-1][0];
}

#define __GSTAMPS_kSeven_sMAX_ 14
inline bint kSeven(std::vector<bint>& points, const size_t s,
		   const int verbose) {
    assert(s>0); points.resize(0); points.reserve(7);
    assert(s<=__GSTAMPS_kSeven_sMAX_);
    const static size_t kSevenLow[__GSTAMPS_kSeven_sMAX_][8]={
{6,1,2,3,4,5,6},
{20,1,2,5,8,9,10},
{52,1,3,7,9,19,24},
{108,1,4,9,16,38,49},
{211,1,4,13,24,56,61},
{388,1,7,11,48,83,115},
{1137,1,7,18,62,104,244,259},
{2001,1,6,28,47,127,412,602},
{3191,1,7,30,86,189,607,920},
{5047,1,6,29,96,246,857,1179},
{7820,1,10,34,153,380,1342,1487},
{11568,1,8,49,127,419,1566,2604},
{17178,1,12,40,223,544,2479,3253},
{24466,1,12,52,225,546,3033,5464},
};
    points.assign(std::begin(kSevenLow[s-1])+1,std::end(kSevenLow[s-1]));
    return kSevenLow[s-1][0];
}

#define __GSTAMPS_kEight_sMAX_ 8
inline bint kEight(std::vector<bint>& points, const size_t s,
		   const int verbose) {
    assert(s>0); points.resize(0); points.reserve(8);
    assert(s<=__GSTAMPS_kEight_sMAX_);
    const static size_t kEightLow[__GSTAMPS_kEight_sMAX_][9]={
{8,1,2,3,4,5,6,7,8},
{32,1,2,5,8,11,14,15,16},
{93,1,3,6,10,24,26,39,41},
{228,1,3,8,19,33,39,92,102},
{524,1,6,8,33,48,77,183,236},
{1007,1,5,18,29,97,170,219,308},
{1911,1,4,17,31,117,209,513,550},
{3485,1,8,13,58,169,295,831,1036},
};
    points.assign(std::begin(kEightLow[s-1])+1,std::end(kEightLow[s-1]));
    return kEightLow[s-1][0];
}

#define __GSTAMPS_sTwo_kMAX_ 24
inline bint sTwo(std::vector<bint>& points, const size_t k,
		 const int verbose) {
    assert(k>0); points.resize(0); points.reserve(k);
    assert(k<=__GSTAMPS_sTwo_kMAX_);
    const static std::vector<bint> sTwoLow[__GSTAMPS_sTwo_kMAX_]={
{2,1},
{4,1,2},
{8,1,3,4},
{12,1,3,5,6},
{16,1,3,5,7,8},
{20,1,2,5,8,9,10},
{26,1,2,5,8,11,12,13},
{32,1,2,5,8,11,14,15,16},
{40,1,3,4,9,11,16,17,19,20},
{46,1,2,3,7,11,15,19,21,22,24},
{54,1,3,4,9,11,16,18,23,24,26,27},
{63,1,3,4,9,11,16,21,23,28,29,31,32},
{72,1,3,4,9,11,16,20,25,27,32,33,35,36},
{80,1,3,4,5,8,14,20,26,32,35,36,37,39,40},
{92,1,3,4,5,8,14,20,26,32,38,41,42,43,45,46},
{104,1,3,4,5,8,14,20,26,32,38,44,47,48,49,51,52},
{116,1,3,4,5,8,14,20,26,32,38,44,50,53,54,55,57,58},
{128,1,3,4,5,8,14,20,26,32,38,44,50,56,59,60,61,63,64},
{140,1,3,4,5,8,14,20,26,32,38,44,50,56,62,65,66,67,69,70},
{152,1,3,4,5,8,14,20,26,32,38,44,50,56,62,68,71,72,73,75,76},
{164,1,3,4,6,10,13,15,21,29,37,45,53,61,67,69,72,76,78,79,81,82},
{180,1,3,4,6,10,13,15,21,29,37,45,53,61,69,75,77,80,84,86,87,89,90},
{196,1,3,4,6,10,13,15,21,29,37,45,53,61,69,77,83,85,88,92,94,95,97,98},
{212,1,3,4,6,10,13,15,21,29,37,45,53,61,69,77,85,91,93,96,100,102,103,105,106},
};
    points.assign(sTwoLow[k-1].begin()+1,sTwoLow[k-1].end());
    return sTwoLow[k-1].front();
}


#define __GSTAMPS_sThree_kMAX_ 15
inline bint sThree(std::vector<bint>& points, const size_t k,
		   const int verbose) {
    assert(k>0); points.resize(0); points.reserve(k);
    assert(k<=__GSTAMPS_sThree_kMAX_);
    const static std::vector<bint> sThreeLow[__GSTAMPS_sThree_kMAX_]={
{3,1},
{7,1,3},
{15,1,4,5},
{24,1,4,7,8},
{36,1,4,6,14,15},
{52,1,3,7,9,19,24},
{70,1,4,5,15,18,27,34},
{93,1,3,6,10,24,26,39,41},
{121,1,3,8,9,14,32,36,51,53},
{154,1,2,6,8,19,28,40,43,91,103},
{186,1,2,3,8,11,26,38,56,69,85,89},
{225,1,3,8,13,15,16,49,53,84,88,108,114},
{271,1,4,6,14,16,20,39,56,79,100,113,122,131},
{323,1,2,4,9,15,27,38,43,46,97,107,127,147,157},
{385,1,4,5,16,18,29,37,43,52,78,98,148,158,178,188},
};
    points.assign(sThreeLow[k-1].begin()+1,sThreeLow[k-1].end());
    return sThreeLow[k-1].front();
}

#define __GSTAMPS_sFour_kMAX_ 12
inline bint sFour(std::vector<bint>& points, const size_t k,
		  const int verbose) {
    assert(k>0); points.resize(0); points.reserve(k);
    assert(k<=__GSTAMPS_sFour_kMAX_);
    const static std::vector<bint> sFourLow[__GSTAMPS_sFour_kMAX_]={
{4,1},
{10,1,3},
{26,1,5,8},
{44,1,3,11,18},
{70,1,3,11,15,32},
{108,1,4,9,16,38,49},
{162,1,4,9,24,35,49,51},
{228,1,3,8,19,33,39,92,102},
{310,1,5,7,22,31,36,83,117,133},
{422,1,4,9,24,26,42,104,115,174,185},
{550,1,4,9,20,34,52,62,137,149,229,242},
{700,1,5,8,20,22,29,45,106,174,240,311,321},
};
    points.assign(sFourLow[k-1].begin()+1,sFourLow[k-1].end());
    return sFourLow[k-1].front();
}

#define __GSTAMPS_sFive_kMAX_ 10
inline bint sFive(std::vector<bint>& points, const size_t k,
		  const int verbose) {
    assert(k>0); points.resize(0); points.reserve(k);
    assert(k<=__GSTAMPS_sFive_kMAX_);
    const static std::vector<bint> sFiveLow[__GSTAMPS_sFive_kMAX_]={
{5,1},
{14,1,4},
{35,1,6,7},
{71,1,4,12,21},
{126,1,4,9,31,51},
{211,1,4,13,24,56,61},
{336,1,4,13,24,30,87,106},
{524,1,6,8,33,48,77,183,236},
{726,1,4,13,18,51,92,163,208,223},
{1016,1,6,8,21,60,93,104,154,378,414},
};
    points.assign(sFiveLow[k-1].begin()+1,sFiveLow[k-1].end());
    return sFiveLow[k-1].front();
}

#define __GSTAMPS_sSix_kMAX_ 9
inline bint sSix(std::vector<bint>& points, const size_t k,
		 const int verbose) {
    assert(k>0); points.resize(0); points.reserve(k);
    assert(k<=__GSTAMPS_sSix_kMAX_);
    const static std::vector<bint> sSixLow[__GSTAMPS_sSix_kMAX_]={
{6,1},
{18,1,4},
{52,1,7,12},
{114,1,4,19,33},
{216,1,7,12,43,52},
{388,1,7,11,48,83,115},
{638,1,4,18,31,104,145,170},
{1007,1,5,18,29,97,170,219,308},
{1545,1,6,10,32,77,114,284,447,471},
};
    points.assign(sSixLow[k-1].begin()+1,sSixLow[k-1].end());
    return sSixLow[k-1].front();
}

// ============================================
// Hybrid methods

// Mrose Divide & Conquer
template<typename List, typename stype_t>
inline bint CutSelect(List& points, const size_t k, const size_t kotwo,
		      const stype_t s, const stype_t sotwo,
		      const int rlevel, const bool approx, const int verbose) {

    assert(kotwo>0); assert(sotwo>0);
    const stype_t smh(s-sotwo); const size_t kmo(k-kotwo);
    assert(smh>0); assert(kmo>0);

    bint m1 = FSelect(points, kmo, smh, rlevel, approx, verbose-1);
    if (verbose>0)
        rangeprint(std::clog << "#[CS] (" << k << ',' << size_t(s) << ")m1["
                   << kmo << ',' << size_t(smh) << "]:" << m1
                   << ", n: ", points) << std::endl;

    std::vector<bint> p2;
    bint m2(m1);
    if ( (kmo!=kotwo) || (sotwo!=smh) ) {
        m2 = FSelect(p2, kotwo, sotwo, rlevel, approx, verbose-1);
    } else {
        p2.assign(points.begin(),points.end());
    }
    ++m1;
    const bint max(m1*(m2+1)-1);
    for( auto& it : p2 )
        points.emplace_back( std::move( it *= m1 ) );

    if (verbose>0) {
        rangeprint(std::clog << "#[CS] (" << k << ',' << size_t(s) << ")m2["
                   << kotwo << ',' << size_t(sotwo) << "]:" << m2
                   << ", n: ", p2) << std::endl;
        std::clog << "#[CS] >= " << max <<std::endl;
    }

    return (approx?max:Cover(points, s, verbose));
}


// Recursive (rlevel) Quadratic exploration, or midpoints only
template<typename List, typename stype_t>
inline bint RecSelect(List& points, const size_t k, const stype_t s,
		      const int rlevel, const bool approx, const int verbose) {
    assert(k>1); assert(s>1);

    bint max(0);
    if (rlevel>0) {
        for(size_t klow(1); klow<k; ++klow) {
            for(stype_t slow(1); slow<s; ++slow) {
                std::vector<bint> p2;
                const bint cm = CutSelect(p2, k, klow, s, slow,
                                          rlevel-1, approx, verbose-1);
                if (verbose>0) std::clog << "#[RS] (" << k << '|' << klow << ','
                                         << size_t(s) << '|' << size_t(slow)
                                         << "): " << cm << std::endl;
                if (cm>max) {
                    points.swap(p2);
                    max = cm;
                }
            }
        }
        if (verbose>0) rangeprint(std::clog << "#[RS] max: " << max
                                  << ", points: ", points) << std::endl;

        if (! approx) {
                // Just try the previous one also
            std::vector<bint> pm;
            bint mpm = FSelect(pm,k-1,s,rlevel-1,approx,verbose-1);
            pm.push_back(mpm+1);
            mpm = _Cover(pm.begin(), pm.end(), s);
            if (verbose>0)
                std::clog << "#[FPM] (" << k-1 << ',' << size_t(s) << "): "
                          << mpm << std::endl;
            if (mpm>max) {
                points.swap(pm);
                max = mpm;
            }
        }
    } else {
            // Chossing midpoints only
        max = CutSelect(points, k, k>>1, s, stype_t(s>>1), 0, approx, verbose-1);
    }
    return max;
}

#include <gstamps_basis.h>

// Switching among different solutions, known extremal first
template<typename stype_t>
inline bint DSelect(std::vector<bint>& points, const size_t k, const stype_t s,
		    const int rlevel, const bool approx, const int verbose) {
    points.resize(0); points.reserve(k);

    if (k == 1) {
	points.push_back(__St_One);
	return bint(s);
    }

    if (s == 1) {
	for(size_t e(1); e<=k; ++e) {
	    points.emplace_back(e);
	}
	return bint(k);
    }

    if (k == 2) {
	points.push_back(__St_One);
	const bint t(s>>1);
	if (s & 0x1) {
	    points.push_back(t+2);
	    return (t*(t+4)+2);
	} else {
	    points.push_back(t+1); // t+2 could workd also ...
	    return t*(t+3);
	}
    }

    if (k == 3) {
	const bint max = kThree(points,s,verbose-1);
	if (verbose>0) rangeprint(std::clog << "#[Fk3] max: " << max
				  << ", points: ", points) << std::endl;
	return max;
    }

    bint max(0);
    if (k == 4) {
	max = kFour(points,s,approx,verbose-1);
	if (verbose>0) rangeprint(std::clog << "#[Fk4] max: " << max
				  << ", points: ", points) << std::endl;
	if (s <= __GSTAMPS_kFour_sKNOWN_) return max;
    }

    if ((s == 2) && (k>4)) {
	if (k<=__GSTAMPS_sTwo_kMAX_) {
	    max = sTwo(points,k,verbose-1);
	    if (verbose>0) rangeprint(std::clog << "#[Fs2] max: " << max
				      << ", points: ", points) << std::endl;
	    return max;
	}
	std::vector<bint> p2;
	const bint ctwo = KloveMossige(p2,k,verbose-1);
	if (verbose>0) rangeprint(std::clog << "#[FKloveMossige] max: "
				  << ctwo << ", points: ", p2) << std::endl;
	if (ctwo>max) {
	    points.swap(p2);
	    max = ctwo;
	}
    }

    if ((s == 3) && (k>4) && (k<=__GSTAMPS_sThree_kMAX_)) {
	max = sThree(points,k,verbose-1);
	if (verbose>0) rangeprint(std::clog << "#[Fs3] max: " << max
				  << ", points: ", points) << std::endl;
	return max;
    }

    if ((s == 4) && (k>4) && (k<=__GSTAMPS_sFour_kMAX_)) {
	max = sFour(points,k,verbose-1);
	if (verbose>0) rangeprint(std::clog << "#[Fs4] max: " << max
				  << ", points: ", points) << std::endl;
	return max;
    }

    if ((s == 5) && (k>4) && (k<=__GSTAMPS_sFive_kMAX_)) {
	max = sFive(points,k,verbose-1);
	if (verbose>0) rangeprint(std::clog << "#[Fs5] max: " << max
				  << ", points: ", points) << std::endl;
	return max;
    }

    if ((s == 6) && (k>4) && (k<=__GSTAMPS_sSix_kMAX_)) {
	max = sSix(points,k,verbose-1);
	if (verbose>0) rangeprint(std::clog << "#[Fs6] max: " << max
				  << ", points: ", points) << std::endl;
	return max;
    }

    if ((k == 5) && (s>6) && (s<=__GSTAMPS_kFive_sMAX_)) {
	max = kFive(points,s,verbose-1);
	if (verbose>0) rangeprint(std::clog << "#[Fk5] max: " << max
				  << ", points: ", points) << std::endl;
	return max;
    }

    if ((k == 6) && (s>6) && (s<=__GSTAMPS_kSix_sMAX_)) {
	max = kSix(points,s,verbose-1);
	if (verbose>0) rangeprint(std::clog << "#[Fk6] max: " << max
				  << ", points: ", points) << std::endl;
	return max;
    }

    if ((k == 7) && (s>6) && (s<=__GSTAMPS_kSeven_sMAX_)) {
	max = kSeven(points,s,verbose-1);
	if (verbose>0) rangeprint(std::clog << "#[Fk7] max: " << max
				  << ", points: ", points) << std::endl;
	return max;
    }

    if ((k == 8) && (s>6) && (s<=__GSTAMPS_kEight_sMAX_)) {
	max = kEight(points,s,verbose-1);
	if (verbose>0) rangeprint(std::clog << "#[Fk8] max: " << max
				  << ", points: ", points) << std::endl;
	return max;
    }

    const stype_t sot(s>>1);
    if ((k<s) && (k>sot)) {
	std::vector<bint> p2;
	const bint mab = AlterBernett(p2,k,s,verbose-1);
	if (verbose>0) rangeprint(std::clog << "#[FAlterBernett] max: " << mab
				  << ", points: ", p2) << std::endl;
	if (mab>max) {
	    points.swap(p2);
	    max = mab;
	}
    }

    if (k>=s) {
	std::vector<bint> p2;
	const bint mdp = BalGreedy(p2,k,s,verbose-1);
	if (verbose>0) rangeprint(std::clog << "#[FBalGreedy] max: " << mdp
				  << ", points: ", p2) << std::endl;
	if (mdp>max) {
	    points.swap(p2);
	    max = mdp;
	}
    }

    std::vector<bint> p2;
    const bint mrs = RecSelect(p2,k,s,rlevel,approx,verbose-1);
    if (verbose>0) rangeprint(std::clog << "#[FSRec] max: " << mrs
			      << ", points: ", p2) << std::endl;
    if (mrs>max) {
	points.swap(p2);
	max = mrs;
    }


    const static StampDB _stampDB;
    std::vector<bint> pDB;
    if (_stampDB(std::pair<size_t,size_t>(k,s),pDB)) {
	if (verbose>0) rangeprint(std::clog << "#[GoodBasis] : "
				  , pDB) << std::endl;
	if (pDB.front()>max) {
	    max = pDB.front();
	    points.assign(std::next(pDB.begin()), pDB.end());
	}
    }


    return max;
}



// Memoization of solutions
template<typename stype_t>
inline bint FSelect(std::vector<bint>& points, const size_t k, const stype_t s,
                    const int rlevel, const bool approx, const int verbose) {
    points.resize(0); points.reserve(k);
    bint max(0);
    static std::map<std::pair<size_t,size_t>,
        std::pair<bint,std::vector<bint>>> memoize;
    std::pair<size_t,size_t> p(k,s);
    if (memoize.count(p)>0) {
        max = memoize[p].first;
        auto& vec(memoize[p].second);
        points.assign(vec.begin(),vec.end());
        if (verbose>0) {
            rangeprint(std::clog << "#[FMM] (" << k << ',' << size_t(s) << "):"
                       << max << ", n: ", points)<<std::endl;
        }
    } else {
        max = DSelect(points, k, s, rlevel, approx, verbose);
        memoize[p]=std::pair<size_t,std::vector<bint>>(max,points);
    }

    return max;
}


// Exhaust additional denominations ...
#ifndef __GSTAMPS_Half_search__
#  define __GSTAMPS_AMX(a,p) (a+__St_One)
#else
	// Search only within the largest half
#  define __GSTAMPS_AMX(a,p) ( ((p+a)>>1) + __St_One )
#endif

template<typename stype_t>
inline bint complement(std::vector<bint>& prescribed,
                       const size_t k, const stype_t s, const int verbose) {

    bint bmax(_Cover(prescribed.begin(),prescribed.end(),s));
    const bint pc(bmax+__St_One);
    const bint amx( __GSTAMPS_AMX(prescribed.back(), bmax) );
    if (verbose>0) {
        std::clog << "#[Cpmt(" << prescribed.size() << ")] amx: " << amx
                  << " bm: " << bmax << " pc: " << pc << " with prescribed: ";
        for(const auto& it: prescribed) std::clog << it << ' ';
        std::clog << std::endl;
    }

    if ( (prescribed.size()>=k) || (amx>pc) ) return bmax;

    std::vector<bint> points; points.reserve(k);
    points.assign(prescribed.begin(), prescribed.end());

    std::vector<bint> bfound; bfound.reserve(k);

    for(bint u(pc); u>=amx; --u) {
        points.push_back(u);
        const bint bu = complement(points, k, s, verbose-1);
        if (bu > bmax) {
            bfound.resize(0);
            bfound.assign(points.begin(), points.end());
            bmax = bu;
            if (verbose>0) {
                std::clog << "#[Cpmt(" << prescribed.size() << ")] max: "
                          << bmax << " with basis: ";
                for(const auto& it: points) std::clog << it << ' ';
                std::clog << std::endl;
            }
        }
        points.resize(prescribed.size());
        if (verbose>0) std::clog << "#[Cpmt(" << prescribed.size() << ")] "
                                 << amx << " <= " << u << " <= " << pc
                                 << " : " << bu << " <= " << bmax << std::endl;
    }

    prescribed.resize(0);
    prescribed.assign(bfound.begin(),bfound.end());

    return bmax;
}


// Exhaust additional denominations ...
template<typename stype_t>
inline bint par_complement(std::vector<bint>& prescribed,
                           const size_t k, const stype_t s, const int verbose) {

    assert(prescribed.size()>=1);
    assert(prescribed.size()<k);

    bint bmax(_Cover(prescribed.begin(),prescribed.end(),s));
    const bint pc(bmax + __St_One);
    const bint amx( __GSTAMPS_AMX(prescribed.back(), pc) );

    std::clog << "#[PCt(" << prescribed.size() << ")]"
              << " amx: " << amx << " pc: " << pc << std::endl;

    std::vector<bint> bfound; bfound.reserve(k);

    const int64_t maxu(pc-amx); // Should not loop more than 2^63 anyway ...

    if (maxu>=0) {

#pragma omp parallel for shared(prescribed,bfound,bmax,amx,pc,k,s,verbose) schedule(dynamic)
        for(int64_t iu=maxu; iu>=0; --iu) {
            std::vector<bint> points; points.reserve(k);
            points.assign(prescribed.begin(), prescribed.end());
            points.push_back(amx+bint(iu));
            const bint bu = complement(points, k, s, verbose-1);
#pragma omp critical
            {
                if (bu > bmax) {
                    bfound.resize(0);
                    bfound.assign(points.begin(), points.end());
                    bmax = bu;
                    if (verbose>0) {
                        std::clog << "#[PCt(" << prescribed.size()
                                  << ")] max: " << bmax << " with basis: ";
                        for(const auto& it: points) std::clog << it << ' ';
                        std::clog << std::endl;
                    }
                }
                if (verbose>0)
                    std::clog << "#[PCt(" << prescribed.size() << ")] " << amx
                              << " <= " << points[prescribed.size()] << " <= "
                              << pc << " : " << bu << " <= " << bmax
                              << std::endl;
            }
        }

        prescribed.resize(0);
        prescribed.assign(bfound.begin(),bfound.end());
    }

    return bmax;
}

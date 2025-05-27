// ==========================================================================
// GStamps: C++ routines for the Global Postage Stamp Problem
// Authors: 
//   L. Colisson, J-G. Dumas, A. Galan, B. Grenet, A. Maignan, D. S. Roche
// ==========================================================================

/****************************************************************
 * GStamps Library, inline implementations
 ****************************************************************/


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
            mpm = _KCover(pm.begin(), pm.end(), s);
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

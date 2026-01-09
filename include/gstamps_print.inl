// ==========================================================================
// GStamps: C++ routines for the Local & the Global Postage Stamp Problem
// Authors:
//   L. Colisson, J-G. Dumas, A. Galan, B. Grenet, A. Maignan, D. S. Roche
// ==========================================================================

/****************************************************************
 * GStamps Library, inline implementations
 ****************************************************************/

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

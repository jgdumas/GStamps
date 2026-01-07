// ==========================================================================
// GStamps: C++ routines for the Local & the Global Postage Stamp Problem
// Authors:
//   L. Colisson, J-G. Dumas, A. Galan, B. Grenet, A. Maignan, D. S. Roche
// ==========================================================================

#ifndef _GSTAMPS_CHRONO_H_
#define _GSTAMPS_CHRONO_H_

#include <iostream>
#include <time.h>


#  ifndef CLOCKTYPE
#    ifdef CLOCK_PROCESS_CPUTIME_ID
    /*  cpu time in the current process */
#      define CLOCKTYPE  CLOCK_PROCESS_CPUTIME_ID
    /* real time in the current process */
// #    define CLOCKTYPE  CLOCK_REALTIME
#    else
    /* this one should be appropriate to
               avoid errors on multiprocessors systems */
#      define CLOCKTYPE  CLOCK_MONOTONIC
#    endif
#  endif
#  define __St_NANO_FACTOR 1.0e9



// Simple chonometer class with start/stop and access/print to the (double) time
class Chrono {
    double ctime;
    struct timespec begin_time,end_time;
public:
    void start() { clock_gettime(CLOCKTYPE, &begin_time); }
    double stop() {
        clock_gettime(CLOCKTYPE, &end_time);
        ctime = difftime(end_time.tv_sec, begin_time.tv_sec);
        return ctime += ((double) (end_time.tv_nsec - begin_time.tv_nsec) )/ __St_NANO_FACTOR;
    }
    const double usertime() {
        return this->ctime;
    }
    const double operator-(const Chrono& c2) {
        return this->ctime-c2.ctime;
    }
    const double operator+(const Chrono& c2) {
        return this->ctime+c2.ctime;
    }
    friend std::ostream& operator<<(std::ostream& o, const Chrono& c) {
        return o << c.ctime << 's';
    }
};


#endif

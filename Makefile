#########################################################################
# GStamps: C++ routines for the Local & the Global Postage Stamp Problem
# Authors:
#   L. Colisson, J-G. Dumas, A. Galan, B. Grenet, A. Maignan, D. S. Roche
#########################################################################


OPTFLAGS = -O3 -ffast-math

#######
# g++ options
OPTFLAGS += -fopenmp
OPTFLAGS += -fprefetch-loop-arrays -funroll-all-loops
OPTFLAGS += -UDEBUG -DNDEBUG
#OPTFLAGS +=  -D__GSTAMPS_EXTENDED_PRECISION
#######


CXXFLAGS += ${OPTFLAGS} -I`pwd`/include/ `pkg-config givaro --cflags`
LOADLIBES+= `pkg-config givaro --libs`

#######

PRG  = range basis
PRG += dynprg supplement complement
PRG += brute
PRG += search

BEN  = fibo alba geom bala
BEN += krange reach srange
BEN += depthrange

#######

EXE = ${PRG} ${BEN}
BIN = ${EXE:%=bin/%}

all: ${BIN}

prg: ${PRG:%=bin/%}

bench: ${BEN:%=bin/%}

VPATH = src:benchmarks

bin/%: %.cpp
	$(LINK.cpp) $^ $(LOADLIBES) $(LDLIBS) -o $@

clean:
	- \rm ${BIN}

range: FDTC.sh ${BIN}
	./$< 5 3 1 3 1 6
	./$< 5 5 0 5 1 6
	./$< 30 3 1 3 20 6
	./$< 50 8 1 3 45 5
	./$< 5 50 0 5 45 4

basis: FDTB.sh ${BIN}
	./$< 15

brute: FDTA.sh ${BIN}
	./$< 10

check: brute basis range

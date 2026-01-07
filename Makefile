#########################################################################
# GStamps: C++ routines for the Global Postage Stamp Problem
# Authors:
#   L. Colisson, J-G. Dumas, A. Galan, B. Grenet, A. Maignan, D. S. Roche
#########################################################################


OPTFLAGS = -O3 -ffast-math

#######
# g++ options
OPTFLAGS += -fopenmp
OPTFLAGS += -fprefetch-loop-arrays -funroll-all-loops
OPTFLAGS += -UDEBUG -DNDEBUG
#######


CXXFLAGS += ${OPTFLAGS} -I`pwd`/include/ `pkg-config givaro --cflags`
LOADLIBES+= `pkg-config givaro --libs`

#######

PRG  = cover basis
PRG += dynprg supplement complement
PRG += brute
PRG += search

BEN  = fibo albe geom
BEN += kcover reach scover
BEN += depthcover

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

cover: FDTC.sh ${BIN}
	./$< 30 3 1 3 20 1
	./$< 50 8 1 3 45 3
	./$< 5 50 0 5 45 3

basis: FDTB.sh ${BIN}
	./$< 15

check: basis cover

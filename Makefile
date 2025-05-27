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

check: FDT.sh ${BIN}
	./$< 20

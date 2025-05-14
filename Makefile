####################################################################
# GStamps: C++ routines for the Global Postage Stamp Problem
# Authors: J-G. Dumas, A. Galan, B. Grenet, A. Maignan, D. S. Roche
####################################################################


OPTFLAGS = -O3 -ffast-math

#######
# g++ options
OPTFLAGS += -fopenmp -march=native -funroll-all-loops
OPTFLAGS += -UDEBUG -DNDEBUG
#######


CXXFLAGS += ${OPTFLAGS} -I`pwd`/include/ `pkg-config givaro --cflags`
LOADLIBES+= `pkg-config givaro --libs`

#######

EXE  = cover basis
EXE += dynprg supplement complement
EXE += brute
EXE += search
EXE += fibo ab geom


SRC=${EXE:%=src/%.cpp}

BIN=${EXE:%=bin/%}

#######

all: ${BIN}

bin/%: src/%.cpp
	$(LINK.cpp) $^ $(LOADLIBES) $(LDLIBS) -o $@

clean:
	- \rm ${BIN}

check: ${BIN}
	./bin/FDT.sh 20

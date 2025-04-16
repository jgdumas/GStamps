####################################################################
# GStamps: C++ routines for the Global Postage Stamp Problem
# Authors: J-G. Dumas, A. Galan, B. Grenet, A. Maignan, D. S. Roche
####################################################################

OPTFLAGS = -O3 -ffast-math

#######

CXXFLAGS += ${OPTFLAGS} -I`pwd`/include/ `pkg-config givaro --cflags`
LOADLIBES+= `pkg-config givaro --libs`

#######

EXE  = cover basis
EXE += brute
EXE += search
EXE += fibo ab geom
EXE += dynprg supplement complement


SRC=${EXE:%=src/%.cpp}

BIN=${EXE:%=bin/%}

#######

all: ${BIN}

bin/%: src/%.cpp
	$(LINK.cpp) $^ $(LOADLIBES) $(LDLIBS) -o $@

clean:
	- \rm ${BIN}

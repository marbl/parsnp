
parsnp_CXXFLAGS = -I$(top_srcdir)/src  -O3 -m64 -fopenmp -L.$(libmuscle)/libMUSCLE/.libs/ -funroll-all-loops -fomit-frame-pointer -ftree-vectorize
parsnp_LDFLAGS = -lgomp -lstdc++ -lpthread -std=gnu++0x -lMUSCLE-3.7

# this lists the binaries to produce, the (non-PHONY, binary) targets in
# the previous manual Makefile
bin_PROGRAMS = parsnp
parsnp_SOURCES = MuscleInterface.cpp MuscleInterface.h parsnp.cpp parsnp.hh LCB.cpp LCB.hh LCR.cpp LCR.hh TMum.cpp TMum.hh Converter.cpp Converter.hh ./ext/iniFile.cpp ./ext/iniFile.h 
bindir = $(top_srcdir)/bin

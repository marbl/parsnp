parsnp_core_CXXFLAGS = -fopenmp -O2 -m64 -funroll-all-loops -fomit-frame-pointer -ftree-vectorize 
parsnp_core_LDFLAGS =  -fopenmp -lstdc++ -lpthread -std=gnu++0x -L${CONDA_PREFIX}/lib -lMUSCLE-3.7 
bin_PROGRAMS = parsnp_core
parsnp_core_SOURCES = MuscleInterface.cpp MuscleInterface.h parsnp.cpp parsnp.hh LCB.cpp LCB.hh LCR.cpp LCR.hh TMum.cpp TMum.hh Converter.cpp Converter.hh ./ext/iniFile.cpp ./ext/iniFile.h 
bindir = $(top_srcdir)/bin
libdir = $(top_srcdir)/lib


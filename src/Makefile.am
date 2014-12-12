parsnp_CXXFLAGS = -Xlinker --no-as-needed -O3 -m64 -fopenmp -funroll-all-loops -fomit-frame-pointer -ftree-vectorize 
parsnp_LDFLAGS = -lgomp -lstdc++ -lpthread -std=gnu++0x -Wl,-rpath,$(libmuscle)/lib -L$(libmuscle)/lib -lMUSCLE-3.7 
bin_PROGRAMS = parsnp
parsnp_SOURCES = MuscleInterface.cpp MuscleInterface.h parsnp.cpp parsnp.hh LCB.cpp LCB.hh LCR.cpp LCR.hh TMum.cpp TMum.hh Converter.cpp Converter.hh ./ext/iniFile.cpp ./ext/iniFile.h 
bindir = $(top_srcdir)/bin
libdir = $(top_srcdir)/lib


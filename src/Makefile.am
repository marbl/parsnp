parsnp_CXXFLAGS = -I$(abs_top_srcdir)/src  -O3 -m64 -fopenmp -L/usr/local/lib -funroll-all-loops -fomit-frame-pointer -ftree-vectorize 
parsnp_LDFLAGS = -lgomp -lstdc++ -lpthread -std=gnu++0x -Wl,-rpath -Wl,$(abs_top_srcdir)/muscle/lib -lMUSCLE-3.7 
bin_PROGRAMS = parsnp
parsnp_SOURCES = MuscleInterface.cpp MuscleInterface.h parsnp.cpp parsnp.hh LCB.cpp LCB.hh LCR.cpp LCR.hh TMum.cpp TMum.hh Converter.cpp Converter.hh ./ext/iniFile.cpp ./ext/iniFile.h 
bindir = $(abs_top_srcdir)/bin
libdir = $(abs_top_srcdir)/lib


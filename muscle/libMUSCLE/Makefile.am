OPTIMIZATION = -O2 -funroll-loops -fomit-frame-pointer -ftree-vectorize 
AM_CXXFLAGS = $(OPTIMIZATION) -DNDEBUG=1 @OPENMP_CXXFLAGS@
AM_LDFLAGS = $(OPTIMIZATION)

LIBMUSCLE_H = \
alpha.h cluster.h clust.h clustsetdf.h clustset.h clustsetmsa.h diaglist.h distcalc.h distfunc.h dpregionlist.h dpreglist.h edgelist.h enumopts.h enums.h estring.h gapscoredimer.h gonnet.h intmath.h msadist.h msa.h muscle.h objscore.h params.h profile.h pwpath.h refine.h scorehistory.h seq.h seqvect.h textfile.h timing.h tree.h types.h unixio.h threadstorage.h 

LIBMUSCLE_SRC = \
aligngivenpath.cpp aligngivenpathsw.cpp aligntwomsas.cpp aligntwoprofs.cpp \
aln.cpp alpha.cpp anchoredpp.cpp anchors.cpp bittraceback.cpp blosumla.cpp \
clust.cpp cluster.cpp clwwt.cpp color.cpp cons.cpp diaglist.cpp \
diffobjscore.cpp diffpaths.cpp difftrees.cpp difftreese.cpp distcalc.cpp \
distfunc.cpp distpwkimura.cpp domuscle.cpp dosp.cpp dpreglist.cpp drawtree.cpp \
edgelist.cpp enumopts.cpp enumtostr.cpp estring.cpp fasta2.cpp fasta.cpp \
fastclust.cpp fastdist.cpp fastdistjones.cpp fastdistkbit.cpp fastdistkmer.cpp \
fastdistmafft.cpp fastdistnuc.cpp fastscorepath2.cpp finddiags.cpp \
finddiagsn.cpp glbalign352.cpp glbalign.cpp glbaligndiag.cpp glbalignle.cpp \
glbalignsimple.cpp glbalignsp.cpp glbalignspn.cpp glbalignss.cpp \
glbalndimer.cpp globals.cpp globalslinux.cpp globalswin32.cpp gonnet.cpp \
gotowt.cpp henikoffweight.cpp henikoffweightpb.cpp html.cpp hydro.cpp \
intmath.cpp local.cpp makerootmsab.cpp makerootmsa.cpp mhack.cpp mpam200.cpp \
msa2.cpp msa.cpp msadistkimura.cpp msf.cpp muscle.cpp muscleout.cpp nucmx.cpp nwdasimple2.cpp nwdasimple.cpp nwdasmall.cpp nwrec.cpp nwsmall.cpp objscore2.cpp objscore.cpp objscoreda.cpp onexception.cpp options.cpp outweights.cpp pam200mafft.cpp params.cpp phy2.cpp phy3.cpp phy4.cpp phy.cpp phyfromclust.cpp phyfromfile.cpp physeq.cpp phytofile.cpp posgap.cpp ppscore.cpp profdb.cpp profile.cpp profilefrommsa.cpp progalign.cpp progress.cpp progressivealign.cpp pwpath.cpp readmx.cpp realigndiffs.cpp realigndiffse.cpp refine.cpp refinehoriz.cpp refinesubfams.cpp refinetree.cpp refinetreee.cpp refinevert.cpp refinew.cpp savebest.cpp scoregaps.cpp scorehistory.cpp scorepp.cpp seq.cpp seqvect.cpp setblosumweights.cpp setgscweights.cpp setnewhandler.cpp spfast.cpp sptest.cpp stabilize.cpp subfam.cpp subfams.cpp sw.cpp termgaps.cpp textfile.cpp threewaywt.cpp traceback.cpp tracebackopt.cpp tracebacksw.cpp treefrommsa.cpp typetostr.cpp upgma2.cpp usage.cpp validateids.cpp vtml2.cpp writescorefile.cpp \
blosum62.cpp globalsother.cpp maketree.cpp scoredist.cpp tomhydro.cpp

library_includedir=$(includedir)/$(GENERIC_LIBRARY_NAME)-$(GENERIC_API_VERSION)/$(GENERIC_LIBRARY_NAME)
library_include_HEADERS = $(LIBMUSCLE_H)

INCLUDES = -I$(top_srcdir)

lib_LTLIBRARIES = libMUSCLE-3.7.la
libMUSCLE_3_7_la_SOURCES = $(LIBMUSCLE_SRC)
libMUSCLE_3_7_la_LDFLAGS= -version-info $(GENERIC_LIBRARY_VERSION) 
libMUSCLE_3_7_la_CXXFLAGS= -D_LIB $(AM_CXXFLAGS)

bin_PROGRAMS = muscle 


muscle_SOURCES = main.cpp
muscle_LDADD = $(top_builddir)/libMUSCLE/libMUSCLE-3.7.la



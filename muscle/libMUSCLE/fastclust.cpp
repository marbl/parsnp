#include "libMUSCLE/muscle.h"
#include "libMUSCLE/seqvect.h"
#include "libMUSCLE/distfunc.h"
#include "libMUSCLE/clust.h"
#include "libMUSCLE/clustsetdf.h"
#include "libMUSCLE/tree.h"
#include "libMUSCLE/clust.h"
#include "libMUSCLE/distcalc.h"
#include <math.h>

namespace muscle {

static void TreeFromSeqVect_NJ(const DistFunc &DF, CLUSTER Cluster, Tree &tree)
	{
    ClustSetDF CSD(DF);

    Clust C;
    C.Create(CSD, Cluster);

    tree.FromClust(C);
	}

static void TreeFromSeqVect_UPGMA(const DistFunc &DF, CLUSTER Cluster, Tree &tree)
	{
	LINKAGE Linkage = LINKAGE_Undefined;
	switch (Cluster)
		{
	case CLUSTER_UPGMA:
		Linkage = LINKAGE_Avg;
		break;
	case CLUSTER_UPGMAMin:
		Linkage = LINKAGE_Min;
		break;
	case CLUSTER_UPGMAMax:
		Linkage = LINKAGE_Max;
		break;
	case CLUSTER_UPGMB:
		Linkage = LINKAGE_Biased;
		break;
	default:
		Quit("TreeFromSeqVect_UPGMA, CLUSTER_%u not supported", Cluster);
		}
	
	DistCalcDF DC;
	DC.Init(DF);
	UPGMA2(DC, tree, Linkage);
	}

static void SaveDF(const SeqVect &v, DistFunc &d, const char *FileName)
	{
	FILE *f = fopen(FileName, "w");
	if (f == 0)
		Quit("Cannot create %s", FileName);

	unsigned n = v.GetSeqCount();
	fprintf(f, "%u\n", n);
	for (unsigned i = 0; i < n; ++i)
		{
		fprintf(f, "%10.10s  ", v.GetSeqName(i));
		for (unsigned j = 0; j < i; ++j)
			fprintf(f, "  %9g", d.GetDist(i, j));
		fprintf(f, "\n");
		}
	fclose(f);
	}

void TreeFromSeqVect(const SeqVect &v, Tree &tree, CLUSTER Cluster,
  DISTANCE Distance, ROOT Root, const char *SaveFileName)
	{
	DistFunc DF;
	DistUnaligned(v, Distance, DF);
	if (SaveFileName != 0)
		SaveDF(v, DF, SaveFileName);
	if (CLUSTER_NeighborJoining == Cluster)
		TreeFromSeqVect_NJ(DF, Cluster, tree);
	else
		TreeFromSeqVect_UPGMA(DF, Cluster, tree);
	FixRoot(tree, Root);
	}
} 

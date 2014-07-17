#include "libMUSCLE/muscle.h"
#include "libMUSCLE/msa.h"
#include "libMUSCLE/tree.h"
#include "libMUSCLE/clust.h"
#include "libMUSCLE/clustsetmsa.h"
#include "libMUSCLE/distcalc.h"

namespace muscle {

static void SaveMSADist(const MSA &msa, MSADist &d, const char *FileName)
	{
	FILE *f = fopen(FileName, "w");
	if (f == 0)
		Quit("Cannot create %s", FileName);

	unsigned n = msa.GetSeqCount();
	for (unsigned i = 0; i < n; ++i)
		{
		fprintf(f, "%10.10s  ", msa.GetSeqName(i));
		for (unsigned j = 0; j < n; ++j)
			fprintf(f, "  %9g", d.ComputeDist(msa, i, j));
		fprintf(f, "\n");
		}
	fclose(f);
	}

static void TreeFromMSA_NJ(const MSA &msa, Tree &tree, CLUSTER Cluster,
  DISTANCE Distance, const char *SaveFileName)
	{
	MSADist MD(Distance);
	ClustSetMSA Set(msa, MD);

	if (SaveFileName != 0)
		SaveMSADist(msa, MD, SaveFileName);

	Clust C;
	C.Create(Set, Cluster);

	tree.FromClust(C);
	}

static void SaveDC(const DistCalcMSA &DC, const char *FileName)
	{
	FILE *f = fopen(FileName, "w");
	if (f == 0)
		Quit("Cannot create %s", FileName);

	unsigned n = DC.GetCount();
	fprintf(f, "%u\n", n);
	float *Dist = new float[n];
	for (unsigned i = 0; i < n; ++i)
		{
		fprintf(f, "%10.10s  ", DC.GetName(i));
		DC.CalcDistRange(i, Dist);
		for (unsigned j = 0; j < i; ++j)
			fprintf(f, "  %9g", Dist[j]);
		fprintf(f, "\n");
		}
	fclose(f);
	}

static void TreeFromMSA_UPGMA(const MSA &msa, Tree &tree, CLUSTER Cluster,
  DISTANCE Distance, const char *SaveFileName)
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
		Quit("TreeFromMSA_UPGMA, CLUSTER_%u not supported", Cluster);
		}
	
	DistCalcMSA DC;
	DC.Init(msa, Distance);
	if (SaveFileName != 0)
		SaveDC(DC, SaveFileName);
	UPGMA2(DC, tree, Linkage);
	}

void TreeFromMSA(const MSA &msa, Tree &tree, CLUSTER Cluster,
  DISTANCE Distance, ROOT Root, const char *SaveFileName)
	{
	if (CLUSTER_NeighborJoining == Cluster)
		TreeFromMSA_NJ(msa, tree, Cluster, Distance, SaveFileName);
	else
		TreeFromMSA_UPGMA(msa, tree, Cluster, Distance, SaveFileName);
	FixRoot(tree, Root);
	}
} 

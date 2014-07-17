#include "libMUSCLE/muscle.h"
#include "libMUSCLE/msa.h"
#include "libMUSCLE/tree.h"
#include "libMUSCLE/profile.h"
#include <stdio.h>

namespace muscle {

#define TRACE	0

void RefineTreeE(MSA &msa, const SeqVect &v, Tree &tree, ProgNode *ProgNodes)
	{
	const unsigned uSeqCount = msa.GetSeqCount();
	if (tree.GetLeafCount() != uSeqCount)
		Quit("Refine tree, tree has different number of nodes");

	if (uSeqCount < 3)
		return;

#if	DEBUG
	ValidateMuscleIds(msa);
	ValidateMuscleIds(tree);
#endif

	const unsigned uNodeCount = tree.GetNodeCount();
	unsigned *uNewNodeIndexToOldNodeIndex= new unsigned[uNodeCount];

	Tree Tree2;
	TreeFromMSA(msa, Tree2, g_Cluster2.get(), g_Distance2.get(), g_Root2.get(), g_pstrDistMxFileName2.get());

#if	DEBUG
	ValidateMuscleIds(Tree2);
#endif

	DiffTreesE(Tree2, tree, uNewNodeIndexToOldNodeIndex);

	unsigned uRoot = Tree2.GetRootNodeIndex();
	if (NODE_CHANGED == uNewNodeIndexToOldNodeIndex[uRoot])
		{
		MSA msa2;
		RealignDiffsE(msa, v, Tree2, tree, uNewNodeIndexToOldNodeIndex, msa2, ProgNodes);
		tree.Copy(Tree2);
		msa.Copy(msa2);
#if	DEBUG
		ValidateMuscleIds(msa2);
#endif
		}

	delete[] uNewNodeIndexToOldNodeIndex;

	SetCurrentAlignment(msa);
	ProgressStepsDone();
	}
} 

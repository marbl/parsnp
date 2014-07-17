#include "libMUSCLE/muscle.h"
#include "libMUSCLE/msa.h"
#include "libMUSCLE/tree.h"
#include "libMUSCLE/profile.h"
#include <stdio.h>

namespace muscle {

void RefineTree(MSA &msa, Tree &tree)
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

	unsigned *IdToDiffsLeafNodeIndex = new unsigned[uSeqCount];
	unsigned uDiffsCount = uSeqCount;
	Tree Tree2;
	for (unsigned uIter = 0; uIter < g_uMaxTreeRefineIters.get(); ++uIter)
		{
		TreeFromMSA(msa, Tree2, g_Cluster2.get(), g_Distance2.get(), g_Root2.get(), g_pstrDistMxFileName2.get());

#if	DEBUG
		ValidateMuscleIds(Tree2);
#endif

		Tree Diffs;
		DiffTrees(Tree2, tree, Diffs, IdToDiffsLeafNodeIndex);

		tree.Copy(Tree2);

		const unsigned uNewDiffsNodeCount = Diffs.GetNodeCount();
		const unsigned uNewDiffsCount = (uNewDiffsNodeCount - 1)/2;

		if (0 == uNewDiffsCount || uNewDiffsCount >= uDiffsCount)
			{
			ProgressStepsDone();
			break;
			}
		uDiffsCount = uNewDiffsCount;

		MSA msa2;
		RealignDiffs(msa, Diffs, IdToDiffsLeafNodeIndex, msa2);

#if	DEBUG
		ValidateMuscleIds(msa2);
#endif

		msa.Copy(msa2);
		SetCurrentAlignment(msa);
		}

	delete[] IdToDiffsLeafNodeIndex;
	}
} 

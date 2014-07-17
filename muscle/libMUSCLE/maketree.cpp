#include "libMUSCLE/muscle.h"
#include "libMUSCLE/msa.h"
#include "libMUSCLE/textfile.h"
#include "libMUSCLE/tree.h"

namespace muscle {

void DoMakeTree()
	{
	if (g_pstrInFileName.get() == 0 || g_pstrOutFileName.get() == 0)
		Quit("-maketree requires -in <msa> and -out <treefile>");

	SetStartTime();

	SetSeqWeightMethod(g_SeqWeight1.get());

	TextFile MSAFile(g_pstrInFileName.get());

	MSA msa;
	msa.FromFile(MSAFile);

	unsigned uSeqCount = msa.GetSeqCount();
	MSA::SetIdCount(uSeqCount);

// Initialize sequence ids.
// From this point on, ids must somehow propogate from here.
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		msa.SetSeqId(uSeqIndex, uSeqIndex);
	SetMuscleInputMSA(msa);

	Progress("%u sequences", uSeqCount);

	Tree tree;
	TreeFromMSA(msa, tree, g_Cluster2.get(), g_Distance2.get(), g_Root2.get());

	TextFile TreeFile(g_pstrOutFileName.get(), true);
	tree.ToFile(TreeFile);

	Progress("Tree created");
	}
} 

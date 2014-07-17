#include "libMUSCLE/muscle.h"
#include "libMUSCLE/textfile.h"
#include "libMUSCLE/seqvect.h"
#include "libMUSCLE/distfunc.h"
#include "libMUSCLE/msa.h"
#include "libMUSCLE/tree.h"
#include "libMUSCLE/clust.h"
#include "libMUSCLE/profile.h"
#include "libMUSCLE/clustsetmsa.h"

namespace muscle {

void Refine()
	{
	SetOutputFileName(g_pstrOutFileName.get());
	SetInputFileName(g_pstrInFileName.get());
	SetStartTime();

	SetMaxIters(g_uMaxIters.get());
	SetSeqWeightMethod(g_SeqWeight1.get());

	TextFile fileIn(g_pstrInFileName.get());
	MSA msa;
	msa.FromFile(fileIn);

	const unsigned uSeqCount = msa.GetSeqCount();
	if (0 == uSeqCount)
		Quit("No sequences in input file");

	ALPHA Alpha = ALPHA_Undefined;
	switch (g_SeqType.get())
		{
	case SEQTYPE_Auto:
		Alpha = msa.GuessAlpha();
		break;

	case SEQTYPE_Protein:
		Alpha = ALPHA_Amino;
		break;

	case SEQTYPE_DNA:
		Alpha = ALPHA_DNA;
		break;

	case SEQTYPE_RNA:
		Alpha = ALPHA_RNA;
		break;

	default:
		Quit("Invalid SeqType");
		}
	SetAlpha(Alpha);
	msa.FixAlpha();

	SetPPScore();
	if (ALPHA_DNA == Alpha || ALPHA_RNA == Alpha)
		SetPPScore(PPSCORE_SPN);

	MSA::SetIdCount(uSeqCount);

// Initialize sequence ids.
// From this point on, ids must somehow propogate from here.
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		msa.SetSeqId(uSeqIndex, uSeqIndex);
	SetMuscleInputMSA(msa);

	Tree GuideTree;
	TreeFromMSA(msa, GuideTree, g_Cluster2.get(), g_Distance2.get(), g_Root2.get());
	SetMuscleTree(GuideTree);

	if (g_bAnchors.get())
		RefineVert(msa, GuideTree, g_uMaxIters.get());
	else
		RefineHoriz(msa, GuideTree, g_uMaxIters.get(), false, false);

	ValidateMuscleIds(msa);
	ValidateMuscleIds(GuideTree);

//	TextFile fileOut(g_pstrOutFileName.get(), true);
//	msa.ToFile(fileOut);
	MuscleOutput(msa);
	}
} 

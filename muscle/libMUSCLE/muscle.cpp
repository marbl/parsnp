#include "libMUSCLE/muscle.h"
#include "libMUSCLE/msa.h"
#include "libMUSCLE/seqvect.h"
#include "libMUSCLE/msa.h"
#include "libMUSCLE/tree.h"
#include "libMUSCLE/profile.h"

namespace muscle {

void MUSCLE(SeqVect &v, MSA &msaOut)
	{
	const unsigned uSeqCount = v.Length();

	if (0 == uSeqCount)
		Quit("No sequences in input file");

	ALPHA Alpha = ALPHA_Undefined;
	switch (g_SeqType.get())
		{
	case SEQTYPE_Auto:
		Alpha = v.GuessAlpha();
		break;

	case SEQTYPE_Protein:
		Alpha = ALPHA_Amino;
		break;

	case SEQTYPE_RNA:
		Alpha = ALPHA_RNA;
		break;

	case SEQTYPE_DNA:
		Alpha = ALPHA_DNA;
		break;

	default:
		Quit("Invalid seq type");
		}
	SetAlpha(Alpha);
	v.FixAlpha();

	if (ALPHA_DNA == Alpha || ALPHA_RNA == Alpha)
		{
		SetPPScore(PPSCORE_SPN);
		g_Distance1.get() = DISTANCE_Kmer4_6;
		}

	unsigned uMaxL = 0;
	unsigned uTotL = 0;
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		unsigned L = v.GetSeq(uSeqIndex).Length();
		uTotL += L;
		if (L > uMaxL)
			uMaxL = L;
		}

	SetIter(1);
	g_bDiags.get() = g_bDiags1.get();
	SetSeqStats(uSeqCount, uMaxL, uTotL/uSeqCount);

	MSA::SetIdCount(uSeqCount);

//// Initialize sequence ids.
//// From this point on, ids must somehow propogate from here.
//	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
//		v.SetSeqId(uSeqIndex, uSeqIndex);

	if (uSeqCount > 1)
		MHackStart(v);

	if (0 == uSeqCount)
		{
		msaOut.Clear();
		return;
		}

	if (1 == uSeqCount && ALPHA_Amino == Alpha)
		{
		const Seq &s = v.GetSeq(0);
		msaOut.FromSeq(s);
		return;
		}

// First iteration
	Tree GuideTree;
	TreeFromSeqVect(v, GuideTree, g_Cluster1.get(), g_Distance1.get(), g_Root1.get());

	SetMuscleTree(GuideTree);

	ProgNode *ProgNodes = 0;
	if (g_bLow.get())
		ProgNodes = ProgressiveAlignE(v, GuideTree, msaOut);
	else
		ProgressiveAlign(v, GuideTree, msaOut);
	SetCurrentAlignment(msaOut);

	if (1 == g_uMaxIters.get() || 2 == uSeqCount)
		{
		MHackEnd(msaOut);
		return;
		}

	g_bDiags.get() = g_bDiags2.get();
	SetIter(2);

	if (g_bLow.get())
		{
		if (0 != g_uMaxTreeRefineIters.get())
			RefineTreeE(msaOut, v, GuideTree, ProgNodes);
		}
	else
		RefineTree(msaOut, GuideTree);

	extern void DeleteProgNode(ProgNode &Node);
	const unsigned uNodeCount = GuideTree.GetNodeCount();
	for (unsigned uNodeIndex = 0; uNodeIndex < uNodeCount; ++uNodeIndex)
		DeleteProgNode(ProgNodes[uNodeIndex]);

	delete[] ProgNodes;
	ProgNodes = 0;

	SetSeqWeightMethod(g_SeqWeight2.get());
	SetMuscleTree(GuideTree);

	if (g_bAnchors.get())
		RefineVert(msaOut, GuideTree, g_uMaxIters.get() - 2);
	else
		RefineHoriz(msaOut, GuideTree, g_uMaxIters.get() - 2, false, false);

	MHackEnd(msaOut);
	}
} 

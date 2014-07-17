#include "libMUSCLE/muscle.h"
#include "libMUSCLE/textfile.h"
#include "libMUSCLE/seqvect.h"
#include "libMUSCLE/distfunc.h"
#include "libMUSCLE/msa.h"
#include "libMUSCLE/tree.h"
#include "libMUSCLE/profile.h"
#include "libMUSCLE/timing.h"

namespace muscle {

static const char g_strUseTreeWarning[] =
"\n******** WARNING ****************\n"
"\nYou specified the -usetree option.\n"
"Note that a good evolutionary tree may NOT be a good\n"
"guide tree for multiple alignment. For more details,\n"
"please refer to the user guide. To disable this\n"
"warning, use -usetree_nowarn <treefilename>.\n\n";

void DoMuscle()
	{
	SetOutputFileName(g_pstrOutFileName.get());
	SetInputFileName(g_pstrInFileName.get());

	SetMaxIters(g_uMaxIters.get());
	SetSeqWeightMethod(g_SeqWeight1.get());

	TextFile fileIn(g_pstrInFileName.get());
	SeqVect v;
	v.FromFASTAFile(fileIn);
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

	case SEQTYPE_DNA:
		Alpha = ALPHA_DNA;
		break;

	case SEQTYPE_RNA:
		Alpha = ALPHA_RNA;
		break;

	default:
		Quit("Invalid seq type");
		}
	SetAlpha(Alpha);
	v.FixAlpha();

//
// AED 21/12/06: Moved matrix loading code inside the PP param function so it gets called for all alignment types
//
	SetPPScore();


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

	SetMuscleSeqVect(v);

	MSA::SetIdCount(uSeqCount);

// Initialize sequence ids.
// From this point on, ids must somehow propogate from here.
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		v.SetSeqId(uSeqIndex, uSeqIndex);

	if (0 == uSeqCount)
		Quit("Input file '%s' has no sequences", g_pstrInFileName.get());
	if (1 == uSeqCount)
		{
		TextFile fileOut(g_pstrOutFileName.get(), true);
		v.ToFile(fileOut);
		return;
		}

	if (uSeqCount > 1)
		MHackStart(v);

// First iteration
	Tree GuideTree;
	if (0 != g_pstrUseTreeFileName.get())
		{
	// Discourage users...
		if (!g_bUseTreeNoWarn.get())
			fprintf(stderr, g_strUseTreeWarning);

	// Read tree from file
		TextFile TreeFile(g_pstrUseTreeFileName.get());
		GuideTree.FromFile(TreeFile);

	// Make sure tree is rooted
		if (!GuideTree.IsRooted())
			Quit("User tree must be rooted");

		if (GuideTree.GetLeafCount() != uSeqCount)
			Quit("User tree does not match input sequences");

		const unsigned uNodeCount = GuideTree.GetNodeCount();
		for (unsigned uNodeIndex = 0; uNodeIndex < uNodeCount; ++uNodeIndex)
			{
			if (!GuideTree.IsLeaf(uNodeIndex))
				continue;
			const char *LeafName = GuideTree.GetLeafName(uNodeIndex);
			unsigned uSeqIndex;
			bool SeqFound = v.FindName(LeafName, &uSeqIndex);
			if (!SeqFound)
				Quit("Label %s in tree does not match sequences", LeafName);
			unsigned uId = v.GetSeqIdFromName(LeafName);
			GuideTree.SetLeafId(uNodeIndex, uId);
			}
		}
	else
		TreeFromSeqVect(v, GuideTree, g_Cluster1.get(), g_Distance1.get(), g_Root1.get(),
		  g_pstrDistMxFileName1.get());

	const char *Tree1 = ValueOpt("Tree1");
	if (0 != Tree1)
		{
		TextFile f(Tree1, true);
		GuideTree.ToFile(f);
		if (g_bClusterOnly.get())
			return;
		}

	SetMuscleTree(GuideTree);
	ValidateMuscleIds(GuideTree);

	MSA msa;
	ProgNode *ProgNodes = 0;
	if (g_bLow.get())
		ProgNodes = ProgressiveAlignE(v, GuideTree, msa);
	else
		ProgressiveAlign(v, GuideTree, msa);
	SetCurrentAlignment(msa);

	if (0 != g_pstrComputeWeightsFileName.get())
		{
		extern void OutWeights(const char *FileName, const MSA &msa);
		SetMSAWeightsMuscle(msa);
		OutWeights(g_pstrComputeWeightsFileName.get(), msa);
		return;
		}

	ValidateMuscleIds(msa);

	if (1 == g_uMaxIters.get() || 2 == uSeqCount)
		{
		//TextFile fileOut(g_pstrOutFileName.get(), true);
		//MHackEnd(msa);
		//msa.ToFile(fileOut);
		MuscleOutput(msa);
		return;
		}

	if (0 == g_pstrUseTreeFileName.get())
		{
		g_bDiags.get() = g_bDiags2.get();
		SetIter(2);

		if (g_bLow.get())
			{
			if (0 != g_uMaxTreeRefineIters.get())
				RefineTreeE(msa, v, GuideTree, ProgNodes);
			}
		else
			RefineTree(msa, GuideTree);

		const char *Tree2 = ValueOpt("Tree2");
		if (0 != Tree2)
			{
			TextFile f(Tree2, true);
			GuideTree.ToFile(f);
			}
		}

	SetSeqWeightMethod(g_SeqWeight2.get());
	SetMuscleTree(GuideTree);

	if (g_bAnchors.get())
		RefineVert(msa, GuideTree, g_uMaxIters.get() - 2);
	else
		RefineHoriz(msa, GuideTree, g_uMaxIters.get() - 2, false, false);

#if	0
// Refining by subfamilies is disabled as it didn't give better
// results. I tried doing this before and after RefineHoriz.
// Should get back to this as it seems like this should work.
	RefineSubfams(msa, GuideTree, g_uMaxIters.get() - 2);
#endif

	ValidateMuscleIds(msa);
	ValidateMuscleIds(GuideTree);

	//TextFile fileOut(g_pstrOutFileName.get(), true);
	//MHackEnd(msa);
	//msa.ToFile(fileOut);
	MuscleOutput(msa);
	}

void Run()
	{
	SetStartTime();
	Log("Started %s\n", GetTimeAsStr());
	for (int i = 0; i < g_argc.get(); ++i)
		Log("%s ", g_argv.get()[i]);
	Log("\n");

#if	TIMING
	TICKS t1 = GetClockTicks();
#endif
	if (g_bRefine.get())
		Refine();
	else if (g_bRefineW.get())
		{
		extern void DoRefineW();
		DoRefineW();
		}
	else if (g_bProfDB.get())
		ProfDB();
	else if (g_bSW.get())
		Local();
	else if (0 != g_pstrSPFileName.get())
		DoSP();
	else if (g_bProfile.get())
		Profile();
	else if (g_bPPScore.get())
		PPScore();
	else if (g_bPAS.get())
		ProgAlignSubFams();
	else if (g_bMakeTree.get())
		{
		extern void DoMakeTree();
		DoMakeTree();
		}
	else
		DoMuscle();

#if	TIMING
	extern TICKS g_ticksDP;
	extern TICKS g_ticksObjScore;
	TICKS t2 = GetClockTicks();
	TICKS TotalTicks = t2 - t1;
	TICKS ticksOther = TotalTicks - g_ticksDP - g_ticksObjScore;
	double dSecs = TicksToSecs(TotalTicks);
	double PctDP = (double) g_ticksDP*100.0/(double) TotalTicks;
	double PctOS = (double) g_ticksObjScore*100.0/(double) TotalTicks;
	double PctOther = (double) ticksOther*100.0/(double) TotalTicks;
	Log("                 Ticks     Secs    Pct\n");
	Log("          ============  =======  =====\n");
	Log("DP        %12ld  %7.2f  %5.1f%%\n",
	  (long) g_ticksDP, TicksToSecs(g_ticksDP), PctDP);
	Log("OS        %12ld  %7.2f  %5.1f%%\n",
	  (long) g_ticksObjScore, TicksToSecs(g_ticksObjScore), PctOS);
	Log("Other     %12ld  %7.2f  %5.1f%%\n",
	  (long) ticksOther, TicksToSecs(ticksOther), PctOther);
	Log("Total     %12ld  %7.2f  100.0%%\n", (long) TotalTicks, dSecs);
#endif

	ListDiagSavings();
	Log("Finished %s\n", GetTimeAsStr());
	}
} 

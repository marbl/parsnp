#include "libMUSCLE/muscle.h"
#include "libMUSCLE/textfile.h"
#include "libMUSCLE/objscore.h"
#include "libMUSCLE/profile.h"
#include "libMUSCLE/enumopts.h"

namespace muscle {

const double DEFAULT_MAX_MB_FRACT = 0.8;

TLS<SCORE> g_scoreCenter(0);
TLS<SCORE> g_scoreGapExtend(0);
TLS<SCORE> g_scoreGapOpen2(MINUS_INFINITY);
TLS<SCORE> g_scoreGapExtend2(MINUS_INFINITY);
TLS<SCORE> g_scoreGapAmbig(0);
TLS<SCORE> g_scoreAmbigFactor(0);

extern SCOREMATRIX VTML_LA;
extern SCOREMATRIX PAM200;
extern SCOREMATRIX PAM200NoCenter;
extern SCOREMATRIX VTML_SP;
extern SCOREMATRIX VTML_SPNoCenter;
extern SCOREMATRIX NUC_SP;

TLS<PTR_SCOREMATRIX> g_ptrScoreMatrix;

TLS<const char *> g_pstrInFileName("-");
TLS<const char *> g_pstrOutFileName("-");
TLS<const char *> g_pstrFASTAOutFileName(0);
TLS<const char *> g_pstrMSFOutFileName(0);
TLS<const char *> g_pstrClwOutFileName(0);
TLS<const char *> g_pstrClwStrictOutFileName(0);
TLS<const char *> g_pstrHTMLOutFileName(0);
TLS<const char *> g_pstrPHYIOutFileName(0);
TLS<const char *> g_pstrPHYSOutFileName(0);
TLS<const char *> g_pstrDistMxFileName1(0);
TLS<const char *> g_pstrDistMxFileName2(0);

TLS<const char *> g_pstrFileName1(0);
TLS<const char *> g_pstrFileName2(0);

TLS<const char *> g_pstrSPFileName(0);
TLS<const char *> g_pstrMatrixFileName(0);

TLS<const char *> g_pstrUseTreeFileName(0);
TLS<bool> g_bUseTreeNoWarn(false);

TLS<const char *> g_pstrComputeWeightsFileName(0);
TLS<const char *> g_pstrScoreFileName(0);	// AED: init these to null to avoid crashes on optimized code

TLS<const char *> g_pstrProf1FileName(0);
TLS<const char *> g_pstrProf2FileName(0);

TLS<unsigned> g_uSmoothWindowLength(7);
TLS<unsigned> g_uAnchorSpacing(32);
TLS<unsigned> g_uMaxTreeRefineIters(1);

TLS<unsigned> g_uRefineWindow(200);
TLS<unsigned> g_uWindowFrom(0);
TLS<unsigned> g_uWindowTo(0);
TLS<unsigned> g_uSaveWindow(uInsane);
TLS<unsigned> g_uWindowOffset(0);

TLS<unsigned> g_uMaxSubFamCount(5);

TLS<unsigned> g_uHydrophobicRunLength(5);
TLS<float> g_dHydroFactor((float) 1.2);

TLS<unsigned> g_uMinDiagLength(24);	// TODO alpha -- should depend on alphabet?
TLS<unsigned> g_uMaxDiagBreak(1);
TLS<unsigned> g_uDiagMargin(5);

TLS<float> g_dSUEFF((float)0.1);

TLS<bool> g_bPrecompiledCenter(true);
TLS<bool> g_bNormalizeCounts(false);
TLS<bool> g_bDiags1(false);
TLS<bool> g_bDiags2(false);
TLS<bool> g_bAnchors(true);
TLS<bool> g_bQuiet(false);
TLS<bool> g_bVerbose(false);
TLS<bool> g_bRefine(false);
TLS<bool> g_bRefineW(false);
TLS<bool> g_bProfDB(false);
TLS<bool> g_bLow(false);
TLS<bool> g_bSW(false);
TLS<bool> g_bClusterOnly(false);
TLS<bool> g_bProfile(false);
TLS<bool> g_bProfileOnStdIn(false);
TLS<bool> g_bAnchoredPP(false);
TLS<bool> g_bPPScore(false);
TLS<bool> g_bBrenner(false);
TLS<bool> g_bDimer(false);
TLS<bool> g_bVersion(false);
TLS<bool> g_bStable(false);
TLS<bool> g_bFASTA(false);
TLS<bool> g_bPAS(false);
TLS<bool> g_bTomHydro(false);
TLS<bool> g_bMakeTree(false);

#if	DEBUG
TLS<bool> g_bCatchExceptions(false);
#else
TLS<bool> g_bCatchExceptions(true);
#endif

TLS<bool> g_bMSF(false);
TLS<bool> g_bAln(false);
TLS<bool> g_bClwStrict(false);
TLS<bool> g_bHTML(false);
TLS<bool> g_bPHYI(false);
TLS<bool> g_bPHYS(false);
//TLS<int> g_scoreGapOpen(1);
//TLS<int> g_scoreGapExtend(0);
TLS<unsigned> g_uMaxIters(8);
TLS<unsigned long> g_ulMaxSecs(0);
TLS<unsigned> g_uMaxMB(16000);

TLS<PPSCORE> g_PPScore(PPSCORE_LE);
TLS<OBJSCORE> g_ObjScore(OBJSCORE_SPM);

TLS<SEQWEIGHT> g_SeqWeight1(SEQWEIGHT_ClustalW);
TLS<SEQWEIGHT> g_SeqWeight2(SEQWEIGHT_ClustalW);

TLS<DISTANCE> g_Distance1(DISTANCE_Kmer6_6);
TLS<DISTANCE> g_Distance2(DISTANCE_PctIdKimura);

TLS<CLUSTER> g_Cluster1(CLUSTER_UPGMB);
TLS<CLUSTER> g_Cluster2(CLUSTER_UPGMB);

TLS<ROOT> g_Root1(ROOT_Pseudo);
TLS<ROOT> g_Root2(ROOT_Pseudo);

TLS<bool> g_bDiags;

TLS<SEQTYPE> g_SeqType(SEQTYPE_Auto);

TLS<TERMGAPS> g_TermGaps(TERMGAPS_Half);

//------------------------------------------------------
// These parameters depending on the chosen prof-prof
// score (g_PPScore), initialized to "Undefined".
TLS<float> g_dSmoothScoreCeil(fInsane);
TLS<float> g_dMinBestColScore(fInsane);
TLS<float> g_dMinSmoothScore(fInsane);
TLS<SCORE> g_scoreGapOpen(fInsane);
//------------------------------------------------------

static unsigned atou(const char *s)
	{
	return (unsigned) atoi(s);
	}

const char *MaxSecsToStr()
	{
	if (0 == g_ulMaxSecs.get())
		return "(No limit)";
	return SecsToStr(g_ulMaxSecs.get());
	}

void ListParams()
	{
	Log("\n");
	Log("%s\n", MUSCLE_LONG_VERSION);
	Log("http://www.drive5.com/muscle\n");
	Log("\n");
	Log("Profile-profile score    %s\n", PPSCOREToStr(g_PPScore.get()));
	Log("Max iterations           %u\n", g_uMaxIters.get());
	Log("Max trees                %u\n", g_uMaxTreeRefineIters.get());
	Log("Max time                 %s\n", MaxSecsToStr());
	Log("Max MB                   %u\n", g_uMaxMB.get());
	Log("Gap open                 %g\n", g_scoreGapOpen.get());
	Log("Gap extend (dimer)       %g\n", g_scoreGapExtend.get());
	Log("Gap ambig factor         %g\n", g_scoreAmbigFactor.get());
	Log("Gap ambig penalty        %g\n", g_scoreGapAmbig.get());
	Log("Center (LE)              %g\n", g_scoreCenter.get());
	Log("Term gaps                %s\n", TERMGAPSToStr(g_TermGaps.get()));

	Log("Smooth window length     %u\n", g_uSmoothWindowLength.get());
	Log("Refine window length     %u\n", g_uRefineWindow.get());
	Log("Min anchor spacing       %u\n", g_uAnchorSpacing.get());
	Log("Min diag length (lambda) %u\n", g_uMinDiagLength.get());
	Log("Diag margin (mu)         %u\n", g_uDiagMargin.get());
	Log("Min diag break           %u\n", g_uMaxDiagBreak.get());
	Log("Hydrophobic window       %u\n", g_uHydrophobicRunLength.get());

	Log("Hydrophobic gap factor   %g\n", g_dHydroFactor.get());
	Log("Smooth score ceiling     %g\n", g_dSmoothScoreCeil.get());
	Log("Min best col score       %g\n", g_dMinBestColScore.get());
	Log("Min anchor score         %g\n", g_dMinSmoothScore.get());
	Log("SUEFF                    %g\n", g_dSUEFF.get());

	Log("Brenner root MSA         %s\n", BoolToStr(g_bBrenner.get()));
	Log("Normalize counts         %s\n", BoolToStr(g_bNormalizeCounts.get()));
	Log("Diagonals (1)            %s\n", BoolToStr(g_bDiags1.get()));
	Log("Diagonals (2)            %s\n", BoolToStr(g_bDiags2.get()));
	Log("Anchors                  %s\n", BoolToStr(g_bAnchors.get()));
	Log("MSF output format        %s\n", BoolToStr(g_bMSF.get()));
	Log("Phylip interleaved       %s\n", BoolToStr(g_bPHYI.get()));
	Log("Phylip sequential        %s\n", BoolToStr(g_bPHYS.get()));
	Log("ClustalW output format   %s\n", BoolToStr(g_bAln.get()));
	Log("Catch exceptions         %s\n", BoolToStr(g_bCatchExceptions.get()));
	Log("Quiet                    %s\n", BoolToStr(g_bQuiet.get()));
	Log("Refine                   %s\n", BoolToStr(g_bRefine.get()));
	Log("ProdfDB                  %s\n", BoolToStr(g_bProfDB.get()));
	Log("Low complexity profiles  %s\n", BoolToStr(g_bLow.get()));

	Log("Objective score          %s\n", OBJSCOREToStr(g_ObjScore.get()));

	Log("Distance method (1)      %s\n", DISTANCEToStr(g_Distance1.get()));
	Log("Clustering method (1)    %s\n", CLUSTERToStr(g_Cluster1.get()));
	Log("Root method (1)          %s\n", ROOTToStr(g_Root1.get()));
	Log("Sequence weighting (1)   %s\n", SEQWEIGHTToStr(g_SeqWeight1.get()));

	Log("Distance method (2)      %s\n", DISTANCEToStr(g_Distance2.get()));
	Log("Clustering method (2)    %s\n", CLUSTERToStr(g_Cluster2.get()));
	Log("Root method (2)          %s\n", ROOTToStr(g_Root2.get()));
	Log("Sequence weighting (2)   %s\n", SEQWEIGHTToStr(g_SeqWeight2.get()));

	Log("\n");
	}

static void SetDefaultsLE()
	{
	g_ptrScoreMatrix.get() = &VTML_LA;

	//g_scoreGapOpen.get() = (SCORE) -3.00;
	//g_scoreCenter.get() = (SCORE) -0.55;
	g_scoreGapOpen.get() = (SCORE) -2.9;
	g_scoreCenter.get() = (SCORE) -0.52;

	g_bNormalizeCounts.get() = true;

	//g_dSmoothScoreCeil.get() = 5.0;
	//g_dMinBestColScore.get() = 4.0;
	//g_dMinSmoothScore.get() = 2.0;
	g_dSmoothScoreCeil.get() = 3.0;
	g_dMinBestColScore.get() = 2.0;
	g_dMinSmoothScore.get() = 1.0;

	g_Distance1.get() = DISTANCE_Kmer6_6;
	g_Distance2.get() = DISTANCE_PctIdKimura;
	}

static void SetDefaultsSP()
	{
	g_ptrScoreMatrix.get() = &PAM200;

	g_scoreGapOpen.get() = -1439;
	g_scoreCenter.get() = 0.0;	// center pre-added into score mx

	g_bNormalizeCounts.get() = false;

	g_dSmoothScoreCeil.get() = 200.0;
	g_dMinBestColScore.get() = 300.0;
	g_dMinSmoothScore.get() = 125.0;

	g_Distance1.get() = DISTANCE_Kmer6_6;
	g_Distance2.get() = DISTANCE_PctIdKimura;
	}

static void SetDefaultsSV()
	{
	g_ptrScoreMatrix.get() = &VTML_SP;

	g_scoreGapOpen.get() = -300;
	g_scoreCenter.get() = 0.0;	// center pre-added into score mx

	g_bNormalizeCounts.get() = false;

	g_dSmoothScoreCeil.get() = 90.0;
	g_dMinBestColScore.get() = 130.0;
	g_dMinSmoothScore.get() = 40.0;

	g_Distance1.get() = DISTANCE_Kmer6_6;
	g_Distance2.get() = DISTANCE_PctIdKimura;
	}

//static void SetDefaultsSPN()
//	{
//	g_ptrScoreMatrix.get() = &NUC_SP;
//
//	g_scoreGapOpen.get() = -400;
//	g_scoreCenter.get() = 0.0;	// center pre-added into score mx
//
//	g_bNormalizeCounts.get() = false;
//
//	g_dSmoothScoreCeil.get() = 999.0;		// disable
//	g_dMinBestColScore.get() = 90;
//	g_dMinSmoothScore.get() = 90;
//
//	g_Distance1.get() = DISTANCE_Kmer4_6;
//	g_Distance2.get() = DISTANCE_PctIdKimura;
//	}

static void SetDefaultsSPN_DNA()
	{
	g_ptrScoreMatrix.get() = &NUC_SP;

	if ( g_scoreGapOpen.get() == fInsane )
		g_scoreGapOpen.get() = -400;
	g_scoreCenter.get() = 0.0;	// center pre-added into score mx
	if ( g_scoreGapExtend.get() == 0 )
		g_scoreGapExtend.get() = 0.0;

	g_bNormalizeCounts.get() = false;

	g_dSmoothScoreCeil.get() = 999.0;		// disable
	g_dMinBestColScore.get() = 90;
	g_dMinSmoothScore.get() = 90;

	g_Distance1.get() = DISTANCE_Kmer4_6;
	g_Distance2.get() = DISTANCE_PctIdKimura;
	}

static void SetDefaultsSPN_RNA()
	{
	g_ptrScoreMatrix.get() = &NUC_SP;

	g_scoreGapOpen.get() = -420;
	g_scoreCenter.get() = -300;	// total center = NUC_EXTEND - 300 
	g_scoreGapExtend.get() = 0.0;

	g_bNormalizeCounts.get() = false;

	g_dSmoothScoreCeil.get() = 999.0;		// disable
	g_dMinBestColScore.get() = 90;
	g_dMinSmoothScore.get() = 90;

	g_Distance1.get() = DISTANCE_Kmer4_6;
	g_Distance2.get() = DISTANCE_PctIdKimura;
	}

static void FlagParam(const char *OptName, bool *ptrParam, bool bValueIfFlagSet)
	{
	bool bIsSet = FlagOpt(OptName);
	if (bIsSet)
		*ptrParam = bValueIfFlagSet;
	}

static void StrParam(const char *OptName, const char **ptrptrParam)
	{
	const char *opt = ValueOpt(OptName);
	if (0 != opt)
		*ptrptrParam = opt;
	}

static void FloatParam(const char *OptName, float *ptrParam)
	{
	const char *opt = ValueOpt(OptName);
	if (0 != opt)
		*ptrParam = (float) atof(opt);
	}

static void UintParam(const char *OptName, unsigned *ptrParam)
	{
	const char *opt = ValueOpt(OptName);
	if (0 != opt)
		*ptrParam = atou(opt);
	}

static void EnumParam(const char *OptName, EnumOpt *Opts, int *Param)
	{
	const char *Value = ValueOpt(OptName);
	if (0 == Value)
		return;

	for (;;)
		{
		if (0 == Opts->pstrOpt)
			Quit("Invalid parameter -%s %s", OptName, Value);
		if (0 == stricmp(Value, Opts->pstrOpt))
			{
			*Param = Opts->iValue;
			return;
			}
		++Opts;
		}
	}

static void SetPPDefaultParams()
	{
	switch (g_PPScore.get())
		{
	case PPSCORE_SP:
		SetDefaultsSP();
		break;

	case PPSCORE_LE:
		SetDefaultsLE();
		break;

	case PPSCORE_SV:
		SetDefaultsSV();
		break;

	case PPSCORE_SPN:
		switch (g_Alpha.get())
			{
		case ALPHA_DNA:
			SetDefaultsSPN_DNA();
			break;
		case ALPHA_RNA:
			SetDefaultsSPN_RNA();
			break;
		default:
			Quit("Invalid alpha %d", g_Alpha.get());
			}
		break;

	default:
		Quit("Invalid g_PPScore.get()");
		}
	}

static void SetPPCommandLineParams()
	{
	FloatParam("GapOpen", &g_scoreGapOpen.get());
	FloatParam("GapOpen2", &g_scoreGapOpen2.get());
	FloatParam("GapExtend", &g_scoreGapExtend.get());
	FloatParam("GapExtend2", &g_scoreGapExtend2.get());
	FloatParam("GapAmbig", &g_scoreAmbigFactor.get());
	FloatParam("Center", &g_scoreCenter.get());
	FloatParam("SmoothScoreCeil", &g_dSmoothScoreCeil.get());
	FloatParam("MinBestColScore", &g_dMinBestColScore.get());
	FloatParam("MinSmoothScore", &g_dMinSmoothScore.get());

	EnumParam("Distance", DISTANCE_Opts, (int *) &g_Distance1.get());
	EnumParam("Distance", DISTANCE_Opts, (int *) &g_Distance2.get());

	EnumParam("Distance1", DISTANCE_Opts, (int *) &g_Distance1.get());
	EnumParam("Distance2", DISTANCE_Opts, (int *) &g_Distance2.get());
	}

void SetPPScore(bool bRespectFlagOpts)
	{
	if (bRespectFlagOpts)
		{
		if (FlagOpt("SP"))
			g_PPScore.get() = PPSCORE_SP;
		else if (FlagOpt("LE"))
			g_PPScore.get() = PPSCORE_LE;
		else if (FlagOpt("SV"))
			g_PPScore.get() = PPSCORE_SV;
		else if (FlagOpt("SPN"))
			g_PPScore.get() = PPSCORE_SPN;
		}

	switch (g_PPScore.get())
		{
	case PPSCORE_LE:
	case PPSCORE_SP:
	case PPSCORE_SV:
		if (ALPHA_RNA == g_Alpha.get() || ALPHA_DNA == g_Alpha.get())
			g_PPScore.get() = PPSCORE_SPN;
		break;
	case PPSCORE_SPN:
		if (ALPHA_Amino == g_Alpha.get())
			g_PPScore.get() = PPSCORE_LE;
		break;
		}

	SetPPDefaultParams();
	SetPPCommandLineParams();

	DYN_PTR_SCOREMATRIX UserMatrix = 0;
	if (0 != g_pstrMatrixFileName.get())
		{
		const char *FileName = g_pstrMatrixFileName.get();
		const char *Path = getenv("MUSCLE_MXPATH");
		if (Path != 0)
			{
			size_t n = strlen(Path) + 1 + strlen(FileName) + 1;
			char *NewFileName = new char[n];
			sprintf(NewFileName, "%s/%s", Path, FileName);
			FileName = NewFileName;
			}
		TextFile File(FileName);
		UserMatrix = ReadMx(File);
// AED 21/12/2006: allow a custom nucleotide substitution matrix (don't force AA alignment)
//		g_Alpha.get() = ALPHA_Amino;
//		g_PPScore.get() = PPSCORE_SP;
		}
	if (0 != UserMatrix)
		g_ptrScoreMatrix.get() = UserMatrix;

// AED 21/12/2006: if a nucleotide matrix was loaded, add the gap extend penalty directly
// to the matrix
	if( 0 != UserMatrix && g_PPScore.get() == PPSCORE_SPN )
	{
		// default gap extend is 30, so add 60
		// command-line gap extend will be negative, so multiply by -2
		float add_score = g_scoreGapExtend.get() == 0 ? 60 : -2 * g_scoreGapExtend.get();
		for( int i = 0; i < 4; ++i )
			for( int j = 0; j < 4; ++j )
				(*UserMatrix)[i][j] += add_score;
	}

	if (g_bVerbose.get())
		ListParams();
	}

void SetPPScore(PPSCORE p)
	{
	g_PPScore.get() = p;
	SetPPScore(true);
	}

static void SetMaxSecs()
	{
	float fMaxHours = 0.0;
	FloatParam("MaxHours", &fMaxHours);
	if (0.0 == fMaxHours)
		return;
	g_ulMaxSecs.get() = (unsigned long) (fMaxHours*60*60);
	}

static bool CanDoLowComplexity()
	{
	if (g_SeqWeight1.get() != SEQWEIGHT_ClustalW)
		return false;
	if (1 == g_uMaxIters.get())
		return true;
	return g_SeqWeight2.get() == SEQWEIGHT_ClustalW;
	}

bool MissingCommand()
	{
	if (strcmp(g_pstrInFileName.get(), "-"))
		return false;
	if (0 != g_pstrFileName1.get())
		return false;
	if (0 != g_pstrSPFileName.get())
		return false;
	return true;
	}

void SetParams()
	{
	SetMaxSecs();

	StrParam("in", &g_pstrInFileName.get());
	StrParam("out", &g_pstrOutFileName.get());

	StrParam("FASTAOut", &g_pstrFASTAOutFileName.get());
	StrParam("ClwOut", &g_pstrClwOutFileName.get());
	StrParam("ClwStrictOut", &g_pstrClwStrictOutFileName.get());
	StrParam("HTMLOut", &g_pstrHTMLOutFileName.get());
	StrParam("PHYIOut", &g_pstrPHYIOutFileName.get());
	StrParam("PHYSOut", &g_pstrPHYSOutFileName.get());
	StrParam("MSFOut", &g_pstrMSFOutFileName.get());

	StrParam("in1", &g_pstrFileName1.get());
	StrParam("in2", &g_pstrFileName2.get());

	StrParam("Matrix", &g_pstrMatrixFileName.get());
	StrParam("SPScore", &g_pstrSPFileName.get());

	StrParam("UseTree_NoWarn", &g_pstrUseTreeFileName.get());
	if (0 != g_pstrUseTreeFileName.get())
		g_bUseTreeNoWarn.get() = true;

	StrParam("UseTree", &g_pstrUseTreeFileName.get());
	StrParam("ComputeWeights", &g_pstrComputeWeightsFileName.get());
	StrParam("ScoreFile", &g_pstrScoreFileName.get());
	StrParam("DistMx1", &g_pstrDistMxFileName1.get());
	StrParam("DistMx2", &g_pstrDistMxFileName2.get());

	FlagParam("Core", &g_bCatchExceptions.get(), false);
	FlagParam("NoCore", &g_bCatchExceptions.get(), true);

	FlagParam("Diags1", &g_bDiags1.get(), true);
	FlagParam("Diags2", &g_bDiags2.get(), true);

	bool Diags = false;
	FlagParam("Diags", &Diags, true);
	if (Diags)
		{
		g_bDiags1.get() = true;
		g_bDiags2.get() = true;
		}

	FlagParam("Anchors", &g_bAnchors.get(), true);
	FlagParam("NoAnchors", &g_bAnchors.get(), false);

	FlagParam("Quiet", &g_bQuiet.get(), true);
	FlagParam("Verbose", &g_bVerbose.get(), true);
	FlagParam("Version", &g_bVersion.get(), true);
	FlagParam("Stable", &g_bStable.get(), true);
	FlagParam("Group", &g_bStable.get(), false);
	FlagParam("Refine", &g_bRefine.get(), true);
	FlagParam("RefineW", &g_bRefineW.get(), true);
	FlagParam("ProfDB", &g_bProfDB.get(), true);
	FlagParam("SW", &g_bSW.get(), true);
	FlagParam("ClusterOnly", &g_bClusterOnly.get(), true);
	FlagParam("Profile", &g_bProfile.get(), true);
	FlagParam("ProfileOnStdIn", &g_bProfileOnStdIn.get(), true);
	FlagParam("AnchoredPP", &g_bAnchoredPP.get(), true);
	FlagParam("PPScore", &g_bPPScore.get(), true);
	FlagParam("Brenner", &g_bBrenner.get(), true);
	FlagParam("Dimer", &g_bDimer.get(), true);

	FlagParam("MSF", &g_bMSF.get(), true);
	FlagParam("PHYI", &g_bPHYI.get(), true);
	FlagParam("PHYS", &g_bPHYS.get(), true);
	FlagParam("clw", &g_bAln.get(), true);
	FlagParam("HTML", &g_bHTML.get(), true);
	FlagParam("FASTA", &g_bFASTA.get(), true);
	FlagParam("PAS", &g_bPAS.get(), true);
	FlagParam("MakeTree", &g_bMakeTree.get(), true);

	bool b = false;
	FlagParam("clwstrict", &b, true);
	if (b)
		{
		g_bAln.get() = true;
		g_bClwStrict.get() = true;
		}

	UintParam("MaxIters", &g_uMaxIters.get());
	UintParam("MaxTrees", &g_uMaxTreeRefineIters.get());
	UintParam("SmoothWindow", &g_uSmoothWindowLength.get());
	UintParam("RefineWindow", &g_uRefineWindow.get());
	UintParam("FromWindow", &g_uWindowFrom.get());
	UintParam("ToWindow", &g_uWindowTo.get());
	UintParam("SaveWindow", &g_uSaveWindow.get());
	UintParam("WindowOffset", &g_uWindowOffset.get());
	UintParam("AnchorSpacing", &g_uAnchorSpacing.get());
	UintParam("DiagLength", &g_uMinDiagLength.get());
	UintParam("DiagMargin", &g_uDiagMargin.get());
	UintParam("DiagBreak", &g_uMaxDiagBreak.get());
	UintParam("MaxSubFam", &g_uMaxSubFamCount.get());

	UintParam("Hydro", &g_uHydrophobicRunLength.get());
	FlagParam("TomHydro", &g_bTomHydro.get(), true);
	if (g_bTomHydro.get())
		g_uHydrophobicRunLength.get() = 0;	

	FloatParam("SUEFF", &g_dSUEFF.get());
	FloatParam("HydroFactor", &g_dHydroFactor.get());

	EnumParam("ObjScore", OBJSCORE_Opts, (int *) &g_ObjScore.get());
	EnumParam("TermGaps", TERMGAPS_Opts, (int *) &g_TermGaps.get());

	EnumParam("Weight", SEQWEIGHT_Opts, (int *) &g_SeqWeight1.get());
	EnumParam("Weight", SEQWEIGHT_Opts, (int *) &g_SeqWeight2.get());
	
	EnumParam("Weight1", SEQWEIGHT_Opts, (int *) &g_SeqWeight1.get());
	EnumParam("Weight2", SEQWEIGHT_Opts, (int *) &g_SeqWeight2.get());

	EnumParam("Cluster", CLUSTER_Opts, (int *) &g_Cluster1.get());
	EnumParam("Cluster", CLUSTER_Opts, (int *) &g_Cluster2.get());

	EnumParam("Cluster1", CLUSTER_Opts, (int *) &g_Cluster1.get());
	EnumParam("Cluster2", CLUSTER_Opts, (int *) &g_Cluster2.get());

	EnumParam("Root1", ROOT_Opts, (int *) &g_Root1.get());
	EnumParam("Root2", ROOT_Opts, (int *) &g_Root2.get());

	EnumParam("SeqType", SEQTYPE_Opts, (int *) &g_SeqType.get());

	g_scoreGapAmbig.get() = g_scoreGapOpen.get()*g_scoreAmbigFactor.get();
	g_bLow.get() = CanDoLowComplexity();

	if (g_bDimer.get())
		g_bPrecompiledCenter.get() = false;

	UintParam("MaxMB", &g_uMaxMB.get());
	if (0 == ValueOpt("MaxMB"))
		g_uMaxMB.get() = (unsigned) (GetRAMSizeMB()*DEFAULT_MAX_MB_FRACT);
	}
} 

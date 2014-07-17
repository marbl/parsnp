#include "libMUSCLE/muscle.h"
#include "libMUSCLE/textfile.h"
#include "libMUSCLE/msa.h"
#include "libMUSCLE/tree.h"
#include "libMUSCLE/profile.h"
#include "libMUSCLE/objscore.h"

namespace muscle {

bool TreeNeededForWeighting(SEQWEIGHT s)
	{
	switch (s)
		{
	case SEQWEIGHT_ClustalW:
	case SEQWEIGHT_ThreeWay:
		return true;
	default:
		return false;
		}
	}

static ProfPos *ProfileFromMSALocal(MSA &msa, Tree &tree)
	{
	const unsigned uSeqCount = msa.GetSeqCount();
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		msa.SetSeqId(uSeqIndex, uSeqIndex);
	if (TreeNeededForWeighting(g_SeqWeight2.get()))
		{
		TreeFromMSA(msa, tree, g_Cluster2.get(), g_Distance2.get(), g_Root1.get());
		SetMuscleTree(tree);
		}
	return ProfileFromMSA(msa);
	}

void SetProfileProfileAlphabet(MSA &msa1, MSA &msa2)
{
	ALPHA Alpha = ALPHA_Undefined;
	switch (g_SeqType.get())
		{
	case SEQTYPE_Auto:
		Alpha = msa1.GuessAlpha();
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

	msa1.FixAlpha();
	msa2.FixAlpha();

	if (ALPHA_DNA == Alpha || ALPHA_RNA == Alpha)
		SetPPScore(PPSCORE_SPN);
}

void ProfileProfile(MSA &msa1, MSA &msa2, MSA &msaOut)
	{

	unsigned uLength1;
	unsigned uLength2;

	uLength1 = msa1.GetColCount();
	uLength2 = msa2.GetColCount();

	Tree tree1;
	Tree tree2;
	ProfPos *Prof1 = ProfileFromMSALocal(msa1, tree1);
	ProfPos *Prof2 = ProfileFromMSALocal(msa2, tree2);

	PWPath Path;
	ProfPos *ProfOut;
	unsigned uLengthOut;
	Progress("Aligning profiles");
	AlignTwoProfs(Prof1, uLength1, 1.0, Prof2, uLength2, 1.0, Path, &ProfOut, &uLengthOut);

	Progress("Building output");
	AlignTwoMSAsGivenPath(Path, msa1, msa2, msaOut);
	delete[] Prof1;
	delete[] Prof2;
	delete[] ProfOut;
	}

// Do profile-profile alignment
void Profile()
	{
	if ( !g_bProfileOnStdIn.get() && (0 == g_pstrFileName1.get() || 0 == g_pstrFileName2.get()))
		Quit("-profile needs -in1 and -in2 or -ProfileOnStdIn");

	SetSeqWeightMethod(g_SeqWeight1.get());

	MSA msa1;
	MSA msa2;
	MSA msaOut;

	if( !g_bProfileOnStdIn.get() )
	{
		TextFile file1(g_pstrFileName1.get());
		TextFile file2(g_pstrFileName2.get());
		msa1.FromFile(file1);
		msa2.FromFile(file2);
	}else{
		TextFile file1("-");
		TextFile file2("-");
		msa1.FromFile(file1);
		msa2.FromFile(file2);
	}

	ALPHA Alpha = ALPHA_Undefined;
	switch (g_SeqType.get())
		{
	case SEQTYPE_Auto:
		Alpha = msa1.GuessAlpha();
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
	msa1.FixAlpha();
	msa2.FixAlpha();
	SetPPScore();

	const unsigned uSeqCount1 = msa1.GetSeqCount();
	const unsigned uSeqCount2 = msa2.GetSeqCount();
	//const unsigned uMaxSeqCount = (uSeqCount1 > uSeqCount2 ? uSeqCount1 : uSeqCount2);
	//MSA::SetIdCount(uMaxSeqCount);
	const unsigned uSumSeqCount = uSeqCount1 + uSeqCount2;
	MSA::SetIdCount(uSumSeqCount);


	SetProfileProfileAlphabet(msa1, msa2);
	if( g_bAnchoredPP.get() )
		AnchoredProfileProfile(msa1, msa2, msaOut);
	else
		ProfileProfile(msa1, msa2, msaOut);

	MuscleOutput(msaOut);
	}
} 

#include "libMUSCLE/muscle.h"
#include "libMUSCLE/textfile.h"
#include "libMUSCLE/msa.h"
#include "libMUSCLE/objscore.h"
#include "libMUSCLE/tree.h"
#include "libMUSCLE/profile.h"

namespace muscle {

void DoSP()
	{
	TextFile f(g_pstrSPFileName.get());

	MSA a;
	a.FromFile(f);

	ALPHA Alpha = ALPHA_Undefined;
	switch (g_SeqType.get())
		{
	case SEQTYPE_Auto:
		Alpha = a.GuessAlpha();
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
	a.FixAlpha();

	SetPPScore();

	const unsigned uSeqCount = a.GetSeqCount();
	if (0 == uSeqCount)
		Quit("No sequences in input file %s", g_pstrSPFileName.get());

	MSA::SetIdCount(uSeqCount);
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		a.SetSeqId(uSeqIndex, uSeqIndex);

	SetSeqWeightMethod(g_SeqWeight1.get());
	Tree tree;
	TreeFromMSA(a, tree, g_Cluster2.get(), g_Distance2.get(), g_Root2.get());
	SetMuscleTree(tree);
	SetMSAWeightsMuscle((MSA &) a);

	SCORE SP = ObjScoreSP(a);

	Log("File=%s;SP=%.4g\n", g_pstrSPFileName.get(), SP);
	fprintf(stderr, "File=%s;SP=%.4g\n", g_pstrSPFileName.get(), SP);
	}
} 

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

void ProfDB()
	{
	SetOutputFileName(g_pstrOutFileName.get());
	SetInputFileName(g_pstrFileName2.get());
	SetStartTime();

	TextFile file1(g_pstrFileName1.get());
	TextFile file2(g_pstrFileName2.get());

	SetMaxIters(g_uMaxIters.get());
	SetSeqWeightMethod(g_SeqWeight1.get());

	TextFile fileIn(g_pstrFileName1.get());
	MSA msa1;
	msa1.FromFile(fileIn);

	const unsigned uSeqCount1 = msa1.GetSeqCount();
	if (0 == uSeqCount1)
		Quit("No sequences in input alignment");

	SeqVect v;
	v.FromFASTAFile(file2);
	const unsigned uSeqCount2 = v.Length();
	if (0 == uSeqCount2)
		Quit("No sequences in input alignment");

	MSA::SetIdCount(uSeqCount1 + uSeqCount2);
	SetProgressDesc("Align sequence database to profile");
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount2; ++uSeqIndex)
		{
		Progress(uSeqIndex, uSeqCount2);
		Seq &s = *(v[uSeqIndex]);
		s.SetId(0);
		MSA msaTmp;
		msaTmp.FromSeq(s);
		MSA msaOut;
		SetProfileProfileAlphabet(msa1, msaTmp);
		ProfileProfile(msa1, msaTmp, msaOut);
		msa1.Copy(msaOut);
		}
	ProgressStepsDone();

	TextFile fileOut(g_pstrOutFileName.get(), true);
	msa1.ToFile(fileOut);
	}
} 

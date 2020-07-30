//Fork of libMUSCLE code from http://sourceforge.net/p/mauve/code/HEAD/tree/muscle/

#include "MuscleInterface.h"
#include "libMUSCLE/muscle.h"
#include "libMUSCLE/params.h"
#include "libMUSCLE/seq.h"
#include "libMUSCLE/seqvect.h"
#include "libMUSCLE/msa.h"
#include "libMUSCLE/tree.h"
#include "libMUSCLE/clust.h"
#include "libMUSCLE/profile.h"
#include "libMUSCLE/distfunc.h"
#include "libMUSCLE/clustsetdf.h"
#include "libMUSCLE/textfile.h"
#include "libMUSCLE/types.h"
#include "libMUSCLE/threadstorage.h"


namespace muscle {
extern void MUSCLE(SeqVect &v, MSA &msaOut);
//extern void RefineW(const MSA &msaIn, MSA &msaOut);
}

using namespace muscle;

MuscleInterface& MuscleInterface::getMuscleInterface()
{
	static MuscleInterface m_ci;
	return m_ci;
}

MuscleInterface::MuscleInterface()
{

}

bool MuscleInterface::CallMuscleFast ( vector< string >& aln_matrix, const vector< string >& seq2aln)
{


//#pragma omp threadprivate (g_SeqType,g_uMaxIters,g_bStable,g_bQuiet,g_SeqWeight1,g_ulMaxSecs, g_bDiags1, g_bDiags2, g_bAnchors,g_bRefine,g_bVerbose)
//	{
	g_SeqType.get() = SEQTYPE_DNA;
	g_uMaxIters.get() = 1;
	g_bStable.get() = true;
	g_bVerbose.get() = false;
	g_bQuiet.get() = true;
	g_SeqWeight1.get() = SEQWEIGHT_ClustalW;
	SetMaxIters(g_uMaxIters.get());
	SetSeqWeightMethod(g_SeqWeight1.get());
	g_ulMaxSecs.get() = 0;
    SeqVect sv;
    const char* seqname = "seq00000";
    for(size_t seqI = 0; seqI < seq2aln.size(); seqI++)
    {
        Seq curseq;
        curseq.SetId(seqI);
        curseq.SetName(seqname);
        curseq.resize(seq2aln.at(seqI).size());
        std::copy(seq2aln[seqI].begin(), seq2aln[seqI].end(), curseq.begin());
        //cout << curseq << endl;
        sv.AppendSeq(curseq);
    }

    MSA msaTmp;
    MUSCLE(sv,msaTmp);

    aln_matrix.clear();
    aln_matrix.resize(msaTmp.GetSeqCount());
    for(size_t seqI = 0; seqI < msaTmp.GetSeqCount(); seqI++)
    {
        unsigned indie = msaTmp.GetSeqIndex(seqI);
        const char* buf = msaTmp.GetSeqBuffer(indie);
        string curseq(buf,msaTmp.GetColCount());
        swap(aln_matrix[seqI],curseq);
    }
//	}
    return true;

}

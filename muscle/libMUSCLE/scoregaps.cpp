#include "libMUSCLE/muscle.h"
#include "libMUSCLE/msa.h"
#include "libMUSCLE/objscore.h"
#include "libMUSCLE/threadstorage.h"

namespace muscle {

#define TRACE	0

struct GAPINFO
	{
	GAPINFO *Next;
	unsigned Start;
	unsigned End;
	};

static TLS<GAPINFO **> g_Gaps;
static TLS<GAPINFO *> g_FreeList;
static TLS<unsigned> g_MaxSeqCount;
static TLS<unsigned> g_MaxColCount;
static TLS<unsigned> g_ColCount;
static TLS<bool *> g_ColDiff;

static GAPINFO *NewGapInfo()
	{
	if (0 == g_FreeList.get())
		{
		const int NEWCOUNT = 256;
		GAPINFO *NewList = new GAPINFO[NEWCOUNT];
		g_FreeList.get() = &NewList[0];
		for (int i = 0; i < NEWCOUNT-1; ++i)
			NewList[i].Next = &NewList[i+1];
		NewList[NEWCOUNT-1].Next = 0;
		}
	GAPINFO *GI = g_FreeList.get();
	g_FreeList.get() = g_FreeList.get()->Next;
	return GI;
	}

static void FreeGapInfo(GAPINFO *GI)
	{
	GI->Next = g_FreeList.get();
	g_FreeList.get() = GI;
	}

// TODO: This could be much faster, no need to look
// at all columns.
static void FindIntersectingGaps(const MSA &msa, unsigned SeqIndex)
	{
	const unsigned ColCount = msa.GetColCount();
	bool InGap = false;
	bool Intersects = false;
	unsigned Start = uInsane;
	for (unsigned Col = 0; Col <= ColCount; ++Col)
		{
		bool Gap = ((Col != ColCount) && msa.IsGap(SeqIndex, Col));
		if (Gap)
			{
			if (!InGap)
				{
				InGap = true;
				Start = Col;
				}
			if (g_ColDiff.get()[Col])
				Intersects = true;
			}
		else if (InGap)
			{
			InGap = false;
			if (Intersects)
				{
				GAPINFO *GI = NewGapInfo();
				GI->Start = Start;
				GI->End = Col - 1;
				GI->Next = g_Gaps.get()[SeqIndex];
				g_Gaps.get()[SeqIndex] = GI;
				}
			Intersects = false;
			}
		}
	}

static SCORE Penalty(unsigned Length, bool Term)
	{
	if (0 == Length)
		return 0;
	SCORE s1 = g_scoreGapOpen.get() + g_scoreGapExtend.get()*(Length - 1);
#if	DOUBLE_AFFINE
	SCORE s2 = g_scoreGapOpen2.get() + g_scoreGapExtend2.get()*(Length - 1);
	if (s1 > s2)
		return s1;
	return s2;
#else
	return s1;
#endif
	}

//static SCORE ScorePair(unsigned Seq1, unsigned Seq2)
//	{
//#if	TRACE
//	{
//	Log("ScorePair(%d,%d)\n", Seq1, Seq2);
//	Log("Gaps seq 1: ");
//	for (GAPINFO *GI = g_Gaps.get()[Seq1]; GI; GI = GI->Next)
//		Log(" %d-%d", GI->Start, GI->End);
//	Log("\n");
//	Log("Gaps seq 2: ");
//	for (GAPINFO *GI = g_Gaps.get()[Seq2]; GI; GI = GI->Next)
//		Log(" %d-%d", GI->Start, GI->End);
//	Log("\n");
//	}
//#endif
//	return 0;
//	}

SCORE ScoreGaps(const MSA &msa, const unsigned DiffCols[], unsigned DiffColCount)
	{
#if	TRACE
	{
	Log("ScoreGaps\n");
	Log("DiffCols ");
	for (unsigned i = 0; i < DiffColCount; ++i)
		Log(" %u", DiffCols[i]);
	Log("\n");
	Log("msa=\n");
	msa.LogMe();
	Log("\n");
	}
#endif
	const unsigned SeqCount = msa.GetSeqCount();
	const unsigned ColCount = msa.GetColCount();
	g_ColCount.get() = ColCount;

	if (SeqCount > g_MaxSeqCount.get())
		{
		delete[] g_Gaps.get();
		g_MaxSeqCount.get() = SeqCount + 256;
		g_Gaps.get() = new GAPINFO *[g_MaxSeqCount.get()];
		}
	memset(g_Gaps.get(), 0, SeqCount*sizeof(GAPINFO *));

	if (ColCount > g_MaxColCount.get())
		{
		delete[] g_ColDiff.get();
		g_MaxColCount.get() = ColCount + 256;
		g_ColDiff.get() = new bool[g_MaxColCount.get()];
		}

	memset(g_ColDiff.get(), 0, g_ColCount.get()*sizeof(bool));
	for (unsigned i = 0; i < DiffColCount; ++i)
		{
		unsigned Col = DiffCols[i];
		assert(Col < ColCount);
		g_ColDiff.get()[Col] = true;
		}

	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		FindIntersectingGaps(msa, SeqIndex);

#if	TRACE
	{
	Log("\n");
	Log("Intersecting gaps:\n");
	Log("      ");
	for (unsigned Col = 0; Col < ColCount; ++Col)
		Log("%c", g_ColDiff.get()[Col] ? '*' : ' ');
	Log("\n");
	Log("      ");
	for (unsigned Col = 0; Col < ColCount; ++Col)
		Log("%d", Col%10);
	Log("\n");
	for (unsigned Seq = 0; Seq < SeqCount; ++Seq)
		{
		Log("%3d:  ", Seq);
		for (unsigned Col = 0; Col < ColCount; ++Col)
			Log("%c", msa.GetChar(Seq, Col));
		Log("  :: ");
		for (GAPINFO *GI = g_Gaps.get()[Seq]; GI; GI = GI->Next)
			Log(" (%d,%d)", GI->Start, GI->End);
		Log("  >%s\n", msa.GetSeqName(Seq));
		}
	Log("\n");
	}
#endif

	SCORE Score = 0;
	for (unsigned Seq1 = 0; Seq1 < SeqCount; ++Seq1)
		{
		const WEIGHT w1 = msa.GetSeqWeight(Seq1);
		for (unsigned Seq2 = Seq1 + 1; Seq2 < SeqCount; ++Seq2)
			{
			const WEIGHT w2 = msa.GetSeqWeight(Seq2);
//			const SCORE Pair = ScorePair(Seq1, Seq2);
			const SCORE Pair = ScoreSeqPairGaps(msa, Seq1, msa, Seq2);
			Score += w1*w2*Pair;
#if	TRACE
			Log("Seq1=%u Seq2=%u ScorePair=%.4g w1=%.4g w2=%.4g Sum=%.4g\n",
			  Seq1, Seq2, Pair, w1, w2, Score);
#endif
			}
		}

	return Score;
	}
} 

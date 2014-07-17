#include "libMUSCLE/muscle.h"
#include "libMUSCLE/msa.h"
#include "libMUSCLE/profile.h"
#include "libMUSCLE/objscore.h"
#include "libMUSCLE/refine.h"
#include "libMUSCLE/tree.h"
#include <vector>
#include <iostream>

namespace muscle {

void MergeBestCols(const SCORE Scores[], const unsigned BestCols[],
  unsigned uBestColCount, unsigned uWindowLength, unsigned AnchorCols[],
  unsigned *ptruAnchorColCount);

void FindBestColsGrade(const SCORE Score[], unsigned uCount,
  double dThreshold, unsigned BestCols[], unsigned *ptruBestColCount);
  
SCORE ScoreSeqPairLetters(const MSA &msa1, unsigned uSeqIndex1,
  const MSA &msa2, unsigned uSeqIndex2, SCORE MatchScore[] )
	{
	const unsigned uColCount = msa1.GetColCount();
	const unsigned uColCount2 = msa2.GetColCount();
	if (uColCount != uColCount2)
		Quit("ScoreSeqPairLetters, different lengths");

#if	TRACE_SEQPAIR
	{
	Log("\n");
	Log("ScoreSeqPairLetters\n");
	MSA msaTmp;
	msaTmp.SetSize(2, uColCount);
	msaTmp.CopySeq(0, msa1, uSeqIndex1);
	msaTmp.CopySeq(1, msa2, uSeqIndex2);
	msaTmp.LogMe();
	}
#endif

	SCORE scoreLetters = 0;
	SCORE scoreGaps = 0;
	bool bGapping1 = false;
	bool bGapping2 = false;

	unsigned uColStart = 0;
	bool bLeftTermGap = false;
	for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
		{
		bool bGap1 = msa1.IsGap(uSeqIndex1, uColIndex);
		bool bGap2 = msa2.IsGap(uSeqIndex2, uColIndex);
		if (!bGap1 || !bGap2)
			{
			if (bGap1 || bGap2)
				bLeftTermGap = true;
			uColStart = uColIndex;
			break;
			}
		}

	unsigned uColEnd = uColCount - 1;
	bool bRightTermGap = false;
	for (int iColIndex = (int) uColCount - 1; iColIndex >= 0; --iColIndex)
		{
		bool bGap1 = msa1.IsGap(uSeqIndex1, iColIndex);
		bool bGap2 = msa2.IsGap(uSeqIndex2, iColIndex);
		if (!bGap1 || !bGap2)
			{
			if (bGap1 || bGap2)
				bRightTermGap = true;
			uColEnd = (unsigned) iColIndex;
			break;
			}
		}

#if	TRACE_SEQPAIR
	Log("LeftTermGap=%d RightTermGap=%d\n", bLeftTermGap, bRightTermGap);
#endif

	for (unsigned uColIndex = uColStart; uColIndex <= uColEnd; ++uColIndex)
		{
		unsigned uLetter1 = msa1.GetLetterEx(uSeqIndex1, uColIndex);
		if (uLetter1 >= g_AlphaSize.get())
			continue;
		unsigned uLetter2 = msa2.GetLetterEx(uSeqIndex2, uColIndex);
		if (uLetter2 >= g_AlphaSize.get())
			continue;

		SCORE scoreMatch = (*g_ptrScoreMatrix.get())[uLetter1][uLetter2];
		scoreLetters += scoreMatch;
		if( MatchScore != NULL )
			MatchScore[uColIndex] = scoreMatch;
		}
	return scoreLetters;
	}

/* a version of ScoreSeqPairGaps that computes a per-residue score */
SCORE ScoreSeqPairGaps(const MSA &msa1, unsigned uSeqIndex1,
  const MSA &msa2, unsigned uSeqIndex2, SCORE MatchScore[] )
	{
	const unsigned uColCount = msa1.GetColCount();
	const unsigned uColCount2 = msa2.GetColCount();
	if (uColCount != uColCount2)
		Quit("ScoreSeqPairGaps, different lengths");

#if	TRACE_SEQPAIR
	{
	Log("\n");
	Log("ScoreSeqPairGaps\n");
	MSA msaTmp;
	msaTmp.SetSize(2, uColCount);
	msaTmp.CopySeq(0, msa1, uSeqIndex1);
	msaTmp.CopySeq(1, msa2, uSeqIndex2);
	msaTmp.LogMe();
	}
#endif

	SCORE scoreGaps = 0;
	bool bGapping1 = false;
	bool bGapping2 = false;

	unsigned uColStart = 0;
	bool bLeftTermGap = false;
	for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
		{
		bool bGap1 = msa1.IsGap(uSeqIndex1, uColIndex);
		bool bGap2 = msa2.IsGap(uSeqIndex2, uColIndex);
		if (!bGap1 || !bGap2)
			{
			if (bGap1 || bGap2)
				bLeftTermGap = true;
			uColStart = uColIndex;
			break;
			}
		}

	unsigned uColEnd = uColCount - 1;
	bool bRightTermGap = false;
	for (int iColIndex = (int) uColCount - 1; iColIndex >= 0; --iColIndex)
		{
		bool bGap1 = msa1.IsGap(uSeqIndex1, iColIndex);
		bool bGap2 = msa2.IsGap(uSeqIndex2, iColIndex);
		if (!bGap1 || !bGap2)
			{
			if (bGap1 || bGap2)
				bRightTermGap = true;
			uColEnd = (unsigned) iColIndex;
			break;
			}
		}

#if	TRACE_SEQPAIR
	Log("LeftTermGap=%d RightTermGap=%d\n", bLeftTermGap, bRightTermGap);
#endif

	unsigned gap_left_col = 0;
	SCORE cur_gap_score = 0;
	for (unsigned uColIndex = uColStart; uColIndex <= uColEnd; ++uColIndex)
		{
		bool bGap1 = msa1.IsGap(uSeqIndex1, uColIndex);
		bool bGap2 = msa2.IsGap(uSeqIndex2, uColIndex);

		if (bGap1 && bGap2)
			continue;

		if (bGap1)
			{
			if (!bGapping1)
				{
#if	TRACE_SEQPAIR
				Log("Gap open seq 1 col %d\n", uColIndex);
#endif
				gap_left_col = uColIndex;
				if (uColIndex == uColStart)
					{
					scoreGaps += TermGapScore(true);
					cur_gap_score += TermGapScore(true);
				}else{
					scoreGaps += g_scoreGapOpen.get();
					cur_gap_score += g_scoreGapOpen.get();
					}
				bGapping1 = true;
				}
			else
				{
				scoreGaps += g_scoreGapExtend.get();
				cur_gap_score += g_scoreGapExtend.get();
				}
			continue;
			}

		else if (bGap2)
			{
			if (!bGapping2)
				{
#if	TRACE_SEQPAIR
				Log("Gap open seq 2 col %d\n", uColIndex);
#endif
				gap_left_col = uColIndex;
				if (uColIndex == uColStart)
					{
					scoreGaps += TermGapScore(true);
					cur_gap_score += TermGapScore(true);
				}else{
					scoreGaps += g_scoreGapOpen.get();
					cur_gap_score += g_scoreGapOpen.get();
					}
				bGapping2 = true;
				}
			else
				{
				scoreGaps += g_scoreGapExtend.get();
				cur_gap_score += g_scoreGapExtend.get();
				}
			continue;
			}

		if( MatchScore != NULL && (bGapping1 || bGapping2) )
		{
			// spread the total gap penalty evenly across all columns
			SCORE per_site_penalty = cur_gap_score / (uColIndex-gap_left_col);
			for( unsigned uGapIndex = gap_left_col; uGapIndex < uColIndex; ++uGapIndex )
				{
				MatchScore[uGapIndex] = per_site_penalty;
				}
			gap_left_col = uInsane;
			cur_gap_score = 0;
		}
		bGapping1 = false;
		bGapping2 = false;
		}

	if (bGapping1 || bGapping2)
		{
		scoreGaps -= g_scoreGapOpen.get();
		scoreGaps += TermGapScore(true);
		cur_gap_score -= g_scoreGapOpen.get();
		cur_gap_score += TermGapScore(true);

		if( MatchScore != NULL )
			{
			// spread the total gap penalty evenly across all columns
			SCORE per_site_penalty = cur_gap_score / (uColCount-gap_left_col);
			for( unsigned uGapIndex = gap_left_col; uGapIndex < uColCount; ++uGapIndex )
				{
				MatchScore[uGapIndex] = per_site_penalty;
				}
			}
		}
	return scoreGaps;
	}




// this is a version of the profile x profile score that computes
// a per-site score suitable for use with anchoring heuristics
SCORE LetterObjScoreXP(const MSA &msa1, const MSA &msa2, SCORE MatchScore[])
	{
	const unsigned uColCount1 = msa1.GetColCount();
	const unsigned uColCount2 = msa2.GetColCount();
	if (uColCount1 != uColCount2)
		Quit("ObjScoreXP, alignment lengths differ %u %u", uColCount1, uColCount2);

	const unsigned uSeqCount1 = msa1.GetSeqCount();
	const unsigned uSeqCount2 = msa2.GetSeqCount();

#if	TRACE
	Log("     Score  Weight  Weight       Total\n");
	Log("----------  ------  ------  ----------\n");
#endif

	SCORE* mmScore = NULL;
	SCORE* ggScore = NULL;
	if( MatchScore != NULL )
		{
		mmScore = new SCORE[uColCount1];
		ggScore = new SCORE[uColCount1];
		memset( MatchScore, 0, sizeof(SCORE)*uColCount1 );
		}

	SCORE scoreTotal = 0;
	unsigned uPairCount = 0;
	for (unsigned uSeqIndex1 = 0; uSeqIndex1 < uSeqCount1; ++uSeqIndex1)
		{
		const WEIGHT w1 = msa1.GetSeqWeight(uSeqIndex1);
		for (unsigned uSeqIndex2 = 0; uSeqIndex2 < uSeqCount2; ++uSeqIndex2)
			{
			if( mmScore != NULL )
				memset( mmScore, 0, sizeof(SCORE)*uColCount1 );
			if( ggScore != NULL )
				memset( ggScore, 0, sizeof(SCORE)*uColCount1 );
			const WEIGHT w2 = msa2.GetSeqWeight(uSeqIndex2);
			const WEIGHT w = w1*w2;
			SCORE scoreLetters = ScoreSeqPairLetters(msa1, uSeqIndex1, msa2, uSeqIndex2, mmScore);
			SCORE scoreGaps = ScoreSeqPairGaps(msa1, uSeqIndex1, msa2, uSeqIndex2, ggScore);
			SCORE scorePair = scoreLetters + scoreGaps;
			scoreTotal += w*scorePair;
			++uPairCount;
			if( MatchScore != NULL )
				for( unsigned uColIndex = 0; uColIndex < uColCount1; ++uColIndex )
					MatchScore[uColIndex] += w*(mmScore[uColIndex]+ggScore[uColIndex]);
#if	TRACE
			Log("%10.2f  %6.3f  %6.3f  %10.2f  >%s >%s\n",
			  scorePair,
			  w1,
			  w2,
			  scorePair*w1*w2,
			  msa1.GetSeqName(uSeqIndex1),
			  msa2.GetSeqName(uSeqIndex2));
#endif
			}
		}
	if (0 == uPairCount)
		Quit("0 == uPairCount");

#if	TRACE
	Log("msa1=\n");
	msa1.LogMe();
	Log("msa2=\n");
	msa2.LogMe();
	Log("XP=%g\n", scoreTotal);
#endif
//	return scoreTotal / uPairCount;
	if(	mmScore != NULL )
		delete[] mmScore;
	if( ggScore != NULL )
		delete[] ggScore;

	return scoreTotal;
	}


// Best col only if all following criteria satisfied:
// (1) Score >= min
// (2) Smoothed score >= min
static void FindBestColsComboPP(unsigned uColCount, const SCORE Score[],
  const SCORE SmoothScore[], double dMinScore, double dMinSmoothScore,
  unsigned BestCols[], unsigned *ptruBestColCount)
	{

	unsigned uBestColCount = 0;
	for (unsigned uIndex = 0; uIndex < uColCount; ++uIndex)
		{
		if (Score[uIndex] < dMinScore)
			continue;
		if (SmoothScore[uIndex] < dMinSmoothScore)
			continue;
		BestCols[uBestColCount] = uIndex;
		++uBestColCount;
		}
	*ptruBestColCount = uBestColCount;
	}


void FindAnchorColsPP(const MSA &msa1, const MSA &msa2, unsigned AnchorCols[],
  unsigned *ptruAnchorColCount)
	{
	const unsigned uColCount = msa1.GetColCount();
	if( uColCount != msa2.GetColCount() )
		{
		*ptruAnchorColCount = 0;
		return;	// the profiles must have equal length to find anchor cols
		}

	SCORE *MatchScore = new SCORE[uColCount];
	SCORE *SmoothScore = new SCORE[uColCount];
	unsigned *BestCols = new unsigned[uColCount];

	LetterObjScoreXP(msa1, msa2, MatchScore);
	g_uSmoothWindowLength.get() = 21;	// this is better for DNA
	g_uAnchorSpacing.get() = 96;
	WindowSmooth(MatchScore, uColCount, g_uSmoothWindowLength.get(), SmoothScore,
	  g_dSmoothScoreCeil.get());


	unsigned uBestColCount;
//	FindBestColsGrade(SmoothScore,uColCount,.85,BestCols,&uBestColCount);
	FindBestColsComboPP(uColCount, MatchScore, SmoothScore, g_dMinBestColScore.get(), g_dMinSmoothScore.get(),
	  BestCols, &uBestColCount);
/*
	std::cerr << "found " << uBestColCount << " anchor cols:\n";
	for( size_t colI = 0; colI < uBestColCount; colI++ )
	{
		if( colI > 0 )
			std::cerr << ", ";
		std::cerr << BestCols[colI];
	}
	std::cerr << std::endl;
*/

#if	TRACE
	ListBestCols(msa, MatchScore, SmoothScore, BestCols, uBestColCount);
#endif

	MergeBestCols(MatchScore, BestCols, uBestColCount, g_uAnchorSpacing.get(), AnchorCols,
	  ptruAnchorColCount);
/*
	std::cerr << "\n\nafter merging, have " << *ptruAnchorColCount << " anchor cols:\n";
	for( size_t colI = 0; colI < *ptruAnchorColCount; colI++ )
	{
		if( colI > 0 )
			std::cerr << ", ";
		std::cerr << AnchorCols[colI];
	}
	std::cerr << std::endl;
*/
	delete[] MatchScore;
	delete[] SmoothScore;
	delete[] BestCols;
	}


void StripGapColumns( MSA& msa )
	{
	unsigned uCurCol = 0;
	for( unsigned uColIndex = 0; uColIndex < msa.GetColCount(); uColIndex++ )
		{
		if( !msa.IsGapColumn(uColIndex) )
			{
			for( unsigned uGapSeq = 0; uGapSeq < msa.GetSeqCount(); uGapSeq++ )
				{
				msa.SetChar(uGapSeq, uCurCol, msa.GetChar(uGapSeq,uColIndex));
				}
			uCurCol++;
			}
		}
	msa.DeleteColumns(uCurCol, msa.GetColCount()-uCurCol);
	}


void PrepareMSAforScoring( MSA& msa )
{
	Tree tree;
	const unsigned uSeqCount = msa.GetSeqCount();
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		msa.SetSeqId(uSeqIndex, uSeqIndex);

	TreeFromMSA(msa, tree, g_Cluster2.get(), g_Distance2.get(), g_Root1.get());
	SetMuscleTree(tree);
	SetMSAWeightsMuscle(msa);
}

// Return true if any changes made
void AnchoredProfileProfile(MSA &msa1, MSA &msa2, MSA &msaOut)
	{

	const unsigned uColCountIn = msa1.GetColCount();
	const unsigned uSeqCountIn = msa1.GetSeqCount() + msa2.GetSeqCount();

	unsigned *AnchorCols = new unsigned[uColCountIn];
	unsigned uAnchorColCount;

	PrepareMSAforScoring(msa1);
	PrepareMSAforScoring(msa2);
	FindAnchorColsPP(msa1, msa2, AnchorCols, &uAnchorColCount);

	const unsigned uRangeCount = uAnchorColCount + 1;
	Range *Ranges = new Range[uRangeCount];

#if	TRACE
	Log("%u ranges\n", uRangeCount);
#endif

	ColsToRanges(AnchorCols, uAnchorColCount, uColCountIn, Ranges);
	ListVertSavings(uColCountIn, uAnchorColCount, Ranges, uRangeCount);

#if	TRACE
	{
	Log("Anchor cols: ");
	for (unsigned i = 0; i < uAnchorColCount; ++i)
		Log(" %u", AnchorCols[i]);
	Log("\n");

	Log("Ranges:\n");
	for (unsigned i = 0; i < uRangeCount; ++i)
		Log("%4u - %4u\n", Ranges[i].m_uBestColLeft, Ranges[i].m_uBestColRight);
	}
#endif

	delete[] AnchorCols;

	msaOut.SetSize(uSeqCountIn, 0);

	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCountIn; ++uSeqIndex)
		{
		const char *ptrName;
		unsigned uId;
		if( uSeqIndex < msa1.GetSeqCount() )
			{
			msa1.SetSeqId(uSeqIndex, uSeqIndex);
			ptrName = msa1.GetSeqName(uSeqIndex);
			}
		else
			{
			msa2.SetSeqId(uSeqIndex-msa1.GetSeqCount(), uSeqIndex);
			ptrName = msa2.GetSeqName(uSeqIndex-msa1.GetSeqCount());
			}
		msaOut.SetSeqName(uSeqIndex, ptrName);
		msaOut.SetSeqId(uSeqIndex, uSeqIndex);
		}

	for (unsigned uRangeIndex = 0; uRangeIndex < uRangeCount; ++uRangeIndex)
		{
		MSA msaRange1;
		MSA msaRange2;
		MSA msaRangeOut;

		const Range &r = Ranges[uRangeIndex];

		const unsigned uFromColIndex = r.m_uBestColLeft;
		const unsigned uRangeColCount = r.m_uBestColRight - uFromColIndex;

		if (0 == uRangeColCount)
			continue;
/*		else if (1 == uRangeColCount)
			{
			MSAFromColRange(msaIn, uFromColIndex, 1, msaRange);
			MSAAppend(msaOut, msaRange);
			continue;
			}
*/
		MSAFromColRange(msa1, uFromColIndex, uRangeColCount, msaRange1);
		MSAFromColRange(msa2, uFromColIndex, uRangeColCount, msaRange2);
		StripGapColumns(msaRange1);
		StripGapColumns(msaRange2);

#if	TRACE
		Log("\n-------------\n");
		Log("Range %u - %u count=%u\n", r.m_uBestColLeft, r.m_uBestColRight, uRangeColCount);
		Log("Before:\n");
		msaRange1.LogMe();
		msaRange2.LogMe();
#endif

		ProfileProfile(msaRange1, msaRange2, msaRangeOut);

#if	TRACE
		Log("After:\n");
		msaRangeOut.LogMe();
#endif
		for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCountIn; ++uSeqIndex)
			msaRangeOut.SetSeqId(uSeqIndex, uSeqIndex);

		MSAAppend(msaOut, msaRangeOut);

#if	TRACE
		Log("msaOut after Cat:\n");
		msaOut.LogMe();
#endif
		}

	delete[] Ranges;
	}

} 

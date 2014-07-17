#include "libMUSCLE/muscle.h"
#include "libMUSCLE/profile.h"
#include "libMUSCLE/pwpath.h"
#include "libMUSCLE/threadstorage.h"

namespace muscle {

struct DP_MEMORY
	{
	unsigned uLength;
	SCORE *GapOpenA;
	SCORE *GapOpenB;
	SCORE *GapCloseA;
	SCORE *GapCloseB;
	SCORE *MPrev;
	SCORE *MCurr;
	SCORE *MWork;
	SCORE *DPrev;
	SCORE *DCurr;
	SCORE *DWork;
	SCORE **ScoreMxB;
	unsigned **SortOrderA;
	unsigned *uDeletePos;
	FCOUNT **FreqsA;
	int **TraceBack;
	};

static TLS<struct DP_MEMORY> DPM;

void FreeDPMemSPN()
	{
	const unsigned uOldLength = DPM.get().uLength;
	if (0 == uOldLength)
		return;

	for (unsigned i = 0; i < uOldLength; ++i)
		{
		delete[] DPM.get().TraceBack[i];
		delete[] DPM.get().FreqsA[i];
		delete[] DPM.get().SortOrderA[i];
		}
	for (unsigned n = 0; n < 4; ++n)
		delete[] DPM.get().ScoreMxB[n];

	delete[] DPM.get().MPrev;
	delete[] DPM.get().MCurr;
	delete[] DPM.get().MWork;
	delete[] DPM.get().DPrev;
	delete[] DPM.get().DCurr;
	delete[] DPM.get().DWork;
	delete[] DPM.get().uDeletePos;
	delete[] DPM.get().GapOpenA;
	delete[] DPM.get().GapOpenB;
	delete[] DPM.get().GapCloseA;
	delete[] DPM.get().GapCloseB;
	delete[] DPM.get().SortOrderA;
	delete[] DPM.get().FreqsA;
	delete[] DPM.get().ScoreMxB;
	delete[] DPM.get().TraceBack;
	}

static void AllocDPMem(unsigned uLengthA, unsigned uLengthB)
	{
// Max prefix length
	unsigned uLength = (uLengthA > uLengthB ? uLengthA : uLengthB) + 1;
	if (uLength < DPM.get().uLength)
		return;

// Add 256 to allow for future expansion and
// round up to next multiple of 32.
	uLength += 256;
	uLength += 32 - uLength%32;

	const unsigned uOldLength = DPM.get().uLength;
	if (uOldLength > 0)
		{
		for (unsigned i = 0; i < uOldLength; ++i)
			{
			delete[] DPM.get().TraceBack[i];
			delete[] DPM.get().FreqsA[i];
			delete[] DPM.get().SortOrderA[i];
			}
		for (unsigned n = 0; n < 4; ++n)
			delete[] DPM.get().ScoreMxB[n];

		delete[] DPM.get().MPrev;
		delete[] DPM.get().MCurr;
		delete[] DPM.get().MWork;
		delete[] DPM.get().DPrev;
		delete[] DPM.get().DCurr;
		delete[] DPM.get().DWork;
		delete[] DPM.get().uDeletePos;
		delete[] DPM.get().GapOpenA;
		delete[] DPM.get().GapOpenB;
		delete[] DPM.get().GapCloseA;
		delete[] DPM.get().GapCloseB;
		delete[] DPM.get().SortOrderA;
		delete[] DPM.get().FreqsA;
		delete[] DPM.get().ScoreMxB;
		delete[] DPM.get().TraceBack;
		}

	DPM.get().uLength = uLength;

	DPM.get().GapOpenA = new SCORE[uLength];
	DPM.get().GapOpenB = new SCORE[uLength];
	DPM.get().GapCloseA = new SCORE[uLength];
	DPM.get().GapCloseB = new SCORE[uLength];

	DPM.get().SortOrderA = new unsigned*[uLength];
	DPM.get().FreqsA = new FCOUNT*[uLength];
	DPM.get().ScoreMxB = new SCORE*[4];
	DPM.get().MPrev = new SCORE[uLength];
	DPM.get().MCurr = new SCORE[uLength];
	DPM.get().MWork = new SCORE[uLength];

	DPM.get().DPrev = new SCORE[uLength];
	DPM.get().DCurr = new SCORE[uLength];
	DPM.get().DWork = new SCORE[uLength];
	DPM.get().uDeletePos = new unsigned[uLength];

	DPM.get().TraceBack = new int*[uLength];

	for (unsigned uLetter = 0; uLetter < 4; ++uLetter)
		DPM.get().ScoreMxB[uLetter] = new SCORE[uLength];

	for (unsigned i = 0; i < uLength; ++i)
		{
		DPM.get().SortOrderA[i] = new unsigned[4];
		DPM.get().FreqsA[i] = new FCOUNT[4];
		DPM.get().TraceBack[i] = new int[uLength];
		}
	}

SCORE GlobalAlignSPN(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, PWPath &Path)
	{
	if (ALPHA_DNA != g_Alpha.get() || ALPHA_RNA == g_Alpha.get())
		Quit("GlobalAlignSPN: must be nucleo");

	const unsigned uPrefixCountA = uLengthA + 1;
	const unsigned uPrefixCountB = uLengthB + 1;

	AllocDPMem(uLengthA, uLengthB);

	SCORE *GapOpenA = DPM.get().GapOpenA;
	SCORE *GapOpenB = DPM.get().GapOpenB;
	SCORE *GapCloseA = DPM.get().GapCloseA;
	SCORE *GapCloseB = DPM.get().GapCloseB;

	unsigned **SortOrderA = DPM.get().SortOrderA;
	FCOUNT **FreqsA = DPM.get().FreqsA;
	SCORE **ScoreMxB = DPM.get().ScoreMxB;
	SCORE *MPrev = DPM.get().MPrev;
	SCORE *MCurr = DPM.get().MCurr;
	SCORE *MWork = DPM.get().MWork;

	SCORE *DPrev = DPM.get().DPrev;
	SCORE *DCurr = DPM.get().DCurr;
	SCORE *DWork = DPM.get().DWork;
	unsigned *uDeletePos = DPM.get().uDeletePos;

	int **TraceBack = DPM.get().TraceBack;

	for (unsigned i = 0; i < uLengthA; ++i)
		{
		GapOpenA[i] = PA[i].m_scoreGapOpen;
		GapCloseA[i] = PA[i].m_scoreGapClose;

		for (unsigned uLetter = 0; uLetter < 4; ++uLetter)
			{
			SortOrderA[i][uLetter] = PA[i].m_uSortOrder[uLetter];
			FreqsA[i][uLetter] = PA[i].m_fcCounts[uLetter];
			}
		}

	for (unsigned j = 0; j < uLengthB; ++j)
		{
		GapOpenB[j] = PB[j].m_scoreGapOpen;
		GapCloseB[j] = PB[j].m_scoreGapClose;
		}

	for (unsigned uLetter = 0; uLetter < 4; ++uLetter)
		{
		for (unsigned j = 0; j < uLengthB; ++j)
			ScoreMxB[uLetter][j] = PB[j].m_AAScores[uLetter];
		}

	for (unsigned i = 0; i < uPrefixCountA; ++i)
		memset(TraceBack[i], 0, uPrefixCountB*sizeof(int));

// Special case for i=0
	unsigned **ptrSortOrderA = SortOrderA;
	FCOUNT **ptrFreqsA = FreqsA;
	assert(ptrSortOrderA == &(SortOrderA[0]));
	assert(ptrFreqsA == &(FreqsA[0]));
	TraceBack[0][0] = 0;

	SCORE scoreSum = 0;
	unsigned *ptrSortOrderAi = SortOrderA[0];
	const unsigned *ptrSortOrderAEnd = ptrSortOrderAi + 4;
	FCOUNT *ptrFreqsAi = FreqsA[0];
	for (; ptrSortOrderAi != ptrSortOrderAEnd; ++ptrSortOrderAi)
		{
		const unsigned uLetter = *ptrSortOrderAi;
		const FCOUNT fcLetter = ptrFreqsAi[uLetter];
		if (0 == fcLetter)
			break;
		scoreSum += fcLetter*ScoreMxB[uLetter][0];
		}
	MPrev[0] = scoreSum - g_scoreCenter.get();

// D(0,0) is -infinity (requires I->D).
	DPrev[0] = MINUS_INFINITY;

	for (unsigned j = 1; j < uLengthB; ++j)
		{
	// Only way to get M(0, j) looks like this:
	//		A	----X
	//		B	XXXXX
	//			0   j
	// So gap-open at j=0, gap-close at j-1.
		SCORE scoreSum = 0;
		unsigned *ptrSortOrderAi = SortOrderA[0];
		const unsigned *ptrSortOrderAEnd = ptrSortOrderAi + 4;
		FCOUNT *ptrFreqsAi = FreqsA[0];
		for (; ptrSortOrderAi != ptrSortOrderAEnd; ++ptrSortOrderAi)
			{
			const unsigned uLetter = *ptrSortOrderAi;
			const FCOUNT fcLetter = ptrFreqsAi[uLetter];
			if (0 == fcLetter)
				break;
			scoreSum += fcLetter*ScoreMxB[uLetter][j];
			}
		MPrev[j] = scoreSum - g_scoreCenter.get() + GapOpenB[0] + GapCloseB[j-1];
		TraceBack[0][j] = -(int) j;

	// Assume no D->I transitions, then can't be a delete if only
	// one letter from A.
		DPrev[j] = MINUS_INFINITY;
		}

	SCORE IPrev_j_1;
	for (unsigned i = 1; i < uLengthA; ++i)
		{
		++ptrSortOrderA;
		++ptrFreqsA;
		assert(ptrSortOrderA == &(SortOrderA[i]));
		assert(ptrFreqsA == &(FreqsA[i]));

		SCORE *ptrMCurr_j = MCurr;
		memset(ptrMCurr_j, 0, uLengthB*sizeof(SCORE));
		const FCOUNT *FreqsAi = *ptrFreqsA;

		const unsigned *SortOrderAi = *ptrSortOrderA;
		const unsigned *ptrSortOrderAiEnd = SortOrderAi + 4;
		const SCORE *ptrMCurrMax = MCurr + uLengthB;
		for (const unsigned *ptrSortOrderAi = SortOrderAi;
		  ptrSortOrderAi != ptrSortOrderAiEnd;
		  ++ptrSortOrderAi)
			{
			const unsigned uLetter = *ptrSortOrderAi;
			SCORE *NSBR_Letter = ScoreMxB[uLetter];
			const FCOUNT fcLetter = FreqsAi[uLetter];
			if (0 == fcLetter)
				break;
			SCORE *ptrNSBR = NSBR_Letter;
			for (SCORE *ptrMCurr = MCurr; ptrMCurr != ptrMCurrMax; ++ptrMCurr)
				*ptrMCurr += fcLetter*(*ptrNSBR++);
			}

		for (unsigned j = 0; j < uLengthB; ++j)
			MCurr[j] -= g_scoreCenter.get();

		ptrMCurr_j = MCurr;
		unsigned *ptrDeletePos = uDeletePos;

	// Special case for j=0
	// Only way to get M(i, 0) looks like this:
	//			0   i
	//		A	XXXXX
	//		B	----X
	// So gap-open at i=0, gap-close at i-1.
		assert(ptrMCurr_j == &(MCurr[0]));
		*ptrMCurr_j += GapOpenA[0] + GapCloseA[i-1];

		++ptrMCurr_j;

		int *ptrTraceBack_ij = TraceBack[i];
		*ptrTraceBack_ij++ = (int) i;

		SCORE *ptrMPrev_j = MPrev;
		SCORE *ptrDPrev = DPrev;
		SCORE d = *ptrDPrev;
		SCORE DNew = *ptrMPrev_j + GapOpenA[i];
		if (DNew > d)
			{
			d = DNew;
			*ptrDeletePos = i;
			}

		SCORE *ptrDCurr = DCurr;

		assert(ptrDCurr == &(DCurr[0]));
		*ptrDCurr = d;

	// Can't have an insert if no letters from B
		IPrev_j_1 = MINUS_INFINITY;

		unsigned uInsertPos;
		const SCORE scoreGapOpenAi = GapOpenA[i];
		const SCORE scoreGapCloseAi_1 = GapCloseA[i-1];

		for (unsigned j = 1; j < uLengthB; ++j)
			{
		// Here, MPrev_j is preserved from previous
		// iteration so with current i,j is M[i-1][j-1]
			SCORE MPrev_j = *ptrMPrev_j;
			SCORE INew = MPrev_j + GapOpenB[j];
			if (INew > IPrev_j_1)
				{
				IPrev_j_1 = INew;
				uInsertPos = j;
				}

			SCORE scoreMax = MPrev_j;

			assert(ptrDPrev == &(DPrev[j-1]));
			SCORE scoreD = *ptrDPrev++ + scoreGapCloseAi_1;
			if (scoreD > scoreMax)
				{
				scoreMax = scoreD;
				assert(ptrDeletePos == &(uDeletePos[j-1]));
				*ptrTraceBack_ij = (int) i - (int) *ptrDeletePos;
				assert(*ptrTraceBack_ij > 0);
				}
			++ptrDeletePos;

			SCORE scoreI = IPrev_j_1 + GapCloseB[j-1];
			if (scoreI > scoreMax)
				{
				scoreMax = scoreI;
				*ptrTraceBack_ij = (int) uInsertPos - (int) j;
				assert(*ptrTraceBack_ij < 0);
				}

			assert(ptrSortOrderA == &(SortOrderA[i]));
			assert(ptrFreqsA == &(FreqsA[i]));

			*ptrMCurr_j += scoreMax;
			assert(ptrMCurr_j == &(MCurr[j]));
			++ptrMCurr_j;

			MPrev_j = *(++ptrMPrev_j);
			assert(ptrDPrev == &(DPrev[j]));
			SCORE d = *ptrDPrev;
			SCORE DNew = MPrev_j + scoreGapOpenAi;
			if (DNew > d)
				{
				d = DNew;
				assert(ptrDeletePos == &uDeletePos[j]);
				*ptrDeletePos = i;
				}
			assert(ptrDCurr + 1 == &(DCurr[j]));
			*(++ptrDCurr) = d;

			++ptrTraceBack_ij;
			}

		Rotate(MPrev, MCurr, MWork);
		Rotate(DPrev, DCurr, DWork);
		}

// Special case for i=uLengthA
	SCORE IPrev = MINUS_INFINITY;

	unsigned uInsertPos;

	for (unsigned j = 1; j < uLengthB; ++j)
		{
		SCORE INew = MPrev[j-1] + GapOpenB[j];
		if (INew > IPrev)
			{
			uInsertPos = j;
			IPrev = INew;
			}
		}

// Special case for i=uLengthA, j=uLengthB
	SCORE scoreMax = MPrev[uLengthB-1];
	int iTraceBack = 0;

	SCORE scoreD = DPrev[uLengthB-1] + GapCloseA[uLengthA-1];
	if (scoreD > scoreMax)
		{
		scoreMax = scoreD;
		iTraceBack = (int) uLengthA - (int) uDeletePos[uLengthB-1];
		}

	SCORE scoreI = IPrev + GapCloseB[uLengthB-1];
	if (scoreI > scoreMax)
		{
		scoreMax = scoreI;
		iTraceBack = (int) uInsertPos - (int) uLengthB;
		}

	TraceBack[uLengthA][uLengthB] = iTraceBack;

	TraceBackToPath(TraceBack, uLengthA, uLengthB, Path);

	return scoreMax;
	}
} 

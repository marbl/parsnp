/***
Conservation value for a column in an MSA is defined as the number
of times the most common letter appears divided by the number of
sequences.
***/

#include "libMUSCLE/muscle.h"
#include "libMUSCLE/msa.h"
#include <math.h>

namespace muscle {

double MSA::GetAvgCons() const
	{
	assert(GetSeqCount() > 0);
	double dSum = 0;
	unsigned uNonGapColCount = 0;
	for (unsigned uColIndex = 0; uColIndex < GetColCount(); ++uColIndex)
		{
		if (!IsGapColumn(uColIndex))
			{
			dSum += GetCons(uColIndex);
			++uNonGapColCount;
			}
		}
	assert(uNonGapColCount > 0);
	double dAvg = dSum / uNonGapColCount;
	assert(dAvg > 0 && dAvg <= 1);
	return dAvg;
	}

double MSA::GetCons(unsigned uColIndex) const
	{
	unsigned Counts[MAX_ALPHA];
	for (unsigned uLetter = 0; uLetter < g_AlphaSize.get(); ++uLetter)
		Counts[uLetter] = 0;

	unsigned uMaxCount = 0;
	for (unsigned uSeqIndex = 0; uSeqIndex < GetSeqCount(); ++uSeqIndex)
		{
		if (IsGap(uSeqIndex, uColIndex))
			continue;
		char c = GetChar(uSeqIndex, uColIndex);
		c = toupper(c);
		if ('X' == c || 'B' == c || 'Z' == c)
			continue;
		unsigned uLetter = GetLetter(uSeqIndex, uColIndex);
		unsigned uCount = Counts[uLetter] + 1;
		if (uCount > uMaxCount)
			uMaxCount = uCount;
		Counts[uLetter] = uCount;
		}

// Cons is undefined for all-gap column
	if (0 == uMaxCount)
		{
//		assert(false);
		return 1;
		}

	double dCons = (double) uMaxCount / (double) GetSeqCount();
	assert(dCons > 0 && dCons <= 1);
	return dCons;
	}

// Perecent identity of a pair of sequences.
// Positions with one or both gapped are ignored.
double MSA::GetPctIdentityPair(unsigned uSeqIndex1, unsigned uSeqIndex2) const
	{
	const unsigned uColCount = GetColCount();
	unsigned uPosCount = 0;
	unsigned uSameCount = 0;
	for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
		{
		const char c1 = GetChar(uSeqIndex1, uColIndex);
		const char c2 = GetChar(uSeqIndex2, uColIndex);
		if (IsGapChar(c1) || IsGapChar(c2))
			continue;
		if (c1 == c2)
			++uSameCount;
		++uPosCount;
		}
	if (0 == uPosCount)
		return 0;
	return (double) uSameCount / (double) uPosCount;
	}

extern unsigned ResidueGroup[];

// Perecent group identity of a pair of sequences.
// Positions with one or both gapped are ignored.
double MSA::GetPctGroupIdentityPair(unsigned uSeqIndex1,
  unsigned uSeqIndex2) const
	{

	const unsigned uColCount = GetColCount();
	unsigned uPosCount = 0;
	unsigned uSameCount = 0;
	for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
		{
		if (IsGap(uSeqIndex1, uColIndex))
			continue;
		if (IsGap(uSeqIndex2, uColIndex))
			continue;
		if (IsWildcard(uSeqIndex1, uColIndex))
			continue;
		if (IsWildcard(uSeqIndex2, uColIndex))
			continue;

		const unsigned uLetter1 = GetLetter(uSeqIndex1, uColIndex);
		const unsigned uLetter2 = GetLetter(uSeqIndex2, uColIndex);
		const unsigned uGroup1 = ResidueGroup[uLetter1];
		const unsigned uGroup2 = ResidueGroup[uLetter2];
		if (uGroup1 == uGroup2)
			++uSameCount;
		++uPosCount;
		}
	if (0 == uPosCount)
		return 0;
	return (double) uSameCount / (double) uPosCount;
	}
} 

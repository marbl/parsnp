#include "libMUSCLE/muscle.h"
#include "libMUSCLE/profile.h"
#include "libMUSCLE/diaglist.h"
#include "libMUSCLE/threadstorage.h"

namespace muscle {

#define TRACE	0

const unsigned KTUP = 5;
const unsigned KTUPS = 6*6*6*6*6;
static TLS<unsigned[KTUPS]> TuplePos;

static char *TupleToStr(int t)
	{
	static TLS<char[7]> s;
	int t1, t2, t3, t4, t5;

	t1 = t%6;
	t2 = (t/6)%6;
	t3 = (t/(6*6))%6;
	t4 = (t/(6*6*6))%6;
	t5 = (t/(6*6*6*6))%6;

	s.get()[4] = '0' + t1;
	s.get()[3] = '0' + t2;
	s.get()[2] = '0' + t3;
	s.get()[1] = '0' + t4;
	s.get()[0] = '0' + t5;
	return s.get();
	}

static unsigned GetTuple(const ProfPos *PP, unsigned uPos)
	{
	const unsigned t0 = PP[uPos].m_uResidueGroup;
	if (RESIDUE_GROUP_MULTIPLE == t0)
		return EMPTY;

	const unsigned t1 = PP[uPos+1].m_uResidueGroup;
	if (RESIDUE_GROUP_MULTIPLE == t1)
		return EMPTY;

	const unsigned t2 = PP[uPos+2].m_uResidueGroup;
	if (RESIDUE_GROUP_MULTIPLE == t2)
		return EMPTY;

	const unsigned t3 = PP[uPos+3].m_uResidueGroup;
	if (RESIDUE_GROUP_MULTIPLE == t3)
		return EMPTY;

	const unsigned t4 = PP[uPos+4].m_uResidueGroup;
	if (RESIDUE_GROUP_MULTIPLE == t4)
		return EMPTY;

	return t0 + t1*6 + t2*6*6 + t3*6*6*6 + t4*6*6*6*6;
	}

void FindDiags(const ProfPos *PX, unsigned uLengthX, const ProfPos *PY,
  unsigned uLengthY, DiagList &DL)
	{
	if (ALPHA_Amino != g_Alpha.get())
		Quit("FindDiags: requires amino acid alphabet");

	DL.Clear();

	if (uLengthX < 12 || uLengthY < 12)
		return;

// Set A to shorter profile, B to longer
	const ProfPos *PA;
	const ProfPos *PB;
	unsigned uLengthA;
	unsigned uLengthB;
	bool bSwap;
	if (uLengthX < uLengthY)
		{
		bSwap = false;
		PA = PX;
		PB = PY;
		uLengthA = uLengthX;
		uLengthB = uLengthY;
		}
	else
		{
		bSwap = true;
		PA = PY;
		PB = PX;
		uLengthA = uLengthY;
		uLengthB = uLengthX;
		}

// Build tuple map for the longer profile, B
	if (uLengthB < KTUP)
		Quit("FindDiags: profile too short");

	memset(TuplePos.get(), EMPTY, sizeof(TuplePos.get()));

	for (unsigned uPos = 0; uPos < uLengthB - KTUP; ++uPos)
		{
		const unsigned uTuple = GetTuple(PB, uPos);
		if (EMPTY == uTuple)
			continue;
		TuplePos.get()[uTuple] = uPos;
		}

// Find matches
	for (unsigned uPosA = 0; uPosA < uLengthA - KTUP; ++uPosA)
		{
		const unsigned uTuple = GetTuple(PA, uPosA);
		if (EMPTY == uTuple)
			continue;
		const unsigned uPosB = TuplePos.get()[uTuple];
		if (EMPTY == uPosB)
			continue;

	// This tuple is found in both profiles
		unsigned uStartPosA = uPosA;
		unsigned uStartPosB = uPosB;

	// Try to extend the match forwards
		unsigned uEndPosA = uPosA + KTUP - 1;
		unsigned uEndPosB = uPosB + KTUP - 1;
		for (;;)
			{
			if (uLengthA - 1 == uEndPosA || uLengthB - 1 == uEndPosB)
				break;
			const unsigned uAAGroupA = PA[uEndPosA+1].m_uResidueGroup;
			if (RESIDUE_GROUP_MULTIPLE == uAAGroupA)
				break;
			const unsigned uAAGroupB = PB[uEndPosB+1].m_uResidueGroup;
			if (RESIDUE_GROUP_MULTIPLE == uAAGroupB)
				break;
			if (uAAGroupA != uAAGroupB)
				break;
			++uEndPosA;
			++uEndPosB;
			}
		uPosA = uEndPosA;

#if	TRACE
		{
		Log("Match: A %4u-%4u   ", uStartPosA, uEndPosA);
		for (unsigned n = uStartPosA; n <= uEndPosA; ++n)
			Log("%c", 'A' + PA[n].m_uResidueGroup);
		Log("\n");
		Log("       B %4u-%4u   ", uStartPosB, uEndPosB);
		for (unsigned n = uStartPosB; n <= uEndPosB; ++n)
			Log("%c", 'A' + PB[n].m_uResidueGroup);
		Log("\n");
		}
#endif

		const unsigned uLength = uEndPosA - uStartPosA + 1;
		assert(uEndPosB - uStartPosB + 1 == uLength);

		if (uLength >= g_uMinDiagLength.get())
			{
			if (bSwap)
				DL.Add(uStartPosB, uStartPosA, uLength);
			else
				DL.Add(uStartPosA, uStartPosB, uLength);
			}
		}
	}
} 

#include "libMUSCLE/muscle.h"
#include "libMUSCLE/msa.h"
#include "libMUSCLE/profile.h"
#include "libMUSCLE/pwpath.h"
#include "libMUSCLE/textfile.h"
#include "libMUSCLE/timing.h"

namespace muscle {

SCORE AlignTwoMSAs(const MSA &msa1, const MSA &msa2, MSA &msaOut, PWPath &Path,
  bool bLockLeft, bool bLockRight)
	{
	const unsigned uLengthA = msa1.GetColCount();
	const unsigned uLengthB = msa2.GetColCount();

	ProfPos *PA = ProfileFromMSA(msa1);
	ProfPos *PB = ProfileFromMSA(msa2);

	if (bLockLeft)
		{
		PA[0].m_scoreGapOpen = MINUS_INFINITY;
		PB[0].m_scoreGapOpen = MINUS_INFINITY;
		}

	if (bLockRight)
		{
		PA[uLengthA-1].m_scoreGapClose = MINUS_INFINITY;
		PB[uLengthB-1].m_scoreGapClose = MINUS_INFINITY;
		}

	float r = (float) uLengthA/ (float) (uLengthB + 1); // +1 to prevent div 0
	if (r < 1)
		r = 1/r;

	SCORE Score = GlobalAlign(PA, uLengthA, PB, uLengthB, Path);

	AlignTwoMSAsGivenPath(Path, msa1, msa2, msaOut);

	delete[] PA;
	delete[] PB;

	return Score;
	}
} 

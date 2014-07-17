#include "libMUSCLE/muscle.h"
#include "libMUSCLE/pwpath.h"
#include "libMUSCLE/timing.h"
#include "libMUSCLE/textfile.h"
#include "libMUSCLE/msa.h"
#include "libMUSCLE/profile.h"

namespace muscle {

#if	VER_3_52

#if	TIMING
TICKS g_ticksDP = 0;
#endif

SCORE GlobalAlign(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, PWPath &Path)
	{
#if	TIMING
	TICKS t1 = GetClockTicks();
#endif
	SCORE Score = 0;
	if (g_bDiags.get())
		Score = GlobalAlignDiags(PA, uLengthA, PB, uLengthB, Path);
	else
		Score = GlobalAlignNoDiags(PA, uLengthA, PB, uLengthB, Path);
#if	TIMING
	TICKS t2 = GetClockTicks();
	g_ticksDP += (t2 - t1);
#endif
	return Score;
	}

SCORE GlobalAlignNoDiags(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, PWPath &Path)
	{
	if (g_bDimer.get())
		return GlobalAlignDimer(PA, uLengthA, PB, uLengthB, Path);

	switch (g_PPScore.get())
		{
	case PPSCORE_LE:
		return GlobalAlignLE(PA, uLengthA, PB, uLengthB, Path);

	case PPSCORE_SP:
	case PPSCORE_SV:
		return GlobalAlignSP(PA, uLengthA, PB, uLengthB, Path);

	case PPSCORE_SPN:
		return GlobalAlignSPN(PA, uLengthA, PB, uLengthB, Path);
		}

	Quit("Invalid PP score (GlobalAlignNoDiags)");
	return 0;
	}

#endif	// VER_3_52
} 

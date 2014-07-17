#include "libMUSCLE/muscle.h"
#include "libMUSCLE/profile.h"

namespace muscle {

extern void TomHydro(ProfPos *Prof, unsigned Length);

// Apply hydrophobicity heuristic to a profile
void Hydro(ProfPos *Prof, unsigned uLength)
	{
	if (ALPHA_Amino != g_Alpha.get())
		return;

	if (g_bTomHydro.get())
		{
		TomHydro(Prof, uLength);
		return;
		}

	if (0 == g_uHydrophobicRunLength.get())
		return;
	if (uLength <= g_uHydrophobicRunLength.get())
		return;


	unsigned uRunLength = 0;
	unsigned L2 = g_uHydrophobicRunLength.get()/2;
	for (unsigned uColIndex = L2; uColIndex < uLength - L2; ++uColIndex)
		{
		ProfPos &PP = Prof[uColIndex];
		bool bHydro = IsHydrophobic(PP.m_fcCounts);
		if (bHydro)
			{
			++uRunLength;
			if (uRunLength >= g_uHydrophobicRunLength.get())
				{
				Prof[uColIndex-L2].m_scoreGapOpen *= (SCORE) g_dHydroFactor.get();
				Prof[uColIndex-L2].m_scoreGapClose *= (SCORE) g_dHydroFactor.get();
				}
			}
		else
			uRunLength = 0;
		}
	}
} 

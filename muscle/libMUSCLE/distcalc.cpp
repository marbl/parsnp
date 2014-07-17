#include "libMUSCLE/muscle.h"
#include "libMUSCLE/distfunc.h"
#include "libMUSCLE/distcalc.h"
#include "libMUSCLE/msa.h"

namespace muscle {

extern double GetScoreDist(const MSA &msa, unsigned SeqIndex1, unsigned SeqIndex2);

void DistCalcDF::Init(const DistFunc &DF)
	{
	m_ptrDF = &DF;
	}

void DistCalcDF::CalcDistRange(unsigned i, dist_t Dist[]) const
	{
	for (unsigned j = 0; j < i; ++j)
		Dist[j] = m_ptrDF->GetDist(i, j);
	}

unsigned DistCalcDF::GetCount() const
	{
	return m_ptrDF->GetCount();
	}

unsigned DistCalcDF::GetId(unsigned i) const
	{
	return m_ptrDF->GetId(i);
	}

const char *DistCalcDF::GetName(unsigned i) const
	{
	return m_ptrDF->GetName(i);
	}

void DistCalcMSA::Init(const MSA &msa, DISTANCE Distance)
	{
	m_ptrMSA = &msa;
	m_Distance = Distance;
	}

void DistCalcMSA::CalcDistRange(unsigned i, dist_t Dist[]) const
	{
	for (unsigned j = 0; j < i; ++j)
		{
		switch (m_Distance)
			{
		case DISTANCE_PctIdKimura:
			{
			const float PctId = (float) m_ptrMSA->GetPctIdentityPair(i, j);
			Dist[j] = (float) KimuraDist(PctId);
			break;
			}
		case DISTANCE_PctIdLog:
			{
			const float PctId = (float) m_ptrMSA->GetPctIdentityPair(i, j);
			Dist[j] = (float) PctIdToMAFFTDist(PctId);
			break;
			}
		case DISTANCE_ScoreDist:
			{
			Dist[j] = (float) GetScoreDist(*m_ptrMSA, i, j);
			continue;
			}
		case DISTANCE_Edit:
			{
			const float PctId = (float) m_ptrMSA->GetPctIdentityPair(i, j);
			if (PctId > 1.0)
				Quit("Internal error, DISTANCE_Edit, pct id=%.3g", PctId);
			Dist[j] = (float) 1.0 - PctId;
			break;
			}
		default:
			Quit("DistCalcMSA: Invalid DISTANCE_%u", m_Distance);
			}
		}
	}

unsigned DistCalcMSA::GetCount() const
	{
	return m_ptrMSA->GetSeqCount();
	}

unsigned DistCalcMSA::GetId(unsigned i) const
	{
	return m_ptrMSA->GetSeqId(i);
	}

const char *DistCalcMSA::GetName(unsigned i) const
	{
	return m_ptrMSA->GetSeqName(i);
	}
} 

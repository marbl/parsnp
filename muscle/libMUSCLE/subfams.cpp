#include "libMUSCLE/muscle.h"
#include "libMUSCLE/distfunc.h"

namespace muscle {

const float INFINITY = float(1e29);
const unsigned NILL = uInsane;

static TLS<float *> ShortestPathEstimate;
static TLS<unsigned *> Predecessor;

static void GetMostDistantPair(DistFunc &DF, unsigned *ptrIndex1, unsigned *ptrIndex2)
	{
	const unsigned uNodeCount = DF.GetCount();
	if (uNodeCount < 2)
		Quit("GetMostDistantPair: < 2 seqs");

	float MaxDist = -1; 
	unsigned Index1 = uInsane;
	unsigned Index2 = uInsane;
	for (unsigned i = 0; i < uNodeCount; ++i)
		{
		for (unsigned j = i + 1; j < uNodeCount; ++j)
			{
			float d = DF.GetDist(i, j);
			if (d > MaxDist)
				{
				MaxDist = d;
				Index1 = i;
				Index2 = j;
				}
			}
		}

	assert(Index1 != uInsane);
	assert(Index2 != uInsane);

	*ptrIndex1 = Index1;
	*ptrIndex2 = Index2;
	}

static void InitializeSingleSource(DistFunc &DF, unsigned uIndex)
	{
	const unsigned uNodeCount = 0;

	for (unsigned i = 0; i < uNodeCount; ++i)
		{
		ShortestPathEstimate.get()[i] = INFINITY;
		Predecessor.get()[i] = NILL;
		}
	ShortestPathEstimate.get()[uIndex] = 0;
	}

static void Relax(DistFunc &DF, unsigned u, unsigned v)
	{
	float w = DF.GetDist(u, v);
	float d = ShortestPathEstimate.get()[u] + w;
	if (ShortestPathEstimate.get()[v] > d)
		{
		ShortestPathEstimate.get()[v] = d;
		Predecessor.get()[v] = u;
		}
	}

void ShortestPath(DistFunc &DF, unsigned uIndex)
	{
	}
} 

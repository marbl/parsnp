#include "libMUSCLE/muscle.h"
#include "libMUSCLE/tree.h"
#include "libMUSCLE/distcalc.h"

namespace muscle {

// UPGMA clustering in O(N^2) time and space.

#define	TRACE	0

#define	MIN(x, y)	((x) < (y) ? (x) : (y))
#define	MAX(x, y)	((x) > (y) ? (x) : (y))
#define	AVG(x, y)	(((x) + (y))/2)

static TLS<unsigned> g_uLeafCount;
static TLS<unsigned> g_uTriangleSize;
static TLS<unsigned> g_uInternalNodeCount;
static TLS<unsigned> g_uInternalNodeIndex;

// Triangular distance matrix is g_Dist.get(), which is allocated
// as a one-dimensional vector of length g_uTriangleSize.get().
// TriangleSubscript(i,j) maps row,column=i,j to the subscript
// into this vector.
// Row / column coordinates are a bit messy.
// Initially they are leaf indexes 0..N-1.
// But each time we create a new node (=new cluster, new subtree),
// we re-use one of the two rows that become available (the children
// of the new node). This saves memory.
// We keep track of this through the g_uNodeIndex.get() vector.
static TLS<dist_t *> g_Dist;

// Distance to nearest neighbor in row i of distance matrix.
// Subscript is distance matrix row.
static TLS<dist_t *> g_MinDist;

// Nearest neighbor to row i of distance matrix.
// Subscript is distance matrix row.
static TLS<unsigned *> g_uNearestNeighbor;

// Node index of row i in distance matrix.
// Node indexes are 0..N-1 for leaves, N..2N-2 for internal nodes.
// Subscript is distance matrix row.
static TLS<unsigned *> g_uNodeIndex;

// The following vectors are defined on internal nodes,
// subscripts are internal node index 0..N-2.
// For g_uLeft.get()/Right, value is the node index 0 .. 2N-2
// because a child can be internal or leaf.
static TLS<unsigned *> g_uLeft;
static TLS<unsigned *> g_uRight;
static TLS<dist_t *> g_Height;
static TLS<dist_t *> g_LeftLength;
static TLS<dist_t *> g_RightLength;

static inline unsigned TriangleSubscript(unsigned uIndex1, unsigned uIndex2)
	{
#if	DEBUG
	if (uIndex1 >= g_uLeafCount.get() || uIndex2 >= g_uLeafCount.get())
		Quit("TriangleSubscript(%u,%u) %u", uIndex1, uIndex2, g_uLeafCount.get());
#endif
	unsigned v;
	if (uIndex1 >= uIndex2)
		v = uIndex2 + (uIndex1*(uIndex1 - 1))/2;
	else
		v = uIndex1 + (uIndex2*(uIndex2 - 1))/2;
	assert(v < (g_uLeafCount.get()*(g_uLeafCount.get() - 1))/2);
	return v;
	}

static void ListState()
	{
	Log("Dist matrix\n");
	Log("     ");
	for (unsigned i = 0; i < g_uLeafCount.get(); ++i)
		{
		if (uInsane == g_uNodeIndex.get()[i])
			continue;
		Log("  %5u", g_uNodeIndex.get()[i]);
		}
	Log("\n");

	for (unsigned i = 0; i < g_uLeafCount.get(); ++i)
		{
		if (uInsane == g_uNodeIndex.get()[i])
			continue;
		Log("%5u  ", g_uNodeIndex.get()[i]);
		for (unsigned j = 0; j < g_uLeafCount.get(); ++j)
			{
			if (uInsane == g_uNodeIndex.get()[j])
				continue;
			if (i == j)
				Log("       ");
			else
				{
				unsigned v = TriangleSubscript(i, j);
				Log("%5.2g  ", g_Dist.get()[v]);
				}
			}
		Log("\n");
		}

	Log("\n");
	Log("    i   Node   NrNb      Dist\n");
	Log("-----  -----  -----  --------\n");
	for (unsigned i = 0; i < g_uLeafCount.get(); ++i)
		{
		if (uInsane == g_uNodeIndex.get()[i])
			continue;
		Log("%5u  %5u  %5u  %8.3f\n",
		  i,
		  g_uNodeIndex.get()[i],
		  g_uNearestNeighbor.get()[i],
		  g_MinDist.get()[i]);
		}

	Log("\n");
	Log(" Node      L      R  Height  LLength  RLength\n");
	Log("-----  -----  -----  ------  -------  -------\n");
	for (unsigned i = 0; i <= g_uInternalNodeIndex.get(); ++i)
		Log("%5u  %5u  %5u  %6.2g  %6.2g  %6.2g\n",
		  i,
		  g_uLeft.get()[i],
		  g_uRight.get()[i],
		  g_Height.get()[i],
		  g_LeftLength.get()[i],
		  g_RightLength.get()[i]);
	}

void UPGMA2(const DistCalc &DC, Tree &tree, LINKAGE Linkage)
	{
	g_uLeafCount.get() = DC.GetCount();

	g_uTriangleSize.get() = (g_uLeafCount.get()*(g_uLeafCount.get() - 1))/2;
	g_uInternalNodeCount.get() = g_uLeafCount.get() - 1;

	g_Dist.get() = new dist_t[g_uTriangleSize.get()];

	g_uNodeIndex.get() = new unsigned[g_uLeafCount.get()];
	g_uNearestNeighbor.get() = new unsigned[g_uLeafCount.get()];
	g_MinDist.get() = new dist_t[g_uLeafCount.get()];
	unsigned *Ids = new unsigned [g_uLeafCount.get()];
	char **Names = new char *[g_uLeafCount.get()];

	g_uLeft.get() = new unsigned[g_uInternalNodeCount.get()];
	g_uRight.get() = new unsigned[g_uInternalNodeCount.get()];
	g_Height.get() = new dist_t[g_uInternalNodeCount.get()];
	g_LeftLength.get() = new dist_t[g_uInternalNodeCount.get()];
	g_RightLength.get() = new dist_t[g_uInternalNodeCount.get()];

	for (unsigned i = 0; i < g_uLeafCount.get(); ++i)
		{
		g_MinDist.get()[i] = BIG_DIST;
		g_uNodeIndex.get()[i] = i;
		g_uNearestNeighbor.get()[i] = uInsane;
		Ids[i] = DC.GetId(i);
		Names[i] = strsave(DC.GetName(i));
		}

	for (unsigned i = 0; i < g_uInternalNodeCount.get(); ++i)
		{
		g_uLeft.get()[i] = uInsane;
		g_uRight.get()[i] = uInsane;
		g_LeftLength.get()[i] = BIG_DIST;
		g_RightLength.get()[i] = BIG_DIST;
		g_Height.get()[i] = BIG_DIST;
		}

// Compute initial NxN triangular distance matrix.
// Store minimum distance for each full (not triangular) row.
// Loop from 1, not 0, because "row" is 0, 1 ... i-1,
// so nothing to do when i=0.
	for (unsigned i = 1; i < g_uLeafCount.get(); ++i)
		{
		dist_t *Row = g_Dist.get() + TriangleSubscript(i, 0);
		DC.CalcDistRange(i, Row);
		for (unsigned j = 0; j < i; ++j)
			{
			const dist_t d = Row[j];
			if (d < g_MinDist.get()[i])
				{
				g_MinDist.get()[i] = d;
				g_uNearestNeighbor.get()[i] = j;
				}
			if (d < g_MinDist.get()[j])
				{
				g_MinDist.get()[j] = d;
				g_uNearestNeighbor.get()[j] = i;
				}
			}
		}

#if	TRACE
	Log("Initial state:\n");
	ListState();
#endif

	for (g_uInternalNodeIndex.get() = 0; g_uInternalNodeIndex.get() < g_uLeafCount.get() - 1;
	  ++g_uInternalNodeIndex.get())
		{
#if	TRACE
		Log("\n");
		Log("Internal node index %5u\n", g_uInternalNodeIndex.get());
		Log("-------------------------\n");
#endif

	// Find nearest neighbors
		unsigned Lmin = uInsane;
		unsigned Rmin = uInsane;
		dist_t dtMinDist = BIG_DIST;
		for (unsigned j = 0; j < g_uLeafCount.get(); ++j)
			{
			if (uInsane == g_uNodeIndex.get()[j])
				continue;

			dist_t d = g_MinDist.get()[j];
			if (d < dtMinDist)
				{
				dtMinDist = d;
				Lmin = j;
				Rmin = g_uNearestNeighbor.get()[j];
				assert(uInsane != Rmin);
				assert(uInsane != g_uNodeIndex.get()[Rmin]);
				}
			}

		assert(Lmin != uInsane);
		assert(Rmin != uInsane);
		assert(dtMinDist != BIG_DIST);

#if	TRACE
		Log("Nearest neighbors Lmin %u[=%u] Rmin %u[=%u] dist %.3g\n",
		  Lmin,
		  g_uNodeIndex.get()[Lmin],
		  Rmin,
		  g_uNodeIndex.get()[Rmin],
		  dtMinDist);
#endif

	// Compute distances to new node
	// New node overwrites row currently assigned to Lmin
		dist_t dtNewMinDist = BIG_DIST;
		unsigned uNewNearestNeighbor = uInsane;
		for (unsigned j = 0; j < g_uLeafCount.get(); ++j)
			{
			if (j == Lmin || j == Rmin)
				continue;
			if (uInsane == g_uNodeIndex.get()[j])
				continue;

			const unsigned vL = TriangleSubscript(Lmin, j);
			const unsigned vR = TriangleSubscript(Rmin, j);
			const dist_t dL = g_Dist.get()[vL];
			const dist_t dR = g_Dist.get()[vR];
			dist_t dtNewDist;

			switch (Linkage)
				{
			case LINKAGE_Avg:
				dtNewDist = AVG(dL, dR);
				break;

			case LINKAGE_Min:
				dtNewDist = MIN(dL, dR);
				break;

			case LINKAGE_Max:
				dtNewDist = MAX(dL, dR);
				break;

			case LINKAGE_Biased:
				dtNewDist = g_dSUEFF.get()*AVG(dL, dR) + (1 - g_dSUEFF.get())*MIN(dL, dR);
				break;

			default:
				Quit("UPGMA2: Invalid LINKAGE_%u", Linkage);
				}

		// Nasty special case.
		// If nearest neighbor of j is Lmin or Rmin, then make the new
		// node (which overwrites the row currently occupied by Lmin)
		// the nearest neighbor. This situation can occur when there are
		// equal distances in the matrix. If we don't make this fix,
		// the nearest neighbor pointer for j would become invalid.
		// (We don't need to test for == Lmin, because in that case
		// the net change needed is zero due to the change in row
		// numbering).
			if (g_uNearestNeighbor.get()[j] == Rmin)
				g_uNearestNeighbor.get()[j] = Lmin;

#if	TRACE
			Log("New dist to %u = (%u/%.3g + %u/%.3g)/2 = %.3g\n",
			  j, Lmin, dL, Rmin, dR, dtNewDist);
#endif
			g_Dist.get()[vL] = dtNewDist;
			if (dtNewDist < dtNewMinDist)
				{
				dtNewMinDist = dtNewDist;
				uNewNearestNeighbor = j;
				}
			}

		assert(g_uInternalNodeIndex.get() < g_uLeafCount.get() - 1 || BIG_DIST != dtNewMinDist);
		assert(g_uInternalNodeIndex.get() < g_uLeafCount.get() - 1 || uInsane != uNewNearestNeighbor);

		const unsigned v = TriangleSubscript(Lmin, Rmin);
		const dist_t dLR = g_Dist.get()[v];
		const dist_t dHeightNew = dLR/2;
		const unsigned uLeft = g_uNodeIndex.get()[Lmin];
		const unsigned uRight = g_uNodeIndex.get()[Rmin];
		const dist_t HeightLeft =
		  uLeft < g_uLeafCount.get() ? 0 : g_Height.get()[uLeft - g_uLeafCount.get()];
		const dist_t HeightRight =
		  uRight < g_uLeafCount.get() ? 0 : g_Height.get()[uRight - g_uLeafCount.get()];

		g_uLeft.get()[g_uInternalNodeIndex.get()] = uLeft;
		g_uRight.get()[g_uInternalNodeIndex.get()] = uRight;
		g_LeftLength.get()[g_uInternalNodeIndex.get()] = dHeightNew - HeightLeft;
		g_RightLength.get()[g_uInternalNodeIndex.get()] = dHeightNew - HeightRight;
		g_Height.get()[g_uInternalNodeIndex.get()] = dHeightNew;

	// Row for left child overwritten by row for new node
		g_uNodeIndex.get()[Lmin] = g_uLeafCount.get() + g_uInternalNodeIndex.get();
		g_uNearestNeighbor.get()[Lmin] = uNewNearestNeighbor;
		g_MinDist.get()[Lmin] = dtNewMinDist;

	// Delete row for right child
		g_uNodeIndex.get()[Rmin] = uInsane;

#if	TRACE
		Log("\nInternalNodeIndex=%u Lmin=%u Rmin=%u\n",
		  g_uInternalNodeIndex.get(), Lmin, Rmin);
		ListState();
#endif
		}

	unsigned uRoot = g_uLeafCount.get() - 2;
	tree.Create(g_uLeafCount.get(), uRoot, g_uLeft.get(), g_uRight.get(), g_LeftLength.get(), g_RightLength.get(),
	  Ids, Names);

#if	TRACE
	tree.LogMe();
#endif

	delete[] g_Dist.get();

	delete[] g_uNodeIndex.get();
	delete[] g_uNearestNeighbor.get();
	delete[] g_MinDist.get();
	delete[] g_Height.get();

	delete[] g_uLeft.get();
	delete[] g_uRight.get();
	delete[] g_LeftLength.get();
	delete[] g_RightLength.get();
	
	for (unsigned i = 0; i < g_uLeafCount.get(); ++i)
		free(Names[i]);
	delete[] Names;
	delete[] Ids;
	}

class DistCalcTest : public DistCalc
	{
	virtual void CalcDistRange(unsigned i, dist_t Dist[]) const
		{
		static dist_t TestDist[5][5] =
			{
			0,  2, 14, 14, 20,
			2,  0, 14, 14, 20,
			14, 14,  0,  4, 20,
			14, 14,  4,  0, 20,
			20, 20, 20, 20,  0,
			};
		for (unsigned j = 0; j < i; ++j)
			Dist[j] = TestDist[i][j];
		}
	virtual unsigned GetCount() const
		{
		return 5;
		}
	virtual unsigned GetId(unsigned i) const
		{
		return i;
		}
	virtual const char *GetName(unsigned i) const
		{
		return "name";
		}
	};

void Test()
	{
	SetListFileName("c:\\tmp\\lobster.log", false);
	DistCalcTest DC;
	Tree tree;
	UPGMA2(DC, tree, LINKAGE_Avg);
	}
} 

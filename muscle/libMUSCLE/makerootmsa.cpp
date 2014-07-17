#include "libMUSCLE/muscle.h"
#include "libMUSCLE/tree.h"
#include "libMUSCLE/seqvect.h"
#include "libMUSCLE/profile.h"
#include "libMUSCLE/msa.h"
#include "libMUSCLE/pwpath.h"
#include "libMUSCLE/estring.h"

namespace muscle {

#define TRACE		0
#define VALIDATE	0

static void PathSeq(const Seq &s, const PWPath &Path, bool bRight, Seq &sOut)
	{
	short *esA;
	short *esB;
	PathToEstrings(Path, &esA, &esB);

	const unsigned uSeqLength = s.Length();
	const unsigned uEdgeCount = Path.GetEdgeCount();

	sOut.Clear();
	sOut.SetName(s.GetName());
	unsigned uPos = 0;
	for (unsigned uEdgeIndex = 0; uEdgeIndex < uEdgeCount; ++uEdgeIndex)
		{
		const PWEdge &Edge = Path.GetEdge(uEdgeIndex);
		char cType = Edge.cType;
		if (bRight)
			{
			if (cType == 'I')
				cType = 'D';
			else if (cType == 'D')
				cType = 'I';
			}
		switch (cType)
			{
		case 'M':
			sOut.AppendChar(s[uPos++]);
			break;
		case 'D':
			sOut.AppendChar('-');
			break;
		case 'I':
			sOut.AppendChar(s[uPos++]);
			break;
		default:
			Quit("PathSeq, invalid edge type %c", cType);
			}
		}
	}

#if	VALIDATE

static void MakeRootSeq(const Seq &s, const Tree &GuideTree, unsigned uLeafNodeIndex,
  const ProgNode Nodes[], Seq &sRoot)
	{
	sRoot.Copy(s);
	unsigned uNodeIndex = uLeafNodeIndex;
	for (;;)
		{
	  	unsigned uParent = GuideTree.GetParent(uNodeIndex);
		if (NULL_NEIGHBOR == uParent)
			break;
		bool bRight = (GuideTree.GetLeft(uParent) == uNodeIndex);
		uNodeIndex = uParent;
		const PWPath &Path = Nodes[uNodeIndex].m_Path;
		Seq sTmp;
		PathSeq(sRoot, Path, bRight, sTmp);
		sRoot.Copy(sTmp);
		}
	}

#endif	// VALIDATE

static short *MakeRootSeqE(const Seq &s, const Tree &GuideTree, unsigned uLeafNodeIndex,
  const ProgNode Nodes[], Seq &sRoot, short *Estring1, short *Estring2)
	{
	short *EstringCurr = Estring1;
	short *EstringNext = Estring2;

	const unsigned uSeqLength = s.Length();
	EstringCurr[0] = uSeqLength;
	EstringCurr[1] = 0;

	unsigned uNodeIndex = uLeafNodeIndex;
	for (;;)
		{
	  	unsigned uParent = GuideTree.GetParent(uNodeIndex);
		if (NULL_NEIGHBOR == uParent)
			break;
		bool bRight = (GuideTree.GetLeft(uParent) == uNodeIndex);
		uNodeIndex = uParent;
		const PWPath &Path = Nodes[uNodeIndex].m_Path;
		const short *EstringNode = bRight ?
		  Nodes[uNodeIndex].m_EstringL : Nodes[uNodeIndex].m_EstringR;

		MulEstrings(EstringCurr, EstringNode, EstringNext);
#if	TRACE
		Log("\n");
		Log("Curr=");
		LogEstring(EstringCurr);
		Log("\n");
		Log("Node=");
		LogEstring(EstringNode);
		Log("\n");
		Log("Prod=");
		LogEstring(EstringNext);
		Log("\n");
#endif
		short *EstringTmp = EstringNext;
		EstringNext = EstringCurr;
		EstringCurr = EstringTmp;
		}
	EstringOp(EstringCurr, s, sRoot);

#if	TRACE
	Log("Root estring=");
	LogEstring(EstringCurr);
	Log("\n");
	Log("Root seq=");
	sRoot.LogMe();
#endif
	return EstringCurr;
	}

static unsigned GetFirstNodeIndex(const Tree &tree)
	{
	if (g_bStable.get())
		return 0;
	return tree.FirstDepthFirstNode();
	}

static unsigned GetNextNodeIndex(const Tree &tree, unsigned uPrevNodeIndex)
	{
	if (g_bStable.get())
		{
		const unsigned uNodeCount = tree.GetNodeCount();
		unsigned uNodeIndex = uPrevNodeIndex;
		for (;;)
			{
			++uNodeIndex;
			if (uNodeIndex >= uNodeCount)
				return NULL_NEIGHBOR;
			if (tree.IsLeaf(uNodeIndex))
				return uNodeIndex;
			}
		}
	unsigned uNodeIndex = uPrevNodeIndex;
	for (;;)
		{
		uNodeIndex = tree.NextDepthFirstNode(uNodeIndex);
		if (NULL_NEIGHBOR == uNodeIndex || tree.IsLeaf(uNodeIndex))
			return uNodeIndex;
		}
	}

void MakeRootMSA(const SeqVect &v, const Tree &GuideTree, ProgNode Nodes[],
  MSA &a)
	{
#if	TRACE
	Log("MakeRootMSA Tree=");
	GuideTree.LogMe();
#endif
	const unsigned uSeqCount = v.GetSeqCount();
	unsigned uColCount = uInsane;
	unsigned uSeqIndex = 0;
	const unsigned uTreeNodeCount = GuideTree.GetNodeCount();
	const unsigned uRootNodeIndex = GuideTree.GetRootNodeIndex();
	const PWPath &RootPath = Nodes[uRootNodeIndex].m_Path;
	const unsigned uRootColCount = RootPath.GetEdgeCount();
	const unsigned uEstringSize = uRootColCount + 1;
	short *Estring1 = new short[uEstringSize];
	short *Estring2 = new short[uEstringSize];
	SetProgressDesc("Root alignment");

	unsigned uTreeNodeIndex = GetFirstNodeIndex(GuideTree);
	do
		{
		Progress(uSeqIndex, uSeqCount);

		unsigned uId = GuideTree.GetLeafId(uTreeNodeIndex);
		const Seq &s = *(v[uId]);

		Seq sRootE;
		short *es = MakeRootSeqE(s, GuideTree, uTreeNodeIndex, Nodes, sRootE,
		  Estring1, Estring2);
		Nodes[uTreeNodeIndex].m_EstringL = EstringNewCopy(es);

#if	VALIDATE
		Seq sRoot;
		MakeRootSeq(s, GuideTree, uTreeNodeIndex, Nodes, sRoot);
		if (!sRoot.Eq(sRootE))
			{
			Log("sRoot=");
			sRoot.LogMe();
			Log("sRootE=");
			sRootE.LogMe();
			Quit("Root seqs differ");
			}
#endif

#if	TRACE
		Log("MakeRootSeq=\n");
		sRoot.LogMe();
#endif
		if (uInsane == uColCount)
			{
			uColCount = sRootE.Length();
			a.SetSize(uSeqCount, uColCount);
			}
		else
			{
			assert(uColCount == sRootE.Length());
			}
		a.SetSeqName(uSeqIndex, s.GetName());
		a.SetSeqId(uSeqIndex, uId);
		for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
			a.SetChar(uSeqIndex, uColIndex, sRootE[uColIndex]);
		++uSeqIndex;

		uTreeNodeIndex = GetNextNodeIndex(GuideTree, uTreeNodeIndex);
		}
	while (NULL_NEIGHBOR != uTreeNodeIndex);

	delete[] Estring1;
	delete[] Estring2;

	ProgressStepsDone();
	assert(uSeqIndex == uSeqCount);
	}
} 

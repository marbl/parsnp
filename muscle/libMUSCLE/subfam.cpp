#include "libMUSCLE/muscle.h"
#include "libMUSCLE/tree.h"
#include "libMUSCLE/textfile.h"	// for test only
#include "libMUSCLE/msa.h"
#include "libMUSCLE/seqvect.h"
#include "libMUSCLE/profile.h"
#ifndef _MSC_VER
#include <unistd.h>	//	for unlink
#endif

namespace muscle {

#define TRACE	0

/***
Find subfamilies from tree by following criteria:

(a) number of leaves <= max,
(b) is monophyletic, i.e. most recent common ancestor is parent
of no more than one subfamily.
***/

static unsigned SubFamRecurse(const Tree &tree, unsigned uNodeIndex, unsigned uMaxLeafCount,
  unsigned SubFams[], unsigned &uSubFamCount)
	{
	if (tree.IsLeaf(uNodeIndex))
		return 1;

	unsigned uLeft = tree.GetLeft(uNodeIndex);
	unsigned uRight = tree.GetRight(uNodeIndex);
	unsigned uLeftCount = SubFamRecurse(tree, uLeft, uMaxLeafCount, SubFams, uSubFamCount);
	unsigned uRightCount = SubFamRecurse(tree, uRight, uMaxLeafCount, SubFams, uSubFamCount);

	unsigned uLeafCount = uLeftCount + uRightCount;
	if (uLeftCount + uRightCount > uMaxLeafCount)
		{
		if (uLeftCount <= uMaxLeafCount)
			SubFams[uSubFamCount++] = uLeft;
		if (uRightCount <= uMaxLeafCount)
			SubFams[uSubFamCount++] = uRight;
		}
	else if (tree.IsRoot(uNodeIndex))
		{
		if (uSubFamCount != 0)
			Quit("Error in SubFamRecurse");
		SubFams[uSubFamCount++] = uNodeIndex;
		}

	return uLeafCount;
	}

void SubFam(const Tree &tree, unsigned uMaxLeafCount, unsigned SubFams[], unsigned *ptruSubFamCount)
	{
	*ptruSubFamCount = 0;
	SubFamRecurse(tree, tree.GetRootNodeIndex(), uMaxLeafCount, SubFams, *ptruSubFamCount);

#if	TRACE
	{
	Log("\n");
	Log("Tree:\n");
	tree.LogMe();
	//void DrawTree(const Tree &tree);
	//DrawTree(tree);
	Log("\n");
	Log("%d subfams:\n", *ptruSubFamCount);
	for (unsigned i = 0; i < *ptruSubFamCount; ++i)
		Log("  %d=%d", i, SubFams[i]);
	Log("\n");
	}
#endif
	}

//unsigned SubFams[9999];
//unsigned uSubFamCount;
//
//static unsigned DistFromRoot(const Tree &tree, unsigned uNodeIndex)
//	{
//	const unsigned uRoot = tree.GetRootNodeIndex();
//	unsigned uDist = 0;
//	while (uNodeIndex != uRoot)
//		{
//		++uDist;
//		uNodeIndex = tree.GetParent(uNodeIndex);
//		}
//	return uDist;
//	}
//
//static void DrawNode(const Tree &tree, unsigned uNodeIndex)
//	{
//	if (!tree.IsLeaf(uNodeIndex))
//		DrawNode(tree, tree.GetLeft(uNodeIndex));
//
//	unsigned uDist = DistFromRoot(tree, uNodeIndex);
//	for (unsigned i = 0; i < 5*uDist; ++i)
//		Log(" ");
//	Log("%d", uNodeIndex);
//	for (unsigned i = 0; i < uSubFamCount; ++i)
//		if (uNodeIndex == SubFams[i])
//			{
//			Log("*");
//			break;
//			}
//	Log("\n");
//
//	if (!tree.IsLeaf(uNodeIndex))
//		DrawNode(tree, tree.GetRight(uNodeIndex));
//	}
//
//static void DrawTree(const Tree &tree)
//	{
//	unsigned uRoot = tree.GetRootNodeIndex();
//	DrawNode(tree, uRoot);
//	}
//
//void TestSubFams(const char *FileName)
//	{
//	Tree tree;
//	TextFile f(FileName);
//	tree.FromFile(f);
//	SubFam(tree, 5, SubFams, &uSubFamCount);
//	DrawTree(tree);
//	}

static void SetInFam(const Tree &tree, unsigned uNodeIndex, bool NodeInSubFam[])
	{
	if (tree.IsLeaf(uNodeIndex))
		return;
	unsigned uLeft = tree.GetLeft(uNodeIndex);
	unsigned uRight = tree.GetRight(uNodeIndex);
	NodeInSubFam[uLeft] = true;
	NodeInSubFam[uRight] = true;

	SetInFam(tree, uLeft, NodeInSubFam);
	SetInFam(tree, uRight, NodeInSubFam);
	}

void AlignSubFam(SeqVect &vAll, const Tree &GuideTree, unsigned uNodeIndex,
  MSA &msaOut)
	{
	const unsigned uSeqCount = vAll.GetSeqCount();

	const char *InTmp = "asf_in.tmp";
	const char *OutTmp = "asf_out.tmp";

	unsigned *Leaves = new unsigned[uSeqCount];
	unsigned uLeafCount;
	GetLeaves(GuideTree, uNodeIndex, Leaves, &uLeafCount);

	SeqVect v;
	for (unsigned i = 0; i < uLeafCount; ++i)
		{
		unsigned uLeafNodeIndex = Leaves[i];
		unsigned uId = GuideTree.GetLeafId(uLeafNodeIndex);
		Seq &s = vAll.GetSeqById(uId);
		v.AppendSeq(s);
		}

#if	TRACE
	{
	Log("Align subfam[node=%d, size=%d] ", uNodeIndex, uLeafCount);
	for (unsigned i = 0; i < uLeafCount; ++i)
		Log(" %s", v.GetSeqName(i));
	Log("\n");
	}
#endif

	TextFile fIn(InTmp, true);

	v.ToFASTAFile(fIn);
	fIn.Close();

	char CmdLine[4096];
	sprintf(CmdLine, "probcons %s > %s 2> /dev/null", InTmp, OutTmp);
//	sprintf(CmdLine, "muscle -in %s -out %s -maxiters 1", InTmp, OutTmp);
	system(CmdLine);

	TextFile fOut(OutTmp);
	msaOut.FromFile(fOut);

	for (unsigned uSeqIndex = 0; uSeqIndex < uLeafCount; ++uSeqIndex)
		{
		const char *Name = msaOut.GetSeqName(uSeqIndex);
		unsigned uId = vAll.GetSeqIdFromName(Name);
		msaOut.SetSeqId(uSeqIndex, uId);
		}

	unlink(InTmp);
	unlink(OutTmp);

	delete[] Leaves;
	}

void ProgAlignSubFams()
	{
	MSA msaOut;

	SetOutputFileName(g_pstrOutFileName.get());
	SetInputFileName(g_pstrInFileName.get());

	SetMaxIters(g_uMaxIters.get());
	SetSeqWeightMethod(g_SeqWeight1.get());

	TextFile fileIn(g_pstrInFileName.get());
	SeqVect v;
	v.FromFASTAFile(fileIn);
	const unsigned uSeqCount = v.Length();

	if (0 == uSeqCount)
		Quit("No sequences in input file");

	ALPHA Alpha = ALPHA_Undefined;
	switch (g_SeqType.get())
		{
	case SEQTYPE_Auto:
		Alpha = v.GuessAlpha();
		break;

	case SEQTYPE_Protein:
		Alpha = ALPHA_Amino;
		break;

	case SEQTYPE_DNA:
		Alpha = ALPHA_DNA;
		break;

	case SEQTYPE_RNA:
		Alpha = ALPHA_RNA;
		break;

	default:
		Quit("Invalid seq type");
		}
	SetAlpha(Alpha);
	v.FixAlpha();

	PTR_SCOREMATRIX UserMatrix = 0;
	if (0 != g_pstrMatrixFileName.get())
		{
		const char *FileName = g_pstrMatrixFileName.get();
		const char *Path = getenv("MUSCLE_MXPATH");
		if (Path != 0)
			{
			size_t n = strlen(Path) + 1 + strlen(FileName) + 1;
			char *NewFileName = new char[n];
			sprintf(NewFileName, "%s/%s", Path, FileName);
			FileName = NewFileName;
			}
		TextFile File(FileName);
		UserMatrix = ReadMx(File);
		g_Alpha = ALPHA_Amino;
		g_PPScore = PPSCORE_SP;
		}

	SetPPScore();

	if (0 != UserMatrix)
		g_ptrScoreMatrix = UserMatrix;

	if (ALPHA_DNA == Alpha || ALPHA_RNA == Alpha)
		{
		SetPPScore(PPSCORE_SPN);
		g_Distance1.get() = DISTANCE_Kmer4_6;
		}

	unsigned uMaxL = 0;
	unsigned uTotL = 0;
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		unsigned L = v.GetSeq(uSeqIndex).Length();
		uTotL += L;
		if (L > uMaxL)
			uMaxL = L;
		}

	SetIter(1);
	g_bDiags.get() = g_bDiags1.get();
	SetSeqStats(uSeqCount, uMaxL, uTotL/uSeqCount);

	SetMuscleSeqVect(v);

	MSA::SetIdCount(uSeqCount);

// Initialize sequence ids.
// From this point on, ids must somehow propogate from here.
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		v.SetSeqId(uSeqIndex, uSeqIndex);

	if (uSeqCount > 1)
		MHackStart(v);

	if (0 == uSeqCount)
		{
		msaOut.Clear();
		return;
		}

	if (1 == uSeqCount && ALPHA_Amino == Alpha)
		{
		const Seq &s = v.GetSeq(0);
		msaOut.FromSeq(s);
		return;
		}

	Tree GuideTree;
	TreeFromSeqVect(v, GuideTree, g_Cluster1.get(), g_Distance1.get(), g_Root1.get());
	SetMuscleTree(GuideTree);

	MSA msa;
	if (g_bLow.get())
		{
		ProgNode *ProgNodes = 0;
		ProgNodes = ProgressiveAlignE(v, GuideTree, msa);
		delete[] ProgNodes;
		}
	else
		ProgressiveAlign(v, GuideTree, msa);
	SetCurrentAlignment(msa);
	TreeFromMSA(msa, GuideTree, g_Cluster2.get(), g_Distance2.get(), g_Root2.get());
	SetMuscleTree(GuideTree);

	unsigned *SubFams = new unsigned[uSeqCount];
	unsigned uSubFamCount;
	SubFam(GuideTree, g_uMaxSubFamCount.get(), SubFams, &uSubFamCount);

	SetProgressDesc("Align node");
	const unsigned uNodeCount = 2*uSeqCount - 1;

	ProgNode *ProgNodes = new ProgNode[uNodeCount];
	bool *NodeIsSubFam = new bool[uNodeCount];
	bool *NodeInSubFam = new bool[uNodeCount];

	for (unsigned i = 0; i < uNodeCount; ++i)
		{
		NodeIsSubFam[i] = false;
		NodeInSubFam[i] = false;
		}

	for (unsigned i = 0; i < uSubFamCount; ++i)
		{
		unsigned uNodeIndex = SubFams[i];
		assert(uNodeIndex < uNodeCount);
		NodeIsSubFam[uNodeIndex] = true;
		SetInFam(GuideTree, uNodeIndex, NodeInSubFam);
		}

	unsigned uJoin = 0;
	unsigned uTreeNodeIndex = GuideTree.FirstDepthFirstNode();
	do
		{
		if (NodeIsSubFam[uTreeNodeIndex])
			{
#if	TRACE
			Log("Node %d: align subfam\n", uTreeNodeIndex);
#endif
			ProgNode &Node = ProgNodes[uTreeNodeIndex];
			AlignSubFam(v, GuideTree, uTreeNodeIndex, Node.m_MSA);
			Node.m_uLength = Node.m_MSA.GetColCount();
			}
		else if (!NodeInSubFam[uTreeNodeIndex])
			{
#if	TRACE
			Log("Node %d: align two subfams\n", uTreeNodeIndex);
#endif
			Progress(uJoin, uSubFamCount - 1);
			++uJoin;

			const unsigned uMergeNodeIndex = uTreeNodeIndex;
			ProgNode &Parent = ProgNodes[uMergeNodeIndex];

			const unsigned uLeft = GuideTree.GetLeft(uTreeNodeIndex);
			const unsigned uRight = GuideTree.GetRight(uTreeNodeIndex);

			ProgNode &Node1 = ProgNodes[uLeft];
			ProgNode &Node2 = ProgNodes[uRight];

			PWPath Path;
			AlignTwoMSAs(Node1.m_MSA, Node2.m_MSA, Parent.m_MSA, Path);
			Parent.m_uLength = Parent.m_MSA.GetColCount();

			Node1.m_MSA.Clear();
			Node2.m_MSA.Clear();
			}
		else
			{
#if	TRACE
			Log("Node %d: in subfam\n", uTreeNodeIndex);
#endif
			;
			}
		uTreeNodeIndex = GuideTree.NextDepthFirstNode(uTreeNodeIndex);
		}
	while (NULL_NEIGHBOR != uTreeNodeIndex);
	ProgressStepsDone();

	unsigned uRootNodeIndex = GuideTree.GetRootNodeIndex();
	ProgNode &RootProgNode = ProgNodes[uRootNodeIndex];

	TextFile fOut(g_pstrOutFileName.get(), true);
	MHackEnd(RootProgNode.m_MSA);
	RootProgNode.m_MSA.ToFile(fOut);

	delete[] NodeInSubFam;
	delete[] NodeIsSubFam;
	delete[] ProgNodes;
	delete[] SubFams;

	ProgNodes = 0;
	NodeInSubFam = 0;
	NodeIsSubFam = 0;
	SubFams = 0;
	}
} 

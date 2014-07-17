#ifndef	MSA_h
#define MSA_h

namespace muscle {

const int MAX_SEQ_NAME = 63;
struct PathEdge;
class TextFile;
class Seq;
class ClusterNode;
class NodeCounts;
class DataBuffer;

class MSA
	{
public:
	MSA();
	virtual ~MSA();

public:
// Ways to create an MSA
	void FromFile(TextFile &File);
	void FromFASTAFile(TextFile &File);
	void FromSeq(const Seq &s);

	void ToFile(TextFile &File) const;
	void ToFASTAFile(TextFile &File) const;
	void ToMSFFile(TextFile &File, const char *ptrComment = 0) const;
	void ToAlnFile(TextFile &File) const;
	void ToHTMLFile(TextFile &File) const;
	void ToPhySequentialFile(TextFile &File) const;
	void ToPhyInterleavedFile(TextFile &File) const;

	void SetSize(unsigned uSeqCount, unsigned uColCount);
	void SetSeqCount(unsigned uSeqCount);
	char GetChar(unsigned uSeqIndex, unsigned uIndex) const;
	unsigned GetLetter(unsigned uSeqIndex, unsigned uIndex) const;
	unsigned GetLetterEx(unsigned uSeqIndex, unsigned uIndex) const;
	const char *GetSeqName(unsigned uSeqIndex) const;
	unsigned GetSeqId(unsigned uSeqIndex) const;
	unsigned GetSeqIndex(unsigned uId) const;
	bool GetSeqIndex(unsigned uId, unsigned *ptruIndex) const;
	double GetOcc(unsigned uColIndex) const;
	void GetFractionalWeightedCounts(unsigned uColIndex, bool bNormalize,
	  FCOUNT fcCounts[], FCOUNT *ptrfcGapStart, FCOUNT *ptrfcGapEnd,
	  FCOUNT *fcGapExtend, FCOUNT *ptrfOcc,
	  FCOUNT *fcLL, FCOUNT *fcLG, FCOUNT *fcGL, FCOUNT *fcGG) const;
	bool IsGap(unsigned uSeqIndex, unsigned uColIndex) const;
	bool IsWildcard(unsigned uSeqIndex, unsigned uColIndex) const;
	bool IsGapColumn(unsigned uColIndex) const;
	bool ColumnHasGap(unsigned uColIndex) const;
	bool IsGapSeq(unsigned uSeqIndex) const;

	void SetChar(unsigned uSeqIndex, unsigned uColIndex, char c);
	void SetSeqName(unsigned uSeqIndex, const char szName[]);
	void SetSeqId(unsigned uSeqIndex, unsigned uId);
	bool HasGap() const;
	bool IsLegalLetter(unsigned uLetter) const;
	void GetSeq(unsigned uSeqIndex, Seq &seq) const;
	void Copy(const MSA &msa);
	double GetCons(unsigned uColIndex) const;
	double GetAvgCons() const;
	double GetPctIdentityPair(unsigned uSeqIndex1, unsigned uSeqIndex2) const;
	bool GetSeqIndex(const char *ptrSeqName, unsigned *ptruSeqIndex) const;
	void DeleteCol(unsigned uColIndex);
	void DeleteColumns(unsigned uColIndex, unsigned uColCount);
	void CopySeq(unsigned uToSeqIndex, const MSA &msaFrom, unsigned uFromSeqIndex);
	void DeleteSeq(unsigned uSeqIndex);
//	void DeleteEmptyCols(bool bProgress = false);
	bool IsEmptyCol(unsigned uColIndex) const;

	WEIGHT GetSeqWeight(unsigned uSeqIndex) const;
	WEIGHT GetTotalSeqWeight() const;
	void SetSeqWeight(unsigned uSeqIndex, WEIGHT w) const;
	void NormalizeWeights(WEIGHT wTotal) const;
	bool WeightsSet() const;

	unsigned GetGCGCheckSum(unsigned uSeqIndex) const;

	ALPHA GuessAlpha() const;
	void FixAlpha();

	unsigned UniqueResidueTypes(unsigned uColIndex) const;

	void UnWeight();

	void GetNodeCounts(unsigned uAlignedColIndex, NodeCounts &Counts) const;
	void ValidateBreakMatrices() const;
	unsigned GetCharCount(unsigned uSeqIndex, unsigned uColIndex) const;
	const char *GetSeqBuffer(unsigned uSeqIndex) const;
	unsigned AlignedColIndexToColIndex(unsigned uAlignedColIndex) const;
	unsigned GetSeqLength(unsigned uSeqIndex) const;
	void GetPWID(unsigned uSeqIndex1, unsigned uSeqIndex2, double *ptrdPWID,
	  unsigned *ptruPosCount) const;

	void GetPairMap(unsigned uSeqIndex1, unsigned uSeqIndex2, int iMap1[],
	  int iMap2[]) const;

	void LogMe() const;
	void ListWeights() const;

	void GapInfoToDataBuffer(DataBuffer &Buffer) const;
	void GapInfoFromDataBuffer(const DataBuffer &Buffer);
	double GetPctGroupIdentityPair(unsigned uSeqIndex1, unsigned uSeqIndex2) const;

	void Clear()
		{
		Free();
		}
	unsigned GetSeqCount() const
		{
		return m_uSeqCount;
		}
	unsigned GetColCount() const
		{
		return m_uColCount;
		}

	static bool SeqsEq(const MSA &a1, unsigned uSeqIndex1, const MSA &a2,
	  unsigned uSeqIndex2);

	static void SetIdCount(unsigned uIdCount);

private:
	friend void SetMSAWeightsMuscle(MSA &msa);
	friend void SetThreeWayWeightsMuscle(MSA &msa);
	void SetHenikoffWeightsPB() const;
	void SetHenikoffWeights() const;
	void SetGSCWeights() const;
	void SetUniformWeights() const;
	void SetClustalWWeights(const Tree &tree);

	void Free();
	void AppendSeq(char *ptrSeq, unsigned uSeqLength, char *ptrLabel);
	void ExpandCache(unsigned uSeqCount, unsigned uColCount);
	void CalcWeights() const;
	void GetNameFromFASTAAnnotationLine(const char szLine[],
	  char szName[], unsigned uBytes);
	void CopyCol(unsigned uFromCol, unsigned uToCol);
	unsigned CalcBLOSUMWeights(ClusterTree &BlosumCluster) const;
	void SetBLOSUMSubtreeWeight(const ClusterNode *ptrNode, double dWeight) const;
	unsigned SetBLOSUMNodeWeight(const ClusterNode *ptrNode, double dMinDist) const;
	void SetSubtreeWeight2(const ClusterNode *ptrNode) const;
	void SetSubtreeGSCWeight(ClusterNode *ptrNode) const;

	void CalcHenikoffWeightsColPB(unsigned uColIndex) const;
	void CalcHenikoffWeightsCol(unsigned uColIndex) const;

private:
	unsigned m_uSeqCount;
	unsigned m_uColCount;
	unsigned m_uCacheSeqLength;
	unsigned m_uCacheSeqCount;
	char **m_szSeqs;
	char **m_szNames;

	static TLS<unsigned> m_uIdCount;

	unsigned *m_IdToSeqIndex;
	unsigned *m_SeqIndexToId;

	WEIGHT *m_Weights;
	};

void SeqVectFromMSA(const MSA &msa, SeqVect &v);
void DeleteGappedCols(MSA &msa);
void MSAFromColRange(const MSA &msaIn, unsigned uFromColIndex, unsigned uColCount,
  MSA &msaOut);
void MSACat(const MSA &msa1, const MSA &msa2, MSA &msaCat);
void MSAAppend(MSA &msa1, const MSA &msa2);
void MSAFromSeqSubset(const MSA &msaIn, const unsigned uSeqIndexes[], unsigned uSeqCount,
  MSA &msaOut);
void AssertMSAEq(const MSA &msa1, const MSA &msa2);
void AssertMSAEqIgnoreCaseAndGaps(const MSA &msa1, const MSA &msa2);
void MSASubsetByIds(const MSA &msaIn, const unsigned Ids[], unsigned uIdCount,
  MSA &msaOut);
void SetMSAWeightsMuscle(MSA &msa);
void SetClustalWWeightsMuscle(MSA &msa);
void SetThreeWayWeightsMuscle(MSA &msa);

} // namespace muscle

#endif	// MSA_h

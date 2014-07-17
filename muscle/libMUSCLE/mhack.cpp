#include "libMUSCLE/muscle.h"
#include "libMUSCLE/seqvect.h"
#include "libMUSCLE/msa.h"
#include "libMUSCLE/threadstorage.h"

namespace muscle {

/***
Methionine hack.
Most proteins start with M.
This results in odd-looking alignments with the terminal Ms aligned followed
immediately by gaps.
Hack this by treating terminal M like X.
***/

static TLS<bool *> M;

void MHackStart(SeqVect &v)
	{
	if (ALPHA_Amino != g_Alpha.get())
		return;

	const unsigned uSeqCount = v.Length();
	M.get() = new bool[uSeqCount];
	memset(M.get(), 0, uSeqCount*sizeof(bool));
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		Seq &s = v.GetSeq(uSeqIndex);
		if (0 == s.Length())
			continue;
		unsigned uId = s.GetId();
		if (s[0] == 'M' || s[0] == 'm')
			{
			M.get()[uId] = true;
			s[0] = 'X';
			}
		}
	}

void MHackEnd(MSA &msa)
	{
	if (ALPHA_Amino != g_Alpha.get())
		return;
	if (0 == M.get())
		return;

	const unsigned uSeqCount = msa.GetSeqCount();
	const unsigned uColCount = msa.GetColCount();
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		unsigned uId = msa.GetSeqId(uSeqIndex);
		if (M.get()[uId])
			{
			for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
				{
				if (!msa.IsGap(uSeqIndex, uColIndex))
					{
					msa.SetChar(uSeqIndex, uColIndex, 'M');
					break;
					}
				}
			}
		}

	delete[] M.get();
	M.get() = 0;
	}
} 

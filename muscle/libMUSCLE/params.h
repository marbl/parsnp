#ifndef params_h
#define params_h
#include "libMUSCLE/threadstorage.h"

namespace muscle {

extern TLS<const char *> g_pstrInFileName;
extern TLS<const char *> g_pstrOutFileName;

extern TLS<const char *> g_pstrFASTAOutFileName;
extern TLS<const char *> g_pstrMSFOutFileName;
extern TLS<const char *> g_pstrClwOutFileName;
extern TLS<const char *> g_pstrClwStrictOutFileName;
extern TLS<const char *> g_pstrHTMLOutFileName;
extern TLS<const char *> g_pstrPHYIOutFileName;
extern TLS<const char *> g_pstrPHYSOutFileName;
extern TLS<const char *> g_pstrDistMxFileName1;
extern TLS<const char *> g_pstrDistMxFileName2;

extern TLS<const char *> g_pstrFileName1;
extern TLS<const char *> g_pstrFileName2;

extern TLS<const char *> g_pstrSPFileName;
extern TLS<const char *> g_pstrMatrixFileName;

extern TLS<const char *> g_pstrUseTreeFileName;
extern TLS<bool> g_bUseTreeNoWarn;

extern TLS<const char *> g_pstrComputeWeightsFileName;
extern TLS<const char *> g_pstrScoreFileName;

extern TLS<SCORE> g_scoreGapOpen;
extern TLS<SCORE> g_scoreCenter;
extern TLS<SCORE> g_scoreGapExtend;
extern TLS<SCORE> g_scoreGapAmbig;

#if	DOUBLE_AFFINE
extern TLS<SCORE> g_scoreGapOpen2;
extern TLS<SCORE> g_scoreGapExtend2;
#endif

extern TLS<unsigned> g_uSmoothWindowLength;
extern TLS<unsigned> g_uAnchorSpacing;
extern TLS<unsigned> g_uMaxTreeRefineIters;

extern TLS<unsigned> g_uMinDiagLength;
extern TLS<unsigned> g_uMaxDiagBreak;
extern TLS<unsigned> g_uDiagMargin;

extern TLS<unsigned> g_uRefineWindow;
extern TLS<unsigned> g_uWindowFrom;
extern TLS<unsigned> g_uWindowTo;
extern TLS<unsigned> g_uSaveWindow;
extern TLS<unsigned> g_uWindowOffset;

extern TLS<unsigned> g_uMaxSubFamCount;

extern TLS<unsigned> g_uHydrophobicRunLength;
extern TLS<float> g_dHydroFactor;

extern TLS<float> g_dSmoothScoreCeil;
extern TLS<float> g_dMinBestColScore;
extern TLS<float> g_dMinSmoothScore;
extern TLS<float> g_dSUEFF;

extern TLS<bool> g_bPrecompiledCenter;
extern TLS<bool> g_bNormalizeCounts;
extern TLS<bool> g_bDiags1;
extern TLS<bool> g_bDiags2;
extern TLS<bool> g_bDiags;
extern TLS<bool> g_bAnchors;
extern TLS<bool> g_bCatchExceptions;

extern TLS<bool> g_bMSF;
extern TLS<bool> g_bAln;
extern TLS<bool> g_bClwStrict;
extern TLS<bool> g_bHTML;
extern TLS<bool> g_bPHYI;
extern TLS<bool> g_bPHYS;

extern TLS<bool> g_bQuiet;
extern TLS<bool> g_bVerbose;
extern TLS<bool> g_bRefine;
extern TLS<bool> g_bRefineW;
extern TLS<bool> g_bRefineX;
extern TLS<bool> g_bLow;
extern TLS<bool> g_bSW;
extern TLS<bool> g_bClusterOnly;
extern TLS<bool> g_bProfile;
extern TLS<bool> g_bProfileOnStdIn;
extern TLS<bool> g_bAnchoredPP;
extern TLS<bool> g_bProfDB;
extern TLS<bool> g_bPPScore;
extern TLS<bool> g_bBrenner;
extern TLS<bool> g_bDimer;
extern TLS<bool> g_bVersion;
extern TLS<bool> g_bStable;
extern TLS<bool> g_bFASTA;
extern TLS<bool> g_bPAS;
extern TLS<bool> g_bTomHydro;
extern TLS<bool> g_bMakeTree;

extern TLS<PPSCORE> g_PPScore;
extern TLS<OBJSCORE> g_ObjScore;

extern TLS<DISTANCE> g_Distance1;
extern TLS<CLUSTER> g_Cluster1;
extern TLS<ROOT> g_Root1;
extern TLS<SEQWEIGHT> g_SeqWeight1;

extern TLS<DISTANCE> g_Distance2;
extern TLS<CLUSTER> g_Cluster2;
extern TLS<ROOT> g_Root2;
extern TLS<SEQWEIGHT> g_SeqWeight2;

extern TLS<unsigned> g_uMaxIters;
extern TLS<unsigned long> g_ulMaxSecs;
extern TLS<unsigned> g_uMaxMB;

extern TLS<SEQTYPE> g_SeqType;
extern TLS<TERMGAPS> g_TermGaps;

} // namespace muscle

#endif // params_h

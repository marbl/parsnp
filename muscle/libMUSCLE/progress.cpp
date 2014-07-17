#include "libMUSCLE/muscle.h"
#include <stdio.h>
#include <time.h>

namespace muscle {

// Functions that provide visible feedback to the user
// that progress is being made.
#ifdef _OPENMP
#ifdef __GNUG__
#warning "progress indicators have not been made parallel-safe yet"
#endif
#endif

static unsigned g_uIter = 0;		// Main MUSCLE iteration 1, 2..
static unsigned g_uLocalMaxIters = 0;	// Max iters
static FILE *g_fProgress = stderr;	// Default to standard error
static char g_strFileName[32];		// File name
static time_t g_tLocalStart;				// Start time
static char g_strDesc[32];			// Description
static bool g_bWipeDesc = false;
static int g_nPrevDescLength;
static unsigned g_uTotalSteps;

double GetCheckMemUseMB()
	{
	unsigned MB = (unsigned) GetMemUseMB();
	if (0 == g_uMaxMB.get() || MB <= g_uMaxMB.get())
		return MB;
	fprintf(stderr, "\n\n*** MAX MEMORY %u MB EXCEEDED***\n", g_uMaxMB.get());
	fprintf(stderr, "Memory allocated so far %u MB, physical RAM %u MB\n",
	  MB, (unsigned) GetRAMSizeMB());
	fprintf(stderr, "Use -maxmb <n> option to increase limit, where <n> is in MB.\n");
	SaveCurrentAlignment();
	exit(EXIT_FatalError);
	return MB;
	}

const char *ElapsedTimeAsStr()
	{
	time_t Now = time(0);
	unsigned long ElapsedSecs = (unsigned long) (Now - g_tLocalStart);
	return SecsToStr(ElapsedSecs);
	}

const char *MemToStr(double MB)
	{
	if (MB < 0)
		return "";

	static char Str[11];
	static double MaxMB = 0;
	static double RAMMB = 0;

	if (RAMMB == 0)
		RAMMB = GetRAMSizeMB();

	if (MB > MaxMB)
		MaxMB = MB;
	double Pct = (MaxMB*100.0)/RAMMB;
	if (Pct > 100)
		Pct = 100;
	sprintf(Str, "%.0f MB(%.0f%%)", MaxMB, Pct);
	return Str;
	}

void SetInputFileName(const char *pstrFileName)
	{
	NameFromPath(pstrFileName, g_strFileName, sizeof(g_strFileName));
	}

void SetSeqStats(unsigned uSeqCount, unsigned uMaxL, unsigned uAvgL)
	{
	if (g_bQuiet.get())
		return;

	fprintf(g_fProgress, "%s %u seqs, max length %u, avg  length %u\n",
	  g_strFileName, uSeqCount, uMaxL, uAvgL);
	if (g_bVerbose.get())
		Log("%u seqs, max length %u, avg  length %u\n",
		  uSeqCount, uMaxL, uAvgL);
	}

void SetStartTime()
	{
	time(&g_tLocalStart);
	}

unsigned long GetStartTime()
	{
	return (unsigned long) g_tLocalStart;
	}

void SetIter(unsigned uIter)
	{
	g_uIter = uIter;
	}

void IncIter()
	{
	++g_uIter;
	}

void SetMaxIters(unsigned uMaxIters)
	{
	g_uLocalMaxIters = uMaxIters;
	}

void SetProgressDesc(const char szDesc[])
	{
	strncpy(g_strDesc, szDesc, sizeof(g_strDesc));
	g_strDesc[sizeof(g_strDesc) - 1] = 0;
	}

static void Wipe(int n)
	{
	for (int i = 0; i < n; ++i)
		fprintf(g_fProgress, " ");
	}

void Progress(const char *szFormat, ...)
	{
	CheckMaxTime();

	if (g_bQuiet.get())
		return;

	double MB = GetCheckMemUseMB();

	char szStr[4096];
	va_list ArgList;
	va_start(ArgList, szFormat);
	vsprintf(szStr, szFormat, ArgList);

	fprintf(g_fProgress, "%8.8s  %12s  %s",
	  ElapsedTimeAsStr(),
	  MemToStr(MB),
	  szStr);

	fprintf(g_fProgress, "\n");
	fflush(g_fProgress);
	}

void Progress(unsigned uStep, unsigned uTotalSteps)
	{
	CheckMaxTime();

	if (g_bQuiet.get())
		return;

	double dPct = ((uStep + 1)*100.0)/uTotalSteps;
	double MB = GetCheckMemUseMB();
	fprintf(g_fProgress, "%8.8s  %12s  Iter %3u  %6.2f%%  %s",
	  ElapsedTimeAsStr(),
	  MemToStr(MB),
	  g_uIter,
	  dPct,
	  g_strDesc);

	if (g_bWipeDesc)
		{
		int n = g_nPrevDescLength - (int) strlen(g_strDesc);
		Wipe(n);
		g_bWipeDesc = false;
		}

	fprintf(g_fProgress, "\r");

	g_uTotalSteps = uTotalSteps;
	}

void ProgressStepsDone()
	{
	CheckMaxTime();

	if (g_bVerbose.get())
		{
		double MB = GetCheckMemUseMB();
		Log("Elapsed time %8.8s  Peak memory use %12s  Iteration %3u %s\n",
		 ElapsedTimeAsStr(),
		 MemToStr(MB),
		 g_uIter,
		 g_strDesc);
		}

	if (g_bQuiet.get())
		return;

	Progress(g_uTotalSteps - 1, g_uTotalSteps);
	fprintf(g_fProgress, "\n");
	g_bWipeDesc = true;
	g_nPrevDescLength = (int) strlen(g_strDesc);
	}
} 

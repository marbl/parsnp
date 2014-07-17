#include "libMUSCLE/muscle.h"
#include "libMUSCLE/msa.h"
#include "libMUSCLE/textfile.h"
#include <time.h>
#include "libMUSCLE/threadstorage.h"

namespace muscle {


static TLS<const char *> pstrOutputFileName;

void SetOutputFileName(const char *out)
	{
	pstrOutputFileName.get() = out;
	}

void SetCurrentAlignment(MSA &msa)
	{
    extern TLS<MSA *>ptrBestMSA;
	ptrBestMSA.get() = &msa;
	}

void SaveCurrentAlignment()
	{
    extern TLS<MSA *>ptrBestMSA;
	static TLS<bool> bCalled(false);
	if (bCalled.get())
		{
		fprintf(stderr,
		  "\nRecursive call to SaveCurrentAlignment, giving up attempt to save.\n");
		exit(EXIT_FatalError);
		}

	if (0 == ptrBestMSA.get())
		{
		fprintf(stderr, "\nAlignment not completed, cannot save.\n");
		Log("Alignment not completed, cannot save.\n");
		exit(EXIT_FatalError);
		}

	if (0 == pstrOutputFileName.get())
		{
		fprintf(stderr, "\nOutput file name not specified, cannot save.\n");
		exit(EXIT_FatalError);
		}

	fprintf(stderr, "\nSaving current alignment ...\n");

	TextFile fileOut(pstrOutputFileName.get(), true);
	ptrBestMSA.get()->ToFASTAFile(fileOut);

	fprintf(stderr, "Current alignment saved to \"%s\".\n", pstrOutputFileName.get());
	Log("Current alignment saved to \"%s\".\n", pstrOutputFileName.get());
	}

void CheckMaxTime()
	{
	if (0 == g_ulMaxSecs.get())
		return;

	time_t Now = time(0);
	time_t ElapsedSecs = Now - GetStartTime();
	if (ElapsedSecs <= (time_t) g_ulMaxSecs.get())
		return;

	Log("Max time %s exceeded, elapsed seconds = %ul\n",
	  MaxSecsToStr(), ElapsedSecs);

	SaveCurrentAlignment();
	exit(EXIT_Success);
	}
} 

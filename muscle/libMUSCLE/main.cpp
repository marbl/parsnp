//@@TODO reconcile /muscle with /muscle3.6

#include "libMUSCLE/muscle.h"
#include <stdio.h>
#ifdef	WIN32
#include <windows.h>	// for SetPriorityClass()
#include <io.h>			// for isatty()
#else
#include <unistd.h>		// for isatty()
#endif

using namespace muscle;

int main(int argc, char **argv)
	{
#if	WIN32
// Multi-tasking does not work well in CPU-bound
// console apps running under Win32.
// Reducing the process priority allows GUI apps
// to run responsively in parallel.
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	g_argc.get() = argc;
	g_argv.get() = argv;

	SetNewHandler();
	SetStartTime();
	ProcessArgVect(argc - 1, argv + 1);
	SetParams();
	SetLogFile();

	//extern void TestSubFams(const char *);
	//TestSubFams(g_pstrInFileName.get());
	//return 0;

	if (g_bVersion.get())
		{
		printf(MUSCLE_LONG_VERSION "\n");
		exit(EXIT_SUCCESS);
		}

	if (!g_bQuiet.get())
		Credits();

	if (MissingCommand() && isatty(0))
		{
		Usage();
		exit(EXIT_SUCCESS);
		}

	if (g_bCatchExceptions.get())
		{
		try
			{
			Run();
			}
		catch (...)
			{
			OnException();
			exit(EXIT_Except);
			}
		}
	else
		Run();

	exit(EXIT_Success);
	}


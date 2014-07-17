#include "libMUSCLE/muscle.h"
#include <stdio.h>

namespace muscle {

static const char szOnExceptionMessage[] =
	{
	"\nFatal error, exception caught.\n"
	};

void OnException()
	{
	fprintf(stderr, szOnExceptionMessage);
	Log(szOnExceptionMessage);
	Log("Finished %s\n", GetTimeAsStr());
	exit(EXIT_Except);
	}
} 

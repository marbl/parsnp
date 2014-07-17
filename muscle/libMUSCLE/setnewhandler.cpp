#include "libMUSCLE/muscle.h"
#include <stdio.h>
#include <new>

namespace muscle {

const int ONE_MB = 1024*1024;
const size_t RESERVE_BYTES = 8*ONE_MB;
static void *EmergencyReserve = 0;

void OnOutOfMemory()
	{
#pragma omp critical 
{
	free(EmergencyReserve);
	fprintf(stderr, "\n*** OUT OF MEMORY ***\n");
	fprintf(stderr, "Memory allocated so far %g MB\n", GetMemUseMB());
    extern TLS<MSA *>ptrBestMSA;
	if (ptrBestMSA.get() == 0)
		fprintf(stderr, "No alignment generated\n");
	else
	SaveCurrentAlignment();
	exit(EXIT_FatalError);
}
	}

void SetNewHandler()
	{
#pragma omp critical 
{
	if(EmergencyReserve == 0)
	{
		EmergencyReserve = malloc(RESERVE_BYTES);
		std::set_new_handler(OnOutOfMemory);
	}
}
	}
} 

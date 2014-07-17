#if	WIN32
#include <intrin.h>

namespace muscle {

typedef unsigned __int64 TICKS;

#pragma warning(disable:4035)
inline TICKS GetClockTicks()
	{
		// aed 16/7/7: get ticks in a way that's portable to x64
		return __rdtsc();
	}

#define	StartTimer()	__int64 t1__ = GetClockTicks()

#define	GetElapsedTicks()	(GetClockTicks() - t1__)

static double TicksToSecs(TICKS t)
	{
	return (__int64) t/2.5e9;
	}

} // namespace muscle

#endif	// WIN32

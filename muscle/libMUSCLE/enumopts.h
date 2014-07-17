#ifndef enumopts_h
#define enumopts_h

namespace muscle {

struct EnumOpt
	{
	const char *pstrOpt;
	int iValue;
	};

#define	s(t)		extern EnumOpt t##_Opts[];
#define c(t, x)		/* empty */
#define e(t)		/* empty */
#include "libMUSCLE/enums.h"	

} // namespace muscle

#endif // enumopts_h

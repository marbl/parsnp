#include "libMUSCLE/muscle.h"
#include "libMUSCLE/enumopts.h"

namespace muscle {

#define	s(t)		EnumOpt t##_Opts[] = {
#define c(t, x)		#x, t##_##x,
#define e(t)		0, 0 };

#include "libMUSCLE/enums.h"
} 

#ifndef Gonnet_h
#define Gonnet_h

namespace muscle {

typedef double t_ROW[20];

const t_ROW *GetGonnetMatrix(unsigned N);
SCORE GetGonnetGapOpen(unsigned N);
SCORE GetGonnetGapExtend(unsigned N);

extern double GonnetLookup[400][400];

} // namespace muscle

#endif	// Gonnet_h

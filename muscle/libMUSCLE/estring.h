#ifndef pathsum_h
#define pathsum_h

namespace muscle {

void PathToEstrings(const PWPath &Path, short **ptresA, short **ptresB);
void EstringsToPath(const short esA[], const short esB[], PWPath &Path);
void MulEstrings(const short es1[], const short es2[], short esp[]);
void EstringOp(const short es[], const Seq &sIn, Seq &sOut);
unsigned EstringOp(const short es[], const Seq &sIn, MSA &a);
void LogEstring(const short es[]);
unsigned LengthEstring(const short es[]);
short *EstringNewCopy(const short es[]);

} // namespace muscle

#endif	// pathsum_h

//Fork of libMUSCLE code from http://sourceforge.net/p/mauve/code/HEAD/tree/muscle/

#ifndef _MuscleInterface_h_
#define _MuscleInterface_h_
#include <iostream>
#include <fstream>
#include <bitset>
#include <string>
#include <vector>
#include <map>
#include <omp.h>
using namespace std;

class MuscleInterface
{
public:
	static MuscleInterface& getMuscleInterface();
	MuscleInterface();
	bool CallMuscleFast( vector< string > & aln_matrix, const vector< string >& seq_table);
};

#endif // _MuscleInterface_h_

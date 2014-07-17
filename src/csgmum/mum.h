/////////////////////////////////////////
// mum.h
// MUM suffix tree search code
/////////////////////////////////////////

#ifndef _MMUM_H_
#define _MMUM_H_

#include "csg.h"

typedef struct UM UM;
struct UM {
 	int UP;
	int EP;
};

typedef struct Mum Mum;

typedef struct SP SP;
struct SP 
{
  ulong *MSP;
  char  *forward; 
};
struct Mum {
	ulong *DSP;
	int LON;
    char *forward;
};

typedef struct IRegion IRegion;
struct IRegion {
	ulong len_region;
	ulong ini_region;
	char *sequence;
	char *rc;
};

void Search_MUM(
		ulong *num_mums, 
		Mum *list_mums, 
		UM *Master, 
		int size, 
		int min_lenght_of_mum, 
		IRegion *regions, 
		int number_of_regions, 
		ulong **SP);

void Find_UM(CSG *csg, const char *s,ulong *SP, UM *Pair);

void Intersect_UM(CSG *csg, UM *Master, UM *Pair, int size, ulong *SP);

void Merge_Pair(UM *Pair, UM *PairRC, int size, int length, SP *SPF, SP *SPR, int pos);
void Merge_Master( UM *Master, UM *MasterRC, int size, int length, SP *SPF, SP *SPR, ulong * MSP, int pos);

int
Find_MUM(
	ulong *num_mums
	, Mum **list_mums
	, int number_of_regions	/* the number of sequences to align */
	, IRegion *regions
	, int min_lenght_of_mum
	, float factor
	, CSG **csg
	, UM ** Master
	, UM ** Pair
);	// return the status

#endif


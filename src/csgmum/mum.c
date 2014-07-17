/////////////////////////////////////////
// mum.c
// MUM suffix tree search code 
/////////////////////////////////////////

//#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include  "types.h"

#include "mum.h"
#include "fastaread.h"
#include "csg.h"


#define Min(a, b) ((a) <= (b) ? (a) : (b))
#define Max(a, b) ((a) >= (b) ? (a) : (b))


/*
 ////////////////////////////////////////////////////////////////////
 // MUM search
 //////////////////////////////////////////////////////////////////////
 */

void
Test_UM(UM *Pair, int k, int l, int MMP, ulong *SP, int pos_label) {
    if (Pair[l].EP == 0) {
        Pair[l].EP = MMP;
        Pair[l].UP = pos_label;
        SP[l] = k - MMP + l;
        return;
    }
    
    if (MMP > Pair[l].UP) {
        if (MMP > Pair[l].EP) {
            Pair[l].UP = Pair[l].EP;
            Pair[l].EP = MMP;
            SP[l] = k - MMP + l;
        } else {
            Pair[l].UP = MMP;
        }
    }
}

void print_state(UM *Pair, int size)
{
    int i;
    printf("    ");
    for (i = 0; i < size; i++) {
        printf("%.2d ", i);
    }
    printf("\nUP  ");
    for (i = 0; i < size; i++) {
        printf("%.2d ", Pair[i].UP);
    }
    printf("\nEP  ");
    for (i = 0; i < size; i++) {
        printf("%.2d ", Pair[i].EP);
    }
    printf("\n");
}
void Merge_Pair(UM *Pair, UM *PairRC, int size, int length,SP *SPF, SP *SPR, int pos )
{
    int i = 0;
	int k = 0;
    
    while (k < size)
    {
        
        
        if ( Pair[k].EP > PairRC[k].EP )
        {
            
            SPF[pos].forward[k] = 1;
            SPR[pos].forward[k] = 0;
        }
        else
        {
            //both false, ties go to forward for now, grab reverse
            Pair[k].EP = PairRC[k].EP;
            Pair[k].UP = PairRC[k].UP;
            SPF[pos].MSP[k] = SPR[pos].MSP[k];
            SPF[pos].forward[k] = 0;
            SPR[pos].forward[k] = 1;
        }
        k++;
	}
    
}
void Merge_Master(UM *Master, UM *MasterRC, int size, int length, SP *SPF, SP *SPR, ulong* MSP, int pos )
{
	int i = 0;
	int k = 0;
    
    
    while (k < size)
    {
        
        if ( Master[k].EP > MasterRC[k].EP )
        {
            MasterRC[k].EP = Master[k].EP;
            MasterRC[k].UP = Master[k].UP;
            MSP[k] = SPF[pos].MSP[k];
            SPF[pos].forward[k] = 1;
        }
        else
        {
            //both false, ties go to forward for now, grab reverse
            Master[k].EP = MasterRC[k].EP;
            //need to update MergedMaster to refect RC position?
            Master[k].UP = MasterRC[k].UP;

            // SP-------UP-------EP>  <EP-------UP-------SP      <SP------UP--------EP
            SPF[pos].MSP[k] = MSP[k];
            SPF[pos].forward[k] = 0;
        }

        k++;
    }

}

void Intersect_UM(CSG *csg, UM *Master, UM *Pair, int size, ulong *SP)
{
    int k = 1;
    int MMP;
    
    Master[0].UP = Max(Master[0].UP, Pair[0].UP);
    //if doesn't exist, EP gets set to 0?
    Master[0].EP = Min(Master[0].EP, Pair[0].EP);
    
    while (k < size) {
        if (Pair[k - 1].UP > Pair[k].UP) {
            Pair[k].UP = Pair[k - 1].UP;
            if (0 && Pair[k].EP == 0)
            {
                Pair[k].EP = Pair[k].UP;
                SP[k] = 0;
            }
            else if (Pair[k].EP < Pair[k].UP) {
                Pair[k].EP = Pair[k].UP;
                SP[k] = SP[k-1] + 1;
            }
        }

        MMP = Pair[k - 1].EP;
        if (MMP > Pair[k].UP ){
            if (MMP >= Pair[k].EP) {
                Pair[k].UP = Pair[k].EP;
                Pair[k].EP = MMP;
                SP[k] = SP[k-1] + 1;
            } else {
                Pair[k].UP = MMP;
            }
        }
        //update UP,EP positions
        
        if (Pair[k].EP != 0 || 1)
        {
            Master[k].UP = Max(Master[k].UP, Pair[k].UP);
            Master[k].EP = Min(Master[k].EP, Pair[k].EP);
        }
        else
        {
            //NO UM/MUM in this genome, conserve previous end position
            Master[k].EP = Master[k].EP;
            Master[k].UP = Max(Master[k].UP, Pair[k].UP);
        }
        Pair[k - 1].EP = Pair[k - 1].UP = 0;
        k++;
    }
    Pair[k - 1].EP = Pair[k - 1].UP = 0;
}

void Find_UM(CSG *csg, const char *s,ulong *SP, UM *Pair)
{
    char *pendent = (char *)s;
    int  lprefix  = 0;
    CSG_Node *state;
    int  lg = 0;
    
    int pos_seq_pivot;
    int pos_seq_act;
    int MMP;
    int num_leaf;
    CSG_Node *next;
    int len_to_next;
    int pos_label;
    int len;
    
    state = NEXT(csg, INITIAL, *pendent);
	while (state == NINDEF)
	{
		pendent++;
        state = NEXT(csg, INITIAL, *pendent);
	}
    len   = LEN_L(state);
    
    while ( *pendent != LAST_SYMBOL )
	{
		if (!lprefix)
		{
			lprefix = (state == INITIAL ? 1 : prefix_comu(csg,  state, (char *) pendent));
		}
        next = NEXT(csg, state, *(pendent + len));
        if ((lprefix > len && !IS_LEAF(state)) || (lprefix == len && next != NINDEF))
		{
            lprefix -= len;
            pendent += len;
            lg      += len;
            
		}
		else
		{
            if (IS_LEAF(state))
			{
                pos_seq_pivot = LABEL(state) - csg->seq;
                pos_seq_act   = pendent - s + lprefix;
                MMP           = pos_seq_pivot + lprefix;
                num_leaf      = pos_seq_pivot - lg;
                pos_label     = MMP - lprefix;
                Test_UM(Pair, pos_seq_act, num_leaf, MMP, SP, pos_label);
                len_to_next = LABEL(NEXT_LEAF(state)) - LABEL(state) + 1;
                while (lprefix >= len_to_next) 
				{
                	state    = NEXT_LEAF(state);
                    
                	pendent += len_to_next - 1;
                	lprefix -= len_to_next - 1;
                	len_to_next = LABEL(NEXT_LEAF(state)) - LABEL(state) + 1;
                } 
			}
            
			next = CLINK(state);
            
			lg   = next == INITIAL ? -1 : LG(next);
		}
        
        if (state == next)
		{
			pendent++;
		}
		state = next;

		len   = LEN_L(state);
        
	}
}



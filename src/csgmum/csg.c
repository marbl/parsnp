/////////////////////////////////////////
// csg.c
// compressed suffix graph structure
/////////////////////////////////////////
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>

#include "csg.h"
#include "fastaread.h"

#define CLONE_EDGES(original, clon) \
do { \
	if (NEXT(csg, nextState, 71) != NINDEF) \
		INS_EDGE(csg, clon, NEXT(csg, nextState, 71)); \
	if (NEXT(csg, nextState, 67) != NINDEF) \
		INS_EDGE(csg, clon, NEXT(csg, nextState, 67)); \
	if (NEXT(csg, nextState, 84) != NINDEF) \
		INS_EDGE(csg, clon, NEXT(csg, nextState, 84)); \
	if (NEXT(csg, nextState, 65) != NINDEF) \
		INS_EDGE(csg, clon, NEXT(csg, nextState, 65)); \
	if (NEXT(csg, nextState, 78) != NINDEF) \
		INS_EDGE(csg, clon, NEXT(csg, nextState, 78)); \
} while (0)

#define CLEAR_EDGES(source) \
do { \
	if (OUT(source) > 1) \
	        free(source->arc); \
	OUT(source) = 0; \
} while (0)


void SWAP_EDGE(CSG *csg, CSG_Node *source, CSG_Node *old, CSG_Node *nnew)
{
	int nout = OUT(source);
	
	if (nout == 1) {
		source->arc = (CSG_Node **)nnew;
		nnew->lg = old->lg;
	} else while (nout--) {
		if (old == source->arc[nout]) {
			source->arc[nout] = nnew;
			nnew->lg = old->lg;
			return;
		}
	}

}

CSG_Node *NEXT(CSG *csg, CSG_Node *source, int target) 
{
	int nout = OUT(source);
        //printf("target %i\n",target);
	if (nout == 1 && target == *LABEL((CSG_Node *)source->arc))
		return (CSG_Node *)source->arc;
	else if (nout == 1)
		return NINDEF;
	while (nout--) {
		if (target == *LABEL(source->arc[nout]))
			return source->arc[nout];
	}
	return NINDEF;	
}

void INS_EDGE(CSG *csg, CSG_Node *source, CSG_Node *target) 
{ 
	CSG_Node **bucket;

	if (*(LABEL(target)) == LAST_SYMBOL) 
		SET_FINAL(source); 
	else {
		if (OUT(source) == 0) {
			source->arc = (CSG_Node **)target;
		} else  {
			if (OUT(source) == 1) {
				bucket = (CSG_Node **)malloc(sizeof(CSG_Node *)
				* 2);
				bucket[0] = (CSG_Node *)source->arc;
				source->arc = bucket;
			} else if (OUT(source) == 2) {
				source->arc = (CSG_Node **)realloc(source->arc 
				, sizeof(CSG_Node *) * 5);
			} 
			source->arc[OUT(source)] = target; 
		}
		OUT(source)++; 
		LG(target) = MAX(LG(target), (IS_INITIAL(source) 
		? 0 :  LG(source) + LEN_L(source))); 
	}
}; 


int GET_SOLID(CSG *csg, CSG_Node *source, CSG_Node *node_t) 
{ 
	int res = (LG(node_t) == (IS_INITIAL(source) ? 0 
	: LG(source) + LEN_L(source) ));

	return res;
}

int emsurto = 0;

CSG *new_CSG(CSG *csg_clone, ulong node_ini, const char *seq, ulong size, int setmem)
{
    CSG *csg;

    if(setmem)
	  csg = csg_clone;
	else
	  csg  = (CSG *) malloc(sizeof(CSG));
	
	csg->seq       = seq;
	csg->size      = size;
  	csg->ini_node  = size;

	if (setmem)
    {
	  csg->node = (CSG_Node *) memset(csg->node,0, sizeof(CSG_Node)*node_ini);

	}	
	else 
    {
	  csg->node = (CSG_Node *) calloc(node_ini,sizeof(CSG_Node));

	}	
	
	CLINK(INITIAL) = INITIAL;
	LEN_L(INITIAL) = 1;
  	csg->last_state = 0;
  	csg->num_nodes  = 1;
  	csg->num_leafs  = 0;

	return csg;
}


void free_CSG(CSG *csg)
{
	int i = 0;
	//int t = 0;
	int j = 0;
	int nout = 0;
	CSG_Node *actual;
	CSG_Node *base = csg->node;

	for (i = 0; i < csg->last_state; i++) {
	  actual = base + i;
	  nout = OUT(actual);

		if (0 && nout==1) {
		  if (((CSG_Node *)actual->arc <= base) 
			  ||((CSG_Node *)actual->arc >= (base + csg->last_state))) {
		    if(actual != NULL)
				free(actual);
				//printf("%s \n","uo");
		  }
	  }	
	  else { 
		  for(j=0;j<nout;j++) {
				CSG_Node *act = actual->arc[j];
				if (0 && ((act <= base) || (act >= (base + csg->last_state)))) {
				  if(act != NULL)
					free(act);
				}
				
			}
		    if(nout>2)
		      {
			if (csg->node[i].arc != NULL)
				free(csg->node[i].arc);
		      }
		}
	}
        if (csg->node != NULL)
   	    free((void *)csg->node);	
        if (csg != NULL)
    	    free((void *)csg);
}


int _conv(char c) {
	switch (c) {
	case 'a': case 'A': return 0;
	case 'c': case 'C': return 1;
	case 'g': case 'G': return 2;
	case 't': case 'T': return 3;
	case 'n': case 'N': return 4;
	case '$': return 5;
	default : printf("error\n"); exit(-1);
	}
}

char conv_(int c) {
	switch (c) {
	case 0: return 'a';
	case 1: return 'c';
	case 2: return 'g';
	case 3: return 't';
	case 4: return 'n';
	case 5: return '$';
	default : return (char)c;
	}
}

int prefix_comu(CSG *csg, CSG_Node *state, const char *pendent)
{
	const char *seq = csg->seq;
	int n     = csg->size;

	const char *s = LABEL(state);
	int l   = (IS_LEAF(state) ? (n - (s - seq)) : LEN_L(state));
	int p   = 0;

	while (*s == *pendent && p < l) {
		s++;
		pendent++;
		p++;
	}
	return p;
}

int class_link(CSG *csg, CSG_Node **state, CSG_Node **parent, int *lprefix)
{
	/* go to position with the class_link*/
	/* return true if the last position is the end of a node */
	const char *pendent = LABEL(*state);
	
	*state = CLINK(*state);
	while  (*lprefix > LEN_L(*state) &&  !IS_LEAF(*state)) {
		*lprefix -= LEN_L(*state);
		pendent  += LEN_L(*state);
		*parent   = *state;
		*state    = NEXT(csg, *state, *pendent);
	}
	return *lprefix == LEN_L(*state);
}

CSG_Node *insert_leaf(CSG *csg, const char *pendent, int split_lprefix)
{
	
	CSG_Node *newState  = INITIAL + ++csg->last_state;
	
	CLINK(newState)  = NINDEF;
 	LABEL(newState)  = pendent;
	LEN_L(newState)  = split_lprefix;
	return newState;
}

CSG_Node *to_clone(CSG *csg, CSG_Node *state, CSG_Node *nextState)
/* state: actual state; nextState: state to clone */
{
	
	CSG_Node *clon         = INITIAL + ++csg->last_state;
	CSG_Node *candidat     = state;
	CSG_Node *parent_state = state;
	int lprefix      = LEN_L(state);
	const char *prefix     = LABEL(nextState);
	int trobat       = TRUE;
	uint lg_clon     = (IS_INITIAL(state)? 0 :  LG(state) + LEN_L(state));  

	if (!IS_LEAF(nextState)) {
		INC_NNODES;
		CLONE_EDGES(nextState, clon);
	} else {
		INC_NLEAFS;
	}
	LABEL(clon)       = LABEL(nextState);
	LEN_L(clon)       = LEN_L(nextState);
	CLINK(clon)       = CLINK(nextState);
	CLINK(nextState)  = clon;
	do {
		if (NEXT(csg, candidat, *prefix) == nextState) {
			SWAP_EDGE(csg, candidat, nextState, clon);
		} else {
		        trobat = FALSE;
		}
		state = candidat;
		class_link(csg,  &candidat, &parent_state, &lprefix);
	} while (trobat && candidat != state);
	LG(clon) = lg_clon; 
	
	return clon;
}


CSG_Node *divide(CSG *csg, CSG_Node *state, int lprefix, CSG_Node *target)
{
	CSG_Node *waiting   = target;
	int waiting_lprefix = LEN_L(state) - lprefix;
	int is_leaf         = IS_LEAF(state);
	//int j;

	/* delete */
	LEN_L(state) = lprefix;
	if (is_leaf) {
		/* the leaf is converted to a node and is created a new leaf */
		DEC_NLEAFS;
		INC_NNODES;
		if (target == NINDEF || LABEL(target)!=(LABEL(state) + lprefix))
		{
			INC_NLEAFS;
			waiting = insert_leaf(csg, LABEL(state) + lprefix, LINDEF);
		}
	} else {
		INC_NNODES;
 		if (target == NINDEF) {
 			waiting = insert_leaf(csg, LABEL(state) + lprefix, waiting_lprefix);
		}
		LG(waiting) = (IS_INITIAL(state) ? 0 
		: LG(state) + LEN_L(state));
		waiting->arc = state->arc;
		OUT(waiting) = OUT(state);
		OUT(state)   = 0;
	}
	INS_EDGE(csg, state, waiting);
	return waiting;
}

CSG_Node *compres_leaf(CSG *csg, CSG_Node *source, CSG_Node *parent_state
, CSG_Node *state, CSG_Node *waiting, int lprefix)
{
	//const char *prefix = LABEL(state);

	LEN_L(source) = LEN_L(source) - lprefix;
	/* update actual state*/
	LEN_L(state) = lprefix;
	/* clear old edge */
	/* insert edge betwen the compress node and split node */
	SWAP_EDGE(csg, source, waiting, state);
	LG(state) = (IS_INITIAL(source)? 0 :  LG(source) + LEN_L(source));  
	/* insert waiting edge */
	DEC_NLEAFS;
	INC_NNODES;
	INS_EDGE(csg, state, waiting);
	return state;
}


CSG_Node *arrange_link(CSG *csg, CSG_Node *state, CSG_Node *waiting, int go_sl)
{
	int lprefix      = LEN_L(state);
	int is_end       = TRUE;
	CSG_Node *source       = state;
	CSG_Node *parent_state = state;
	const char *wprefix    = LABEL(waiting);
	CSG_Node *next;
	CSG_Node *split;
	CSG_Node *node_act = state;
	
	if (go_sl && waiting != NINDEF)
		is_end  = class_link(csg,  &state, &parent_state, &lprefix);
	while (!is_end) {
		if (IS_LEAF(state)) {
			state = compres_leaf(csg,  source, parent_state, state
			, waiting, lprefix);
			source = state;
			node_act = state;
  		} else {
			split = divide(csg,  state, lprefix, NINDEF);
   			CLINK(waiting) = split;
   			waiting = split;
			if (state == node_act)
				node_act = waiting;
		}
		is_end  = class_link(csg,  &state, &parent_state, &lprefix);
	}
	/* search for another leaf in the same t-path */
	next = NEXT(csg, state, *wprefix);
	while (CLINK(waiting) == NINDEF) {
		if (state != source && IS_LEAF(next)) {
			if (!GET_SOLID(csg, state, next)) {
				SWAP_EDGE(csg, state, next, waiting);
				is_end  = class_link(csg,  &state
    				, &parent_state, &lprefix);
				next = NEXT(csg, state, *wprefix);
			} else {
				SWAP_EDGE(csg, source, waiting, next);
				waiting = next;
				csg->last_state--;
				DEC_NLEAFS;
			}
		} else  {
			CLINK(waiting) = (next == NINDEF) || (next == waiting) 
			? INITIAL : next;
		}
	}
	return node_act;
}

CSG_Node *add_suffixes(CSG *csg, CSG_Node *state, CSG_Node *waiting
, const char **pendent)
{
	/* return a actual state (the last state after do the operations */
	CSG_Node *parent   = state;
	CSG_Node *source   = state;
	int lprefix        = LEN_L(state);
	const char *wprefix      = LABEL(waiting);
	int is_end;
	CSG_Node *next;
	CSG_Node *split;

	is_end  = class_link(csg,  &state, &parent, &lprefix) ;
	next    = NEXT(csg, state, *wprefix);
	
	while (!(is_end && (next != NINDEF || HAS_FINAL(state)))) {	
		if (!is_end) {
			split = divide(csg,  state, lprefix, NINDEF);
			arrange_link(csg,  state, split, FALSE);
		}
		INS_EDGE(csg, state, waiting);
		source = state;
		is_end  = class_link(csg,  &state, &parent, &lprefix) ;
		next    = NEXT(csg, state, *wprefix);
	}
	if (source == INITIAL && **pendent != LAST_SYMBOL) {
		CLINK(waiting) = INITIAL;
		(*pendent)++;
	} else if (**pendent != LAST_SYMBOL) {
 		if (GET_SOLID(csg, state, next))
 			state = next;
		else 
			state = to_clone(csg,  state, next);
		CLINK(waiting) = state;
	}
	return state;
}

				
CSG_Node *divide_and_compress(CSG *csg, CSG_Node *state, CSG_Node *cloned
, int lprefix)
{ 
	CSG_Node *waiting = NINDEF;
	//int j;

	if (!IS_LEAF(cloned)) {
 		waiting = divide(csg,  cloned, lprefix, NINDEF);
 		LEN_L(state) = lprefix;
		CLEAR_EDGES(state);
 		INS_EDGE(csg, state, waiting);
	} else {
		waiting = divide(csg,  state, lprefix, NINDEF);
	}
	return waiting;
}

void build_CSG(CSG *csg, const char *s, long int size, int debug)
{
	CSG_Node *state = INITIAL;
	int lprefix     = 1;
	const char *pendent   = s;
	CSG_Node *nextState;
	int offset;
	int cloned = FALSE;
	CSG_Node *node_cloned = NINDEF;

	int go;
	CSG_Node *waiting;

	while (*pendent != LAST_SYMBOL) {
		lprefix = state == INITIAL ? 1 
		: prefix_comu(csg,  state, pendent);
		offset  = state == INITIAL ? 0 : lprefix;
		if (lprefix == LEN_L(state)) {
			go        = *(pendent + offset);
			nextState = NEXT(csg, state, go);
			pendent  += offset;
			cloned = FALSE;
			if (nextState == NINDEF) {
				/* INSERT SUFIXES */
				INC_NLEAFS;
				waiting = insert_leaf(csg, pendent, LINDEF);
				INS_EDGE(csg, state, waiting);
				nextState =  add_suffixes(csg,  state, waiting, &pendent);
			} else if (!GET_SOLID(csg, state, nextState) 
			|| cloned) {
				node_cloned = nextState;
				nextState   = to_clone(csg,  state, nextState);
				cloned      = TRUE;
			}
		} else {
			if (!cloned)
				waiting = divide(csg,  state, lprefix, NINDEF);
			else {
				waiting = divide_and_compress(csg,  state
    				, node_cloned, lprefix);
        			cloned = FALSE;	
     			}
			state = arrange_link(csg,  state, waiting, TRUE);
			offset  = state == INITIAL ? 0 : lprefix;
			pendent += offset;
			INC_NLEAFS;
			waiting = insert_leaf(csg, pendent, LINDEF);
   			INS_EDGE(csg, state, waiting);
			nextState = add_suffixes(csg,  state, waiting, &pendent);
		}
		state = nextState;

	}
	// eliminate the last leaf "$"
	DEC_NLEAFS;

}


void find_leaves(CSG *csg)
{
	const char *seq  = csg->seq;
	ulong size = csg->size;
    const char *pendent = seq;
    CSG_Node *state = INITIAL;
    int  nleaf = -1;
    CSG_Node *last_leaf = NINDEF;
    CSG_Node *next;
	int shucount = 0;
	int smallest = 0;

    while (*pendent != LAST_SYMBOL) 
	{
		if (IS_LEAF(state)) 
		{
			if (++nleaf) 
				SET_NEXT_LEAF(last_leaf, state);
			last_leaf = state;
			
			
			if (smallest == 0 || state->lg < smallest )
			{
				shucount = 0;
				smallest = state->lg;
				shucount+=1;
			}
			else if( smallest == state->lg)
			{
				shucount+=1;
			}
		
			state = CLINK(state);
			if (state == INITIAL) 
				pendent++;
			
		}
		else if(1)
		{
			if (smallest == 0 || state->lg < smallest )
			{
				shucount = 0;
				smallest = state->lg;
				shucount+=1;
			}
			else if( smallest == state->lg)
			{
				shucount+=1;
			}

		}
		if (state != INITIAL) 
		{
			pendent += LEN_L(state);
		}
		
		next = NEXT(csg, state, *pendent);
		state = next;
    }

	state = (CSG_Node *)malloc(sizeof(CSG_Node));
	state->prefix  = seq + size;
	state->lprefix = -1;
	state->s        = INITIAL;
	
	SET_NEXT_LEAF(last_leaf, state);
	SET_NEXT_LEAF(state, INITIAL);
			
}

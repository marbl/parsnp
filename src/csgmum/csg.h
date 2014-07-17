/////////////////////////////////////////
// csg.c
// compressed suffix graph structure
/////////////////////////////////////////

#ifndef _CSG_H_
#define _CSG_H_

#include "types.h"

/* constants */
#define NINDEF NULL 
#define LINDEF -1
#define MAXALPHA 5 
#define INITIAL csg->node
#define FINAL 2
#define FALSE 0
#define TRUE 1
#define LAST_SYMBOL 5

#define MAX(x,y) ((x) > (y) ? (x) : (y)) 

#define LG(state)                       ((state)->lg)
#define OUT(state)                      (state)->nout
#define HAS_FINAL(state)                ((state)->is_final)
#define SET_FINAL(state)                ((state)->is_final = 1)
#define SET_NOFINAL(state)              ((state)->is_final = 0)

#define LABEL(state)			((state)->prefix)
#define LEN_L(state)			((state)->lprefix)
#define CLINK(source)			(source)->s
#define IS_LEAF(state)			(LEN_L(state) == LINDEF)
#define NSOLID(source)			(source)->in_solid
#define INC_NNODES			csg->num_nodes++
#define INC_NLEAFS			csg->num_leafs++
#define DEC_NLEAFS			csg->num_leafs--

#define IS_INITIAL(state)		(state == INITIAL)
#define NEXT_LEAF(source)     (IS_LEAF(source)?(CSG_Node *)source->arc : NINDEF)
#define SET_NEXT_LEAF(source, next) (source->arc = (CSG_Node **)next)

typedef struct CSG_Node CSG_Node;
struct CSG_Node {
	CSG_Node **arc;		/* Transition Vector */
	CSG_Node *s;			/* suffix link */
	const char *prefix;
	int   lprefix;
	uint  nout : 3;
	uint is_final: 1;
	uint lg;
};

typedef struct struct_csg CSG;
struct struct_csg {
	CSG_Node *node;
	int last_state;
	ulong ini_node;
	int num_nodes;
	int num_leafs;
	const char *seq;
	ulong size;
};

CSG_Node *NEXT(CSG *csg, CSG_Node *source, int target) ;
int GET_SOLID(CSG *csg, CSG_Node *source, CSG_Node *node_t);

int _conv(char c);
char conv_(int c);

CSG *new_CSG(CSG * csg_clone,ulong node_ini, const char *seq, ulong size, int setmem);
void free_CSG(CSG *);
void build_CSG(CSG *csg, const char *s, long int size, int out);
int prefix_comu(CSG *csg, CSG_Node *state, const char *pendent);
void find_leaves(CSG *csg);

#endif


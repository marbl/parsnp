//#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

//#include <sys/types.h>
#include <fcntl.h>

#if defined(_WIN32)
#include <io.h>
#include <share.h>
#  pragma warning(disable : 4996) 
#endif
#include <sys/stat.h>
#include "fastaread.h"

#define READBUFF(c) \
	do {							\
		if (!reads) {					\
			reads = read(fpin, buf, sizeof(buf));	\
			bufp = buf;				\
		}						\
		if (--reads < 0)				\
			*c = 0;					\
		else					      	\
	        	*c = (unsigned char) *bufp++;           \
} while (0)

char Sf(char c) {
	switch (c) {
	case 'a': case 'A': return 0; break;
	case 'c': case 'C': return 1; break;
	case 'g': case 'G': return 2; break;
	case 'n': case 'N': return 4; break;
	case 't': case 'T': case 'u': case 'U': return 3; break;
	default: return LAST_SYMBOL;
	}
}

char _Sf(char c) {
	switch (c) {
	case 0: return 'a'; break;
	case 1: return 'c'; break;
	case 2: return 'g'; break;
	case 3: return 't'; break;
	case 4: return 'n'; break;
	default: return '$';
	}
}

char *to_string(char *str, int len)
{
	char *result = (char *)calloc(sizeof(char), len + 1);
	int i = 0;
	
	for (i = 0; i < len; i++)
		result[i] = _Sf(str[i]);
	return result;
}

char *to_code(char *str, int len)
{
	char *result = (char *)calloc(sizeof(char), len + 1);
	int i = 0;
	
	for (i = 0; i < len; i++)
		result[i] = Sf(str[i]);
	return result;
}

/*
/////////////////////////////////////////////////////////////////
// FILTER(C):
// 	We read the input file in BUFSIZ blocks, and when we have
//      read the block then load another from the disk.
//	The BUFSIZ is maximal optimum.
//      This function process the actual char and filter any non bp
//      char.
//	return only a know char or 5 if is the end.
//////////////////////////////////////////////////////////////////
*/
void
filter(char *c, int ini, int fpin)
{
	static char buf[BUFSIZ];
	static char *bufp = buf;
	static int reads = 0;
	char stop = 0;

	if (ini)
		reads = 0;

	while (!stop) {
		READBUFF(c);
		if (!*c) {
			stop = 1;
			break;
		}
		switch (*c) {
			case 'a': case 'A':
		        case 'c': case 'C':
			case 'g': case 'G':
			case 'n': case 'N':
			case 't': case 'T':
			case 'u': case 'U':
				stop = 1;
				break;
			case '\n': case '\t': case ' ':
				break;
			case '>':
				do {
					READBUFF(c);
				} while (*c && (*c != '\n'));
				break;
			default:
				break;
		}
	}
	if (stop && *c) {
		if (*c == 'A') *c = 'a';
		else if (*c == 'C') *c = 'c';
		else if (*c == 'G') *c = 'g';
		else if (*c == 'T') *c = 't';
		else if (*c == 'N') *c = 'n';
	}	
}


/* llegim la sequencia */
unsigned char *get_seq(char *name,long int *size) {
	long i_b = 0;
	char c;
	unsigned char *seq;
	struct stat atbuf;
	int fpin;

	stat(name, &atbuf);
	//*size = atbuf.st_size - 1;
	*size = atbuf.st_size;

	/* inicialitza taula char2int */
	if ((fpin = open(name, O_RDONLY)) == -1) {
		printf("I can't open the input file\n");
		exit(-1);
	}

	if (!(seq = (unsigned char *)malloc(*size + 2))) {
		printf("size NOK \n");
		exit(0);
	}

	filter(&c, 1, fpin);
	while (c) {
		  seq[i_b++] = Sf(c);
		  filter(&c, 0, fpin);
	}
	seq[i_b]   = Sf('$');
	seq[i_b+1] = '\0';
	close(fpin);
	*size = i_b;
	return seq;
}


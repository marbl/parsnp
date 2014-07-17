
#ifndef _FASTA_READ_H_
#define _FASTA_READ_H_

#define LAST_SYMBOL 5

unsigned char *get_seq(char *name,long int *size);
char *to_string(char *str, int len);
char *to_code(char *str, int len);

#endif

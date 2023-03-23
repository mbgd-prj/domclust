#ifndef _SSTACK_H
#define SSTACK_ERROR -999999555

typedef struct sStack {
	char *data;
	int size;
	int datasize;
	int currpos;
} sStack;

sStack *create_sStack(int stacksize, int datasize);
char *pop_sStack(sStack *st);
char *get_sStackIdx(sStack *st, int idx);

#define _SSTACK_H
#endif

#include <stdio.h>
#include <stdlib.h>
#include "sstack.h"

sStack *create_sStack(int stacksize, int datasize)
{
	sStack *stack;
	if ((stack = (sStack*) malloc(sizeof(sStack))) == NULL) {
		fprintf(stderr, "Can't allocate memory: stack\n");
		exit(1);
	}
	if ((stack->data = (char*) malloc(datasize * stacksize)) == NULL) {
		fprintf(stderr, "Can't allocate memory: stack.data\n");
		exit(1);
	}
	stack->currpos = 0;
	stack->size = stacksize;
	stack->datasize = datasize;
	return(stack);
}
free_sStack(sStack *stack)
{
	free(stack->data);
	free(stack);
}
reset_sStack(sStack *stack)
{
	stack->currpos = 0;
}
push_sStack(sStack *stack, char *data)
{
	if (stack->currpos >= stack->size) {
		fprintf(stderr, "stack overflows\n");
		return(-1);
	}
	memcpy(&(stack->data[stack->currpos * stack->datasize]), data, stack->datasize);
	stack->currpos++;
	return (0);
}
char *pop_sStack(sStack *stack)
{
	if (stack->currpos <= 0) {
		fprintf(stderr, "stack underflows\n");
		return(NULL);
	}
	--stack->currpos;
	return(&(stack->data[stack->currpos * stack->datasize]));
}
char *get_sStackIdx(sStack *stack, int idx)
{
	if (idx < -stack->currpos  || idx >= stack->currpos) {
/**
		fprintf(stderr, "stack: acceess point out of range: %d,%d\n",
				idx, stack->currpos);
**/
		return(NULL);
	}
	if (idx >= 0) {
		return(&(stack->data[idx * stack->datasize]));
	} else {
		return(&(stack->data[(stack->currpos + idx) * stack->datasize]));
	}
}

#ifdef STACK_DEBUGMAIN
main(int argc, char **argv)
{
	sStack *st;
	int i, *p;
	st = create_sStack(200, sizeof(int));
	i = 10; push_sStack(st, (char*) &i);
	i = 50; push_sStack(st, (char*) &i);
	i = 0; push_sStack(st, (char*) &i);
	i = 20; push_sStack(st, (char*) &i);
	i = 5; push_sStack(st, (char*) &i);
	p = (int*) pop_sStack(st); printf("%d\n", *p);
	p = (int*) pop_sStack(st); printf("%d\n", *p);
	p = (int*) get_sStackIdx(st, -1); printf("get>-1>%d\n", *p);
	p = (int*) get_sStackIdx(st, 0); printf("get>1>%d\n", *p);
	p = (int*) pop_sStack(st); printf("%d\n", *p);
}
#endif

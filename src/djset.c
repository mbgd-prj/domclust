/*
 * djset.c: disjoint set data structure 
 * Copyright (c) 2000-2008, Ikuo Uchiyama
 * All rights reserved.
 */
#include <stdlib.h>
#include <stdio.h>
#include "djset.h"

DJset *create_DJset(int nel)
{
	DJset * djset;
	int i;
	if ( (djset = (DJset *) malloc(sizeof (DJset))) == NULL ) {
		fprintf(stderr, "Can't allocate memory (DJset)\n");
		exit(1);
	}
	if ( (djset->data = (DJsetData *) malloc(
			sizeof (DJsetData) * nel)) == NULL ) {
		fprintf(stderr, "Can't allocate memory (DJsetData)\n");
		exit(1);
	}
	for (i = 0; i < nel; i++) {
		djset->data[i].parent = i;
		djset->data[i].rank = 0;
	}
	djset->numelem = nel;
	return djset;
}
free_DJset(DJset *djset)
{
	free(djset->data);
	free(djset);
}

DJset_merge(DJset *djset, int i, int j)
{
	if (i >= djset->numelem || j >= djset->numelem) {
		fprintf(stderr, "DJset: Too large number: %d,%d,%d\n",
			i,j,djset->numelem);
		exit(1);
	}
	DJset_link(djset, DJset_find_repr(djset,i), DJset_find_repr(djset,j));
}
DJset_link(DJset *djset, int i, int j)
{
	if (i == j) {
		return(0);
	}
	if (djset->data[i].rank == djset->data[j].rank) {
		if (i < j) {
			djset->data[j].parent = i;
			djset->data[i].rank ++;
		} else {
			djset->data[i].parent = j;
			djset->data[j].rank ++;
		}
	} else if (djset->data[i].rank > djset->data[j].rank) {
		djset->data[j].parent = i;
	} else {
		djset->data[i].parent = j;
	}
}

DJset_find_repr(DJset *djset, int i)
{
	if (i >= djset->numelem) {
		fprintf(stderr, "DJset: Too large number: %d,%d\n",
			i, djset->numelem);
abort();
		exit(1);
	}
	if (djset->data[i].parent != i) {
		djset->data[i].parent = DJset_find_repr(
			djset, djset->data[i].parent);
	}
	return(djset->data[i].parent);
}


#ifdef  DJSET_MAIN
#include <string.h>
main(int argc, char ** argv)
{
	char filename[120];
	FILE *fp;
	DJset *djset;
	char buf[BUFSIZ];
	int i, j;
	int nel = 20;

	if (argc > 1 && strcmp(argv[1], "-")!=0) {
		strcpy(filename, argv[1]);
		if ((fp = fopen(filename, "r")) == NULL) {
			fprintf(stderr, "Can't open %s\n", filename);
			exit(1);
		}
	} else {
		fp = stdin;
	}
	if (argc > 2) {
		nel = atoi(argv[2]);
	}
	djset = create_DJset(nel);
	while (fgets(buf, BUFSIZ, fp) != NULL) {
		sscanf(buf, "%d%d", &i, &j);
		DJset_merge(djset, i,j);
	}
	for (i = 0; i < nel; i++) {
		printf("%d %d\n", i, DJset_find_repr(djset, i));
	}
}
#endif

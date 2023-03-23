/*
 * DomClust: Hierarchical Clustering for Orthologous Domain Classification
 * Copyright (c) 2000-2007, Ikuo Uchiyama
 * All rights reserved.
 */

#include "seqreg.h"
#include <stdio.h>
#include <stdlib.h>

cmpr_poslist(PosList_t *poslist1, PosList_t *poslist2)
{
	int tmp;
	if (tmp = (int) poslist1->pos - (int) poslist2->pos) {
		return tmp;
	} else {
		return strcmp(poslist1->type, poslist2->type);
	}
}
PosList *createPosList(int size)
{
	PosList *pol;
	if ((pol = (PosList *) malloc(sizeof(PosList))) == NULL) {
		fprintf(stderr, "Can't alloc memory\n");
		exit(1);
	}
	if ((pol->poslist = (PosList_t *) malloc(sizeof(PosList_t)*size))==NULL) {
		fprintf(stderr, "Can't alloc memory\n");
		exit(1);
	}
	pol->num = 0;
	pol->segnum = 0;
	pol->maxsize = size;
	return pol;
}

initPosList(PosList **poslist)
{
	if (! *poslist) {
		*poslist = createPosList(8);
	} else {
		initPosList0(*poslist);
	}
}
initPosList0(PosList *pl)
{
	pl->num = pl->segnum = 0;
}
addPosList(PosList *poslist, Region *reg, char *str)
{
	PosList_t *pl;
	if (poslist->num > poslist->maxsize) {
		fprintf(stderr, "poslist overflows\n");
		exit(1);
	}
	pl = &(poslist->poslist[poslist->num++]);
	strcpy(pl->type, str);
	strcat(pl->type, "f");
	pl->pos = reg->from;
	pl->seg = poslist->segnum;

	poslist->num++; pl++;
	strcpy(pl->type, str);
	strcat(pl->type, "t");
	pl->pos = reg->to;
	pl->seg = poslist->segnum++;
}
sortPosList(PosList *poslist)
{
	qsort(poslist->poslist, poslist->num, sizeof(PosList_t), cmpr_poslist);
}
printPosList(PosList *poslist)
{
	int i;
	for (i = 0; i < poslist->num; i++) {
		printf(">%d,%s,%d\n",
			poslist->poslist[i].seg,
			poslist->poslist[i].type,
			poslist->poslist[i].pos);
	}
}
freePosList(PosList *poslist)
{
	free(poslist);
}

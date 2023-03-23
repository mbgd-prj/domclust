/*
 * sgraph.c: simple graph data structure for single linkage clustering
 * Copyright (c) 2000-2007, Ikuo Uchiyama
 * All rights reserved.
 */

#include <stdio.h>
#include <stdlib.h>
#include "sgraph.h"
#include "plist.h"
#include "hash.h"

/***
typedef struct sGraphNode{
	int id;
	pList *nlist;
	int cluster;
} sGraphNode;

typedef struct sGraph{
	sGraphNode *nodes;
	Hash *nodehash;
	int nodenum;
	int clustnum;
} sGraph;
sGraph *create_sGraph(int nodenum, int edgenum);
***/

sGraph *create_sGraph(int nodenum, int edgenum)
{
	sGraph *sg;
	if ((sg = (sGraph *) malloc(sizeof(sGraph)))==NULL) {
		fprintf(stderr, "Can't alloc memory\n");
		exit(1);
	}
	if ((sg->nodes = (sGraphNode*) malloc(sizeof(sGraphNode)*nodenum))
			== NULL) {
		fprintf(stderr, "Can't allocate memory\n");
		exit(1);
	}
	sg->nodehash = Hcreate(nodenum * 120);
	sg->nodenum = 0;
	sg->clustnum = 0;
	return sg;
}
free_sGraph(sGraph *sg)
{
	Hdestroy(sg->nodehash);
	free(sg->nodes);
	free(sg);
}
sGraphNode *sgraph_hashNode(sGraph *sg, int nodeid)
{
	sGraphNode *node;
	HENTRY hent;

	node = &(sg->nodes[sg->nodenum]);
	hent.key= nodeid; hent.datum = (void*) node;
	if (HIsearch(sg->nodehash, &hent, ENTER) == 1) {
		/* already entered */
		node = (sGraphNode *) hent.datum;
	} else {
		node->id = nodeid; node->nlist = create_pList();
		node->cluster = 0;
		sg->nodenum++;
	}
	return node;
}
sgraph_addEdge(sGraph *sg, int a, int b)
{
	sGraphNode *n1,*n2;
	HENTRY hent;

	n1 = sgraph_hashNode(sg, a);
	n2 = sgraph_hashNode(sg, b);

#ifdef AHO
	n1 = &(sg->nodes[sg->nodenum]);
	hent.key= a; hent.datum = (void*) n1;
	if (HIsearch(sg->nodehash, &hent, ENTER) == 1) {
		/* already entered */
		n1 = (sGraphNode *) hent.datum;
	} else {
		n1->id = a; n1->nlist = create_pList();
		n1->cluster = 0;
		sg->nodenum++;
	}

	n2 = &(sg->nodes[sg->nodenum]);
	hent.key= b; hent.datum=(void*) n2;
	if (HIsearch(sg->nodehash, &hent, ENTER) == 1) {
		/* already entered */
		n2 = (sGraphNode *) hent.datum;
	} else {
		n2->id = b; n2->nlist = create_pList();
		n2->cluster = 0;
		sg->nodenum++;
	}
#endif
	pushList(n1->nlist, n2);
	pushList(n2->nlist, n1);
}
sgraph_slink(sGraph *sg)
{
	sGraphNode *n1;
	int i;
	for (i = 0; i < sg->nodenum; i++) {
		n1 = &(sg->nodes[i]);
		if (n1->cluster) {
			continue;
		}
		sg->clustnum++;
		sgraph_slink0(sg, n1);
	}
}
sgraph_slink0(sGraph *sg, sGraphNode *n1)
{
	sGraphNode *n2;
	listIter iter;
	if (n1->cluster) {
		return;
	}
	n1->cluster = sg->clustnum;
	setListIter(&iter, n1->nlist, 1);
	while (n2 = (sGraphNode*) getListIter(&iter)) {
		sgraph_slink0(sg, n2);
	}
}
sgraph_print_cluster(sGraph *sg)
{
	int i;
	for (i = 0; i < sg->nodenum; i++) {
		printf("%d, %d\n", sg->nodes[i].id, sg->nodes[i].cluster);
	}
}

#ifdef SLINK_MAIN
#include "namehash.h"
#define MAXNODE 50000
#define MAXEDGE 800000
main(int argc, char **argv)
{
	sGraph *sg;
	int i;
	FILE *fp;
	int id1, id2, newid, nid;
	NameHash *nhash;
	char filename[120];
	char buf[BUFSIZ];
	char name1[20], name2[20], *namep;
	char *names[MAXNODE];
	double dist, cutoff = 0.0;

	nhash = initNames(MAXNODE);
	sg = create_sGraph(MAXNODE, MAXEDGE);

	if (argc > 1) {
		strcpy(filename, argv[1]);
		if (strcmp(filename, "-")==0) {
			fp = stdin;
		} else if ((fp = fopen(filename, "r")) == NULL) {
			fprintf(stderr, "Can't open %s\n", filename);
			exit(1);
		}
		if (argc > 2) {
			cutoff = atof(argv[2]);
		}
	} else {
		fp = stdin;
	}
	newid = 1;
	while (fgets(buf, BUFSIZ, fp) != NULL) {
		sscanf(buf, "%s%s%lf",name1,name2,&dist);
		if ( (id1 = getNameID(nhash, name1)) < 0 ) {
			id1 = newid++;
			namep = addName(nhash, name1, id1);
			names[id1] = namep;
		}
		if ( (id2 = getNameID(nhash, name2)) < 0 ) {
			id2 = newid++;
			namep = addName(nhash, name2, id2);
			names[id2] = namep;
		}
		if (! cutoff || dist < cutoff) {
			sgraph_addEdge(sg, id1, id2);
		} else {
			sgraph_hashNode(sg, id1);
			sgraph_hashNode(sg, id2);
		}
	}
	sgraph_slink(sg);

/*
	sgraph_print_cluster(sg);
*/
	for (i = 0; i < sg->nodenum; i++) {
		nid = sg->nodes[i].id;
		printf("%s %d %d\n", names[nid],
			nid, sg->nodes[i].cluster);
	}

	free_sGraph(sg);
}
#endif

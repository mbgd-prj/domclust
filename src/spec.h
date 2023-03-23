/*
 * DomClust: Hierarchical Clustering for Orthologous Domain Classification
 * Copyright (c) 2000-2007, Ikuo Uchiyama
 * All rights reserved.
 */

#ifndef _SPEC_H_
#define _SPEC_H_

#define SPFLAGSIZ 800
#define MAXSP (SPFLAGSIZ * sizeof(char) * 8)
#define MAXSPBUF MAXSP*4
#define MAXTAXNUM 20


int SPnum;
/*
int bitcnt[256] = {
	0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,
	4,5,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,
	4,5,5,6,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,
	4,5,4,5,5,6,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,
	4,5,5,6,5,6,6,7,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,
	4,5,3,4,4,5,4,5,5,6,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,
	4,5,5,6,4,5,5,6,5,6,6,7,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,
	4,5,4,5,5,6,4,5,5,6,5,6,6,7,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
	4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8
};
double wbitcnt[SPFLAGSIZ][256];
*/
NameHash *SpHash;
char *SPnames[MAXSP];
double SPweights[MAXSP];
typedef unsigned char specFlag[SPFLAGSIZ];
typedef unsigned char *specFlagP;

#define FLAG_TAXOR 1
#define FLAG_TAXAND 2

typedef struct SPTreeNode {
	int spid;
	int parent;
	int child, sibling;
	double weight;
	specFlag spflag;
	char flag;
	char *name;
} SPTreeNode;

typedef struct SPTree {
	SPTreeNode node[MAXSPBUF];
	int nodenum;
} SPTree;
SPTree spTree;

struct {
	specFlag inGroup, outGroup;
	specFlag spMask;
	specFlag meta;
	specFlag taxQuery;
	specFlag ignore;
	specFlag unknown;
	specFlag partial;
} SPflags;

char *getSPname(), *getTaxName();
double spFlagCntW(), spFlagCntW_All(), spFlagANDcntW(), getSPweight(), sptree_MatchFlagsCntW(),
	sptree_spFlagCountTaxOrW();
int parse_spinfo(char *);
int readSPfile(char *);
#endif

/*
 * DomClust: Hierarchical Clustering for Orthologous Domain Classification
 * Copyright (c) 2000-2007, Ikuo Uchiyama
 * All rights reserved.
 */
#include "domclust.h"
#include "vararray.h"
#include "util.h"
#include "seqreg.h"
#include "hitlist.h"
#include <assert.h>

#define DUMMY 99999999
/*
#define MAXNBR 15000
*/

#define MAXNEWDOM 5 /* 1 central, 2 left and 2 right nodes */

static int DEBUGFLAG = 0;
static int DEBUGSTAT = 0;

#define IsMatchedSeg(reg) ((reg)==D_MATCH1||(reg)==D_MATCH2)
static DomIdx Idx1[3] = {D_LEFT1, D_MATCH1, D_RIGHT1};
static DomIdx Idx2[3] = {D_LEFT2, D_MATCH2, D_RIGHT2};
static DomIdx Idx3[3] = {D_THIRD, D_THIRD, D_THIRD};
/***
static int NewNodeIdx[7] = {0, 1, 2, 3, 4, 0, 0};
***/
static int segLen[8];

typedef enum {LEFT, RIGHT, UP, DOWN} Direction;


int call_domCluster_core(Edge *bestedge, char **args);
int createNewSideNode(NodeSet *nodes, Edge *best, Node *n, SeqPos brkpnt, int newlen, int nodeflag, Node **newnodes);
calNewAli_Normal( Region *ali12, Region *ali21, Region *ali13, Region *ali31,
	Region *ali23, Region *ali32, Region *brk1, Region *brk2,
	NewSeqInfo *newseq, Region *newreg, Region *newconsreg,
	Region *newali, Region *newali3,
	char *hitflag, HitList *hitlist,
	Count cnt1, Count cnt2, Count cnt3, Node *n1, Node *n2, Node *n3,
	Edge *e1, Edge *e2);
calNewAliM( Region *ali12, Region *ali21,
	Region *ali12_2, Region *ali21_2,
	Region *brk1, Region *brk2,
	NewSeqInfo *newseq, Region *newali, char *hitflag,
	Count cnt1, Count cnt2, Node *n1, Node *n2, Edge *e1, HitList *hitlist);
calNewAliSelf3( Region *ali12, Region *ali21,
	Region *ali11_1, Region *ali11_2,
	Region *brk1, Region *brk2,
	NewSeqInfo *newseq, Region *newali, char *hitflag,
	Count cnt1, Count cnt2, Node *n1, Node *n2, Edge *e1,
	HitList *hitlist, int swapflag);
calNewAliSelf( Region *ali12, Region *ali21, Region *ali13, Region *ali31,
	Region *brk1, Region *brk2, NewSeqInfo *newseq,
	Region *newali, Region *newali3, char *hitflag,
	Count cnt1, Count cnt2, Node *n1, Node *n2, Edge *e1,
	HitList *hitlist);

getHitSeg(int seg){
	if (seg==D_MATCH1||seg==D_MATCH2) {
		return 0;
	} else {
		return seg;
	}
}

int calNewDist_NEW(int nn, Count cnt1, Count cnt2, Dist dist1, Dist dist2,
		Dist score1, Dist score2, ConnCount conn1, ConnCount conn2,
		Dist *newdist, Dist *newscore, ConnCount *newconn);

calNewAli_NEW(
    /* aliXY: aligment region on sequence X against sequence Y */
        Region *ali12,          /* region on seq1 aligned against seq2 */
        Region *ali21,          /* region on seq2 aligned against seq1 */
        Region *ali13,          /* region on seq1 aligned against seq3 */
        Region *ali31,          /* region on seq3 aligned against seq1 */
        Region *ali23,          /* region on seq2 aligned against seq3 */
        Region *ali32,          /* region on seq3 aligned against seq2 */
        Region *brk1, Region *brk2,     /* break points on seq1 and seq2 */
        NewSeqInfo *newseq,
        Region *newreg,         /* new sequence by merging seq1 and seq2 */
        Region *newconsreg,     /* aligned reg. (seq1 vs seq2) on newreg */
        Region *newali,         /* output: aligned reg on newreg against seq3 */
        Region *newali3,        /* output: aligned reg on seq3 against newreg */
        char *hitflag,
        HitList *hitlist,
        Count cnt1, Count cnt2, Count cnt3, Node *n1, Node *n2, Node *n3,
        Edge *e1, Edge *e2, EdgeType edgeType)
{
	if (edgeType != MULTI_EDGE && (n1 == n3 || n2 == n3)) {
		edgeType = SELF_EDGE;
	}
/*
printf("ET=%d\n",edgeType);
*/

	clear_segLen(segLen);
if (ali12->to > n1->len){
printf("ali12>>>%d,%d\n",ali12->to,n1->len);
}

	if (edgeType == SELF_EDGE) {
		Region *ali11_1, *ali11_2, *tmp_brk1, *tmp_brk2, *tmp;
		Node *nn1, *nn2;
		Edge *ee1;
		int swap;
if (Opt.DEBUG) {
printf("SELF_3_EDGE\n");
}
		if (n1 == n3) {
			ali11_1 = ali13; ali11_2 = ali31;
			nn1 = n1; nn2 = n2;
			tmp_brk1 = brk1; tmp_brk2 = brk2;
			ee1= e1;
			swap = 0;
		} else {
			tmp = ali12; ali12 = ali21; ali21 = tmp;
			ali11_1 = ali23; ali11_2 = ali32;
			nn1 = n2; nn2 = n1;
			tmp_brk1 = brk2; tmp_brk2 = brk1;
			ee1= e2;
			swap = 1;
		}
		calNewAliSelf3(ali12, ali21, ali11_1, ali11_2,
			tmp_brk1, tmp_brk2,
			newseq, newali, hitflag,
			nn2->cnt, nn2->cnt, nn1,nn2,ee1, hitlist, swap);
	} else if (edgeType == MULTI_EDGE) {
	/** calculate new alignment regions (newali and newali3) **/
if (Opt.DEBUG) {
printf("MULTI_EDGE\n");
}
		assert (e1 == e2);
		calNewAliM(ali12, ali21, ali13, ali31, brk1, brk2,
			newseq, newali, hitflag,
			n1->cnt, n2->cnt, n1,n2,e1, hitlist);
	} else if (n1 == n2) {
if (Opt.DEBUG) {
printf("SELF_EDGE\n");
}
		if (!ali13){
			ali13 = ali23;
			ali31 = ali32;
		}
		calNewAliSelf(ali12, ali21, ali13, ali31,
			brk1, brk2,
			newseq, newali,
			newali3, hitflag,
			n1->cnt, n2->cnt, n1,n2,e1, hitlist);
	} else {
	        calNewAli_Normal(ali12, ali21, ali13, ali31, ali23, ali32,
			brk1, brk2,
			newseq, &(newseq->seq), &(newseq->aliM),
			newali, newali3, hitflag, hitlist,
			n1->cnt, n2->cnt, n3->cnt, n1,n2,n3,e1,e2);
	}
}

calNewAli_Normal(
    /* aliXY: aligment region on sequence X against sequence Y */
	Region *ali12,		/* region on seq1 aligned against seq2 */
	Region *ali21,		/* region on seq2 aligned against seq1 */
	Region *ali13,		/* region on seq1 aligned against seq3 */
	Region *ali31,		/* region on seq3 aligned against seq1 */
	Region *ali23,		/* region on seq2 aligned against seq3 */
	Region *ali32,		/* region on seq3 aligned against seq2 */
	Region *brk1, Region *brk2,	/* break points on seq1 and seq2 */
	NewSeqInfo *newseq,
	Region *newreg,		/* new sequence by merging seq1 and seq2 */
	Region *newconsreg,	/* aligned reg. (seq1 vs seq2) on newreg */
	Region *newali,		/* output: aligned reg on newreg against seq3 */
	Region *newali3,	/* output: aligned reg on seq3 against newreg */
	char *hitflag,
	HitList *hitlist,
	Count cnt1, Count cnt2, Count cnt3, Node *n1, Node *n2, Node *n3,
	Edge *e1, Edge *e2)
{
	Region brk3, tmp1, tmp2, tmp3, tmp4;
	int i, j;
	int numh1, numh2;
	static PosList *poslist1, *poslist2, *poslist31, *poslist32;
	Region brk1v, brk2v;
	HitList_t orighit;

	initPosList(&poslist1);
	initPosList(&poslist2);
	initPosList(&poslist31);
	initPosList(&poslist32);

	if (ali13 && ali31) {
		addPosList(poslist1, brk1, "b1");
		addPosList(poslist1, ali13, "a1");
		sortPosList(poslist1);

		addPosList(poslist31, ali31, "a1");
		copyReg(&brk1v, brk1); transformReg(&brk1v, ali13, ali31);
		addPosList(poslist31, &brk1v, "b1");
		sortPosList(poslist31);
	}
	if (ali23 && ali32) {
		addPosList(poslist2, brk2, "b1");
		addPosList(poslist2, ali23, "a1");
		sortPosList(poslist2);
		addPosList(poslist32, ali32, "a1");
		copyReg(&brk2v, brk2); transformReg(&brk2v, ali23, ali32);
		addPosList(poslist32, &brk2v, "b1");
		sortPosList(poslist32);
	}

if (Opt.DEBUG) {
	printf("---\n");
	printPosList(poslist1);
	printf("---\n");
	printPosList(poslist31);
	printf("---\n");
	printf("---\n");
	printPosList(poslist2);
	printf("---\n");
	printPosList(poslist32);
	printf("---\n");
}

	initHitList(hitlist);
	orighit.dist = orighit.score = orighit.connect = 0;
	if (ali13 && ali31) {
		createOrigHit(&orighit, ali31, e1, 1);
		calNewAli_common(poslist1, poslist31, Idx1, Idx3,
			n1, n3, newseq, 3, hitlist, &orighit, segLen,
			FlagHitList_NEW);
		postCheckHitList(hitlist,&orighit,segLen);
	}
	if (ali23 && ali32) {
		createOrigHit(&orighit, ali32, e2, 2);
		calNewAli_common(poslist2, poslist32, Idx2, Idx3,
			n2, n3, newseq, 3, hitlist, &orighit, segLen,
			FlagHitList_ADD);
		postCheckHitList(hitlist,&orighit,segLen);
	}
/*
	postCheckHitList(hitlist,NULL,NULL);
*/
}

/** for multi-edge **/
calNewAliM(
    /* aliXY: aligment region on sequence X against sequence Y */
	Region *ali12,		/* region on seq1 aligned against seq2 */
	Region *ali21,		/* region on seq2 aligned against seq1 */
	Region *ali12_2,	/* the other reg on seq1 aligned against seq2 */
	Region *ali21_2,	/* the other reg on seq2 aligned against seq1 */
	Region *brk1, Region *brk2,	/* break points on seq1 and seq2 */
	NewSeqInfo *newseq,
	Region *newali,		/* output: aligned reg on newreg against seq3 */
	char *hitflag,
	Count cnt1, Count cnt2, Node *n1, Node *n2, Edge *e1,
	HitList *hitlist)
{
	int i, j, k;
	int r1, r2;
	char flag;
	Region brk1v, brk2v, tmp1, tmp2, tmp3;
	SeqPos from, to, pos1,pos2,prevpos1,prevpos2;
	SeqPos bk1,bk2,prevbk1,prevbk2;
	static PosList *poslist1, *poslist2;
	int numh;
	HitList_t orighit;

	/* the break points in seq1,2 are mapped onto the other seq */
	copyReg(&brk1v, brk1); transformReg(&brk1v, ali12_2, ali21_2);
	copyReg(&brk2v, brk2); transformReg(&brk2v, ali21_2, ali12_2);

	initPosList(&poslist1);
	initPosList(&poslist2);

	addPosList(poslist1, brk1, "b1");
	addPosList(poslist1, &brk2v, "b2");
	addPosList(poslist1, ali12_2, "a1");
	sortPosList(poslist1);

	addPosList(poslist2, brk2, "b2");
	addPosList(poslist2, &brk1v, "b1");
	addPosList(poslist2, ali21_2, "a1");
	sortPosList(poslist2);

 
	initHitList(hitlist);
	createOrigHit(&orighit, ali12_2, e1, 1);
	calNewAli_common(poslist1, poslist2, Idx1, Idx2,
		n1, n2, newseq, 3, hitlist, &orighit, segLen,
		FlagHitList_NEW);
	postCheckHitList(hitlist,&hitlist,segLen);
}

calNewAliSelf3(
    /* aliXY: aligment region on sequence X against sequence Y */
	Region *ali12,		/* region on seq1 aligned against seq2 */
	Region *ali21,		/* region on seq2 aligned against seq1 */
	Region *ali11_1,	/* region on seq1 aligned itself (reg1)*/
	Region *ali11_2,	/* region on seq1 aligned itself (reg2)*/
	Region *brk1, Region *brk2,	/* break points on seq1 and seq2 */
	NewSeqInfo *newseq,
	Region *newali,		/* output: aligned reg on newreg against seq3 */
	char *hitflag,
	Count cnt1, Count cnt2, Node *n1, Node *n2, Edge *e1,
	HitList *hitlist, int swapflag)
{
	int i, j, k;
	int r1, r2, reg1, reg2;
	int *idx = Idx1;
	Region brk1v, brk2v, tmp1, tmp2, tmp3;
	int numh;
	SeqPos from, to, pos1,pos2,prevpos1,prevpos2;
	SeqPos bk1,bk2,prevbk1,prevbk2;
	static PosList *poslist1, *poslist2;
	HitList_t orighit;

	if (swapflag) {
		idx = Idx2;
	}

	initPosList(&poslist1);
	initPosList(&poslist2);
/*
printf("Self3, swap=%d\n",swapflag);
*/

	addPosList(poslist1, brk1, "b1");
	if (overlapCheck(brk1,ali11_2) ==  0) {
		copyReg(&brk1v, brk1); transformReg(&brk1v, ali11_2, ali11_1);
		addPosList(poslist1, &brk1v, "b2");
	}
	addPosList(poslist1, ali11_1, "a1");
	sortPosList(poslist1);

	addPosList(poslist2, brk1, "b1");
	if (overlapCheck(brk1,ali11_1) ==  0) {
		copyReg(&brk2v, brk1); transformReg(&brk2v, ali11_1, ali11_2);
		addPosList(poslist2, &brk2v, "b2");
	}
	addPosList(poslist2, ali11_2, "a2");
	sortPosList(poslist2);

/*
	printPosList(poslist1);
	printf("---\n");
	printPosList(poslist2);
*/

	initHitList(hitlist);
	createOrigHit(&orighit, ali11_1, e1, 1);
	calNewAli_common(poslist1, NULL, idx, NULL,
		n1, NULL, newseq, 1, hitlist, &orighit, segLen,
		FlagHitList_NEW);

	calNewAli_common(poslist2, NULL, idx, NULL,
		n1, NULL, newseq, 2, hitlist, &orighit, segLen,
		FlagHitList_NEW);

	postCheckHitList(hitlist, &orighit, segLen);

/*
	printHitList(hitlist, 3);
*/


/*
for (k = 1; k <= *numhit; k++) {
int l1 = regLen(&(hit[k-1].reg1));
int l2 = regLen(&(hit[k-1].reg2));
	printf(">>%d,%d,%d,%d\n",pos1,prevpos1,bk1,prevbk1);
	printf(">>%d,%d,%d,%d\n",pos2,prevpos2,bk2,prevbk2);
	printf("%d,%d,%d,(%d,%d),(%d,%d),%d,%d\n",
		k-1,hit[k-1].seg1,hit[k-1].seg2,
		hit[k-1].reg1.from,hit[k-1].reg1.to,
		hit[k-1].reg2.from,hit[k-1].reg2.to,bk1,bk2);
assert( ( (l1 > l2) ? ((double) l2 / l1) : ((double) l1 / l2) ) >= 0.7 );

if( ( (l1 > l2) ? ((double) l2 / l1) : ((double) l1 / l2) ) < 0.8 ){
fprintf(stderr, "Warning: unbalanced alignment: %d,%d\n",l1,l2);
}
}
*/

}

calNewAliSelf(
    /* aliXY: aligment region on sequence X against sequence Y */
	Region *ali12,		/* 1st region of seq1 */
	Region *ali21,		/* 2nd region of seq1 */
	Region *ali13,	/* region on seq1 aligned against seq3 */
	Region *ali31,	/* region on seq3 aligned against seq3 */
	Region *brk1, Region *brk2,	/* break points on seq1 and seq2 */
	NewSeqInfo *newseq,
	Region *newali,		/* output: aligned reg on newreg against seq3 */
	Region *newali3,	/* output: aligned reg on seq3 against newreg */
	char *hitflag,
	Count cnt1, Count cnt2, Node *n1, Node *n2, Edge *e1,
	HitList *hitlist)
{
	Region *tmp, brk3_1, brk3_2;
	static PosList *poslist1;
	SeqPos from, to, pos1,pos2,prevpos1,prevpos2;
	SeqPos bk1,bk2,prevbk1,prevbk2;
	int r1, r2, reg1, reg2;
	int i,j,k,flag, seg, prevseg;
	HitList_t *hit = hitlist->hit;
	HitList_t orighit;

	int *idx = Idx1;
	int len1, len2;

	if(  ali12->from < ali21->from) {
	} else {
		tmp = ali21; ali21 = ali12; ali12 = tmp;
		tmp = brk1; brk1 = brk2; brk2 = tmp;
	}
	copyReg(&tmp, brk1); transformReg(&brk3_1, ali13, ali31);
	copyReg(&tmp, brk2); transformReg(&brk3_2, ali13, ali31);

	initPosList(&poslist1);
	if (ali12->to > ali21->from) {
		/** overlapping alignment **/
		idx = Idx2;
	}
	len1 = regLen(ali12); len2 = regLen(ali21);
	len1 = len1 < len2 ? len1 : len2;

/*
	addPosList(poslist1, brk1, "b1");
	addPosList(poslist1, brk2, "b2");
*/
	addPosList(poslist1, ali12, "b1");
	addPosList(poslist1, ali21, "b2");
	addPosList(poslist1, ali13, "a1");
	sortPosList(poslist1);

/*
	printf("--\n");
	printPosList(poslist1);
*/

	initHitList(hitlist);
	createOrigHit(&orighit, ali13, e1, 1);

	calNewAli_common(poslist1, NULL, Idx1, Idx2,
		n1, n2, newseq, 3, hitlist, &orighit, segLen,
		FlagHitList_NEW);
/*
	copyReg(&orighit, ali13);
	orighit.dist = e1->dist;
	orighit.score = e1->score;
	orighit.connect = e1->connect;
*/
	postCheckHitList(hitlist, &orighit, segLen);

	for (i = 0; i < 5; i++){
		hitflag[i] = 0;
	}
	for (i = 0, j = 0; i < hitlist->numhit; i++) {
		if ( hit[i].seg1 == D_THIRD  || hit[i].seg2 == D_THIRD ) {
			printf("ovlpreg\n");
		} else if ( aliLenCheck(&(hit[i].reg1),0,0) ) {
			int seg = hit[i].seg1;
			int offset = hit[i].offset1;
			if (j != i) {
				copyHitList(&hit[j], &hit[i]);
			}
			j++;
			if (seg != 0 && hit[i].seg2 == 0) {
				seg = 0;
				copyReg(&newali[seg], &hit[i].reg2);
				offset = hit[i].offset2;
			} else {
				copyReg(&newali[seg], &hit[i].reg1);
			}
			copyReg(&tmp, &newali[seg]);
			shiftReg(&tmp, offset);
/*
printReg1(&hit[i].reg1);
printReg1(&hit[i].reg2);
				
printf("seg=%d, %d, (%d,%d)\n",seg, &newali[seg],hit[i].reg2.from,hit[i].reg1.to);
printReg1(&newali[seg]);
*/
			transformReg(&tmp, ali13,ali31);
/*
printReg1(ali13);printReg1(ali31);
printReg1(&tmp);
*/
			ovlpReg(&tmp, ali31, &newali3[seg],NULL,1);
/*
printReg1(&newali[seg]); printReg1(&newali3[seg]);
printReg1(n1->newreg);printReg1(&(newseq->seq));
*/
			if (regLen(&(newali[seg])) >= Opt.minlen &&
				regLen(&newali3[seg]) >= Opt.minlen) {
				reg1++;
				hitflag[seg] = 1;
			} else {
				hitflag[seg] = 0;
			}
		} else {
		}
	}
/*
	printHitList(hit,3);
*/
}

/* split the given alignment at given break points and calculate new
	alignment positions */
calNewAli_common(
	PosList *poslist1, PosList *poslist2, int *idx1, int *idx2,
	Node *n1, Node *n2, NewSeqInfo *newseq,
	int segflag,
	HitList *hitlist, HitList_t *orighit, int *segLen,
	FlagHitList flag_hitlist)
{
	int i, j, k;
	int r1, r2;
	char flag;
	Region tmp;
	SeqPos pos1,pos2, prevpos1,prevpos2;
	SeqPos bk1,bk2,prevbk1,prevbk2;
	Region *newreg1, *newreg2;
	Region hitreg1, hitreg2;
	int segnum1, segnum2;
	int ovflag;

	if (flag_hitlist == FlagHitList_NEW) {
		initHitList(hitlist);
	}

	pos1 = pos2 = INFPOS;
	newreg1 = &(n1->newreg);

/*
printf("==##==\n");
printNode(n1);
printNode(n2);
printReg1(&(n1->newreg));
printReg1(&(n2->newreg));
printf("==##==\n");
*/

	if (n1==n2) {
		newreg2 = &(n1->newreg2);
	} else if (n2) {
		newreg2 = &(n2->newreg);
	}
	if (idx1) {
		segnum1 = idx1[0];
	}
	if (idx2) {
		segnum2 = idx2[0];
	}

	r1 = r2 = k = flag = 0;
	bk1 = bk2 = prevbk1 = prevbk2 = 0;

	for (i = 0; i < poslist1->num; i++) {
		pos1 = poslist1->poslist[i].pos;
		if (! definedPos(pos1)) {
			/** assuming the type is 'from' **/
			pos1 = 0;
		} else if (poslist1->poslist[i].type[2]=='f' && pos1 > 0) {
			pos1--;
		}
		if (poslist2) {
			pos2 = poslist2->poslist[i].pos;
			if (! definedPos(pos2)) {
				/** assuming the type is 'from' **/
				pos2 = 0;
			} else if (poslist2->poslist[i].type[2]=='f' && pos2 > 0) {
				pos2--;
			}
		} else {
			pos2 = pos1;
		}

		if (strncmp(poslist1->poslist[i].type, "b1", 2)==0) {
			/* breakpoint */
			r1++;
			bk1 = pos1;
			if (poslist1->poslist[i].type[2]=='f' && bk1 > 0) {
				bk1--;
			}
			segLen[idx1[r1-1]] = bk1 - prevbk1 + 1;
		} else if (poslist1->poslist[i].type[0] == 'a') {
			/* beginning (flag==1) or ending (flag==2) of
				the alignment */
			flag++;
		}
		if (poslist2) {
			if (strncmp(poslist2->poslist[i].type, "b2", 2)==0) {
				/* breakpoint */
				r2++;
				bk2 = pos2;
				if (poslist2->poslist[i].type[2]=='f' && bk2 > 0) {
					bk2--;
				}
			}
		} else if (strncmp(poslist1->poslist[i].type, "b2", 2)==0) {
			r2++;
			bk2 = pos2;
			if (poslist1->poslist[i].type[2]=='f' && bk2 > 0) {
				bk2--;
			}
			if (idx2) {
				segLen[idx2[r2-1]] = bk2 - prevbk2 + 1;
			}
		} else if (! idx1 && idx2) {
			r2 = r1;
			bk2 = bk1;
		}
/*
printf("rr>>>%d,%d, %d,%d  %d\n",r1,segnum1, r2, segnum2, flag);
*/
		if (flag == 0) {
			continue;
		}

		if (k < 1) {
			/* skip the first pos */
			k++;
		} else {
			setSeqReg(&hitreg1,
				prevpos1-prevbk1, pos1-prevbk1);
if (Opt.DEBUG){
printf("%d\n",prevbk1);
printReg1(&hitreg1);
printReg1(newreg1);
printReg1(&newseq->seq);
}
			ovflag = 0;
			if (IsMatchedSeg(segnum1)) {
				shiftReg(&hitreg1, prevbk1);
				ovflag = mapReg(
					&hitreg1, newreg1, &(newseq->seq));
			}
if (Opt.DEBUG){
printf("%d\n",segnum1);
printReg1(&hitreg1);
printReg1(newreg1);
printReg1(&(newseq->seq));
printf("----\n");
}
{ int ll = regLen(&(newseq->seq));
if (IsMatchedSeg(segnum1)&& ovflag == 0 && hitreg1.to > ll){
printf("ll>???:%d,%d,%d,%d,%d\n",pos1,prevpos1,prevbk1,hitreg1.to,ll);
printReg1(&hitreg1);
printReg1(newreg1);
printReg1(&(newseq->seq));
}
}
			if (ovflag != 0) {
				/** no overlap -- skip **/
			} else if (segflag == 1) {
				addHit(hitlist, segnum1, &hitreg1, prevbk1,
					0, NULL, 0, orighit);
			} else if (segflag == 2) {
				addHit(hitlist, 0, NULL, 0,
					segnum1, &hitreg1, prevbk1, orighit);
			} else if (segflag == 3) {
				setSeqReg(&hitreg2,
					prevpos2-prevbk2, pos2-prevbk2);
				ovflag = 0;
				if ( IsMatchedSeg(segnum2) ) {
				   shiftReg(&hitreg2, prevbk2);
/** boundary should be 'newreg' rather than 'break' **/
				   ovflag = mapReg(
					&hitreg2, newreg2, &(newseq->seq));
				}
/*
printReg1(&hitreg2);
printReg1(newreg2);
printReg1(&(newseq->seq));
printf("----\n");
*/
				if (ovflag != 0) {
					/** no overlap -- skip **/
				} else {
				    addHit(hitlist, segnum1, &hitreg1,
						prevbk1, segnum2,
						&hitreg2, prevbk2, orighit);
				}
			}
			k++;
		}

		if (idx1) {
			segnum1 = idx1[r1];
		}
		if (idx2) {
			segnum2 = idx2[r2];
		}
		prevpos1 = pos1+1;
		prevpos2 = pos2+1;
		prevbk1 = bk1;
		prevbk2 = bk2;
		if (flag == 2) break;
	}
	if (! segLen[idx1[r1]]) {
		segLen[idx1[r1]] = n1->len - bk1 + 1;
	}
	if (idx2 && ! segLen[idx2[r2]]) {
		segLen[idx2[r2]] = n2->len - bk2 + 1;
	}
}
aliLenCheck(Region *ali, int len1, int len2)
{
	int cutoff = overlapCutoff(len1, len2);
	return (regLen(ali) >= cutoff);
}
createOrigHit(HitList_t *orighit, Region *ali, Edge *e, int nnum)
{
	copyReg(&orighit->reg1, ali);
	orighit->dist = e->dist;
	orighit->score = e->score;
	orighit->connect = e->connect;
	if (nnum==1) {
		orighit->cnt = e->node1 ? e->node1->cnt : 0;
	} else {
		orighit->cnt = e->node2 ? e->node2->cnt : 0;
	}
}

clear_segLen(int *segLen) {
	int i;
	for (i = 0;  i <= 7; i++) {
		segLen[i] = 0;
	}
}

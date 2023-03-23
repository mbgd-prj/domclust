/*
 * DomClust: Hierarchical Clustering for Orthologous Domain Classification
 * Copyright (c) 2000-2007, Ikuo Uchiyama
 * All rights reserved.
 */

/**
	Domain cutting procedure by Heger and Holm
**/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "domclust.h"
#include "vararray.h"
#include "util.h"
#include <assert.h>
#include "domcut.h"

#define MAXSEQLEN 5000
#define INI_ALINUM 100
#define ALIPOSNUM (MAXNBR*2)

typedef enum {Left, Right} Side;

double Rtab[MAXSEQLEN];
double Stab[MAXSEQLEN];
int pos[MAXNBR];

Region *Align[MAXNBR];
int AliPos[ALIPOSNUM];
int DomSep[MAXNBR];
char domcutFlag[MAXNBR];
double R(), S();
/*
typedef struct {
	SeqPos pos1, pos2;
} PosPair;
*/

int addPosArray(DomCut *domcut, Region *ali1, Region *ali2, char flag);
int calc_score(Region *ali, SeqPos pos, int seqlen, double *Sscore, double *Rscore);

int AliNum, PosNum, DomNum;

int num_cmpr(int *a, int *b) {
	return *a - *b;
}
int ali_cmpr(Region **a, Region **b) {
	int tmp =  (*a)->from - (*b)->from;
	return tmp ? tmp : (*a)->to - (*b)->to;
}
int cmpr_pospair(PosPair *a, PosPair *b) {
	return ((a->pos1 != b->pos1) ?
		(a->pos1 - b->pos1) : (a->pos2 - b->pos2));
}

int value_equal(int a1, int a2)
{
	return (a1 == a2);
}
align_equal(Region *a1, Region *a2)
{
	return (a1->from == a2->from && a1->to == a2->to);
}
copyPosPair(PosPair *pp1, PosPair *pp2)
{
	pp1->pos1 = pp2->pos1;
	pp1->pos2 = pp2->pos2;
}


init_domcut()
{
	static x = 0;
	if (x) return 0;
	make_Rtab();
	make_Stab();
	x = 1;
}


DomCut *createDomCut(NodeSet *nodes)
{
	DomCut *domcut;
	if ( (domcut = (DomCut *) malloc(sizeof(DomCut))) == NULL) {
		fprintf(stderr, "Can't allocate memory\n");
		exit(1); 
	}
	init_domcut();
	domcut->AliPosL = createVarArray(INI_ALINUM, sizeof(PosPair));
	domcut->AliPosR = createVarArray(INI_ALINUM, sizeof(PosPair));
	domcut->nodes = nodes;
	return domcut;
}
initDomCut(DomCut *domcut, Region *ali12, Region *ali21)
{
	clearArray(domcut->AliPosL); clearArray(domcut->AliPosR);
	domcut->maxposL.pos1 = ali12->from; domcut->maxposR.pos1 = ali12->to;
	domcut->maxposL.pos2 = ali21->from; domcut->maxposR.pos2 = ali21->to;
	domcut->ali12 = ali12;
	domcut->ali21 = ali21;
	addPosArray(domcut, ali12, ali21, 3);
	domcut->maxscoreL = domcut->maxscoreR = -1;
}

addPos(DomCut *domcut, Region *ali1, Region *ali2)
{
	char flag = 0;
	if (checkDiff(domcut->ali12->from, ali1->from, Left)
	 && checkDiff(domcut->ali21->from, ali2->from, Left)) {
		flag |= 1;
	}
	if (checkDiff(domcut->ali12->to, ali1->to, Right)
	 && checkDiff(domcut->ali21->to, ali2->to, Right)) {
		flag |= 2;
	}
	addPosArray(domcut, ali1, ali2, flag);
}
addPosArray(DomCut *domcut, Region *ali1, Region *ali2,
				char flag)
{
	PosPair pair;
	if (flag & 1) {
		pair.pos1 = ali1->from;
		pair.pos2 = ali2->from;
		addArray(domcut->AliPosL, &pair);
	}
	if (flag & 2) {
		pair.pos1 = ali1->to;
		pair.pos2 = ali2->to;
		addArray(domcut->AliPosR, &pair);
	}
}
find_maxcuts(DomCut *domcutObj, Node *n1, Node *n2, Node3List *nlist, int nbrnum)
{
	sortArray(domcutObj->AliPosL, cmpr_pospair);
	find_cut(domcutObj->AliPosL, n1->len, n2->len, nlist, nbrnum,
			&(domcutObj->maxscoreL), &(domcutObj->maxposL));

	sortArray(domcutObj->AliPosR, cmpr_pospair);
	find_cut(domcutObj->AliPosR, n1->len, n2->len, nlist, nbrnum,
			&(domcutObj->maxscoreR), &(domcutObj->maxposR));
}
assignNewSeq( DomCut *domcutObj,
	int cutflagL1, int cutflagR1, int cutflagL2, int cutflagR2,
	NewSeqInfo *newseq, Node *n1, Node *n2,
	Region *ali12, Region *ali21)
{
	int newreg_len1, newreg_len2, newlen;
	int newalilen;
	int minlen1, minlen2, maxlen1, maxlen2;
	SeqPos from1, to1, from2, to2;

	from1 = domcutObj->maxposL.pos1; to1 = domcutObj->maxposR.pos1;
	from2 = domcutObj->maxposL.pos2; to2 = domcutObj->maxposR.pos2;

	setSeqReg(&(newseq->boundary1), from1, to1);
	setSeqReg(&(newseq->boundary2), from2, to2);
	newreg_len1 = (int) regLen(&(newseq->boundary1));
	newreg_len2 = (int) regLen(&(newseq->boundary2));

	/** ** creating the new sequence ** **/
	newlen = meanValue(newreg_len1, n1->cnt, newreg_len2, n2->cnt);
	setSeqReg(&(newseq->seq), 1, newlen);

	/** conserved (aligned) region on the new sequence,
		each end of which is used as an anchor point **/
	newalilen = meanValue((int) regLen(ali12), n1->cnt,
		(int) regLen(ali21), n2->cnt);
	newseq->aliM.from = meanValue(
		(int) (ali12->from - newseq->boundary1.from), n1->cnt,
		(int) (ali21->from - newseq->boundary2.from), n2->cnt);
	if (newseq->aliM.from < 0) {
		newalilen += newseq->aliM.from;
		newseq->aliM.from = 0;
	}
	newseq->aliM.to = newseq->aliM.from + newalilen - 1;
if (newseq->aliM.from < 0 || newseq->aliM.to < 0) {
printf("aliM>>\n");
printReg(ali12);
printf("%d,%d,%d,%d\n",from1,to1,from2,to2);
printReg(&(newseq->aliM));
}

        /** determine the length range of the new sequence **/
/*
        minlen1 = n1->minlen * newreg_len1 / n1->len;
        minlen2 = n2->minlen * newreg_len2 / n2->len;
        newseq->minlen = (SeqPos)min(minlen1, minlen2);
        maxlen1 = n1->maxlen * newreg_len1 / n1->len;
        maxlen2 = n2->maxlen * newreg_len2 / n2->len;
        newseq->maxlen = (SeqPos)max(maxlen1, maxlen2);
*/

        /** determine the break point **/
        /*  cut the alignment end only when some other segments are matched
                outside of the current region OR the overhang region
                is sufficiently long as an independent domain [Opt.minlen2] */
        newseq->break1.from = newseq->break2.from = INFPOS;
        newseq->break1.to = newseq->break2.to = SUPPOS;

    if (Opt.nobreak){
	/* do nothing */
    } else {
	if ( (cutflagL1 && Opt.minlen < from1) || Opt.minlen2 < from1) {
		newseq->break1.from = from1;
	}
	if ( (cutflagR1 && Opt.minlen < n1->len - to1 + 1)
			|| Opt.minlen2 < n1->len - to1 + 1) {
		newseq->break1.to = to1;
	}
	if ( (cutflagL2 && Opt.minlen < from2) || Opt.minlen2 < from2) {
		newseq->break2.from = from2;
	}
	if ( (cutflagR2 && Opt.minlen < n2->len - to2 + 1)
			|| Opt.minlen2 < n2->len - to2 + 1) {
		newseq->break2.to = to2;
	}
    }
/*
    printf("%d,%d,%d,%d\n",
	newseq->break1.from,newseq->break1.to,
	newseq->break2.from,newseq->break2.to);
*/

	/* save break points in each child node : ??? */
	copyReg(&(n1->brk), &(newseq->break1));
	copyReg(&(n1->newreg), &(newseq->boundary1));
	if (n1 == n2) {
		/* self */
		copyReg(&(n1->brk2), &(newseq->break2));
		copyReg(&(n1->newreg2), &(newseq->boundary2));
	} else {
		copyReg(&(n2->brk), &(newseq->break2));
		copyReg(&(n2->newreg), &(newseq->boundary2));
	}

}

checkDiff(int pos1, int pos2, Side side)
{
	static int maxDiffIn = 1, maxDiffOut = 500;
	if (side == Left) {
		return checkDiff_sub(pos1, pos2, maxDiffOut, maxDiffIn);
	} else {
		return checkDiff_sub(pos1, pos2, maxDiffIn, maxDiffOut);
	}
}
checkDiff_sub(int pos1, int pos2, int diff1, int diff2)
{
	return (pos1 >= pos2 && pos1 - pos2 <= diff1) ||
		(pos1 < pos2 && pos2 - pos1 <= diff2);
}
find_cut( varArray *AliPos, int seqlen1, int seqlen2,
		Node3List *nlist, int nbrnum,
		double *ret_maxscore, PosPair *ret_maxpos)
{
	arrayIter a_iter;
	PosPair *pp, prevpos;
	Region *ali = 0;
	int i;
	double Sscore, Rscore, score;
	PosPair maxpos;
	double maxscore = -1.0;

	setArrayIter(&a_iter, AliPos);
	prevpos.pos1 = prevpos.pos2 = -1;
	while (pp = (PosPair *) getArrayIter(&a_iter)) {
		if (pp->pos1 == prevpos.pos1 && pp->pos2 == prevpos.pos2) {
			continue;
		}
		Sscore = 0.0, Rscore = 0.0;
		for (i = 0; i < nbrnum; i++) {
			calc_score(nlist[i].ali13, pp->pos1, seqlen1,
					&Sscore, &Rscore);
			calc_score(nlist[i].ali23, pp->pos2, seqlen2,
					&Sscore, &Rscore);

		}
		score = Rscore + Sscore;

if (Opt.DEBUG){
	printf("pp>>%d,%d,%lf\n",pp->pos1,pp->pos2,score);
}
		if (score > maxscore) {
			maxscore = score;
			copyPosPair(&maxpos, pp);
/*
printf("upd:(%d,%d): %lf,%lf  %lf\n",maxpos.pos1,maxpos.pos2,
	Rscore,Sscore,maxscore);
*/
		}

		copyPosPair(&prevpos, pp);
	}
	*ret_maxscore = maxscore;
	if (maxscore > 0) {
		copyPosPair(ret_maxpos, &maxpos);
	}
}

calc_score(Region *ali, SeqPos pos, int seqlen, double *Sscore, double *Rscore)
{
	if (!ali) {
		return 0;
	}
	if (ali->from <= pos && pos <= ali->to) {
		*Sscore += S(pos - ali->from + 1)
			+ S(ali->to - pos + 1)
			- S(ali->to - ali->from + 1);
/*
printf("%d,%d,%d: %lf,%lf,%lf\n",ali->from,pos,ali->to,
	S(pos - ali->from + 1),
	S(ali->to - pos + 1) ,S(ali->to - ali->from + 1));
*/

	} else if (pos <= ali->from) {
		*Rscore += -R(pos);
	} else if (ali->to <= pos) {
		*Rscore += -R(seqlen - pos);
	}
	if (Opt.DEBUG){
	printf("ali>>%d,%d,%d,%lf,%lf\n",pos,ali->from,ali->to,*Sscore,*Rscore);
	}
}



addDom(int pos) {
	int i, tmp;
	for (i = DomNum-1; i >= 0; i--) {
/*
printf("%d: %d,%d\n",i,DomSep[i],pos);
*/
		if (pos > DomSep[i]) {
			DomSep[i+1] = pos;
			break;
		} else {
			DomSep[i+1] = DomSep[i];
		}
	}
	DomNum++;
/*
	for (i = 0; i < DomNum; i++) {
		printf("dom=%d\n",DomSep[i]);
	}
*/
}

/****
domcut(NodeSet *nodes, Edge *edge, Node3List *nlist, int nbrnum)
{
	Edge *e1, *e2;
	Node *n1 = edge->node1, *n2 = edge->node2, *n3;
	Region *ali12 = edge->ali1, *ali21 = edge->ali2;
	Region *ali13 = NULL, *ali23 = NULL;
	int i;
	
	for (i = 0; i < nbrnum; i++) {
		e1 = nlist[i].e1;
		e2 = nlist[i].e2;
		n3 = getNode(nodes, nlist[i].n3);
		if (e1 && e2) {
			printf("%d,%d\n",nlist[i].ali13->from,
				nlist[i].ali13->to);
			printf("%d,%d\n",nlist[i].ali23->from,
				nlist[i].ali23->to);
		}
	}
}
****/

double calc_R(int len) {
	double a = 0.06;
	double a2 = 0.007;
	double b = 0.5;

	return(log(b * a * exp(-a * len) + (1-b) * a2 * exp(-a2 * len)));
}
double R(int len) {
	if (len >= MAXSEQLEN) {
		return calc_R(len);
	}
	return Rtab[len];
}
make_Rtab() {
	int len;
	for (len = 0; len < MAXSEQLEN; len++) {
		Rtab[len] = calc_R(len);
	}
}

double calc_S(int len) {
	double p0 = 0.004;
	double A = 1 - p0;
	double tc = 8.578;
	double w = 76.411;
	double z = ((double)len - tc) / w;
	return(log( 1 - (p0 + A * exp(-exp(-z) - z + 1)) ));
}
double S(int len) {
	if (len >= MAXSEQLEN) {
		return calc_S(len);
	}
	return Stab[len];
}
make_Stab() {
	int len;
	for (len = 0; len < MAXSEQLEN; len++) {
		Stab[len] = calc_S(len);
	}
}


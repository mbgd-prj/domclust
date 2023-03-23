#include <stdio.h>
#include "domclust.h"
#include "plist.h"
#include <sys/types.h>
#include <time.h>
#include "neighbor.h"

typedef struct ScoreElem {
	Edge *edge;
	Dist scoreF, scoreB;
	struct ScoreElem *prev, *next;
} ScoreElem;

static ScoreElem *nbrScores;
static NodeSet *Nodes;

/*
static Alloc_Object *scObj;
*/
Dist gapcheck();

printScel(const ScoreElem *scel)
{
	printf("{%d,%d}",scel->edge->node1->id,scel->edge->node2->id);
}

checkNeighbor2(EdgeSet *edges, NodeSet *nodes)
{
	int i;
	time_t tt, tt2;
	Dist origscore;
	Edge *edge;
	char *nodeflag;
	enum {LEFT =1, RIGHT=2} tmpdir;
	char segflag;

/*
	printf("Start\n");
	time(&tt);
*/

	Nodes = nodes;
	if ( (nbrScores = (ScoreElem *)
		    calloc(edges->edgenum, sizeof(ScoreElem))) == NULL) {
		fprintf(stderr, "Can't alloc memory\n");
		exit(1);
	}

/*
	scObj = init_alloc_object_with_freelist(sizeof(ScoreElem),EDGE_BLKSIZ);
*/

	/* forward */
	checkNeighbor22(edges, nodes, 1);
	/* backward */
	checkNeighbor22(edges, nodes, -1);

	for (i = 0; i < edges->edgenum; i++) {
		origscore = nbrScores[i].edge->score;
/*
printf("%d,%lf,%lf\n",i,nbrScores[i].scoreF,origscore);
*/
		edge = nbrScores[i].edge;
		segflag = 0;

		if (nbrScores[i].scoreF > origscore
			/* bi-directional path check */
			&& nbrScores[i].prev->next == &nbrScores[i]) {

			edge->score +=
			  (nbrScores[i].scoreF - origscore) * Opt.nbrScoreRatio;
			if (Opt.nbrScoreLim) {
				int maxscore = origscore + Opt.nbrScoreLim;
/*
				int maxscore = origscore * Opt.nbrScoreLim;
*/
				if (edge->score > maxscore) {
					edge->score = maxscore;
				}
			}
/*
			pushList(edge->left, nbrScores[i].prev->edge);
*/

/*
			printEdge(nbrScores[i].prev->edge);
			printEdge(edge);
*/

/*
printf("====\n");
printEdge(nbrScores[i].prev->edge);
printEdge(edge);
*/
			addEdgeNeighbor(nbrScores[i].prev->edge, edge);
/*
			if (Opt.sim) {
				edge->dist = edge->score;
			}
*/


/*
printf("F>>%s,%s,%lf\n",nbrScores[i].edge->node1->name,nbrScores[i].edge->node2->name,nbrScores[i].edge->score);
traceBack(&nbrScores[i],1);
*/

			segflag = 1;
		}
/*
else {
printf(">F>>%s,%s,%lf\n",nbrScores[i].edge->node1->name,nbrScores[i].edge->node2->name,nbrScores[i].edge->score);
}
*/
		if (nbrScores[i].scoreB > origscore
			/* bi-directional path check */
			&& nbrScores[i].next->prev == &nbrScores[i]) {
			edge->score +=
			  (nbrScores[i].scoreB - origscore) * Opt.nbrScoreRatio;
			if (Opt.nbrScoreLim) {
				int maxscore = origscore + Opt.nbrScoreLim;
/*
				int maxscore = origscore * Opt.nbrScoreLim;
*/
				if (edge->score > maxscore) {
					edge->score = maxscore;
				}
			}
/*
			pushList(edge->right, nbrScores[i].next->edge);
*/
/***
printf("##%d,%d\n",edge->id, nbrScores[i].next->edge->id);
***/
/*
			addEdgeNeighbor(edge, nbrScores[i].next->edge);
*/
/*
			if (Opt.sim) {
				edge->dist = edge->score;
			}
*/

/*
printf("B>>%s,%s,%lf\n",nbrScores[i].edge->node1->name,nbrScores[i].edge->node2->name,nbrScores[i].edge->score);
traceBack(&nbrScores[i],-1);
*/
			segflag = 1;
		}
		if (Opt.neighbor == 3 && segflag) {
			printf("%s %s\n",edge->node1->name,edge->node2->name);
		}

/*
printf("%lf\n",nbrScores[i].edge->score);
*/
/*
		printf(">>%lf\n",nbrScores[i].score);
*/
	}


	free(nbrScores);
	if (Opt.neighbor == 3) {
		exit(0);
	}



/*
	printf("End\n");
	time(&tt2);
	printf("Time: %d\n", tt2 - tt);
*/

	if (Opt.neighbor == 2) {
		Edge *edge;

		/* print relationships to which neighbor effects are added */
		for (i = 0; i < edges->edgenum; i++) {
			edge = (Edge *) get_objdata_idx(edges->edgeobj, i);
			printf("%s %s %d %d %d %d %.1lf %.1lf\n",
				edge->node1->name, edge->node2->name,
				edge->ali1->from, edge->ali1->to,
				edge->ali2->from, edge->ali2->to,
				edge->dist, edge->score);
		}
		exit(0);
	}

}
traceBack(ScoreElem *scel, int dir)
{
	if (! scel) return 0;
	if (dir > 0) {
		printf("## %lf : ", scel->scoreF); printEdge(scel->edge);
		traceBack(scel->prev,dir);
	} else if (dir > 0) {
		printf("## %lf : ", scel->scoreB); printEdge(scel->edge);
		traceBack(scel->next,dir);
	}
}

ScoreElem *allocNbrScore(EdgeID id)
{
	return &nbrScores[id];
}

checkNeighbor22(EdgeSet *edges, NodeSet *nodes, int dir)
{
	int n, nn;
	listIter *eiter;
	Edge *e;
	pList *sclF, *sclB; 
	listIter sclIterF, sclIterB, *sclIter, sclIter2;
	ScoreElem **scp, *scpF, *scpB, *scelIn, *scp2, *maxscp;
	Dist score, maxscore, gap;
	NodeID en2;
	int matchDir;
	enum status {INS_BEFR, INS_AFTR, SKIP} status;
int flag;

	sclF = create_pList();
	sclB = create_pList();

	for (n = 0; n < nodes->nodenum; n++) {
		if (dir > 0) {
			nn = n;
		} else {
			nn = (nodes->nodenum - 1 - n);
		}
		setListIter(&sclIterF, sclF, 1);
		setListIter(&sclIterB, sclB, 1);
		eiter = NULL;
		e = getEdgeByNode(edges, nn, &eiter);

		if (! e) continue;
		scpF = (ScoreElem *) getListIter(&sclIterF);
		scpB = (ScoreElem *) getListIter(&sclIterB);

		while (e) {
			if (e->node1->id != nn) {
				/* nn == e->node2->id */
				e = getEdgeByNode(edges, nn, &eiter);
				continue;
			}

			matchDir = ((int)e->node1->dir) * ((int)e->node2->dir);
			if (matchDir > 0) {
				scp = &scpF;
				sclIter = &sclIterF;
			} else {
				scp = &scpB;
				sclIter = &sclIterB;
			}

			if (dir > 0) {
				while (*scp && ((int) (*scp)->edge->node1->id
						< nn - Opt.dpWin)) {
					delCurrElem(sclIter);
					(*scp) = (ScoreElem *)
						getListIter(sclIter);
				}
			} else {
				while (*scp && ((int) (*scp)->edge->node1->id
						> nn + Opt.dpWin)) {
					delCurrElem(sclIter);
					(*scp) = (ScoreElem *)
						getListIter(sclIter);
				}
			}

			en2 = e->node2->id;
			if (! *scp || en2 <= (*scp)->edge->node2->id) {
				if (! *scp) {
					status = INS_AFTR;
				} else if (en2 == (*scp)->edge->node2->id) {
					NodeID en1 = e->node1->id;
					if (en1 < (*scp)->edge->node1->id) {
						status = INS_BEFR;
					} else if (en1>(*scp)->edge->node1->id){
						status = SKIP;
					} else {
						int ov;
						ov = ovlpReg(e->ali1,e->ali2,
							NULL,NULL,Opt.minovlp);
						if (ov < 0) {
							status = INS_BEFR;
						} else {
							status = SKIP;
						}
					}
				} else {
					status = INS_BEFR;
				}
			} else {
					status = SKIP;
			}

			if (status == INS_BEFR || status == INS_AFTR) {
				scelIn = allocNbrScore(e->id);
				scelIn->edge = e;
				if (status == INS_BEFR) {
					insBeforeCurrElem(sclIter,scelIn);
				} else {
					insAfterCurrElem(sclIter,scelIn);
				}
				e = getEdgeByNode(edges, nn, &eiter);
			} else {
				*scp = (ScoreElem *) getListIter(sclIter);
			}
		}
		geneDP2(sclF, nn, dir, dir);
		geneDP2(sclB, nn, -dir, dir);
/*
printf("%s,%d,%d\n",getNode(nodes, nn)->name,sclF->numelem,sclB->numelem);
*/

		setListIter(&sclIterF, sclF, 1);
		setListIter(&sclIterB, sclB, 1);

/*
printf("***************\n");
if (sclF) {printList(sclF,printScel); putchar('\n');}
if (sclB) {printList(sclB,printScel); putchar('\n');}
*/
	}
}

geneDP2(pList *scl, NodeID nn, int dirH, int dirV)
{
	listIter sclIter, sclIter2;
	Dist score, maxscore, gap;
	ScoreElem *scp, *scp2, *maxscp;
	Edge *e;
	NodeID en2;
int flag;

	setListIter(&sclIter, scl, dirH);
	while (scp = (ScoreElem *) getListIter(&sclIter)) {
		if (scp->edge->node1->id != nn) continue;
		e = scp->edge; en2 = e->node2->id;

		copyListIter(&sclIter, &sclIter2);
		setListIterRev(&sclIter2);
		maxscore = 0;
		/** discard current elem **/ 
		getListIter(&sclIter2);

flag = 0;
		while (scp2 = (ScoreElem *) getListIter(&sclIter2)) {
/*
printf("({%s,%s},{%s,%s},%.0lf)",e->node1->name,e->node2->name,scp2->edge->node1->name,scp2->edge->node2->name,maxscore);
if (dirH < 0 && dirV < 0){
printEdge(scp2->edge);
}
*/
			if ((dirH > 0 &&
				  scp2->edge->node2->id < en2 - Opt.dpWin)
			 || (dirH < 0 &&
				  scp2->edge->node2->id > en2 + Opt.dpWin)) {
				break;
			}
			if ((gap = gapcheck(scp2->edge, e, dirH * dirV, dirV)) < 0){ 
				continue;
			}
			if (dirV > 0) {
				score = scp2->scoreF - gap;
			} else {
				score = scp2->scoreB - gap;
			}
/*
if (strcmp(scp2->edge->node1->name,"bsu:FLIJ")==0) printf("%lf\n",score);
*/
/*
			score = scp2->score - gap;
*/
			if (score > maxscore) {
				maxscore = score;
				maxscp = scp2;
			}
flag=1;
		}
/*
if(flag)
printf("\n");
*/

		if (dirV > 0) {
			if (maxscore > 0) {
/*
printf("**");
printEdge(e);
printf("%lf,%lf\n",maxscore,e->score);
*/
				scp->scoreF = maxscore * Opt.nbrScoreDecay + e->score;
				scp->prev = maxscp;
			} else {
				scp->scoreF = e->score;
			}
		} else {
			if (maxscore > 0) {
				scp->scoreB = maxscore * Opt.nbrScoreDecay + e->score;
				scp->next = maxscp;
			} else {
				scp->scoreB = e->score;
			}
		}
/*
		scp->score = maxscore + e->score;
*/
	}
}

Dist gapcheck(Edge *e1, Edge *e2, int dirMatch, int dirV)
	/* dirMatch > 0 ... direct match, dirMatch < 0 ... inverse match */
	/* dirV > 0 ... forward DP, dirV < 0 ... backward DP */
{
	int ndiff1, ndiff2;
	Dist gap = 0;
	int ov;
	int flag = 0;
	int len1 = 0, len2 = 0;
	NodeID n1, n2;

	ndiff1 = (int) e1->node1->id - e2->node1->id;
	ndiff2 = (int) e1->node2->id - e2->node2->id;

	if (e1->node1->id < e2->node1->id) {
		len1 = e1->node1->len - e1->ali1->to;
		for (n1 = e1->node1->id + 1; n1 < e2->node1->id; n1++) {
			len1 += getNode(Nodes, n1)->len;
		}
		len1 += e2->ali1->from;
	} else if (e1->node1->id > e2->node1->id) {
		len1 = e2->node1->len - e2->ali1->to;
		for (n1 = e2->node1->id + 1; n1 < e1->node1->id; n1++) {
			len1 += getNode(Nodes, n1)->len;
		}
		len1 += e1->ali1->from;
	} else {
		if (e1->ali1->to < e2->ali1->from) {
			len1 += e2->ali1->from - e1->ali1->to;
		} else if (e2->ali1->to < e1->ali1->from) {
			len1 += e1->ali1->from - e2->ali1->to;
		}
	}
	if (e1->node2->id < e2->node2->id) {
		len2 = e1->node2->len - e1->ali2->to;
		for (n2 = e1->node2->id + 1; n2 < e2->node2->id; n2++) {
			len2 += getNode(Nodes, n2)->len;
		}
		len2 += e2->ali2->from;
	} else if (e1->node2->id > e2->node2->id) {
		len2 = e2->node2->len - e2->ali2->to;
		for (n2 = e2->node2->id + 1; n2 < e1->node2->id; n2++) {
			len2 += getNode(Nodes, n2)->len;
		}
		len2 += e1->ali2->from;
	} else {
		if (e1->ali2->to < e2->ali2->from) {
			len2 += e2->ali2->from - e1->ali2->to;
		} else if (e2->ali2->to < e1->ali2->from) {
			len2 += e1->ali2->from - e2->ali2->to;
		}
	}

	if (! (ndiff1 == 1 && ndiff2 == 1)) {
		gap = Opt.iniGap + (Dist) abs(len1 - len2) * Opt.extGap +
				(Dist) (len1 + len2) * Opt.skipGap;
	}
/*
if (strcmp(e1->node1->name,"bha:BH3780")==0 || strcmp(e1->node2->name,"bha:BH3780")==0  || strcmp(e2->node1->name,"bha:BH3780")==0 || strcmp(e2->node2->name,"bha:BH3780")==0) {
*/

/*
if (strcmp(e1->node1->name,"mth:MTH43")==0 || strcmp(e1->node2->name,"mth:MTH43")==0  || strcmp(e2->node1->name,"mth:MTH43")==0 || strcmp(e2->node2->name,"mth:MTH43")==0) {
printf("**%s,%s,%s,%s,%lf,%d,%d\n",e1->node1->name,e1->node2->name,e2->node1->name,e2->node2->name,gap,len1,len2);
}
*/

	if (dirMatch > 0) {
		if (dirV * ndiff1 <= 0 && dirV * ndiff2 <= 0) {
			flag = 1;
			if (ndiff1 == 0) {
				ov = ovlpReg(e1->ali1,e2->ali1,NULL,NULL,
						Opt.minovlp);
				if (dirV * e1->node1->dir * ov >= 0) {
					flag = 0;
				}
			}
			if (ndiff2 == 0) {
				ov = ovlpReg(e1->ali2,e2->ali2,NULL,NULL,
						Opt.minovlp);
				if (dirV * e1->node2->dir * ov >= 0) {
					flag = 0;
				}
			}
		}
	} else {
		if (dirV * ndiff1 <= 0 && dirV * ndiff2 >= 0) {
			flag = 1;
			if (ndiff1 == 0) {
				ov = ovlpReg(e1->ali1,e2->ali1,NULL,NULL,
						Opt.minovlp);
				if (dirV * e1->node1->dir * ov >= 0) {
					flag = 0;
				}
			}
			if (ndiff2 == 0) {
				ov = ovlpReg(e1->ali2,e2->ali2,NULL,NULL,
						Opt.minovlp);
				if (dirV * e1->node2->dir * ov <= 0) {
					flag = 0;
				}
			}
		}
	}

	if (flag) {
		return gap;
	} else {
		return (Dist) -1;
	}
}

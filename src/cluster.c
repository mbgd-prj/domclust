/*
 * DomClust: Hierarchical Clustering for Orthologous Domain Classification
 * Copyright (c) 2000-2007, Ikuo Uchiyama
 * All rights reserved.
 */
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "domclust.h"
#include "vararray.h"
#include "util.h"
#include "seqreg.h"

#ifdef WITH_CALNEWALI
#include "hitlist.h"
#endif

#ifdef WITH_DOMCUT
#include "domcut.h"
static DomCut *domcutObj;
#endif

#define DUMMY 99999999

/** for checkAliOvlp **/
#define ALI_OVRATIO 0.9
#define ALI_MAXOVDIFF 0.01

#define MAXNEWDOM 5 /* 1 central, 2 left and 2 right nodes */

static int DEBUGFLAG = 0;
static int DEBUGSTAT = 0;

typedef enum {LEFT, RIGHT, UP, DOWN} Direction;


int call_domCluster_core(Edge *bestedge, char **args);
int createNewSideNode(NodeSet *nodes, Edge *best, Node *n, SeqPos brkpnt,
	int newlen, Region *newconsreg,
	NodeFlag nodeflag, Node **newnodes, char truncFlag);

int print_seqpos(SeqPos *datum)
{
	printf("%d\n", *datum);
}
int cmpr_edge_nodeid(Edge *e1, Edge *e2)
{
	if (e1->node1->id != e2->node1->id) {
		return e1->node1->id - e2->node1->id;
	} else{
		return e1->node2->id - e2->node2->id;
	}
}
int cmpr_nodeid(Node *n1, Node *n2)
{
	printf("%s,%s\n",n1->name,n2->name);
	return (n1->id - n2->id);
}

/*****
	DomClust main loop
*****/
domCluster(EdgeSet *edges, NodeSet *nodes)
{
	Edge *bestedge;
	pList *nbre;
	int cnt = 0;
	char *dfs_func_args[3];
	char *DEBUG_entry = Opt.DEBUG_ent;
	if (Opt.neighbor) {
		dfs_func_args[0] = (char *) edges;
		dfs_func_args[1] = (char *) nodes;
		dfs_func_args[2] = NULL;
	}
#ifdef WITH_DOMCUT
	if (Opt.domcut) {
		domcutObj = createDomCut(nodes);
	}
#endif

	while (1) {

		bestedge = getBestEdge(edges);

		if (! bestedge ) {
			break;
		}
/*
		if (isNeighbor(bestedge->node1,bestedge->node2)) {
		}
*/

		if (DEBUG_entry) {
			if (strcmp(bestedge->node1->name,DEBUG_entry)==0
					&& ! isFlankNode(bestedge->node1)) {
				Opt.DEBUG = Opt.DEBUG_val;
			} else if (strcmp(bestedge->node2->name,DEBUG_entry)==0
					&& ! isFlankNode(bestedge->node2)) {
				Opt.DEBUG = Opt.DEBUG_val;
				DEBUG_entry = bestedge->node1->name;
			} else {
				Opt.DEBUG = 0;
			}
		}
		if (Opt.DEBUG & DBG_basic) {
			printf("==========\n");
			printf("Best: "); printEdge(bestedge);
		}
		if (Opt.VERBOSE) {
			if (++cnt % Opt.verbose_step2 == 0) {
				fprintf(stderr, "clustering %d\n", cnt);
			}
		}

		if (BETTER(Opt.cutoff, MEASURE(bestedge))) {
			break;
		}

		if (isNeighbor(bestedge->node1,bestedge->node2)){
			/** discard weak similarity between neighboring nodes
			     to avoid incorporating spurious duplications **/
			if (BETTER(Opt.cutoff2, MEASURE(bestedge))) {
				delEdge(edges, bestedge);
				continue;
			}
		}

		domCluster_core(edges, nodes, bestedge);

#ifdef WITH_NEIGHBOR
		if (Opt.neighbor) {
			edgeNeighbor_df_search(bestedge, 1,
				&call_domCluster_core, dfs_func_args,
				NULL, NULL, NULL, NULL, NULL, NULL);
		}
#endif
	}
	if (Opt.VERBOSE) {
		fprintf(stderr, "Done\n");
	}
	Opt.DEBUG = Opt.DEBUG_val;	/* reset the DEBUG flag */
}

/*****
	DomClust core routine
*****/
domCluster_core(EdgeSet *edges, NodeSet *nodes, Edge *bestedge)
{
	pList oldEdges, newEdges; /* temporary lists maintaining pointers to
					the edges before and after merging
					the two nodes, respectively. */
	pList newEdges2;
	pList elist1, elist2;
	Node3List *nlist;
	NewSeqInfo newseq;
	Node *newnodes[MAXNEWDOM];
	int nbrnum;
	int i;
	int self_match = 0;

 	init_pList(&oldEdges); init_pList(&newEdges); init_pList(&newEdges2);

	if (bestedge->node1 == bestedge->node2) {
		self_match = 1;
	}
#ifdef LARGE
	if (Opt.revMatch) {
		if (bestedge->dir < 0) {
/*
			revStrReg(bestedge->ali2, bestedge->node2->len);
			bestedge->node2->dir *= -1;
*/
		}
	}
#endif

	/* create a list of the third sequences homologous to either
		of the sequences that are connected by the bestedge */
	nbrnum = create_nlist(nodes, edges, bestedge, &nlist);


	/* determine the break points using the nlist */

	checkAliOvlp(nodes, bestedge, nlist, nbrnum, &newseq);

	/* merge the sequences connected by the bestedge
		and split it into domains at the break points */
	createNewNodes(nodes, newnodes, bestedge, &newseq);

	/* update distances(scores) and aligned regions
		between each new node and the third sequences adjacent to it
		that are listed in the nlist */
	updateDist(nodes, edges, nlist, nbrnum, bestedge, newnodes, &newseq,
			&oldEdges, &newEdges, &newEdges2);


	clearFlagNode(bestedge->node1);
	if (bestedge->node1 != bestedge->node2) {
		clearFlagNode(bestedge->node2);
	}

	/* update edge indices */ 
	addNewEdgeIndices(edges, &newEdges);
	addNewEdgeIndices(edges, &newEdges2);
	deleteOldEdges(edges, bestedge, &oldEdges);
	edgeSelectFlag(bestedge);

 	clearList(&oldEdges);
 	clearList(&newEdges);
 	clearList(&newEdges2);
}

/* for depth first search */
call_domCluster_core(Edge *edge, char **args)
{
	EdgeSet *edges = (EdgeSet *) args[0];
	NodeSet *nodes = (NodeSet *) args[1];
	if (! deleted(edge) && ! edgeSelected(edge)) {
		if (Opt.DEBUG & DBG_basic) {
			printf("BestNbr: "); printEdge(edge);
		}
		delBinData(edges->bin, (double)MEASURE(edge), edge->binelem);
		domCluster_core(edges, nodes, edge);
	}
	return 0;
}

/*
             n3
            /  \
           e1  e2
          /      \
	n1--best--n2
*/
int updateDist(NodeSet *nodes, EdgeSet *edges, Node3List *nlist, int nbrnum,
	Edge *best, Node **newnodes, NewSeqInfo *newseq, pList *oldEdges,
	pList *newEdges, pList *newEdges2)
{
	Node *n1, *n2, *n3 = NULL, *n3_2 = NULL, *newn;
	register int i, j;
	int cnt1, cnt2;
	int offset, offset3, offset4;
	Edge *newedge;
	Edge *e1, *e2;
	Edge *newedges[MAXNEWDOM];
	Region *ali13 = NULL, *ali31 = NULL, *ali23 = NULL, *ali32 = NULL;
	Region ali12M, ali21M;
	Region *ali12 = &ali12M, *ali21 = &ali21M;
	Region *brk1 = &(newseq->break1), *brk2 = &(newseq->break2);

	int ovlp;
	Region newali[MAXNEWDOM+1], newali3[MAXNEWDOM+1];
	char hitflag[MAXNEWDOM];
	listIter iter1, iter2;
	EdgeType edgeType;
	NodeID prevn3 = 0;
	Node3List *nlp;

#ifdef WITH_CALNEWALI
	HitList hitlist, hitlist_merged;
	initHitList(&hitlist_merged);
#endif

	n1 = best->node1; n2 = best->node2;
	cnt1 = n1->cnt; cnt2 = n2->cnt;
	StrReg2Reg(ali12, best->ali1);
	StrReg2Reg(ali21, best->ali2);

    	for (j = 0; j < nbrnum; j++) {
		e1 = nlist[j].e1;
		e2 = nlist[j].e2;
		nlp = &(nlist[j]);
/*
		n3 = getNode(nodes, nlist[i].n3);
*/
	
		if (Opt.DEBUG & DBG_cluster) {
			printf("#prevn3=%d, nlist[%d].n3=%d\n", prevn3, j, nlist[j].n3);
			if (! prevn3 || prevn3==DUMMY || prevn3 <= nlist[j].n3) {
			} else {
				printf("assertion failed!!\n");
				fflush(stdout);
			}
			
		}
		assert(! prevn3 || prevn3==DUMMY || prevn3 <= nlist[j].n3);
		prevn3 = nlist[j].n3;
	
		ali13 = ali31 = ali23 = ali32 = NULL;
		n3 = n3_2 = NULL;
		edgeType = NORMAL_EDGE;
	
		for (i = 0; i < 5; i++) {
			hitflag[i] = 0;
		}
		if (e1) {
			ali13 = nlp->ali13;
			ali31 = nlp->ali31;
			if (e1->node1 == n1) {
				n3 = e1->node2;
			} else {
				n3 = e1->node1;
			}
			checkHit(ali13, brk1, n1->len,
				&hitflag[1], &hitflag[0], &hitflag[3], 1);
	
			if (Opt.DEBUG & DBG_basic) {
				printf("Edge1: "); printEdge(e1);
			}
		}
		if (e2) {
			ali23 = nlp->ali23;
			ali32 = nlp->ali32;
			if (e2->node1 == n2) {
				n3_2 = e2->node2;
			} else {
				n3_2 = e2->node1;
			}
			checkHit(ali23, brk2, n2->len,
				&hitflag[2], &hitflag[0], &hitflag[4], 2);
	
			if (! n3) {
				n3 = n3_2;
			} else if (n3 == n2 && n3_2 == n1) {
		       /* a different alignment of the same sequence pair */
	
				edgeType = MULTI_EDGE;
	/*
			} else if (n3 != n3_2) {
				fprintf(stderr, "ERROR: mismatch node3: %d,%d\n",n3->id,n3_2->id);
				exit(1);
	*/
			}
			if (Opt.DEBUG & DBG_basic) {
				printf("Edge2: "); printEdge(e2);
			}
		}

		if (nlist[j].flag && hitflag[0]) {
			/* to avoid spurious overlap */
			hitflag[0] = 0;
		}
	
		if (Opt.DEBUG & DBG_basic) {
			printf("hit=%d,%d,%d,%d,%d\n",
				hitflag[0],hitflag[1],hitflag[2],
				hitflag[3],hitflag[4]);
		}
	
#ifdef WITH_CALNEWALI
	    if (! Opt.NEW) {
#endif
/*
if (ali13 && ! ali31) {
printf("ali13=%d,%d\n",ali13,ali31);
}
if (ali23 && ! ali32) {
printf("ali23=%d,%d\n",ali23,ali32);
}
*/
		calNewAli(ali12, ali21, ali13, ali31, ali23, ali32,
			brk1, brk2,
				newseq, &(newseq->seq), &(newseq->aliM),
				newali, newali3, hitflag,
				n1->cnt, n2->cnt, n3->cnt,
				n1,n2,n3,e1,e2);
		post_calNewAli(n1, n2, n3, n3_2, newnodes,
			edges, best, e1, e2, newedges,
			newali, newali3, hitflag, newEdges);
#ifdef WITH_CALNEWALI
	    } else {
		/* experimental */
		calNewAli_NEW(ali12, ali21, ali13, ali31, ali23, ali32,
			brk1, brk2,
			newseq, &(newseq->seq), &(newseq->aliM),
			newali, newali3, hitflag, &hitlist,
			n1->cnt, n2->cnt, n3->cnt,
			n1,n2,n3,e1,e2, edgeType);
		mergeHitList(&hitlist_merged, &hitlist);
		postCheckHitList(&hitlist_merged,NULL,NULL);

		if (j == nbrnum - 1 || nlist[j+1].n3 != n3->id) {
		    mergeHitSeg(&hitlist_merged, n1, n2);

		    addEdgeFromHitList(&hitlist_merged, edges, newnodes,
		      best, e1, e2, n3, newedges, newEdges, newEdges2);
		    initHitList(&hitlist_merged);
	    	    if (Opt.neighbor) {
			/** add neighboring links to the new edges **/
			addLinksForNewEdge(newedges, e1, e2);
	    	    }
		}
	    }
#endif
	
		/** mark the old edges to be deleted **/
		if (e1 && ! deleted(e1)) {
			setEdgeFlag(e1, EDGE_CLEARED);
			pushList(oldEdges, e1);
		}
		if (e2 && ! deleted(e2)) {
			setEdgeFlag(e2, EDGE_CLEARED);
			pushList(oldEdges, e2);
		}
    	}
}

#ifdef WITH_CALNEWALI
addEdgeFromHitList(HitList *hitlist, EdgeSet *edges, Node **newnodes,
	Edge *best, Edge *e1, Edge *e2, Node *n3,
	Edge **newedges, pList *newEdges, pList *newEdges2)
{
	int i;
	HitList_t *hit = hitlist->hit;
	Dist dist1, dist2, newdist, score1, score2, newscore;
	Count cnt1, cnt2;
	ConnCount conn1, conn2, newconn; 
	signed char newedir; 
	Edge *newedge;
	Node *n1 = best->node1, *n2 = best->node2;
	Node *nn1, *nn2;
	int seg1, seg2;

	cnt1 = n1->cnt; cnt2 = n2->cnt;
	conn1 = e1 ? e1->connect : 0;
	conn2 = e2 ? e2->connect : 0;
	dist1 = e1 ? e1->dist : Opt.missdist;
	dist2 = e2 ? e2->dist : Opt.missdist;
	score1 = e1 ? e1->score : Opt.missscore;
	score2 = e2 ? e2->score : Opt.missscore;

	for (i = 0; i < hitlist->numhit; i++) {
/*
		calNewDist(hit[i].seg1, cnt1, cnt2, dist1, dist2,
			score1, score2,
			conn1, conn2, &newdist, &newscore, &newconn);
*/
		seg1 = getHitSeg(hit[i].seg1);
		seg2 = getHitSeg(hit[i].seg2);

		nn1 = ( (seg1 == D_THIRD) ? n3 : newnodes[seg1]);
		nn2 = ( (seg2 == D_THIRD) ? n3 : newnodes[seg2]);
		newscore = hit[i].score;
		newdist = hit[i].dist;

		if (! nn1 || ! nn2) continue;

		newedge = addEdgeWithScore(edges, nn1, nn2,
				&(hit[i].reg1), &(hit[i].reg2),
				newdist, newscore, newconn, newedir);
		if (newedge == NULL) {
			continue;
		}

		if(Opt.DEBUG & DBG_basic){
			printf("NewEdge: %d\n",newedge->id);
			printEdge(newedge);
			printf("SEG: %d\n",seg1);
		}

		newedges[seg1] = newedge;
		if (seg1 != D_THIRD && seg2 != D_THIRD) {
			addElemSortedListInv(newEdges2, newedge, cmpr_edge_nodeid);
		} else {
			addElemSortedListInv(newEdges, newedge, cmpr_edge_nodeid);
		}
	}
}
#endif

deleteOldEdges(EdgeSet *edges, Edge *bestedge, pList *oldEdges)
{
	Edge *e;

	while (e = (Edge *) shiftList(oldEdges)) {
		if (deleted(e) && e != bestedge) {
			delEdge(edges, e);
		}
	}
}
addNewEdgeIndices(EdgeSet *edges, pList *newEdges)
{
	Edge *e;
	/** newEdges shold be orderd by node IDs of the 3rd node **/
	while (e = (Edge *) shiftList(newEdges)) {
		addEdgeIndex(edges, e);
	}
}

create_nlist(NodeSet *nodes, EdgeSet *edges, Edge *bestedge,
			Node3List **nlist_ret)
	/** nlist: list of the 3rd nodes connected to at least **/
	/**        one of the best pair; ordered by node ids.   **/
	/** Return: the number of nodes just added. **/
{
	Edge *e1 = NULL, *e2 = NULL;
	NodeID node1, node2, node3_1, node3_2, prev1, prev2;
	Node *n1 = bestedge->node1, *n2 = bestedge->node2, *n3;
	Region ali12M, ali21M;
	Region *ali12 = &ali12M, *ali21 = &ali21M;
	listIter *iter1 = NULL, *iter2 = NULL;
	int nm = 0;
	Node3List *nlp;
	int i;
	int dir1, dir2;
	static varArray *nlistArray;
	Node3List *nlist;

	if (nlistArray == NULL) {
		nlistArray = createVarArray(MAXNBR, sizeof(Node3List));
	} else {
		clearArray(nlistArray);
	}

	node1 = bestedge->node1->id; node2 = bestedge->node2->id;
	prev1 = prev2 = DUMMY;
	StrReg2Reg(ali12, bestedge->ali1);
	StrReg2Reg(ali21, bestedge->ali2);

	/* getEdgeByNode() returns each of edges incident to a given node
	   by successive calls; Orders should be according to node IDs */

	e1 = getEdgeByNode(edges, node1, &iter1);
	if (node1 == node2) {
		/* self match */
	} else {
		e2 = getEdgeByNode(edges, node2, &iter2);
	}

	while (e1 && e2) {
		if (e1 == bestedge) {
			/* skip the best edge itself */
			e1 = getEdgeByNode(edges, node1, &iter1);
			continue;
		} else if (e2 == bestedge) {
			/* skip the best edge itself */
			e2 = getEdgeByNode(edges, node2, &iter2);
			continue;
		}
		node3_1 = otherNodeID(e1, node1);
		node3_2 = otherNodeID(e2, node2);
		if (node3_1 == node2) {
			/* duplicated edge */
printf("DUP\n");
printNode(getNode(nodes, node1));
printNode(getNode(nodes, node2));
putchar('\n');
printEdge(e1);
printEdge(e2);
/*
*/
			add_nlist(nlistArray,e1,e1,DUMMY);
			e1 = getEdgeByNode(edges, node1, &iter1);
		} else if (node3_2 == node1) {
			/* duplicated edge: do nothing
				(already processed via node1 */
			e2 = getEdgeByNode(edges, node2, &iter2);
		} else if (node3_1 == node3_2) {
			add_nlist(nlistArray,e1,e2,node3_1);
			prev1 = node3_1; prev2 = node3_2;
			e1 = getEdgeByNode(edges, node1, &iter1);
			e2 = getEdgeByNode(edges, node2, &iter2);
		} else if (node3_1 < node3_2) {
			add_nlist(nlistArray,e1,NULL,node3_1);
			prev1 = node3_1;
			e1 = getEdgeByNode(edges, node1, &iter1);
		} else {
			add_nlist(nlistArray,NULL,e2,node3_2);
			prev2 = node3_2;
			e2 = getEdgeByNode(edges, node2, &iter2);
		}
	}
	while (e1) {
		node3_1 = otherNodeID(e1, node1);
		if (e1 == bestedge) {
			/* skip the best edge itself */
		} else if (node3_1 == node2) {
			/* duplicated edge */

			add_nlist(nlistArray,e1,e1,DUMMY);
/*
		} else if (! (node3_1 == prev1 || node3_1 == node2)) {
*/
		} else {
			add_nlist(nlistArray,e1,NULL,node3_1);
		}
		prev1 = node3_1;
		e1 = getEdgeByNode(edges, node1, &iter1);
	}
	while (e2) {
		node3_2 = otherNodeID(e2, node2);
		if (e2 == bestedge) {
			/* skip the best edge itself */
		} else if (node3_2 == node1) {
			/* duplicated edge: do nothing
				(already processed via node1) */
/*
		} else if (! (node3_2 == prev2 || node3_2 == node1)) {
*/
		} else {
			add_nlist(nlistArray,NULL,e2,node3_2);
		}
		prev2 = node3_2;
		e2 = getEdgeByNode(edges, node2, &iter2);
	}
	nm = arraySize(nlistArray);
	nlist = (Node3List*) getRawArray(nlistArray);

	for (i = 0; i < nm; i++) {
		nlp = &(nlist[i]);
		e1 = nlist[i].e1; e2 = nlist[i].e2;
		if (e1) {
			StrReg2Reg(&(nlp->ali_e11),e1->ali1);
			StrReg2Reg(&(nlp->ali_e12),e1->ali2);
			if (bestedge->node1 == e1->node2) {
				nlp->ali13 = &(nlp->ali_e12);
				nlp->ali31 = &(nlp->ali_e11);
			} else {
				nlp->ali13 = &(nlp->ali_e11);
				nlp->ali31 = &(nlp->ali_e12);
			}
		} else {
			nlist[i].ali13 = nlist[i].ali31 = NULL;
		}
		if (e2) {
			StrReg2Reg(&(nlp->ali_e22), e2->ali2);
			StrReg2Reg(&(nlp->ali_e21), e2->ali1);
			if (bestedge->node2 == e2->node2) {
				nlp->ali23 = &(nlp->ali_e22);
				nlp->ali32 = &(nlp->ali_e21);
			} else {
				nlp->ali23 = &(nlp->ali_e21);
				nlp->ali32 = &(nlp->ali_e22);
			}
		} else {
			nlist[i].ali23 = nlist[i].ali32 = NULL;
		}
		n3 = getNode(nodes, nlist[i].n3);
		nlist[i].flag = 0;

		/** check to avoid spurious duplication **/
		/**
			hom(n1,n2) hom(n1,n3) but not hom(n2,n3),
			AND neighbor(n2,n3),
			AND not overlap(ali12, ali13)
			=> probably n2 and n3 are different domains.
		**/
		if (e1 && ! e2 && (dir1=isNeighbor(n2,n3))) {
			dir2=overlapCheck(nlist[i].ali13, ali12);
			if (dir2 != 0) {
				if (Opt.DEBUG&DBG_basic){
					printf("spurious up:\n");
					printEdge(bestedge);
					printEdge(e1); printNode(n2);
					printf("\n");
				}
				nlist[i].flag = 1;
			}
		} else if (! e1 && e2 && (dir1=isNeighbor(n1,n3))) {
			dir2 = overlapCheck(nlist[i].ali23, ali21);
			if (dir2 != 0) {
				if (Opt.DEBUG&DBG_basic){
					printf("spurious up:\n");
					printEdge(bestedge);
					printEdge(e2); printNode(n1); printf("\n");
				}
				nlist[i].flag = 1;
			}
		}
	}
	*nlist_ret = nlist;
	return nm;
}
add_nlist(varArray *narray, Edge *e1, Edge *e2, NodeID n3id)
{
	Node3List nlp;
	nlp.e1 = e1; nlp.e2 = e2; nlp.n3 = n3id;
	if (Opt.DEBUG && DBG_cluster) {
		printf("node3: %d\n",n3id);
	}
	addArray(narray, &nlp);
}
to_skip_nlist(Node3List *nlp)
{
	return (nlp->flag);
}


/*
	-------->cutL  cutR<----
	extL<--------------->extR

	      <----best---->

	      --------------> diffR
	      ----------------------->
*/

/* determine the break points using the nlist */
checkAliOvlp(NodeSet *nodes, Edge *edge, Node3List *nlist, int nbrnum, NewSeqInfo *newseq)
{
	listIter iter1, iter2;
	Edge *e1, *e2;
	Node *n1 = edge->node1, *n2 = edge->node2;

	Region ali12M, ali21M;
	Region *ali12, *ali21;
	Region *ali13 = NULL, *ali23 = NULL;
	Region *ali31 = NULL, *ali32 = NULL;
	Region trali13, trali23;
	Region brk1, brk2;
	int newreg_len1, newreg_len2;
	SeqPos newlen, newalilen;
	Node *n3;
	char mchflag;
	int ov1, ov2;
	int i;
	
	int cutL1 = 1, cutL2 = 1;
	int cutR1 = n1->len, cutR2 = n2->len;
	int numcutL1, numcutL2, numcutR1, numcutR2;

	int minextL1, minextL2, minextR1, minextR2;
	int maxextL1, maxextL2, maxextR1, maxextR2;


	int flag = 0, updateflag;
	int skipflag, skipflag1, skipflag2, nbrflag;
	int minlen1, minlen2, maxlen1, maxlen2;
	int self_match = 0;
	int tmp_extL1, tmp_extR1, tmp_extL2, tmp_extR2;
	int nobreak_flag = 0;

	if (n1 == n2) {
		self_match = 1;
	}

/*
	ali12 = edge->ali1; ali21 = edge->ali2;
*/
	StrReg2Reg(&ali12M, edge->ali1);
	StrReg2Reg(&ali21M, edge->ali2);
	ali12 = &ali12M; ali21 = &ali21M;

	minextL1 = ali12->from; minextL2 = ali21->from;
	minextR1 = ali12->to; minextR2 = ali21->to;
	maxextL1 = ali12->from; maxextL2 = ali21->from;
	maxextR1 = ali12->to; maxextR2 = ali21->to;

	copyReg(&brk1, ali12);
	copyReg(&brk2, ali21);
#ifdef WITH_DOMCUT
	if (Opt.domcut) {
		initDomCut(domcutObj, ali12, ali21);
	}
#endif
    if (Opt.nobreak) {
	nobreak_flag = 3;
    } else {
	/** partial sequence (metagenome data) should not split non-partial data **/
	if (partialOnlyCluster(n1->spflag)) {
		nobreak_flag |= 2;
	}
	if (partialOnlyCluster(n2->spflag)) {
		nobreak_flag |= 1;
	}
/*
	if (metagenomeOnlyCluster(n1->spflag)) {
		nobreak_flag |= 2;
	}
	if (metagenomeOnlyCluster(n2->spflag)) {
		nobreak_flag |= 1;
	}
*/
/*
if (strcmp(n1->name, "syn:SLL0816")==0) {
	printf("1) %s: %d\n",n2->name, nobreak_flag);
} else if (strcmp(n2->name, "syn:SLL0816")==0) {
	printf("2) %s: %d\n",n1->name, nobreak_flag);
}
*/
    }

    if (Opt.noextend || nobreak_flag==3) {
	/** do nothing **/
    } else {
	while (1) {
	    updateflag = 0;
	    if (Opt.DEBUG & DBG_basic){
		printf("[%d,%d]; [%d,%d]\n",
			brk1.from,brk1.to,brk2.from,brk2.to);
	    }
	    for (i = 0; i < nbrnum; i++) {
		e1 = nlist[i].e1; e2 = nlist[i].e2;
		if (nlist[i].n3 != DUMMY) {
			n3 = getNode(nodes, nlist[i].n3);
			if (partialOnlyCluster(n3->spflag)) {
				continue;
			}
/*
			if (metagenomeOnlyCluster(n3->spflag)) {
				continue;
			}
*/
		}

/* spflg check and cutoff2 check:
	to prevent a probable noise to cut the region */

		if (e1) {
			ali13 = nlist[i].ali13;
			ali31 = nlist[i].ali31;
			ov1 = overlapCheck(ali13, &brk1);
		}
		if (e2) {
			ali23 = nlist[i].ali23;
			ali32 = nlist[i].ali32;
			ov2 = overlapCheck(ali23, &brk2);

		}

		if (e1 && e2) {
			/* both ali13 and ali23 are found */

			if (! ov1 && ! ov2) {
			/* overlap(ali13,ali12) and overlap(ali23,ali21) */

	if ( (ali13->from>=brk1.from || ali23->from>=brk2.from)
	  && (ali13->to<=brk1.to || ali23->to<=brk2.to)) {
		/* cannot extend either direction */
		continue;
	}
	if (! alignmentConsistency(ali12, ali21, ali13, ali31, ali23, ali32)) {
		/* no alignment consistency */
		continue;
	}
			/**  extension **/
if (ali13->to > n1->len) {
fprintf(stderr, "%d-%d, %d, %s(%d)\n", ali13->from, ali13->to, n1->len, n1->name, n1->id);
}
if (ali23->to > n2->len) {
fprintf(stderr, "%d-%d, %d, %s(%d)\n", ali23->from, ali23->to, n2->len, n2->name, n2->id);
}
				assert(ali13->to <= n1->len);
				assert(ali23->to <= n2->len);

				transformReg2(ali13, &trali23, ali12, ali21);
				transformReg2(ali23, &trali13, ali21, ali12);

				tmp_extL1 = max(ali13->from,trali13.from);
				if (minextL1 > tmp_extL1) {
					minextL1 = tmp_extL1; updateflag = 1;
				}
				tmp_extL2 = max(ali23->from,trali23.from);
				if (minextL2 > tmp_extL2) {
					minextL2 = tmp_extL2; updateflag = 1;
				}
				tmp_extR1 = min(ali13->to,trali13.to);
				if (minextR1 < tmp_extR1) {
					minextR1 = tmp_extR1; updateflag = 1;
				}
				tmp_extR2 = min(ali23->to,trali23.to);
				if (minextR2 < tmp_extR2) {
					minextR2 = tmp_extR2; updateflag = 1;
				}

#ifdef WITH_DOMCUT
				if (Opt.domcut) {
					Region tmpr1, tmpr2;
					tmpr1.from = tmp_extL1;
					tmpr1.to = tmp_extR1;
					tmpr2.from = tmp_extL2;
					tmpr2.to = tmp_extR2;
					addPos(domcutObj, &tmpr1, &tmpr2);
				}
#endif
if (Opt.DEBUG & DBG_basic){
	printf("(%d,%d) (%d,%d)\n",
		tmp_extL1,tmp_extR1,tmp_extL2,tmp_extR2);
}

			}
		}
	    } /* the end of the for loop*/

	    if (! updateflag) {
		break;
	    } else {
		brk1.from= minextL1; brk1.to=minextR1;
		brk2.from= minextL2; brk2.to=minextR2;
	    }
	}
    }

#ifdef WITH_DOMCUT
	if (Opt.domcut) {
		find_maxcuts(domcutObj, n1, n2, nlist, nbrnum);

		brk1.from = domcutObj->maxposL.pos1;
		brk1.to = domcutObj->maxposR.pos1;
		brk2.from = domcutObj->maxposL.pos2;
		brk2.to = domcutObj->maxposR.pos2;
	}
#endif

	numcutL1 = numcutL2 = numcutR1 = numcutR2 = 0;
	for (i = 0; i < nbrnum; i++) {
		Dist cutoff2 = Opt.cutoff2;
		e1 = nlist[i].e1; e2 = nlist[i].e2;
		n3 = getNode(nodes, nlist[i].n3);
		updateflag = 0;
		nbrflag = 0;
		skipflag = 0;
		if (Opt.cutoff2_score_per_len && Opt.cutoff2_score_per_len * n3->len > cutoff2) {
			cutoff2 = (Dist) (Opt.cutoff2_score_per_len * n3->len);
		}
		if (e1) {
			ali13 = nlist[i].ali13;
			ali31 = nlist[i].ali31;
			ov1 = overlapCheck(ali13, &brk1);
			skipflag1 = skipflag2 = 0;
			if (ov1 != 0 && isNeighbor(n2, n3)) {
/*
printf("nnn1:%d\n",isNeighbor(n2,n3));
*/
				/* should not merge neighborhood nodes */
			} else if (ov1 != 0) {
				if ( BETTER(Opt.cutoff, MEASURE(e1))) {
					continue;
				}
				if (! covCheck(ali31, n3->len)) {
					skipflag1 = 1;
				}
				if ( cutoff2 > e1->score ) {
					skipflag2 = 1;
				}
				if ( Opt.cutoff2dist &&
						Opt.cutoff2dist < e1->dist ) {
					skipflag2 = 1;
				}
				if (skipflag1 || skipflag2) {
					skipflag = 1;
				}
			}
		}
		if (!skipflag && e2) {
			ali23 = nlist[i].ali23;
			ali32 = nlist[i].ali32;
			ov2 = overlapCheck(ali23, &brk2);
			skipflag1 = skipflag2 = 0;
			if (ov2 != 0 && isNeighbor(n1, n3)) {
/*
printf("nnn2:%d\n",isNeighbor(n1,n3));
*/
				/* should not merge neighborhood nodes */
			} else if (ov2 != 0) {
				if ( BETTER(Opt.cutoff, MEASURE(e2))) {
					continue;
				}
				if (! covCheck(ali32, n3->len)) {
					skipflag1 = 1;
				}
				if ( cutoff2 > e2->score ) {
					skipflag2 = 1;
				}
				if ( Opt.cutoff2dist &&
						Opt.cutoff2dist < e2->dist ) {
					skipflag2 = 1;
				}
				if (skipflag1 || skipflag2) {
					skipflag = 1;
				}
			}
		}
		if (skipflag) {
			continue;
		}
		if (e1 && e2) {
			if (ov1 != ov2) {
				if (ov1 < 0) {
					numcutL1 += n3->cnt;
					if (cutL1 < ali13->to) {
						cutL1 = ali13->to;
						updateflag = 1;
					}
				} else if (ov1 > 0) {
					numcutR1 += n3->cnt;
					if (cutR1 > ali13->from) {
						cutR1 = ali13->from;
						updateflag = 1;
					}
				}
				if (ov2 < 0) {
					numcutL2 += n3->cnt;
					if (cutL2 < ali23->to) {
						cutL2 = ali23->to;
						updateflag = 1;
					}
				} else if (ov2 > 0) {
					numcutR2 += n3->cnt;
					if (cutR2 > ali23->from) {
						cutR2 = ali23->from;
						updateflag = 1;
					}
				}
			}
		} else if (e1) {
			if (ov1 < 0) {
				numcutL1 += n3->cnt;
				if (cutL1 < ali13->to) {
					cutL1 = ali13->to;
					updateflag = 1;
				}
			} else if (ov1 > 0) {
				numcutR1 += n3->cnt;
				if (cutR1 > ali13->from) {
					cutR1 = ali13->from;
					updateflag = 1;
				}
			}
		} else if (e2) {
			if (ov2 < 0) {
				numcutL2 += n3->cnt;
				if (cutL2 < ali23->to) {
					cutL2 = ali23->to;
					updateflag = 1;
				}
			} else if (ov2 > 0) {
				numcutR2 += n3->cnt;
				if (cutR2 > ali23->from) {
					cutR2 = ali23->from;
					updateflag = 1;
				}
			}
		}
if (Opt.DEBUG & DBG_basic) {
	if (updateflag)
		printf("{%d,%d} {%d,%d}\n", cutL1,cutR1,cutL2,cutR2);
}
	}

	/** set a new region on each sequence **/
	copyReg(&newseq->boundary1, &brk1);
	copyReg(&newseq->boundary2, &brk2);
	newreg_len1 = (int) regLen(&(newseq->boundary1));
	newreg_len2 = (int) regLen(&(newseq->boundary2));

	/** determine the break point **/
	/*  cut the alignment end only when some other segments are matched
		outside of the current region OR the overhang region
		is sufficiently long as an independent domain [Opt.minlen2] */

	newseq->break1.from = newseq->break2.from = INFPOS;
	newseq->break1.to = newseq->break2.to = SUPPOS;

    if ( (nobreak_flag==3)
		|| (newreg_len1 < Opt.minlen && newreg_len2 < Opt.minlen)
			) {
		/* do nothing */
	brk1.from = 1; brk1.to = n1->len; brk2.from = 1; brk2.to = n2->len;
    } else {
	int lenL1, lenR1, lenL2, lenR2;

	/* checking the length of overhang regions */
	/* for this check, we consider the boundary of each segment
		is consreg */

	if (nobreak_flag & 1) {
		brk1.from = 1; brk1.to = n1->len;
	} else {
		if (n2->truncFlag & 0b10) {
			/* do not cut N-term */
			brk1.from = 1;
		} else {
			lenL1 = brk1.from - n1->consreg.from;
       			if ( (numcutL1 >= Opt.mincutcnt && Opt.minlen <= lenL1)
					|| Opt.minlen2 <= lenL1 ) {
				newseq->break1.from = brk1.from;
			} else {
				/* do not cut */
				brk1.from = 1;
			}
		}
		if (n2->truncFlag & 0b01) {
			/* do not cut C-term */
			brk1.to = n1->len;
		} else {
			lenR1 = n1->consreg.to - brk1.to;
			if ( (numcutR1 >= Opt.mincutcnt && Opt.minlen <= lenR1)
					|| Opt.minlen2 <= lenR1) {
				newseq->break1.to = brk1.to;
			} else {
				/* do not cut */
				brk1.to = n1->len;
			}
		}
/*
if (strcmp(n1->name, "syn:SLL0816")==0) {
	printf("%s, %d, %d; %d, %d, %d, %d; %d\n", n2->name, brk1.from, brk1.to, numcutL1, numcutR1, lenL1, lenR1, n1->len);
}
*/
	}
	if (nobreak_flag & 2) {
		brk2.from = 1; brk2.to = n2->len;
	} else {
		if (n1->truncFlag & 0b10) {
			/* do not cut N-term */
			brk2.from = 1;
		} else {
			lenL2 = brk2.from - n2->consreg.from;
			if ( (numcutL2 >= Opt.mincutcnt && Opt.minlen <= lenL2)
					|| Opt.minlen2 <= lenL2) {
				newseq->break2.from = brk2.from;
			} else {
				/* do not cut */
				brk2.from = 1;
			}
		}
		if (n1->truncFlag & 0b01) {
			/* do not cut C-term */
			brk2.to = n2->len;
		} else {
			lenR2 = n2->consreg.to - brk2.to;
			if ( (numcutR2 >= Opt.mincutcnt && Opt.minlen <= lenR2) 
					|| Opt.minlen2 <= lenR2) {
				newseq->break2.to = brk2.to;
			} else {
				/* do not cut */
				brk2.to = n2->len;
			}
		}
/*
if (strcmp(n2->name, "syn:SLL0816")==0) {
	printf("%s, %d, %d; %d, %d, %d, %d: %d\n", n1->name, brk2.from, brk2.to, numcutL2, numcutR2, lenL2, lenR2, n2->len);
}
*/
	}
    }


if (! Opt.nobreak) {
  if (newseq->boundary1.to > n1->len){
    fprintf(stderr, "1>%d,%d\n",newseq->boundary1.to, n1->len);
    printf("1>%d,%d\n",newseq->boundary1.to, n1->len);
    printEdge(edge);
    newseq->boundary1.to = n1->len;
    abort();
  }
  if (newseq->boundary2.to > n2->len){
    fprintf(stderr, "2>%d,%d\n",newseq->boundary2.to, n2->len);
    printf("2>%d,%d\n",newseq->boundary2.to, n2->len);
    printEdge(edge);
    newseq->boundary2.to = n2->len;
    abort();
  }
}

	/** aligned segment on the new sequence (aliM1==aliM2),
		each end of which is used as an anchor point **/
	/* the length of alignment is the average length of both seq. */
	newalilen = meanValue((int) regLen(ali12), n1->cnt,
		(int) regLen(ali21), n2->cnt);
	/* we add the longer unaligned segments to both sides of
		aligned region to keep alignment infomation **/
	newseq->aliM.from = 1 + max(ali12->from - brk1.from,
					ali21->from - brk2.from);
	newseq->aliM.to = newseq->aliM.from + newalilen - 1;

     /** ** creating the new sequence ** **/
	/* taking the maximum length to keep the alignment info. */
	newlen = newseq->aliM.to +
		max((int)brk1.to - ali12->to, (int)brk2.to - ali21->to);
	setSeqReg(&(newseq->seq), 1, newlen);

  	newseq->meanlen = meanValue(newreg_len1, n1->cnt, newreg_len2, n2->cnt);
  	newseq->minlen = (SeqPos) min(n1->minlen * newreg_len1 / n1->len,
				n2->minlen * newreg_len2 / n2->len);

	/* consreg is used as anchor points on the newseq */
	newseq->consreg.from = 1 + max(newseq->boundary1.from - brk1.from,
			newseq->boundary2.from - brk2.from);
	newseq->consreg.to = newseq->consreg.from + newseq->meanlen - 1;

	/* save break points in each child node */
	copyReg(&(n1->brk), &(newseq->break1));
	/* newreg is used as anchor points on the old seqs */
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
createNewNodes(NodeSet *nodes, Node **newnodes, Edge *best,
	NewSeqInfo *newseq)
{
	int i;
	Region *break1 = &(newseq->break1), *break2 = &(newseq->break2);

	int newcnt = best->node1->cnt + best->node2->cnt;
	int newlen = regLen(&(newseq->seq));
	int newlen1, newlen2, newlen3, newlen4;
	int newtotlen;
	Region newconsreg;
	int newminlen, newmaxlen;
	specFlag newspflag;
	Node *n1 = best->node1;
	Node *n2 = best->node2;

	newlen = regLen(&(newseq->seq));
	spFlagOR(n1->spflag, n2->spflag, newspflag);

	for (i = 0; i < MAXNEWDOM; i++) {
		newnodes[i] = NULL;
	}

	newlen1 = break1->from - 1;
	newlen2 = break2->from - 1;
	newlen3 = n1->len - break1->to + 1;
	newlen4 = n2->len - break2->to + 1;

	newtotlen = newlen;
	if (newlen1 > 0 && newlen2 > 0) {
		newtotlen += min( newlen1, newlen2 );
	}
	if (newlen3 > 0 && newlen4 > 0) {
		newtotlen += min( newlen3, newlen4 );
	}

	newnodes[0] = addNode(nodes, n1->name, newcnt, newlen,
		&(newseq->consreg), newtotlen,
		best, NODE_INT1|NODE_INT2, newspflag,
		n1->truncFlag & n2->truncFlag); /* trunc if both are trunc */

	newconsreg.from = n1->consreg.from; newconsreg.to = newlen1;
	createNewSideNode(nodes, best, n1, break1->from,
			newlen1, &newconsreg,
			NODE_INT1, &newnodes[1], n1->truncFlag & 0b10); /* reset C-term truncFlag */
	if (n1 != n2) {
		/* we eliminate the case of self match
			because in such a case node2 and node3 are same */
		newconsreg.from = n2->consreg.from; newconsreg.to = newlen2;
		createNewSideNode(nodes, best, n2, break2->from,
			newlen2, &newconsreg,
			NODE_INT2, &newnodes[2], n2->truncFlag & 0b10);
	}
	newconsreg.from = 1; newconsreg.to = n1->consreg.to - break1->to;
	createNewSideNode(nodes, best, n1, break1->to,
			newlen3, &newconsreg,
			NODE_INT1, &newnodes[3], n1->truncFlag & 0b01); /* reset N-term truncFlag */
	newconsreg.from = 1; newconsreg.to = n2->consreg.to - break2->to;
	createNewSideNode(nodes, best, n2, break2->to,
			newlen4, &newconsreg,
			NODE_INT2, &newnodes[4], n2->truncFlag & 0b01);

	addLinksForNewNode(newnodes, best);
if (Opt.DEBUG & DBG_basic){
	for (i = 0; i < MAXNEWDOM; i++) {
		if (newnodes[i]) {
			printf("newnode[%d]=%d\n",i,newnodes[i]->id);
if (i == 0) DEBUGSTAT = newnodes[i]->id;
		}
	}
}
}
createNewSideNode(NodeSet *nodes, Edge *best, Node *n, SeqPos brkpnt,
		int newlen, Region *newconsreg,
		NodeFlag nodeflag, Node **newnodes, char truncFlag)
{
	int newminlen;
	int newmeanlen;
	*newnodes = NULL;
	if (definedPos(brkpnt) && regLen(newconsreg) >= Opt.minlen) {
	    *newnodes = addNode(nodes, n->name, n->cnt,
			newlen,	/* len (max) */
			newconsreg,
			n->totlen,
			best, nodeflag, n->spflag, truncFlag); /* todo: set truncFlag appropriately */
	}
}

/*** defined in util.h
#define min(x, y) ((x) < (y) ? (x) : (y))
#define max(x, y) ((x) > (y) ? (x) : (y))
***/

overlapCheck(Region *reg1, Region *reg2)
{
	int len1 = regLen(reg1);
	int len2 = regLen(reg2);
	int cutofflen = overlapCutoff(len1, len2);

	return ovlpReg(reg1, reg2, NULL, NULL, cutofflen);
}
overlapRegion(Region *reg1, Region *reg2, Region *andreg, Region *orreg)
{
	int len1 = regLen(reg1);
	int len2 = regLen(reg2);
	int cutofflen = overlapCutoff(len1, len2);
	return ovlpReg(reg1, reg2, andreg, orreg, cutofflen);
}
overlapCutoff(int len1, int len2)
{
	register int minlen, maxlen;
	if (len1 < len2) {
		minlen = len1; maxlen = len2;
	} else {
		minlen = len2; maxlen = len1;
	}
	minlen = minlen * Opt.ovlpratio;
	maxlen = maxlen * Opt.ovlpratio2;
	maxlen = max(Opt.minovlp, maxlen);
	/* minimum length */
	minlen = max(MINLEN_LIMIT, minlen);

	return min(minlen, maxlen);

}
/** checking the coverage criterion for each aligned segments */
covCheck(Region *ali12, int len1)
{
	int minlen, maxlen, alilen;
	alilen = regLen(ali12);
	if ( alilen >=  len1 * Opt.coverage) {
		return 1;
	} else {
		return 0;
	}
}

alignmentConsistency(Region *ali12, Region *ali21,
	Region *ali13, Region *ali31, Region *ali23, Region *ali32)
{
	Region trali13;
	int status1, status2;
	int minov = overlapCutoff(regLen(ali31), regLen(ali32));

	transformReg2(ali23, &trali13, ali21, ali12);
	if (cmprOvlpStatusReg(ali13, &trali13, ali31, ali32, minov)) {
		return 1;
	} else {
		return 0;
	}
}

checkHit(Region *reg, Region *brk, int seqlen,
	char *hitL, char *hitM, char *hitR, int n)
{
        int maxfrom, minto;
	Region tmpreg;
	int ovlp_stat, lenL = 0, lenC = 0, lenR = 0;

	if (Opt.DEBUG & DBG_basic ) {
		printf("reg:        "); printReg(reg);
		printf("breakpoint: "); printReg(brk);
	}
	
/**
	if ( (ovlp_stat = ovlpReg(reg, brk, &tmpreg, NULL, 0)) == 0 ) {
		lenC = regLen(&tmpreg);
		if (reg->to > brk->to ) {
			lenR = reg->to - brk->to;
		}
		if (reg->from < brk->from ) {
			lenL = brk->from - reg->from;
		}
	} else if (ovlp_stat > 0) {
		lenR = regLen(reg);
	} else {
		lenL = regLen(reg);
	}
**/

	/* check if reg is overlapped with the left region */
	if (reg->from < brk->from) {
		tmpreg.from = 1;
		tmpreg.to = brk->from - 1;
		if (overlapCheck(&tmpreg, reg) == 0) {
			*hitL = n;
		}
	}
	/* check if reg is overlapped with the middle region */
        tmpreg.from = (int) (brk->from >= 1 ? brk->from : 1);
        tmpreg.to =  (int) (brk->to <= seqlen ? brk->to : seqlen);

	/* hitM = 0 (none), 1 or 2 (one), 3 (both) */
        if (overlapCheck(&tmpreg, reg) == 0) {
		*hitM += n;
	}
	/* check if reg is overlapped with the right region */
	if (brk->to < reg->to) {
		tmpreg.from = brk->to + 1;
		tmpreg.to = seqlen;
		if (overlapCheck(&tmpreg, reg) == 0) {
			*hitR = n;
		}
	}
}

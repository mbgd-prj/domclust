/*
 * DomClust: Hierarchical Clustering for Orthologous Domain Classification
 * Copyright (c) 2000-2007, Ikuo Uchiyama
 * All rights reserved.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "domclust.h"
#include "memalloc.h"
#include "plist.h"
#include "neighbor.h"
#include "vararray.h"
#include "util.h"

#define BLKSIZ 5000
#define isNode(func) (strcmp(func->classname,"Node")==0)
#define fDepth 3	/* findDepth */

static Alloc_Object *nbrObj;

/**** ------------------------------------------ ****/
/***
	base class for generic neighbor objects (node neighbor or edge neighbor)
***/

/** a set of generic object access methods **/
typedef struct {
	char *classname;
	pList * (*getLeft)();
	pList * (*getRight)();
	int (*getID)();
	NodeID *(*getIDptr)();
	void (*setLeft)();
	void (*setRight)();
	char (*getDir)();
	void (*setMark)();
	void (*unsetMark)();
	int (*Marked)();
	void (*setMark2)();
	void (*unsetMark2)();
	int (*Marked2)();
	void (*printObj)();
	int (*delcheck)();
} NbrFunc;

typedef struct NeighborIter {
	char *node;
	listIter *iter;
	NbrFunc *func;
} NeighborIter;


/** public methods **/
Neighbor *createNeighbor(char *node1, char *node2);
int _addNeighbor(char *leftnode, char *rightnode, int dir1, int dir2,
		NbrFunc *func);
int _addNeighborCnt(char *leftnode, char *rightnode, int dir1, int dir2,
		int count, NbrFunc *func);
int _addNeighborRestore(char *leftnode, char *rightnode, int dir1, int dir2,
		int count, char status, NbrFunc *func);
int _mergeNeighborList(char *node1, char *node2, char *newnode, int dir,
		NbrFunc *func);
int _restoreNeighborList(char *node1, char *node2, NbrFunc *func);
char *_getNeighbor(char *node, int dir, int *nbrdir, NbrFunc *func);
void _printNeighborList(pList *list, char *orignode, NbrFunc *func);
Neighbor *_addElemNeighborList(char *node1, char *node2, int dir1, int dir2, int cnt, char status, NbrFunc *func);
int _getMaxNeighbor(pList *list, char *orignode, char nodelCheck,
	char **ret_maxnode, Neighbor **ret_nbr);
int _neighborDeleted(Neighbor *nbr, NbrFunc *nbrf);

/** private methods **/
NeighborIter *createNeighborIter(char *node, int dir, NbrFunc *func);
char *getNeighborIter(NeighborIter *iter, Neighbor **retnbr);
char *_getNeighborIter(NeighborIter *iter, Neighbor **retnbr, char nodelCheck);
void delCurrElemNeighbor(NeighborIter *iter);
char *getNeighborIterNoCheck(NeighborIter *iter, Neighbor **retnbr);

/* ----------------------------------- */
/** internal function and variable for addElemSortedList */
static char *__orig_node;
static NbrFunc *__func;
#define otherNbrNode(nbr) (nbr->node1 != __orig_node ? nbr->node1 : nbr->node2)
int cmpr_nbr(Neighbor *nbr1, Neighbor *nbr2)
{
	char *n1 = otherNbrNode(nbr1);
	char *n2 = otherNbrNode(nbr2);
	if ( __func->getID(n1) == __func->getID(n2) ) {
		return __func->getDir(n1) - __func->getDir(n2);
	} else {
		return __func->getID(n1) - __func->getID(n2);
	}
}
/* ----------------------------------- */


initNeighbor()
{
	nbrObj = init_alloc_object_with_freelist(sizeof(Neighbor), BLKSIZ);
}

Neighbor *createNeighbor0()
{
	Neighbor *nbr = (Neighbor *) memalloc(nbrObj);
	nbr->count = nbr->status = nbr->dir = 0;
	return nbr;
}
Neighbor *createNeighbor(char *node1, char *node2)
{
	Neighbor *nbr;
	if (! nbrObj) initNeighbor();
	nbr = (Neighbor *) memalloc(nbrObj);
	nbr->node1 = node1;
	nbr->node2 = node2;
	nbr->count = 1;
	nbr->status = nbr->dir = 0;
	return nbr;
}

int _addNeighbor(char *leftnode, char *rightnode, int dir1, int dir2,
	NbrFunc *func)
{
	_addNeighborCnt(leftnode, rightnode, dir1, dir2, 1, func);
}
int _addNeighborCnt(char *leftnode, char *rightnode, int dir1, int dir2,
	int count, NbrFunc *func)
{
if(leftnode==rightnode && ((Node*)leftnode)->id==25249){
printf("L=R:");
printNode((Node*)leftnode);
printNode((Node*)rightnode);
abort();
}
	_addNeighborRestore(leftnode, rightnode, dir1, dir2, count, 0, func);
}
int _addNeighborRestore(char *leftnode, char *rightnode, int dir1, int dir2,
	int count, char status, NbrFunc *func)
{
	_addElemNeighborList(leftnode, rightnode, dir1, dir2,
		count, status, func);
}

Neighbor *_addElemNeighborList(char *leftnode, char *rightnode, int dir1, int dir2,
	int count, char status, NbrFunc *func)
{
	Neighbor *nbr, *retnbr;
	pList *np1, *np2;

	nbr = createNeighbor(leftnode, rightnode);
	nbr->count = 0;

#ifdef WITH_DEBUG
/*
if (! isNode(func)){
	func->printObj(leftnode); func->printObj(rightnode);
	printf("%d\n",nbr);
}
*/
#endif

	__func = func;
	__orig_node = leftnode;

	if (dir1 >= 0) {
		if (! func->getRight(leftnode)) {
			func->setRight(leftnode, create_pList());
		}
		np1 = func->getRight(leftnode);
	} else {
		if (! func->getLeft(leftnode)) {
			func->setLeft(leftnode, create_pList());
		}
		np1 = func->getLeft(leftnode);
	}
	retnbr = NULL;
	findElemSortedListInv(np1, nbr, cmpr_nbr, (char**)&retnbr, LIST_INSERT);

	if (retnbr) {
	    /* A neighborhood relationship between leftnode and rightnode
			is already registered:
			increment the count and return */
		retnbr->count += count;
		return retnbr;
	}


	__orig_node = rightnode;
	if (dir2 >= 0) {
		if (! func->getLeft(rightnode)) {
			func->setLeft(rightnode, create_pList());
		}
		np2 = func->getLeft(rightnode);
	} else {
		if (! func->getRight(rightnode)) {
			func->setRight(rightnode, create_pList());
		}
		np2 = func->getRight(rightnode);
	}
	findElemSortedListInv(np2, nbr, cmpr_nbr, (char**)&retnbr, LIST_INSERT);

	nbr->dir = dir1 * dir2;
#ifdef WITH_DEBUG
/*
if (isNode(func)){
	printf("%d: %d,%d %d,%s,%s %d,%d\n",nbr, dir1,dir2, nbr->dir,rightnode->name,leftnode->name,((Node*)rightnode)->id,((Node*)leftnode)->id);
}
*/
#endif
	nbr->count = count;
	nbr->status = status;
	return nbr;
}

int _isNeighbor(char *node1, char *node2, NbrFunc *func)
{
	pList *nbr1L, *nbr1R, *nbr2L, *nbr2R, *nbrL, *nbrR;
	char *orignode, *nbrnode, *nd;
	int num1, num2, dir;
	Neighbor *retnbr;
	NeighborIter iter;

	nbr1L = func->getLeft(node1); nbr1R = func->getRight(node1);
	nbr2L = func->getLeft(node2); nbr2R = func->getRight(node2);
	num1 = numelemList(nbr1L) + numelemList(nbr1R);
	num2 = numelemList(nbr2L) + numelemList(nbr2R);
	if (num1 == 0 && num2 == 0) {
		return 0;
	}
	if (num1 <= num2) {
		orignode = node1; nbrnode = node2; dir = 1;
	} else {
		orignode = node2; nbrnode = node1; dir = -1;
	}
	if (initNeighborIter(&iter, orignode, 1, func) >= 0) {
	    while (nd = getNeighborIterNoCheck(&iter, &retnbr)) {
		if (func->getID(nd) == func->getID(nbrnode)) {
			/* node1->node2 => -1; node2->node1 => 1 */
			return -dir;
		} else if (func->getID(nd) > func->getID(nbrnode)) {
			break;
		}
	    }
	}
	if (initNeighborIter(&iter, orignode, -1, func) >= 0) {
	    while (nd = getNeighborIterNoCheck(&iter, &retnbr)) {
		if (func->getID(nd) == func->getID(nbrnode)) {
			/* node1->node2 => -1; node2->node1 => 1 */
			return dir;
		} else if (func->getID(nd) > func->getID(nbrnode)) {
			break;
		}
	    }
	}
	return 0;
}

_checkOrder(char *orignode, pList *p, NbrFunc *func)
{
	listIter iter;
	Neighbor *nbr;
	char *nbrnode;
	int previd =0;
	setListIter(&iter, p, 1);
	while (nbr = (Neighbor *) getListIter(&iter)) {
		if (orignode == nbr->node1) {
			nbrnode = nbr->node2;
		} else {
			nbrnode = nbr->node1;
		}
		if (previd && previd > func->getID(nbrnode)) {
			return 1;
		}
		previd = func->getID(nbrnode);
	}
	return 0;
}

_mergeNeighborList(char *node1, char *node2, char *newnode, int dir,
		NbrFunc *func)
{
	NeighborIter *iter1, *iter2;
	pList *p1 = NULL, *p2 = NULL, *new = create_pList();
	char *n1 = NULL, *n2 = NULL;
	int dir1, dir2;
	Neighbor *nbr1, *nbr2;
	char flag = 0;

	if (dir < 0) {
		func->setLeft(newnode, new);
	} else {
		func->setRight(newnode, new);
	}
#ifdef WITH_DEBUG
/*
if (isNode(func) && node1 && node2){
	printf("%s,%s,%d,%d\n",node1->name,node2->name,newnode,dir);
}
*/
#endif

	if (node1) {
		if (iter1 = createNeighborIter(node1, dir, func)) {
			n1 = getNeighborIter(iter1, &nbr1);
		}
	}
	if (node2) {
		if (iter2 = createNeighborIter(node2, dir, func)) {
			n2 = getNeighborIter(iter2, &nbr2);
		}
	}

	while (n1 || n2) {
		flag = 1;

		/** ignore mutual nbr relationship b/w node 1 and 2 **/
		if (n1 == node2){
			n1 = getNeighborIter(iter1, &nbr1);
			continue;
		}
		if (n2 == node1){
			n2 = getNeighborIter(iter2, &nbr2);
			continue;
		}

		/** Delete obsolete (=marked) neighbors (delayed process) **/
		if (n1 && nbr1->status) {

			if (Opt.DEBUG & DBG_nbrconst) {
				printf("[##%d]",func->getID(n1));
			}

			n1 = getNeighborIter(iter1, &nbr1);
			continue;
		}
		if (n2 && nbr2->status) {
			if (Opt.DEBUG & DBG_nbrconst) {
				printf("(##%d)",func->getID(n2));
			}
			n2 = getNeighborIter(iter2, &nbr2);
			continue;
		}
		if (n1 && n2 && func->getID(n1) == func->getID(n2)) {
			/** same nodes --- merge **/
			int newcnt = nbr1->count + nbr2->count;

			/** mark these neighbors as obsolete one **/
			nbr1->status = nbr2->status = 1;
			if (nbr1->dir != nbr2->dir) {
				/* Loop!! do not add neighobors **/
			} else if (dir < 0) {
				/* n1 is the left hand side of the newnode */
				_addNeighborCnt(n1,newnode,nbr1->dir,1,newcnt,func);
			} else {
				_addNeighborCnt(newnode,n1,1,nbr1->dir,newcnt,func);
			}

			if (Opt.DEBUG & DBG_nbrconst) {
				printf("[#%d]",func->getID(n1));
			}

			n1 = getNeighborIter(iter1, &nbr1);
			n2 = getNeighborIter(iter2, &nbr2);
		} else if (! n2 || (n1 && func->getID(n1) <= func->getID(n2))){
			if (Opt.DEBUG & DBG_nbrconst) {
				printf("[%d]",func->getID(n1));
			}

			nbr1->status = 1;
			if (dir < 0) {
				_addNeighborCnt(n1,newnode,nbr1->dir,1,nbr1->count,func);
			} else {
				_addNeighborCnt(newnode,n1,1,nbr1->dir,nbr1->count,func);
			}
			n1 = getNeighborIter(iter1, &nbr1);
		} else {
			if (Opt.DEBUG & DBG_nbrconst) {
				printf("(%d)",func->getID(n2));
			}

			nbr2->status = 1;
			if (dir < 0) {
				_addNeighborCnt(n2,newnode,nbr2->dir,1,nbr2->count,func);
			} else {
				_addNeighborCnt(newnode,n2,1,nbr2->dir,nbr2->count,func);
			}

			n2 = getNeighborIter(iter2, &nbr2);
		}
	}
	if ( (Opt.DEBUG & DBG_nbrconst) && flag) {
		putchar('\n');
	}

	if (_checkOrder(newnode, new, func)) {
		func->printObj(newnode);
		_printNeighborList(new,newnode,func);
		printf(": ERROR\n");
	}

}
/*
	Restore neighbors of a child node (node1) with delete flag
		after the parent node (node0) is deleted
*/
_restoreNeighborList(char *node0, char *node1, NbrFunc *func)
{
	Neighbor *nbr;
	listIter iter;
	static char *prev_node0 = NULL;
	static pList *deleteNbr = NULL;
	if (! deleteNbr) {
		deleteNbr = create_pList();
	}

	if (prev_node0 != node0) {
		setListIter(&iter, deleteNbr, 1);
		while (nbr = (Neighbor *) getListIter(&iter)) {
			if (Opt.DEBUG & DBG_nbrrestore){
				printf("==>DEL: ");
				func->printObj(nbr->node1);
				func->printObj(nbr->node2);
				printf(" %d\n",nbr->status);
			}

			nbr->status = 2;
		}
		freeList(deleteNbr);
		deleteNbr = create_pList();
		prev_node0 = node0;
	}
	_restoreNeighborList0(node0, node1, 1, func, deleteNbr);
	_restoreNeighborList0(node0, node1, -1, func, deleteNbr);

/****
	while (nbr = (Neighbor *) shiftList(deleteNbr)) {
		nbr->status = 2;
	}
***/
}

#define NID(n) ((n) ? (((Node*)n)->id) : 0)
#define NSTAT(nbr) ((nbr) ? ((nbr)->status) : -1)

_restoreNeighborList0(char *node0, char *node1, int dir, NbrFunc *func,
	pList *deleteNbr)
{
	NeighborIter *iter1, *iter2;
	char *n1 = NULL, *n2 = NULL;
	Neighbor *nbr1, *nbr2;
	HENTRY ent;
	Hash *hash = NULL;

/* node0: old root to be deleted */
/* node1: new root */
/* node0<->n1, node1<->n2 */

/** dir > 0 :  node(01) -> nbr ; dir < 0 :  nbr <- node(01) **/


	if (Opt.DEBUG & DBG_nbrrestore){
		printf("### %d,%d,%d\n",func->getID(node0),func->getID(node1),dir);
		printf("L=%d,",numelemList( ((Node*)node1)->left ) );
		printf("R=%d\n",numelemList( ((Node*)node1)->right ) );
	}

	if (node0 && (iter1 = createNeighborIter(node0, dir, func))) {
		n1 = getNeighborIterNoCheck(iter1, &nbr1);
	} else {
		return 0;
	}
	if (node1 && (iter2 = createNeighborIter(node1, dir, func))) {
		n2 = getNeighborIterNoCheck(iter2, &nbr2);
	} else {
		return 0;
	}
	while (n1 || n2) {
		if (Opt.DEBUG & DBG_nbrrestore){
			printf("n>%d,%d,%d,%d\n",NID(n1),NID(n2),NSTAT(nbr1),NSTAT(nbr2));
		}

		if (n1 && n2 && func->getID(n1) == func->getID(n2)
				&& nbr1->dir == nbr2->dir) {
			/* a node adjacent to both node0 and node1 */
			if (nbr1->status != 1) {
				Node *c;
				/* CASE I: nbr1 is not deleted */
				nbr2->status = 0;

		if (Opt.DEBUG & DBG_nbrrestore){
			printf(">>>>restore (case I)");
			func->printObj(nbr2->node1); printf(" ");
			func->printObj(nbr2->node2); printf("\n");
		}

				pushList(deleteNbr, nbr1);

			} else if (nbr2->status == 1) {
				/* CASE II: both nbr1 and nbr2 are deleted */
				if (hash == NULL) {
					hash = Hcreate(27449);
				}
				/* save the neighbor of the parent */
				ent.key = (char *) (func->getID(n1));
				/* to be used for counting*/
				ent.datum = (char *) nbr2;
		if (Opt.DEBUG & DBG_nbrrestore){
			printf("enter>%p\n",ent.key);
		}

				HIsearch(hash, &ent, ENTER);

		if (Opt.DEBUG & DBG_nbrrestore){
			printf(">>>>restore (case II) store: ");
			func->printObj(nbr2->node1); printf(" ");
			func->printObj(nbr2->node2); printf("\n");
		}
			}
			n1 = getNeighborIterNoCheck(iter1, &nbr1);
			n2 = getNeighborIterNoCheck(iter2, &nbr2);
		} else {
			if (! n2 || (n1 && func->getID(n1) <= func->getID(n2))){
			    int cnt = 0;
			    if (hash && nbr1->status != 2) {
				/* CASE II: add new link */
				Neighbor *nbr0 = NULL;
				char *child[2];
				int i;
		if (Opt.DEBUG & DBG_nbrrestore){
			printf("CHECK: %d\n",((Node*)n1)->id);
		}
				if (cnt = hashFindInTree(hash,(Node*)n1,
						nbr1->dir,fDepth)) {
		/* add a new relationship node1(child)<->n1(parent_neighbor) */
		/* direction is same as nbr1->dir */

		if (Opt.DEBUG & DBG_nbrrestore){
			printf("new neighbor: %d:%d %d\n",func->getID(node1), func->getID(n1), cnt);
		}
					if (dir > 0) {
						nbr0 = _addElemNeighborList(
					 	    node1, n1, 1, nbr1->dir,
						    cnt,nbr1->status,func);
					} else {
						nbr0 = _addElemNeighborList(
						    n1, node1, nbr1->dir, 1,
						    cnt,nbr1->status,func);
					}
					ent.key = (char *) (func->getID(n1));
					ent.datum = (char *) nbr0;

					/* assert( nbr0 != 0 ); */

					HIsearch(hash, &ent, ENTER);

		if (Opt.DEBUG & DBG_nbrrestore){
			printf(">>>>restore (case II) add: ");
			func->printObj(node1); printf(" ");
			func->printObj(n1); printf("\n");
		}
				}
				pushList(deleteNbr, nbr1);
			    } else {	/* ! hash */
				if (Opt.DEBUG & DBG_nbrrestore){
					printf("no hash\n");
				}
			    }
			    n1 = getNeighborIterNoCheck(iter1, &nbr1);
			} else {
			    n2 = getNeighborIterNoCheck(iter2, &nbr2);
			}
		}
	}
	if (hash) {
		Hdestroy(hash);
	}
	if (Opt.DEBUG & DBG_nbrrestore){
		printf("Done\n");
	}
}

findTreeNeighbor(Node *n1, Node *n2, int dir1, int nbrdir, int count,
		NbrFunc *func)
{
	NeighborIter *iter;
	Neighbor *nbr;
	Node *nn;
	Hash *h;
	HENTRY ent;
	int cnt = 0;

 	iter = (NeighborIter *) createNeighborIter((char*)n1, dir1, func);

	h = Hcreate(27449);
	while (nn = (Node *) getNeighborIter(iter, &nbr)) {
		ent.key = (char *) n1->id;
		ent.datum = (char *) nbr; /* to be used for counting*/
		HIsearch(h, &ent, ENTER);
	}
	cnt = hashFindInTree(h, n2, nbrdir, fDepth);
	if (cnt){
		if (dir1 > 0) {
			_addElemNeighborList((char*)n1, (char*)n2, 1, nbrdir, cnt, 0, func);
		} else {
			_addElemNeighborList((char*)n2, (char*)n1, 1, nbrdir, cnt, 0, func);
		}
	}
	free(iter);
	Hdestroy(h);
}
hashFindInTree(Hash *hash, Node *node, int dir, int depth)
{
	HENTRY ent, *e;
	Neighbor *nbr0;
	int cnt = 0;
	HashIter iter;

	if (node == NULL) return 0;
	if (depth <= 0) return 0;

/*
	initHashIter(&iter, hash);
	while (e = nextHashIter(&iter)) {
		nbr0 = (Neighbor *) e->datum;
		printf("(%d,%d,%d)",(int)e->key,
			((Node*)nbr0->node1)->id,((Node*)nbr0->node2)->id);
	}
	printf("\n");
*/

	ent.key = (char *) (node->id);
	if (HIsearch(hash, &ent, FIND) > 0) {
		nbr0 = (Neighbor *) ent.datum;

		if (Opt.DEBUG & DBG_nbrrestore){
			printf("Found: "); printNode(node); printf("\n");
		}

		if (dir == nbr0->dir) {
			cnt += nbr0->count;
		}
	} else if (node->child) {
		cnt += hashFindInTree(hash, node->child->node1,dir,depth-1);
		cnt += hashFindInTree(hash, node->child->node2,dir,depth-1);
	}
	return cnt;
}

char *_getNeighbor(char *node, int dir, int *nbrdir, NbrFunc *func)
{
	Neighbor *nbr = NULL;
	if (dir < 0) {
		pList *lnode = (pList *) func->getLeft(node);
		if (lnode) {
			nbr = (Neighbor *) getListIdx(lnode, 0);
		}
	} else {
		pList *rnode = (pList *) func->getRight(node);
		if (rnode) {
			nbr = (Neighbor *) getListIdx(rnode, 0);
		}
	}
	if (nbr) {
		if (nbrdir) {
			*nbrdir = dir * nbr->dir;
		}
		if (nbr->node1 == node) {
			return nbr->node2;
		} else {
			return nbr->node1;
		}
	} else {
		return NULL;
	}
}
#ifdef WITH_NEIGHBOR
int _df_search(char *node, int dir, NbrFunc *nbrf,
		int (*func_before)(), char **args_before,
		int (*func_after)(), char **args_after,
		int (*func_nodecut)(), char **args_nodecut,
		int (*func_nbrcut)(), char **args_nbrcut)
{
	pList *inilist = create_pList();
	pList *visit = create_pList();
	listIter *iter;
	char *n;
	_df_search0(node, -dir, nbrf, inilist, NULL,NULL,NULL,NULL,
			func_nodecut,args_nodecut,func_nbrcut,args_nbrcut);
	iter = createListIter(inilist, 1);
	while (n = (char *) getListIter(iter)) {
		nbrf->unsetMark(n);
	}

	setListIter(iter, inilist, 1);

	while (n = (char *) getListIter(iter)) {
		_df_search0(n, dir, nbrf, visit,
			func_before, args_before, func_after, args_after,
			func_nodecut,args_nodecut,func_nbrcut, args_nbrcut);
	}

	iter = createListIter(visit, 1);
	while (n = (char *) getListIter(iter)) {
		nbrf->unsetMark(n);
	}
	freeList(visit); free(visit);
	freeList(inilist); free(inilist);
}
int _df_search0(char *node, int dir, NbrFunc *nbrf, pList *visit,
		int (*func_before)(), char **args_before,
		int (*func_after)(),  char **args_after,
		int (*func_nodecut)(),  char **args_nodecut,
		int (*func_nbrcut)(),  char **args_nbrcut)
{
	pList *pl;
	listIter *iter;
	Neighbor *nbr;
	int nextdir;
	char *nbrnode;
	int ee; /*for debug */

	if (! node) {
		return 0;
	} else if (nbrf->Marked(node)) {
		return 0;
	}
	nbrf->setMark(node);
	pushList(visit, node);
/*
if (isNode(nbrf)){
printList(visit,printNode);
printNode(node);
printf("OK\n");
}
*/

	if (dir < 0) {
		pl = nbrf->getLeft(node);
	} else {
		pl = nbrf->getRight(node);
	}
	if (!pl) {
		return 0;
	}
/*
if (isNode(nbrf)){
printf("##\n");
printNode(node); putchar('\n');
}
*/

/*
if (strcmp(nbrf->classname,"Edge")==0){
nbrf->printObj(node);
}
*/
	iter = createListIter(pl, 1);
	if (func_nodecut) {
		int ret;
		ret = (args_nodecut ? (*func_nodecut)(node, args_nodecut)
			: (*func_nodecut)(node));
		if (ret == -1) return -1;
	}
/*
if (isNode(nbrf)){
printf("NodeCutOK\n");
}
*/
	if (func_before) {
		args_before ? (*func_before)(node, args_before) : (*func_before)(node);
	}
	while (nbr = (Neighbor *) getListIter(iter)) {
		int ret;
/*
		if (nbr->status) {
			continue;
		}
*/

/*
if (isNode(nbrf)){
nbrf->printObj(node); putchar('\n');
}
*/
/*
		if (_neighborDeleted(nbr,nbrf)==2) {
*/
/*
		if ((!isNode(nbrf) && (ee=_neighborDeleted(nbr,nbrf))) ||
			(isNode(nbrf) && _neighborDeleted(nbr,nbrf)==2)) {
*/

		if (_neighborDeleted(nbr,nbrf)) {
/*
if (! isNode(nbrf)){
printf("Del: ");
nbrf->printObj(nbr->node1); putchar(' ');
nbrf->printObj(nbr->node2); printf("%d,%d\n",nbr,ee);
}
*/
			continue;
		}
/*
if (! isNode(nbrf)){
printf("OK: ");
nbrf->printObj(nbr->node1); putchar(' ');
nbrf->printObj(nbr->node2); printf("%d\n",nbr);
}
*/
		if (func_nbrcut) {
			ret = (args_nbrcut ?
				(*func_nbrcut)(nbr, node, args_nbrcut)
					: (*func_nbrcut)(nbr, node));
			if (ret == -1) {
				continue;
			}
		}
/*
if (isNode(nbrf)){
printf("NbrCutOK\n");
}
*/
		nextdir = dir * nbr->dir;
if (Opt.DEBUG == 2) {
    if (! isNode(nbrf) && dir > 0) {
	nbrf->printObj(nbr->node1);
	nbrf->printObj(nbr->node2);
	printf(">>>%d\n",nextdir);
    }
}
/*
if (isNode(nbrf)) {
printf("%s\n",nbrf->classname);
printf("DD>%d,%d\n",dir,nbr->dir);
printf(">%d,%s,%s\n",nbr->dir,nbr->node1->name,nbr->node2->name);
nbrf->printObj(nbr->node1);
nbrf->printObj(nbr->node2);
nbrf->printObj(node);
}
*/
		if (nbr->node1 == node && nbr->node2 != node) {
			nbrnode = nbr->node2;
		} else if (nbr->node2 == node && nbr->node1 != node) {
			nbrnode = nbr->node1;
		} else {
			continue;
		}
/*
if (isNode(nbrf)) {
	printf("****\n");
	printNode(nbrnode);
	printf("****\n");
}
*/
		_df_search0(nbrnode, nextdir, nbrf, visit,
			func_before, args_before, func_after, args_after,
			func_nodecut,args_nodecut,func_nbrcut, args_nbrcut);
	}
	if (func_after) {
		args_after ? (*func_after)(node, args_after)
			: (*func_after)(node);
	}
	return 0;
}
#endif
#ifdef WITH_NEIGHBOROUT
Hash *gaphash;
varArray *countArray;
pList *gapsrch_nbrlist;
pList *gapsrch_visit;
int gap_search(char *node, int dir, NbrFunc *nbrf,
	int *out_maxid, int *out_maxcnt, int *out_nextdir)
{
	pList *visit;
	listIter *iter;
	pList *nbrpl;
	Neighbor *nbr;
	int cnt, maxcnt;
	int id, maxid;
	int nextdir;
	HENTRY entry;
	char *nn;

	if ( (gaphash = Hcreate(200000)) == NULL ) {
		exit(0);
	}
	countArray = createVarArray(2000, sizeof(int));
	gapsrch_nbrlist = create_pList();

	gap_search0(node, dir, nbrf, 0);

	maxid = 0;
	maxcnt = 0;
	iter = createListIter(gapsrch_nbrlist, 1);
	while (nn = (char *) getListIter(iter)) {
		id =  nbrf->getID(nn);
		entry.key = (char *) id;
/*
printf("ID>%d\n",id);
*/
		if (HIsearch(gaphash, &entry, FIND)) {
			cnt = * ((int*)entry.datum);
/*
printf("CNT>>>%d\n",cnt);
*/
			if (maxcnt < cnt) {
				maxcnt = cnt;
				maxid = id;
				nextdir = 1;
			}
		} else {
			entry.key = (char *) (- id);
/*
printf("id>>>%d\n",entry.key);
printf("Find>>%d,%d\n", entry.key,HIsearch(gaphash, &entry, FIND));
*/
			if (HIsearch(gaphash, &entry, FIND)) {
				cnt = * ((int*)entry.datum);
/*
printf("CNT2>>>%d\n",cnt);
*/
				if (maxcnt < cnt) {
					maxcnt = cnt;
					maxid = id;
					nextdir = -1;
				}
			}
		}
	}
/*
printf("=====\n");
*/
	freeListIter(iter);
	freeList(gapsrch_nbrlist);
	free(gapsrch_nbrlist);
	freeArray(countArray);
	Hdestroy(gaphash);
	*out_maxcnt = maxcnt;
	*out_maxid = maxid;
	*out_nextdir = nextdir;
	return maxid;
}
int gap_search0(char *node, int dir, NbrFunc *nbrf, int pathlen)
{
	HENTRY entry;
	Neighbor *nbr;
	int nextdir;
	char *nbrnode;
	pList *nbrpl, *pl;
	int cnt;
	listIter *iter;
	char *n;
	int eee;

	if (pathlen > Opt.nbrConnGap) {
		return 0;
	}
	if (! node) {
		return 0;
	} else if (nbrf->Marked(node)) {
		return 0;
	}
	if (pathlen > 0) {
		nbrf->setMark(node);
		pushList(gapsrch_visit, node);
	}
	if (dir < 0) {
		nbrpl = nbrf->getLeft(node);
	} else {
		nbrpl = nbrf->getRight(node);
	}
	if (!nbrpl) {
		return 0;
	}
	iter = createListIter(nbrpl, 1);
	while (nbr = (Neighbor *) getListIter(iter)) {
/*
		if (_neighborDeleted(nbr,nbrf)==2) {
			continue;
		}
*/
/*
nbrf->printObj(nbr->node1);putchar(' ');
nbrf->printObj(nbr->node2);putchar('\n');
*/
		if ((eee=_neighborDeleted(nbr,nbrf))==2) {
/*
printf("Cont: %d\n", eee);
*/
			continue;
		}
/*
printf("ok\n");
*/

		nextdir = dir * nbr->dir;
		if (nbr->node1 == node && nbr->node2 != node) {
			nbrnode = nbr->node2;
		} else if (nbr->node2 == node && nbr->node1 != node) {
			nbrnode = nbr->node1;
		} else {
			continue;
		}
		entry.key = (char *) (nbrf->getID(nbrnode) * nextdir);
/*
printf(">%d,%d,%d\n",nbrf->getID(nbrnode),entry.key,nextdir);
*/
/*
		nbrf->printObj(nbrnode);
*/
		if (HIsearch(gaphash, &entry, FIND)) {
/*
			pl = (pList*)entry.datum;
			pushList(pl, nbrf->getIDptr(nbrnode));
*/
			cnt = *((int*)entry.datum); cnt+=nbr->count;
/*
printf(">%d,%d\n",*((int*)entry.datum),cnt);
*/
			*((int*)entry.datum) = cnt;
		} else {
/*
			pl = create_pList();
			pushList(pl, nbrf->getIDptr(nbrnode));
			entry.datum = (char *) pl;
*/
			cnt = nbr->count;
			addArray(countArray, &cnt);
			entry.datum = getArrayItem(countArray, countArray->num-1);
/*
printf("II>%d,%d,%d,%d\n",*((int*)entry.datum),entry.key,pathlen,cnt);
*/
			HIsearch(gaphash, &entry, ENTER);
/*
printf("Find>>%d,%d\n", entry.key,HIsearch(gaphash, &entry, FIND));
*/
		}
		if (pathlen == 0) {
			gapsrch_visit = create_pList();
			pushList(gapsrch_nbrlist, nbrnode);
		}

		gap_search0(nbrnode, nextdir, nbrf, pathlen+1);

		if (pathlen == 0) {
			listIter *iter0 = createListIter(gapsrch_visit, 1);
			while (n = (char *) getListIter(iter0)) {
/*
printf(">>>>%d\n",n);
*/
				nbrf->unsetMark(n);
			}
			freeList(gapsrch_visit); free(gapsrch_visit);
			gapsrch_visit = NULL;
			freeListIter(iter0);
		}
	}
	freeListIter(iter);
}
#endif

#ifdef WITH_NEIGHBOR
/* for topological sorting */
typedef struct {
	pList *list;
} topoVis;
_toposort_func(char *node, char **args)
{
	topoVis *tvis = (topoVis *) args[0];
	pushList(tvis->list, node);
}
pList *_toposort(char *node, NbrFunc *nbrf,
		int (*func_nodecut)(), int (*func_nbrcut)())
{
	topoVis topo;
	char *toposort_func_args[2];
	toposort_func_args[0] = (char*) &topo;
	toposort_func_args[1] = NULL;
	topo.list = create_pList();
	_df_search(node, 1, nbrf, NULL, NULL, _toposort_func,
		toposort_func_args, func_nodecut, NULL, func_nbrcut, NULL);
	return topo.list;
}
#endif

NeighborIter *createNeighborIter(char *node, int dir, NbrFunc *func)
{
	NeighborIter *nbriter;
	if ( (nbriter = (NeighborIter *) malloc(sizeof(NeighborIter))) == NULL ) {
		fprintf(stderr, "Can't alloc memory\n");
		exit(1);
	}
	if (initNeighborIter(nbriter, node, dir, func) < 0) {
		free(nbriter);
		return NULL;
	}
	return nbriter;
}
initNeighborIter(NeighborIter *nbriter, char *node, int dir, NbrFunc *func)
{
	pList *nbrlist;
	if (dir < 0) {
		if (! (nbrlist = func->getLeft(node)) ) {
			return -1;
		}
	} else {
		if (! (nbrlist = func->getRight(node)) ) {
			return -1;
		}
	}
	nbriter->node = node;
	nbriter->iter = createListIter(nbrlist, 1);
	nbriter->func = func;
	return 0;
}
char *getNeighborIter(NeighborIter *iter, Neighbor **retnbr)
{
	return _getNeighborIter(iter, retnbr, 1);
}
char *getNeighborIterNoCheck(NeighborIter *iter, Neighbor **retnbr)
{
	return _getNeighborIter(iter, retnbr, 0);
}
char *_getNeighborIter(NeighborIter *iter, Neighbor **retnbr, char nodelCheck)
{
	Neighbor *nbr = NULL;
	*retnbr = NULL;
 	while (nbr = (Neighbor *) getListIter(iter->iter)) {
		*retnbr = nbr;
		if (nodelCheck && _neighborDeleted(nbr, iter->func)) {
			continue;
		} else if (nbr->node1 == iter->node) {
			return nbr->node2;
		} else {
			return nbr->node1;
		}
	}
	return NULL;
}
void delCurrElemNeighbor(NeighborIter *iter)
{
	if (iter->iter->ptr) {
		((Neighbor *)iter->iter->ptr->prev->datum)->status = 1;
	} else {
		((Neighbor *)iter->iter->plist->last->datum)->status = 1;
	}
}


void _printNeighborList(pList *list, char *orignode, NbrFunc *func)
{
	listIter iter;
	Neighbor *nbr;
	char *node2;
	int i = 0;

	setListIter(&iter, list, 1);
	putchar('[');
	while (nbr = (Neighbor *) getListIter(&iter)) {
		if (nbr->node1 == orignode) {
			node2 = nbr->node2;
		} else {
			node2 = nbr->node1;
		}

/*
if (nbr->count != orignode->cnt) continue;
*/
		if (i++) putchar(',');
if (nbr->status) printf("#");
		func->printObj(node2);
		printf("<%d>",nbr->count);
	}
	putchar(']');
}
void _dumpNeighborList(FILE *ofp, pList *list, char *header,
	char *currnode, int dir, NbrFunc *func)
{
	listIter iter;
	Neighbor *nbr;
	int id1, id2;
	int currid = func->getID(currnode), otherid;
	setListIter(&iter, list, 1);
	while (nbr = (Neighbor *) getListIter(&iter)) {
		id1 = func->getID(nbr->node1);
		id2 = func->getID(nbr->node2);
		if (currid == id1) {
			otherid = id2;
		} else if (currid == id2) {
			otherid = id1;
		} else {
			fprintf(stderr, "dump error\n");
			abort();
		}
		if (currid <= otherid) {
			fprintf(ofp, "%s %d %d", header, currid, otherid);
			fprintf(ofp, " %d %d %d %d\n", nbr->count,
				dir, dir * nbr->dir, nbr->status);
		}
	}
}
int _getMaxNeighbor(pList *list, char *orignode, char nodelCheck,
		char **ret_maxnode, Neighbor **ret_nbr)
{
	listIter iter;
	Neighbor *nbr, *maxnbr = NULL;
	char *node2, *maxnode = NULL;
	int i = 0;
	int maxcount = 0;

	setListIter(&iter, list, 1);
	while (nbr = (Neighbor *) getListIter(&iter)) {
		if (nbr->node1 == orignode) {
			node2 = nbr->node2;
		} else {
			node2 = nbr->node1;
		}
		/** valid for Node neighbor only **/
/*
		if (nbr->status) {
		if (nodelCheck && nbr->status) {
			continue;
		}
*/
		if (nodelCheck && ! isRoot(node2)) {
			continue;
		}
		if (maxcount < nbr->count) {
			maxcount = nbr->count;
			maxnode = node2;
			maxnbr = nbr;
		}
	}
	if (ret_maxnode) *ret_maxnode = maxnode;
	if (ret_nbr) *ret_nbr = maxnbr;
	
	return maxcount;
}
int _getLargeNeighbors(pList *nlist, double ovlpcut, double inclcut,
	int mincut, char *orignode, int nodelCheck, pList *retlist)
{
	listIter iter;
	Neighbor *nbr;
	char *node2;
	int flag;
	int cnt1, cnt2;
	setListIter(&iter, nlist, 1);

	while (nbr = (Neighbor *) getListIter(&iter)) {
		if (nbr->node1 == orignode) {
			node2 = nbr->node2;
		} else {
			node2 = nbr->node1;
		}
		/** valid for Node neighbor only **/
/*
		if (nbr->status) {
			continue;
		}
*/
		if (nodelCheck && ! isRoot(node2)) {
			continue;
		}
		flag = 0;

		cnt1 = ((Node *)nbr->node1)->cnt;
		cnt2 = ((Node *)nbr->node2)->cnt;

		if (ovlpcut) {
			if ( nbr->count >= rint(max(cnt1,cnt2) * ovlpcut) ) {
				flag = 1;
			}
		}
		if (inclcut) {
			if ( nbr->count >= rint(min(cnt1,cnt2) * inclcut) ) {
				flag = 1;
			}
		}
		if (mincut) {
			if ( nbr->count > min(cnt1,cnt2) - mincut ) {
				flag = 1;
			}
		}
		if (flag) {
			pushList(retlist, nbr);
		}
	}
}
_neighborDeleted(Neighbor *nbr, NbrFunc *nbrf)
{
#ifdef WITH_DEBUG
/*
printf("DelCheck\n");
nbrf->printObj(nbr->node1); printf("%d\n", nbrf->delcheck(nbr->node1));
nbrf->printObj(nbr->node2); printf("%d\n", nbrf->delcheck(nbr->node2));
*/
#endif
	if (nbrf->delcheck(nbr->node1) || nbrf->delcheck(nbr->node2)) {
		return 2;
	}
	if (nbr->status) {
		return 1;
	}
	return 0;
}

/**** ------------------------------------------ ****/
/***
	Node Neighbor
**/

/** node access methods **/
pList *getNodeLeft (Node *node){ return node->left; }
pList *getNodeRight (Node *node){ return node->right; }
int getNodeID(Node *node) {return node->id; }
NodeID *getNodeIDptr(Node *node) {return &(node->id); }
char getNodeDir(Node *node) {return node->dir; }
void setNodeLeft (Node *node, pList *plist){ node->left = plist; }
void setNodeRight (Node *node, pList *plist){ node->right = plist; }
void setNodeMark (Node *node){ node->flag |= NODE_TMPMARK; }
void unsetNodeMark (Node *node){ node->flag &= ~ NODE_TMPMARK; }
void setNodeMark2 (Node *node){ node->flag |= NODE_TMPMARK2; }
void unsetNodeMark2 (Node *node){ node->flag |= NODE_TMPMARK2; }
/*
int nodeMarked (Node *node){ return (node->flag & NODE_TMPMARK); }
*/
int nodeMarked2 (Node *node){ return (node->flag & NODE_TMPMARK2); }

NbrFunc nbrFuncNode = {
	"Node",
	getNodeLeft,
	getNodeRight,
	getNodeID,
	getNodeIDptr,
	setNodeLeft,
	setNodeRight,
	getNodeDir,
	setNodeMark,
	unsetNodeMark,
	nodeMarked,
	setNodeMark2,
	unsetNodeMark2,
	nodeMarked2,
	printNode,
	isNotRoot,
};

int addNeighbor(Node *leftnode, Node *rightnode, int dir1, int dir2)
{
	return _addNeighbor((char *) leftnode, (char *) rightnode, dir1, dir2,
		&nbrFuncNode);
}
int addNeighborCnt(Node *leftnode, Node *rightnode, int dir1,
		int dir2, Count count)
{
	return _addNeighborCnt((char *) leftnode, (char *) rightnode,
		dir1, dir2, (int) count, &nbrFuncNode);
}
/** for restoring dumped data */
int addNeighborRestore(Node *leftnode, Node *rightnode, int dir1, int dir2,
		Count count, char status)
{
	return _addNeighborRestore((char *) leftnode, (char *) rightnode,
		dir1, dir2, (int) count, status, &nbrFuncNode);
}
int isNeighbor(Node *node1, Node *node2)
{
	return _isNeighbor((char *) node1, (char *) node2, &nbrFuncNode);
}
int mergeNeighborList(Node *node1, Node *node2, Node *newnode, int dir)
{
	int ret;
#ifdef WITH_DEBUG
/*
if (newnode && newnode->id==16490 ||
	newnode && newnode->id==18826 ||
	newnode && newnode->id==22408 ||
	newnode && newnode->id==21386 ||
	newnode && newnode->id==22926){
printf("ID>>%d\n",newnode->id);
if (node1) {
	printf("1)%d,%d ",node1->id,dir);
	if (dir>0){
		if (node1->right) printNeighborList(node1->right,node1);
	} else {
		if (node1->left) printNeighborList(node1->left,node1);
	}
	printf("\n");
}
if (node2) {
	printf("2)%d,%d ",node2->id,dir);
	if (dir>0){
		if (node2->right) printNeighborList(node2->right,node2);
	} else {
		if (node2->left) printNeighborList(node2->left,node2);
	}
	printf("\n");
}
}
**/
#endif
	ret = _mergeNeighborList((char *) node1, (char *) node2,
		(char *) newnode, dir, &nbrFuncNode);
#ifdef WITH_DEBUG
/**
if (newnode && newnode->id==16490 ||
	newnode && newnode->id==18826 ||
	newnode && newnode->id==22408 ||
	newnode && newnode->id==21386 ||
	newnode && newnode->id==22926){
printf("MERGE: ");
if (dir>0){
printf("%d", newnode->right);
	if (newnode->right) printNeighborList(newnode->right,newnode);
} else {
printf("%d", newnode->left);
	if (newnode->left) printNeighborList(newnode->left,newnode);
}
printf("\n");
}
**/
#endif
	return ret;
}
int restoreNeighborList(Node *node0, Node *node1, int dir)
{
	return _restoreNeighborList((char *) node0, (char *) node1, &nbrFuncNode);
}
Node *getNeighborNode(Node *node, int dir, int *nbrdir)
{
	return (Node *) _getNeighbor((char *)node, dir, nbrdir, &nbrFuncNode);
}
void printNeighborList(pList *list, Node *orignode)
{
	_printNeighborList(list, (char *)orignode, &nbrFuncNode);
}
void dumpNeighborList(FILE *ofp, pList *list, char *header, Node *currnode,
		int dir)
{
	_dumpNeighborList(ofp, list, header, (char*)currnode,
				dir, &nbrFuncNode);
}
void addElemNeighborList(Node *node1, Node *node2, int dir1, int dir2,
		int cnt)
{
	_addElemNeighborList((char*)node1, (char*)node2, dir1, dir2,
					cnt, 0, &nbrFuncNode);
}
int neighborDeleted(Neighbor *nbr) {
	_neighborDeleted(nbr, &nbrFuncNode);
}
int getMaxNeighbor(pList *list, Node *orignode, char nodelCheck,
		Node **ret_maxnode, Neighbor **ret_nbr)
{
	return _getMaxNeighbor(list, (char *) orignode, nodelCheck,
		(char **) ret_maxnode, ret_nbr);
}
int getLargeNeighbors(pList *nlist, double ovlpcut, double inclcut,
		int mincut, Node *orignode, int nodelCheck, pList *retlist)
{
	return _getLargeNeighbors(nlist, ovlpcut, inclcut, mincut,
		(char *) orignode, nodelCheck, retlist);
}
#ifdef WITH_NEIGHBOR
int nodeNeighbor_df_search(Node *node, int dir,
	int (*func_before)(), char **args_before,
	int (*func_after)(), char **args_after,
	int (*func_nodecut)(), char **args_nodecut,
	int (*func_nbrcut)(), char **args_nbrcut)
{
	_df_search((char *) node, dir, &nbrFuncNode, func_before, args_before,
		func_after, args_after,
		func_nodecut, args_nodecut, func_nbrcut, args_nbrcut);
}
pList *nodeNeighbor_toposort(Node *node, int (*func_nodecut)(),
			int (*func_nbrcut)())
{
	return _toposort((char*)node, &nbrFuncNode, func_nodecut, func_nbrcut);
}
#endif
#ifdef WITH_NEIGHBOROUT
pList *nodeNeighbor_gapsearch(NodeSet *nodes, Node *node)
{
	int id, cnt, nextdir;
	pList *pl  = create_pList();
	Node *nn = node;
	Node *nnn;
	int nodecnt = node->cnt;
	int dir;
	pushList(pl, node);
	dir = 1;
	while ( 1 ) {
		gap_search((char*)node, dir, &nbrFuncNode, &id, &cnt, &nextdir);
		nnn = getNode(nodes, id);

/*
printf("#L>>>%d,%s:%d,%d,%d,%d\n", cnt,nnn->name,id,nodecnt,nnn->cnt, nextdir);
*/
/*
		if (cnt > 1 && cnt >= (double) nodecnt * Opt.nbrConnRatio
			|| cnt >= (double) node->cnt * Opt.nbrConnRatio) {
*/
		if (cnt > 1 && cnt >= (double) nodecnt * Opt.nbrConnRatio
			&& cnt >= (double) nnn->cnt * Opt.nbrConnRatio) {
			node = getNode(nodes, id);
			if (nbrFuncNode.Marked2(node)) {
/*
printf("---- %s\n", node->name);
*/
				break;
			}
/*
printf(">>%s\n", node->name);
*/
			pushList(pl, node);
			nbrFuncNode.setMark2(node);
		} else {
/*
printf("==== %s %d %d %d\n", node->name, cnt, node->cnt, nodecnt);
*/
			break;
		}
		dir = nextdir;
	}
	node = nn;
	dir = -1;
	while ( 1 ) {
		gap_search((char*)node, dir, &nbrFuncNode, &id, &cnt, &nextdir);
		nnn = getNode(nodes, id);
/*
printf("#R>>>%d,%s:%d,%d,%d,%d\n", cnt,nnn->name,id,nodecnt,nnn->cnt, nextdir);
*/
/*
		if (cnt > 1 && cnt >= nodecnt * Opt.nbrConnRatio
			|| cnt >= node->cnt * Opt.nbrConnRatio) {
*/
		if (cnt > 1 && cnt >= nodecnt * Opt.nbrConnRatio
			&& cnt >= nnn->cnt * Opt.nbrConnRatio) {
			node = getNode(nodes, id);
			if (nbrFuncNode.Marked2(node)) {
/*
printf("---- %s\n", node->name);
*/
				break;
			}
			unshiftList(pl, node);
			nbrFuncNode.setMark2(node);
		} else {
/*
printf("==== %s %d %d %d\n", node->name, cnt, node->cnt, nodecnt);
*/
			break;
		}
		dir = nextdir;
	}
	return pl;
}
#endif

#ifdef WITH_NEIGHBOR
/**** ------------------------------------------ ****/
/***
	Edge Neighbor
***/

/** edge access methods **/
pList *getEdgeLeft (Edge *edge){ return edge->left; }
pList *getEdgeRight (Edge *edge){ return edge->right; }
int getEdgeID(Edge *edge) {return edge->id; }
EdgeID* getEdgeIDptr(Edge *edge) {return (EdgeID*)&(edge->id); }
char getEdgeDir(Edge *edge) {return (char)1; }	/** don't use edge->dir **/
void setEdgeLeft (Edge *edge, pList *plist){ edge->left = plist; }
void setEdgeRight (Edge *edge, pList *plist){ edge->right = plist; }
void setEdgeMark (Edge *edge){ edge->flag |= EDGE_TMPMARK; }
void unsetEdgeMark (Edge *edge){ edge->flag &= ~ EDGE_TMPMARK; }
void setEdgeMark2 (Node *edge){ edge->flag |= EDGE_TMPMARK2; }
void unsetEdgeMark2 (Node *edge){ edge->flag |= EDGE_TMPMARK2; }
/*
int edgeMarked (Edge *edge){ return (edge->flag & EDGE_TMPMARK); }
*/
int edgeMarked2 (Edge *edge){ return (edge->flag & EDGE_TMPMARK2); }

NbrFunc nbrFuncEdge  = {
	"Edge",
	getEdgeLeft,
	getEdgeRight,
	getEdgeID,
	getEdgeIDptr,
	setEdgeLeft,
	setEdgeRight,
	getEdgeDir,
	setEdgeMark,
	unsetEdgeMark,
	edgeMarked,
	setEdgeMark2,
	unsetEdgeMark2,
	edgeMarked2,
	printEdge,
	deleted,
};

int addEdgeNeighbor(Edge *leftnode, Edge *rightnode)
{
	return _addNeighbor((char *) leftnode, (char *) rightnode, 1, 1,
		&nbrFuncEdge);
}
int addEdgeNeighborRestore(Edge *leftnode, Edge *rightnode,
		ConnCount count, char status)
{
	return _addNeighborCnt((char *) leftnode, (char *) rightnode,
		1, 1, (int) count, &nbrFuncEdge);
}
int addEdgeNeighborCnt(Edge *leftnode, Edge *rightnode, ConnCount count)
{
	return _addNeighborCnt((char *) leftnode, (char *) rightnode,
		1, 1, (int) count, &nbrFuncEdge);
}
int mergeEdgeNeighborList(Edge *node1, Edge *node2, Edge *newnode, int dir)
{
	return _mergeNeighborList((char *) node1, (char *) node2,
		(char *) newnode, dir, &nbrFuncEdge);
}
int restoreEdgeNeighborList(Edge *node0, Edge *node1, int dir)
{
	return _restoreNeighborList((char *) node0, (char *) node1, &nbrFuncEdge);
}
Edge *getNeighborEdge(Edge *node, int dir, int *nbrdir)
{
	return (Edge *) _getNeighbor((char *)node, dir, nbrdir, &nbrFuncEdge);
}
void printEdgeNeighborList(pList *list, Edge *orignode)
{
	_printNeighborList(list, (char *)orignode, &nbrFuncEdge);
}
void dumpEdgeNeighborList(FILE *ofp, pList *list, char *header, char *curredge,
		int dir)
{
	_dumpNeighborList(ofp, list, header, curredge, &nbrFuncEdge, dir);
}
int getMaxEdgeNeighbor(pList *list, Edge *orignode, char nodelCheck,
		Node **ret_maxnode, Neighbor **ret_nbr)
{
	return _getMaxNeighbor(list, (char *) orignode, nodelCheck,
		(char **) ret_maxnode, ret_nbr);
}
int edgeNeighbor_df_search(Edge *edge, int dir,
	int (*func_before)(), char **args_before,
	int (*func_after)(), char **args_after,
	int (*func_nodecut)(), char **args_nodecut,
	int (*func_nbrcut)(), char **args_nbrcut)
{
	_df_search((char *) edge, dir, &nbrFuncEdge, func_before, args_before,
		func_after, args_after,
		func_nodecut, args_nodecut, func_nbrcut, args_nbrcut);
}
pList *edgeNeighbor_toposort(Edge *edge, int (*func_nodecut)(),
			int (*func_nbrcut)())
{
	return _toposort(edge, &nbrFuncNode, func_nodecut, func_nbrcut);
}
#endif

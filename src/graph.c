/*
 * DomClust: Hierarchical Clustering for Orthologous Domain Classification
 * Copyright (c) 2000-2007, Ikuo Uchiyama
 * All rights reserved.
 */

#include "domclust.h"
#include "memalloc.h"
#include "plist.h"
#include "neighbor.h"

#include <assert.h>

pList *array2list();
pList *array2list_with_free();
int setEdge(EdgeSet *, Edge *, Node*, Node*, Region*, Region*, Dist);

typedef struct {
	int id;
	SeqPos from;
	Edge *edge;
} EdgePtr;

typedef int (*SortFunc) (const void *, const void *);

idcmp(EdgePtr *a, EdgePtr *b)
{
	if (a->id == b->id) {
		return a->from - b->from;
	} else {
		return a->id - b->id;
	}
}
idcmp_score(EdgePtr *a, EdgePtr *b)
{
	if (a->id == b->id) {
		return b->edge->score - a->edge->score;
	} else {
		return a->id - b->id;
	}
}

NodeSet *createNodes()
{
	NodeSet *nodes;

	if ((nodes = (NodeSet *) malloc(sizeof(NodeSet))) == 0) {
		allocError("NodeSet");
	}
	nodes->nodeobj = init_alloc_object(sizeof(Node), NODE_BLKSIZ);
	nodes->intobj = init_alloc_object(sizeof(int), INT_BLKSIZ);
	nodes->lenobj = init_alloc_object(sizeof(SeqPos), INT_BLKSIZ);
	nodes->idobj = init_alloc_object(sizeof(DomID), INT_BLKSIZ);
	nodes->nodenum = 0;
	return nodes;
}
setLeafNodeNum(NodeSet *nodes)
{
	nodes->leafnum = nodes->nodenum;
}
setTotalNodeNum(NodeSet *nodes)
{
	nodes->total_nodenum = nodes->nodenum;
}
TotalNodeNum(NodeSet *nodes)
{
	return nodes->total_nodenum;
}
getRootInit(NodeSetIter *iter, NodeSet *nodes, int dir, int withLeaf)
{
	iter->nodes = nodes;
	iter->dir = dir;
	if (dir > 0) {
		if (withLeaf) {
			iter->i = 0;
		} else {
			iter->i = nodes->leafnum;
		}
		iter->end = nodes->total_nodenum - 1;
	} else {
		iter->i = nodes->total_nodenum - 1;
		if (withLeaf) {
			iter->end = 0;
		} else {
			iter->end = nodes->leafnum;
		}
	}
}
Node *getRootNext_all(NodeSetIter *iter)
{
	return getRootNext_sub(iter, 1);
}
Node *getRootNext(NodeSetIter *iter)
{
	return getRootNext_sub(iter, 0);
}
#define IS_ROOT(x,allflag) ( (allflag) ? isAllRoot(x) : isRoot(x) )
Node *getRootNext_sub(NodeSetIter *iter, int allflag)
{
	Node *node = NULL;
	do {
		if ( (iter->dir > 0 && iter->i > iter->end) ||
			(iter->dir < 0 && iter->i < iter->end) ) {
			break;
		}
		if (! (node = getNode(iter->nodes, iter->i)) ) {
			break;
		}
		if (iter->dir > 0) {
			++iter->i;
		} else {
			--iter->i;
		}
	} while (! IS_ROOT(node,allflag));
	if (! IS_ROOT(node,allflag)) {
		/* No root node remains */
		return NULL;
	}
	return node;
}

addLinksForNewNode(Node **newnodes, Edge *edge)
{
	int i;
	int self_match = 0;
	pList *left1, *left2, *right1, *right2;

	if (edge->node1 == edge->node2) {
		self_match = 1;
	}

	/* add new node to the old nodes as a parent */
	edge->node1->parent = newnodes[0];
	if (! self_match) {
		edge->node2->parent = newnodes[0];
	}

	if (self_match) {
		if (newnodes[1]) {
			mergeNeighborList(edge->node1, NULL, newnodes[1], -1);
		} else {
			mergeNeighborList(edge->node1, NULL, newnodes[0], -1);
		}
		if (newnodes[4]) {
			mergeNeighborList(edge->node1, NULL, newnodes[4], 1);
		} else {
			mergeNeighborList(edge->node1, NULL, newnodes[0], 1);
		}
	} else {
		if (newnodes[1]) {
			mergeNeighborList(edge->node1, NULL, newnodes[1], -1);
		} else if (newnodes[2]) {
			mergeNeighborList(edge->node1, NULL, newnodes[0], -1);
		} else {
			mergeNeighborList(edge->node1, edge->node2, newnodes[0], -1);
		}
		if (newnodes[2]) {
			mergeNeighborList(NULL, edge->node2, newnodes[2], -1);
		} else if (newnodes[1]) {
			mergeNeighborList(NULL, edge->node2, newnodes[0], -1);
		}

		if (newnodes[3]) {
			mergeNeighborList(edge->node1, NULL, newnodes[3], 1);
		} else if (newnodes[4]) {
			mergeNeighborList(edge->node1, NULL, newnodes[0], 1);
		} else {
			mergeNeighborList(edge->node1, edge->node2, newnodes[0], 1);
		}
		if (newnodes[4]) {
			mergeNeighborList(NULL, edge->node2, newnodes[4], 1);
		} else if (newnodes[3]) {
			mergeNeighborList(NULL, edge->node2, newnodes[0], 1);
		}
	}

	if (newnodes[1]) {
		addNeighborCnt(newnodes[1], newnodes[0], 1, 1,edge->node1->cnt);
		edge->node1->parentL = newnodes[1];
	}
	if (newnodes[2]) {
		addNeighborCnt(newnodes[2], newnodes[0], 1, 1,edge->node2->cnt);
		edge->node2->parentL = newnodes[2];
	}
	if (newnodes[3]) {
		addNeighborCnt(newnodes[0], newnodes[3], 1, 1,edge->node1->cnt);
		if (self_match) {
			edge->node1->parentM = newnodes[3];
		} else {
			edge->node1->parentR = newnodes[3];
		}
	}
	if (newnodes[4]) {
		addNeighborCnt(newnodes[0], newnodes[4], 1, 1,edge->node2->cnt);
		edge->node2->parentR = newnodes[4];
	}
}
addLinksForNewEdge(Edge **newedges, Edge *edge1, Edge *edge2)
{
#ifdef WITH_NEIGHBOR
	if (newedges[1]) {
		mergeEdgeNeighborList(edge1, NULL, newedges[1], -1);
	} else if (newedges[2] && newedges[0]) {
		mergeEdgeNeighborList(edge1, NULL, newedges[0], -1);
	} else if (newedges[0]) {
		mergeEdgeNeighborList(edge1, edge2, newedges[0], -1);
	}
	if (newedges[2]) {
		mergeEdgeNeighborList(NULL, edge2, newedges[2], -1);
	} else if (newedges[1] && newedges[0]) {
		mergeEdgeNeighborList(NULL, edge2, newedges[0], -1);
	}

	if (newedges[3]) {
		mergeEdgeNeighborList(edge1, NULL, newedges[3], 1);
	} else if (newedges[4] && newedges[0]) {
		mergeEdgeNeighborList(edge1, NULL, newedges[0], 1);
	} else if (newedges[0]) {
		mergeEdgeNeighborList(edge1, edge2, newedges[0], 1);
	}
	if (newedges[4]) {
		mergeEdgeNeighborList(NULL, edge2, newedges[4], 1);
	} else if (newedges[3] && newedges[0]) {
		mergeEdgeNeighborList(NULL, edge2, newedges[0], 1);
	}

	if (newedges[1] && newedges[0]) {
		addEdgeNeighborCnt(newedges[1], newedges[0], edge1->connect);
	}
	if (newedges[2] && newedges[0]) {
		addEdgeNeighborCnt(newedges[2], newedges[0], edge2->connect);
	}
	if (newedges[3] && newedges[0]) {
		addEdgeNeighborCnt(newedges[0], newedges[3], edge1->connect);
	}
	if (newedges[4] && newedges[0]) {
		addEdgeNeighborCnt(newedges[0], newedges[4], edge2->connect);
	}
#endif
}

Node *addNode(NodeSet *nodes, char *name, int cnt, int len,
		Region *consreg,
		int totlen,
		Edge *child, NodeFlag flag, specFlag spflag, char truncFlag)
{
	Node *node = (Node *)memalloc(nodes->nodeobj);
	node->id = nodes->nodenum;
	setNode(nodes, node, name, cnt, len, child);
	copyReg(&(node->consreg), consreg);
	node->totlen = totlen;

	nodes->nodenum++;
	node->parent = node->parentL = node->parentR = NULL;
	node->parentM = NULL;
	node->left = node->right = NULL;
	node->domains = NULL;
	node->flag = flag;
	memcpy(node->spflag, spflag, sizeof(specFlag));
	node->dir = 1; /** ?? **/
	resetReg(&node->brk); resetReg(&node->brk2);
	resetReg(&node->newreg); resetReg(&node->newreg2);
	node->truncFlag = truncFlag;

	return node;
}

setNode(NodeSet *nodes, Node *node, char *name, int cnt, int len, Edge *child)
{
	if (name) {
		node->name = name;
	}
	node->child = child;
	node->cnt = (Count) cnt;
	node->len = (SeqPos) len;
}
Node *dupNode(NodeSet *nodes, Node *node)
{
	Node *newnode = addNode(nodes, node->name, node->cnt, node->len,
		&(node->consreg),
		node->totlen,
		node->child, node->flag, node->spflag, node->truncFlag);
	return newnode;
}
Node *dupNodeID(NodeSet *nodes, NodeID nodeid)
{
	return dupNode(nodes, getNode(nodes, nodeid));
}
setFlagNode(Node *node, NodeFlag flag)
{
	node->flag |= flag;
}
testFlagNode(Node *node, NodeFlag flag)
{
	return (node->flag & flag);
}
unsetFlagNode(Node *node, NodeFlag flag)
{
	node->flag &= ~flag;
}
clearFlagNode(Node *node)
{
	if (node->flag & NODE_CLEARED) {
		printf("ERROR: duplicated node: %s %d\n",node->name,node->id);
		abort();
	}
	setFlagNode(node, NODE_CLEARED);
}
printFlagNode(Node *node) {
	int i;
	int MAXIDX = 16;
	int flag = node->flag;
	int val;
	for (i = 0; i < MAXIDX; i++) {
		val = flag % 2;
		printf("%d ", val);
		flag >>= 1;
	}
	printf("\n");
}
makeRootNode(Node *node)
{
/*
printf("MakeRoot:");printNode(node);printf("\n");
*/
	unsetFlagNode(node, NODE_CLEARED);
}
/* make root node for outgroup mode */
makeRootNode_Out(Node *node)
{
	unsetFlagNode(node, NODE_DELETED);
	unsetFlagNode(node, NODE_OUTGRP);
	unsetFlagNode(node, NODE_CLEARED);
}
unclearFlagNode(Node *node)
{
	unsetFlagNode(node, NODE_CLEARED);
}

Node *getNode(NodeSet *nodes, int nodeid)
{
	Node *node = (Node *) get_objdata_idx(nodes->nodeobj, nodeid);
	return node;
}

deleteRootNode(Node *node)
{
	node->flag |= NODE_DELETED;
}
deleteRootNode_restoreDom(Node *node)
{
	deleteRootNode(node);
	resetNodeDomainAll(node);
}
nodeDeleted0(Node *node)
{
	return ( (node->flag & NODE_DELETED) && ! (node->flag & NODE_OUTGRP) );
}
nodeVisited_flag(NodeFlag flag)
{
	/* root or deleted or merged */
	return (! (flag & NODE_CLEARED));
}
nodeVisited(Node *node)
{
	return (nodeVisited_flag(node->flag));
}
nodeMerged(Node *node)
{
	return (node->flag & NODE_MERGED);
}
setMergeDir(Node *node,char dir)
{
	if (dir == 'R') {
		node->flag |= NODE_MERGED2;
	}
}
getMergeDir(Node *node) {
	if (node->flag & NODE_MERGED2) {
		/* right node */
		return 1;
	} else {
		return 0;
	}
}
isIntNode(Node *node)
{
	return (node->flag & NODE_INTERNAL);
}
isFlankNode1(Node *node)
{
	return ((node->flag & NODE_INT1) && ! (node->flag & NODE_INT2));
}
isFlankNode2(Node *node)
{
	return (! (node->flag & NODE_INT1) && (node->flag & NODE_INT2));
}
isFlankNode(Node *node)
{
	register NodeFlag flag = node->flag & NODE_INTERNAL;
	return (flag && flag != NODE_INTERNAL);
}
isLeaf(Node *node)
{
	return (! isIntNode(node));
}
isNotRoot(Node *node)
{
	return (! isRoot(node) );
}
isRoot_flag(NodeFlag flag)
{
	int isAllRoot_flag(NodeFlag);
/*
	return (! (flag & (NODE_CLEARED|NODE_DELETED|NODE_MERGED|NODE_OUTGRP)) );
*/
	return (! (flag & (NODE_CLEARED|NODE_DELETED|NODE_MERGED|NODE_OUTGRP)) ||
			(isAllRoot_flag(flag) &&  (flag & NODE_CENTER_ROOT)) );
}
isRoot(Node *node)
{
	if (node==NULL) return 0;
	return (isRoot_flag(node->flag));
/*
	return (isRoot0(node->flag) || 
		(isAllRoot(node) &&  (node->flag & NODE_CENTER_ROOT)) );
*/
}
/* without center root */
isRootRaw_flag(NodeFlag flag)
{
	return (! (flag & (NODE_CLEARED|NODE_DELETED|NODE_MERGED|NODE_OUTGRP)) );
}
isRootRaw(Node *node)
{
	if (node==NULL) return 0;
	return (isRootRaw_flag(node->flag));
}

isOutRoot(Node *node)
{
	return(isRoot(node) || isOuterRoot(node));
}
/* Root + InRoot */
isAllRoot_flag(NodeFlag flag)
{
	return (! (flag & (NODE_CLEARED|NODE_DELETED|NODE_MERGED)) );
}
isAllRoot(Node *node)
{
	if (node==NULL) return 0;
	return (isAllRoot_flag(node->flag));
}
/* InRoot: root of ingroup subcluster */
isInRoot(Node *node)
{
	if (node==NULL) return 0;
	return (isAllRoot(node) && ! isRoot(node));
}
isAllInRoot(Node *node)
{
	/* in-root under an out-root plus independent in-root (=out-root) */

	if (node==NULL) return 0;

	if (isAllRoot(node)) {
		if (node->child) {
			if ( (! node->child->node1 || ! nodeVisited(node->child->node1)) ||
				(! node->child->node2 || ! nodeVisited(node->child->node2)) ) {
				return 1;
			}
		} else {
			return 1;
		}
	}
	return 0;
}
isDeletedOuterRoot_flag(NodeFlag flag)
{
	return ((flag & NODE_OUTERROOT_DEL) != 0);
}
isDeletedOuterRoot(Node *node)
{
	return (isDeletedOuterRoot_flag(node->flag));
}
isOuterRoot_flag(NodeFlag flag)
{
	return ( ((flag & NODE_OUTERROOT) != 0) && (flag & NODE_OUTERROOT_DEL) == 0 );
}
isOuterRoot(Node *node)
{
	return ( isOuterRoot_flag(node->flag) );
}

setOutGrpNode(Node *node)
{
	node->flag |= NODE_OUTGRP;
}
unsetOutGrpNode(Node *node)
{
	node->flag &= ~NODE_OUTGRP;
}
isOutGrpNode(Node *node)
{
	return (node->flag & NODE_OUTGRP);
}

/* InGrpRoot -- currently not used */
setInGrpRoot(Node *node)
{
	node->flag |= NODE_INGRP;
}
isInGrpRoot(Node *node)
{
	return (node->flag & NODE_INGRP);
}
makeCenterRoot(Node *node)
{
	makeRootNode(node);
	node->flag |= NODE_CENTER_ROOT;
}
resetCenterRoot_flag(NodeFlag flag)
{
	flag &= ~NODE_CENTER_ROOT;
}
resetCenterRoot(Node *node)
{
	resetCenterRoot_flag(node->flag);
}
isCenterRoot_flag(NodeFlag flag)
{
	return (isRoot_flag(flag) && flag & NODE_CENTER_ROOT);
/*
	return (isAllRoot_flag(flag) && flag & NODE_CENTER_ROOT);
*/
/*
	return (flag & NODE_CENTER_ROOT);
*/
}
isCenterRoot(Node *node)
{
	return (isCenterRoot_flag(node->flag));
}
nodeDefinedParentL(Node *node)
{
	return (node->parentL && node->brk.from != INFPOS);
}
nodeDefinedParentC(Node *node)
{
	return (node->parent && node->brk.from != node->brk.to);
}
nodeDefinedParentR(Node *node)
{
	return (node->parentR && node->brk.to != SUPPOS);
}
nodeExistParentL(Node *node)
{
	return (node->parentL != NULL);
}
nodeExistParentC(Node *node)
{
	return (node->parent != NULL);
}
nodeExistParentR(Node *node)
{
	return (node->parentR != NULL);
}
/* outer Root */
setDeletedOuterRoot(Node *node)
{
	node->flag  |= NODE_OUTERROOT_DEL;
}
setOuterRoot(Node *node)
{
	node->flag |= NODE_OUTERROOT;
}
nodeOuterVisited(Node *node)
{
	return ( (node->flag & NODE_OUTERROOT_DEL) != 0 );
}

nodeMarked(Node *node)
{
	return ((node->flag & NODE_TMPMARK) != 0);
}
markNode(Node *node)
{
	node->flag |= NODE_TMPMARK;
}

clearAllMarks_sub(pList *nodelist, NodeFlag flag)
{
	listIter iter;
	Node *n;
	setListIter(&iter, nodelist, 1);
	while (n = (Node *) getListIter(&iter)) {
		unsetFlagNode(n, flag);
	}
}
clearAllMarks(pList *nodelist)
{
	clearAllMarks_sub(nodelist, (NodeFlag) NODE_TMPMARK);
}

getChilds(Node *node, Node **child1, Node **child2)
{
	*child1 = *child2 = NULL;
	if (!node->child || isLeaf(node)) {
	} else if (isFlankNode1(node)) {
		*child1 = node->child->node1;
	} else if (isFlankNode2(node)) {
		*child2 = node->child->node2;
	} else if (isIntNode(node)) {
		*child1 = node->child->node1;
		*child2 = node->child->node2;
	} else {
		printf(">>>????\n");
	}
}
specFlagP getSpFlag(Node *node) {
	return node->spflag;
}
nodeLen(Node *node)
{
	return regLen(&(node->consreg));
}
checkNodeLen(Node *node)
{
	int nodelen = nodeLen(node);
	if (nodelen >= Opt.minlen2) {
		return 2;
	} else if (nodelen >= Opt.minlen) {
		return 1;
	} else {
		return 0;
	}
}

EdgeSet *createEdges(int nodenum)
{
	register int i;
	int size;
	EdgeSet *edges;
	if ((edges = (EdgeSet *) malloc(sizeof(EdgeSet))) == 0) {
		allocError("EdgeSet");
	}
	edges->edgeobj = init_alloc_object_with_freelist(sizeof(Edge), EDGE_BLKSIZ);
	edges->idobj = init_alloc_object(sizeof(DomID), EDGE_BLKSIZ);
	edges->lenobj = init_alloc_object(sizeof(SeqPos), EDGE_BLKSIZ);
	edges->regobj = init_alloc_object(sizeof(StrRegion), EDGE_BLKSIZ);
	edges->edgenum = 0;
/*
	size = nodenum / 5000;
	size = (size > AVE_REL) ? size : AVE_REL;
*/
	edges->tmplist = create2DVarArray(nodenum, AVE_REL,
				sizeof(EdgePtr));
	return edges;
}
Edge *addEdgeWithScore(EdgeSet *edges, Node* node1, Node* node2, 
		Region *ali1, Region *ali2, Dist dist, Dist score,
		ConnCount connect, signed char dir)
{
	Edge *edge = addEdge(edges, node1, node2, ali1, ali2, dist, connect, dir);
	if (edge) {
		edge->score = score;
	}
	return edge;
}
Edge *addEdge(EdgeSet *edges, Node* node1, Node* node2, 
		Region *ali1, Region *ali2, Dist dist,
		ConnCount connect, signed char dir)
{
	Edge *edge;
	if (node1 == node2) {
		if (ali1->from == ali2->from && ali1->to == ali2->to) {
/*
			fprintf(stderr, "Warning: edges between identical sequences\n"); 
*/
			return 0;
		} else if (ali1->from > ali2->from) {
			Node *tmpn;
			Region *tmpr;
			tmpn = node1; node1 = node2; node2 = tmpn;
			tmpr = ali1; ali1 = ali2; ali2 = tmpr;
		}
	}

	edge = (Edge *) memalloc(edges->edgeobj);
	if (! edge){
		allocError("edge");
	}

/*
assert(ali1->to <= node1->len);
assert(ali2->to <= node2->len);
*/
	setEdge(edges, edge, node1, node2, ali1, ali2, dist);
	if (edges->tmplist) {
		/* initial construction. */ 
		addEdgePtr(edges, edge);
	} else {
		/* updating. node1 is the node newly created */ 
/*
		addEdgeIndex(edges, edge, node1->id, node2->id);
*/
	}
#ifdef LARGE
	edge->dir = dir;
#endif
#ifdef EXPERIMENTAL
	edge->connect = connect;
#endif
	edges->edgenum++;
#ifdef WITH_NEIGHBOR
	edge->id = edges->edgenum;
	edge->left = edge->right = NULL;
#endif
/*
	if (Opt.neighbor) {
		edge->left = create_pList(); edge->right = create_pList();
	} else {
		edge->left = edge->right = NULL;
	}
*/
	return edge;
}


setEdge(EdgeSet *edges, Edge *edge, Node* node1, Node* node2,
		Region *ali1, Region *ali2, Dist dist)
{
	register int i;
	Node *tmpnode;
	Region *tmpali;

	if (node1->id > node2->id) { /* swap */
		tmpnode = node1; node1 = node2; node2 = tmpnode;
		tmpali = ali1; ali1 = ali2; ali2 = tmpali;
	}

	edge->node1 = node1; edge->node2 = node2;
	edge->ali1 = (StrRegion *) memalloc(edges->regobj);
	edge->ali2 = (StrRegion *) memalloc(edges->regobj);
	Reg2StrReg(edge->ali1, ali1);
	Reg2StrReg(edge->ali2, ali2);

	edge->dist = dist;
	edge->flag = 0;
	edge->ordelem1 = edge->ordelem2 = NULL;
}
setEdgeFlag(Edge *edge, EdgeFlag flag)
{
	edge->flag |= flag;
}

delEdge(EdgeSet *edges, Edge *edge)
{

	if (Opt.DEBUG & DBG_cluster) {
		printf("Del: "); printEdge(edge);
	}

	if (isFree(edge->binelem)) {
		/* probably the case of deleting the best edge */
	} else {
		delBinData(edges->bin, (double) MEASURE(edge), edge->binelem);
	}
	edge->binelem = NULL;

/*
if (edge->ordelem1->datum != edge) {
	printf("ERROR!!!: "); printEdge(edge);
	printEdge((Edge*) edge->ordelem1->datum);
}
*/

	delelemList(&(edges->ordlist[edge->node1->id]), edge->ordelem1);
/*
execList(&(edges->ordlist[edge->node1->id]),printEdge2); putchar('\n');
*/

	delelemList(&(edges->ordlist[edge->node2->id]), edge->ordelem2);
	free_object_data(edges->edgeobj, edge);
}

addEdgePtr(EdgeSet *edges, Edge *edge)
{
	NodeID id1, id2;
	EdgePtr ep;
	id1 = edge->node1->id; id2 = edge->node2->id;

	ep.edge = edge;

	ep.id = id2;
	ep.from = edge->ali2->from;
	add2DArray(edges->tmplist, id1, &ep);

	ep.id = id1;
	ep.from = edge->ali1->from;
	add2DArray(edges->tmplist, id2, &ep);

/*	adding the edge to the bin is deferred until createIndex().
	addEdgeToBin(edges->bin, edge);
*/
}

addEdgeIndex(EdgeSet *edges, Edge *edge)
{
    /** node3_id < new_id **/
	NodeID node3 = edge->node1->id, new = edge->node2->id;
	addEdgeToBin(edges->bin, edge);

   /* Creating the sorted list of the edge index: */ 
   /* newnode should have the largest ID, so we can simply push the new edge */
	pushList(&(edges->ordlist[node3]), edge);
	edge->ordelem1 = getListElemIdx(&(edges->ordlist[node3]), -1);

   /* node3 should also have the largest ID among those ever pushed */
	if (node3 < new) {
		pushList(&(edges->ordlist[new]), edge);
		edge->ordelem2 = getListElemIdx(&(edges->ordlist[new]), -1);
	}
}

Better(double a, double b, int simflag)
{
	if (simflag) {
		return (a > b);
	} else {
		return (a < b);
	}
}
double DistRatio(double a, double b, int simflag)
{
	if (simflag) {
		return (b / a);
	} else {
		return (a / b);
	}
}

Edge *getBestEdge(EdgeSet *edges)
{
	return (Edge *) getBestData(edges->bin, Opt.sim);
}
#ifdef EXPERIMENTAL
Dist getEdgeScoreSum(Edge *e)
{
/*
	return (e->score - Opt.missscore) * e->node1->cnt * e->node2->cnt
			+ Opt.missscore * e->connect;
*/
	return (e->score - Opt.missscore) * e->connect
		+ Opt.missscore * e->node1->cnt * e->node2->cnt;
}
#endif

/* Return each edge incident to a given node by successive calls
	with the same arguments. See slink0() for an example. */

Edge *getEdgeByNode(EdgeSet *edges, NodeID node, listIter **iter)
{
	Edge *e;
	if (*iter == NULL) {
		*iter = createListIter(&(edges->ordlist[node]), 1);
	}
	e = getEdgeByNode0(edges, *iter);
	if (e) {
		return e;
	} else {
		*iter = NULL;
		return NULL;
	}
}
Edge *getEdgeByNode0(EdgeSet *edges, listIter *iter)
{
	Edge *edge;

	do {
		edge = (Edge *) getListIter(iter);
	} while (edge && deleted(edge));

	return edge;
}

printNodeName(Node *node)
{
	if (isName(node->name)) {
		printf("%s",node->name);
	} else {
		printf("NODE");
	}
}
void printNode(Node *node)
{
	printf("%s(%d)",node->name,node->id);
}
void printNode2(Node *node)
{
	printf("%s(%d)[%d,%d,%d]",node->name,node->id,node->minlen,node->len,node->meanlen);
}

printEdgeERR(Edge *edge)
{
	fprintf(stderr, "%s(%d) %s(%d) %lf\n",
		edge->node1->name, edge->node1->id,
		edge->node2->name, edge->node2->id, (double)MEASURE(edge));
}
void printEdge(Edge *edge)
{
	printf("%s[%d](%d,[%d:%d],%d) %s[%d](%d,[%d,%d],%d) {%d,%d}{%d,%d} %lf",
		edge->node1->name, edge->node1->id, edge->node1->len,
		edge->node1->consreg.from,edge->node1->consreg.to,
		edge->node1->cnt,
		edge->node2->name, edge->node2->id, edge->node2->len,
		edge->node2->consreg.from,edge->node2->consreg.to,
		edge->node2->cnt,
		edge->ali1->from,edge->ali1->to,
		edge->ali2->from,edge->ali2->to,(double)MEASURE(edge));
#ifdef WITH_NEIGHBOR
	printf(" <%d>", ,edge->id);
#endif
	printf("\n");
}
void printEdge2(Edge *edge)
{
	printf("(%d,%d)", edge->node1->id, edge->node2->id);
}

createOrdList(SimGraph *SimG)
{
	int i;
	Edge *edge;

	sort2DArrayVerbose(SimG->edges->tmplist, idcmp_score);
	SimG->edges->ordlist = array2list(SimG->edges->tmplist);
	free2DArray(SimG->edges->tmplist);
	SimG->edges->tmplist = NULL;
}
sort2DArrayVerbose(varArray *vararray, int (*func)(const void *, const void *))
{
	int i;
	varArray **var = (varArray **) vararray->array;
#pragma omp parallel for
	for (i = 0; i < vararray->num; i++) {
		qsort(var[i]->array, var[i]->num, var[i]->datasize, func);
		if (Opt.VERBOSE) {
			if ((i+1) % Opt.verbose_step2 == 0) {
				fprintf(stderr, "sorting %d\n", i+1);
			}
		}
	}
}

createIndex(SimGraph *SimG)
{
	unsigned long long i;
	EdgeSet *edges = SimG->edges;
	Edge *edge;
	Dist getMaxDist(), getMinDist();
	ObjDataIter objIter;

	if (Opt.sim) {
		Dist maxscore = getMaxDist(edges);
		edges->bin = createBin(Opt.missscore, maxscore + 2,
			Opt.distscale, Opt.logdistscale);
	} else {
		Dist mindist = getMinDist(edges);
		edges->bin = createBin(mindist, Opt.missdist + 2,
			Opt.distscale, Opt.logdistscale);
	}
/*
	for (i = 0; i < edges->edgenum; i++) {
		edge = (Edge *) get_objdata_idx(edges->edgeobj, i);
		addEdgeToBin(edges->bin, edge);
		if (Opt.VERBOSE) {
			if ((i+1) % Opt.verbose_step1 == 0) {
				fprintf(stderr, "indexing %d\n", i+1);
			}
		}
	}
*/
	initObjIter(&objIter, edges->edgeobj);
	i = 0;
	while (edge = (Edge*) nextObjIter(&objIter)) {
		if (edge->flag & EDGE_DELETED) {
			/** initial deletion */
			continue;
		}
		addEdgeToBin(edges->bin, edge);
		if (Opt.VERBOSE) {
			if ((++i) % Opt.verbose_step1 == 0) {
				fprintf(stderr, "indexing %lld\n", i);
			}
		}
	}
}
Dist getMaxDist(EdgeSet *edges)
{
	Dist maxdist = -999999, s;	/* actually score rather than dist */
	Edge *edge;
	ObjDataIter objIter;
	initObjIter(&objIter, edges->edgeobj);
	while (edge = (Edge*) nextObjIter(&objIter)) {
		if (edge->flag & EDGE_DELETED) {
			/** initial deletion */
			continue;
		}
		if (maxdist < (s=MEASURE(edge))) {
			maxdist = s;
		}
	}
	return maxdist;
}
Dist getMinDist(EdgeSet *edges)
{
	Dist mindist = 999999, s;	/* actually score rather than dist */
	Edge *edge;
	ObjDataIter objIter;
	initObjIter(&objIter, edges->edgeobj);
	while (edge = (Edge*) nextObjIter(&objIter)) {
		if (edge->flag & EDGE_DELETED) {
			/** initial deletion */
			continue;
		}
		if (mindist > (s=MEASURE(edge))) {
			mindist = s;
		}
	}
	return mindist;
}

/*** routine for single linkage clustering ***/

static EdgeSet *Edges;
static NodeSet *Nodes;
static Cluster *cluster;

clidcmp(Cluster *a, Cluster *b)
{
	if (a->clid == b->clid) {
		return (a->id - b->id);
	} else {
		return (a->clid - b->clid);
	}
}
Cluster *slink(EdgeSet *edges, NodeSet *nodes)
{
	Edge *edge;
	Node *node;
	NodeID i, id;
	int clnum = 1;
	char *name;
	
	Edges = edges; Nodes = nodes;

	if ((cluster = (Cluster *) malloc(Nodes->nodenum * sizeof(Cluster)))
			== NULL) {
		allocError("clnums for scluster");
	}
	for (id = 0; id < Nodes->nodenum; id++) {
		cluster[id].id = id;
		cluster[id].clid = 0;
	}
	for (id = 0; id < Nodes->nodenum; id++) {
		if (cluster[id].clid) continue;
		cluster[id].clid = clnum++;
		slink0(id, clnum);
	}

	qsort(cluster, Nodes->nodenum, sizeof(Cluster), (SortFunc) clidcmp);
/*
	for (i = 0; i < Nodes->nodenum; i++) {
		id = nodenums[i];
		node = getNode(Nodes,id);
		printf("%s %d %d\n", node->name, id, clnums[id]);
	}
*/
	return cluster;
}
slink0(NodeID node1, int clnum)
{
	Edge *edge;
	listIter *iter = NULL;
	NodeID node2;
	while (edge = getEdgeByNode(Edges, node1, &iter)) {
#ifdef WITH_DEBUG
/*** DEBUG
		printf(">%d,%d,%d\n",node1,pair->id1,pair->id2);
		printf(">%d,%d, %s,%s\n",
			pair->edge->node1->id, pair->edge->node2->id,
			pair->edge->node1->name, pair->edge->node2->name);
***/
#endif
		node2 = otherNodeID(edge, node1);
		if (! cluster[node2].clid) {
			cluster[node2].clid = clnum;
			slink0(node2, clnum);
		}
	}
}

printCluster(Cluster *cluster, NodeSet *nodes)
{
	int i;
	Node *node;
	for (i = 0; i < nodes->nodenum; i++) {
		node = getNode(nodes,i);
		printf("%s %d %d\n", node->name, cluster[i].id, cluster[i].clid);
	}
}
addEdgeToBin(Bin *bin, Edge *edge)
{
	edge->binelem =
		addBinRetElem(bin, (double) MEASURE(edge), edge, NULL);
}

deleted(Edge *edge)
{
	if (! edge || edge == (Edge *) -1 || edge->flag & EDGE_CLEARED || edge->flag & EDGE_DELETED) {
		return 1;
	} else {
		return 0;
	}
}
edgeMarked(Edge *edge)
{
	return (edge->flag & EDGE_TMPMARK);
}
edgeSelectFlag(Edge *edge)
{
	edge->flag |= EDGE_SELECTED;
}
edgeSelected(Edge *edge)
{
	return (edge->flag & EDGE_SELECTED);
}
edgeSkipFlag(Edge *edge)
{
	edge->flag |= EDGE_SKIP;
}
edgeSkip(Edge *edge)
{
	return (edge->flag & EDGE_SKIP);
}
#define ARRAYSIZE_FACTOR 5
pList *array2list(varArray *ar)
{
	return(array2list_with_free(ar, NULL));
}
pList *array2list_with_free(varArray *ar, Alloc_Object *edgeobj)
{
	int i, j;
	varArray **ar2;
	EdgePtr *p, *prev_p;
	pList *plist;
	if ( (plist = create_pList_array(
			(int)((double)arraySize(ar) * ARRAYSIZE_FACTOR))) == NULL) {
		allocError("plist");
	}
	for (i = 0; i < arraySize(ar); i++) {
		ar2 = (varArray **) getArrayItem(ar, i);
		prev_p = NULL;
		for (j = 0; p = (EdgePtr *) getArrayItem(*ar2, j); j++) {
			if (p->edge->flag & EDGE_DELETED) {
				continue;
			}
			if (prev_p && prev_p->id == p->id) {
				/* skip duplicated edge (edges should be sorted by id and score */
				p->edge->flag |= EDGE_DELETED; /* initial deletion */
				if (edgeobj) {
					free_object_data(edgeobj, p->edge);
				}
				continue;
/*
				fprintf(stderr, "ERROR: Duplicated edge found:\n");
				printEdgeERR(prev_p->edge);
				printEdgeERR(p->edge);
				exit(1);
*/
			}
			prev_p = p;

			pushList(&(plist[i]), p->edge);
			if (! p->edge->ordelem1) {
				p->edge->ordelem1=getListElemIdx(&plist[i], -1);
			} else if (! p->edge->ordelem2) {
				p->edge->ordelem2=getListElemIdx(&plist[i], -1);
			} else {
				fprintf(stderr, "ERROR: ordelem\n");
			}
		}
	}
	return plist;
}

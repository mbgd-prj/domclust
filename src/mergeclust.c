/*
 * DomClust: Hierarchical Clustering for Orthologous Domain Classification
 * Copyright (c) 2000-2007, Ikuo Uchiyama
 * All rights reserved.
 */

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "domclust.h"
#include "util.h"
#include "neighbor.h"
#include "hash.h"
#include "plist.h"
#include "djset.h"

#define MAXSUBGRP 20000

int concatNeighbor(Node *nodeS, Node *nodeL, char dir);

cmpr_node_by_id(Domain **dom1, Domain **dom2)
{
	return (*dom1)->leaf->id - (*dom2)->leaf->id;
}
cmpr_node_by_id_pos(Domain **dom1, Domain **dom2)
{
	if ((*dom1)->leaf->id == (*dom2)->leaf->id) {
		return (*dom1)->from - (*dom2)->from;
	} else {
		return (*dom1)->leaf->id - (*dom2)->leaf->id;
	}
}
/*
cmpr_cinfo_by_homclust(ClusterInfo **cinfo1, ClusterInfo **cinfo2)
*/
cmpr_cinfo_by_homclust(const void *_cinfo1, const void *_cinfo2)
{
	ClusterInfo *cinfo1 = *(ClusterInfo**) _cinfo1;
	ClusterInfo *cinfo2 = *(ClusterInfo**) _cinfo2;

	if (cinfo1->homclust == cinfo2->homclust) {
		return cinfo1->outerclust - cinfo2->outerclust;
	} else {
		return cinfo1->homclust - cinfo2->homclust;
	}
}

Hash *create_cinfoHash(pList *clustRoots) {
        listIter iter;
        ClusterInfo *cInfo;
        HENTRY hent;
        Hash *clustInfoHash;
        int cinfonum = numelemList(clustRoots);
        clustInfoHash = Hcreate(cinfonum * 50);
        
        setListIter(&iter, clustRoots, 1);
        while (cInfo = (ClusterInfo*) getListIter(&iter)) {
                hent.key = (char *) (cInfo->root->node->id);
                hent.datum = (char *) cInfo;
                HIsearch(clustInfoHash, &hent, ENTER);
        }       
        return clustInfoHash;
}               

cinfoHsearch(Hash *cinfoHash, int nodeid, ClusterInfo **ci)
{
	HENTRY hent;
	hent.key = (char *) nodeid;
	if (HIsearch(cinfoHash, &hent, FIND)){
		*ci = (ClusterInfo *) hent.datum;
		return 1;
	} else {
		return 0;
	}
}
createClusterInfo_All(pList *treenodes, ClusterInfoData *cinfoData)
{
	ClusterInfo *cinfo, *ci;
	TreeNode *treenode;
	Node *node;
	listIter iter;
	int clustnum = numelemList(treenodes);
	int cinfo_size = clustnum * 1.2;
	int i;
	ClusterInfo *addClusterInfo(ClusterInfoData *, TreeNode *);
/*
	HENTRY hent;
*/
	int cinfoHash_flag = 1;

	if ( (cinfo = (ClusterInfo *)
		malloc(cinfo_size * sizeof(ClusterInfo)) ) == NULL) {
		fprintf(stderr, "Can't allocate memory\n");
		exit(1);
	}
	cinfoData->cinfo = cinfo;
	cinfoData->clustnum = 0;
	cinfoData->cinfo_link = NULL;
	setListIter(&iter, treenodes, 1);
	if (cinfoHash_flag) {
		cinfoData->cinfoHash = Hcreate(cinfo_size * 50);
	}

	i = 0;
	while (treenode = (TreeNode *) getListIter(&iter)) {
/*
		if (! isRoot(treenode->node)) {
*/
		if (! isOutTreeRoot(treenode)) {
			Node *n = treenode->node;
			fprintf(stderr, ">> invalid root: %s,%d,%d\n",n->name,n->id, treenode->flag);
			printf(">> invalid root: %s,%d,%d\n",n->name,n->id, treenode->flag);

			continue;
		}
/*
		cinfo[i].root = cinfo[i].origroot = treenode;
		cinfo[i].clusters = create_pList();
		pushList(cinfo[i].clusters, &cinfo[i]);
		cinfo[i].members = create_pList();
		cinfo[i].homclust = -1;
		cinfo[i].outgroup = cinfo[i].outgroup2 = NULL;
		cinfo[i].outer_root = NULL;
*/
		addClusterInfo(cinfoData, treenode);

		i++;
	}
	clustnum = i;
/*
	if (cinfoHash_flag) {
		cinfoData->cinfoHash = Hcreate(cinfo_size * 50);
		for (i = 0; i < clustnum; i++) {
			hent.key = (char *) cinfo[i].root->node->id;
			hent.datum = (char *) &(cinfo[i]);
			HIsearch(cinfoData->cinfoHash, &hent, ENTER);
			addCinfoHash(cinfoData, &(cinfo[i]));
		}
	}
*/
#ifdef WITH_DEBUG
/** for debug
	for (i = 0; i < cinfo_size; i++) {
		cinfoHsearch(*cinfoHash, cinfo[i].root->node->id, &ci);
		printTreeNode(cinfo[i].root);
		printTreeNode(ci->root);
		printf("--\n");
		printf("%d,%s; %d,%s\n", &cinfo[i], cinfo[i].root->node->name,
			ci, ci->root->node->name);
	}
**/
#endif
	cinfoData->clustnum = clustnum;
	if ( (cinfoData->cinfo_link = (ClusterInfo **)
			malloc(cinfoData->clustnum * sizeof(ClusterInfo *)) ) == NULL) {
		fprintf(stderr, "Can't allocate memory\n");
		exit(1);
	}
	for (i = 0; i < cinfoData->clustnum; i++) {
		cinfoData->cinfo_link[i] = &(cinfoData->cinfo[i]);
	}
	return clustnum;
}
ClusterInfo *addClusterInfo(ClusterInfoData *cinfoData, TreeNode *treenode)
{
	ClusterInfo *cinfo = &(cinfoData->cinfo[ cinfoData->clustnum ]);
	cinfoData->clustnum++;
	setClusterInfo(cinfo, treenode);
	if (cinfoData->cinfoHash) {
		addCinfoHash(cinfoData, cinfo);
	}
	return cinfo;
}
setClusterInfo(ClusterInfo *cinfo, TreeNode *treenode)
{
	cinfo->root = cinfo->origroot = treenode;

	cinfo->clusters = create_pList();
/*
	pushList(cinfo->clusters, treenode);
*/
	pushList(cinfo->clusters, cinfo);
	cinfo->members = create_pList();
	cinfo->homclust = cinfo->outerclust = -1;
	cinfo->outgroup = cinfo->outgroup2 = NULL;
	cinfo->outer_root = NULL;
/*
	cinfo->clusters_outer = NULL;
*/
	collectCluster(treenode, cinfo->members);
/*
	copySPFlag(treenode->node->spflag, cinfo->spflag);
*/
	copySPFlag(treenode->spflag, cinfo->spflag);
}
addCinfoHash(ClusterInfoData *cinfoData, ClusterInfo *cinfo) {
	HENTRY hent;
	hent.key = (char *) cinfo->root->node->id;
	hent.datum = (char *) cinfo;
	HIsearch(cinfoData->cinfoHash, &hent, ENTER);
}
cinfoMeanLen(ClusterInfo *cinfo)
{
	listIter iter;
	Domain *dom;
	int len = 0, cnt = 0;

	setListIter(&iter, cinfo->members, 1);
	while (dom = (Domain *) getListIter(&iter)) {
		len += (dom->to - dom->from + 1);
		cnt++;
	}
	if (len == 0) return 0;
	return (len / cnt);
}
printCinfo(ClusterInfo *cinfo)
{
	listIter iter;
	Domain *dom;
	printf("ClustID:%d (%d/%d)\n",cinfo->clustid, cinfo->homclust, cinfo->outerclust);
	printf("Root:"); printTreeNode(cinfo->root); printf("\n");
	if (cinfo->outer_root) {
		printf("OuterRoot:"); printTreeNode(cinfo->outer_root); printf("\n");
	}
	setListIter(&iter, cinfo->members, 1);
	while (dom = (Domain *) getListIter(&iter)) {
		printf("%s %d %d %d\n",dom->leaf->name,dom->num,dom->from,dom->to);
	}
	printf("--\n");

if (cinfo->outgroup && cinfo->clusters) {
ClusterInfo *cinfo0;
	setListIter(&iter, cinfo->clusters, 1);
	while (cinfo0 = (ClusterInfo *) getListIter(&iter)) {
		printf("ROOT>");printTreeNode(cinfo0->root);printf("\n");
	}
}
	printf("//\n");
}
/* renum domains */
renumCinfo(ClusterInfo *cinfo, int clustnum)
{
	int i;
	pList *alldom = create_pList();
	Domain *dom, *prevdom;
	int domn;
	for (i = 0; i < clustnum; i++) {
		pushListAll(alldom, cinfo[i].members);
/*
printf("cluster=%d:%d; ",i, cinfo[i].clustid);
printTreeNode(cinfo[i].root);
printList(cinfo[i].members,printDomain);
printf("\n");
*/
	}
	sortList(alldom, cmpr_node_by_id_pos);
	prevdom = NULL;
	while (dom = (Domain *) shiftList(alldom)) {
		if (! prevdom || dom->leaf->id != prevdom->leaf->id) {
			domn = 1;
		} else {
			domn++;
		}
		dom->num= domn;

		if (Opt.DEBUG & DBG_domain) {
			if (Opt.DEBUG_ent == NULL || strcmp(dom->leaf->name, Opt.DEBUG_ent) == 0) {
				printf(">renum_domnum>");printDomain(dom);printf(":%p\n",dom);
			}
		}

		prevdom = dom;
	}
	free_pList(alldom);
}
/** merge cluster1 into cluster2: cinfo1=>null, cinfo2=>merged_list **/
mergeClusters(ClusterInfo *cinfo1, ClusterInfo *cinfo2, NodeSet *nodes)
{
	listIter iter1, iter2;
	pList *merged_list;
	Domain *d1, *d2;
	int mergecnt = 0, cnt1 = 0, cnt2 = 0;
	int mergeflag = 0;

#ifdef WITH_DEBUG
/***
	if (! checkMergeClusters(cinfo1, cinfo2)) {
		return 0;
	}
***/
#endif

	sortList(cinfo1->members, cmpr_node_by_id);
	sortList(cinfo2->members, cmpr_node_by_id);

	/* Check merge condition */
	mergeClustList_count(cinfo1->members, cinfo2->members,
			&cnt1, &cnt2, &mergecnt);
/*
printf(">%s,%s,%d,%d,%d\n",cinfo1->root->node->name,cinfo2->root->node->name,
cnt1,cnt2,mergecnt);
*/

	if (Opt.adjOvlpRatio &&
			mergecnt >= rint(max(cnt1, cnt2) * Opt.adjOvlpRatio)) {
		mergeflag = 1;
	} else if (Opt.adjInclRatio || Opt.mincutcnt) {
		if ( cnt1 <= cnt2 && (
				(mergecnt >= rint(cnt1 * Opt.adjInclRatio)) ||
				(cnt1 - mergecnt < Opt.mincutcnt)) ) {
			if (cinfoMeanLen(cinfo1) <= Opt.minlen2) {
				mergeflag = 1;
			}
		} else if ( cnt2 <= cnt1 && (
				(mergecnt >= rint(cnt2 * Opt.adjInclRatio)) ||
				(cnt2 - mergecnt < Opt.mincutcnt)) ) {
			if (cinfoMeanLen(cinfo2) <= Opt.minlen2) {
				mergeflag = 1;
			}
		}
	}
	if (! mergeflag) {
		/* NG: do not update members */
		return 0;
	}

	/* OK: update members */
	merged_list = create_pList();
	mergeClustList(cinfo1->members, cinfo2->members, merged_list);
		
	if (mergeflag) {
		Node *root1, *root2;
		Node *newnode;
		TreeNode *newtreenode;
		root1 = cinfo1->root->node;
		root2 = cinfo2->root->node;

		newtreenode = dupTreeNode(nodes, cinfo2->root);
		newnode = newtreenode->node;

		/* OK: update members */
		freeList(cinfo1->members);
		freeList(cinfo2->members);
		cinfo2->members = merged_list;

/*
printList(cinfo1->clusters, printTreeNode);
printList(cinfo2->clusters, printTreeNode);
*/
		concatList(cinfo2->clusters, cinfo1->clusters);
		free(cinfo1->clusters);
		cinfo1->clusters = NULL;
		cinfo1->root = newtreenode;
/*
printf("merge>>>>\n");
printTreeNode(cinfo2->root);printf(" ");printTreeNode( cinfo2->origroot);printf("\n");
printTreeNode(cinfo1->root);printf(" ");printTreeNode( cinfo1->origroot);printf("\n");
*/
	
/*
		if (cinfo2->clusters_outer) {
			if (cinfo1->clusters_outer) {
				concatList(cinfo2->clusters_outer, cinfo1->clusters_outer);
				free(cinfo1->clusters_outer);
				cinfo1->clusters_outer = NULL;
			}
			cinfo1->outer_root = cinfo2->outer_root;
		}
*/

		mergeNeighborList(root1, root2, newnode, 1);
		mergeNeighborList(root1, root2, newnode, -1);

		spFlagOR(cinfo1->spflag, cinfo2->spflag, cinfo2->spflag);
		setTreeNode(cinfo2->root, newnode);

		setListIter(&iter1, merged_list, 1);
		while (d1 = (Domain *) getListIter(&iter1)) {
			d1->root = newnode;
		}

		return 1;
	}
}
reset_mergedDomains_all(ClusterInfo *cinfo)
{
	listIter iter;
	setListIter(&iter, cinfo->clusters, 1);
	ClusterInfo *cinfo0;
	TreeNode *tnode;
/*
	while (tnode = (TreeNode *) getListIter(&iter)) {
		reset_mergedDomains(tnode);
	}
*/
	while (cinfo0 = (ClusterInfo *) getListIter(&iter)) {
		reset_mergedDomains(cinfo0->origroot);
	}
}
reset_mergedDomains(TreeNode *tnode)
{
	int i = 0;
	static int maxiter = 10000;
int tmpflag = 0;
	if (tnode == NULL) {
	} else if (ClustTree_isLeaf(tnode)) {
/*
if (tnode->dom->num<0){
printf("RESET>");printDomain(tnode->dom);printf("\n");
tmpflag = 1;
}
*/
		/* renum merged domain which has domnum = -1 */
		while (tnode->dom->num == -1) {
			if (tnode->dom == (Domain*)tnode->dom->root || i++ > maxiter) {
				/* to avoid an infinite loop */
				printf("Too many loops??\n");
				exit(1);
			}
			tnode->dom = (Domain *) tnode->dom->root;
		}
/*
if (tmpflag) {
printDomain(tnode->dom);printf("\n");
}
*/
	} else {
		reset_mergedDomains(tnode->child1);
		reset_mergedDomains(tnode->child2);
	}
}
#ifdef WITH_DEBUG
/***
checkMergeClusters(ClusterInfo *cinfo1, ClusterInfo *cinfo2)
{
	int l1, l2;
	if ((l1=cinfoMeanLen(cinfo1)) > Opt.minlen2 &&
		(l2=cinfoMeanLen(cinfo2)) > Opt.minlen2) {
		return 0;
	}
	return 1;
}
***/
#endif

mergeClustList_count(pList *clist1, pList *clist2,
			int *cnt1, int *cnt2, int *mergecnt)
{
	_mergeClustList(clist1, clist2, NULL, cnt1, cnt2, mergecnt);
}
mergeClustList(pList *clist1, pList *clist2, pList *merged_list)
{
	_mergeClustList(clist1, clist2, merged_list, NULL, NULL, NULL);
}

_mergeClustList(pList *clist1, pList *clist2, pList *merged_list,
			int *cnt1, int *cnt2, int *mergecnt)
{
	listIter iter1, iter2;
	Domain *d1, *d2;
	int ign1, ign2;

	setListIter(&iter1, clist1, 1);
	d1 = (Domain *) getListIter(&iter1);
	setListIter(&iter2, clist2, 1);
	d2 = (Domain *) getListIter(&iter2);

	while (d1 || d2) {
		if (! d2 || (d1 && cmpr_node_by_id(&d1,&d2) < 0)) {
			if (merged_list) pushList(merged_list, d1);
			if (cnt1 && !(metagenomeOnlyCluster(d1->leaf->spflag))){
				++ *cnt1;
			}
			d1 = (Domain *) getListIter(&iter1);
		} else if (! d1 || (d2 && cmpr_node_by_id(&d2,&d1) < 0)) {
			if (merged_list) pushList(merged_list, d2);
			if (cnt2 && !(metagenomeOnlyCluster(d2->leaf->spflag))){
				++ *cnt2;
			}
			d2 = (Domain *) getListIter(&iter2);
		} else {
			ign1 = (metagenomeOnlyCluster(d1->leaf->spflag)) ? 1: 0;
			ign2 = (metagenomeOnlyCluster(d2->leaf->spflag)) ? 1: 0;

			if (mergecnt && !ign1 && !ign2) {
				++ *mergecnt;
			}
			if (cnt1 && ! ign1) {
				++ *cnt1;
			}
			if (cnt2 && ! ign2) {
				++ *cnt2;
			}
			/** concatenate domains **/
			if (merged_list) {
				if (d1->to <= d2->from) {
					d1->to = d2->to;
				} else if (d1->from >= d2->to) {
					d1->from = d2->from;
				} else {
				}
				pushList(merged_list, d1);

				deleteDomain(d2);
				/* save retained domain
					(for reset_mergedDomains) */
/*
printf(">>>DDDOM:");
printDomain(d1);
printNode(d1->root);
printf("; ");
printDomain(d2);
printNode(d2->root);
printf("\n");
*/
				d2->num = -1;
				d2->root = (Node *) d1;
			}
			d1 = (Domain *) getListIter(&iter1);
			d2 = (Domain *) getListIter(&iter2);
		}
	}
}

int checkOverlapClustersAll(pList *treenodes, ClusterInfoData *cinfoData, NodeSet *nodeset)
{
	int mergecnt;
	ClusterInfo *clust;
	Node *node;
	listIter iter;
	ClusterInfo *cinfo;
	int clustnum;
	int i, j;
	Hash *cinfoHash = cinfoData->cinfoHash;

/*
	Hash *tmp_cinfoHash = NULL;
	if (cinfoHash == NULL) {
		cinfoHash = &tmp_cinfoHash;
	}
*/

/*
	clustnum = createClusterInfo(treenodes, &cinfo, cinfoHash);
*/
	cinfo = cinfoData->cinfo;
	clustnum = cinfoData->clustnum;

	/** ascending order of size **/
	for (i = clustnum-1; i >= 0; i--) {
/*
	printf("checkOverlapCluster:cinfo: i=%d,root=%d\n",i,cinfo[i].root);
*/
		checkOverlapCluster(&cinfo[i], cinfoHash, nodeset);
	}
/*
	if (tmp_cinfoHash) {
		/ * destroy temporal hash if (cinfoHash == NULL) * /
		Hdestroy(tmp_cinfoHash);
	}
*/
	renumCinfo(cinfo, clustnum);	/* renum domains */
	for (i = 0; i < clustnum; i++) {
		if (numelemList(cinfo[i].members)==0) {
			continue;
		}
		reset_mergedDomains_all(&cinfo[i]);
/*
		assignSubgroups_cinfo(&cinfo[i]);
*/
	}
	return clustnum;
}

skipNeighborCheckFlag(Node *node)
{
	/*  <---L (internal)---><---CenterRoot---><---R (internal)---> */
	/* skip the check of merging neighorhood
 * 		if there is a flanking unvisited (=internal node) domain on either side */

	int skip_flag = 0;
	if (isCenterRoot(node)) {
		if (nodeDefinedParentL(node) && ! nodeVisited(node->parentL)) {
			skip_flag |= 1; /* skip left */
		}
		if (nodeDefinedParentR(node) && ! nodeVisited(node->parentR)) {
			skip_flag |= 2; /* skip right */
		}
/*
printf("Center: %d: ", skip_flag);
printNode(node);
printf("\n");
*/
	}
	return skip_flag;
}
checkOverlapCluster(ClusterInfo *cinfo, Hash *cinfoHash, NodeSet *nodeset)
{
	Neighbor *nbr;
	pList *nbrlist, *nbrlistL, *nbrlistR;
	pList *merged_list;
	Node *nbrnode;
	Node *node = cinfo->root->node;
	Node *newnode;
	ClusterInfo *nbrcinfo;
	int skip_flag = 0;

	nbrlistL = create_pList();
	nbrlistR = create_pList();
	nbrlist = create_pList();

	skip_flag = skipNeighborCheckFlag(node);

	/* create left side neighbor list */
	if (node->left && ! (skip_flag & 1)) {
		getLargeNeighbors(node->left,Opt.adjOvlpRatio,Opt.adjInclRatio,
				Opt.mincutcnt, node, 1,nbrlistL);
		while (nbr = (Neighbor *) popList(nbrlistL)) {
			nbrnode = ((Node*)nbr->node1 == node) ?
				(Node*)nbr->node2 : (Node*)nbr->node1;
			int skip_flag2 = skipNeighborCheckFlag(nbrnode);
			if (skip_flag2 & 2) { /* skip right of the counterpart */
				continue;
			}
			pushList(nbrlist, nbr);
		}
	}
	/* create right side neighbor list */
	if (node->right && skip_flag & 2) {
		getLargeNeighbors(node->right,Opt.adjOvlpRatio,Opt.adjInclRatio,
				Opt.mincutcnt, node, 1,nbrlistR);
		while (nbr = (Neighbor *) popList(nbrlistR)) {
			nbrnode = ((Node*)nbr->node1 == node) ?
				(Node*)nbr->node2 : (Node*)nbr->node1;
			int skip_flag2 = skipNeighborCheckFlag(nbrnode);
			if (skip_flag2 & 1) { /* skip left of the counterpart */
				continue;
			}
			pushList(nbrlist, nbr);
		}
	}
	free_pList(nbrlistL);
	free_pList(nbrlistR);

	while (nbr = (Neighbor *) popList(nbrlist)) {
		nbrnode = ((Node*)nbr->node1 == node) ?
			(Node*)nbr->node2 : (Node*)nbr->node1;
/*
if (nbrnode->id==3075055){
printf("rr=%d\n",isRoot(nbrnode));
}
*/
		if (cinfoHsearch(cinfoHash, nbrnode->id, &nbrcinfo)==0){
			fprintf(stderr, "node id not found: %d\n",nbrnode->id);
			fprintf(stderr, "(%d,%s)\n", node->id,node->name);
			continue;
		}

		/** merging (deleted) node: find undeleted node **/
		while (! nbrcinfo->clusters && nbrcinfo->root) {
			int origid = nbrcinfo->root->node->id;
			cinfoHsearch(cinfoHash, origid, &nbrcinfo);
			if(origid==nbrcinfo->root->node->id) break;
		}
		if (! nbrcinfo->clusters) {
			continue;
		}
		if (nbrcinfo->root->node->id == node->id) {
			/** already merged ? **/
			continue;
		}
		if (mergeClusters(cinfo, nbrcinfo, nodeset)) {
			HENTRY hent;
			Node *newnode = cinfo->root->node;
			hent.key = (char *) (newnode->id);
			hent.datum = (char *) nbrcinfo;
			HIsearch(cinfoHash, &hent, ENTER);

			break;
		}
	}
	free_pList(nbrlist);
}
/* check outer root for NEW_OUTGROUP_MODE */
checkOuterRoot(Node *node)
{
	if (IS_NEW_OUTGROUP_MODE(Opt.outgroupMode)) {
		if (isOuterRoot(node)) {
			return 1;
		} else {
			return 0;
		}
	}
	/* always OK unless NEW_OUTGROUP_MODE */
	return 1;
}

checkIncludedClustersAll(NodeSet *nset)
{
	int mergecnt;
	Node *node;
	listIter iter;
	int i;

	do {
		/** ascending order of size **/
		mergecnt = 0;
		for (i = nset->nodenum-1; i>=nset->leafnum; i--){
			node=getNode(nset,i);
			if (isRoot(node) && node->child) {
				mergecnt += checkIncludedCluster(node);
			}
		}
	} while (mergecnt);
}
checkIncludedCluster(Node *node)
{
	Node *nbrnode;
	Neighbor *nbr = NULL;
	int mincnt;
	int mergecnt = 0;
	pList *nbrlist;

	mincnt = node->cnt;
	nbrlist = create_pList();

	if (node->left) {
		getLargeNeighbors(node->left, 0.0, 1.0, 0,
				node,1,nbrlist);
	}
	while (nbr = (Neighbor *) popList(nbrlist)) {
		nbrnode = ((Node *) nbr->node1 == node) ?
			(Node *) nbr->node2 : (Node *) nbr->node1;

		if (concatNeighbor(node, nbrnode, 'L')) {
			mergecnt++;
		}
	}
	if (node->right) {
		getLargeNeighbors(node->right, 0.0, 1.0, 0,
				node,1,nbrlist);
	}
	while (nbr = (Neighbor *) popList(nbrlist)) {
		nbrnode = ((Node *) nbr->node1 == node) ?
			(Node *) nbr->node2 : (Node *) nbr->node1;
		if (node->cnt != nbr->count) continue;
		if (concatNeighbor(node, nbrnode, 'R')) {
			mergecnt++;
		}
	}
	free_pList(nbrlist);
	return mergecnt;
}
concatNeighbor(Node *nodeS, Node *nodeL, char dir)
{
	if(nodeS->cnt > nodeL->cnt) {
		/* swap*/
		Node *tmp = nodeS;
		nodeS = nodeL; nodeL = tmp;
		dir = (dir == 'R') ? 'L' : 'R';
	}

	if (checkOuterRoot(nodeS) == 0) {
		/* merged node should be outer_root if NEW_OUTGROUP_MODE */
		return 0;
	}

	if (checkNodeLen(nodeS)==2) {
		return 0;
	}

/*
printf("Merge==>");printNode(nodeS);printf(" ");printNode(nodeL);printf("%c; %d,%d\n",dir,nodeS->cnt,nodeL->cnt);
*/

	if (Opt.DEBUG & DBG_nbrrestore) {
		printf("DEL: ");
		printNode(nodeS); printf("=%c=>",dir);
		printNode(nodeL); putchar('\n');
	}

	nodeL->len += nodeS->len;
	nodeS->len = 0;
	/* mark the smaller side node as 'MERGED' */
	setFlagNode(nodeS, NODE_MERGED);

	deleteNode(nodeS);

	/* mark the smaller node according to the side */
	setMergeDir(nodeS, dir);
	return 1;
}

struct assignInfo {
	int max_subgrp;
	int excheck;
	DJset *djset;
};

mergeHomgroups_cinfo(ClusterInfoData *cinfoData, pList *clustRootsHom)
{
	_mergeHighergroups_cinfo(cinfoData, clustRootsHom, 0);
}
mergeOutergroups_cinfo(ClusterInfoData *cinfoData, pList *clustRootsHom)
{
	_mergeHighergroups_cinfo(cinfoData, clustRootsHom, 1);
}
_mergeHighergroups_cinfo(ClusterInfoData *cinfoData, pList *clustRootsHom, int outer_mode)
{
	listIter iter;
	TreeNode *tnode;
	SubClusterInfo **clustRootsHom_array;
	int num_homclust = numelemList(clustRootsHom);
	int i, idx, repr_idx;
	HENTRY hent;
	DJset *djset;
	ClusterInfo *cinfo = cinfoData->cinfo;
	ClusterInfo *rootinfo;
	int cinfonum = cinfoData->clustnum;
	int cn = (cinfonum % 2 == 0) ? cinfonum+1 : cinfonum;
	Hash *homClustHash = Hcreate(cn * 53);

	if ( (clustRootsHom_array = (SubClusterInfo**) list2array(clustRootsHom)) == NULL ) {
		fprintf(stderr, "Can't alloc memory for clustRootsHom\n"); exit(1);
	}
	for (i = 0; i < num_homclust; i++) {
		setListIter(&iter, clustRootsHom_array[i]->members, 1);
		/* store root node of each orthoCluster as a key and the ClusterInfo of homCluster that the orthoCluster belonging to as a datum in a hash */
		while (tnode = (TreeNode*) getListIter(&iter)) {
/*
printf("homclust_hash:%d:",i); printTreeNode(tnode); printf(" "); printTreeNode(clustRootsHom_array[i]->root); printf("\n");
*/
			hent.key = (char *) (tnode->node->id);
			hent.datum = (char *) &clustRootsHom_array[i];
			HIsearch(homClustHash, &hent, ENTER);
		}
	}

	/* merging homologous groups */
	djset = create_DJset(num_homclust);
	for (i = 0; i < cinfonum; i++) {
		if (! cinfo[i].clusters) continue;
		setListIter(&iter, cinfo[i].clusters, 1);
		repr_idx = -1;
		while (rootinfo = (ClusterInfo *) getListIter(&iter)) {
			tnode = rootinfo->origroot;
			hent.key = (char *) (tnode->node->id);
			idx = -1;
			if (HIsearch(homClustHash, &hent, FIND)==0) {
				printf("Data not found in hash %d,%s\n", tnode->node->id,tnode->node->name);
				fprintf(stderr, "Data not found in hash %d,%s\n", tnode->node->id,tnode->node->name);
				continue;
			}
			idx = (SubClusterInfo**) hent.datum - clustRootsHom_array;
			if (repr_idx < 0) {
				repr_idx = idx;
				/* homclust/outer group index */
				if (outer_mode) {
					cinfo[i].outerclust = idx;
				} else {
					cinfo[i].homclust = idx;
				}
			} else {
				DJset_merge(djset, repr_idx, idx);
			}
		}
	}
	for (i = 0; i < cinfonum; i++) {
		if (! cinfo[i].clusters) continue;
/*
		printf("hom[%d]=repr[%d]=%d\n", i, idx, repr_idx);
*/
		if (outer_mode) {
			idx =  cinfo[i].outerclust;
			if (idx >= 0) {
				repr_idx = DJset_find_repr(djset,idx);
				cinfo[i].outerclust = repr_idx;
			}
		} else {
			idx =  cinfo[i].homclust;
			repr_idx = DJset_find_repr(djset,idx);
			cinfo[i].homclust = repr_idx;
		}
		if (idx >= 0 && idx != repr_idx) {
			Edge *e1 = clustRootsHom_array[repr_idx]->root->node->child;
			Edge *e2 = clustRootsHom_array[idx]->root->node->child;
			if (! e1 ||
				( e2 && BETTER(MEASURE(e1),MEASURE(e2)) ) ) {
				clustRootsHom_array[repr_idx]->root
					= clustRootsHom_array[idx]->root;
			}
		}
	}
	if (! outer_mode) {
		for (i = 0; i < cinfonum; i++) {
			if (! cinfo[i].clusters) continue;
			repr_idx = cinfo[i].homclust;
			cinfo[i].homclustRoot = clustRootsHom_array[repr_idx]->root->node;
		}
	}
	free(clustRootsHom_array);
}
sortCinfoData_byHomClust(ClusterInfoData *cinfoData)
{
	int i;
	qsort(cinfoData->cinfo_link, cinfoData->clustnum, sizeof(ClusterInfo*), cmpr_cinfo_by_homclust);
/*
	for (i = 0; i < cinfoData->clustnum; i++) {
		printCinfo(cinfoData->cinfo_link[i]);
	}
*/
}

assignSubgroups_cinfo_all(ClusterInfoData *cinfoData)
{
	int i;
	ClusterInfo * cinfo = cinfoData->cinfo;
	int clustnum = cinfoData->clustnum;

	for (i = 0; i < clustnum; i++) {
		if (numelemList(cinfo[i].members)==0) {
			continue;
		}
		assignSubgroups_cinfo(&cinfo[i]);
	}
}
assignSubgroups_cinfo(ClusterInfo *cinfo)
{
	listIter iter, iter2, iter3;
	TreeNode *tnode;
	Domain *dom;
	pList *ingrps;
	SubClusterInfo *outgrp;
	int subid = 1, max_subgrp = 0;
	int origclustnum = 0;
	int i;
	SubClusterInfo **subgrp_space;
	pList *tmp_inList, *tmp_outList;
	DJset *djset;
	struct assignInfo ainfo;
	ClusterInfo *cinfo0;


	ingrps = create_pList();
	outgrp = (IS_OUTGROUP_MODE(Opt.outgroupMode)) ? createSubClusterInfo(NULL) : NULL;
	tmp_inList = create_pList();
	tmp_outList = create_pList();

	djset = create_DJset(MAXSUBGRP);


	ainfo.max_subgrp=0;
	ainfo.djset = djset;

	setListIter(&iter, cinfo->clusters, 1);
/*
	while (tnode = (TreeNode *) getListIter(&iter)) {
*/
	while (cinfo0 = (ClusterInfo *) getListIter(&iter)) {
		tnode = cinfo0->origroot;
/*
printf("assignSub:root");printTreeNode(tnode);printf("\n");
*/
		/* excheck should be skipped at the first time */
		ainfo.excheck = ( (origclustnum++ == 0) ? 0 : 1 );
		if (IS_METAGENOME_MODE(Opt.outgroupMode)) {
			/** OBSOLETE **/
			collectCluster_outgroup_meta(tnode, ingrps, outgrp);
		} else {
			collectCluster_outgroup(tnode, ingrps, outgrp);
		}
	}
/*
printf("assignSub:ingrpsize="); printTreeNode(tnode); printf("%d\n",numelemList(ingrps));
*/
	assignSubgroups_sub(ingrps, outgrp, &ainfo, tmp_inList, tmp_outList);

	/* ingrp */
	subgrp_space = (SubClusterInfo **) list2array(ingrps);
	max_subgrp = numelemList(ingrps);

	for (i = 0; i < max_subgrp; i++) {
		clearList(subgrp_space[i]->members);
	}
	
	/* tmp_inList is a unique ingroup domain list */
	setListIter(&iter, tmp_inList, 1);
	while (dom = (Domain *) getListIter(&iter)) {
		/* reassign subgroup id according to the slink clustering */
		subid = DJset_find_repr(djset,dom->subgrp);
/*
printf("reassignSub:dom=");printDomain(dom);printf(":%d\n",subid);
*/
		pushList(subgrp_space[subid-1]->members, dom);
	}
	cinfo->ingroups = create_pList();
	subid = 1;
	for (i = 0; i < max_subgrp; i++) {
		int repr_idx = DJset_find_repr(djset, i + 1);
/*
printf("repr=%d,%d;%d\n",repr_idx,i,numelemList(subgrp_space[i]->members));
*/
		if (i != repr_idx - 1) {

/** 090729 ???? */

if (subgrp_space[i]->root == NULL) {
printf("Warning: root variable has not been set: %d\n", i);
	continue;
} else if (subgrp_space[i]->root->node == NULL) {
printf("Warning: root node is not properly set: %d\n", i);
	continue;
}
			Edge *e1 = subgrp_space[i]->root->node->child,
				*e0 = subgrp_space[repr_idx-1]->root->node->child;
			if ( e0 && e1 && BETTER(MEASURE(e0), MEASURE(e1)) ) {
/***
printNode(subgrp_space[i]->root->node);
printNode(subgrp_space[repr_idx-1]->root->node);
printf("%d:%f,%d:%f\n",repr_idx,MEASURE(e0),i+1,MEASURE(e1));
***/

				/* replace with the larger distance */
				subgrp_space[repr_idx-1]->root = subgrp_space[i]->root;
			}
			if (numelemList(subgrp_space[i]->members) > 0) {
				/** never come here **/
				printf("error: %d,%d\n",i,repr_idx);
			}
		} else if (numelemList(subgrp_space[i]->members) > 0) {
			/* reassign subgroup id for each domain */
			subgrp_space[i]->clustid = subid;
			setListIter(&iter, subgrp_space[i]->members, 1);
			while (dom = (Domain *) getListIter(&iter)) {
				dom->subgrp = subid;
			}
			pushList(cinfo->ingroups, subgrp_space[i]);
			subid++;
		}
	}

	/* outgrp */
	cinfo->outgroup = createSubClusterInfo(NULL);
	setListIter(&iter, tmp_outList, 1);
	while (dom = (Domain *) getListIter(&iter)) {
		pushList(cinfo->outgroup->members, dom);
	}
	free_DJset(djset);
	free(subgrp_space);
}

assignSubgroups_sub(pList *ingrps, SubClusterInfo *outgrp,
	struct assignInfo *ainfo,
	pList *tmp_inList, pList *tmp_outList)
{
	listIter iter2, iter3;
	SubClusterInfo *subInfo;
	Domain *dom;
	int subid = 0;

	if (ingrps) {
		setListIter(&iter2, ingrps, 1);
		/** assgin subid for each subcluster **/
		while (subInfo = (SubClusterInfo *) getListIter(&iter2)) {
			setListIter(&iter3, subInfo->members, 1);
			subid++;

			while (dom = (Domain *) getListIter(&iter3)) {
				if (! dom->subgrp) {
					if (tmp_inList) {
						/* newly classified domains */
						pushList(tmp_inList, dom);
					}
					dom->subgrp = subid;
				} else if (dom->subgrp > 0) {
					if (subid != dom->subgrp){
						DJset_merge(ainfo->djset, subid, dom->subgrp);
					}
				}
			}
		}
		ainfo->max_subgrp = subid;
	}
	/** outgroup **/
	if (outgrp) {
	    setListIter(&iter2, outgrp->members, 1);
	    while (dom = (Domain *) getListIter(&iter2)) {
		if (! dom->subgrp) {
			if (tmp_outList) {
				/* newly classified domains */
				pushList(tmp_outList, dom);
			}
			dom->subgrp = -1;
		}
	    }
	}
}
SubClusterInfo *createSubClusterInfo(TreeNode *root)
{
	return createSubClusterInfo_general(root, SubCluster);
}
SubClusterInfo *createHomClusterInfo(TreeNode *root)
{
	return createSubClusterInfo_general(root, HomCluster);
}
SubClusterInfo *createOuterClusterInfo(TreeNode *root)
{
	return createSubClusterInfo_general(root, OuterCluster);
}
SubClusterInfo *createSubClusterInfo_general(TreeNode *root, ClusterLevel level)
{
	SubClusterInfo * retInfo = createSubClusterInfoN(1);
	retInfo->root = root;
	retInfo->level = level;
	return(retInfo);
}
SubClusterInfo *createSubClusterInfoN(int nel)
{
	int i;
	SubClusterInfo *retInfo;
	if ( (retInfo = (SubClusterInfo *)
		malloc(sizeof(SubClusterInfo) * nel) ) == NULL) {
		fprintf(stderr, "Can't allocate memory: SubClusterInfo\n");
		return NULL;
	}
	for (i = 0; i < nel; i++) {
		retInfo[i].members = create_pList();
		retInfo[i].root = NULL;
	}
	return retInfo;
}

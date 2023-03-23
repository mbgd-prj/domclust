/*
 * DomClust: Hierarchical Clustering for Orthologous Domain Classification
 * Copyright (c) 2000-2007, Ikuo Uchiyama
 * All rights reserved.
 */
#include <stdio.h>
#include <math.h>
#include "domclust.h"
#include "util.h"
#include "spec.h"
#include "plist.h"


/* collect cluster members by traversing the tree */
collectCluster(TreeNode *root, pList *nodelist)
{
/*
	printf("collect:root");printTreeNode(root);printf("\n");
*/
	collectCluster_sub(root, nodelist);
	clearDomainMark(nodelist);
}
collectCluster_sub(TreeNode *tnode, pList *nodelist)
{
	Domain *dom;
	if (tnode == NULL) {
	} else if (ClustTree_isLeaf(tnode)) {
		/* leaf node */
/*
	printf("collect:Dom");printDomain(tnode->dom);printf(":%d\n",tnode);
*/
		pushList(nodelist, tnode->dom);
	} else {
		collectCluster_sub(tnode->child1, nodelist);
		collectCluster_sub(tnode->child2, nodelist);
	}
}
collectCluster_outgroup(TreeNode *root, pList *subclstList,
	SubClusterInfo *outgrpSubClust)
{
	listIter iter;

	SubClusterInfo *ingrp;
/*
 	printf("OutRoot");printNode(root->node);printf("\n");
*/
	collectCluster_outgroup_sub(root, subclstList, outgrpSubClust, 0);
	if (subclstList) {
		setListIter(&iter, subclstList, 1);
		while (ingrp = (SubClusterInfo *) getListIter(&iter)) {
			clearDomainMark(ingrp->members);
		}
	}
	if (outgrpSubClust) {
		clearDomainMark(outgrpSubClust->members);
	}
}
collectCluster_outgroup_sub(TreeNode *tnode, pList *subclstList,
			SubClusterInfo *outgrpSubClust, int outRootFlag)
{
	Domain *dom;
	int domn;
	int flg;
	OutGrpStat outchk;
	Node *node;
	int cut;
OutGrpStat outchk_node;


	if (tnode == NULL) {
		return 0;
	}
	node = tnode->node;

	if (! ClustTree_isLeaf(tnode)) {

		outchk = outgroupCheck_TreeNode(tnode);
/*
printf("collect:Internal:");printNode(node);printf(":%d\n",outchk);
*/

/*
outchk_node = outgroupCheck(tnode->node);
if (outchk != outchk_node) {
printf("###@###");
}
printf("outchk=%d,%d: ",outchk,outchk_node);
printNode(node); printf("\n");
outchk = outchk_node;
*/


/****
		if (! outRootFlag) {
			if (isTreeRoot(tnode)) {
				outRootFlag = 1;
printf("root: "); printNode(tnode->node);printf("\n");
			} else {
				cut = checkOutGrpRootTreeNode(tnode, &outchk);
				if (cut == 0) {
printf("CUTCUT:");
printNode(tnode->node);
printf("\n");
					outRootFlag = 1;
					makeRootTreeNode(tnode);
				}
			}
		} else if (isTreeRoot(tnode)) {
		}
*/


/*
		if (outchk == In_In) {
*/
		if (outchk == In_In || outchk == In_In_Horiz) {
/*
printf("collect:internal:IN_IN:");printTreeNode(tnode);printf("\n");
*/
			/* both are ingroup nodes */
/*`
			if (isAllRoot(node) || ! nodeVisited(node)) {
*/
/*
if (isAllInTreeRoot(tnode) != isAllInRoot(node)) {
	printf("AllInRoot: %d,%d;", isAllInTreeRoot(tnode), isAllInRoot(node));
	printf("flag=%d,%d;",tnode->flag, node->flag);
	printTreeNode(tnode); printf("\n");
}
*/
/*
			if (duplicationCheck_inGrp_TreeNode(tnode)==0) {
			if (ingroupFlagCheck(tnode->spflag) == 2) {
			if (isAllInRoot(tnode->node)) {
*/
			if (isAllInTreeRoot(tnode)) {
/*
printf("collect:internal:InRoot:");printNode(node);printf("\n");
*/
				SubClusterInfo *ingrp = createSubClusterInfo(tnode);
/*
printf("collectOutGrp:AllInRoot>>");printTreeNode(tnode);printf(":%d\n",tnode);
*/
				collectCluster_sub(tnode, ingrp->members);
				pushList(subclstList, ingrp);
			} else {
				collectCluster_outgroup_sub(tnode->child1,
					subclstList, outgrpSubClust,
					outRootFlag);
				collectCluster_outgroup_sub(tnode->child2,
					subclstList, outgrpSubClust,
					outRootFlag);
			}
		} else if (outchk == Out_Out) {
			/* both are outgroup nodes */
			collectCluster_sub(tnode->child1, outgrpSubClust->members);
			collectCluster_sub(tnode->child2, outgrpSubClust->members);
		} else if (outchk == Out_Both) {
			/* node1 is outgroup nodes */
			collectCluster_outgroup_sub(tnode->child2, subclstList,
					outgrpSubClust, outRootFlag);
			collectCluster_sub(tnode->child1, outgrpSubClust->members);
		} else if (outchk == Both_Out) {
			/* node2 is outgroup nodes */
			collectCluster_outgroup_sub(tnode->child1, subclstList,
					outgrpSubClust, outRootFlag);
			collectCluster_sub(tnode->child2, outgrpSubClust->members);
		} else if (outchk == Both_Both) {
			/** outgrp spec in a subcluster: possible HGT **/
			/** never come here unless -Ohorizweight is specified*/
			SubClusterInfo *ingrp = createSubClusterInfo(tnode);
			collectCluster_sub(tnode, ingrp->members);
			pushList(subclstList, ingrp);
		}
	} else {
		/* leaf node */
		if (tnode->dom) {
			dom = tnode->dom;

			flg = outgroupFlagCheck(node->spflag);
			if (flg == 0) {
				/* ingroup */
				SubClusterInfo *ingrp = createSubClusterInfo(tnode);
				pushList(ingrp->members, dom);
				pushList(subclstList, ingrp);
			} else {
				/* outgroup */
				pushList(outgrpSubClust->members, dom);
			}
		}
	}
}

/** OBSOLETE!! **/
collectCluster_outgroup_meta_all(pList *clustRoots, pList **newClustRoots)
{
	
	listIter iter, iter2;
	TreeNode *tnode;
	SubClusterInfo *subInfo;
	pList *subclstList = create_pList();
	char delflag;

	setListIter(&iter, clustRoots, 1);
	*newClustRoots = create_pList();
Node *tmpnode;
	while (tnode = (TreeNode *) getListIter(&iter)) {
		collectCluster_outgroup_meta(tnode, subclstList);
		delflag = 1;
		setListIter(&iter2, subclstList, 1);
		while (subInfo = (SubClusterInfo *) getListIter(&iter2)) {
			pushList(*newClustRoots, subInfo->root);
			if (subInfo->root == tnode) {
				delflag = 0;
				continue;
			}
			setFlagNode(subInfo->root->node, NODE_TMPMARK);
		}
		if (delflag) {
			deleteRootNode(tnode->node);
			phyloCut_without_check(tnode->node, NULL);
		}
		while (subInfo = (SubClusterInfo *) popList(subclstList)) {
			setFlagNode(subInfo->root->node, NODE_TMPMARK);
		}
	}
}
collectCluster_outgroup_meta(TreeNode *root, pList *subclstList)
{
	listIter iter;

	SubClusterInfo *ingrp;
	SubClusterInfo *outgrpSubClust = createSubClusterInfo(NULL);
	SubClusterInfo *new_outgrpSubClust = NULL;
	collectCluster_outgroup_meta_sub(root, subclstList, outgrpSubClust, &new_outgrpSubClust, 1);
	setListIter(&iter, subclstList, 1);

	while (ingrp = (SubClusterInfo *) getListIter(&iter)) {
		clearDomainMark(ingrp->members);
	}

	if (numelemList(outgrpSubClust->members) > 0) {
		pushList(subclstList, outgrpSubClust);
	}
}
/* MetaGenomeMode: outgroup == unk-only (meta) node */
collectCluster_outgroup_meta_sub(TreeNode *tnode, pList *subclstList,
			SubClusterInfo *outgrpSubClust,
			SubClusterInfo **new_outgrpSubClust,
			int outRootFlag)
{
	Domain *dom;
	int domn;
	int flg;
	OutGrpStat outchk;
	Node *node;
/*
	pList *new_outgrp = outgrp;
*/

	if (tnode == NULL) {
		return 0;
	}
	node = tnode->node;

	if (! ClustTree_isLeaf(tnode)) {

		outchk = outgroupCheck_TreeNode(tnode);
/*
outchk = outgroupCheck(tnode->node);
*/

/*
if (node){
  printNode(node);
  if (node->child){
	printf(" %f", node->child->dist);
  }
  printf(": %d\n", outchk);
}
*/

		if (outchk == In_In) {
			/* both are ingroup nodes */
			if (isAllInRoot(node)) {
				SubClusterInfo *ingrp = createSubClusterInfo(tnode);
				collectCluster_sub(tnode, ingrp->members);
				pushList(subclstList, ingrp);
				*new_outgrpSubClust = NULL;
				return 1;
			} else {
				collectCluster_outgroup_meta_sub(tnode->child1,
					subclstList, outgrpSubClust,
					new_outgrpSubClust,
					outRootFlag);
				collectCluster_outgroup_meta_sub(tnode->child2,
					subclstList, outgrpSubClust,
					new_outgrpSubClust,
					outRootFlag);
				*new_outgrpSubClust = NULL;
			}
		} else if (outchk == Out_Out) {
			/* both are outgroup (unk only) nodes */
			if (*new_outgrpSubClust == NULL) {
				*new_outgrpSubClust = createSubClusterInfo(tnode);
				pushList(subclstList, *new_outgrpSubClust);
			}
			collectCluster_sub(tnode->child1, (*new_outgrpSubClust)->members);
			collectCluster_sub(tnode->child2, (*new_outgrpSubClust)->members);
		} else if (outchk == Out_Both) {
			/* node1 is outgroup (unk only) nodes */
			int ret = collectCluster_outgroup_meta_sub(tnode->child2,
					subclstList, outgrpSubClust,
					new_outgrpSubClust,
					outRootFlag);
			if (ret == 1) {
				/* sister is a single OG -- put into this OG */
				*new_outgrpSubClust = (SubClusterInfo*) getListIdx(subclstList, -1);
				(*new_outgrpSubClust)->root = tnode;
			    } else {
				/* multiple OGs -- put into another OG */
				if (*new_outgrpSubClust == NULL) {
					*new_outgrpSubClust = createSubClusterInfo(tnode->child1);
					pushList(subclstList, *new_outgrpSubClust);
				}
			}
			collectCluster_sub(tnode->child1, (*new_outgrpSubClust)->members);
			return ret;
		} else if (outchk == Both_Out) {
			/* node2 is outgroup (meta) nodes */
			int ret = collectCluster_outgroup_meta_sub(tnode->child1,
					subclstList, outgrpSubClust,
					new_outgrpSubClust,
					outRootFlag);
			if (ret == 1) {
				/* sister is a single OG -- put into this OG */
				*new_outgrpSubClust = (SubClusterInfo*) getListIdx(subclstList, -1);
				(*new_outgrpSubClust)->root = tnode;
			    } else {
				/* multiple OGs -- put into another OG */
				if (*new_outgrpSubClust == NULL) {
					*new_outgrpSubClust = createSubClusterInfo(tnode->child2);
					pushList(subclstList, *new_outgrpSubClust);
				}
			}
			collectCluster_sub(tnode->child2, (*new_outgrpSubClust)->members);
			return ret;
		} else if (outchk == Both_Both) {
			/** outgrp spec in a subcluster: possible HGT **/
			/** never come here unless -Ohorizweight is specified*/
			SubClusterInfo *ingrp = createSubClusterInfo(tnode);
			collectCluster_sub(tnode, ingrp->members);
			pushList(subclstList, ingrp);
		}
	} else {
		/* leaf node */

		if (tnode->dom) {
			dom = tnode->dom;

			flg = outgroupFlagCheck(node->spflag);
			if (flg == 0) {
				/* ingroup */
				SubClusterInfo *ingrp = createSubClusterInfo(tnode);
				pushList(ingrp->members, dom);
				pushList(subclstList, ingrp);
				return 1;
			} else {
				/* outgroup */
				pushList(outgrpSubClust->members, dom);
			}
		}
	}
	return 0;
}

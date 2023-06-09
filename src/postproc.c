/*
 ;* DomClust: Hierarchical Clustering for Orthologous Domain Classification
 * Copyright (c) 2000-2007, Ikuo Uchiyama
 * All rights reserved.
 */

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <limits.h>
#include "domclust.h"
#include "util.h"
#include "neighbor.h"
#include "spec.h"

enum {MODE_NODE, MODE_TREENODE} TreeNodeMode;

cmpr_node_by_size(Node **a, Node **b)
{
	return (*b)->cnt - (*a)->cnt;
}
cmpr_treenode_by_size(TreeNode **a, TreeNode **b)
{
	return ( (*b)->node->cnt == (*a)->node->cnt ) ?
		(*b)->node->id - (*a)->node->id :
		( (*b)->node->cnt - (*a)->node->cnt );
}
cmpr_list_by_size(pList **a, pList **b)
{
	return (*b)->numelem - (*a)->numelem;
}

pList *createRootNodeList(NodeSet *nodes);

outputResults(SimGraph *SimG, int argc, char **argv)
{
	NodeSet *nodes = SimG->nodes;
	EdgeSet *edges = SimG->edges;
	pList *origClustRoots;
	pList *clustRoots=NULL, *clustRootsHom=NULL, *outerRoots=NULL;
	ClusterInfoData cInfoData;
	int sortflag = 0;

/*
	ClusterInfo *cInfo;
	int i, clustnum;
	Hash *cinfoHash = NULL;
*/

	setTotalNodeNum(nodes);

	if (Opt.outstyle == DUMP) {
		dumpGraph(SimG, argc, argv);
		exit(0);
	}

	if (Opt.spmaskStr) {
		checkSkipMaskAll(nodes);
	}

	/* save original roots (homolog clusters) */

	if (Opt.homClustOut || IS_OUTGROUP_MODE(Opt.outgroupMode)) {
		origClustRoots = createRootNodeList(nodes);
		sortList(origClustRoots, cmpr_node_by_size);
	}


	/** phylogenetic tree cutting procedure **/
        phyloCutAll(nodes, edges);


	deleteOutgroupClusters(nodes);


	/** reexamine neighboring clusters **/
	if (Opt.adjInclRatio >= ONE) {
		checkIncludedClustersAll(nodes);
	}

	if (Opt.delete_small) {
		/** delete small clusters before domain definition **/
		checkClusterSize(nodes);
	}
	/** domain definition **/
	checkDomains(nodes);

	if (Opt.homClustOut || IS_NEW_OUTGROUP_MODE(Opt.outgroupMode)) {
		convClustTree_fromHom(nodes, origClustRoots, &clustRootsHom, &clustRoots);
		sortList(clustRootsHom, cmpr_list_by_size);
	} else {
		clustRoots = convClustTree(nodes);
	}

	clearDomainMarkAll(nodes);



/**
{
 listIter iter;
 TreeNode *treenode;
 setListIter(&iter, clustRoots, 1);
 while (treenode = (TreeNode *) getListIter(&iter)) {
	if (treenode->node->id==94527||treenode->node->id==93118) {
		printf("##>>"); printNode(treenode->node);printf("\n");
	}
	if (! isRoot(treenode->node)) {
printf("isnot_root: "); printTreeNode(treenode); printf("\n");
	}
 }
}
**/


/**
	clustRoots = convClustTree(nodes);

print_clustList(clustRoots, "00>>>")

	if (IS_METAGENOME_MODE(Opt.outgroupMode) && Opt.outputSubClustTree) {
		pList *newClustRoots;
		collectCluster_outgroup_meta_all(clustRoots, &newClustRoots);
		free_pList(clustRoots);
		clustRoots = newClustRoots;

print_clustList(clustRoots, "01>>>")

		checkDomains(nodes);
	}
**/

/*
	if (Opt.homClustOut) {
		convClustTree_fromHom(nodes, origClustRoots, &clustRootsHom, NULL);
		sortList(clustRootsHom, cmpr_list_by_size);
	}
*/


if (clustRootsHom) {
	listIter iter;
	SubClusterInfo *hom_cInfo;
        setListIter(&iter, clustRootsHom, 1);
        while (hom_cInfo = (SubClusterInfo*) getListIter(&iter)) {
                TreeNode *root = hom_cInfo->root;
/*
                if (root->node->id == 92988) {
                        printNode(root->node);
                        printf("; %d,%d\n", root->node->child->node1, root->node->child->node2);
                }
*/
        }
}

	if (IS_NEW_OUTGROUP_MODE(Opt.outgroupMode)) {
		addRootNode(clustRootsHom, clustRoots);
	}
	createClusterInfo_All(clustRoots, &cInfoData);
	if (IS_NEW_OUTGROUP_MODE(Opt.outgroupMode)) {
		addOuterGroup(clustRootsHom, cInfoData.cinfoHash, &outerRoots);
	}

	sortList(clustRoots, cmpr_treenode_by_size);

	if (Opt.adjOvlpRatio || (Opt.mincutcnt > 1) ||
			(Opt.adjInclRatio && (Opt.adjInclRatio < ONE)) ) {
		/* merging adjacent clusters */
		checkOverlapClustersAll(clustRoots, &cInfoData, nodes);
	}

	assignSubgroups_cinfo_all(&cInfoData);

	assign_clusterid(&cInfoData);


	if (Opt.homClustOut) {
		/** assign and merge homologous groups **/
		mergeHomgroups_cinfo(&cInfoData, clustRootsHom);
		sortflag = 1;
	}
	if (IS_NEW_OUTGROUP_MODE(Opt.outgroupMode)) {
		mergeOutergroups_cinfo(&cInfoData, outerRoots);
		sortflag = 1;
	}
	if (sortflag) {
		sortCinfoData_byHomClust(&cInfoData);
	}
		
	if (Opt.taxMapOut != 3) {
		if (IS_NEW_OUTGROUP_MODE(Opt.outgroupMode) && Opt.outputOuterClustTree) {
			outputOuterClusterInfo(outerRoots, nodes);
		} else {
			outputClusterInfo(&cInfoData, nodes);
		}
	}
	if (Opt.taxMapOut > 1) {
		mapTaxInfoAll(&cInfoData);
	}
}
/** for debug **/
print_clustList(pList *clustList, char *string)
{
	listIter iter;
	TreeNode *tn;
	setListIter(&iter, clustList, 1);
	while (tn = (TreeNode*)getListIter(&iter)) {
		if (string) printf(string);
		printNode(tn->node); putchar('\n');
	}
}

assign_clusterid(ClusterInfoData *cInfoData)
{
	int i, clsize;
	double spcnt;
	int clustid = 0;
	ClusterInfo *cInfo = cInfoData->cinfo;
	int clustnum = cInfoData->clustnum;

	for (i = 0; i < clustnum; i++) {
		clsize = numelemList(cInfo[i].members);
		spcnt = spFlagCntW_All(cInfo[i].spflag);
		if ( (clsize == 0) ||
			(spcnt < Opt.minsp || clsize < Opt.minent) ||
			/* small clusters */
			(Opt.outgroupMode &&
				outgroupFlagCheck(cInfo[i].spflag)==2) ) {
			/* only outgroup species */

			cInfo[i].clustid = -1;
			continue;
		}
		cInfo[i].clustid = ++clustid;
	}
}


#ifdef EXPERIMENTAL
checkConnect(Node *node) {
	float cutcnt;
	int totcnt = node->child->node1->cnt * node->child->node2->cnt;
	if (Opt.chkConnect > 1) {
		/** count **/
		cutcnt = Opt.chkConnect;
		if (cutcnt > totcnt) {
			cutcnt = totcnt;
		}
	} else {
		/** ratio **/
		cutcnt = ((float)node->child->node1->cnt *
			node->child->node2->cnt) * Opt.chkConnect;
		cutcnt = (cutcnt > MAX_CONN)  ? MAX_CONN : cutcnt;
	}
/*
printf("cutcnt>>>%f,%d (%d,%d)\n",cutcnt,(int) node->child->connect, (int) node->child->node1->cnt, (int) node->child->node2->cnt);
*/
	if (node->child->connect < cutcnt) {
		/* cut */
		return 0;
	} else {
		return 1;
	}
}
#endif

#define DEBUG_ON(x) if(x){ tmpflag=1; } else { tmpflag=0; }
int tmpflag;

phyloCutAll(NodeSet *nodes, EdgeSet *edges)
{
	int i;
	Node *node;
	NodeSetIter iter;

	/* decreasing order:
		note that new root nodes may be created after the cut */
/*
	getRootInit(&iter, nodes, -1, 0);
*/
	getRootInit(&iter, nodes, -1, 1);

/*
	while (node = getRootNext(&iter)) {
*/
	while (node = getRootNext_all(&iter)) {

#define TMPDBG
#ifdef TMPDBG	/* for DEBUG */
DEBUG_ON(node->id == 12391)
#endif

		/** cut non-recursively
		   further cut should be done in the succeeding loops **/
		if ( checkUnvisitedParent(node) == 0) {
			/* a node in a tree under an undeleted root */
			continue;
		}
		phyloCut(node);
	}
}

phyloCut_balancetest(Node *node)
{
	double dup_cutoff = 3.0;
	Node *largeNode;
	if (node->child == NULL || isFlankNode(node)) {
		return 0;
	}
	if (node->child->node1->cnt > node->child->node2->cnt) {
		largeNode = node->child->node1;
	} else {
		largeNode = node->child->node2;
	}
	if (largeNode->cnt > spFlagCnt(largeNode->spflag) * dup_cutoff) {
		return 1;
	}
}
phyloCut_disttest(Node *node)
{
	int cut1 = 0, cut2 = 0;
	Node *childNode;
	if (node->child == NULL) {
		return 0;
	}
	if (isFlankNode(node)) {
		/* never come here? */
		return 0;
	} else {
		if (node->child->node1->child) {
			cut1 = phyloCut_disttest_sub(node, node->child->node1);
		}
		if (node->child->node1->child) {
			cut2 = phyloCut_disttest_sub(node, node->child->node2);
		}
	}
	return (cut1 || cut2);
}
phyloCut_disttest_sub(Node *node, Node *child)
{
	/* test mode */
	while (isFlankNode(child)) {
		if (isFlankNode1(child)) {
			child = child->child->node1;
		} else {
			child = child->child->node2;
		}
	}
	if (! child->child) {
		return 0;
	}
/*
printf(">>node=%s,%d; %d,%d\n", node->name,node->id, node, child);
*/
	double distr = DISTRATIO( MEASURE(child->child),
			MEASURE(node->child) );
/*
printf(">dist=%lf,%lf\n",MEASURE(child->child),MEASURE(node->child));
*/
	if (distr <= Opt.distdiffcut){
		return 1;
	} else {
		return 0;
	}
}
/* cut parent and make node as a new root */
postPhyloCut(Node *node, Node *parent)
{
	int stat;
	int outgroupflag = 0;

/*
printf("postPhyloCut:");printNode(node);printf(" ");
if(parent) printNode(parent);
printf(":%d\n",nodeDeleted0(node));
*/

if (nodeDeleted0(node)){
	printf("phylocut: ????\n");
	if (node->child) {
		printf("%d\n",isRoot(node));
		printf("%d\n",isRoot(node->child->node1));
		printf("%d\n",isRoot(node->child->node2));
	}
	return 0;
}

	if (! parent) {
		return 0;
	}
	/** When there is a parent, we have entered this node after
		cutting this parent node. So we must cancel
		the domain boundaries created during the cluster
		merging process.
	**/

#ifdef TMPDBG	/* for DEBUG */
if (tmpflag) {
printf("******\n");
printNode(node);
printNode(parent);
printf("\n");
printf("(%d,%d);(%d,%d)\n",parent->left,parent->right,node->left,node->right);
}
#endif

	/* outgroupflag: the node is located between root and ingroup-root */
	if (isOutGrpNode(parent)) {
		setOutGrpNode(node);
		outgroupflag = 1;
	}
	if (! outgroupflag) {
		if (node->parent == parent && checkNodeLen(node) == 2) {
/*
printf("makeCenter:");printNode(node);printf("\n");
*/
			makeCenterRoot(node);
		} else {
			restoreNeighborList(parent, node);
			restoreBreak(node, parent);
		}
	}

/*
if (node->id==81441) {
	printf("pp>%d,%d,%d\n",node->parent, node->parentL, node->parentR);
	if (node->parent) {
		printf("parent=");printNode(node->parent);printf("\n");
	}
	if (node->parentL) {
		printf("parentL=");printNode(node->parentL);printf("\n");
	}
	if (node->parentR) {
		printf("parentR=");printNode(node->parentR);printf("\n");
	}
}
*/

	if ( (stat=checkUnvisitedParent(node)) == 0) {
		/** do not proceed further if there remains an unvisited parent **/
		/** never come here unless the node is an outgroup node **/
		if (outgroupflag) {
			/** set the node as an ingroup root to prevent further cut 
				because the flanking node is in a single cluster **/
/*
printf(">>>outgrpflag\n");
printNode(node);
printf("\n");
printf("makeroot0:");printNode(node);printf(" ");printNode(parent);printf("\n");
*/
			makeRootNode(node);
/*
		} else {
			if ( checkUnvisitedCenterParent(node) ) {
printf(">>>CENTER:");printNode(node);printf("<\n");
				makeCenterRoot(node);
			}
*/
		}
		return 0;
	}


	/* ### There is no unvisited parent */

	/* delete the flanking root node if exist,
		which is probably retained as a long domain */
	if (! outgroupflag) {
		if (nodeDefinedParentL(node) && isRoot(node->parentL)) {
			deleteNode(node->parentL);
			restoreBreak(node, node->parentL);
		}
		if (nodeDefinedParentR(node) && isRoot(node->parentR)) {
			deleteNode(node->parentR);
			restoreBreak(node, node->parentR);
		}
	}
/*
if (node->parentL && node->parentL->id==10869) {
printf("L******%d,%d\n",parent->id,outgroupflag);
} else if (node->parentR && node->parentR->id==10869) {
printf("R******%d,%d\n",parent->id,outgroupflag);
}
*/

	/* if the parent is not a root then make this node as a root */
/***
	if (! node->parent || ! isRoot(node->parent)) {
		makeRootNode(node);
	}
***/
	makeRootNode(node);
/*
printf("makeroot:");printNode(node);printf(" ");printNode(parent);printf(";%d,%d,%d\n",isRoot(node),isAllRoot(node),isOuterRoot(node));
printf(">postPhyloCut:>makeroot:"); printNode(node); printf(" "); printNode(parent); printf(" %d\n",node->flag);
*/

/*
printFlagNode(node);
printf(": isroot=%d\n",isRoot(node));
*/

	return 0;
}
phyloCut(Node *node)
{
	int cut = 0;
	int cutO = 0;
	int outroot_flag = 0;
	OutGrpStat outchk = In_In;
/*
printf("phyloCut:");printNode(node);printf("\n");
*/

	if (node->child) {
		if (isFlankNode(node)) {
			Node *nbrNode = NULL, *child = NULL;
	
			if (IS_NEW_OUTGROUP_MODE(Opt.outgroupMode)) {
				if (checkDeletedOuterParent(node)) {
					if (checkNodeLen(node) == 2) {
						/* temporary set as an outer root node */
						/* delete the outerRoot flag if child node is not a valid ortholog node */
						setOuterRoot(node);
/*
	printf("outerRoot(flank)>>");printNode(node);printf("\n");
*/
					} else {
/*
	printf("DelOuterRoot(flank)>>");printNode(node);printf("\n");
*/
						setDeletedOuterRoot(node);
					}
				}
			}

			/* nbrNode = central node (newnode0) */
			if (isFlankNode1(node)) {
				child = node->child->node1;
			} else if (isFlankNode2(node)) {
				child = node->child->node2;
			} else {
				fprintf(stderr, "???? node type mismatch\n");
			}
			if (isOutGrpNode(node)) {
				postPhyloCut(child, node);
				if (isInRoot(node)) {
					deleteNode(node);
				}
				return 0;
			}
			/* find the neighboring (central) node */
			if (node == child->parentL) {
				/* right neighbor */
				nbrNode = (Node *) getNeighborNode(
							node, 1, NULL);
			} else if (node == child->parentR) {
				/* left neighbor */
				nbrNode = (Node *) getNeighborNode(
							node, -1, NULL);
			} else {
				/** ???? never come here **/
				/** node == child->parent || child->parentM */
				fprintf(stderr, "ERROR?? a central node is treated as a flanking node\n");
				assert(0);
				makeRootNode(child);
				deleteRootNode(node);
				return 0;
			}
			if (nbrNode && nbrNode->flag & NODE_DELETED) {
				/* the central node is already cut */
				/* never come here ? */
				fprintf(stderr, "Central node is already cut:%p\n",nbrNode);
				deleteNode(node);
				return postPhyloCut(child, node);
			} else {
/*
				if (IS_NEW_OUTGROUP_MODE(Opt.outgroupMode) &&  ! checkDeletedOuterParent(node)) {
printf("####>>>>OO>>");printNode(node);printf("\n");
					makeRootNode(node);
				} else
*/
				if (checkNodeLen(node) < 2 && checkOuterRoot(node)) {
					/* short segment */
					deleteNode(node);
					checkIncludedCluster(node);
				} else {
					/* having sufficient length */
					/* make this node as a root */
					makeRootNode(node);
				}
				return 0;
			}
			/* never come here */
		}

		if (Opt.metagenomeMode) {
			cut = phyloCheck_forMeta(node);
		} else {
	/** Test for the phylogenetic tree cutting procedure **/
			cut = duplicationCheck(node);
		}
/*
if (node->child->node1) printNode(node->child->node1);
if (node->child->node2) printNode(node->child->node2);
printf("Cut: "); printNode(node); printf(">>>cut=%d\n",cut);
*/
#ifdef EXPERIMENTAL
/*
printf("cut: ");
printEdge(node->child);
printf("\n");
*/
/*
		if ((Opt.chkConnect && ! checkConnect(node))) {
			cut = 1;
		}
*/
#endif
		if (cut == 1) {
			/** immediately cut this node **/
		} else if ( cut == 2 ) {
#ifdef EXPERIMENTAL
			if ((Opt.chkConnect && ! checkConnect(node)) ||
				(Opt.sumcut &&
			  	Opt.sumcut > getEdgeScoreSum(node->child))) {

				/* weak cut condition (connect/sumcut) */
			} else
#endif
			if (Opt.distdiffcut) {
				/* weak cut condition (distdiffcut) */
				int cut1 = phyloCut_disttest(node);
int cut2 = 0;
/*
				int cut2 = phyloCut_balancetest(node);
*/
				if (cut1 == 0 && cut2 == 0) {
					cut = 0;
				} else {
					/* cut = 2 */
				}
			} else {
				cut = 0;
			}
		}
		if (IS_NEW_OUTGROUP_MODE(Opt.outgroupMode)) {
			int outchk = outgroupCheckRaw(node);
			cutO = checkOutGrpCutNode(outchk, node);
/*
printf("cut=%d,cutO=%d\n",cut,cutO);
if (node->id==165784) {
	printf(">>");printNode(node);printf("\n");
	if (node->parent) {printf("parent>>");printNode(node->parent);printf(":%d,%d\n",isOuterRoot(node->parent),isDeletedOuterRoot(node->parent));}
	if (node->parentL) {printf("parentL>>");printNode(node->parentL);printf(":%d,%d\n",isOuterRoot(node->parentL),isDeletedOuterRoot(node->parentL));}
	if (node->parentR) {printf("parentR>>");printNode(node->parentR);printf(":%d,%d\n",isOuterRoot(node->parentR,isDeletedOuterRoot(node->parentR)));}
}
*/

			if (checkDeletedOuterParent(node)) {
/*
				if (cut == 1 || cutO == 1) {
printf("%d,%d,%d\n",cut,cutO,outchk);
if (node->id==171084) {
*/
				if ((cut == 1 || cutO == 1) && ! isInGrpNode(outchk)) {
					/* cutO == 1: in+out, in+out -- not outer root */
					setDeletedOuterRoot(node);
/*
if (node->id==14334) {
	printf("DelOuterRoot(internal)>>");printNode(node);printf("\n");
	if (node->child) {printEdge(node->child);}
}
*/
				} else {
					/* set outer root if all parents are "visited" (= above root) */
					setOuterRoot(node);
/*
if (node->parent && node->parent->id==14434) {
	printf("outerRoot(internal)>>");printNode(node);printf("\n");
}
*/
				}
/*
			} else if (isCenterRoot(node)) {
				setOuterRoot(node);
*/
			}
		}

		if (IS_OUTGROUP_MODE(Opt.outgroupMode)) {
			cutO = checkOutGrpRoot(node, &outchk);
			if (cutO) {
				cut = 1;
			} else {
				cut = post_checkOutGrpRoot(node, outchk);
			}

			if (outchk == Out_Both || outchk == Both_Out) {
				/* out,in+out || in+out,out */
				if (cut == 0) {
					/* possible outgroup root */
					/* force to cut to find ingroup root */
					cut = 3;
				}
			}
/*
printf("CutO: "); printNode(node); printf(">>>cutO=%d,outchk=%d\n",cut,outchk);
*/
		}
		if (cut) {
			cutNode(node, cut, outchk);
		}
	}
	return cut;
}
phyloCut_without_check(Node *node, Node *parent)
{
	Edge *child = node->child;
/*
	if (nodeDeleted0(node)) {
		return 0;
	}
*/
	if (testFlagNode(node, NODE_TMPMARK)) {
/*
printf("OK");
printNode(node);
printf("\n");
*/
		return 0;
	}
/*
printf("phylocut:");
printNode(node);
printf("\n");
*/
	if (child) {
		if (isFlankNode(node)) {
			Node *nbrNode = NULL, *childNode = NULL;
			if (isFlankNode1(node)) {
				childNode = node->child->node1;
			} else if (isFlankNode2(node)) {
				childNode = node->child->node2;
			}
			if (node == childNode->parentL) {
				nbrNode = (Node *) getNeighborNode(node, 1, NULL);
			} else if (node == childNode->parentR) {
				nbrNode = (Node *) getNeighborNode(node, -1, NULL);
			}
			if (nbrNode && nbrNode->flag & NODE_DELETED) {
				fprintf(stderr, "????\n");
				deleteNode(node);
			} else {
				if (checkNodeLen(node) < 2) {
					/* short segment */
					deleteNode(node);
					checkIncludedCluster(node);
				} else {
					/* having sufficient length */
					/* make this node as a root */
					makeRootNode(node);
				}
/*
				return 0;
*/
			}
			phyloCut_without_check(childNode, node);
		} else {
/*
			if (parent) {
printf("RESTORE::");
printNode(node);
printf("\n");
				restoreNeighborList(parent, node);
				restoreBreak(node, parent);
			}
*/
			unsetOutGrpNode(node);
			cutNode(node, 1, In_In);

			phyloCut_without_check(child->node1, node);
			phyloCut_without_check(child->node2, node);
		}
	}
}
cutNode(Node *node, int cut, OutGrpStat outchk)
{
	int outroot_flag = 0;
	int skip_child = 0;
	if (Opt.DEBUG & DBG_nbrrestore) {
		printf("Cut: "); printNode(node); putchar('\n');
		print_specFlag(node->child->node1->spflag);
		print_specFlag(node->child->node2->spflag);
	}
	if (! isOutGrpNode(node) && cut == 3) {
	    /** this is an outgroup root: do not delete **/
		outroot_flag = 1;
	    /* temporarily set the flag to propagate it to children */
		setOutGrpNode(node);
	} else if (isCenterRoot(node) && cut == 3) {
	} else {
/*
printf("unroot: ");printNode(node);printf("\n");
*/
		deleteRootNode(node);
	}
	if (outchk == Out_Both) {
		skip_child = 1;
	} else if (outchk == Both_Out) {
		skip_child = 2;
	}
	if (skip_child != 1) {
		postPhyloCut(node->child->node1, node);
	}
	if (skip_child != 2 &&
		(node->child->node1 != node->child->node2)) {
		/* NOT self_match */
		postPhyloCut(node->child->node2, node);
	}
	if (outroot_flag) {
		/* reset the flag in the outgroup root */
		unsetOutGrpNode(node);
	}
}

/** metagenomeMode **/
phyloCheck_forMeta(Node *node)
{
	specFlagP spflag1, spflag2;
	int cut = 0;
	if (metagenomeOnlyCluster(node->spflag)) {
		return 0;
	} else if (node->child) {
		if (isFlankNode(node)) {
			if (isFlankNode1(node)) {
				cut = phyloCheck_forMeta(node->child->node1);
			} else if (isFlankNode2(node)) {
				cut = phyloCheck_forMeta(node->child->node2);
			}
		} else {
			spflag1 = node->child->node1->spflag; 
			spflag2 = node->child->node2->spflag; 
			if (metagenomeOnlyCluster(spflag1)) {
				cut = phyloCheck_forMeta(node->child->node2);
			} else if (metagenomeOnlyCluster(spflag2)) {
				cut = phyloCheck_forMeta(node->child->node1);
			} else {
				cut = duplicationCheck(node);
			}
		}
	} else {
		/** leaf **/
		return 0;
	}
	return(cut);
}

post_checkOutGrpRoot(Node *node, OutGrpStat outchk) {
	int cut;
	if (outchk==In_In || outchk==In_In_Horiz) {
		if ( cut==0 && Opt.in_cutoff &&
			BETTER(Opt.in_cutoff, MEASURE(node->child)) ) {
			/* cutoff for ingroup */
			cut = 1;
		} else {
		/* in, in (or in,in+out_horiz) */
			cut = duplicationCheck_inGrp(node);
		/* force to check weak condition */
			if (cut==0 && outchk==In_In_Horiz) {
				cut = 2;
			}
/*
printf("ingrp_root:"); printNode(node); printf(":%d\n",cut);
*/
		}
	} else {
		cut = duplicationCheck(node);
	}
	if ( cut == 2 ) {
		if (Opt.distdiffcut) {
			/* weak cut condition (distdiffcut) */
			if (phyloCut_disttest(node) == 0) {
				cut = 0;
			}
		} else {
			cut = 0;
		}
	}
/*
printf("postCheckOutGrpRoot:cut=%d\n",cut);
*/
	/* try to set out-root at the node that will be split into multiple ingroup clusters
		like  (((in1,in1), (in2, in2)), (in3, in3))
		(i.e. a paralogous node containing only ingroup species)
	*/
	if ( (outchk==In_In || outchk==In_In_Horiz) && cut) {
		cut = 3;
	}

	return cut;
}
checkOutGrpRoot(Node *node, OutGrpStat *ret_outchk) {
	OutGrpStat outchk;
	int cut;
	outchk = outgroupCheckRaw(node);
/*
printf("Cut00:");printNode(node);printf("%d,%d\n", outchk, isOutGrpNode(node));
*/
	if (ret_outchk) {
		*ret_outchk = outchk;
	}
	cut = checkOutGrpCutNode(outchk, node);
	if (cut == 2) {
		/* check children */
		if (node->child) {
			if (isFlankNode1(node)) {
				cut = checkOutGrpRoot(node->child->node1,NULL);
			} else if (isFlankNode2(node)) {
				cut = checkOutGrpRoot(node->child->node2,NULL);
			} else {
				int cut1, cut2;
				cut1 = checkOutGrpRoot(node->child->node1,NULL);
				cut2 = checkOutGrpRoot(node->child->node2,NULL);
				if (cut1 || cut2) {
					cut = 1;
				} else {
					cut = 0;
				}
			}
		} else {
			cut = 0;
		}
	}
	return cut;
}
checkOutGrpRoot_forTreeNode_All(TreeNode *tnode, OutGrpStat *ret_outchk)
{
	int cut = checkOutGrpRoot_forTreeNode(tnode, ret_outchk);
	if (cut == 0) {
		cut = post_checkOutGrpRoot(tnode->node, *ret_outchk);
	}
	return cut;
}
checkOutGrpRoot_forTreeNode(TreeNode *tnode, OutGrpStat *ret_outchk)
{
	OutGrpStat outchk;
OutGrpStat outchk_node;
	int cut;

	if (isLeaf(tnode->node)) {
		return 0;
	}

	outchk = outgroupCheckRaw_TreeNode(tnode);

/*
outchk_node = outgroupCheckRaw(tnode->node);
printTreeNode(tnode);printf("\n");
if (outchk != outchk_node) {
printf("ooout>%d,%d\n",outchk,outchk_node);
} else {
printf("same=ooout>%d,%d\n",outchk,outchk_node);
}
outchk = outchk_node;
*/

	if (ret_outchk) {
		*ret_outchk = outchk;
	}
	cut = checkOutGrpCutNode(outchk, tnode->node);
/*
printf("checkOutGrpRoot:cut=%d,%d\n",cut,outchk);
*/
	if (cut == 2) {
		int cut1 = 0, cut2 = 0;
		/* check children */
		if (tnode->child1) {
			cut1 = checkOutGrpRoot_forTreeNode(tnode->child1,NULL);
		}
		if (tnode->child2) {
			cut2 = checkOutGrpRoot_forTreeNode(tnode->child2,NULL);
		}
		if (cut1 || cut2) {
			cut = 1;
		} else {
			cut = 0;
		}
	}
	return cut;
}
/* basic rule for cutting outgruoup nodes */
checkOutGrpCutNode(OutGrpStat outchk, Node *node)
{
	if (outchk == Out_Out) {
		/* out,out -> not cut */
		return 0;
	} else if (outchk == In_In || outchk == In_In_Horiz) {
/*
	} else if (outchk == In_In) {
*/
		/* in,in -> not cut */
		return 0;
/*
	} else if (outchk == In_In_Horiz) {
		return duplicationCheck_inGrp(node);
*/
	} else if (outchk == Both_Both) {
		/* in+out,in+out -> cut */
		return 1;
	} else {
		/* outchk == Out_Both || Both_Out */
		/* out,in+out || in+out,out -> need to check children */
		return 2;
	}
}
isInGrpNode(OutGrpStat outchk) {
	if (outchk == In_In || outchk == In_In_Horiz) {
		return 1;
	} else {
		return 0;
	}
}
/**
checkOutRoot(Node *node)
{
	if (node->child) {
		if (isFlankNode(node)) {
		} else {
		}
	}
}
**/

#ifdef OOOO
cancelOutRoot(Node *node)
{
	if (isRoot(node)) {
		/* cancel this root node */
/*
printNode(node);
printf(":CANCEL\n");
*/
		deleteRootNode(node);
	} else {
/*
printNode(node);printf(">>parent:%d,%d,%d\n",node->parentL,node->parent,node->parentR);
if (node->parent){
printNode(node->parent);printf(">>cantest\n");
}
*/
		if (nodeDefinedParentL(node)) cancelOutRoot(node->parentL);
		if (nodeDefinedParentC(node)) cancelOutRoot(node->parent);
		if (nodeDefinedParentR(node)) cancelOutRoot(node->parentR);
	}
}
#endif

/* skip nodes that contain only masked species */
checkSkipMaskAll(NodeSet *nodes)
{
	Node *node;
	NodeSetIter iter;

	getRootInit(&iter, nodes, -1, 0);
	while (node = getRootNext(&iter)) {
		if (! nodeDeleted0(node)) {
			checkSkipMask(node);
		}
	}
}
checkSkipMask(Node *node)
{
	Node *child1, *child2;
	int maskchk = spFlagMaskCheck(node->spflag, SPflags.spMask);
#ifdef WITH_DEBUG
/***
printf("Match:%d\n",maskchk);
print_specFlag(node->spflag);
print_specFlag(SPflags.spMask);
printf("===\n");
**/
#endif
	if (maskchk == 2) {
		/** all organisms are matched with specified ones **/
		return 0;
	} else if (maskchk == 0) {
		/** no organism is matched with specified ones **/
		deleteRootNode(node);
		return 0;
	}
#ifdef AA
	if (node->child) {
		child1 = node->child->node1;
		child2 = node->child->node2;
		if (! isFlankNode(node)) {
			if (spFlagMaskCheck(child1->spflag, SPflags.spMask) == 0) {
				/* child1 should be deleted */
				unsetFlagNode(node, NODE_INT1);
			} else if (spFlagMaskCheck(child2->spflag, SPflags.spMask) == 0) {
				/* child2 should be deleted */
				unsetFlagNode(node, NODE_INT2);
			}
		}
	}
#endif
	return 0;
}

/* delete node and recover neighbor relationships */
deleteNode(Node *node)
{
	listIter *iterL, *iterR;
	Neighbor *nbrL, *nbrR;
	Node *nodeL, *nodeR;
	int cnt;
	char status;

	if (nodeDeleted0(node)) {
		return 0;
	}

if (Opt.DEBUG & DBG_nbrrestore) {
	printf("DelFlank: "); printNode(node); printf("\n");
}

	/***
		 +<-------newlink------->+
		 |     nbrL       nbrR   |
		nodeL --+-- node --+-- nodeR
			   delete
	**/

	deleteRootNode(node);
	if (node->left && node->right) {
		iterL = createListIter(node->left, 1);
		iterR = createListIter(node->right, 1);
		while (nbrL = (Neighbor *)getListIter(iterL)) {
			nodeL = (Node *) (((Node *) nbrL->node1 == node)
					? nbrL->node2 : nbrL->node1);

			setListIter(iterR, node->right, 1);

			while (nbrR = (Neighbor *)getListIter(iterR)) {
				nodeR = (Node *) (
					    ((Node *) nbrR->node1 == node) ?
					    nbrR->node2 : nbrR->node1);

				cnt = (nbrL->count < nbrR->count) ?
					nbrL->count : nbrR->count;
			/***
				   nbrL       nbrR
				L -{1}- curr -{-1}- R
				1         1        -1
			**/
				if (nbrL->status == 0 && nbrR->status == 0){
					status = 0;
				} else if (nbrL->status == 2 || nbrR->status == 2){
					status = 2;
				} else {
					status = 1;
				}

if (Opt.DEBUG & DBG_nbrrestore) {
	printf("newnode: ");
	printNode(nodeL); printNode(nodeR);
	printf(" %d\n", status);
}
				addNeighborRestore(nodeL,nodeR,
					nbrL->dir,nbrR->dir,cnt, status);
			}
		}
	}
	return 1;
}
pList *createRootNodeList(NodeSet *nodes)
{
	NodeSetIter iter;
	Node *node;
	pList *rootNodes = create_pList();

	getRootInit(&iter, nodes, -1, 1);
	while (node = getRootNext(&iter)) {
		pushList(rootNodes, node);
	}
	return rootNodes;
}

deleteOutgroupClusters(NodeSet *nodes)
{
	NodeSetIter iter;
	Node *node;

	getRootInit(&iter, nodes, -1, 0);
	while (node = getRootNext(&iter)) {
		if (outgroupFlagCheck(node->spflag) == 2) {
		/* the set contains only the outgroup organisms */
			deleteRootNode(node);
		}
	}
}

typedef struct {
	int (*getChilds)();
	specFlagP (*getSpFlag)();
	int (*isLeaf)();
} NodeFunc;

int getChilds(), isLeaf();
int getChilds_TreeNode(), ClustTree_isLeaf();
specFlagP getSpFlag(), getSpFlag_TreeNode();

OutGrpStat _outgroupCheckRaw();
int _countGenesWithSPflags();
NodeFunc nodeFunc = {
	getChilds, getSpFlag, isLeaf
};

NodeFunc treeNodeFunc = {
	getChilds_TreeNode, getSpFlag_TreeNode, ClustTree_isLeaf
};

/* Classifying out/in-group patterns of both subtrees */
/* Do not discriminate In_In and In_In_Horiz */
OutGrpStat outgroupCheck(Node *node)
{
	OutGrpStat outchk = outgroupCheckRaw(node);
	if (outchk == In_In_Horiz) {
		return In_In;
	} else {
		return outchk;
	}
}
OutGrpStat outgroupCheck_TreeNode(TreeNode *tnode)
{
	OutGrpStat outchk = outgroupCheckRaw_TreeNode(tnode);
	if (outchk == In_In_Horiz) {
		return In_In;
	} else {
		return outchk;
	}
}

OutGrpStat outgroupCheckRaw(Node *node)
{
	return _outgroupCheckRaw((char*) node, &nodeFunc);
}
OutGrpStat outgroupCheckRaw_TreeNode(TreeNode *tnode)
{
	return _outgroupCheckRaw((char*) tnode, &treeNodeFunc);
}

OutGrpStat _outgroupCheckRaw(char *node, NodeFunc *nfunc)
{
	specFlagP spflag1, spflag2;
	int check1, check2;

	char *child1, *child2;

	specFlagP spflag = nfunc->getSpFlag(node);

	if (nfunc->isLeaf(node)) {
		check1 = outgroupFlagCheck(spflag);
		if (check1 == 0) {
			return In_In;
		} else {
			return Out_Out;
		}
	}
	nfunc->getChilds(node, &child1, &child2);

	if (! child1) {
		return _outgroupCheckRaw(child2, nfunc);
	} else if (! child2) {
		return _outgroupCheckRaw(child1, nfunc);
	}

	spflag1 = nfunc->getSpFlag(child1);
	spflag2 = nfunc->getSpFlag(child2);
	check1 = outgroupFlagCheck(spflag1);
	check2 = outgroupFlagCheck(spflag2);
/*
	printf("check=%d,%d\n",check1,check2);
*/

	if (check1 == 0 && check2 == 0) {
		/* in, in : check */
		return In_In;
	} else if (check1 == 2 && check2 == 2) {
		/* out, out : ignore */
		return Out_Out;
	} else if (check1 == 2) {
		/* out, in+out : add left, check right */
		return Out_Both;
	} else if (check2 == 2) {
		/* in+out, out : add right, check left */
		return Both_Out;
	} else {
		/* in, in+out or in+out, in+out : cut */

/*
printf("spflag===\n");
print_specFlag(spflag1);
print_specFlag(spflag2);
*/
		if (Opt.horizweight > 0) {
			double outcnt1,outcnt2,outcntM, incnt1,incnt2,incntM;
			double in_del, out_del, horiz_O2I, horiz_I2O, dupLoss;
			double delcnt1, delcnt2;

			specFlag spflagM; /* merged set */
			specFlag spflagM_outgrp, spflagM_ingrp;

			spFlagOR(spflag1,spflag2,spflagM);

			outcnt1 = spFlagANDcntW(spflag1, SPflags.outGroup);
			incnt1 = spFlagANDcntW(spflag1, SPflags.inGroup);
			outcnt2 = spFlagANDcntW(spflag2, SPflags.outGroup);
			incnt2 = spFlagANDcntW(spflag2, SPflags.inGroup);
			outcntM = spFlagANDcntW(spflagM, SPflags.outGroup);
			incntM = spFlagANDcntW(spflagM, SPflags.inGroup);
			in_del = (incntM - incnt1) + (incntM - incnt2);
			out_del = (outcntM - outcnt1) + (outcntM - outcnt2);
/*
			horiz = outcntM;
*/

			spFlagAND(spflagM, SPflags.outGroup, spflagM_outgrp);
			spFlagAND(spflagM, SPflags.inGroup, spflagM_ingrp);
	
			/* HGT from ingroup to outgroup species; count num_spec for each ougroup cluster (count_one=0) */
			horiz_I2O = _countGenesWithSPflags(node, spflagM_outgrp, nfunc, 0);
			/* HGT from outgroup to ingroup species; count 1 for each ingroup cluster (count_one=1) */
			horiz_O2I = _countGenesWithSPflags(node, spflagM_ingrp, nfunc, 1);
			dupLoss = (in_del + out_del) * Opt.horizweight;

			if ((in_del + out_del) < horiz_I2O / Opt.horizweight) {
				return Both_Both;
			} else {
				return In_In_Horiz;
			}
		} else if (Opt.horizweight < 0) {
			return In_In_Horiz;
		}
		return Both_Both;
		
	}
}

_countGenesWithSPflags(char *node, specFlag spflag, NodeFunc *nfunc, int count_one)
{
	int count = 0;
	specFlagP node_spflag = nfunc->getSpFlag(node);
	char *child1, *child2;

	if (spFlagInclude(spflag, node_spflag)) {
		if (count_one) {
			return 1;
		} else {
			return spFlagCnt(node_spflag);
		}
	} else if (! nfunc->isLeaf(node)) {
		nfunc->getChilds(node, &child1, &child2);
		if (child1) {
			count += _countGenesWithSPflags(child1, spflag, nfunc, count_one);
		}
		if (child2) {
			count += _countGenesWithSPflags(child2, spflag, nfunc, count_one);
		}
	}
	return count;
}
countGenesWithSPflags(Node *node, specFlag spflag)
{
	int count = 0;
	if (spFlagInclude(spflag, node->spflag)) {
		return spFlagCnt(node->spflag);
	} else if (node->child) {
		if (isFlankNode1(node)) {
			count = countGenesWithSPflags(node->child->node1, spflag);
		} else if (isFlankNode2(node)) {
			count = countGenesWithSPflags(node->child->node2, spflag);
		} else {
			count = countGenesWithSPflags(node->child->node1, spflag);
			count += countGenesWithSPflags(node->child->node2, spflag);
		}
	}
	return count;
}

#ifdef AAA
/*	0 .. no outgroup species,
	1 .. there is a species belonging to the outgroup,
	2 .. all species are belonging to the outgroup */
outgroupFlagCheck(specFlag spflag)
{
	return spFlagMaskCheck(spflag, SPflags.outGroup);
}
/*	0 .. no ingroup species,
	1 .. there is a species belonging to the ingroup,
	2 .. all species are belonging to the ingroup */
ingroupFlagCheck(specFlag spflag)
{
	return spFlagMaskCheck(spflag, SPflags.inGroup);
}
#endif

#ifdef AAA
spFlagMaskCheck(specFlag spflag, specFlag spFlagPat)
{
	specFlag matched;
	int matched_spcnt, spcnt;

	spFlagAND(spflag, spFlagPat, matched);
	matched_spcnt = spFlagCnt(matched);
	if (matched_spcnt) {
		spcnt = spFlagCnt(spflag);
		if (matched_spcnt == spcnt) {
		/* all organisms in the set are matched */
			return 2;
		} else {
		/* there exists at least one matched organism */
			return 1;
		}
	}
	return 0;
}
#endif

int duplicationCheck_All(Node *node)
{
	int cut;
	cut = duplicationCheck(node);
/*
printf("checkAll:"); printNode(node);printf("%d<<f\n",cut);
*/
	if (cut == 1) {
		/** immediately cut this node **/
	} else if ( cut == 2 ) {
#ifdef EXPERIMENTAL
		if ((Opt.chkConnect && ! checkConnect(node)) ||
			(Opt.sumcut &&
			Opt.sumcut > getEdgeScoreSum(node->child))) {

			/* weak cut condition (connect/sumcut) */
		} else
#endif		  
		if (Opt.distdiffcut) {
			/* weak cut condition (distdiffcut) */
			int cut1 = phyloCut_disttest(node);
			if (cut1 == 0) {
				cut = 0;
			} else {
				/* cut = 2 */
			}
		} else {
			cut = 0;
		}
	}
	return cut;
}

/* check the phylogenetic tree cutting criterion */
duplicationCheck(Node *node)
{
	return duplicationCheck0(node->spflag, node->child->node1->spflag,
		node->child->node2->spflag);
}
/* check duplications only among ingroup species */
duplicationCheck_inGrp(Node *node)
{
	duplicationCheck_inGrp_flag(node->spflag, node->child->node1->spflag, node->child->node2->spflag);
}
duplicationCheck_inGrp_TreeNode(TreeNode *tnode)
{
	duplicationCheck_inGrp_flag(tnode->spflag, tnode->child1->spflag, tnode->child2->spflag);
}
duplicationCheck_inGrp_flag(specFlag spflag, specFlag spflag_ch1, specFlag spflag_ch2)
{
	specFlag spflag0, spflag1, spflag2;
	spFlagAND(spflag, SPflags.inGroup, spflag0);
	spFlagAND(spflag_ch1, SPflags.inGroup, spflag1);
	spFlagAND(spflag_ch2, SPflags.inGroup, spflag2);
	return duplicationCheck0(spflag0, spflag1, spflag2);
}
duplicationCheck0(specFlag spflag0, specFlag spflag1, specFlag spflag2)
{
	int mchcnt, spcnt1, spcnt2, mchcnt2;
	double wt_mchcnt, wt_spcnt1, wt_spcnt2;

	/** Either one is an unknown-only cluster -- cut **/
/**
	if (unknownOnlyCluster(spflag1) || unknownOnlyCluster(spflag2)) {
		return 1;
	}
**/

	mchcnt = spFlagANDcnt(spflag1, spflag2);
	spcnt1 = spFlagCnt(spflag1);
	spcnt2 = spFlagCnt(spflag2);

	if (mchcnt * spcnt1 * spcnt2 == 1) {
		/** Do not cut this node **/
		return 0;
	}

	wt_mchcnt = spFlagANDcntW(spflag1, spflag2);
	wt_spcnt1 = spFlagCntW(spflag1);
	wt_spcnt2 = spFlagCntW(spflag2);
/*
printf("mchcnt,spcnt1,spcnt2>>%d,%d,%d\n",mchcnt,spcnt1,spcnt2);
printf("wt_mchcnt,wt_spcnt1,wt_spcnt2>>%lf,%lf,%lf\n",wt_mchcnt,wt_spcnt1,wt_spcnt2);
printf(">>>%lf,%lf\n",wt_spcnt1*Opt.phylocutratio,wt_spcnt2*Opt.phylocutratio);
*/


	/** cut node when two subgroups share organisms with a ratio of
		more than Opt.phylocutratio: do not cut when the
		ratio is just Opt.phylocutratio */
	if ((wt_mchcnt > (double) wt_spcnt1 * Opt.phylocutratio
		 || wt_mchcnt > (double) wt_spcnt2 * Opt.phylocutratio)
			 && (mchcnt * spcnt1 * spcnt2 != 1) ) {
		/** Immediately cut this node **/
		return 1;
	}
	if (Opt.treecheck) {
		/* Mapping the nodes of gene tree onto the species tree */
		int spnode0, spnode1, spnode2;
		spnode0 = sptree_MatchFlags(spflag0);
		spnode1 = sptree_MatchFlags(spflag1);
		spnode2 = sptree_MatchFlags(spflag2);
		if (spnode0 == spnode1 || spnode0 == spnode2) {
			/* duplication !! */
#ifdef WITH_DEBUG
/*
	printf("==\n");
	print_specFlag(spflag0);
	print_specFlag(spflag1);
	print_specFlag(spflag2);
	printf("%d,%d,%d\n",spnode0,spnode1,spnode2);
*/
#endif
			/* sum of the weights of the duplicated species */
			wt_mchcnt =
			    sptree_MatchFlagsCntW(spflag1,spflag2,spnode0);

			if ((wt_mchcnt > (double) wt_spcnt1 * Opt.phylocutratio
			 || wt_mchcnt > (double) wt_spcnt2 * Opt.phylocutratio)
		 		&& (mchcnt * spcnt1 * spcnt2 != 1) ) {
				/** Cut this node **/
				return 1;
			}
		}
	}
	if (! Opt.phylocutratio2 && ! Opt.treeBalanceRatio) {
		return 2;
	}
	if ( Opt.phylocutratio2 &&
		(wt_mchcnt > (double) wt_spcnt1 * Opt.phylocutratio2
		 || wt_mchcnt > (double) wt_spcnt2 * Opt.phylocutratio2) ) {
		return 2;
	}
	if ( Opt.treeBalanceRatio &&
		(wt_spcnt1 < wt_spcnt2 * Opt.treeBalanceRatio || wt_spcnt2 < wt_spcnt1 * Opt.treeBalanceRatio) ) {
		return 2;
	}
	/** Do not cut this node **/
	return 0;
}

restoreBreak(Node *node, Node *parent)
{
	if (! parent) return(0);
/*
printf("RESTORE:");printNode(node);printf(" ");printNode(parent);printf(":%d\n");
*/
	if (parent == node->parentL) {
/*
printf("delL>>");printNode(node->parentL);printf(" ");printf("%d,%d\n",isOuterRoot(node->parentL),isDeletedOuterRoot(node->parentL));printf("\n");;
*/
		node->brk.from = INFPOS;
		if (isOuterRoot(node->parentL)) {
			/* reset temporary outerRootFlag at the flanking node */
			setDeletedOuterRoot(node->parentL);
		}
/*
		node->parentL = NULL;
*/
	} else if (parent == node->parentR) {
/*
printf("delR>>");printNode(node->parentR);printf(" ");printf("%d,%d\n",isOuterRoot(node->parentR),isDeletedOuterRoot(node->parentR));printf("\n");;
*/
		if (node->parentM) {
			node->brk2.to = SUPPOS;
		} else {
			node->brk.to = SUPPOS;
		}
		if (isOuterRoot(node->parentR)) {
			/* reset temporary outerRootFlag at the flanking node */
			setDeletedOuterRoot(node->parentR);
		}
/*
		node->parentR = NULL;
*/
	} else if (parent == node->parentM) {
		node->brk.to = SUPPOS;
		node->brk2.from = INFPOS;
/*
		node->parentM = NULL;
*/
	} else if (parent == node->parent) {
		if(definedPos(node->brk.from)) {
/*
			node->brk.from = INFPOS;
*/
			node->brk.from = SUPPOS;
		} else if (definedPos(node->brk.to)){
			node->brk.to = INFPOS;
/*
			node->brk.to = SUPPOS;
*/
		}
		if (parent->child->node1 == parent->child->node2) {
			if(definedPos(node->brk2.from)) {
				node->brk2.from = INFPOS;
			} else if(definedPos(node->brk2.to)) {
				node->brk2.to = SUPPOS;
			}
		}
	/* the centeral outerRoot node is a valid orholog node; should not be deleted */
/*
		if (isOuterRoot(node->parent)) {
			setDeletedOuterRoot(node->parent);
		}
*/
/*
		if (! nodeMerged(node->parent)) {
			node->parent = NULL;
		}
*/
	} else {
		/*** error **/
		if (nodeDefinedParentC(node) && nodeDefinedParentR(node) && nodeDefinedParentL(node)  && node->parentM!=NULL) {
			fprintf(stderr, ">>ERROR %d: %d,%d,%d\n",
				parent->id,node->parent!=NULL,node->parentL!=NULL,node->parentR!=NULL);
#ifdef WITH_DEBUG
			printf("??? %s,%d %d,%d,%d,%d,%d\n",node->name,node->id,node->parent,node->parentL,node->parentR,parent,parent->id);
			printf("??? %d %d\n", parent->child->node1, parent->child->node2);
#endif
			printNode(node); putchar(' ');
			if (nodeDefinedParentC(node)){
				printNode(node->parent); putchar('\n');
			}
			if (nodeDefinedParentL(node)){
				putchar('L');
				printNode(node->parentL); putchar('\n');
			}
			if (nodeDefinedParentL(node)){
				putchar('R');
				printNode(node->parentR); putchar('\n');
			}
		}
	}
}
checkUnvisitedParent(Node *node)
{
	if ( (! nodeDefinedParentC(node) || nodeVisited(node->parent)) &&
		(! nodeDefinedParentL(node) || nodeVisited(node->parentL)) &&
		(! nodeDefinedParentR(node) || nodeVisited(node->parentR)) ) {

		/** there remains no unvisited parent **/
		return 1;
	}
/*
if (node->parent != NULL && ! nodeVisited(node->parent)) {printNode(node->parent);printf(":M\n");}
if (node->parentL != NULL && ! nodeVisited(node->parentL)) {printNode(node->parentL);printf(":L\n");}
if (node->parentR != NULL && ! nodeVisited(node->parentR)) {printNode(node->parentR);printf(":R\n");}
printNode(node); printf(": unvisited=yes\n");
*/
	/** there remains unvisited parent **/
	return 0;
}
check1nvisitedCenterParent(Node *node)
{
	if (! nodeDefinedParentC(node) || nodeVisited(node->parent))  {
		/* visited */
		return 1;
	} else {
		return 0;
	}
}
checkDeletedOuterParent(Node *node)
{

/*
printf("checkDelOuterParent>>");
if (node->parent) {printf("M:");printNode(node->parent);}
if (node->parentL) {printf("L:");printNode(node->parentL);}
if (node->parentR) {printf("R:");printNode(node->parentR);}
printf("\n");
*/
	/* check existence of parents whether or not they are deleted during defining clustering roots */
/*
	if ( (! nodeDefinedParentC(node) || nodeOuterVisited(node->parent)) &&
		(! nodeDefinedParentL(node) || nodeOuterVisited(node->parentL)) &&
		(! nodeDefinedParentR(node)|| nodeOuterVisited(node->parentR)) ) {
*/
	if ( (! nodeExistParentC(node) || nodeOuterVisited(node->parent)) &&
		(! nodeExistParentL(node) || nodeOuterVisited(node->parentL)) &&
		(! nodeExistParentR(node)|| nodeOuterVisited(node->parentR)) ) {

		/** there remains no unvisited parent **/
		return 1;
	}
}

checkClusterSize(NodeSet *nodes)
{
	NodeSetIter iter;
	Node *node;
	double spcnt;
	int i;

	getRootInit(&iter, nodes, -1, 0);
	while (node = getRootNext(&iter)) {
		if (node->child) {
			spcnt = spFlagCntW_All(node->spflag);

			if (spcnt < Opt.minsp || node->cnt < Opt.minent) {
				deleteNode(node);
			}
		}
	}
}

/***
 addRootNode: find additional root nodes by checking trees (represented in TreeNode) collected in clustRootsHom
***/

static TreeNode *curr_outerroot;
addRootNode(pList *clustRootsHom, pList *clustRoots)
{
	listIter iter;
	SubClusterInfo *homclustInfo;

	setListIter(&iter, clustRootsHom, 1);
	while (homclustInfo = (SubClusterInfo *) getListIter(&iter)) {
		curr_outerroot = homclustInfo->root;
		addRootNode_sub(homclustInfo->root, clustRoots, homclustInfo);
	}
}
addRootNode_sub(TreeNode *tnode, pList *clustRoots, SubClusterInfo *homclustInfo)
{
	OutGrpStat outchk = None;
	int find_root = 0;
	int find_root_ch1 = 0, find_root_ch2 = 0;
	int cut;
/*
printf("add_RootNode:"); printNode(tnode->node); printf(":%d,%d\n",isOuterRoot(tnode->node),isRoot(tnode->node));
*/

	if (! tnode) {
		return 0;
	}

	if (isTreeRoot(tnode)) {
		find_root = 1;
		/* obtain outchk */
/*
	} else if (isOuterRoot(tnode->node)) {
		addRootNode_sub(tnode->child1);
		addRootNode_sub(tnode->child2);
*/
	} else if (ClustTree_isLeaf(tnode)) {
		if (outgroupFlagCheck(tnode->node->spflag) == 0) {
				makeRootTreeNode_Out(tnode);
				pushList(clustRoots, tnode);
				pushList(homclustInfo->members, tnode);
				return 1;
/*
printf("NEWROOT_LEAF:"); printNode(tnode->node); printf(":%d,%d\n",isOuterTreeRoot(tnode),isTreeRoot(tnode));
*/
		} 
	} else {
		/* no root is found among descendants */
		cut = checkOutGrpRoot_forTreeNode_All(tnode, &outchk);
		if (cut == 0) {
			/* check if children have a root rootnode */ 
			find_root_ch1 = addRootNode_check(tnode->child1, clustRoots, homclustInfo);
			find_root_ch2 = addRootNode_check(tnode->child2, clustRoots, homclustInfo);
/*
printf("cut=%d,find_ch1=%d,find_ch2=%d\n",cut,find_root_ch1, find_root_ch2);
*/
			if (! find_root_ch1 && ! find_root_ch2) {
				/* new root node */
				if (outgroupFlagCheck(tnode->node->spflag) != 2) {
					/* not outgroup node */
/*
printf("NEWROOT:"); printNode(tnode->node); printf(":%d:%d,%d\n",tnode,isOuterTreeRoot(tnode),isTreeRoot(tnode));
*/
					makeRootTreeNode_Out(tnode);
					pushList(clustRoots, tnode);
					pushList(homclustInfo->members, tnode);
					resetTreeNodeDomain(tnode, curr_outerroot);
					/* check to find ingroup root */
/*
					if (outchk == In_In || outchk == In_In_Horiz)  {
						find_root_ch1 = addRootNode_In(tnode->child1, clustRoots, homclustInfo);
						find_root_ch2 = addRootNode_In(tnode->child2, clustRoots, homclustInfo);
					}
*/
				}
				find_root = 1;
			}
		} else {
			find_root_ch1 = addRootNode_sub(tnode->child1, clustRoots, homclustInfo);
			find_root_ch2 = addRootNode_sub(tnode->child2, clustRoots, homclustInfo);
		}
	}
	if (find_root) {
/*
		addRootNode_In(tnode, clustRoots, homclustInfo);
printf("findroot:");printTreeNode(tnode);printf("\n");
*/
		if (! ClustTree_isLeaf(tnode)) {
			if (outchk == None) {
				checkOutGrpRoot_forTreeNode_All(tnode, &outchk);
			}
			if (outchk == Both_Out) {
				addRootNode_In(tnode->child1, clustRoots, homclustInfo);
			} else if (outchk == Out_Both) {
				addRootNode_In(tnode->child2, clustRoots, homclustInfo);
			}
		}
		return 1;
	}
	return (find_root_ch1 || find_root_ch2);
}
addRootNode_check(TreeNode *tnode, pList *clustRoots, SubClusterInfo *homclustInfo)
{
	int find_root_ch1 = 0, find_root_ch2 = 0;

	if (! tnode) {
	} else if (tnode->treenodeFlag & TREENODE_HASROOT) { 
		return 1;
	} else if (isTreeRoot(tnode)) {
		return 1;
	} else if (! ClustTree_isLeaf(tnode)) {
		find_root_ch1 = addRootNode_check(tnode->child1, clustRoots, homclustInfo);
		find_root_ch2 = addRootNode_check(tnode->child2, clustRoots, homclustInfo);
	}
	if (find_root_ch1 || find_root_ch2) {
		tnode->treenodeFlag |= TREENODE_HASROOT;
		return 1;
	}
	return 0;
}
addRootNode_In(TreeNode *tnode, pList *clustRoots, SubClusterInfo *homclustInfo)
{
	int cut;
	OutGrpStat outchk = None;
	int find_root_ch1 = 0, find_root_ch2 = 0;
/*
if (tnode) {
printf("check_add_in_node:"); printTreeNode(tnode); printf(":%d,%d\n",isOuterRoot(tnode->node),isRoot(tnode->node));
}
*/
	if (! tnode) {
		return 0;
/*
	} else if (isInTreeRoot(tnode)) {
printf("already root");printTreeNode(tnode);printf("\n");
		return 1;
*/
	} else if (! ClustTree_isLeaf(tnode)) {
		cut = checkOutGrpRoot_forTreeNode_All(tnode, &outchk);
		if (cut==0) {
			cut = post_checkOutGrpRoot(tnode->node, outchk);
			if (cut == 0 && (outchk == Both_Out || outchk == Out_Both)) {
				cut = 3;
			}
		}
		if (cut==0) {
			/* find inroot among children */
			if (tnode->child1){
				find_root_ch1 = addRootNode_In_check(tnode->child1, clustRoots, homclustInfo);
			}
			if (tnode->child2){
				find_root_ch2 = addRootNode_In_check(tnode->child2, clustRoots, homclustInfo);
			}
			if (! find_root_ch1 && ! find_root_ch2) {
				/* add root if no root is found among children */
				makeRootTreeNodeIn(tnode);
/*
printf("NEW_INROOT:"); printNode(tnode->node); printf(":%d:%d,%d,%d\n",tnode,isOuterTreeRoot(tnode),isTreeRoot(tnode),isAllInTreeRoot(tnode));
*/
			}
		} else {
			addRootNode_In(tnode->child1, clustRoots, homclustInfo);
			addRootNode_In(tnode->child2, clustRoots, homclustInfo);
		}
	}
/*
printf("done:");printTreeNode(tnode);printf("\n");
*/
	return (find_root_ch1 || find_root_ch2);
}
addRootNode_In_check(TreeNode *tnode, pList *clustRoots, SubClusterInfo *homclustInfo)
{
	int find_root_ch1 = 0, find_root_ch2 = 0;
	if (isAllTreeRoot(tnode)) {
		return 1;
	} else if (! ClustTree_isLeaf(tnode)) {
		find_root_ch1 = addRootNode_In_check(tnode->child1, clustRoots, homclustInfo);
		find_root_ch2 = addRootNode_In_check(tnode->child2, clustRoots, homclustInfo);
	}
	return (find_root_ch1 || find_root_ch2);
}

resetTreeNodeDomain(TreeNode *tnode, TreeNode *rootnode)
{
	Domain *dom;
	if (tnode == NULL) {
	} else if (ClustTree_isLeaf(tnode)) {
		getDomainNoMark(tnode->node->domains, rootnode->node, &dom);
		if (dom) {
/*
printf("resetTreeNodeDomain:");
printTreeNode(tnode);
printDomain(dom);
printf("\n");
*/
			tnode->dom = dom;
		}
	} else {
		resetTreeNodeDomain(tnode->child1, rootnode);
		resetTreeNodeDomain(tnode->child2, rootnode);
	}
}

/***
 addOuterGroup: find root node descendant from a given outer-root and assign outer-root each root node
***/

/* for (new) outgroup mode */
typedef struct {
	TreeNode *root;
	Hash *hash;
	pList *outergroup_mem;
	pList *outerRoots;
	SubClusterInfo *outerInfo;
} OuterGroupInfo;

addOuterGroup(pList *clustRootsHom, Hash *cinfoHash, pList **outerRoots)
{
	listIter iter;
	SubClusterInfo *hom_cInfo;
	TreeNode *root;
	static OuterGroupInfo outerGroupInfo;

	outerGroupInfo.hash = cinfoHash;
	outerGroupInfo.outergroup_mem = create_pList();
	*outerRoots = outerGroupInfo.outerRoots = create_pList();

	setListIter(&iter, clustRootsHom, 1);
/*
	while (hom_cInfo = (SubClusterInfo*) getListIter(&iter)) {
		root = hom_cInfo->root;
		if (root->node->id == 92988) {
			printf(">>>>92988"); printNode(root->node); printf("; %d,%d\n", root->node->child->node1, root->node->child->node2); }
	}
*/
	setListIter(&iter, clustRootsHom, 1);
	while (hom_cInfo = (SubClusterInfo*) getListIter(&iter)) {
		root = hom_cInfo->root;
/*
printf("addOuter:");printTreeNode(root);printf("\n");
*/
		clearList(outerGroupInfo.outergroup_mem);
		outerGroupInfo.root = NULL;
		addOuterGroup_sub(root, &outerGroupInfo);
	}
	free_pList(outerGroupInfo.outergroup_mem);
}
addOuterGroup_sub(TreeNode *tnode, OuterGroupInfo *ogrpInfo)
{
	int setRoot_flag = 0;
	int check_flag = 0;
	int push_flag = 0;
	OutGrpStat outchk;
	ClusterInfo *cInfo;
	listIter iter;
	TreeNode *n;
	pList *ogrpList;
	SubClusterInfo *outerInfo;
int i;

/*
if (tnode) {
  Node *node = tnode->node;
  if (strcmp(node->name, "mth:MTH802")==0 || strcmp(node->name, "pab:PAB1674")==0){
	printf(">>##>>");printNode(node);printf("\n");
  }
}
printf("###");printTreeNode(tnode);printf(":%d\n",isOuterTreeRoot(tnode));
*/
	if (! tnode) {
		/** ??? */
		return(0);
	}
/*
fprintf(stderr, "tnode:%d, %d, %d,%d,%d\n",tnode, tnode->node->id, isOuterTreeRoot(tnode), isDeletedOuterRoot_flag(tnode->flag), isTreeRoot(tnode));
*/
	if (isOuterTreeRoot(tnode) || (checkTreeRoot(tnode) && ogrpInfo->root==NULL)) {
		/* a root node of the common outgroup orthoogous group including multiple ortholog group */
		/* ( (I1,I2,O1), (I3, O2), O4 )  */
/*
	printf("AddOuter:OuterRoot: ");printTreeNode(tnode);printf(":%d\n",tnode->flag);
*/
		outerInfo = createOuterClusterInfo(tnode);
		pushList(ogrpInfo->outerRoots, outerInfo);
		ogrpInfo->outerInfo = outerInfo;
		if (ogrpInfo->root == NULL) {
			ogrpInfo->root = tnode;
			setRoot_flag = 1;
		} else {
			/* should not come here (outerRoot should not be nested) */
			printf("outer root is already defined: ");
			printNode(ogrpInfo->root->node);
			printf(" ");
			printNode(tnode->node);
			printf("\n");
		}
	}

/*
	if (isTreeRoot(tnode)) {
*/
	if (checkTreeRoot(tnode)) {

/*
	printf("AddOuter:Root>"); printTreeNode(tnode); printf(":%d,%d,%d\n",tnode->node->id, cInfo,tnode);
*/
		/* unique root node for outgroup mode */
/* do not add root node to the outerRoots list
		outerInfo = createSubClusterInfo(tnode);
		pushList(ogrpInfo->outerRoots, outerInfo);
*/

		if (cinfoHsearch(ogrpInfo->hash, tnode->node->id, &cInfo) == 0) {
			printf("Node not found: "); printTreeNode(tnode); printf("\n");
			printf("STATUS: %d,%d,%d\n",
				isTreeRootRaw(tnode), isCenterRoot_flag(tnode->flag), (tnode->parent && (tnode->node->parent == tnode->parent->node)));
		}

/*
printf("CINFOSERCH>>>>");printTreeNode(tnode);printf("\n");
fflush(stdout);
printf("Found>>");printTreeNode(cInfo->root);printf("\n");
fflush(stdout);
*/

		setListIter(&iter, ogrpInfo->outergroup_mem, 1);
		ogrpList = create_pList();

		while (n = (TreeNode *) getListIter(&iter)) {
			/* collect all domains contained in each internal node in outergroup_mem into ogrpList */
/*
if (n) {
printf("NNNN>node>>"); printNode(n->node);
	if (n->dom) {
		printf(" : "); printDomain(n->dom);
	}
printf("\n");
}
*/
			ClustTree_alldom(n, ogrpList);
		}

/*
		if (cInfo->outgroup == NULL) {
			cInfo->outgroup = createSubClusterInfo(NULL);
		}
*/
		if (cInfo->outer_root == NULL) {
			if (ogrpInfo->root) {
				cInfo->outer_root = ogrpInfo->root;
			} else {
				cInfo->outer_root = tnode;
			}
			/* add this rootnode to the current outerClusterInfo */
			pushList(ogrpInfo->outerInfo->members, tnode);
/*
fprintf(stderr, "out>>>%d\n",cInfo->outer_root);
printf("cinfo=>%d, outer_root=>%d;",cInfo, cInfo->outer_root);printTreeNode(cInfo->outer_root);printf("\n");
			if (cInfo->clusters_outer == NULL) {
				cInfo->clusters_outer = create_pList();
			}
			pushList(cInfo->clusters_outer, cInfo->outer_root);
*/
		} else {
			if (! (cInfo->outer_root == ogrpInfo->root || cInfo->outer_root == tnode)) {
				if ((ogrpInfo->root && cInfo->outer_root->node != ogrpInfo->root->node) ||
						(tnode && cInfo->outer_root->node != tnode->node)) {
					printf("null?????????????????");
					printTreeNode(cInfo->outer_root); printf(" != ");
					if (ogrpInfo->root) {
						printTreeNode(ogrpInfo->root);
					} else if (tnode) {
						printTreeNode(tnode);
					} else {
						printf("null");
					}
					printf("\n");
				}
			}
			if (setRoot_flag && ogrpInfo->root) {
				ogrpInfo->root = NULL;
			}
			return(0);
		}

		/* additional outgroup genes */
		if (cInfo->outgroup2 == NULL) {
			cInfo->outgroup2 = createSubClusterInfo(NULL);
		}

		/* adding outgroup members */
		pushListAll(cInfo->outgroup2->members, ogrpList);
/*
if (numelemList(cInfo->outgroup2->members)>0) {
printf("DOM:\n");
printList(ogrpList,printDomain);
printf("\n");
}
*/

		assignSubgroups_sub(NULL, cInfo->outgroup2, NULL, NULL, NULL);
/*
		setListIter(&iter, ogrpInfo->outergroup_mem, 1);
		while (n = (TreeNode *) getListIter(&iter)) {
			assignSubgroups_sub(NULL, cInfo->outgroup2, NULL, NULL, NULL);
		}
*/

		free_pList(ogrpList);

	} else if (ClustTree_isLeaf(tnode)) {
	} else if (ClustTree_isIgnoredLeaf(tnode)) {
	} else {
#ifdef AAA
		if (isOuterTreeRoot(tnode)) {
			/* a root node of the common outgroup orthoogous group including multiple ortholog group */
			/* ( (I1,I2,O1), (I3, O2), O4 )  */
/*
	printf("AddOuter:OuterRoot: ");printTreeNode(tnode);printf("\n");
*/
			outerInfo = createOuterClusterInfo(tnode);
			pushList(ogrpInfo->outerRoots, outerInfo);
			ogrpInfo->outerInfo = outerInfo;
			if (ogrpInfo->root == NULL) {
				ogrpInfo->root = tnode;
				setRoot_flag = 1;
			} else {
				/* should not come here (outerRoot should not be nested) */
				printf("outer root is already defined: ");
				printNode(ogrpInfo->root->node);
				printf(" ");
				printNode(tnode->node);
				printf("\n");
			}
		}
#endif
		if (ogrpInfo->root) {
			outchk = outgroupCheck_TreeNode(tnode);
/*
printf("spec=");print_specFlag(tnode->node->spflag);printf("\n");
printf("outchk>>%d\n",outchk);
*/

			if (outchk == Out_Both && tnode->child1 && tnode->child2) {
				pushList(ogrpInfo->outergroup_mem, tnode->child1);
				push_flag = 1;
				check_flag = 2;
			} else if (outchk == Both_Out && tnode->child1 && tnode->child2) {
				pushList(ogrpInfo->outergroup_mem, tnode->child2);
				push_flag = 1;
				check_flag = 1;
			} else if ((outchk == Both_Both)  ||
					(outchk == In_In || outchk == In_In_Horiz)) {
				if (tnode->child1) {
					check_flag += 1;
				}
				if (tnode->child2) {
					check_flag += 2;
				}
			}
		} else {
			if (tnode->child1) {
				check_flag += 1;
			}
			if (tnode->child2) {
				check_flag += 2;
			}
		}
/*
if (numelemList(ogrpInfo->outergroup_mem) > 0) {
	printf("OUTERgroup=");
	printList(ogrpInfo->outergroup_mem, printTreeNode);
	printf("\n");
}
*/
		if (check_flag & 1) {
			addOuterGroup_sub(tnode->child1, ogrpInfo);
		}
		if (check_flag & 2) {
			addOuterGroup_sub(tnode->child2, ogrpInfo);
		}
		if (push_flag) {
			popList(ogrpInfo->outergroup_mem);
		}
/*
fprintf(stderr, "00>>>%d,%d,%d\n",ogrpInfo->root,setRoot_flag,ogrpInfo);
if (tnode->node) {
	fprintf(stderr, "NID>>%d\n",tnode->node->id);
}
*/
	}
	if (setRoot_flag && ogrpInfo->root) {
/*
for (i = 0; i < 10; i++) {
fprintf(stderr, "%d>>%d\n",i,ogrpInfo->root);
}
*/

		ogrpInfo->root = NULL;
/*
*/
	}
}

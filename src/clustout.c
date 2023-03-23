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
#include "spec.h"
#include "plist.h"

#define MAXLINE 8000
#define WITH_POSITION 1

/*
static char head[MAXLINE];
static Node *rootnode, *save_rootnode;
*/

static struct {
	char head[MAXLINE];
	Node *rootnode, *save_rootnode, *currRoot;
	int root_flag;
	Hash *cinfoHash, *outergrpHash;
} statData;

int outputHier0(Node *node, int lev, signed char dir);
int clearNodeOutMarks(pList *nodelist, NodeFlag flag);
int outputNewick0(TreeNode *node, Dist prevdist);
int printDomInfo0(Domain *dom, char flag);
static pList **clustTab;

#ifdef WITH_NEIGHBOROUT
outputNeighborNode(NodeSet *nodes, Node *node)
{
	Node *nbrnode, *nn, *n;
	Neighbor *nbr;
	signed char dir;
	int flag = 0;
	static int posclustnum, clustnum;
	listIter iter;
	pList *nbrlist;
	int id, cnt;

	if (! isRoot(node) || (Opt.minent > 1 && ! node->child)) {
		return 0;
	}

	nbrlist = nodeNeighbor_gapsearch(nodes, node);
	setListIter(&iter, nbrlist, 1);
	while (n = (Node *) getListIter(&iter)) {
		if (! isRoot(n) || (Opt.minent > 1 && ! n->child)) {
			continue;
		}
		if (! flag) {
			printf("PosCluster %d\n", ++posclustnum);
			clustnum = 0;
			flag = 1;
		}
		statData.rootnode = n;
		printf("Cluster %d\n", ++clustnum);
		outputHier(n);
		printf("//\n");
		/* clear printed node */
		setFlagNode(n, NODE_DELETED);
/*
		if (checkNbrClust(n, dir, &nbrnode, &nbr) == 0) {
			break;
		}
		n = nbrnode;
*/
	}
}
#endif

static struct {
	int hom_clnum;
	int outer_clnum;
	int ortho_clnum;
	int count_clusters;
} Counter;

outputClusterInfo(ClusterInfoData *cInfoData, NodeSet *nodes)
{
	Domain *dom;
	int i, clsize, spcnt;
	int prev_homclust = -1;
	int prev_outerclust = -1;
	int ortho_clustid;
	char str_clustid[20];
	int subclustid;
	SubClusterInfo *subInfo;
	static char buf[64];
	ClusterInfo *cInfo;
	int clustnum = cInfoData->clustnum;

	statData.cinfoHash = cInfoData->cinfoHash;
	

	if (Opt.outstyle == DOMAINOUT) {
		return 0;
	}

	if (Opt.outstyle == CLUSTTAB) {
		outputSpNamesForClustTab();
		alloc_clustTabList();
	}

	Counter.ortho_clnum = 1;
	if (Opt.outstyle == TREEDETAIL) {
		outputHeader();
	}
	for (i = 0; i < clustnum; i++) {
 		cInfo = cInfoData->cinfo_link[i];
		if (cInfo->clustid < 0) continue;
/***
		clsize = numelemList(cInfo->members);
		if (clsize == 0) {
			continue;
		}
		spcnt = (int) spFlagCntW_All(cInfo->spflag);
		if (spcnt < Opt.minsp || clsize < Opt.minent) {
			continue;
		}
		f (Opt.outgroupMode && outgroupFlagCheck(cInfo->spflag)==2){
			continue;
		}
***/

		ortho_clustid = cInfo->clustid;
		str_clustid[0] = '\0';
		if (Opt.homClustOut && cInfo->homclust != prev_homclust) {
			++Counter.hom_clnum;
			prev_homclust = cInfo->homclust;
			if (Opt.outstyle != TABOUT && Opt.outstyle != CLUSTTAB) {
				int scoreout_flag = 0;
				sprintf(buf, "HomCluster %d", Counter.hom_clnum);
				if (Opt.outstyle == TREEDETAIL) {
					printNodeScoreFmt(cInfo->homclustRoot,
						stdout, buf, 2);
				} else {
					printf("%s", buf);
					if (Opt.outputScore &&
						cInfo->homclustRoot->child) {
						printNodeScore(
							cInfo->homclustRoot,
							Opt.outputScoreFp, buf);
					}
				}
				printf("\n");
			}
		}
		if (IS_NEW_OUTGROUP_MODE(Opt.outgroupMode) && cInfo->outerclust != prev_outerclust) {
			++Counter.outer_clnum;
			prev_outerclust = cInfo->outerclust;
			if (Opt.outstyle != TABOUT && Opt.outstyle != CLUSTTAB) {
				printf("OuterCluster %d\n", Counter.outer_clnum);
			}
		}
		if (Opt.outstyle != TABOUT && Opt.outstyle != CLUSTTAB) {
			sprintf(buf, "Cluster %d", ortho_clustid);
			printf("%s", buf);
			if (Opt.outputScore && cInfo->root->node->child) {
				printNodeScore(cInfo->root->node, Opt.outputScoreFp, buf);
			}
			printf("\n");
		} else {
			/** Opt.outstyle == TABOUT **/
			if (Opt.homClustOut) {
				sprintf(str_clustid, "%d|", Counter.hom_clnum);
			}
			sprintf(&str_clustid[strlen(str_clustid)],
					"%d", ortho_clustid);
		}
		switch (Opt.outstyle) {
		case TREEOUT:
		case TREEDETAIL:
		case GRAPHOUT:
		case NEWICK:
		case NEWICK_DIST:
		case HIER_DETAILOUT:
			if (IS_OUTGROUP_MODE(Opt.outgroupMode) && Opt.outputSubClustTree) {
				outputSubClusterInfo(cInfo, nodes);
				if (Opt.outstyle != TABOUT) {
					putchar('\n');
				}
			} else {
				outputClusterInfoMulti(
					cInfo, clustnum, nodes);
			}
			if (Opt.outputScoreFp && Opt.outgroupMode) {
				subclustid = 1;
				while (subInfo = (SubClusterInfo *)
						shiftList(cInfo->ingroups)) {
					sprintf(buf, "SubCluster %d", subclustid);
					printNodeScore(subInfo->root->node,
						Opt.outputScoreFp, buf);
					subclustid++;
				}
			}
			break;
		default:
			/** default output **/
			/** output members **/
			if (! Opt.outgroupMode || Opt.outputSubClustTree) {
				outputClustInfoLines(cInfo->members, str_clustid, 0, NULL);

			} else {    /** outgroupMode **/
				subclustid = 1;
				while (subInfo = (SubClusterInfo *)
					shiftList(cInfo->ingroups)) {
					outputClustInfoLines(subInfo->members,
							str_clustid, subclustid,
							subInfo->root);
					subclustid++;
				}
				/** OutGroup **/
				if (IS_OUTGROUP_MODE(Opt.outgroupMode)) {
					outputClustInfoLines(
						cInfo->outgroup->members,
						str_clustid, -1,
						cInfo->outgroup->root);
					if (cInfo->outgroup2 && numelemList(cInfo->outgroup2->members) > 0) {
					    outputClustInfoLines(
						cInfo->outgroup2->members,
						str_clustid, -2,
						cInfo->outgroup->root);
					}
				}
			}
			break;
		}
		if (Opt.taxMapOut == 1) {
			printf("======\n");
			mapTaxInfo_Cinfo(cInfo);
			printf("======\n\n");
		}
	}
}
outputClustInfoLines(pList *cInfoList, char *str_clustid, int subclustid,
		TreeNode *root)
{
	Domain *dom;
	static char buf[64];
	if (Opt.outstyle != TABOUT && Opt.outstyle != CLUSTTAB) {
		if (subclustid > 0) {
			sprintf(buf, "SubCluster %d", subclustid);
			printf("%s",buf);
			if (Opt.outputScore && root && root->node->child) {
				printNodeScore(root->node,Opt.outputScoreFp,buf);
			}
			putchar('\n');
		} else if (subclustid < 0) {
			if (subclustid == -1) {
				printf("OutGroup\n");
			} else {
				printf("OutGroup2\n");
			}
			sprintf(&str_clustid[strlen(str_clustid)], ".0");
		}
	} else {
		if (subclustid > 0) {
			sprintf(buf, "%s.%d ", str_clustid, subclustid);
		} else if (subclustid < 0) {
			sprintf(buf, "%s.0 ", str_clustid);
		} else {
			sprintf(buf, "%s ", str_clustid);
		}
	}
	if (Opt.outstyle == CLUSTTAB) {
		outputClustTab_members(cInfoList, buf);
	} else {
	    while (dom = (Domain *) shiftList(cInfoList)) {
		if (Opt.outstyle == TABOUT) {
/*
			if (subclustid > 0) {
				printf("%s.%d ", str_clustid, subclustid);
			} else if (subclustid < 0) {
				printf("%s.0 ", str_clustid);
			} else {
				printf("%s ", str_clustid);
			}
*/
			printf("%s ", buf);
		} else {
		}
		printDomInfo(dom);
		putchar('\n');
	    }
	}
	if (Opt.outstyle != TABOUT) {
		printf("\n");
	}
}
outputClusterInfoMulti(ClusterInfo *cInfo, int clustnum, NodeSet *nodes)
{
	listIter iter;
	TreeNode *tnode;
	int tn = 0;
	pList *nodelist;
	int testflag = 1;
	int clnum;
	ClusterInfo *cinfo0;
/*
	pList *clusters_List;
*/

/*
	if (IS_NEW_OUTGROUP_MODE(Opt.outgroupMode) && numelemList(cInfo->clusters_outer)) {
		clusters_List = cInfo->clusters_outer;
	} else {
		clusters_List = cInfo->clusters;
	}
*/

/*
	clnum = numelemList(cInfo->clusters);
	setListIter(&iter, cInfo->clusters, 1);
*/
	clnum = numelemList(cInfo->clusters);
/*
printf("clist>>>>%d,clnum=%d\n",cInfo->clusters, clnum);
*/
	setListIter(&iter, cInfo->clusters, 1);
	if (testflag) {
		nodelist = create_pList();
/*
		while (tnode = (TreeNode *) getListIter(&iter)) {
*/
		while (cinfo0 = (ClusterInfo *) getListIter(&iter)) {
			if (IS_NEW_OUTGROUP_MODE(Opt.outgroupMode) && cinfo0->outer_root) {
				tnode = cinfo0->outer_root;
			} else {
				tnode = cinfo0->origroot;
			}
			testTreeOverlap(tnode, nodelist);
		}
	}

	if (IS_NEW_OUTGROUP_MODE(Opt.outgroupMode) && cInfo->outgroup2 && numelemList(cInfo->outgroup2->members)>0) {
		listIter iter_ogrp;
		Domain *dom;
		HENTRY hent;
		int grpsiz = numelemList(cInfo->outgroup2->members);
		if (grpsiz % 2 == 0) grpsiz++;
		/* cInfo->outgroup2 contains outer-group nodes */
		setListIter(&iter_ogrp, cInfo->outgroup2->members, 1);
		statData.outergrpHash = Hcreate(grpsiz * 17);
		/** create hash for searching outer node */
		while(dom = (Domain*) getListIter(&iter_ogrp)) {
 			hent.key = (char*) dom->leaf->id;
			HIsearch(statData.outergrpHash, &hent, ENTER);
		}
	}


	setListIter(&iter, cInfo->clusters, 1);
/*
	while (tnode = (TreeNode *) getListIter(&iter)) {
*/
	while (cinfo0 = (ClusterInfo *) getListIter(&iter)) {

		statData.root_flag = 0;
		if (IS_NEW_OUTGROUP_MODE(Opt.outgroupMode) && cinfo0->outer_root) {
			tnode = cinfo0->outer_root;
		} else {
			tnode = cinfo0->origroot;
			statData.root_flag = 1;
		}
/*
printf(">>>>>>");printTreeNode(tnode);printf(" ");printTreeNode(cinfo0->origroot);printf(" ");printTreeNode(cinfo0->root);printf(">");printTreeNode(cInfo->origroot);printf("\n");
printf("root>>");printNode(cinfo0->origroot->node);printf("%d,%d\n",cinfo0->homclust,cinfo0->outerclust);
printf("outer_root>>");printTreeNode(tnode);printf("%d,%d\n",isOuterTreeRoot(tnode),isTreeRoot(tnode));
*/
		statData.currRoot = cinfo0->origroot->node;
		outputCluster_sub(nodes, tnode, cInfo->root, 0);
		if (Opt.outstyle != TABOUT) {
			putchar('\n');
		}
		if (++tn < clnum){
			if (Opt.outstyle != NEWICK &&
					Opt.outstyle != NEWICK_DIST) {
				printf("--");
			}
			printf("\n");
		}
	}
	if (Opt.outstyle == TREEDETAIL) {
		printf("//\n");
	}
	if (testflag) {
		clearNodeOutMarks(nodelist, NODE_TMPMARK);
		clearNodeOutMarks(nodelist, NODE_DUPMARK);
		free_pList(nodelist);
	}
	if (IS_NEW_OUTGROUP_MODE(Opt.outgroupMode) && cInfo->outgroup2 && cInfo->outgroup2->members) {
		listIter *iter;
		Domain *dom;
		setListIter(&iter, cInfo->outgroup2->members);
		while (dom = (Domain*) getListIter(&iter)) {
			printf("# Additional outgroup: ");
			printDomInfo(dom);
			printf("\n");
		}
		printf("\n");
		if (statData.outergrpHash) {
			Hdestroy(statData.outergrpHash);
			statData.outergrpHash = NULL;
		}
	}
}

outputSubClusterInfo(ClusterInfo *cInfo, NodeSet *nodes)
{
	SubClusterInfo *subInfo;
	while (subInfo = (SubClusterInfo *) shiftList(cInfo->ingroups)) {
		printf("SubCluster %d\n", subInfo->clustid);
		outputCluster_sub(nodes, subInfo->root, cInfo->root, 0);
		if (Opt.outstyle != TABOUT) {
			putchar('\n');
		}
	}
}
outputOuterClusterInfo(pList *outerRoots, NodeSet *nodes)
{
	listIter iter;
	SubClusterInfo *outerInfo;
	setListIter(&iter, outerRoots, 1);
	while (outerInfo = (SubClusterInfo *) getListIter(&iter)) {
		if (outgroupFlagCheck(outerInfo->root->node->spflag)==2) {
			continue;
		}
		if (outerInfo->level == OuterCluster) {
			printf("OuterCluster\n");
		} else {
			printf("Cluster\n");
		}
		outputCluster_sub(NULL, outerInfo->root, outerInfo->root, 0);
		printf("\n");
	}
}

/**
outputClusterData(pList *roots, NodeSet *nodes)
{
	int i;
	TreeNode *tnode;
	listIter iter;
	int spcnt;

	if (Opt.outstyle == DOMAINOUT) {
		return;
	}
	if (Opt.outstyle == CLUSTTAB) {
		outputSpNamesForClustTab();
		alloc_clustTabList();
	}

	setListIter(&iter, roots, 1);
	while (tnode = (TreeNode *) getListIter(&iter)) {
		Node *node = tnode->node;

		spcnt = (int) spFlagCntW_All(node->spflag);
/ *
printf(">spcnt=%d,cnt=%d; %d,%d\n", spcnt, node->cnt,Opt.minsp,Opt.minent);
* /
		if (spcnt < Opt.minsp || node->cnt<Opt.minent) {
			continue;
		}
		if (Opt.outgroupMode && outgroupFlagCheck(node->spflag)==2){
			/ * only outgroup species * /
			continue;
		}
		if (Opt.outstyle != NEIGHBOR) {
			if (Opt.homClustOut && Counter.count_clusters == 0) {
				printf("HomCluster %d\n", ++Counter.hom_clnum);
			}
			if (Opt.outstyle != TABOUT && Opt.outstyle != CLUSTTAB){
				printf("Cluster %d", ++Counter.ortho_clnum);
				putchar('\n');
			}
			++Counter.count_clusters;
		}
		outputCluster_sub(nodes, tnode, tnode, Counter.ortho_clnum);
		if (Opt.outstyle != TABOUT) {
			putchar('\n');
		}
	}
	if (Opt.outstyle == CLUSTTAB) {
		free_clustTabList();
	}
}
**/
clearNodeOutMarks(pList *nodelist, NodeFlag flag)
{
/*
	listIter iter;
	Node *n;
	setListIter(&iter, nodelist, 1);
	while (n = (Node *) getListIter(&iter)) {
		unsetFlagNode(n, flag);
	}
*/
	/* move to graph.c */
	clearAllMarks_sub(nodelist, flag);
}
testTreeOverlap(TreeNode *tnode, pList *nodelist)
{
	Node *child1 = NULL, *child2 = NULL;
	
	if (testFlagNode(tnode->node, NODE_TMPMARK)) {
		/* an overlapping node */
		setFlagNode(tnode->node, NODE_DUPMARK);
	} else {
		setFlagNode(tnode->node, NODE_TMPMARK);
		pushList(nodelist, tnode->node);
	}
	if (tnode->child1) {
		testTreeOverlap(tnode->child1, nodelist);
	}
	if (tnode->child2) {
		testTreeOverlap(tnode->child2, nodelist);
	}
/*
	getChilds(node, &child1, &child2);
	if (child1) {
		testTreeOverlap(child1, nodelist);
	}
	if (child2) {
		testTreeOverlap(child2, nodelist);
	}
*/
}
outputCluster_sub(NodeSet *nodes, TreeNode *tnode, TreeNode *root, int clustnum)
{
	if  (root) {
		statData.rootnode = root->node;
	}
	switch (Opt.outstyle) {
	case TREEOUT:
		statData.head[0] = '\0';
		outputTree(tnode);
		break;
	case TREEDETAIL:
		outputTreeDetail(tnode);
		break;
	case GRAPHOUT:
		outputGraph(tnode->node);
		break;
#ifdef NEIGHBOROUT
	case NEIGHBOR:
		outputNeighborNode(nodes, tnode->node);
		return 0;
		break;
#endif
	case NEWICK:
	case NEWICK_DIST:
		outputNewick(tnode);
		break;
	case TABOUT:
		outputSimpleTab(tnode, clustnum);
		break;
	case NORMALOUT:
		outputNormal(tnode);
		break;
	case CLUSTTAB:
		outputClustTab(tnode, clustnum);
		break;
	case HIER_DETAILOUT:
/*
	case NORMALOLD:
*/
	default:
		outputHier(tnode->node);
		break;
	}
}

outputNormal(TreeNode *tnode)
{
	if (IS_OUTGROUP_MODE(Opt.outgroupMode)) {
		outputNormal_outgroup(tnode);
	} else {
		outputNormal_normal(tnode);
	}
}
outputNormal_normal(TreeNode *tnode)
{
	pList *members = create_pList();
	Domain *dom;
	collectCluster(tnode, members);
	while (dom = (Domain *) shiftList(members)) {
		printDomInfo(dom);
		printf("\n");
	}
	freeList(members);
}
outputNormal_outgroup(TreeNode *tnode)
{
	pList *ingroups = create_pList();
	SubClusterInfo *ingroup;
	SubClusterInfo *outgroup = NULL;
	listIter iter;
	Domain *dom;
	int subclustid = 1;
	if (IS_OUTGROUP_MODE(Opt.outgroupMode)) {
		outgroup = createSubClusterInfo(NULL);
	}

	collectCluster_outgroup(tnode, ingroups, outgroup);
	while (ingroup = (SubClusterInfo *) shiftList(ingroups)) {
		printf("SubCluster %d\n", subclustid++);
		while (dom = (Domain *) shiftList(ingroup->members)) {
			printDomInfo(dom); printf("\n");
		}
		printf("\n");
	}
	if (IS_OUTGROUP_MODE(Opt.outgroupMode)) {
		printf("OutGroup\n");
		setListIter(&iter, outgroup->members, 1);

		while (dom = (Domain *) getListIter(&iter)) {
			printDomInfo(dom); printf("\n");
		}
	}
}
alloc_clustTabList()
{
	int i;
	if ( (clustTab = (pList **) malloc(sizeof(pList*) * SPnum))==NULL ){
		allocError("clusters");
	}
	for (i = 0; i < SPnum; i++) {
		if ( (clustTab[i] = create_pList())==NULL ) {
			allocError("clustTab");
		}
	}
}
free_clustTabList()
{
	int i;
	for (i = 0; i < SPnum; i++) {
		free(clustTab[i]);
	}
	free(clustTab);
}
outputSpNamesForClustTab()
{
	int i;
	char *spname;
	printf("#id");
	for (i = 0; i < SPnum; i++) {
		spname = getSPname(i);
		printf("\t%s",spname);
	}
	printf("\n");
}
outputClustTab(TreeNode *tnode, int clustnum)
{
	pList *members = create_pList();
	char str_clustnum[15];
	sprintf(str_clustnum, "%d", clustnum);

	collectCluster(tnode, members);
	outputClustTab_members(members, str_clustnum);
}
outputClustTab_members(pList *members, char *str_clustnum)
{
	Domain *dom;
	Node *lnode;
	char spname[SPNAMELEN];
	char *p, *q;
	int i, j, spid;
	while (dom = (Domain *) shiftList(members)) {
		lnode = dom->leaf;
		for (p = lnode->name, q = spname; *p && *p != ':'; p++, q++) {
			*q = *p;
		}
		*q = '\0';
		spid = getSPid(spname);
		pushList(clustTab[spid], dom);
	}
	printf("%s",str_clustnum);
	for (i = 0; i < SPnum; i++) {
		j = 0;
		printf("\t");
		while (dom = (Domain *) shiftList(clustTab[i])) {
			if (j++ > 0) printf(" ");
			printDomName(dom);
		}
	}
	freeList(members);
}
outputSimpleTab(TreeNode *tnode, int clustnum)
{
	pList *members = create_pList();
	Domain *dom;
	collectCluster(tnode, members);
	while (dom = (Domain *) shiftList(members)) {
		printDomInfo(dom);
		printf("\n");
	}
	freeList(members);
}
printNodeInfo(Node *node, Node *root)
{
	Domain *dom;
	int domn;
	if (node->domains) {
		domn = getDomainNoMark(node->domains, root, &dom);
		if (domn) {
			printf("%s", node->name);
			if (numelemList(node->domains) > 1) {
				printf("(%d)", domn);
			}
			printf(" %d %d",(int)dom->from, (int)dom->to);
		} else {
			printf("(%s deleted)",node->name);
/*
			printf(" %d %d",1,(int)node->len);
*/
		}
	} else {
		printf("%s %d %d",node->name,1,(int)node->len);
	}
}
printDomName(Domain *dom)
{
	printDomInfo0(dom, 0);
}
printDomInfo(Domain *dom)
{
	printDomInfo0(dom, WITH_POSITION);
}
printDomInfo0(Domain *dom, char flag)
{
	Node *node = dom->leaf, *root = dom->root;
	int domn;
	printf("%s", node->name);
	if (node->domains) {
		if (dom->num) {
			if (Opt.outstyle == TABOUT) {
				printf(" %d", dom->num);
			} else if (numelemList(node->domains) > 1) {
				printf("(%d)", dom->num);
			}
			if (flag & WITH_POSITION) {
			    printf(" %d %d",(int)dom->from, (int)dom->to);
			}
/*
printf("FOUND_DOMNUM:%d,%d\n",dom->num,dom);
*/
		} else {
			domn = getDomainNoMark(node->domains, root, &dom);
/*
printf("NO_DOMNUM:%d,%d,%d\n",domn,dom->num,dom);
*/
			if (domn) {
				if (numelemList(node->domains) > 1) {
					printf("(%d)", domn);
				}
				if (flag & WITH_POSITION) {
				    printf(" %d %d",(int)dom->from, (int)dom->to);
				}
			} else if (flag & WITH_POSITION) {
				printf(" %d %d",1,(int)node->len);
			}
		}
	} else if (flag & WITH_POSITION) {
		printf(" %d %d",1,(int)node->len);
	}
}
printNodeScore(Node *node, FILE *fp, char *string)
{
	int format = 0;
	if (fp == NULL) {
		fp = stdout;
		format = 1;
	}
	printNodeScoreFmt(node, fp, string, format);
}
printNodeScoreFmt(Node *node, FILE *fp, char *string, int format)
{
	if (format == 0 || format == 2) {
		if (string) {
			fprintf(fp, "%s", string);
		}
	}
	if (node->child) {
		if (format==0) {
			fprintf(fp, ":score=%f:dist=%f\n",
					node->child->score, node->child->dist);
		} else if (format==1) {
			fprintf(fp, " [score=%f dist=%f]",
					node->child->score, node->child->dist);
		} else if (format==2) {
			fprintf(fp, " : %f %f",
				node->child->score, node->child->dist);
		}
	}
}


outputHier(Node *node)
{
	outputHier0(node, 1, (signed char) 1);
}
outputHier0(Node *node, int lev, signed char dir)
{
	int i;
	Edge *child;
	char indent[200];
	Node *maxnode;
	int maxcnt;
	Domain *dom;

	indent[0] = '\0';
	if (Opt.outstyle == HIER_DETAILOUT) {
		for (i = 0; i < lev; i++) {
			strcat(indent, "  ");
		}
	}
	if (isFlankNode1(node)) {
		if (Opt.outstyle == HIER_DETAILOUT) {
			printf("%s - %d (%s,%d) %d/%d (flanking node 1)\n",
				indent, lev, node->name,node->id,
				nodeLen(node),node->len);
		}
		outputHier0(node->child->node1, lev, dir);
	} else if (isFlankNode2(node)) {
		if (Opt.outstyle == HIER_DETAILOUT) {
			printf("%s - %d (%s,%d) %d/%d (flanking node 2)\n",
				indent, lev, node->name,node->id,
				nodeLen(node),node->len);
		}
#ifdef LARGE
		dir *= node->child->dir;
#endif
		outputHier0(node->child->node2, lev, dir);
	} else if (isIntNode(node)) {
		child = node->child;
		if (Opt.outstyle == HIER_DETAILOUT) {
			ConnCount conn = 0;
#ifdef EXPERIMENTAL
			conn = child->connect;
#endif
			printf("%s", indent);
			printf("[%d](%s,%d) [%d:%d] %d/%d %.2f %d  ",
				lev,node->name,node->id,
				node->newreg.from,node->newreg.to,
				nodeLen(node), node->len,
				MEASURE(child),conn);
			print_specFlag(node->spflag);

			if (node->left) {
				printf("%s", indent);
				printf("    Left: ");
/*
				printList(node->left, printNode);
				printNeighborList(node->left, node);
*/
				maxcnt = getMaxNeighbor(node->left,node,1,
					&maxnode,NULL);
				if (maxcnt && maxcnt >=
					(double) node->cnt * Opt.nbrConnRatio) {
					printNode(maxnode);
					printf(" [%d]",maxcnt);
				}
				putchar('\n');
			}
			if (node->right) {
				printf("%s", indent);
				printf("    Right: ");
/*
				printList(node->right, printNode);
				printNeighborList(node->right, node);
*/
				maxcnt = getMaxNeighbor(node->right,node,1,
					&maxnode,NULL);
				if (maxcnt && maxcnt >=
					(double) node->cnt * Opt.nbrConnRatio) {
					printNode(maxnode);
					printf(" [%d]",maxcnt);
				}
				putchar('\n');
			}
		}
		outputHier0(child->node1, lev+1, dir);
#ifdef LARGE
		dir *= node->child->dir;
#endif
		outputHier0(child->node2, lev+1, dir);
	} else {
		printf("%s", indent);
		printNodeInfo(node,statData.rootnode);
		if (Opt.revMatch) {
			printf(" %d", dir);
		}
		putchar('\n');
	}
}
outputGraph(Node *node)
{
	Node *maxnode = NULL;
	enum {Leaf,NonLeaf} type = NonLeaf;
	int maxcnt;

	if (isFlankNode1(node)) {
		if (isIntNode(node->child->node1)) {
			printf("%d %d\n", node->id, node->child->node1->id);
		} else {
			printf("%d %s\n", node->id, node->child->node1->name);
		}
		outputGraph(node->child->node1);
	} else if (isFlankNode2(node)) {
		if (isIntNode(node->child->node2)) {
			printf("%d %d\n", node->id, node->child->node2->id);
		} else {
			printf("%d %s\n", node->id, node->child->node2->name);
		}
		outputGraph(node->child->node2);
	} else if (isIntNode(node)) {
		if (isIntNode(node->child->node1)) {
			printf("%d %d\n", node->id, node->child->node1->id);
		} else {
			printf("%d %s\n", node->id, node->child->node1->name);
		}
		if (isIntNode(node->child->node2)) {
			printf("%d %d\n", node->id, node->child->node2->id);
		} else {
			printf("%d %s\n", node->id, node->child->node2->name);
		}
		outputGraph(node->child->node1);
		outputGraph(node->child->node2);
	} else {
		type = Leaf;
	}

	if (node->left) {
		maxcnt = getMaxNeighbor(node->left,node,0,&maxnode,NULL);
		if (maxcnt && maxcnt >=
				rint(node->cnt * Opt.nbrConnRatio)) {
			if (type == NonLeaf) {
				printf("%d %d L\n", node->id, maxnode->id);
			} else {
				if (isLeaf(maxnode) && ! isRoot(maxnode)) {
					printf("%s %s L\n",node->name,maxnode->name);
				}
			}
		}
	}
	if (node->right) {
		maxcnt = getMaxNeighbor(node->right,node,0,&maxnode,NULL);
		if (maxcnt && maxcnt >=
				rint(node->cnt * Opt.nbrConnRatio)) {
			if (type == NonLeaf) {
				printf("%d %d R\n", node->id, maxnode->id);
			} else {
				if (isLeaf(maxnode) && ! isRoot(maxnode)) {
					printf("%s %s R\n",node->name,maxnode->name);
				}
			}
		}
	}
}

outputTree(TreeNode *tnode)
{
	if (tnode->child1 && ! tnode->child2) {
		tnode = tnode->child1;
	} else if (tnode->child2 && ! tnode->child1) {
		tnode = tnode->child2;
	}
	outputTree_sub(tnode, 0, 1, 1);
}
outputTree_sub(TreeNode *tnode, int lev, int branch_dir, int gene_dir)
{

	Edge *child;
	Domain *dom;
	Node *node;
	char markchar;
	int setroot_flag = 0;
	int outer_flag = 0;
	static char clname_buf[100];
	char *clname = "[undefined cluster]";

	if (tnode == NULL) {
		printf("%s%c%s", statData.head, '+', "- (null)\n");
		return -1;
	}

 	node = tnode->node;
	markchar = ( (node->flag & NODE_DUPMARK) ? '*' : '+');
	dom = tnode->dom;

/*
if (tnode->child1||tnode->child2) {
printf("tnode=%d,child1=%d,child2=%d\n",tnode,tnode->child1,tnode->child2);
}
*/
	if (IS_NEW_OUTGROUP_MODE(Opt.outgroupMode) && ! Opt.outputOuterClustTree) {
		if (isTreeRoot(tnode)) {
			if (node != statData.currRoot && ! statData.root_flag) {
				ClusterInfo *cinfo0;
				int clustid = -999;
				Node *prev_rnode;
				Node *rnode = node;
/*
printNode(node); printf(" "); printNode(statData.currRoot);printf("\n");
*/
				/* different root node */
				if (cinfoHsearch(statData.cinfoHash, rnode->id, &cinfo0)) {
					clustid = cinfo0->clustid;
					while (clustid < 0) {
						prev_rnode = rnode;
						rnode = cinfo0->root->node;
						if (rnode == prev_rnode) {
/*
							fprintf(stderr, "Cannot found root: %d\n",node->id);
*/
							break;
						}
						if (cinfoHsearch(statData.cinfoHash, rnode->id, &cinfo0)) {
							clustid = cinfo0->clustid;
						}
					}
				}
				if (clustid >= 0) {
					clname = clname_buf;
					sprintf(clname, "[cluster %d]", clustid);
				}
				printf("%s%c- %s\n", statData.head, markchar, clname);
				return 0;
			} else if (! statData.root_flag) {
				statData.root_flag = 1;
				setroot_flag = 1;
			}
		}
	}

	if (ClustTree_isLeaf(tnode)) {
		/** leaf **/
		printf("%s%c%s", statData.head, markchar, "- ");
		printDomInfo(tnode->dom);
#ifdef LARGE
		printf(" %d", gene_dir);
#endif

		if (IS_OUTGROUP_MODE(Opt.outgroupMode)) {
			if (dom) {
			    if (dom->subgrp >= 0) {
				if (dom){
					printf(" [I%d]", dom->subgrp);
/*
					printf("; %d,%d", tnode,dom);
*/
					if(outgroupFlagCheck(node->spflag)) {
						printf(" *");
					}
				}
			    } else {
				HENTRY hent;
				if (statData.outergrpHash) {
					hent.key= (char*) node->id;
					if (HIsearch(statData.outergrpHash, &hent, FIND)) {
						outer_flag = 1;
					} else {
						outer_flag = 0;
					}
				}
				if (outer_flag) {
					printf(" [OO]");
				} else if (statData.root_flag) {
					printf(" [O]");
				} else {
					printf(" [o]");
				}
			    }
			} else {
			}
		}
		putchar('\n');
	} else {
		child = node->child;
		if (child) {
			statData.head[lev*2] = '\0';
			strcat(statData.head, (lev > 0 && branch_dir == -1 ? "| " : "  "));
			outputTree_sub(tnode->child1, lev+1, 1, gene_dir);
			statData.head[lev*2] = '\0';
			printf("%s%c%s", statData.head, markchar, "-| ");
/*
printf(">>>%d,%d[%d]<<<",isOuterRoot(node),isDeletedOuterRoot(node),node->id);
*/
			if (Opt.decinum) {
				printf("%.*f\n", Opt.decinum,
					MEASURE(child));
			} else {
				printf("%.1f\n", MEASURE(child));
			}

#ifdef LARGE
			gene_dir *= node->child->dir;
#endif
	
			strcat(statData.head, (lev > 0 && branch_dir == 1 ? "| " : "  "));
			outputTree_sub(tnode->child2, lev+1, -1, gene_dir);
		} else {
			printf("%s%c%s[%s]", statData.head, markchar,
				"- ", node->name);
			printf(" %d %d",1,(int)node->len);
			putchar('\n');
		}
	}
	if (setroot_flag) {
		statData.root_flag = 0;
	}
	return 0;
}
outputTreeDetail(TreeNode *tnode)
{
	if (tnode->child1 && ! tnode->child2) {
		tnode = tnode->child1;
	} else if (tnode->child2 && ! tnode->child1) {
		tnode = tnode->child2;
	}
	outputTreeDetail_sub(tnode, NULL);
}
outputTreeDetail_sub(TreeNode *tnode, TreeNode *parent)
{

	Edge *child;
	Domain *dom;
	Node *node;
	char markchar;
	int pid = -1;

	if (tnode == NULL) return(0);

 	node = tnode->node;
	markchar = ( (node->flag & NODE_DUPMARK) ? '*' : '+');
	dom = tnode->dom;
	
	if (isTreeRoot(tnode)) {
		statData.save_rootnode = statData.rootnode;
		statData.rootnode = node;
	}
	if (parent != NULL && parent->node != NULL) {
		pid = parent->node->id;
	}

	if (ClustTree_isLeaf(tnode)) {
		/** leaf **/
		printf("L %d %d %s ", node->id, pid, node->name);
		printf("%d %d ", tnode->dom->from, tnode->dom->to);
		printf("%d", (numelemList(node->domains) > 1) ?  dom->num : 0);

/*
		if (IS_OUTGROUP_MODE(Opt.outgroupMode)) {
			if (dom) {
			    if (dom->subgrp >= 0) {
				if (dom){
					printf(" [I%d]", dom->subgrp);
					if(outgroupFlagCheck(node->spflag)) {
						printf(" *");
					}
				}
			    } else {
				printf(" [O]");
			    }
			} else {
			}
		}
*/
		putchar('\n');
	} else {
		printf("I %d %d %s ", node->id, pid, node->name);
		child = node->child;
		if (child) {
			if (Opt.decinum) {
				printf("%.*f %.*f", Opt.decinum, child->score,
						Opt.decinum, child->dist);
			} else {
				printf("%.1f %.1f", child->score, child->dist);
			}
			putchar('\n');
	
			outputTreeDetail_sub(tnode->child1, tnode);
			outputTreeDetail_sub(tnode->child2, tnode);
		} else {
/*
			printf("%s%c%s[%s]", statData.head, markchar,
				"- ", node->name);
			printf(" %d %d",1,(int)node->len);
			putchar('\n');
*/
		}
	}
	if (isTreeRoot(tnode)) {
		statData.rootnode = statData.save_rootnode;
	}
}
outputHeader()
{
	printf("# simtype=%d\n",Opt.sim);
	printf("# cutoff=%f\n",Opt.cutoff);
	if (Opt.sim) { /* similarity */
		printf("# missdist=%f\n",Opt.missscore);
	} else {
		printf("# missdist=%f\n",Opt.missdist);
	}
	printf("# phylocutratio=%f\n",Opt.phylocutratio);
	printf("# distscale=%d\n",Opt.distscale);
}

static char namebuf[NAMELEN];
outputNewick(TreeNode *tnode)
{
	if (tnode->child1 && ! tnode->child2) {
		tnode = tnode->child1;
	} else if (tnode->child2 && ! tnode->child1) {
		tnode = tnode->child2;
	}
	outputNewick0(tnode, ((Dist) -1.0) );
	printf(";");
}
outputNewick0(TreeNode *tnode, Dist prevdist)
{
	Dist dist = 0;
	Dist currdist = 0;

	if (tnode == NULL) {
	} else if (ClustTree_isLeaf(tnode)) {
		if (numelemList(tnode->node->domains) > 1) {
			sprintf(namebuf, "%s#%d", tnode->node->name,
						tnode->dom->num);
		} else {
			strcpy(namebuf, tnode->node->name);
		}
		convname(namebuf);
		if (Opt.outstyle == NEWICK_DIST) {
			printf("%s", namebuf);
			if (! Opt.sim && prevdist >= 0) {
				dist = prevdist;
				if (Opt.decinum) {
					printf(":%.*f ", Opt.decinum,
						(float)dist);
				} else {
					printf(":%.1f ",  (float)dist);
				}
			}
		} else {
			printf("%s", namebuf);
		}
	} else {
		currdist = tnode->node->child->dist / 2;
		printf("(");
		outputNewick0(tnode->child1, currdist);
		printf(", ");
		outputNewick0(tnode->child2, currdist);
		printf(")");
		if (Opt.outstyle == NEWICK_DIST) {
			if (! Opt.sim && prevdist >= 0) {
				dist = prevdist - currdist;
				if (Opt.decinum) {
					printf(":%.*f ", Opt.decinum,
						(float)dist);
				} else {
					printf(":%.1f ",  (float)dist);
				}
			}
		}
	}
}
convname(char *name)
{
	register char *p;
	for (p = name; *p; p++) {
		if (*p == ':') {
			*p = '_';
		}
	}
}


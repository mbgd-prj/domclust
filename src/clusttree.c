#include<stdio.h>
#include "domclust.h"

#define BLKSIZ 5000

static Alloc_Object *clustTreeObj;

initClustTreeNode()
{
	if (clustTreeObj == NULL) {
		clustTreeObj = init_alloc_object(sizeof(TreeNode),BLKSIZ);
	}
}
TreeNode *createTreeNode(Node *node)
{
	TreeNode *tnode = (TreeNode *) memalloc(clustTreeObj);
	if (! tnode) {
		allocError("TreeNode");
	}
	tnode->node = node;
	tnode->child1 = tnode->child2 = NULL;
	tnode->dom = NULL;
	tnode->flag = node->flag;
/*
printf("createTreeNode:"); printTreeNode(tnode);printf(" "); printNode(node); printf(":%d",node->flag); printf("\n");
*/
	clearSPflag(tnode->spflag);
	return tnode;
}
TreeNode *createTreeNodeUnroot(Node *node)
{
	TreeNode *retnode = createTreeNode(node);
	resetCenterRoot_flag(retnode->flag);
	return retnode;
}
setTreeNode(TreeNode *treenode, Node *node)
{
	treenode->node = node;
}
ClustTree_addChild(TreeNode *parent, TreeNode *child, int childnum)
{
	if (childnum==1) {
		parent->child1 = child;
	} else {
		parent->child2 = child;
	}
	child->parent = parent;
}
ClustTree_isLeaf(TreeNode *tnode)
{
	if (! tnode) return 0;
	return ((tnode->dom) ? 1 : 0);
}
ClustTree_isIgnoredLeaf(TreeNode *tnode)
{
	return (tnode->dom == 0 && tnode->node->child == NULL);
}

ClustTree_alldom(TreeNode *tnode, pList *list)
{

/*
printf("ALLDOM:%d,%d\n", tnode, tnode->dom);
*/

	if (ClustTree_isLeaf(tnode)) {
/*
printf("Dom>>%d,%d:",tnode->dom->leaf->name,tnode->dom->leaf->id);
printDomain(tnode->dom);
printf("child>>%d>>%d\n",tnode->child1, tnode->child2);
printf("\n");
*/
		pushList(list, tnode->dom);
	} else {
		if (tnode->child1) {
			ClustTree_alldom(tnode->child1, list);
		}
		if (tnode->child2) {
			ClustTree_alldom(tnode->child2, list);
		}
	}
}

static Node *rootnode;

pList *convClustTree(NodeSet *nodes)
{
	Node *node;
	TreeNode *ret_tnode, *tnode;
	NodeSetIter iter;
	pList *root_treeNodes;
	int type;

	initClustTreeNode();
	getRootInit(&iter, nodes, -1, 1);
	root_treeNodes = create_pList();
	while (node = getRootNext(&iter)) {
		rootnode = node;
		type = convClustTree_sub(node, &ret_tnode);
		if (ret_tnode) {
			if (node != ret_tnode->node){
				if (type == 1 || type == 2) {
					tnode = createTreeNode(node);
					ClustTree_addChild(tnode, ret_tnode, type);
					pushList(root_treeNodes, tnode); 
				} else {
					if (outgroupFlagCheck(ret_tnode->node->spflag) != 2) {
/*
						makeRootNode_Out(ret_tnode->node);
*/
						makeRootTreeNode_Out(ret_tnode);
						resetNodeDomainAll(ret_tnode->node);
						pushList(root_treeNodes, ret_tnode); 
					}
					deleteRootNode(node);
				}
			} else {
				pushList(root_treeNodes, ret_tnode); 
			}
		}
/*
		printf("Cluster\n");
		printClustTree(ret_tnode);
		printf("\n");
*/
	}
	assignSpecFlag(root_treeNodes);
	return root_treeNodes;
}


int tmpcount;
pList *root_treeNodes;
Node *root00;

convClustTree_fromHom(NodeSet *nodes, pList *origRoots, pList **allClust_hom, pList **allClust_plain)
{
	Node *node;
	TreeNode *tnode;
/*
	NodeSetIter iter;
*/
	listIter iter;
	int type;
	SubClusterInfo *homInfo;
	SubClusterInfo *createHomClusterInfo();

	initClustTreeNode();
	setListIter(&iter, origRoots, 1); /* origRoots = HomClustRoots */
	*allClust_hom = create_pList();
	if (allClust_plain) {
		*allClust_plain = create_pList();
	}

	while (node = (Node *) getListIter(&iter)) {
/*
root00 = node;
*/
		root_treeNodes = create_pList();
		rootnode = NULL;
		findOrthoRoot(node, &tnode, NULL); /* orthoroots are gathered into root_treeNodes */
		if (tnode == NULL) {
			continue;
		}
		homInfo = createHomClusterInfo(tnode);
		homInfo->members = root_treeNodes;
		pushList(*allClust_hom, homInfo);
		if (allClust_plain) {
			pushListAll(*allClust_plain, root_treeNodes);
		}
		assignSpecFlag_sub(tnode);
	}
	clearAllMarks(origRoots);
/*
	if (allClust_plain) {
		assignSpecFlag(*allClust_plain);
	}
*/
}
makeClustTreeHom(pList *clustRoots, pList *origRoots, pList **allClust_hom)
{
	
	listIter iter;
	Node *node;
	setListIter(&iter, origRoots, 1);
	while (node = (Node *) getListIter(&iter)) {
	}
}

findOrthoRoot(Node *node, TreeNode **treeNode, Node *parent)
{
/*
	TreeNode *treeNode;
*/
	TreeNode *retnode;
	Domain *dom;
	int type;
	Node *save_root;

/*
printf("findOrthoRoot:");printNode(node);printf("\n");
printf("t>%d\n",*treeNode);
*/

/*
	if (IS_OUTGROUP_MODE(Opt.outgroupMode)) {
		int outchk = outgroupCheck(node);
		if (outchk == Both_Out || outchk == Out_Both) {
printf("root:::"); printNode(node); printf("\n");
			rootnode = node;
		}
	}
*/

	/* hom root */
	if (nodeMarked(node)) {
		/* skip duplicate: already visited */
	} else if (nodeMerged(node)) {
		/* merged node: ignore */
/*
printf("Merged:");printNode(node);printf("\n");
*/
	} else if ( isRootRaw(node) || checkCenterRoot(node, parent) ) {
		save_root = rootnode;
		rootnode = node;
		markNode(node);
/*
printf("rootnode:");printNode(rootnode);printf(":%d,%d\n",isRootRaw(node),checkCenterRoot(node));
*/
		type = convClustTree_sub(node, &retnode);
/*
printf("ROOT>");printNode(node);printf(";%d\n",retnode);
*/

/*
if (retnode->node->id==3075055){
printNode(node);
printNode(retnode);
printf("%d,%d\n",type,retnode->node->id);
}
*/
		if (retnode) {
			if (node != retnode->node){
/*
printf("type=%d\n",type);
printNode(node);putchar(' ');
printNode(retnode->node);putchar('\n');
*/
				if (type == 1 || type == 2) {
					/* flanking node: create a new node as a parent of retnode and add it to the root list  */
/*
printf("treeNode:1:"); printNode(node); printf("\n");
*/
					*treeNode = createTreeNode(node);
					ClustTree_addChild(*treeNode, retnode, type);
					pushList(root_treeNodes, *treeNode); 
				} else {
					/* one of the child nodes is deleted */
					*treeNode = retnode;
/*
printf("treeNode:1.1:"); printNode(node); printf("\n");
*/
					if (outgroupFlagCheck(retnode->node->spflag) != 2) {
 						/* set retnode as a root unless it contains only outgroup species */
/*
						makeRootNode_Out(retnode->node);
*/
						makeRootTreeNode_Out(retnode);
						resetNodeDomainAll(retnode->node);
						pushList(root_treeNodes, retnode); 
					}
					deleteRootNode(node);
/*
printf("makeNewRoot>>");printTreeNode(retnode);printf("\n");
print_specFlag(retnode->node->spflag);
printf("outgroupFlagCheck=%d\n",outgroupFlagCheck(retnode->spflag));
*/
				}
			} else {
				/* retnode=current_node: node is a bifurcating node: just add retnode to the root list */
				*treeNode = retnode;
/*
printf("treeNode:1.2:"); printNode(node); printf("\n");
*/
				pushList(root_treeNodes, retnode); 
			}
		} else {
/*
printf("delete_rootnode:"); printNode(node); printf("\n");
*/
			/** deleted cluster **/
			deleteRootNode(node);
		}
		rootnode = save_root;
	} else {
		if ( isOuterRoot(node) ) {
/*
if (rootnode != NULL) {
printf("ROOT IS NOT NULL!!!! ");
printNode(rootnode); printf("\n");
}
*/
			rootnode = node;
/*
printf("NEW ROOT IS "); printNode(rootnode); printf("\n");
*/
		}

		if (isFlankNode1(node)) {
			findOrthoRoot(node->child->node1, treeNode, node);
/*
printf("01 treenode=>%d\n",*treeNode);
*/
		} else if (isFlankNode2(node)) {
			findOrthoRoot(node->child->node2, treeNode, node);
/*
printf("02 treenode=>%d\n",*treeNode);
*/
		} else if (isIntNode(node)) {
			TreeNode *tnode_child1 = NULL;
			TreeNode *tnode_child2 = NULL;
	
			findOrthoRoot(node->child->node1, &tnode_child1, node);
			findOrthoRoot(node->child->node2, &tnode_child2, node);

			if (tnode_child1 && tnode_child2) {
				*treeNode = createTreeNodeUnroot(node);
/*
printf("treeNode:2:"); printNode(node); printf("\n");
*/
				ClustTree_addChild(*treeNode, tnode_child1, 1);
				ClustTree_addChild(*treeNode, tnode_child2, 2);
			} else if (tnode_child1) {
				*treeNode = tnode_child1;
			} else if (tnode_child2) {
				*treeNode = tnode_child2;
			} else {
				*treeNode = NULL;
			}
		} else {
			/** unclassified leaf node (probably outer group node) **/
			if (rootnode && node->domains) {
/*
printNode(node);
printf("; root=");printNode(rootnode);printf("\n");
printf(">>DOMTEST:%d\n", numelemList(node->domains));
*/
				getDomainMark(node->domains, rootnode, &dom);
				if (dom) {
					*treeNode = createTreeNodeUnroot(node);
/* 
printf("treenode:3:"); printNode(node); printf("\n");
*/
					(*treeNode)->dom = dom;
/* do not add as root
printf("DOMAIN::"); printDomain(dom); printf("\n");
					pushList(root_treeNodes, *treeNode); 
*/
				} else {
				/* domain is not defined */
					*treeNode = NULL;
/*
printf("nodom\n");
printf("root=%d\n",rootnode);
if (rootnode) {
	printNode(node);printf(" ");printNode(rootnode);printf("\n");
}
*/
				}
			} else {
				*treeNode = NULL;
/*
				*treeNode = createTreeNode(node);
printf("no_domains ");
printNode(node); printf("\n");
*/
			}
/*
printf("2treenode=>%d\n",*treeNode);
*/
		}
		if ( isOuterRoot(node) ) {
			rootnode = NULL;
		}
/*
if (*treeNode == NULL) {
printf("NULL!!\n");
} else {
printf("NOT NULL!!\n");
}
*/
	}

}
checkCenterRoot(Node *node, Node *parent)
{
	return (isCenterRoot(node) && parent == node->parent);
}

/*
 * converting an original Node strcture to a new TreeNode structure (a strictly bifurcating tree)
 * node: original Node
 * retnode: new TreeNode that should be connected to the parent of the current node
 */
convClustTree_sub(Node *node, TreeNode **retnode)
{
	TreeNode *treeNode;
	TreeNode *retnode1, *retnode2;
	Domain *dom;

	if (isFlankNode1(node)) {
		/* skip a flanking node in TreeNode structure */
		convClustTree_sub(node->child->node1, &retnode1);
		*retnode = retnode1;
		return 1;
	} else if (isFlankNode2(node)) {
		convClustTree_sub(node->child->node2, &retnode2);
		*retnode = retnode2;
		return 2;
	} else if (isIntNode(node)) {
		convClustTree_sub(node->child->node1, &retnode1);
		convClustTree_sub(node->child->node2, &retnode2);
		if (retnode1 && retnode2) {
			/* bifurcated node: create a new TreeNode */
			treeNode = createTreeNodeUnroot(node);
/*
printf("treenode:4:"); printNode(node); printf("\n");
*/
			ClustTree_addChild(treeNode, retnode1, 1);
			ClustTree_addChild(treeNode, retnode2, 2);
			*retnode = treeNode;
		} else if (! retnode1 && retnode2) {
			*retnode = retnode2;
		} else if (retnode1 && ! retnode2) {
			*retnode = retnode1;
		} else {
			/* both children are deleted */
			*retnode = NULL;
		}
		return 3;
	} else {
		/** leaf **/
		if (node->domains) {
/*
if (isRoot(rootnode)){
 printf("ROOT: "); printNode(rootnode); printf("\n");
}
if (isOuterRoot(rootnode)) {
 printf("OUTER ROOT\n");
}
*/
			getDomainMark(node->domains, rootnode, &dom);

			if (dom) {
/*
printf("DOM:"); printDomain(dom); printf("\n");
*/
				treeNode = createTreeNodeUnroot(node);
				treeNode->dom = dom;
				*retnode = treeNode;
			} else {
				/* domain not found: the leaf is deleted */
/*
printf("DELETED:"); printNode(node); printf(" "); printNode(rootnode); printf("\n");
*/
				*retnode = NULL;
			}
		}
		return 0;
	}
}
assignSpecFlag(pList *rootnodes)
{
	listIter iter;
	TreeNode *root;
	setListIter(&iter, rootnodes, 1);
	while (root = (TreeNode*) getListIter(&iter)) {
		assignSpecFlag_sub(root);
	}
}
assignSpecFlag_sub(TreeNode *tnode)
{
	specFlagP spflag1, spflag2;
	specFlag nullSpecFlag;

	/* initialize spflags */
	clearSPflag(nullSpecFlag);
	spflag1 = nullSpecFlag;
	spflag2 = nullSpecFlag;

	if (! tnode->child1 && ! tnode->child2) {
		copySPFlag(tnode->node->spflag, tnode->spflag);
		return(0);
	}
	if (tnode->child1) {
		assignSpecFlag_sub(tnode->child1);
		spflag1 = tnode->child1->spflag;
	}
	if (tnode->child2) {
		assignSpecFlag_sub(tnode->child2);
		spflag2 = tnode->child2->spflag;
	}
	spFlagOR(spflag1, spflag2, tnode->spflag);
/*
printf("assignSpecflag:");
printTreeNode(tnode);printf("%d,%d\n",tnode->flag,tnode->node->flag);
print_specFlag(tnode->spflag);
*/
}

printClustTree(TreeNode *tnode)
{
	if (ClustTree_isLeaf(tnode)) {
		printDomain(tnode->dom);
	} else {
		printClustTree(tnode->child1);
		printClustTree(tnode->child2);
	}
}

getChilds_TreeNode(TreeNode *tnode, TreeNode **child1, TreeNode **child2)
{
	*child1 = tnode->child1;
	*child2 = tnode->child2;
}
specFlagP getSpFlag_TreeNode(TreeNode *tnode)
{
	return(tnode->spflag);
}
void printTreeNode(TreeNode *tnode)
{
	printNode(tnode->node);
}
TreeNode *dupTreeNode(NodeSet *nodes, TreeNode *tnode)
{
	TreeNode *new_tnode;
	Node *newnode = dupNode(nodes, tnode->node);
	new_tnode = createTreeNode(newnode);
	copySPFlag(tnode->spflag, new_tnode->spflag);
	return new_tnode;
}
setFlagTreeNode(TreeNode *tnode, NodeFlag flag)
{
	tnode->flag |= flag;
}
unsetFlagTreeNode(TreeNode *tnode, NodeFlag flag)
{
	tnode->flag &= ~flag;
}
makeRootTreeNodeIn(TreeNode *tnode)
{
	unsetFlagTreeNode(tnode, NODE_DELETED);
	unsetFlagTreeNode(tnode, NODE_CLEARED);
	setFlagTreeNode(tnode, NODE_OUTGRP);
}
makeRootTreeNode_Out(TreeNode *tnode)
{
	unsetFlagTreeNode(tnode, NODE_DELETED);
	unsetFlagTreeNode(tnode, NODE_OUTGRP);
	unsetFlagTreeNode(tnode, NODE_CLEARED);
}


isTreeRootRaw(TreeNode *tnode)
{
	if (tnode==NULL) return 0;
	return (isRootRaw_flag(tnode->flag));
}
isTreeRoot(TreeNode *tnode)
{
	if (tnode==NULL) return 0;
	return (isRoot_flag(tnode->flag));
}
/* root or center-root (descendant of center parent)*/
checkTreeRoot(TreeNode *tnode)
{
/*
if (! isTreeRootRaw(tnode) && isCenterRoot_flag(tnode->flag)) {
	if (tnode->parent) {
		printf("checkCenter>"); printTreeNode(tnode);printf("%d,%d\n",tnode->node->parent, tnode->parent->node);
	} else {
		printf("NULL:: tnode->parent=%d\n",tnode->parent);
	}
}
*/
/*
if (tnode->parent) {
printf("checkCenterRoot>>");printTreeNode(tnode);printf("%d,%d,%d,%d,%d\n",isCenterRoot(tnode->node),tnode->parent->node,tnode->node->parent,tnode->node->parentL,tnode->node->parentR);
printTreeNode(tnode);printTreeNode(tnode->parent);printNode(tnode->node->parent);
printf("\n");
}
*/
/*
	return (isTreeRootRaw(tnode) ||
		(isCenterRoot_flag(tnode->flag) && tnode->parent && (tnode->node->parent == tnode->parent->node)));
*/
	return (isTreeRootRaw(tnode) || isCenterRoot_flag(tnode->flag));
}
isOutTreeRoot(TreeNode *tnode)
{
	return(isTreeRoot(tnode) || isOuterTreeRoot(tnode));
}
isOuterTreeRoot(TreeNode *tnode)
{
	if (tnode==NULL) return 0;
	return (isOuterRoot_flag(tnode->flag));
}
isAllTreeRoot(TreeNode *tnode)
{
	if (tnode==NULL) return 0;
	return (isAllRoot_flag(tnode->flag));
}
isInTreeRoot(TreeNode *tnode)
{
        if (tnode==NULL) return 0;
        return (isAllTreeRoot(tnode) && ! isTreeRoot(tnode));
}
isAllInTreeRoot(TreeNode *tnode)
{
        /* in-root under an out-root plus independent in-root (=out-root) */

        if (tnode==NULL) return 0;

        if (isAllTreeRoot(tnode)) {
                if (! ClustTree_isLeaf(tnode)) {
/*
if (tnode->child1) {
printf("c1=%d:",nodeVisited_flag(tnode->child1->flag));
}
if (tnode->child2) {
printf("c2=%d:",nodeVisited_flag(tnode->child2->flag));
}
printf("\n");
*/
                        if ( (! tnode->child1 || ! nodeVisited_flag(tnode->child1->flag)) ||
                                (! tnode->child2 || ! nodeVisited_flag(tnode->child2->flag)) ) {
                                return 1;
                        }
                } else {
                        return 1;
                }
        }
        return 0;
}


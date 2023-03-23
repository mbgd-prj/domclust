/*
 * DomClust: Hierarchical Clustering for Orthologous Domain Classification
 * Copyright (c) 2000-2007, Ikuo Uchiyama
 * All rights reserved.
 */
#include <stdio.h>
#include <strings.h>
#include "domclust.h"
#include "util.h"

static Alloc_Object *domobj;
static int checkNodeDomains_sub(Node *node, SeqPos from, SeqPos to,
		SeqPos savefrom, SeqPos saveto,
		Region *anchor, SeqPos seqlen);
static int _getDomainMark(pList *domains, Node *rootnode, Domain **dom, char mark);
Domain *createDomain();

createDomainObj()
{
	domobj = init_alloc_object(sizeof(Domain), NODE_BLKSIZ);
}
Domain *createDomain(Node *node)
{
	Domain *dom;
	dom = (Domain *) memalloc(domobj);
	dom->root = NULL;
	dom->leaf = node;
	dom->from = 1;
	dom->to = node->len;
	dom->mark = 0;
	dom->subgrp = 0;
	dom->num = 0;
/*
  printf("createDom:");printDomain(dom);printf(":%d\n",dom);
*/
	return dom;
}
deleteDomain(Domain *deldom)
{
	Node *leaf = deldom->leaf;
	Domain *dom;
	listIter iter;
	if (leaf==NULL) return(0);
	setListIter(&iter, leaf->domains, 1);
	while (dom = (Domain *) getListIter(&iter)) {
		if (dom == deldom) {
			delCurrElem(&iter);
		}
	}
}
getDomain(pList *domains, Node *rootnode, Domain **dom)
{
	return getDomainMark(domains, rootnode, dom);
}
getDomainNoMark(pList *domains, Node *rootnode, Domain **dom)
{
	return _getDomainMark(domains, rootnode, dom, 0);
}
getDomainMark(pList *domains, Node *rootnode, Domain **dom)
{
	return _getDomainMark(domains, rootnode, dom, 1);
}
_getDomainMark(pList *domains, Node *rootnode, Domain **dom, char mark)
{
	listIter iter;
	Domain *dm;
	int domn = 0, i = 1;
	if (domains == NULL) {
		return 0;
	}
	*dom = NULL;
	setListIter(&iter, domains, 1);
/*
printList(domains, printDomain);
*/
	while (dm = (Domain *) getListIter(&iter)) {
/*
if (dm->root==rootnode){
}
printNode(dm->root); printNode(rootnode); printf("::%d,%d\n",mark,dm->mark);
*/
		if (! dm->num) {
			dm->num = i;
		}
		if (dm->root == rootnode && dm->mark == 0) {
			*dom = dm;
			dm->mark = mark;
			domn = i;
/*
printf("#@#");printDomain(dm);printf(" ");printNode(dm->root); printNode(rootnode); printf("::%d,%d;%d\n",mark,dm->mark,dm);
*/
			break;
		}
		i++;
	}
	/** domn might be 0 when the domain has been discarded **/

	return domn;
}
clearDomainMark(pList *domains)
{
	listIter iter;
	Domain *dm;
	if (domains == NULL) {
		return 0;
	}
	setListIter(&iter, domains, 1);
	while (dm = (Domain *) getListIter(&iter)) {
		dm->mark = 0;
	}
}
clearDomainMarkAll(NodeSet *nodes)
{
	int i;
	Node *node;
	for (i = 0; i < nodes->leafnum; i++) {
		node = getNode(nodes, i);
		clearDomainMark(node->domains);
	}
}

static Node *leafnode; 
static int domid;
static int save_from;

int tmpflag;

resetNodeDomainAll(Node *node)
{
	if (isLeaf(node)) {
		checkNodeDomains(node);
	} else if (isFlankNode1(node)) {
		resetNodeDomainAll(node->child->node1);
	} else if (isFlankNode2(node)) {
		resetNodeDomainAll(node->child->node2);
	} else {
		resetNodeDomainAll(node->child->node1);
		resetNodeDomainAll(node->child->node2);
	}
}

/* For each sequence, determine the domain boundaries */
checkDomains(NodeSet *nodes)
{
	Node *node;
	int i;

	if (Opt.DEBUG & DBG_domain) {
		printf("## checkDomains\n");
	}
	createDomainObj();
	for (i = 0; i < nodes->leafnum; i++) {
		node = getNode(nodes, i);
		if (Opt.DEBUG_ent) {
			if (strcmp(node->name, Opt.DEBUG_ent) == 0) {
				Opt.DEBUG = Opt.DEBUG_val;
			} else {
				Opt.DEBUG = 0;
			}
		}
		checkNodeDomains(node);

/*
		checkNodeDomains(node, (SeqPos) 1, node->len,
				(SeqPos) 1, node->len, NULL, node->len);
*/
	}
	Opt.DEBUG = Opt.DEBUG_val;
}

checkNodeDomains(Node *node)
{
	leafnode = node;
	domid = 1;
	if (leafnode->domains) {
		clearList(leafnode->domains);
	} else {
		leafnode->domains = (pList *) create_pList();
	}
	save_from = 0;
/*
printf("Leaf>");printNode(node);printf("\n");
*/
	checkNodeDomains_sub(node, (SeqPos) 1, node->len, (SeqPos) 1, node->len, NULL, node->len);

/*
    if (strcmp(node->name, "hpy:HP0286")==0) {
    if (strcmp(node->name, "cje:CJ1005C")==0) {
    if (strcmp(node->name, "hpy:HP1521")==0) {
    if (strcmp(node->name, "hpy:HP0289")==0) {
    if (strcmp(node->name, "sce:MBGD_120330")==0) {
    if (strcmp(node->name, "sce:YNR051C")==0) {
    if (strcmp(node->name, "ape:APE1507")==0) {
    if (strcmp(node->name, "hpy:HP1521")==0) {
    if (strcmp(node->name, "hpy:HP1512")==0) {
    if (strcmp(node->name, "hpy:HP0527")==0) {
    if (strcmp(node->name, "sce:MBGD_120453")==0) {
*/
	if (Opt.DEBUG & DBG_domain) {
/*
		if (Opt.DEBUG_ent == NULL || strcmp(node->name, Opt.DEBUG_ent) == 0) {
*/
			listIter iter;
			Domain *dom;
			setListIter(&iter, node->domains, 1);
			while (dom = (Domain *) getListIter(&iter)) {
				printDomain(dom); printf(" root=");
				printNode(dom->root);
				printf(",%p; ",dom);
				printf("\n");
			}
/*
		}
*/
    	}
}

/* from, to: anchor points on on the current sequence */
/* savefrom, saveto: saved domain boundary on the current sequence */
checkNodeDomains_sub(Node *node, SeqPos from, SeqPos to,
		SeqPos savefrom, SeqPos saveto, Region *anchor, SeqPos seqlen)
{
	char flag = 0;
	SeqPos from1, to1;
	SeqPos savefrom1, saveto1;
	SeqPos len = to - from + 1;
	Domain *dom;
	Edge *e, *child;
	SeqPos alilen;
	Region *reg, thisNode, currReg;
	int parentCheck = 1;

	/* currReg: coordinate on the examined sequence */
	/* thisNode: coordinate on this node */
	/*     they are aligned to each other */
	if (anchor != NULL) {
		copyReg(&thisNode, anchor);
	} else {
		thisNode.from = 1; thisNode.to = node->len;
	}
	currReg.from = from; currReg.to = to;
	
/*
if (node && node->id==110799) {
printf("Dom==");printNode(node);printf("\n");
printf("is_outroot=%d\n", isOutRoot(node));
}
printf("dom>>");printNode(node);printf("%d,%d,%d; %d,%d\n",isRoot(node),isOuterRoot(node),isCenterRoot(node), savefrom, saveto);
*/

	if (nodeMerged(node)) {
		if (Opt.DEBUG & DBG_domain) {
			printf("merged>>");printNode(node);printf("%d,%d;%d,%d\n",isRoot(node),isOuterRoot(node),savefrom,saveto);
		}

		if (getMergeDir(node) == 1) {
			/* merge this node into the right (next) node */
			if (save_from == 0) {
				save_from = savefrom;
			}
		} else {
			/* merge this node into the left (previous) node */
			dom = (Domain *) getListIdx(leafnode->domains, -1);
			if (dom) {
				/** merged region **/
				dom->to = saveto;
			}
		}
		parentCheck = 0;
	} else if (isOutRoot(node)) {
		/* root or outer-root */

		if (Opt.DEBUG & DBG_domain) {
			printf("ROOT>>");printNode(node);printf("%d,%d,%d; %d,%d\n",isRoot(node),isOuterRoot(node),isCenterRoot(node), savefrom, saveto);
		}

	    if (isCenterRoot(node)) {
			/* root is placed on the center domain: skip ascending the center parent */
			parentCheck = 2;
/*
printf("ccc:%d,%d,%d\n",node->parentL,node->parent,node->parentR);
*/
	    } else {
		int tmplen;
/*
	if (isRoot(node)) {
*/

		if (Opt.outstyle == DOMAINOUT) {
			printf("%s(%d) %d %d %d %d\n",
				leafnode->name, domid++,
				node->id, from, to, node->len);
		}
		if (save_from > 0) {
			savefrom = save_from;
		}
		savefrom = max(savefrom, 1);
		saveto   = min(saveto, seqlen);
		tmplen = saveto - savefrom + 1;
			
		if (tmplen < Opt.min_minlen && tmplen < seqlen * Opt.ovlpratio2) {
			/** discard short domains **/
/*
printf(">short: "); printNode(node); printf("\n");
*/
			return 0;
		}

		dom = createDomain(node);

		if (save_from > 0) {
			save_from = 0;
		}
		dom->from = savefrom;
		dom->to = saveto;
		dom->root = node;
		dom->leaf = leafnode;
		pushList(leafnode->domains,dom);
/*
printf(">dom>>>");printDomain(dom);printf(":%d,%d\n",savefrom,saveto);
*/

		parentCheck = 0;
	    }

/*
	} else if (nodeDeleted0(node)) {
		/ * node that have been already cut--do nothing * /
		printf(">deleted\n");
		parentCheck == 0;
*/
	}
/*
printf("brk=(%d,%d)\n",node->brk.from,node->brk.to);
if (! definedPos(node->brk.from) && ! definedPos(node->brk.to) && node->parent && (node->parentL || node->parentR)) {
printf("????????????:"); printNode(node); printf("\n");
}
*/

	if (parentCheck) {
		/***** ------ParentL------- *****/
/*
printf("%d,%d,%d;(%d,%d)\n",node->parentL,node->parent,node->parentR,node->brk.from,node->brk.to);
		if (node->parentL) {
*/
		if (nodeDefinedParentL(node)) {
			if (definedPos(node->brk.from)) {
				saveto1 = transformSeqPos(node->brk.from-1, &thisNode, &currReg);
			} else if (definedPos(node->brk.to)) {
				saveto1 = transformSeqPos(node->brk.to, &thisNode, &currReg);
			} else{
				saveto1 = saveto;
			}
			from1 = from;
			to1 = transformSeqPos(node->newreg.from-1, &thisNode, &currReg);

			if (Opt.DEBUG & DBG_domain) {
				printf(">Left>>>");printNode(node);printf("; %d,%d\n",savefrom1,saveto1);
			}

#ifdef WITH_DEBUG
/*
printf(">>%d,%d  %d,%d,%d,%d  (%d,%d),(%d,%d)<\n",from1,to1,from,to,len,node->len,node->brk.from,node->brk.to, node->newreg.from,node->newreg.to);
printf("## L %d,%d <- %d,%d,%d,%d %s\n",from1,to1,from,to,node->brk.from,node->brk.to, leafnode->name);
*/
#endif
			checkNodeDomains_sub(node->parentL, from1, to1,
				savefrom, saveto1, NULL, seqlen);
			flag = 1;
		}

		/***** ------ParentC------- *****/
/*
		if (node->parent) {
*/
		if (nodeDefinedParentC(node)) {
			if (Opt.domBoundary == B_BREAK) {
				/* break point */
				reg = &(node->brk);
			} else {
				/* aligned region */
				reg = &(node->newreg);
			}
			if (definedPos(node->brk.from)) {
				savefrom1 = transformSeqPos(node->brk.from, &thisNode, &currReg);
			} else {
				savefrom1 = savefrom;
			}

			if (definedPos(node->brk.to)) {
				saveto1 = transformSeqPos(node->brk.to, &thisNode, &currReg);
			} else {
				saveto1 = saveto;
			}

			/* mapping: node->newreg <-> parent->consreg */
			from1 = transformSeqPos(node->newreg.from, &thisNode, &currReg);
			to1 = transformSeqPos(node->newreg.to, &thisNode, &currReg);

			if (parentCheck == 2) {
				if (Opt.DEBUG & DBG_domain) {
					printf(">CenterRoot>>>");printNode(node);printf("; %d,%d\n",savefrom1,saveto1);
				}
				outputDomain(node, from1, to1, savefrom1, saveto1, seqlen);
			} else {
				if (Opt.DEBUG & DBG_domain) {
					printf(">Center>>>");printNode(node);printf("; %d,%d\n",savefrom1,saveto1);
				}

				checkNodeDomains_sub(node->parent, from1, to1,
					savefrom1, saveto1,
					&(node->parent->consreg), seqlen);


			    if (definedPos(node->newreg2.from)
					&& definedPos(node->newreg2.to)) {
				int from2, to2;
				Region *reg2;

				if (Opt.domBoundary == 1) {
					/* break point */
					reg2 = &(node->brk2);
				} else 	{
					/* aligned region */
					reg2 = &(node->newreg2);
				}
				
				from2 = transformSeqPos(node->newreg2.from, &thisNode, &currReg);
				to2 = transformSeqPos(node->newreg2.to, &thisNode, &currReg);

if (node->parentM){
				checkNodeDomains_sub(node->parentM, to1+1, from2-1,
					to1+1, from2-1,
					NULL, seqlen);
}
				checkNodeDomains_sub(node->parent, from2, to2,
					from2, to2,
					&(node->parent->consreg), seqlen);
			    }
			}
#ifdef WITH_DEBUG
/*
printf("## > %d,%d %s,%s\n", alilen, len, leafnode->name, node->name);
printf("## C %d,%d <- %d,%d,%d,%d %d/%d %s\n",from1,to1,from,to,node->brk.from,node->brk.to, node->len, alilen, leafnode->name);
*/
#endif
			flag = 1;
		}

		/***** ------ParentR------- *****/
/*
		if (node->parentR) {
*/
		if (nodeDefinedParentR(node)) {
			if (definedPos(node->brk.to)) {
				savefrom1 = transformSeqPos(node->brk.to+1,&thisNode,&currReg);
			} else {
				savefrom1 = savefrom;
			}
			from1 = transformSeqPos(node->newreg.to+1,&thisNode,&currReg);
			to1 = to;

			if (Opt.DEBUG & DBG_domain) {
				printf(">Right>>>");printNode(node);printf("; %d,%d\n",savefrom1,saveto1);
			}

#ifdef WITH_DEBUG
/*
printf("## R %d,%d <- %d,%d,%d,%d %s; ",from1,to1,from,to,node->brk.from,node->brk.to, leafnode->name);
printf(">>%d,%d\n",node->brk.from,node->brk.to);
printNode(node);printf("\n");
*/
#endif
			checkNodeDomains_sub(node->parentR, from1, to1,
					savefrom1, saveto, NULL, seqlen);
			flag = 1;
		}
/*
		if (! flag) {
			/ * never come here <<== NOT TRUE when OuterRoot is considered * /
			fprintf(stderr, "Warning: No parent at unrooted node: %s (%d)\n",
				node->name,node->id);
		}
*/
	}
}

outputDomain(Node *node, SeqPos from, SeqPos to, SeqPos savefrom, SeqPos saveto, SeqPos seqlen)
{
	int tmplen;
	Domain *dom;
/*
	printf("root>>");printNode(node);printf("%d,%d; \n",isRoot(node),isOuterRoot(node), savefrom, saveto);
*/
	if (Opt.outstyle == DOMAINOUT) {
		printf("%s(%d) %d %d %d %d\n",
			leafnode->name, domid++,
			node->id, from, to, node->len);
	}
	if (save_from > 0) {
		savefrom = save_from;
	}
	savefrom = max(savefrom, 1);
	saveto   = min(saveto, seqlen);
	tmplen = saveto - savefrom + 1;
		
	if (tmplen < Opt.min_minlen && tmplen < seqlen * Opt.ovlpratio2) {
		/** discard short domains **/
/*
printf(">short: "); printNode(node); printf("\n");
*/
		return 0;
	}

	dom = createDomain(node);

	if (save_from > 0) {
		save_from = 0;
	}
	dom->from = savefrom;
	dom->to = saveto;
	dom->root = node;
	dom->leaf = leafnode;
	pushList(leafnode->domains,dom);
/*
printf(">dom>>>");printDomain(dom);printf(":%d,%d\n",savefrom,saveto);
*/
}

void printDomain(Domain *dom)
{
	printNode(dom->leaf);
	printf(" (%d) %d %d",dom->num, dom->from, dom->to);
}

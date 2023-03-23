#include <stdio.h>
#include "domclust.h"
#include "sstack.h"
#include "hash.h"
#include "namehash.h"

#define STACKSIZE 5000

static sStack *taxNodes;
typedef struct {
/*
	int taxnode_id;
*/
	specFlag spFlag_cmpr;
	Node* genenode;
} TaxNodeData;
static Node *root;

/** buffer **/
typedef struct {
	TreeNode *tnode;
	int domn;
	TaxNodeData taxNode1;
	TaxNodeData taxNode2;
	char flag;
} TaxMapOutBuf;
static varArray *taxMapOutData;

mapTaxInfoAll(ClusterInfoData *cinfoData)
{
	listIter iter;
	int i;
	FILE *fp;
	ClusterInfo *cinfo = cinfoData->cinfo;
	int clustnum = cinfoData->clustnum;

	fp = Opt.taxMapOutFp;
	if (fp == NULL) fp = stdout;
	for (i = 0; i < clustnum; i++) {
		if (cinfo[i].clustid < 0) {
			continue;
		}
		fprintf(fp, "#clustid=%d\n",cinfo[i].clustid);
		mapTaxInfo_Cinfo(&cinfo[i]);
	}
	free_sStack(taxNodes);
}
mapTaxInfo_Cinfo(ClusterInfo *cinfo)
{
	listIter iter;
	setListIter(&iter, cinfo->clusters, 1);
	TreeNode *tnode;
	int i;
	arrayIter arrayIter;
	TaxMapOutBuf *taxinfo;
	Hash *taxMapHash;
	NameHash *nhash;
	HENTRY ent;
	int arysize;
	char strbuf[512];
	
	FILE *fp = Opt.taxMapOutFp;

	taxMapOutData = createVarArray(1000, sizeof(TaxMapOutBuf));
	while (tnode = (TreeNode *) getListIter(&iter)) {
		root = tnode->node;
		mapTaxInfo(tnode);
	}
	arysize = arraySize(taxMapOutData);
	if (arysize == 0) {
		return(0);
	}
	
	arysize = (arysize < 31) ? 31 : arysize;
	
	nhash = initNames(arysize + 3);
	taxMapHash = Hcreate(arysize * 7 - 6);
	setArrayIter(&arrayIter, taxMapOutData);
	/* eliminate redundant data */
	while (taxinfo = (TaxMapOutBuf *) getArrayIter(&arrayIter)) {
		int score = taxinfo_score(taxinfo);
		getNameFromTaxInfo(taxinfo, strbuf);
		ent.key = addName(nhash, strbuf, 1);
		ent.datum = (char*)taxinfo;
		if (Hsearch(taxMapHash, &ent, ENTER)) {
			TaxMapOutBuf *tmpbuf = (TaxMapOutBuf*) ent.datum;
			if (tmpbuf != taxinfo) {
				int maxScore = taxinfo_score(tmpbuf);
				if (score > maxScore) {
					memcpy(tmpbuf, taxinfo, sizeof(TaxMapOutBuf)); /* copy taxinfo to hashed data */
				}
				/* delete */
				taxinfo->flag = 0;
			}
		}
	}
	setArrayIter(&arrayIter, taxMapOutData);
	while (taxinfo = (TaxMapOutBuf *) getArrayIter(&arrayIter)) {
		if (taxinfo->flag == 0) continue;
		outputTaxInfo(fp,  taxinfo);
	}
	freeArray(taxMapOutData);
	Hdestroy(taxMapHash);
	freeNames(nhash);
}
mapTaxInfo(TreeNode *tnode)
{
	if (taxNodes == NULL) {
		taxNodes = create_sStack(STACKSIZE, sizeof(TaxNodeData));
	}
	mapTaxInfo_sub(tnode);
	reset_sStack(taxNodes);		
}
mapTaxInfo_sub(TreeNode *tnode)
{
	Node *node = tnode->node;
	TreeNode *child1 = tnode->child1;
	TreeNode *child2 = tnode->child2;
	int skipflag1 = 0, skipflag2 = 0;
	TaxNodeData spnode;
/*
	specFlag spFlag_cmpr, spFlag_taxQueryCompl;;
*/

/*
	if (taxQueryCluster(tnode->spflag)
*/
	if ( (taxQueryCluster(tnode->spflag) && unknownOnlyCluster(tnode->spflag))
			|| spFlagCntAll(tnode->spflag)==1 ) {
		/* output result */
		outputTaxMap(tnode, taxNodes, Opt.taxMapOutFp);
		return(0);
	} else if (taxKnownCluster(tnode->spflag)) {
		/* no query taxon: do nothing */
		return(0);
	}
	spFlagANDNOT(tnode->spflag, SPflags.unknown, spnode.spFlag_cmpr);
/*
	spnode.taxnode_id = sptree_MatchFlags(spFlag_cmpr);
*/
	spnode.genenode = node;
	push_sStack(taxNodes, &spnode);

	/* set skipflag if the other subtree contains only taxquery clusters
	 * --sp2
	 * |
         * |  +- unk
	 * +- +
         *    +- sp1
	 */

	if (! child1 || unknownOnlyCluster(child1->spflag)) {
		skipflag1 = 1;
	}
	if (! child2 || unknownOnlyCluster(child2->spflag)) {
		skipflag2 = 1;
	}

	if (skipflag2 && child1) {
		if (child2) mapTaxInfo_sub(child2);
		pop_sStack(taxNodes);
		mapTaxInfo_sub(child1);
	} else if (skipflag1 && child2) {
		if (child1) mapTaxInfo_sub(child1);
		pop_sStack(taxNodes);
		mapTaxInfo_sub(child2);
	} else {
		if (child1) mapTaxInfo_sub(child1);
		if (child2) mapTaxInfo_sub(child2);
		pop_sStack(taxNodes);
	}
}

outputTaxMap(TreeNode *tnode, sStack *taxNodes, FILE *fp)
{
	int *p;
	TaxNodeData *spnode, *prev_spnode;
	spnode = (TaxNodeData*) get_sStackIdx(taxNodes, -1);
	prev_spnode = (TaxNodeData*) get_sStackIdx(taxNodes, -2);

	if (fp == NULL) {
		fp = stdout;
	}
	outputTaxMap_sub(tnode, spnode, prev_spnode, fp);
}
outputTaxMap_sub(TreeNode *tnode, TaxNodeData *taxNode1, TaxNodeData *taxNode2, FILE *fp)
{
	Domain dom;
	int domn = 0;
	TaxMapOutBuf taxInfo;
	if (ClustTree_isLeaf(tnode) && taxQueryCluster(tnode->spflag)) {
/*
		domn = getDomainNoMark(tnode->node->domains,root,&dom);
		if (numelemList(tnode->node->domains) <= 1) {
			domn = 0;
		}
*/
		if (numelemList(tnode->node->domains) > 1) {
			domn = tnode->dom->num;
		}

		copyTaxNode(&(taxInfo.taxNode1), taxNode1);
		copyTaxNode(&(taxInfo.taxNode2), taxNode2);
		taxInfo.domn = domn;
		taxInfo.tnode = tnode;
		taxInfo.flag = 1;

		addArray(taxMapOutData, &taxInfo);
/*

		fprintf(fp, "%s", tnode->node->name);
		if (domn) {
			fprintf(fp, "(%d)", domn);
		}
		fprintf(fp, "\t");
printf("ff>>%d\n", tnode->spflag);
		taxnode_output(fp, taxNode1, tnode->spflag);
		fprintf(fp, "\t");
		taxnode_output(fp, taxNode2, tnode->spflag);
		fprintf(fp, "\n");
*/
	} else {
		if (tnode->child1) {
			outputTaxMap_sub(tnode->child1, taxNode1, taxNode2, fp);
		}
		if (tnode->child2) {
			outputTaxMap_sub(tnode->child2, taxNode1, taxNode2, fp);
		}
	}
}
copyTaxNode(TaxNodeData *taxData1, TaxNodeData *taxData2) {
	if (! taxData1) {
		fprintf(stderr, "taxdata1 is null\n");
	} else if (taxData2) {
		memcpy(taxData1, taxData2, sizeof(TaxNodeData));
	} else {
		taxData1->genenode = NULL;
	}
}
getNameFromTaxInfo(TaxMapOutBuf *taxInfo, char *strbuf) {
	sprintf(strbuf, "%s:%d", taxInfo->tnode->node->name, taxInfo->domn);
	
}
outputTaxInfo(FILE *fp,  TaxMapOutBuf *taxInfo) {
	if (fp == NULL) {
		fp = stdout;
	}
	fprintf(fp, "%s", taxInfo->tnode->node->name);
	if (taxInfo->domn) {
		fprintf(fp, "(%d)", taxInfo->domn);
	}
	fprintf(fp, "\t");

	taxnode_output(fp, &(taxInfo->taxNode1), taxInfo->tnode->spflag);
	fprintf(fp, "\t");
	taxnode_output(fp, &(taxInfo->taxNode2), taxInfo->tnode->spflag);
	fprintf(fp, "\n");
}
taxnode_output(FILE *fp, TaxNodeData *taxNode, specFlag selfFlag)
{
	char *name;
	specFlag spFlag_cmpr;
	static int taxnode_id;
	if (taxNode->genenode == NULL) {
		fprintf(fp, "null");
	} else {
		if (selfFlag) {
			spFlagANDNOT(taxNode->spFlag_cmpr, selfFlag, spFlag_cmpr);
			taxnode_id = sptree_MatchFlags(spFlag_cmpr);
		} else {
			taxnode_id = sptree_MatchFlags(taxNode->spFlag_cmpr);
		}
		name = getTaxName(taxnode_id);
		if (name == NULL) {
			fprintf(stderr, "Warning: taxonomhy mapping failed: taxid=%d\n",taxnode_id);
		}
		fprintf(fp, "%s(%.1f,%.1f,%d)", name,
			taxNode->genenode->child->score,
			taxNode->genenode->child->dist,
			spFlagCnt(taxNode->genenode->spflag)-1);
	}
}
taxinfo_score(TaxMapOutBuf *taxOut)
{
	return ( taxnode_score(&(taxOut->taxNode1)) + taxnode_score(&(taxOut->taxNode2)) );
}
taxnode_score(TaxNodeData *taxNode)
{
	int spcnt;
	if (taxNode->genenode==NULL) {
		return 0;
	}
	spcnt = spFlagCnt(taxNode->genenode->spflag)-1;
	return(spcnt);
}

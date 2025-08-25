/*
 * DomClust: Hierarchical Clustering for Orthologous Domain Classification
 * Copyright (c) 2000-2007, Ikuo Uchiyama
 * All rights reserved.
 */
#ifndef _DOMCLUST_
#define _DOMCLUST_
#include "memalloc.h"
#include "vararray.h"
#include "plist.h"
#include "bin.h"
#include "seqreg.h"
#include "namehash.h"
#include "spec.h"
#include<limits.h>

#define DOMCLUST_VERSION "1.2.8h"

#define MAXDOM 256
#define todigit(c) (c - '0')
#define BETTER(a,b) (Better((double) a,(double) b,Opt.sim))
#define DISTRATIO(a,b) ( DistRatio((double) a,(double) b, Opt.sim) )
#define MEASURE(e) (Opt.sim ? (e)->score : (e)->dist)

#define BADVALUE 99999999
#define DISTSCALE 40
#define MAXDIST 270
#define UNDEFINED 99999999
#define MISSRATIO 0.95
#define MINSCORE 60
#define MAXSCORE 50000
#define MAXUCHR 255
#define MAXNBR 2500
#define MINLEN_LIMIT 10

#define INIT_NODENUM 200000
#define NODE_BLKSIZ 40000
#define EDGE_BLKSIZ 120000
#define INT_BLKSIZ 100000
#define AVE_REL 50
#define ONE .99999999

#define NAMELEN 200
#define SPNAMELEN 40

#define otherNodeNoCheck(e,n) (e->node1->id == n ? e->node2 : e->node1)
#define otherNodeCheck(e,n) \
  ( e->node1->id == n ? e->node2 : ((e->node2->id == n) ? e->node1 : NULL) )
#define otherNode(e,n) otherNodeCheck(e,n)
#define otherNodeID(e,n) otherNodeNoCheck(e,n)->id
#define otherNodeAli(e,n) (e->node1->id == n ? e->ali2 : e->ali1);

#define NODE_CLEARED ((NodeFlag) 1)
#define NODE_INT1 ((NodeFlag) 2)
#define NODE_INT2 ((NodeFlag) 4)
#define NODE_INTERNAL (NODE_INT1|NODE_INT2)
#define NODE_DELETED ((NodeFlag) 8)
#define NODE_MERGED ((NodeFlag) 16)	/* the cluster is included in the adjacent cluster */
#define NODE_MERGED2 ((NodeFlag) 32) /* sequences once split are merged during the clustering */
#define NODE_TMPMARK ((NodeFlag) 128)
#define NODE_TMPMARK2 ((NodeFlag) 64)
#define NODE_DUPMARK ((NodeFlag) 256)
#define NODE_OUTGRP ((NodeFlag) 512)	/* marking outgroup nodes */
#define NODE_INGRP ((NodeFlag) 1024)	/* marking ingroup root nodes */
#define NODE_OUTERROOT ((NodeFlag) 2048)	/* marking outer root nodes */
#define NODE_OUTERROOT_DEL ((NodeFlag) 4096)	/* deleted outer root nodes */
#define NODE_CENTER_ROOT ((NodeFlag) 8192) /* root is placed only on the center domain and suppress ascending only center parent; the flag is needed only for center because root can be placed on the other domains using flanking node */

#define TREENODE_HASROOT ((NodeFlag) 1)

#define EDGE_CLEARED 1
#define EDGE_SELECTED 2
#define EDGE_SKIP 4
#define EDGE_DELETED 8
#define EDGE_TMPMARK 128
#define EDGE_TMPMARK2 64

#define IS_OUTGROUP_MODE(x)  (x == OutGroupMode || x == OldOutGroupMode)
#define IS_NEW_OUTGROUP_MODE(x)  (x == OutGroupMode)
#define IS_METAGENOME_MODE(x)  (x == MetaGenomeMode)

/*
#define WITH_NEIGHBOR
*/

enum outStyle {
	NORMALOUT, TREEOUT, NEWICK, NEWICK_DIST, HIER_DETAILOUT,
		CLUSTTAB, GRAPHOUT, NEIGHBOR, DOMAINOUT, TABOUT, DUMP,
		TREEDETAIL
};
typedef enum {NORMAL_EDGE, MULTI_EDGE, SELF_EDGE, SELF_MATCH} EdgeType;
typedef enum {B_BREAK, B_ALIREG} BoundaryOutput;
typedef enum {D_MATCH, D_LEFT1, D_LEFT2, D_RIGHT1, D_RIGHT2,
                D_MATCH1, D_MATCH2, D_THIRD} DomIdx;
typedef enum {DBG_basic=1, DBG_cluster=2, DBG_nbrconst=4, DBG_nbrrestore=8, DBG_domain=16}
	DBG_Mode;
typedef enum OutGrpStat {
	None = -1, In_In, Out_Out, Out_Both, Both_Out, Both_Both, In_In_Horiz,
} OutGrpStat;
typedef enum OutGrpModeType {
	Off, OutGroupMode, OldOutGroupMode, MetaGenomeMode,
} OutGrpModeType;
/* for distinguishing different level of SubClusterInfo */
typedef enum {
	SubCluster, OuterCluster, HomCluster,
} ClusterLevel;

typedef unsigned char DomID;
typedef unsigned short ConnCount;
typedef unsigned short Count;
typedef unsigned int NodeID;
typedef unsigned int EdgeID;
typedef unsigned int ClustID;

#define MAX_CONN USHRT_MAX

/*
typedef double Dist; 
*/
typedef float Dist; 
typedef unsigned short NodeFlag; 
typedef unsigned char EdgeFlag; 

/* SeqRegion */
/*
#define regLen(reg) ((reg)->to - (reg)->from + 1)
#define checkOvlpReg(reg1,reg2) \
        ((reg1)->from < (reg2)->to && (reg2)-> from < (reg1)->to)

typedef struct {
	SeqPos from, to;
} Region;
*/

/* Graph */
typedef struct Node {
	char *name;
	NodeID id;
	NodeFlag flag;
	signed char dir;
	pList *domains;
	SeqPos from, len;
	SeqPos meanlen;
	SeqPos minlen;
	SeqPos totlen;
	struct Node *parent, *parentL, *parentR;
	pList *left, *right;
	Region consreg;		/* conserved region within the seq */

	Region brk, newreg;

	struct Node *parentM;
	Region brk2, newreg2;	/* for self match */

	struct Edge *child;
	Count cnt;
	Dist meandist;
	specFlag spflag;
	char truncFlag;
} Node;
typedef struct {
	Alloc_Object *nodeobj, *intobj, *lenobj, *idobj;
	int leafnum;
	int total_nodenum;	/* leaf nodes + internal nodes */
	int nodenum;		/* extra nodes may be added */
} NodeSet;
typedef struct {
	NodeSet *nodes;
	int i;
	int dir;
	int end;
} NodeSetIter;

typedef struct Edge {
	Node *node1, *node2;
	StrRegion *ali1, *ali2;
	Dist dist, score;
	listElem *binelem;
	listElem *ordelem1, *ordelem2;
	EdgeID id;
#ifdef WITH_NEIGHBOR
	pList *left, *right;
#endif
	EdgeFlag flag;
#ifdef LARGE
	signed char dir;	/* reverse match on nucleotide comparison */
#endif
#ifdef EXPERIMENTAL
	ConnCount connect;	/* number of true connection */
#endif
} Edge;
typedef struct {
	Alloc_Object *edgeobj, *idobj, *lenobj, *regobj;
	unsigned long long edgenum;
	int *bestdist;
	Edge **bestedge;
	varArray *tmplist;
	pList *ordlist;
	int *nodeindex;
	Bin *bin;

/*
	int (*better)(double a,double b);
	char sim;
*/
} EdgeSet;

typedef enum {
	n3flag_skip = 1
} Node3Flag;
typedef struct {
	NodeID n3;
	Edge *e1, *e2;
	Region *ali13, *ali31, *ali23, *ali32;
	Region ali_e11, ali_e12, ali_e21, ali_e22;
	char flag;
} Node3List;

typedef struct {
	NodeID id1, id2;
	Edge *edge;
	listElem *ordelem;
} pairList;

typedef struct {
	NodeID id;
	ClustID clid;
} Cluster;

typedef struct {
	Node *root;
	Node *leaf;
	int num;
	SeqPos from;
	SeqPos to;
	char mark;
	int subgrp;
} Domain;

typedef struct TreeNode {
	Node *node;
	struct TreeNode *parent, *child1, *child2;
	Domain *dom;
	specFlag spflag;
	NodeFlag flag;
	NodeFlag treenodeFlag;
} TreeNode;

typedef struct {
	NodeSet *nodes;
	EdgeSet *edges;
	NameHash *nhash;
} SimGraph;

typedef struct {
	Region break1, break2, boundary1, boundary2;
	Region seq, aliM;
	Region consreg;
	SeqPos meanlen;
	SeqPos minlen;
/*
	SeqPos maxlen;
*/
} NewSeqInfo;
typedef struct {
	TreeNode *root;
	pList *members;
	int clustid;
	ClusterLevel level;
} SubClusterInfo;
typedef struct {
	TreeNode *root;
	TreeNode *origroot;
	TreeNode *outer_root;
	pList *clusters;
	pList *clusters_outer;
	pList *members;
	pList *ingroups;
	SubClusterInfo *outgroup;
	SubClusterInfo *outgroup2;
	specFlag spflag;
	int clustid;
	int outerclust;
	int homclust;
	Node *homclustRoot;
} ClusterInfo;
typedef struct {
	ClusterInfo *cinfo;
	ClusterInfo **cinfo_link;
	Hash *cinfoHash;
	int clustnum;
} ClusterInfoData;

SubClusterInfo *createSubClusterInfo(), *createSubClusterInfo_general(), *createSubClusterInfoN(), *createOuterClusterInfo();

typedef struct {
	Dist cutoff, cutoff2, in_cutoff, sumcut, cutoff2dist;
	int distscale;
	int logdistscale;
	Dist missdist, missscore;
	float missratio, cutoff2ratio;
	float cutoff2_score_per_len;
	DBG_Mode DEBUG;
	int  DEBUG_val;
	char *DEBUG_ent;
	int VERBOSE;
	int verbose_step1, verbose_step2;
	int minsp, minent;	/* minimum num. spececies/genes */
	int delete_small;	/* delete small grp. before dom. split */
	float ovlpratio, ovlpratio2;
	float coverage, covfilt;
	int mincutcnt;
	int minovlp, minovlp2;
	int minlen, minlen2, min_minlen;
	int min_alilen;
	enum outStyle outstyle;
	int sim;
	int nobreak;
	double phylocutratio, phylocutratio2;
	double distdiffcut;
	double treeBalanceRatio;
	int neighbor;
	int dpWin;
	int decinum;
/*
	double chkInclClst, chkOvlpClst;
*/
	double adjInclRatio, adjOvlpRatio;
	char chkSpFlg;
	BoundaryOutput domBoundary;
	int domcut;
	float chkConnect;
	Dist iniGap, extGap, skipGap;
	double nbrScoreRatio,nbrScoreDecay,nbrScoreLim;
	double nbrConnRatio;
	int nbrConnGap;
	int phyloMode;
	int skipErrEnt;
	int revMatch;
	char *dumpfile;
	int treecheck;
	int ignoreTrunc;
	int noextend; /*tmp*/

	char *outgroupStr;
	char *spmaskStr;
	char outgroupRev;	/* specifying ingroup */
	char outgroupMode;
	char metagenomeMode;
	double horizweight; 
	char *metaSpecPref;
	char *metaSpecList;
	char *partialSpecList;
	char *taxMapSpecList;
	char outputSubClustTree;
	char outputOuterClustTree;

	char homClustOut;
	char outputScore;
	FILE *outputScoreFp;
	char taxMapOut;
	FILE *taxMapOutFp;
	int noReplaceSpTreeLeafName;

	char *genes;
	char NEW;
	char blastIn;
} Options;

Options Opt;

NodeSet *createNodes();
Node *getNode();
Node *dupNode(), *dupNodeID();
TreeNode *dupTreeNode();
/*
Node *addNode(NodeSet *nodes, char *name, int cnt, int len,
	int minlen,int maxlen, Edge *child, char flag, specFlag spflag);
*/
Node *addNode(NodeSet *nodes, char *name, int cnt, int len,
	Region *consreg, int totlen,
	Edge *child, NodeFlag flag, specFlag spflag, char truncFlag);
EdgeSet *createEdges();
Edge *getBestEdge();
Edge *getEdgeByNode(), *getEdgeByNode0();

Edge *addEdge(EdgeSet *, Node*, Node*, Region *, Region *, Dist,
		ConnCount, signed char);
Edge *addEdgeWithScore(EdgeSet *, Node*, Node*, Region *, Region *, Dist, Dist,
		ConnCount, signed char);
Node *getRootNext();

int calNewDist(int nn, Count cnt1, Count cnt2, Dist dist1, Dist dist2,
                Dist score1, Dist score2, ConnCount conn1, ConnCount conn2,
                Dist *newdist, Dist *newscore, ConnCount *newconn);
OutGrpStat outgroupCheck(Node *node), outgroupCheckRaw(Node *node),
	outgroupCheck_TreeNode(TreeNode *tnode), outgroupCheckRaw_TreeNode(TreeNode *tnode);

int restoreGraph(char *, SimGraph *, int, char**);
pairList *getPairByNode(), *getPairByNode0();
Cluster *slink();
Region *transformReg();
char *addName();
char *getName();
void printNode();
void printEdge();
void printEdge2();
void printDomain();
void printTreeNode();
int deleted();
int isNotRoot();
int nodeMarked();
int edgeMarked();
int cmpr_nodeid();
Dist getEdgeScoreSum();
double DistRatio();
Node *getRootNext(NodeSetIter *);
Node *getRootNext_all(NodeSetIter *);
Node *getRootNext_sub(NodeSetIter *, int);
pList *convClustTree(NodeSet *);
#endif

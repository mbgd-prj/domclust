/*
 * DomClust: Hierarchical Clustering for Orthologous Domain Classification
 * Copyright (c) 2000-2007, Ikuo Uchiyama
 * All rights reserved.
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include "domclust.h"
#include "readfile.h"

#ifdef WITH_NEIGHBOR
#include "genedp.h"
#endif

/*
#define SWAP_INT(a,b) {int t = a; a = b; b = t;}
#define SWAP_USHORT(a,b) {unsigned short t = a; a = b; b = t;}
*/

char filename[100] = "stdin";
char genefilename[100];
char *geneclustfile;
char *sptreefile;
char *restorefile;
SimGraph SimG;
int getargs(int argc, char **argv);

/** setting default parameters */
defaultOpt() {
	Opt.distscale = DISTSCALE;
	Opt.logdistscale = 0;
	Opt.DEBUG = 0;
	Opt.VERBOSE = 0;
	Opt.verbose_step1 = 0;
	Opt.verbose_step2 = 0;
	Opt.minsp = 2;
	Opt.minent = 2;
	Opt.delete_small = 1;
	Opt.outstyle = NORMALOUT;

	Opt.minlen = 40;
	Opt.min_minlen = 40;
	Opt.minlen2 = 400; /* used in checkIncludedCluster() */

	/* used for overlap condition defined in overlapCheck() */
	/* ($ovlen > Opt.ovlpratio * $len1 && $ovlen > Opt.minlen)
	    || (($ovlen > Opt.ovlpratio2 * $len2) && $ovlen > Opt.minovlp)
		($len1 < $len2) */
	Opt.minovlp = 50;
	Opt.ovlpratio = 0.6;
	Opt.ovlpratio2 = 0.3;

	Opt.sim = 1;	/* 1..similarity, 0..distance; default=similarity */
	Opt.phylocutratio = 0.5;
	Opt.phylocutratio2 = 0;	/* checking distdiffcut */
	Opt.mincutcnt = 1;	/* minimum number of seq. for dom. cut */

	Opt.min_alilen = 250;

	Opt.distdiffcut = 0;
	Opt.domBoundary = B_BREAK;
	Opt.neighbor = 0;
	Opt.dpWin = 5;
	Opt.iniGap = 20.0;
	Opt.extGap = 0.2;
	Opt.skipGap = 0.1;

	Opt.nbrScoreRatio = 0.3;
	Opt.nbrScoreDecay = 1;
	Opt.nbrConnRatio = 0.5;

	Opt.adjInclRatio = 1.0;
	Opt.adjOvlpRatio = 0.0;

	Opt.noextend=0;
	Opt.blastIn=0;
}

main(int argc, char **argv)
{
	FILE *fp;
	Dist dist, score, scorecorr, distcorr;
	char name1[NAMELEN], name2[NAMELEN], lenstr[NAMELEN];
	int id1, id2, domnum;
	static char buf[BUFSIZ];
	int dir = 0;
	Node *node, *node1, *node2;
	Edge *edge;
	Cluster *sclust;
	int i;
	int from1, to1, from2, to2;
	Region ali1, ali2;
	int tmplen;
	int errflg;
	int scannum;
	SelFile selfile;
	SelData seldata;
	int minlen, alilen;
	unsigned long long cnt = 0;
	int spnum;

	defaultOpt();
	getargs(argc, argv);
	if (restorefile) {
		restoreGraph(restorefile, &SimG, argc, argv);
		spnum = getSPnum();
	}
	if (Opt.sim) {
		if (! Opt.cutoff) {
			Opt.cutoff = MINSCORE;
		}
		if (! Opt.missscore) {
			if (Opt.missratio) {
				Opt.missscore = Opt.cutoff * Opt.missratio;
			} else {
				Opt.missscore = Opt.cutoff * MISSRATIO;
			}
		} 
		if (! Opt.missdist) {
			if (Opt.missratio) {
				Opt.missdist = MAXDIST / Opt.missratio;
			} else {
				Opt.missdist = MAXDIST / MISSRATIO;
			}
		}
	} else {
		if (! Opt.cutoff) {
			Opt.cutoff = MAXDIST;
		}
		if (! Opt.missdist) {
			if (Opt.missratio) {
				Opt.missdist = Opt.cutoff / Opt.missratio;
			} else {
				Opt.missdist = Opt.cutoff / MISSRATIO;
			}
		}
		if (! Opt.missscore) {
			if (Opt.missratio) {
				Opt.missscore = MINSCORE * Opt.missratio;
			} else {
				Opt.missscore = MINSCORE * MISSRATIO;
			}
		}
	}
	if (Opt.cutoff <= 1) {
		Opt.decinum = 2;
		if (Opt.distscale < 100) {
			Opt.distscale = 100;
		}
	}
	if (Opt.taxMapOut) {
		if (! sptreefile) {
			fprintf(stderr, "%s: FATAL: taxonomy mapping requires sptree file that should be specified as -t<sptreefile>\n", argv[0]);
			exit(1);
		}
	}
	if (Opt.noReplaceSpTreeLeafName) {
		unset_ReplaceSpTreeLeafName();
	}
	if (sptreefile) {
		if (*sptreefile == '(') {
			parse_spinfo(sptreefile);
		} else {
			readSPfile(sptreefile);
		}
	} else if (Opt.treecheck) {
		fprintf(stderr, "No tree is specified "
			"-- the treecheck option is canceled\n");
		Opt.treecheck = 0;
	}

	if (restorefile) {
		setSPnum(spnum);
		setIgnoreSpec(Opt.metaSpecPref, Opt.metaSpecList, Opt.taxMapSpecList, Opt.partialSpecList);
		setOutGroup(Opt.outgroupStr, Opt.outgroupRev);
		setSpMask(Opt.spmaskStr);
		preproc_for_SpecInfo();
		if (Opt.taxMapOut) {
			/* ignore species absent in sptree */
			ignoreAbsentSpec();
		}
		outputResults(&SimG, argc, argv);
		exit(0);
	}

	SimG.nodes = createNodes();
	SimG.nhash = initNames(INIT_NODENUM);

	if (! *genefilename) {
		FILE *fp;
		strcpy(genefilename, filename);
		read_genefile(filename, &SimG, &fp);
		open_seldata(filename, &selfile, fp);
	} else {
		read_genefile(genefilename, &SimG, NULL);
		open_seldata(filename, &selfile, NULL);
	}
	if (geneclustfile) {
		FILE *fp;
		if ( (fp = fopen(geneclustfile, "r")) == NULL )  {
			fprintf(stderr, "Can't open file: %s\n", geneclustfile);
			exit(1);
		}
		read_geneclustfile(geneclustfile, &SimG);
	}
	if (Opt.blastIn) {
		set_blastIn(&selfile);
	}

	setIgnoreSpec(Opt.metaSpecPref, Opt.metaSpecList, Opt.taxMapSpecList, Opt.partialSpecList);
	setOutGroup(Opt.outgroupStr, Opt.outgroupRev);
	setSpMask(Opt.spmaskStr);
	preproc_for_SpecInfo();

	if (! SimG.nodes->nodenum) {
		fprintf(stderr, "No node read (Is the gene file correct?)\n");
	}

	SimG.edges = createEdges(SimG.nodes->nodenum);

	while ((scannum = read_seldata(&selfile, &seldata)) > 0) {
		if (scannum < MIN_HOMDATA) {
			/** read error **/
			continue;
		}
		if (Opt.VERBOSE) {
			if (++cnt % Opt.verbose_step1 == 0) {
				fprintf(stderr, "reading %lld\n", cnt);
			}
		}
		if (Opt.sim) {
			if (seldata.score < Opt.missscore) continue;
		} else {
			if (seldata.dist > Opt.missdist) continue;
		}


		id1 = getNameID(SimG.nhash, seldata.name1);
		id2 = getNameID(SimG.nhash, seldata.name2);

			/** swap seq1 and seq2 */
/*
		if (strcmp(seldata.name1, seldata.name2) > 0) {
printf("swap:%s,%s; %d,%d\n",seldata.name1,seldata.name2,id1,id2);
#ifdef LARGE
			SWAP_INT(seldata.from1, seldata.from2);
			SWAP_INT(seldata.to1, seldata.to2);
#else
			SWAP_USHORT(seldata.from1, seldata.from2);
			SWAP_USHORT(seldata.to1, seldata.to2);
#endif
			SWAP_INT(id1, id2);
printf("swap:%s,%s; %d,%d\n",seldata.name1,seldata.name2,id1,id2);
		}
*/

		errflg = 0;
		if ((node1 = getNode(SimG.nodes, id1)) == NULL) {
			fprintf(stderr, "not found: %s\n", seldata.name1);
			errflg = 1;
		}
		if ((node2 = getNode(SimG.nodes, id2)) == NULL) {
			fprintf(stderr, "not found: %s\n", seldata.name2);
			errflg = 1;
		}
		if (errflg) {
			if (Opt.skipErrEnt) {
				continue;
			} else {
				exit(1);
			}
		}

		if (node1->len < node2->len) {
			minlen = node1->len;
			alilen = (int)seldata.to1 - seldata.from1 + 1;
		} else {
			minlen = node2->len;
			alilen = (int)seldata.to2 - seldata.from2 + 1;
		}
if (seldata.to1 > node1->len) {
fprintf(stderr, "%s %d %d %d\n",node1->name,node1->len,seldata.from1,seldata.to1);
}
if (seldata.to2 > node2->len) {
fprintf(stderr, "%s %d %d %d\n",node2->name,node2->len,seldata.from2,seldata.to2);
}
		if (Opt.covfilt) {
			if ((float) alilen / minlen < Opt.covfilt) {
				continue;
			}
		}

		dist = (Dist) seldata.dist; score = (Dist) seldata.score;
#ifdef LARGE
		dir = seldata.dir;
#endif
		if (! Opt.sim) {
			if (alilen < Opt.min_alilen) {
			    minlen = (minlen < Opt.min_alilen)
				? minlen : Opt.min_alilen;
			    dist = (dist * alilen + (Dist) Opt.missdist *
					(minlen - alilen))/minlen;
			}
		}
		getDistCorr(node1, node2, &scorecorr, &distcorr);
		score += scorecorr; dist += distcorr;

		setSeqReg(&ali1, (SeqPos)seldata.from1, (SeqPos)seldata.to1);
		setSeqReg(&ali2, (SeqPos)seldata.from2, (SeqPos)seldata.to2);

		addEdgeWithScore(SimG.edges, node1, node2, &ali1, &ali2,
			(Dist) dist, (Dist) score, (ConnCount) 1, (char)dir);

	}
	if (SimG.edges->edgenum == 0) {
		fprintf(stderr, "No data\n");
		exit(0);
	}
	if (Opt.VERBOSE) {
		fprintf(stderr, "Done (%lld)\n", cnt);
	}
	if (SimG.edges->edgenum == 0) {
		fprintf(stderr, "No edge read (Is the homology file correct?)\n");
		exit(0);
	}
	setLeafNodeNum(SimG.nodes);

	if (Opt.VERBOSE) fprintf(stderr, "Sorting\n");

	createOrdList(&SimG);

#ifdef WITH_NEIGHBOR
	if (Opt.neighbor) {
		checkNeighbor2(SimG.edges, SimG.nodes);
	}
#endif
	if (Opt.VERBOSE) fprintf(stderr, "Indexing\n");
	createIndex(&SimG);
	if (Opt.VERBOSE) fprintf(stderr, "Done\n");

#ifdef WITH_DEBUG
/*
	sclust = slink(Edges, Nodes);
	printCluster(sclust, Nodes);
	free(sclust);
*/
#endif
	domCluster(SimG.edges, SimG.nodes, &Opt);
	outputResults(&SimG, argc, argv);
	exit(0);
}
/* a simple distance correction to avoid unnecessary split occurring when all sequences are identical */
getDistCorr(Node *node1, Node *node2, Dist *scoreCorr, Dist *distCorr)
{
	int flg1, flg2;
	Dist corr = (Dist) 1 / Opt.distscale;
	*scoreCorr = *distCorr = 0.0;
	/* intraspecies paralog should be connected first (add score) */
	if (spFlagEqual(node1->spflag, node2->spflag)) {
		*scoreCorr = corr;
		*distCorr = -corr;
		return(0);
	}
	/* outgroup species should be connected last (subtract score) */
	if (IS_OUTGROUP_MODE(Opt.outgroupMode)) {
		flg1 = outgroupFlagCheck(node1->spflag);
		flg2 = outgroupFlagCheck(node2->spflag);
		if ((flg1 != 0 && flg2 == 0) || (flg1 == 0 && flg2 != 0)) {
			*scoreCorr = -corr;
			*distCorr = corr;
		}
	}
}

getargs(int argc, char **argv)
{
	int i, nn = 1;
	char *p;
	char *scoreout_filename = NULL;
	char *taxinfo_filename = NULL;
	if (argc == 1 && isatty(0) == 1 && ! restorefile) {
		usage();
		fprintf(stderr, "  -h for more help\n");
		exit(0);
	}
	for (i = 1; i < argc; i++) {
		p = argv[i];
		if (*p == '-') {
			switch (*++p) {
			case 'a':
			/* parameters for merging adjacent clusters */
				p++;
				if (*p == 'i') {
					Opt.adjInclRatio = atof(++p);
				} else if (*p == 'o') {
					Opt.adjOvlpRatio = atof(++p);
				}
				break;
			case 'c':
			/* cutoff score or distance */
				if (isdigit(*++p)) {
					Opt.cutoff = (Dist) atof(p);
				} else if (*p=='i') {
					if (isdigit(*++p)) {
						Opt.in_cutoff = (Dist) atof(p);
					}
				} else if (*p == 'd') {
					Opt.distdiffcut = atof(++p);
				} else if (*p == 'b') {
					Opt.treeBalanceRatio = atof(++p);
				}
				break;
			case 'C':
			/* cutoff for splitting domains */
				if (*(++p)=='d') {
					Opt.cutoff2dist = atof(++p);
				} else if (*p=='r') {
					Opt.cutoff2ratio = atof(++p);
				} else if (*p=='l') {
					Opt.cutoff2_score_per_len = atof(++p);
				} else {
					Opt.cutoff2 = (Dist) atof(p);
				}
				break;
			case 's':
				Opt.sumcut = atoi(++p);
				break;
			case 'd':
				Opt.sim = 0;
				if (*(++p)) {
					Opt.cutoff = atoi(p);
				}
				break;
			case 'l':
				Opt.logdistscale = 1;
				break;
			case 'm':
			/* set missscore or missdist */
				if (*(++p)=='r') {
					Opt.missratio = atof(++p);
				} else {
					if (Opt.sim) {
						Opt.missscore = (Dist) atof(p);
					} else {
						Opt.missdist = (Dist) atof(p);
					}
				}
				break;
			case 'V':
			/* minimum alignment coverage for domain split */
				Opt.coverage = atof(++p);
				break;
			case 'v':
				Opt.VERBOSE = 1;
				int verbose_step1 = atoi(++p);
				if (verbose_step1) {
					Opt.verbose_step1 = verbose_step1;
				} else if (! Opt.verbose_step1) {
					Opt.verbose_step1 = 100000;
				}
				if (! Opt.verbose_step2) {
					Opt.verbose_step2 = 5000;
				}
				break;
			case 'D':
			/* debug flag */
				++p;
				if (*p=='\0') {
				} else if (isdigit(*p)) {
					Opt.DEBUG = Opt.DEBUG_val = atoi(p);
				} else {
					Opt.DEBUG_ent = p;
				}
				if (! Opt.DEBUG) {
					Opt.DEBUG = Opt.DEBUG_val = 1;
				}
				break;
			case 'o':
			/* output format */
				Opt.outstyle = atoi(++p);
				break;
			case 'n':
			/* minimum size of clusters to be output */
				p++;
				if (isdigit(*p)) {
					Opt.minsp = atoi(p);
					Opt.minent = Opt.minsp;
				} else if (*p == 's') {
					Opt.minsp = atoi(++p);
					if (Opt.minent<Opt.minsp){
						Opt.minent = Opt.minsp;
					}
				} else if (*p == 'e') {
					Opt.minent = atoi(++p);
					if (Opt.minent<Opt.minsp){
						Opt.minsp = Opt.minent;
					}
				} else if (*p == 'n') {
					Opt.delete_small = 0;
				}
				break;
			case 'S':
			/* use similarity measure */
				Opt.sim = 1;
				if (*(++p)) {
					Opt.cutoff = atoi(p);
				}
				break;
			case 'p':
			/* ratio for phylogenetic tree cutting */
				++p;
				if (*p == 'p') {
					Opt.phylocutratio2 = atof(++p);
				} else {
					Opt.phylocutratio = atof(p);
				}
				break;
			case 'H':
				if (*++p=='O') {
			/* homology-orthology hierarchical output */
					Opt.homClustOut = 1;
				} else {
			/* homologous rather than orthologous clustering */
					Opt.phylocutratio = 2.0;
				}
				break;
#ifdef WITH_NEIGHBOR
			case 'N':
				Opt.neighbor = 1;
				if (*++p == 'o') {
					Opt.neighbor = 2;
					if (isdigit(*++p)) {
						Opt.neighbor = atoi(p);
					}
				}
				Opt.logdistscale = 1;
				Opt.distscale = 2000;
				break;
#endif
			case 'R':
			/* restore from dump file */
				restorefile = ++p;
				break;
			case 't':
			/* specifying related genomes with a tree file */
				if (*(p+1)) {
					sptreefile = ++p;
				} else if (*argv[i+1]!='-'){
					sptreefile = argv[++i];
				}
				break;
			case 'g':
				Opt.genes = ++p;
				break;
			case 'O':
			case '-':
			/* additional options */
			/* many of them are experimental */
				++p;
				if (strncmp(p, "minovlp=", 8)==0) {
					Opt.minovlp = atoi(p+8);
				} else if (strncmp(p, "ovlpratio=", 10)==0) {
					Opt.ovlpratio = atof(p+10);
				} else if (strncmp(p, "ovlpratio2=", 11)==0) {
					Opt.ovlpratio2 = atof(p+11);
				} else if (strncmp(p, "minlen=", 7)==0) {
					Opt.minlen = atoi(p+7);
				} else if (strncmp(p, "minlen2=", 8)==0) {
					Opt.minlen2 = atoi(p+8);
				} else if (strncmp(p, "min_alilen=", 11)==0) {
					Opt.min_alilen = atof(p+11);
				} else if (strncmp(p, "covfilt=", 8)==0) {
					Opt.covfilt = atof(p+8);
				} else if (strncmp(p, "distscale=",10)==0) {
					Opt.distscale = atoi(p+10);
				} else if (strncmp(p, "dpwin=", 6)==0) {
					Opt.dpWin = atoi(p+6);
				} else if (strncmp(p, "nbrScoreRatio=", 14)==0) {
					Opt.nbrScoreRatio = atof(p+14);
				} else if (strncmp(p, "nbrScoreDecay=", 14)==0) {
					Opt.nbrScoreDecay = atof(p+14);
				} else if (strncmp(p, "nbrScoreLim=", 12)==0) {
					Opt.nbrScoreLim = atof(p+12);
				} else if (strncmp(p, "iniGap=", 7)==0) {
					Opt.iniGap = (Dist)atof(p+7);
				} else if (strncmp(p, "extGap=", 7)==0) {
					Opt.extGap = (Dist)atof(p+7);
				} else if (strncmp(p, "skipGap=", 8)==0) {
					Opt.skipGap = (Dist)atof(p+8);
				} else if (strncmp(p, "chkInclClst=", 12)==0) {
					Opt.adjInclRatio = atof(p+12);
				} else if (strncmp(p, "adjInclRatio=", 13)==0) {
					Opt.adjInclRatio = atof(p+13);
				} else if (strncmp(p, "chkOvlpClst=", 12)==0) {
					Opt.adjOvlpRatio = atof(p+12);
				} else if (strncmp(p, "adjOvlpRatio=", 13)==0) {
					Opt.adjOvlpRatio = atof(p+13);
				} else if (strncmp(p, "nochkInclClst", 13)==0) {
					Opt.adjInclRatio = 0;
				} else if (strncmp(p, "chkConnect", 10)==0) {
					if (p[10] == '=') {
						Opt.chkConnect = atof(p+11);
					} else {
						Opt.chkConnect = 2.0;
					}
				} else if (strncmp(p, "phyloMode=", 10)==0) {
					Opt.phyloMode = atoi(p+10);
				} else if (strncmp(p, "treecheck", 9)==0) {
					Opt.treecheck = 1;
				} else if (strncmp(p, "skipErrEnt", 10)==0) {
					Opt.skipErrEnt = 1;
				} else if (strncmp(p, "revMatch", 8)==0) {
					Opt.revMatch = 1;
				} else if (strncmp(p, "domBoundary", 11)==0) {
					Opt.domBoundary = B_ALIREG;
				} else if (strncmp(p, "mincutcnt=", 10)==0) {
					Opt.mincutcnt = atoi(p+10);
				} else if (strncmp(p, "domcut", 6)==0) {
					Opt.domcut = 1;
				} else if (strncmp(p, "nobreak", 7)==0) {
					Opt.nobreak = 1;
					Opt.noextend = 1;
				} else if (strncmp(p, "nbrConnRatio=", 13)==0) {
					Opt.nbrConnRatio = atof(p+13);
				} else if (strncmp(p, "outgroup=", 9)==0) {
					Opt.outgroupStr = p+9;
					Opt.outgroupRev = 0;
					if (! Opt.outgroupMode)
						Opt.outgroupMode = OutGroupMode;
				} else if (strncmp(p, "ingroup=", 8)==0) {
					Opt.outgroupStr = p+8;
					Opt.outgroupRev = 1;
					if (! Opt.outgroupMode)
						Opt.outgroupMode = OutGroupMode;
				} else if (strncmp(p, "homClustOut",11)==0) {
					Opt.homClustOut = 1;
				} else if (strncmp(p, "outputScore",11)==0) {
					Opt.outputScore = 1;
					if (p[11] == '=') {
						scoreout_filename = &p[12];
					}
				} else if (strncmp(p, "spmask=", 7)==0) {
					Opt.spmaskStr = p+7;
				} else if (strncmp(p, "outgroupMode=",13)==0){
					Opt.outgroupMode = atoi(p+13);
				} else if (strncmp(p, "outputSubClustTree",18)==0){
					Opt.outputSubClustTree = 1;
					if (strcmp(&p[18],"=0")==0) {
						Opt.outputSubClustTree = 0;
					}
				} else if (strncmp(p, "outputOuterClustTree",18)==0){
					Opt.outputOuterClustTree = 1;
				} else if (strncmp(p, "horizweight=",12)==0){
					Opt.horizweight = atof(p+12);
				} else if (strncmp(p, "ignoreTrunc", 11)==0) {
					Opt.ignoreTrunc = 1;
				} else if (strncmp(p, "noextend", 8)==0) {
					Opt.noextend = 1;
				} else if (strncmp(p, "metaPref=", 9)==0) {
					Opt.metaSpecPref = p+9;
					Opt.metagenomeMode = 1;
/*
					Opt.outputSubClustTree = 1;
*/
				} else if (strncmp(p, "meta=", 5)==0) {
					Opt.metaSpecList = p+5;
					Opt.metagenomeMode = 1;
/*
					Opt.outputSubClustTree = 1;
*/
				} else if (strncmp(p, "partial=", 8)==0) {
					Opt.partialSpecList = p+8;
				} else if (strncmp(p, "taxMapSpec=", 11)==0) {
					if (! Opt.taxMapOut) {
						Opt.taxMapOut = 1;
					}
					Opt.taxMapSpecList = p+11;
				} else if (strncmp(p, "taxMapOut", 9)==0) {
					Opt.taxMapOut = 1;
					if (p[9] == '=') {
						if (strcmp(&p[10],"2") == 0) {
							Opt.taxMapOut = 2;
						} else if (strcmp(&p[10],"3") == 0) {
							Opt.taxMapOut = 3;
						} else {
							Opt.taxMapOut = 2;
							taxinfo_filename = &p[10];
						}
					}
				} else if (strncmp(p, "noReplaceSpTreeLeafName", 23) == 0) {
					Opt.noReplaceSpTreeLeafName = 1;
				} else if (strncmp(p, "verbose_step1=", 14)==0) {
					Opt.verbose_step1 = atoi(p+14);
				} else if (strncmp(p, "verbose_step2=", 14)==0) {
					Opt.verbose_step2 = atoi(p+14);
				} else if (strncmp(p, "NEW", 3)==0) {
					/** version flag: temporal option **/
					Opt.NEW = 1;
				} else if (strncmp(p, "blastIn", 7)==0) {
					Opt.blastIn = 1;
				} else if (strncmp(p, "geneClustFile=", 14)==0) {
					geneclustfile =  p+14;
				}
				break;
			case 'h':
				help();
				exit(0);
			case '\0':
				if (nn == 1) {
					strcpy(filename, "stdin");
				} else {
					strcpy(genefilename, "stdin");
				}
				break;
			default:
				break;
			}
		} else {
			switch (nn) {
			case 1:
				strcpy(filename, p);
				break;
			case 2:
				strcpy(genefilename, p);
				break;
			}
			nn++;
		}
	}
	if (! Opt.cutoff2) {
		if (Opt.cutoff2ratio) {
			Opt.cutoff2 = Opt.cutoff * Opt.cutoff2ratio;
		} else if (Opt.sim) {
			Opt.cutoff2 = Opt.cutoff;
		}
	}
#ifdef WITH_NEIGHBOR
	if (Opt.neighbor && ! Opt.nbrConnGap) {
		Opt.nbrConnGap = Opt.dpWin;
	}
#endif
	if (scoreout_filename) {
		if ( (Opt.outputScoreFp = fopen(scoreout_filename, "w")) == NULL) {
			fprintf(stderr, "Can't open %s\n", scoreout_filename);
		}
	}
	if (taxinfo_filename) {
		if ( (Opt.taxMapOutFp = fopen(taxinfo_filename, "w")) == NULL) {
			fprintf(stderr, "Can't open %s\n", taxinfo_filename);
		}
	}
}
usage(void)
{
	fprintf(stderr, "DomClust ver.%s\n", DOMCLUST_VERSION);
	fprintf(stderr, "Usage: domclust [options] homfile genetab\n");
}
help()
{
	usage();
	fprintf(stderr, "\
    -S     use similarity as a measure of relatedness [on]\n\
    -d     use distance (or disimilarity) as a measure of relatedness\n\
    -c#    cutoff score/distance (can also be spcified as -S# or -d#) [60]\n\
    -ci#   cutoff score/distance among ingroup organisms\n\
    -m#    score/distance for missing relationships (m<c)\n\
    -mr#   specify a missing score as a ratio to c (0<mr<1) [0.95]\n\
    -C#    cutoff score for domain split (c<=C)\n\
    -V#    alignment coverage for domain split (0<=V<=1)\n\
    -n#    minimum # of organisms in clusters to be output [2]\n\
    -ne#   minimum # of entries in clusters to be output [2]\n\
    -p#    ratio of phylogenetic pattern overlap for tree cutting [0.5]\n\
    -H     homology clustering (i.e. skip the tree cutting)\n\
    -HO    hierachical combination of homology/orthology clustering\n\
    -ai#   member overlap for absorbing adjacent small clusters (0<=ai<=1)\n\
    -ao#   member overlap for merging adjacent clusters (0<=ao<=ai)\n\
    -t<fn> use a tree file for weighting related genomes\n\
    -R<fn> restore from dump file\n\
    -o#    output format (default: 0:SimpleList)\n\
            1:Tree, 2:Newick, 3:Newick with length (when -d is specified),\n\
            4:Hierarchical, 5:ClusterTable, 6:Graph, 9:Table, 10:Dump\n\
    -Ooutgroup=sp1,sp2,..  treat sp1,sp2,.. as outgroup\n\
    -Ohorizweight=#        relative weight for horiz. transfer (0<=x<=1)\n\
    -OoutputScore=#        output score/distance at the root of each cluster\n\
    -Ometa=sp1,..          the specified genomes are treated as metagenomes\n\
    -OtaxMapOut=#          output taxonomy mapping of metagenomic data\n\
    -OtaxMapSpec=sp1,..    target species for taxonomy mapping [=meta]\n\
");
}

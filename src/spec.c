/*
 * DomClust: Hierarchical Clustering for Orthologous Domain Classification
 * Copyright (c) 2000-2007, Ikuo Uchiyama
 * All rights reserved.
 */

#include <stdio.h>
#include <string.h>
#include "domclust.h"
#include "spec.h"
#include "namehash.h"

static int spflagsiz = SPFLAGSIZ;
static unsigned char bitmask[8] = {1,2,4,8,16,32,64,128};

static int bitcnt[256] = {
	0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,
	4,5,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,
	4,5,5,6,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,
	4,5,4,5,5,6,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,
	4,5,5,6,5,6,6,7,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,
	4,5,3,4,4,5,4,5,5,6,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,
	4,5,5,6,4,5,5,6,5,6,6,7,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,
	4,5,4,5,5,6,4,5,5,6,5,6,6,7,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
	4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8
};
static double wt_bitcnt[SPFLAGSIZ][256];
static int taxFlag[MAXTAXNUM][SPFLAGSIZ][256];
int _getSPid(char *sp, char add);
static int replaceLeafName = 1;

/*
static char *IgnoreSpecName;
static int IgnoreSpecNameLen;
*/

getSPid(char *sp)
{
	return _getSPid(sp, 1);
}
findSPid(char *sp)
{
	return _getSPid(sp, 0);
}
_getSPid(char *sp, char add)
{
	int id;
	if (! SpHash) {
		SpHash = initNames(MAXSP*2);
	}
	if ( (id = getNameID(SpHash, sp)) < 0 ) {
	    if (add) {
		/* Not found in the hash. Add this name and assign a new ID */
		id = SPnum++;
		if (SPnum >= MAXSP) {
			fprintf(stderr,
				"Too many species. Raise SPFLAGSIZ value. (%d,%d)\n",SPnum,(int)MAXSP);
			exit(1);
		}
		setSPnum(SPnum);
		SPnames[id] = addName(SpHash, sp, id);
		SPweights[id] = 1.0;
	    } else {
		return -1;
	    }
	}
	return id;
}
char *getSPname(int id)
{
	if (id < SPnum) {
		return SPnames[id];
	} else {
		return NULL;
	}
}

char unset_ReplaceSpTreeLeafName() {
	replaceLeafName = 0;
}
char *getTaxName(int id)
{
	if (replaceLeafName && sptree_isLeaf(id) && spTree.node[id].name == 0) {
		/* leaaf node: the name of the node should be stored in the parent node */
		int parent = spTree.node[id].parent;
		if (sptree_countChild(parent) == 1) {
			return(spTree.node[parent].name);
		}
	}
	if (spTree.node[id].name) {
		return(spTree.node[id].name);
	} else if (spTree.node[id].spid >= 0) {
		return(getSPname(spTree.node[id].spid));
	}
}
/****
setIgnoreSPname(char *name, char *splist)
{
	char *spname;
	IgnoreSpecName = name;
	IgnoreSpecNameLen = strlen(name);
}
check_ignoreSPname(char *name)
{
	if (IgnoreSpecNameLen &&
		strncmp(name, IgnoreSpecName, IgnoreSpecNameLen)==0) {
		return 1;
	} else if (getNameID(ignrSpecHash, name)) {
		return 1;
	}
	return 0;
}
****/

unknownOnlyCluster(specFlag spflag)
{
	/** all species are ignored (unknown) **/
	return (spFlagMaskCheck(spflag, SPflags.unknown)==2);
}
knownOnlyCluster(specFlag spflag)
{
	return (spFlagMaskCheck(spflag, SPflags.unknown)==0);
}
metagenomeOnlyCluster(specFlag spflag)
{
	return (spFlagMaskCheck(spflag, SPflags.meta)==2);
}
nonMetagenomeOnlyCluster(specFlag spflag)
{
	return (spFlagMaskCheck(spflag, SPflags.meta)==0);
}
partialOnlyCluster(specFlag spflag)
{
	return (spFlagMaskCheck(spflag, SPflags.partial)==2);
}
taxQueryCluster(specFlag spflag)
{

	return (spFlagMaskCheck(spflag, SPflags.taxQuery)==2);
}
taxKnownCluster(specFlag spflag)
{
	return (spFlagMaskCheck(spflag, SPflags.taxQuery)==0);
}

/*      0 .. no outgroup species,
	1 .. there is a species belonging to the outgroup,
	2 .. all species are belonging to the outgroup */
outgroupFlagCheck(specFlag spflag)
{
	return spFlagMaskCheck(spflag, SPflags.outGroup);
}
/*      0 .. no ingroup species,
	1 .. there is a species belonging to the ingroup,
	2 .. all species are belonging to the ingroup */
ingroupFlagCheck(specFlag spflag)
{
	return spFlagMaskCheck(spflag, SPflags.inGroup);
}


spFlagMaskCheck(specFlag spflag, specFlag spFlagPat)
{
	specFlag matched;
	int matched_spcnt, spcnt;

	spFlagAND(spflag, spFlagPat, matched);
	matched_spcnt = spFlagCntAll(matched);
	if (matched_spcnt) {
		spcnt = spFlagCntAll(spflag);
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

setSPnum(int spnum)
{
	SPnum = spnum;
	spflagsiz = spnum / 8 + 1;
}
getSPnum()
{
	return SPnum;
}
setSPweight(int id, double w)
{
	SPweights[id] = w;
}
double getSPweight(int id)
{
	return SPweights[id];
}

setSpMask(char *str)
{
	setSPflagBySPname(str, SPflags.spMask);
}
setOutGroup(char *str, char rev)
{
	if (rev) {
		setSPflagBySPname(str, SPflags.inGroup);
		/* reverse the inGroup flags to obtain the outGroup flags */
		spFlagCOMPL(SPflags.inGroup, SPflags.outGroup);
	} else {
		setSPflagBySPname(str, SPflags.outGroup);
		spFlagCOMPL(SPflags.outGroup, SPflags.inGroup);
	}
}
setIgnoreSpec(char *metaPref, char *metaList, char *taxQueryList, char *partialList)
{
/*
	int i;
*/
	/* metagenome */
	setSpecList(metaPref, metaList, SPflags.meta);
	/* taxmap query species */
	if (taxQueryList) {
		setSpecList(NULL, taxQueryList, SPflags.taxQuery);
	} else {
		copySPFlag(SPflags.meta, SPflags.taxQuery);
	}
	if (partialList) {
		setSpecList(NULL, partialList, SPflags.partial);
	} else {
		copySPFlag(SPflags.meta, SPflags.partial);
	}
	setIgnoreSpec0(SPflags.meta);
/**
	// meta genome species should be ignored in taxon count 
	copySPFlag(SPflags.meta, SPflags.ignore);

	// both metagenomes and taxQuery seq are ignored (unknown taxon) 
	copySPFlag(SPflags.meta, SPflags.unknown);
	
	for (i = 0; i < SPnum; i++) {
		if (getSPflag(SPflags.ignore, i)) {
			SPweights[i] = 0.0;
		}
	}
*/
}
setIgnoreSpec0(specFlag ignore_spflag)
{
	copySPFlag(ignore_spflag, SPflags.ignore);
	copySPFlag(ignore_spflag, SPflags.unknown);
	int i;
	for (i = 0; i < SPnum; i++) {
		if (getSPflag(SPflags.ignore, i)) {
			SPweights[i] = 0.0;
		}
	}
}
ignoreAbsentSpec()
{
	specFlag absentSpec, newIgnoreFlag;
	if (spTree.nodenum == 0) {
		/* do nothing */
		return(0);
	}
	/* species absent in the root node */
	spFlagCOMPL(spTree.node[1].spflag, absentSpec);
	/* add absentSpec to ignore spec */
	spFlagOR(SPflags.ignore, absentSpec, newIgnoreFlag);
	setIgnoreSpec0(newIgnoreFlag);
}

setSpecList(char *str, char *splist, specFlagP spFlag)
{
	int i;
	int mchlen;
	char *spname;
	NameHash *specHash;

	if (str == 0 && splist == 0) return 0;

	if (splist) {
		char *strtok_in = splist;
		specHash = initNames(MAXSP);
		while (spname = strtok(strtok_in, ",")) {
			addName(specHash, spname, 1);
			strtok_in = NULL;
		}
	}
	if (str) {
		mchlen = strlen(str);
	}
	for (i = 0; i < SPnum; i++) {
		if (splist) {
			if (getNameID(specHash, SPnames[i]) >= 0) {
				addSPflag(spFlag, i);
				continue;
			}
		}
		if (str) {
			if (strncmp(SPnames[i], str, mchlen) == 0) {
				addSPflag(spFlag, i);
			}
		}
	}
	return 0;
/*
	setSPflagBySPname(str, SPflags.ignore);
*/
}
setSPflagBySPname(char *str, specFlagP readGrpP)
{
	char *sp_ptr, *sp_next;
	int id, flag;
	static char *Delim = ":, ";
	static char sp_buf[SPNAMELEN];

	if (! str) return 0;

	sp_ptr = str;
	flag = 1;
	while ( flag ) {
		sp_next = strpbrk(sp_ptr, Delim);
		if (sp_next == NULL) {
			flag = 0;
			strcpy(sp_buf, sp_ptr);
		} else {
			strncpy(sp_buf, sp_ptr, sp_next - sp_ptr);
			sp_buf[sp_next - sp_ptr] = '\0';
		}
		if ( (id = _getSPid(sp_buf, 0)) < 0 ) {
			fprintf(stderr, "Not found: %s\n",sp_buf);
		} else {
			addSPflag(readGrpP, id);
		}
		sp_ptr = sp_next; sp_ptr++;
	}
	return 0;
}
spFlagInGroup(specFlag spflag, specFlag newflag)
{
	spFlagAND(spflag, SPflags.inGroup, newflag);
}
spFlagOutGroup(specFlag spflag, specFlag newflag)
{
	spFlagAND(spflag, SPflags.outGroup, newflag);
}

/** preproc_for_SpecInfo():
	this routine should be called before the clustering process
	but after reading the homology data (or reading all of
	the species names) **/
preproc_for_SpecInfo()
{
	sptree_add_species();
	create_wt_bitcnt();
}
create_wt_bitcnt()
{
	int i, j, jj, k;
	for (i = 0; i < spflagsiz; i++) {
		for (j = 0; j < 256; j++) {
			wt_bitcnt[i][j] = 0;
			jj = j;
			for (k = 0; k < 8 && jj>0; k++) {
				if (jj % 2 == 1) {
					wt_bitcnt[i][j] += SPweights[i*8+k];
				}
				jj /= 2;
			}
		}
	}
}
clearSPflag(specFlag spflag)
{
	bzero(spflag, SPFLAGSIZ);
}
setSPflag(specFlag spflag, int spnum)
{
	clearSPflag(spflag);
	spflag[spnum / 8] = bitmask[spnum % 8];
}
addSPflag(specFlag spflag, int spnum)
{
	spflag[spnum / 8] |= bitmask[spnum % 8];
}
getSPflag(specFlag spflag, int spnum)
{
	 return((spflag[spnum / 8] & bitmask[spnum % 8]) != 0);
}
copySPFlag(specFlag spflag1, specFlag spflag2)
{
	int i;
	for (i = 0; i < spflagsiz; i++) {
		spflag2[i] = spflag1[i];
	}
}

spFlagOR(specFlag spflag1, specFlag spflag2, specFlag newflag)
{
	register int i;
	for (i = 0; i < spflagsiz; i++) {
		newflag[i] = spflag1[i] | spflag2[i];
	}
}
spFlagAND(specFlag spflag1, specFlag spflag2, specFlag newflag)
{
	register int i;
	for (i = 0; i < spflagsiz; i++) {
		newflag[i] = spflag1[i] & spflag2[i];
	}
}
spFlagXOR(specFlag spflag1, specFlag spflag2, specFlag newflag)
{
	register int i;
	for (i = 0; i < spflagsiz; i++) {
		newflag[i] = spflag1[i] ^ spflag2[i];
	}
}
/* for flag clear */
spFlagANDNOT(specFlag spflag1, specFlag spflag2, specFlag newflag)
{
	specFlag not_spflag2;
	spFlagCOMPL(spflag2, not_spflag2);
	spFlagAND(spflag1, not_spflag2, newflag);
}
spFlagCOMPL(specFlag spflag, specFlag newflag)
{
	int i;
	for (i = 0; i < spflagsiz; i++) {
		newflag[i] = ~ spflag[i];
		newflag[i] &= ~ SPflags.ignore[i]; /* turn off ignore spec */
	}
}
spFlagANDcnt(specFlag spflag1, specFlag spflag2)
{
	/** ignore spec are not counted **/
	register int i;
	register int cnt = 0;
	for (i = 0; i < spflagsiz; i++) {
		cnt += bitcnt[
			spflag1[i] & spflag2[i] & (~ SPflags.ignore[i]) ];
	}
	return cnt;
}
spFlagANDcntAll(specFlag spflag1, specFlag spflag2)
{
	/** count all spec **/
	register int i;
	register int cnt = 0;
	for (i = 0; i < spflagsiz; i++) {
		cnt += bitcnt[spflag1[i] & spflag2[i]];
	}
	return cnt;
}
double spFlagANDcntW(specFlag spflag1, specFlag spflag2)
{
	register int i;
	double cnt = 0;
	specFlag tmpflag;
	spFlagAND(spflag1, spflag2, tmpflag);
	return spFlagCntW(tmpflag);
}
spFlagCnt(specFlag spflag)
{
	/** ignore spec are not counted **/
	register int i;
	register int cnt = 0;
	for (i = 0; i < spflagsiz; i++) {
		cnt += bitcnt[ spflag[i] & (~ SPflags.ignore[i]) ];
	}
	return cnt;
}
spFlagCntAll(specFlag spflag)
{
	/** count all spec **/
	register int i;
	register int cnt = 0;
	for (i = 0; i < spflagsiz; i++) {
		cnt += bitcnt[spflag[i]];
	}
	return cnt;
}
double spFlagCntW(specFlag spflag)
{
	register int i;
	register double cnt = 0;
	double cnt2 = 0;
	for (i = 0; i < spflagsiz; i++) {
		cnt += wt_bitcnt[i][spflag[i]];
	}
	cnt2 = sptree_spFlagCountTaxOrW(spflag);
	return cnt + cnt2;
}
double spFlagCntW_All(specFlag spflag)
{
	double cnt = spFlagCntW(spflag);
	if (spFlagANDcntAll(spflag, SPflags.ignore)) {
		/* add ignored species */ 
		cnt += 1.0;
	}
	return cnt;
}
/* spflag1 >= (includes) spflag2 */
spFlagInclude(specFlag spflag1, specFlag spflag2)
{
/*
printf ("%d,%d\n",spFlagCnt(spflag2), spFlagANDcnt(spflag1, spflag2));
*/
	return (spFlagCnt(spflag2) == spFlagANDcnt(spflag1, spflag2));
}
spFlagEqual(specFlag spflag1, specFlag spflag2)
{
	int i;
	for (i = 0; i < spflagsiz; i++) {
		if (spflag1[i] != spflag2[i]) {
			break;
		}
	}
	return (i == spflagsiz);
}

print_specFlag(specFlag spflag)
{
	int i;
	for (i = 0; i < SPnum; i++) {
		printf("%d", (spflag[i / 8] & bitmask[i % 8]) != 0);
	}
	putchar('\n');
}
dump_specFlag(FILE *ofp, specFlag spflag)
{
	int i;
	for (i = 0; i < spflagsiz; i++) {
		fprintf(ofp, " %d", spflag[i]);
	}
}
restore_specFlag(char *str, specFlag spflag)
{
	int i;
	for (i = 0; i < spflagsiz; i++) {
		spflag[i] = (unsigned char) strtol(str, &str, 10);
		if (! str) break;
	}
}

/*
 * DomClust: Hierarchical Clustering for Orthologous Domain Classification
 * Copyright (c) 2000-2007, Ikuo Uchiyama
 * All rights reserved.
 */

#include "domclust.h"
#include "vararray.h"
#include "util.h"
#include "seqreg.h"
#include "hitlist.h"
#include <assert.h>

initHitList(HitList *hitlist) {
	hitlist->numhit = 0;
}

addHit(HitList *hitlist, int seg1, Region *reg1, int offset1,
			int seg2, Region *reg2, int offset2,
			HitList_t *orighit)
{
	HitList_t *hit = &hitlist->hit[hitlist->numhit];
	int reglen = 0, orig_reglen;
	if (reg1) {
		hit->seg1 = seg1;
		copyReg(&hit->reg1, reg1);
		hit->offset1 = offset1;
	}
	if (reg2) {
		hit->seg2 = seg2;
		copyReg(&hit->reg2, reg2);
		hit->offset2 = offset2;
	}
	if (reg2) {
		reglen = regLen(reg2);
	} else {
		reglen = regLen(reg1);
	}
	orig_reglen = regLen(&(orighit->reg1));

/***
	if (orig_reglen - reglen > Opt.minlen){
		hit->score = orighit->score * reglen / orig_reglen;
	} else {
		hit->score = orighit->score;
	}
***/
	hit->score = orighit->score;

/**
if(Opt.DEBUG){
printf("reglen>>>%d,%d: score=%lf>%lf\n",
	orig_reglen,reglen,orighit->score,hit->score);
}
**/
	hit->dist = orighit->dist;
	hit->cnt = orighit->cnt;
	hit->connect =orighit->connect;
	hit->flag = 0;
	if (++hitlist->numhit >= MAXHIT) {
		fprintf(stderr, "Too many hit: raise MAXHIT (>%d)\n",
			hitlist->numhit);
		exit (1);
	}
}

postCheckHitList(HitList *hitlist, HitList_t *orighit, int *segLen) {
	int i, j;
	HitList_t *hit = hitlist->hit;
	int orig_reglen = 0, reglen = 0;
	int slen1, slen2;
	if (orighit) {
		orig_reglen = regLen(&orighit->reg1);
	}

	for (i = 0, j = 0; i < hitlist->numhit; i++) {
/***
		if (orighit && hit[i].score==0) {
			reglen = regLen(&hit[i].reg2);
			if (reglen < orig_reglen) {
				hit[i].score = orighit->score * reglen / orig_reglen;
			} else {
				hit[i].score = orighit->score;
			}
			hit[i].dist = orighit->dist;
			hit[i].connect = orighit->connect;
			hit[i].cnt = orighit->cnt;
		}
***/

/**
printf("len>%d,%d  %d,%d\n",
hit[i].reg1.to-hit[i].reg1.from+1,orig_reglen,hit[i].seg1,hit[i].seg2);
**/

		slen1 = segLen ? segLen[hit[i].seg1] : 0;
		slen2 = segLen ? segLen[hit[i].seg2] : 0;

		if (definedReg(&hit[i].reg1) && definedReg(&hit[i].reg2) &&
	  	   aliLenCheck(&hit[i].reg1, orig_reglen, slen1)
	   	&& aliLenCheck(&hit[i].reg2, orig_reglen, slen2) ){
			if (j != i) {
				copyHitList(&hit[j], &hit[i]);
			}
			j++;
		} else {
			/* skip; */
		}
	}
	hitlist->numhit = j;
}
copyHitList(HitList_t *hit1, HitList_t *hit2)
{
	copyReg(&(hit1->reg1), &(hit2->reg1));
	copyReg(&(hit1->reg2), &(hit2->reg2));
	hit1->seg1 = hit2->seg1;
	hit1->seg2 = hit2->seg2;
	hit1->score = hit2->score;
	hit1->dist = hit2->dist;
	hit1->cnt = hit2->cnt;
	hit1->connect = hit2->connect;
	hit1->flag = hit2->flag;
}

mergeHitList(HitList *hitlist1, HitList *hitlist2)
{
	int i, j;
	int flag;
	HitList_t *hit1 = hitlist1->hit, *hit2 = hitlist2->hit;
	for (i = 0; i < hitlist2->numhit; i++) {
/*
printf(">>>>");
printHitList_t(&hit2[i],3);
*/
		if (! definedReg(&(hit2[i].reg1)) || ! definedReg(&(hit2[i].reg2))) {
				continue;
		}
		flag = 0;
		for (j = 0; j < hitlist1->numhit; j++) {
			if (isMergeableHit(&hit1[j], &hit2[i], FlagSameSeq)) {
/*
printf("Merge>>>>%d,%d\n",i,j);
printHitList_t(&hit1[j],3);
printHitList_t(&hit2[i],3);
printf("//\n");
*/
				mergeOvlpHitReg(&hit1[j], &hit2[i]);
				flag = 1;
				break;
			}
		}
		if (!flag) {
			copyHitList(&hit1[hitlist1->numhit], &hit2[i]);
			if (++hitlist1->numhit >= MAXHIT) {
				fprintf(stderr, "Too many hit: raise MAXHIT\n");
				exit (1);
			}
		}
	}
}
diagDiff(HitList_t *hit1, HitList_t *hit2)
{
	int diff = abs( (hit1->reg1.from - hit2->reg1.from) - 
			(hit1->reg2.from - hit2->reg2.from) );
	return diff;
}
isMergeableHit(HitList_t *hit1, HitList_t *hit2, FlagMergeMode mflag)
{
	if (ovlpHitList(hit1, hit2, mflag)) {
		int diff = diagDiff(hit1, hit2);
		int len1 = regLen(&hit1->reg1);
		int len2 = regLen(&hit1->reg2);
		if (diff < len1 * .6  && diff < len2 * .6) {
			return 1;
		}
	}
	return 0;
}
isMergeableHit2(HitList_t *hit1, HitList_t *hit2, FlagMergeMode mflag)
{
	if (ovlpHitList(hit1, hit2, mflag)) {
		return 1;
	}
	return 0;
}
mergeOvlpHitReg(HitList_t *hit1, HitList_t *hit2)
{
	Region andreg, orreg;
	int ovlen1,ovlen2;
	int reglen1, reglen2;
	double ovlp1, ovlp2;
	Dist score;
	if (ovlpReg(&hit1->reg1, &hit2->reg1, &andreg, &orreg, Opt.minovlp)==0){
		copyReg(&hit1->reg1, &orreg);
		resetReg(&hit2->reg1);
	}
	ovlen1 = regLen(&andreg); if (ovlen1 < 0) ovlen1 = 0;
	if (ovlpReg(&hit1->reg2, &hit2->reg2, &andreg, &orreg, Opt.minovlp)==0){
		copyReg(&hit1->reg2, &orreg);
		resetReg(&hit2->reg2);
	}
	ovlen2 = regLen(&andreg); if (ovlen2 < 0) ovlen2 = 0;
	reglen1 = regLen(&(hit1->reg1)) + regLen(&(hit1->reg2));
	reglen2 = regLen(&(hit2->reg1)) + regLen(&(hit2->reg2));
	ovlp1 = (double) ovlen1 * 2 / reglen1;
	ovlp2 = (double) ovlen2 * 2 / reglen2;

	score = (hit1->score * ovlp1 + hit2->score * ovlp2) / 2;
	hit1->score = score + hit1->score * (1-ovlp1) + hit2->score * (1-ovlp2);
	hit2->dist = (hit1->dist + hit2->dist) / 2;
}

/***
mergeHitReg2(HitList_t *hit1, HitList_t *hit2, Count cnt1, Count cnt2)
{
	Region andreg, orreg;
	int ovlen1,ovlen2;
	int reglen1, reglen2;
	double ovlp1, ovlp2;
	if (ovlpReg(&hit1->reg1, &hit2->reg1, &andreg, &orreg, Opt.minovlp)==0){
		copyReg(&hit1->reg1, &orreg);
		resetReg(&hit2->reg1);
	}
	ovlen1 = regLen(&andreg);
	if (ovlpReg(&hit1->reg2, &hit2->reg2, &andreg, &orreg, Opt.minovlp)==0){
		copyReg(&hit1->reg2, &orreg);
		resetReg(&hit2->reg2);
	}
	ovlen2 = regLen(&andreg);
	reglen1 = regLen(&(hit1->reg1)) + regLen(&(hit1->reg2));
	reglen2 = regLen(&(hit2->reg1)) + regLen(&(hit2->reg2));
	ovlp1 = (double) ovlen1 * 2 / reglen1;
	ovlp2 = (double) ovlen2 * 2 / reglen2;

	hit1->score = (hit1->score * cnt1 + hit2->score * cnt2) / (cnt1+cnt2);
	hit1->dist = (hit1->dist * cnt1 + hit2->dist * cnt2 ) / (cnt1+cnt2);
	hit1->flag = FlagMerged;
	hit2->flag = FlagDeleted;
	hit1->connect = ((int)hit1->connect + hit2->connect > MAXUCHR)
                        ? MAXUCHR : hit1->connect + hit2->connect;
}
***/

ovlpHitList(HitList_t *hit1, HitList_t *hit2, FlagMergeMode mflag)
{
/*
printf(">>>\n");
printHitList_t(hit1,3);
printHitList_t(hit2,3);
printf(">>>\n");
*/
	if ( ((mflag != FlagSameSeq || hit1->seg1 == hit2->seg1) &&
		overlapCheck(&hit1->reg1, &hit2->reg1) == 0) &&
 	     ((mflag != FlagSameSeq || hit1->seg2 == hit2->seg2) &&
		overlapCheck(&hit1->reg2, &hit2->reg2) == 0) ) {
		return 1;
	} else {
		return 0;
	}
}
printHitList(HitList *hitlist)
{
	int i;
	HitList_t *hit = hitlist->hit;
	printf("## hitlist hitnum=%d\n",hitlist->numhit);
	for (i = 0; i < hitlist->numhit; i++) {
		printHitList_t(&hit[i], 3);
	}
	printf("//\n");
}

printHitList_t(HitList_t *hit, int flag)
{
	if (flag & 1) {
		printf("%d:(%d,%d)",
			hit->seg1, hit->reg1.from, hit->reg1.to);
	}
	if (flag & 2) {
		if (flag & 1) {
			printf(" ");
		}
		printf("%d:(%d,%d)",
			hit->seg2, hit->reg2.from, hit->reg2.to);
	}
	printf(" %f,%d",hit->score,hit->cnt);
	printf("\n");
/*
	{
	 int l1 = regLen(&(hit->reg1)), l2 = regLen(&(hit->reg2));
	 if( ( (l1 > l2) ? ((double) l2 / l1) : ((double) l1 / l2) ) < 0.5 ){
	   fprintf(stderr, "Warning: unbalanced alignment: %d,%d\n",l1,l2);
 	   abort();
	 }
	}
*/
}

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

#define MAXHOMSEG 500
static int seg5[MAXHOMSEG], seg6[MAXHOMSEG];
static int matchpair[MAXHOMSEG][2];
calNewDist2(HitList_t *hit1, HitList_t *hit2, Count cnt1, Count cnt2,
		HitList_t *newhit);

mergeHitSeg(HitList *hitlist, Node *n1, Node *n2)
{
	int i, j;
	int n5, n6;
	int flghit1, flghit2;
	int delnum;
	HitList_t *hit = hitlist->hit;
	HitList_t misshit;
	Count cnt1, cnt2;

	misshit.dist = Opt.missdist;
	misshit.score = Opt.missscore;
	cnt1 = n1->cnt; cnt2 = n2->cnt;

	/** seg: D_MATCH=0, D_MATCH1=5, D_MATCH2=6 **/

/*
	for (i = 0; i < MAXHOMSEG; i++) {
		seg5[i] = seg6[i] = -1;
		matchpair[i][0] = matchpair[i][1] = -1;
	}
*/
	n5 = n6 = 0;
	flghit1 = flghit2 = 0;
	for (i = 0; i < hitlist->numhit; i++) {
		if ( hit[i].seg1 == D_MATCH1 ) {
			if (n5 >= MAXHOMSEG) {
				fprintf(stderr, "Too many HOMSEG: raise MAXHOMSEG > %d\n",n5);
				exit(1);
			}
			seg5[n5++] = i;
			flghit1 = 1;
		} else if ( hit[i].seg1 == D_MATCH2 ) {
			if (n6 >= MAXHOMSEG) {
				fprintf(stderr, "Too many HOMSEG: raise MAXHOMSEG > %d\n",n5);
				exit(1);
			}
			seg6[n6++] = i;
			flghit2 = 1;
		} else if ( hit[i].seg1 == D_LEFT1 || hit[i].seg1 == D_RIGHT1){
			flghit1 = 1;
		} else if ( hit[i].seg1 == D_LEFT2 || hit[i].seg1 == D_RIGHT2){
			flghit2 = 1;
		}
	}
	if (n5 > 0 && n6 > 0) {
	    for (i = 0; i < n5; i++) {
		for (j = 0; j < n6; j++) {
			if (isMergeableHit2(&hit[seg5[i]], &hit[seg6[j]],
					FlagDifferentSeq)) {
				calNewDist2(&hit[seg5[i]],&hit[seg6[j]],
					cnt1,cnt2,&hit[seg5[i]]);
				hit[seg5[i]].flag = FlagMerged;
				hit[seg6[j]].flag = FlagDeleted;
				ovlpReg(&(hit[seg5[i]].reg1),
					&(hit[seg6[j]].reg1),
					NULL,&(hit[seg5[i]].reg1),0);
				ovlpReg(&(hit[seg5[i]].reg2),
					&(hit[seg6[j]].reg2),
					NULL,&(hit[seg5[i]].reg2),0);
			}
		}
	    }
	}
	delnum = 0;
	for (i = 0, j = 0; i < hitlist->numhit; i++) {
		if (hit[i].flag == FlagDeleted) {
			delnum++;
			continue;
		} else if (hit[i].flag == FlagMerged) {
			hit[i].seg1 = 0;
		} else {
/*
printf("flghit>>>>%d,%d\n",flghit1,flghit2);
*/
			if ( hit[i].seg1 == D_MATCH1 && ! flghit2 ) {
				/** match against only seg1 **/
				calNewDist2(&hit[i], &misshit,
						cnt1,cnt2,&hit[i]);
				hit[i].seg1 = D_MATCH;
			} else if (hit[i].seg1 == D_MATCH2 && ! flghit1){
				/** match against only seg2 **/
				calNewDist2(&misshit, &hit[i],
						cnt1,cnt2,&hit[i]);
				hit[i].seg1 = D_MATCH;
			}
		}
		if (i > j) {
			copyHitList(&hit[j], &hit[i]);
		}
		j++;
	}
	hitlist->numhit -= delnum;
	mergeOvlpSeg(hitlist);
}
mergeOvlpSeg(HitList *hitlist, Node *n1, Node *n2)
{
	int i, j;
	int ov1;
	int delnum = 0;
	for (i = 0; i < hitlist->numhit; i++) {
		for (j = i + 1; j < hitlist->numhit; j++) {
/*
			if ((ovlpReg(&(hitlist->hit[i].reg1),
				&(hitlist->hit[j].reg1), NULL, NULL, 0) == 0) &&
			    (ovlpReg(&(hitlist->hit[i].reg2),
				&(hitlist->hit[j].reg2), NULL, NULL, 0) == 0)) {
*/
			if (ovlpHitList(&(hitlist->hit[i]), &(hitlist->hit[j]), FlagSameSeq)) {

				if (hitlist->hit[i].score > hitlist->hit[j].score) {
					hitlist->hit[j].flag = FlagDeleted;
				} else {
					hitlist->hit[i].flag = FlagDeleted;
				}
				delnum++;
			}
		}
	}
	hitlist->numhit -= delnum;
}
calNewDist2(HitList_t *hit1, HitList_t *hit2, Count cnt1, Count cnt2,
		HitList_t *newhit)
{
	newhit->dist =
	(Dist) (hit1->dist * cnt1 + hit2->dist * cnt2)
		/ (Dist) (cnt1 + cnt2);
	newhit->score =
	(Dist) (hit1->score * cnt1 + hit2->score * cnt2)
		/ (Dist) (cnt1 + cnt2);
	newhit->connect = ((int)hit1->connect + hit2->connect > MAXUCHR)
			? MAXUCHR : hit1->connect + hit2->connect;
}

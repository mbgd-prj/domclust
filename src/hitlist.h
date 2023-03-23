/*
 * DomClust: Hierarchical Clustering for Orthologous Domain Classification
 * Copyright (c) 2000-2007, Ikuo Uchiyama
 * All rights reserved.
 */

#ifndef _HITLIST_H
#define _HITLIST_H

#define MAXHIT 800
typedef enum {
	FlagHitList_NEW, FlagHitList_ADD, FlagHitList_MOD,
} FlagHitList;

typedef enum {
	FlagSameSeq, FlagDifferentSeq
} FlagMergeMode;

typedef enum {
	FlagNotMerged, FlagMerged, FlagDeleted
} FlagMergeStatus;

typedef struct {
	int seg1, seg2;
	Region reg1, reg2;
	int offset1, offset2;
	Dist score, dist;
	ConnCount connect;
	Count cnt;
	char flag;
} HitList_t;
typedef struct {
	HitList_t hit[MAXHIT];
	int numhit;
} HitList;
#endif

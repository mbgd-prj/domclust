/*
 * DomClust: Hierarchical Clustering for Orthologous Domain Classification
 * Copyright (c) 2000-2007, Ikuo Uchiyama
 * All rights reserved.
 */

#ifndef _SEQREG_H

#define regLen(reg) ((int) (reg)->to - (reg)->from + 1)
#define checkOvlpReg(reg1,reg2) \
        ((reg1)->from < (reg2)->to && (reg2)-> from < (reg1)->to)

/*
#define BIGPOS 32765
*/
#ifdef LARGE
#define BIGPOS 2147483647
#else
#define BIGPOS 65535
#endif

/*
#define BIGPOS 99999999
*/
#define INFPOS -BIGPOS
#define SUPPOS BIGPOS
#define definedPos(x) ((x) > -BIGPOS && (x) < BIGPOS)

#define printReg1(reg) {printf(">%s: ", #reg); printReg(reg);}
#define printReg2(reg, msg) {printf(">%s: ", (msg)); printReg(reg);}

/* sequece positions stored in Edge; use ushort to reduce memory usage */
#ifdef LARGE
typedef int StrSeqPos;
#else
typedef unsigned short StrSeqPos;
#endif
typedef struct {
        StrSeqPos from, to;
} StrRegion;

typedef int SeqPos;
typedef struct {
        SeqPos from, to;
} Region;

typedef struct {
	char type[5];
	short seg;
        SeqPos pos;
} PosList_t;
typedef struct {
	PosList_t *poslist;
	int num;
	int segnum;
	int maxsize;
} PosList;

double ovlpRatio(Region *, Region *, Region *);
double gapRatio(Region *, Region *);
int setSeqPos(Region *, SeqPos, SeqPos);

PosList* createPosList(int );

#define _SEQREG_H
#endif

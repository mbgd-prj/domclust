/*
 * DomClust: Hierarchical Clustering for Orthologous Domain Classification
 * Copyright (c) 2000-2007, Ikuo Uchiyama
 * All rights reserved.
 */

#ifndef _DOMCUT_H_
#define _DOMCUT_H_

typedef struct {
	SeqPos pos1, pos2;
} PosPair;

typedef struct {
	varArray *AliPosL, *AliPosR;
	double maxscoreL, maxscoreR;
	PosPair maxposL, maxposR;
	Region *ali12, *ali21;
	NodeSet *nodes;
} DomCut;
DomCut *createDomCut();
#endif

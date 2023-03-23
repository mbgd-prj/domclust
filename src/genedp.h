/*
 * DomClust: Hierarchical Clustering for Orthologous Domain Classification
 * Copyright (c) 2000-2007, Ikuo Uchiyama
 * All rights reserved.
 */

#ifndef _GENEDP_H_
#include "domclust.h"
#include "hash.h"

typedef struct GeneDP {
	int win;
	Dist *score, *Vscore;
	Dist iniGap;
	Dist extGap;
	Hash *hash;
} GeneDP;

GeneDP *createGeneDP();
Dist geneDP();
Edge *EdgeHashSearch();
Hash *createEdgeHash();

#define _GENEDP_H_
#endif

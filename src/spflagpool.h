/*
 * DomClust: Hierarchical Clustering for Orthologous Domain Classification
 * Copyright (c) 2000-2007, Ikuo Uchiyama
 * All rights reserved.
 */

#ifndef _SPFLAGPOOL_H_
#define _SPFLAGPOOL_H_

#include "namehash.h"
#include "spec.h"

void initSpecFlagPool(void);
specFlagP internSpecFlag(const unsigned char *buf);

extern specFlagP emptySpecFlag;

#endif

/*
 * djset.h
 * Copyright (c) 2000-2008, Ikuo Uchiyama
 * All rights reserved.
 */
#ifndef _DJSET_H_
typedef struct DJsetData {
	int parent;
	int rank;
} DJsetData;

typedef struct DJset {
	DJsetData *data;
	int numelem;
} DJset;
DJset *create_DJset();
#define _DJSET_H_
#endif

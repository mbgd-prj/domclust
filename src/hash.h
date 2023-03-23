/*
 * hash.h
 * Copyright (c) 2000-2007, Ikuo Uchiyama
 * All rights reserved.
 */

#ifndef _HASH_H_

typedef enum {FIND, ENTER} ACTION;

typedef unsigned long hashsize_t;

typedef struct {
	char *key;
	char *datum;
} HENTRY;

typedef struct HashDataBlock {
	HENTRY **datablock;
	int blksiz;
	int recnum;
} HashDataBlock;
typedef struct {
	HENTRY **table;
	hashsize_t hashsize;
	HashDataBlock *hashdatablock;
} Hash;
Hash *Hcreate();
typedef struct {
	char *key;
	int datumidx;
} SHash;
typedef struct {
	Hash *hash;
	int idx;
} HashIter;


HENTRY *nextHashIter();

#define hashf(x,i,hashsize) \
	(((x) % hashsize + (1 + (x) % (hashsize-2)) * (i)) % hashsize)


#define _HASH_H_
#endif

/*
 * DomClust: Hierarchical Clustering for Orthologous Domain Classification
 * Copyright (c) 2000-2007, Ikuo Uchiyama
 * All rights reserved.
 */

#include <stdlib.h>
#include <string.h>
#include "namehash.h"
#include "spflagpool.h"

#define POOL_HASHSIZ (1 << 20)   /* 2^20 = 1M buckets */

typedef struct PoolEntry {
    unsigned char *data;
    struct PoolEntry *next;
} PoolEntry;

static PoolEntry *pool_table[POOL_HASHSIZ];

specFlagP emptySpecFlag = NULL;

static unsigned int hash_spflag(const unsigned char *buf)
{
    /* FNV-1a over all SPFLAGSIZ bytes */
    unsigned int h = 2166136261u;
    int i;
    for (i = 0; i < SPFLAGSIZ; i++) {
        h ^= buf[i];
        h *= 16777619u;
    }
    return h & (POOL_HASHSIZ - 1);
}

specFlagP internSpecFlag(const unsigned char *buf)
{
    unsigned int h = hash_spflag(buf);
    PoolEntry *e;
    for (e = pool_table[h]; e; e = e->next) {
        if (memcmp(e->data, buf, SPFLAGSIZ) == 0)
            return e->data;
    }
    e = (PoolEntry *) malloc(sizeof(PoolEntry));
    e->data = (unsigned char *) malloc(SPFLAGSIZ);
    memcpy(e->data, buf, SPFLAGSIZ);
    e->next = pool_table[h];
    pool_table[h] = e;
    return e->data;
}

void initSpecFlagPool(void)
{
    static unsigned char zerobuf[SPFLAGSIZ]; /* zero-initialized by C */
    emptySpecFlag = internSpecFlag(zerobuf);
}

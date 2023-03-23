/*
 * spindex.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gdbm.h>
#include <errno.h>
#include "readfile.h"

#define VERSION "1.0.2"
#define FILENAMESIZE 1024
#define MAXNAME 25

/*
 *
 */
void print_usage(void) {
    printf("Spindex ver. %s\n", VERSION);
    printf("Usage :: spindex [-D dir] [-m | -d] spec1 [spec2 ... specN]\n");
    printf("         Default output directory is $DIR_SPINDEX.\n");
    printf("         Not defined $DIR_SPINDEX, then use $MBGD_HOME/database.work/bldp\n");
    printf("\n");
    printf("         default :: create spindex for spec.\n");
    printf("                    e.g. spindex eco     <-- create 'spindex.eco'\n");
    printf("         -m      :: merge spindex.spec to spindex.\n");
    printf("                    e.g. spindex -m eco  <-- merge 'spindex.eco' to 'spindex'\n");
    printf("         -d      :: dump spindex.\n");
    printf("                    e.g. spindex -d ''   <-- dump 'spindex'\n");
    printf("                         spindex -d eco  <-- dump 'spindex.eco'\n");
    printf("         -h      :: print this message.\n");

    return;
}

/*
 *
 */
void store_spindex(
    GDBM_FILE dbh,
    char *sp1,
    char *sp2,
    char *file_bldp,
    long pos_fp0,
    long pos_fp)
{
    char buf_gdbm_key[1024];
    char buf_gdbm_val[1024];
    datum key, val;

    sprintf(buf_gdbm_key, "%s:%s", sp1, sp2);
    sprintf(buf_gdbm_val, "%s:%ld:%ld", file_bldp, pos_fp0, pos_fp);

    key.dptr = buf_gdbm_key;
    key.dsize = strlen(buf_gdbm_key);
    val.dptr = buf_gdbm_val;
    val.dsize = strlen(buf_gdbm_val);
    gdbm_store(dbh, key, val, GDBM_REPLACE);

    return;
}

/*
 *
 */
int build_spindex(
    char *dir,
    char *spec,
    int  mode)
{
    int sta = 1;
    char path_bldp[FILENAMESIZE];
    char file_bldp[FILENAMESIZE];
    char path_gdbm[FILENAMESIZE];
    char file_gdbm[FILENAMESIZE];
    FILE *fp;
    GDBM_FILE dbh;
    int rsize;
    int save_errno;
    HomData buf_homdata;
    int size_homdata;
    long pos_fp, pos_fp0;

    /* */
    char sp1[MAXNAME], sp2[MAXNAME];
    char prev_sp1[MAXNAME], prev_sp2[MAXNAME];

    /* */
    strcpy(prev_sp1, "");
    strcpy(prev_sp2, "");
    size_homdata = sizeof(HomData);

    /* */
    sprintf(file_bldp, "blastdpres.%s", spec);
    sprintf(path_bldp, "%s/%s", dir, file_bldp);
    fp = fopen(path_bldp, "r");
    if (fp == NULL) {
        save_errno = errno;
        fprintf(stderr, "Can not fopen %s(%d)\n", path_bldp, save_errno);
        return -1;
    }

    /* */
    sprintf(file_gdbm, "spindex.%s", spec);
    sprintf(path_gdbm, "%s/%s", dir, file_gdbm);
    dbh = gdbm_open(path_gdbm, 0, mode, 0640, 0);
    if (dbh == NULL) {
        save_errno = errno;
        fprintf(stderr, "Can not gdbm_open %s(%d)\n", path_gdbm, save_errno);
        fclose(fp);
        return 0;
    }

    /* */
    pos_fp0 = ftell(fp);
    int i;
    for (;;) {
        pos_fp = ftell(fp);
        rsize = fread((char *)&buf_homdata, size_homdata, 1, fp);
        if (rsize == 0) {
            store_spindex(dbh, prev_sp1, prev_sp2, file_bldp, pos_fp0, pos_fp);
            break;
        }
        memset(sp1, 0, sizeof(sp1));
        for (i = 0; buf_homdata.name1[i] != ':'; i++) {
            sp1[i] = buf_homdata.name1[i];
        }
        memset(sp2, 0, sizeof(sp2));
        for (i = 0; buf_homdata.name2[i] != ':'; i++) {
            sp2[i] = buf_homdata.name2[i];
        }

        if ((pos_fp0 != pos_fp) && ((strcmp(prev_sp1, sp1) != 0) || (strcmp(prev_sp2, sp2) != 0))) {
            store_spindex(dbh, prev_sp1, prev_sp2, file_bldp, pos_fp0, pos_fp);
            pos_fp0 = pos_fp;
        }

        strcpy(prev_sp1, sp1);
        strcpy(prev_sp2, sp2);
    }

    /* */
    gdbm_close(dbh);

    /* */
    fclose(fp);
    
    return sta;
}

/*
 *
 */
int merge_spindex(
    char *dir,
    char *spec,
    int  mode)
{
    int sta = 1;
    char path_gdbm0[FILENAMESIZE];
    char file_gdbm0[FILENAMESIZE];
    char path_gdbm[FILENAMESIZE];
    char file_gdbm[FILENAMESIZE];
    FILE *fp;
    GDBM_FILE dbh0, dbh;
    int save_errno;
    datum key0, key, val;

    /* */
    sprintf(file_gdbm0, "spindex");
    sprintf(path_gdbm0, "%s/%s", dir, file_gdbm0);
    dbh0 = gdbm_open(path_gdbm0, 0, mode, 0640, 0);
    if (dbh0 == NULL) {
        save_errno = errno;
        fprintf(stderr, "Can not gdbm_open %s(%d)\n", path_gdbm0, save_errno);
        return 0;
    }

    /* */
    sprintf(file_gdbm, "spindex.%s", spec);
    sprintf(path_gdbm, "%s/%s", dir, file_gdbm);
    dbh = gdbm_open(path_gdbm, 0, GDBM_READER, 0440, 0);
    if (dbh == NULL) {
        save_errno = errno;
        fprintf(stderr, "Can not gdbm_open %s(%d)\n", path_gdbm, save_errno);
        gdbm_close(dbh0);
        return 0;
    }

    /* */
    key0.dptr = NULL;
    for (key = gdbm_firstkey(dbh); key.dptr != NULL; key = gdbm_nextkey(dbh, key)) {
        if (key0.dptr != NULL) {
            free(key0.dptr);
            key0.dptr = NULL;
        }

        val = gdbm_fetch(dbh, key);
        gdbm_store(dbh0, key, val, GDBM_REPLACE);

        free(val.dptr);
        key0.dptr = key.dptr;
    }
    if (key0.dptr != NULL) {
        free(key0.dptr);
        key0.dptr = NULL;
    }

    /* */
    gdbm_close(dbh);
    gdbm_close(dbh0);
    
    return sta;
}

/*
 *
 */
int dump_spindex(
    char *dir,
    char *spec,
    int  mode)
{
    int sta = 1;
    char path_gdbm[FILENAMESIZE];
    char file_gdbm[FILENAMESIZE];
    GDBM_FILE dbh;
    int save_errno;
    datum key0, key, val;

    /* */
    sprintf(file_gdbm, "spindex");
    if (strcmp(spec, "") != 0) {
        sprintf(file_gdbm, "spindex.%s", spec);
    }
    sprintf(path_gdbm, "%s/%s", dir, file_gdbm);
    dbh = gdbm_open(path_gdbm, 0, mode, 0440, 0);
    if (dbh == NULL) {
        save_errno = errno;
        fprintf(stderr, "Can not gdbm_open %s(%d)\n", path_gdbm, save_errno);
        return 0;
    }

    /* */
    key0.dptr = NULL;
    for (key = gdbm_firstkey(dbh); key.dptr != NULL; key = gdbm_nextkey(dbh, key)) {
        if (key0.dptr != NULL) {
            free(key0.dptr);
            key0.dptr = NULL;
        }

        val = gdbm_fetch(dbh, key);
        printf("%-.*s :: %-.*s\n", key.dsize, key.dptr,
                                   val.dsize, val.dptr);

        free(val.dptr);
    }
    if (key0.dptr != NULL) {
        free(key0.dptr);
        key0.dptr = NULL;
    }

    /* */
    gdbm_close(dbh);

    return sta;
}

/*
 *
 */
int main(
    int argc,
    char **argv)
{
    char dir[FILENAMESIZE];
    int mode;
    int i;
    int argshift = 0;

    /* */
    if ((2 <= argc) && (strcmp(argv[1], "-h") == 0)) {
        print_usage();
        exit(0);
    }

    /* */
    if ((2 <= argc) && (strcmp(argv[1], "-D") == 0)) {
	strcpy(dir, argv[2]);
	argv += 2;
	argc -= 2;
    } else if (getenv("DIR_SPINDEX") == NULL) {
        if (getenv("MBGD_HOME") == NULL) {
            printf("ERROR :: Please set 'MBGD_HOME'\n");
            exit(-1);
        }
        sprintf(dir, "%s/database.work/bldp", getenv("MBGD_HOME"));
    }
    else {
        sprintf(dir, "%s", getenv("DIR_SPINDEX"));
    }
    printf("Directory for spindex :: %s\n", dir);

    if ((2 <= argc) && (strcmp(argv[1], "-d") == 0)) {
        mode = GDBM_READER;
        for (i = 2; i < argc; i++) {
            printf("DUMP :: [%s]\n", argv[i]);
            dump_spindex(dir, argv[i], mode);
        }
    }
    else if ((2 <= argc) && (strcmp(argv[1], "-m") == 0)) {
        mode = GDBM_WRCREAT;
        for (i = 2; i < argc; i++) {
            printf("MERGE :: [%s]\n", argv[i]);
            merge_spindex(dir, argv[i], mode);
        }
    }
    else {
        mode = GDBM_WRCREAT;
        for (i = 1; i < argc; i++) {
            printf("CREATE :: [%s]\n", argv[i]);
            build_spindex(dir, argv[i], mode);
        }
    }

    exit(0);
}

/* eof */

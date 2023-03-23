/*
 * select.c
 * Copyright (c) 2000-2007, Ikuo Uchiyama
 * All rights reserved.
 */

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <dirent.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "memalloc.h"
#include "readfile.h"
#include "namehash.h"
#include "vararray.h"

#include <gdbm.h>

#define VERSION "1.0.3f"

#define MAXNAME 25
#define MAXSPEC 50000

#define MAXGENE 5000000
#define HASHSIZE 100000000
/*
#define HASHSIZE MAXGENE * 2
#define HASHSIZE2 MAXGENE * 20
*/
#define GENEHASHSIZE 800000

#define NAMELEN 64
#define MAXFILES 40480
#define NAMESPACE 500000000
#define GENENAMESPACE 500000

#define BLKSIZ 10000

#define SUFF_PAMSORT ".pam.sort.pack"
#define SUFF_SCORESORT ".score.sort.pack"
#define SUFF_DEFAULT ".pack"

#define FILENAMESIZE 1024
#define INDEX_FILE "spindex"
#define PREFIX "blastdpres"
#define SUFFIX "pack"
#define DIRINFO_FILE "dirinfo"
#define DELIM ':'

#define MAX_DIRNAME 15

char *filename[MAXFILES];
int numfiles;
char *Dirname = ".";
char *Prefix = PREFIX;
char *Suffix = SUFFIX;
HomData *homdata;
/*
char filereg[1000000];
*/
Hash *genehash;
int check_genenames;
int orthoflag;
int paircheck = 1;
int next;
char rank;
char disttype = 'p'; /* pam */
char plasmid_check;
char no_paralog_check;
char update_mode;
NameHash *nhash;
Alloc_Object *genes_obj, *local_gobj;
char verbose;

typedef struct {
	char file[FILENAMESIZE];
	long long begpos, endpos;
	int uniqcheck;
} FilePtr;
typedef struct {
	char *dirname;
	int uniqcheck;
	char *idxfile;
	char *prefix;
	char *suffix;
	char delim;
} DirInfo;
DirInfo *DirList[MAX_DIRNAME];
typedef struct {
	HomData *homdata;
	int numelem;
	char *homdata_reg;
	size_t alloc_length;
} HomDataRead;

static int sort_byname(HomData *a, HomData *b)
{
	int retval = strcmp(a->name1, b->name1);
	if (retval == 0) {
		retval = strcmp(a->name2, b->name2);
		if (retval == 0) {
			retval = a->score - b->score;
		}
	}
	return retval;
}
static int sort_byscore(HomData *a, HomData *b)
{
	return (b->score - a->score);
}
static int sort_bydist(HomData *a, HomData *b)
{
	return (a->pam - b->pam);
}

typedef struct {
	int score;
} Gene;

char buf[BUFSIZ];
char *spec[MAXSPEC];
int specnum = 0;
char *q_spec[MAXSPEC];
int q_specnum = 0;
char namespace[NAMESPACE];
char genenamespace[GENENAMESPACE];
int pamcut = 99999;
int scorecut = 0;
double identcut;
double evalcut = 1.0;
char print_all = 0;
int binout = 0;
int binin = 1;
int output_self = 0;

float dist;

Gene *createGene();
int currID = 1;
void destroyHomDataRead(HomDataRead *homd);

main(int argc, char **argv)
{

	get_envval();
	getargs(argc, argv);
	read_homfile_all();
	return(0);
}

get_envval()
{
	char *dirnames = getenv("BLASTDPDIR");
	if (dirnames) {
		Dirname = dirnames;
	}
}
getargs(int argc, char **argv)
{
	int i;
	char *specstr = NULL, *q_specstr = NULL, *p;
	int genenameflag = 0;
	if (argc == 1 && isatty(0) == 1) {
		usage();
		fprintf(stderr, "  -h for more help\n");
		exit(0);
	}

	for (i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			switch (*++argv[i]) {
			case 'b':
				binout = 1;
				printf("%c", (unsigned) BIN_MAGIC);
				break;
			case 'B':
				binin = 0;
				break;
			case 'P':
				pamcut = atoi(++argv[i]);
				break;
			case 'S':
				scorecut = atoi(++argv[i]);
				break;
			case 'E':
				evalcut = atof(++argv[i]);
				break;
			case 'I':
				identcut = atof(++argv[i]);
				break;
			case 's':
				specstr = ++argv[i];
				break;
			case 'q':
				q_specstr = ++argv[i];
				break;
			case 'u':
				update_mode = 1;
				break;
			case 'g':
				genenameflag = 1;
				check_genenames = 1;
				break;
			case 'G':
				genenameflag = 1;
				check_genenames = 2;
				break;
			case 'r':
				rank = *++argv[i];
				if (! rank){
					rank = 's';
				}
				break;
			case 'd':
				disttype = *++argv[i];
				break;
			case 'o':
				orthoflag = atoi(++argv[i]);
				if (orthoflag == 0) {
					orthoflag = 1;
				}
				break;
			case 'D':
				Dirname = ++argv[i];
				break;
			case 'f':
				argv[i]++;
				if (*argv[i]=='t') {
					Suffix = ++argv[i];
				} else if (*argv[i]=='h') {
					Prefix = ++argv[i];
				}
				break;
			case 'n':
				next = atoi(++argv[i]);
				break;
			case 'a':
				print_all = 1;
				break;
			case 'p':
				no_paralog_check = 1;
				break;
			case 'y':
				output_self = 1;
				break;
			case 'm':
				plasmid_check = 1;
				break;
			case 'v':
				verbose = 1;
				break;
			case 'h':
				usage_detail();
				exit(0);
			default:
				break;
			}
		} else {
			if (genenameflag) {
				strcpy(genenamespace, argv[i]);
				genenameflag = 0;
			} else if (numfiles < MAXFILES) {
				filename[numfiles++] = argv[i];
			}
		}
	}
	if (specstr) {
		for (p = strtok(specstr, ","); p; p = strtok(NULL, ",")) {
			spec[specnum] = p;
			if (++specnum > MAXSPEC) {
				fprintf(stderr, "too many species: %d\n", specnum);
			}
		}
	}
	if (q_specstr) {
		for (p = strtok(q_specstr, ","); p; p = strtok(NULL, ",")) {
			q_spec[q_specnum++] = p;
		}
	}
	makeDirList(Dirname, DirList);
	if (genenamespace[0]) {
		HENTRY ent;
		genehash = Hcreate(GENEHASHSIZE);
		for (p = strtok(genenamespace, ","); p; p = strtok(NULL,",")) {
			ent.key = p;
			ent.datum = NULL;
			Hsearch(genehash, &ent, ENTER);
		}
	}
}
usage() {
	fprintf(stderr, "Select ver.%s\n", VERSION);
	fprintf(stderr, "Usage: select -ssp1,sp2,.. [options]\n");
}
usage_detail() {
	usage();
	fprintf(stderr, "Options:\n"
	"	-s<splist>  list of species to compare\n"
	"	-S#  Score cutoff (delim=,)\n"
	"	-P#  PAM cutoff\n"
	"	-I#  Identity cutoff\n"
	"	-E#  E-value cutoff\n"
	"	-g <genelist>  genes list (delim=,): the listed genes vs others\n"
	"	-G <genelist>  genes list (delim=,): compare among the listed genes\n"
	"	-D<dirlist>  directory list (delim=:)\n"
	"	-r{s|d}  sort by score/PAM\n");
}

read_homfile_all()
{
	int i, j;
	int ord;
/*
	char filenamebuf[1024];
*/
	FilePtr fileptr;
	Hash *tot_scorehash;

	/* storing names and total maximum scores of all genes */
	nhash = initNames(MAXGENE);
	/* storing the maximal score among the interspecies homologs */
	genes_obj = (Alloc_Object *)
			init_alloc_object(sizeof(Gene), BLKSIZ);

	if (numfiles) {
		for (i = 0; i < numfiles; i++) {
/*
			sprintf(fileptr.file, "%s/%s",Dirname, filename[i]);
*/
			if (get_filename_simple(
				&fileptr, DirList, filename[i]) != 0) {
				read_homfile(&fileptr, 0, 1);
			}
		}
	} else if (specnum && q_specnum) {
		for (i = 0; i < q_specnum; i++) {
			for (j = 0; j < specnum; j++) {
				if ( (ord = get_filename(&fileptr, DirList, q_spec[i], spec[j])) != 0 ) {
					read_homfile(&fileptr, 0, ord);
				}
			}
			if (update_mode) {
				for (j = i; j < q_specnum; j++) {
					if (i == j && no_paralog_check) {
						continue;
					}
					if ( (ord = get_filename(&fileptr, DirList, q_spec[i], q_spec[j])) != 0 ) {
						read_homfile(&fileptr, 0, ord);
					}
				}
			}
		}
	} else if (specnum) {
		for (i = 0; i < specnum; i++) {
			for (j = i+1; j < specnum; j++) {
				if ( (ord = get_filename(&fileptr, DirList, spec[i], spec[j])) != 0 ) {
					read_homfile(&fileptr, 0, ord);
				}
			}
		}
		if (! no_paralog_check) {
		    for (i = 0; i < specnum; i++) {
			if ( (ord = get_filename(&fileptr, DirList, spec[i], spec[i])) != 0 ) {
				read_homfile(&fileptr, 1, ord);
			}
		    }
		}
	}
}
read_homfile(FilePtr *fptr, int samespec, int order)
{
	struct stat buf;
	int i, j, numelem;
	HomData *h;
	char *nptr = namespace;
	Gene *gptr = NULL;
	HENTRY ent, *e;
	int srchflg = 0;
	Hash *pairhash = NULL, *scorehash = NULL;
	int score, bscore1, bscore2;
	int nspacelen;
	HomDataRead homdataRead;
	HomData *homdata;
	
	homdataRead.homdata_reg = NULL;

	if (binin) {
		read_homfile_bin(fptr, &homdataRead);
	} else {
		read_homfile_ascii(fptr, &homdataRead);
	}
	homdata = homdataRead.homdata;
	numelem = homdataRead.numelem;

	if (rank=='s') {
		qsort(homdata, numelem, sizeof(HomData),
			(int (*) (const void*,const void*)) sort_byscore);
	} else if (rank=='d') {
		qsort(homdata, numelem, sizeof(HomData),
			(int (*) (const void*,const void*)) sort_bydist);
	}
	if (fptr->uniqcheck) {
		pairhash = Hcreate(HASHSIZE);
	}
	if (orthoflag) {
		scorehash = Hcreate(HASHSIZE);
		if (local_gobj) {
			free_object(local_gobj);
		}
		local_gobj = (Alloc_Object *)
			init_alloc_object(sizeof(Gene), BLKSIZ);
		gptr = createGene(local_gobj);
	}

	nptr = namespace;
	nspacelen = 0;

	for (j = 0; j < numelem; j++) {
		h = &homdata[j];
		if (check_genenames) {
			int flg1, flg2;
			ent.key=h->name1;
			flg1 = Hsearch(genehash, &ent, FIND);
			ent.key=h->name2;
			flg2 = Hsearch(genehash, &ent, FIND);
			if (check_genenames == 1 &&
				! (flg1 > 0 || flg2 > 0)) {
				continue;
			} else if (check_genenames == 2 &&
				! (flg1 > 0 && flg2 > 0)) {
				continue;
			}
		}
		if (h->pam > pamcut || h->score < scorecut
				|| h->eval > evalcut
				|| h->ident < identcut) {
			continue;
		}
		if (! output_self && strcmp(h->name1, h->name2) == 0) {
			continue;
		}
		if (disttype == 'i') {
			dist = 100 - h->ident;
		} else {
			dist = h->pam;
		}

		if (orthoflag) {
			int flg1, flg2;
			if(rank=='d') {
				score = dist;
			} else {
				score = h->score;
			}

			if (samespec) {
				/* the total maximum score */
				bscore1 = getTotalHash(h->name1);
				bscore2 = getTotalHash(h->name2);
			} else {
				/* maximum score in this genome pair */
				bscore1 = addHash(scorehash, h->name1, &gptr, score);
				bscore2 = addHash(scorehash, h->name2, &gptr, score);
			}

			if (rank=='d') {
				flg1 = ((bscore1+0.1)*100/(score+0.1)<next);
				flg2 = ((bscore2+0.1)*100/(score+0.1)<next);
			} else {
				flg1 = ((score+0.1)*100/(bscore1+0.1)<next);
				flg2 = ((score+0.1)*100/(bscore2+0.1)<next);
			}
			if (orthoflag == 1 && (flg1 || flg2)) {
				continue;
			} else if (orthoflag == 2 && (flg1 && flg2)) {
				continue;
			}
		}

		 if (fptr->uniqcheck) {
			if (strcmp(h->name1, h->name2) < 0) {
				ent.key = strcpy(nptr, h->name1);
				strcat(ent.key, h->name2);
			} else {
				ent.key = strcpy(nptr, h->name2);
				strcat(ent.key, h->name1);
			}
			srchflg = Hsearch(pairhash, &ent, ENTER);
		}

		if (! fptr->uniqcheck || srchflg == 0) {
			/* newly enterd */

			if (print_all) {
				print_all_homdata(h, order);
			} else {
				print_homdata(h, order);
			}

			nspacelen += strlen(nptr) + 1;
			nptr += strlen(nptr) + 1;

		} else if (srchflg == 1) {
			/* duplicated data: do nothing */
		} else {
			fprintf(stderr, "Hash table is full\n");
		}

	}

	if (fptr->uniqcheck) {
		Hdestroy(pairhash);
	}
	if (orthoflag) {
		Hdestroy(scorehash);
	}
	destroyHomDataRead(&homdataRead);
}
addHash(Hash *scorehash, char *name, Gene **gptr, int score)
{
	char *nptr;
	HENTRY ent, *e;
	int srchflg;
	int bscore;
	int nscore;
	Gene *gptr0;

	if ((nscore = getNameID(nhash, name)) <0) {
		/* new gene name */
		nptr = addName(nhash, name, score);
	} else {
		if ( (rank=='d' && nscore > score) ||
			(rank=='s' && nscore < score)) {
			/* maximum score */
			resetNameID(nhash, name, score);
		}
		nptr = getName(nhash, name);
	}
	ent.key = nptr;
	ent.datum = (void *) *gptr;
	if ((srchflg = Hsearch(scorehash, &ent, ENTER))
			== 0) {
		/* store the best */
		(*gptr)->score = score;
		(*gptr) = createGene(local_gobj);
		bscore = score;
	} else if (srchflg == 1) {
		/* found the best score */
		bscore = ((Gene *) ent.datum)->score;
	} else {
		fprintf(stderr, "Hash table is full\n");
		exit(1);
	}
	return bscore;
}
getTotalHash(char *name)
{
	Gene *gptr0;
	int nscore;
	if ((nscore = getNameID(nhash, name)) <0) {
		return 0;
	} else {
		return nscore;
	}
}
/*
read_homfile_bin(FilePtr *fptr, HomData **homdata, int *numelem)
*/
read_homfile_bin(FilePtr *fptr, HomDataRead *homd)
{
	struct stat buf;
	int fd;
	char tmpbuf[HOMFILE_NAMELEN + 2];
/*
	static char *homdata_reg;
	size_t alloc_length;
*/
	off_t alloc_begpos, data_begin_offset;
	int pagesize = getpagesize();


	if (homd->homdata_reg) {
		munmap((char*) homd->homdata_reg, homd->alloc_length);
	}
	if ((fd = open(fptr->file, O_RDONLY)) < 0) {
		fprintf(stderr, "Can't open file: %s\n", fptr->file);
		return 0;
	}

	fstat(fd, &buf);
	if (fptr->begpos >= 0) {
		/* adjust alloc_begpos to a multiple of the pagesize for mmap */
		data_begin_offset = fptr->begpos % pagesize;
		alloc_begpos = fptr->begpos - data_begin_offset;
		homd->alloc_length = fptr->endpos  - alloc_begpos;
	} else if (fptr->begpos == -1)  {
		/** entier file **/
		alloc_begpos = 0;
		homd->alloc_length = buf.st_size;
		data_begin_offset = 0;
	} else {
		/* error */
		fprintf(stderr, "Illegal position: %s: %ld\n",
				fptr->file, fptr->begpos);
		exit(1);
	}

	if ((homd->homdata_reg = mmap((caddr_t) 0, homd->alloc_length,
		(PROT_READ | PROT_WRITE), MAP_PRIVATE, fd, (off_t) alloc_begpos))
			== MAP_FAILED) {
		fprintf(stderr, "memory allocation error\n");
		exit(1);
	}
	if (fd >= 2) {
		close(fd);
	}
/*
	homd->numelem = (fptr->endpos - fptr->begpos) / sizeof(HomData);
*/
	if (fptr->begpos >= 0) {
		homd->homdata = (HomData*) (homd->homdata_reg + data_begin_offset);
		homd->numelem = (fptr->endpos - fptr->begpos) / sizeof(HomData);
	} else {
		homd->homdata = (HomData*) homd->homdata_reg;
		homd->numelem = buf.st_size / sizeof(HomData);
	}
	strncpy(tmpbuf, homd->homdata->name1, HOMFILE_NAMELEN+1);
	if (strlen(tmpbuf) > HOMFILE_NAMELEN) {
		fprintf(stderr, "format error? <%s..>\n", tmpbuf);
	}

}
/*
read_homfile_ascii(FilePtr *fptr, HomData **homdata, int *numelem)
*/
read_homfile_ascii(FilePtr *fptr, HomDataRead *homd)
{
	char name1[MAXNAME],name2[MAXNAME];
/*
	int from1, to1, len1, from2, to2, len2, alilen, bscore, pam, score;
	float ident, eval, mpam, spam;
*/
	int from1, to1, len1, from2, to2;
	float ident, eval, pam, score;
	FILE *fp;
	static char buf[BUFSIZ];
	varArray *arry;
	HomData homdata;
	off_t recpos;

/*
	homd->numelem = 0;
*/
	if ((fp = fopen(fptr->file, "r")) == NULL) {
		fprintf(stderr, "Can't open file: %s\n", fptr->file);
		exit(1);
	}
	arry = createVarArray(10000, sizeof(HomData));
	if (fptr->begpos >= 0) {
		fseek(fp, fptr->begpos, SEEK_SET);
	} else if (fptr->begpos == -1) {
		/* entire file */
	} else {
		/* error */
		fprintf(stderr, "Illegal position: %s: %ld\n",
				fptr->file, fptr->begpos);
		exit(1);
	}
	while (fgets(buf, BUFSIZ, fp) != NULL) {
/**
		fscanf(fp, "%s%s%d%d%d%d%d%d%d%f%d%f%d%f%f%d",
			name1,name2,&from1,&to1,&len1,&from2,&to2,&len2,
			&alilen,&ident,&bscore,&eval,&pam,&mpam,&spam,&score);
**/
		sscanf(buf, "%s%s%d%d%d%d%f%f%f%f",
			name1,name2,&from1,&to1,&from2,&to2,
			&ident,&eval,&score,&pam);
		strcpy(homdata.name1,name1); strcpy(homdata.name2,name2);
		homdata.from1 = from1; homdata.to1 = to1;
		homdata.from2 = from2; homdata.to2 = to2;
		homdata.ident = ident; homdata.eval = eval;
		homdata.score = score; homdata.pam = pam;
		addArray(arry, &homdata);
/*
		homd->numelem++;
*/
		if ((recpos = ftell(fp)) >= fptr->endpos) {
			break;
		}
	}
	homd->homdata = (HomData *) arry->array;
	homd->numelem = arraySize(arry);
	fclose(fp);
}
void destroyHomDataRead(HomDataRead *homd) {
	if (homd->homdata_reg) {
		munmap((char*) homd->homdata_reg, homd->alloc_length);
	}
}
/****
****/

print_homdata(HomData *homdata, int order)
{
/*
	printf("%s %s %d %d %d %d %d %d %d %d\n",
		homdata->name1, homdata->name2,
		homdata->from1, homdata->to1, homdata->len1,
		homdata->from2, homdata->to2, homdata->len2,
		homdata->pam, homdata->score);
*/
	float dist = (disttype=='i') ? (100-homdata->ident) : homdata->pam;

	if (plasmid_check) {
		remove_plasmid_suffix(homdata->name1);
		remove_plasmid_suffix(homdata->name2);
	}
	if (binout) {
		SelData seldata;
	    if (order > 0) {
		strcpy(seldata.name1,homdata->name1);
		strcpy(seldata.name2,homdata->name2);
		seldata.from1 = homdata->from1;
		seldata.to1 = homdata->to1;
		seldata.from2 = homdata->from2;
		seldata.to2 = homdata->to2;
	    } else {
		strcpy(seldata.name1,homdata->name2);
		strcpy(seldata.name2,homdata->name1);
		seldata.from1 = homdata->from2;
		seldata.to1 = homdata->to2;
		seldata.from2 = homdata->from1;
		seldata.to2 = homdata->to1;
	    }
		seldata.dist = dist;
		seldata.score = homdata->score;
		fwrite(&seldata, sizeof(SelData), 1, stdout);
	} else if (order > 0) {
		printf("%s %s %d %d %d %d %d %d\n",
		homdata->name1, homdata->name2, homdata->from1, homdata->to1,
		homdata->from2, homdata->to2,
		(int)dist, (int)homdata->score);
	} else {
		printf("%s %s %d %d %d %d %d %d\n",
		homdata->name2, homdata->name1, homdata->from2, homdata->to2,
		homdata->from1, homdata->to1,
		(int)dist, (int)homdata->score);
	}
}
print_all_homdata(HomData *homdata, int order)
{
	float dist = (disttype=='i') ? (100-homdata->ident) : homdata->pam;
	if (plasmid_check) {
		remove_plasmid_suffix(homdata->name1);
		remove_plasmid_suffix(homdata->name2);
	}
	if (order > 0) {
	    printf("%s %s %d %d %d %d %.2f %.2g %d %d\n",
		homdata->name1, homdata->name2, homdata->from1, homdata->to1,
		homdata->from2, homdata->to2, homdata->ident, homdata->eval,
		(int) dist, (int)homdata->score);
	} else {
	    printf("%s %s %d %d %d %d %.2f %.2g %d %d\n",
		homdata->name2, homdata->name1, homdata->from2, homdata->to2,
		homdata->from1, homdata->to1, homdata->ident, homdata->eval,
		(int) dist, (int)homdata->score);
	}
}

get_filename_simple(FilePtr *fptr, DirInfo *dirList, char *filename)
{
	FILE *fp;
	int i;
	for (i = 0; dirList[i].dirname; i++) {
		sprintf(fptr->file, "%s/%s",dirList[i].dirname, filename);
		if ((fp = fopen(fptr->file, "r")) != NULL) {
			fclose(fp);
			fptr->begpos = -1; 
			fptr->uniqcheck = dirList[i].uniqcheck;
			return 1;
		}
	}
	return 0;
}

get_filename(FilePtr *fptr, DirInfo *dirList, char *spec1, char *spec2)
{
	int ret;
	int i;
	for (i = 0; dirList[i].dirname; i++) {
/*
fprintf(stderr, "dirList: %s\n",dirList[i].dirname);
*/
		ret = get_filename_sub(fptr, &dirList[i], spec1, spec2);
		if (ret != 0) {

/*
printf(">>>%s,%s: %s, %s, %d, %d, %d\n", spec1,spec2, dirList[i].dirname,
	fptr->file,fptr->begpos, fptr->endpos, dirList[i].uniqcheck);
*/

			return ret;
		}
	}
	return 0;
}

get_filename_sub(FilePtr *fptr, DirInfo *dirList, char *spec1, char *spec2)
{
	FILE *fp;
	int ord;
	char indexfile[FILENAMESIZE];
	FilePtr tmp_fptr;

	fptr->begpos = -1;
	fptr->uniqcheck = 1;
	fptr->uniqcheck = dirList->uniqcheck;

	if (dirList->suffix && dirList->suffix[0]) {
		sprintf(fptr->file, "%s/%s.%s-%s.%s", dirList->dirname,
			dirList->prefix, spec1, spec2, dirList->suffix);
	} else {
		sprintf(fptr->file, "%s/%s.%s-%s",
			dirList->dirname, dirList->prefix, spec1, spec2);
	}
	if ((fp = fopen(fptr->file, "r")) != NULL) {
		fclose(fp);
		return 1;
	}
	/* reverse */
	if (dirList->suffix && dirList->suffix[0]) {
		sprintf(fptr->file, "%s/%s.%s-%s.%s", dirList->dirname,
			dirList->prefix, spec2, spec1, dirList->suffix);
	} else {
		sprintf(fptr->file, "%s/%s.%s-%s",
			dirList->dirname, dirList->prefix, spec2, spec1);
	}
	if ((fp = fopen(fptr->file, "r")) != NULL) {
		fclose(fp);
		return -1;
	}

	sprintf(indexfile, "%s/%s", dirList->dirname, dirList->idxfile);
	ord = get_file_from_index(spec1, spec2, indexfile, &tmp_fptr, dirList);
	if (verbose) {
		fprintf(stderr, "%s; %s: ord=%d\n", fptr->file,tmp_fptr.file, ord);
	}
	if (ord) {
		sprintf(fptr->file, "%s/%s", dirList->dirname, tmp_fptr.file);
		fptr->begpos = tmp_fptr.begpos;
		fptr->endpos = tmp_fptr.endpos;

	}
	return ord;
}
get_file_from_index(char *spec1, char *spec2, char *idxfile, FilePtr *fptr,
			DirInfo *dirList)
{
	GDBM_FILE dbf;
	int ord = 0;

	if ( (dbf = gdbm_open(idxfile, 0, GDBM_READER, 0744, 0)) == NULL ) {
/*
		fprintf(stderr, "index file %s open error\n", idxfile);
*/
		return 0;
	}
	if (get_file_from_index0(spec1, spec2, fptr, dbf, dirList)) {
		ord = 1;
	} else if (get_file_from_index0(spec2, spec1, fptr, dbf, dirList)) {
		ord = -1;
	}
/*
fprintf(stderr, "ord=%d\n",ord);
*/
	gdbm_close(dbf);
	return ord;
}
get_file_from_index0(char *spec1, char *spec2, FilePtr *fptr,
			GDBM_FILE dbf, DirInfo *dirList) {
	datum key, cont;
	char spbuf[1024];
	char databuf[1024];
	char *p;
	char delim = dirList->delim;

	sprintf(spbuf, "%s%c%s", spec1, delim, spec2);
	key.dptr = spbuf;
	key.dsize = strlen(spbuf);
	if (gdbm_exists(dbf, key)) {
		cont = gdbm_fetch(dbf, key);
		strncpy(databuf, cont.dptr, cont.dsize);
		databuf[cont.dsize] = '\0';
		if ((p = strtok(databuf, ":")) == NULL) return 0;
		strcpy(fptr->file, p);
		if ((p = strtok(NULL, ":")) == NULL) return 0;
		fptr->begpos = atoll(p);
		if ((p = strtok(NULL, ":")) == NULL) return 0;
		fptr->endpos = atoll(p);
		if (verbose) {
			fprintf(stderr, ">>%s - %s, %s, %lld, %lld<<<\n",spec1, spec2, fptr->file, fptr->begpos, fptr->endpos);
		}
		return 1;
	} else {
		if (verbose) {
			fprintf(stderr, "key not found: %s\n", spbuf);
		}
	}
	return 0;
}

Gene *createGene(Alloc_Object *gobj)
{
	return (Gene *) memalloc(gobj);
}

remove_plasmid_suffix(char *name)
{
	char *p;
	if (p = strstr(name, "_p:")) {
		while (*(p+2)) {
			*p = *(p+2);
			p++;
		}
		*p = '\0';
	}
}

makeDirList(char *Dirname, DirInfo *DirList) {
	char *p;
	int i = 0;
	FILE *fp;
	char tmpfile[FILENAMESIZE];
	char buf[BUFSIZ];
	char tmpDirname[FILENAMESIZE];
	strcpy(tmpDirname, Dirname);
	for (p = strtok(tmpDirname, ":"); p; p = strtok(NULL, ":")) {
		if (i+1 > MAX_DIRNAME) {
			fprintf(stderr, "DIRNAME_LIST overflows\n");
			exit(1);
		}
		DirList[i].dirname = strdup(p);
		initDirList(&DirList[i]);

		sprintf(tmpfile, "%s/%s", p, DIRINFO_FILE);
		if ((fp = fopen(tmpfile, "r")) != NULL) {
			while (fgets(buf, sizeof(buf), fp) != NULL) {
				chomp(buf);
				if (strncmp(buf, "uniqcheck=",10)==0) {
					DirList[i].uniqcheck = atoi(&buf[10]);
				} else if (strncmp(buf, "idxfile=", 8) == 0) {
					DirList[i].idxfile = strdup(&buf[8]);
				} else if (strncmp(buf, "prefix=", 7) == 0) {
					DirList[i].prefix = strdup(&buf[7]);
				} else if (strncmp(buf, "suffix=", 7) == 0) {
					DirList[i].suffix = strdup(&buf[7]);
				} else if (strncmp(buf, "delim=", 6) == 0) {
					DirList[i].delim = buf[6];
				}
			}
		}
		i++;
	}
}
initDirList(DirInfo *dlist) {
	dlist->uniqcheck = 1;
	dlist->idxfile = INDEX_FILE;
	dlist->prefix = Prefix;
	dlist->suffix = Suffix;
	dlist->delim = DELIM;
}
chomp(char *str) {
	int len = strlen(str);
	if (str[len-1] == '\n') {
		str[len-1] = '\0';
	}
}

/*
 * DomClust: Hierarchical Clustering for Orthologous Domain Classification
 * Copyright (c) 2000-2007, Ikuo Uchiyama
 * All rights reserved.
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "domclust.h"
#include "readfile.h"

read_seldata(SelFile *selfile, SelData *seldata) {
	int retval = selfile->read(selfile->fp, seldata);
	if (selfile->rsize && retval && selfile->rsize != retval) {
		fprintf(stderr, "read error: %d,%d\n", selfile->rsize, retval);
		retval = -1;
	} else {
		selfile->rsize = retval;
	}
	return retval;
}

FILE *open_seldata(char *filename, SelFile *selfile, FILE *infp)
{
	FILE *fp;
	unsigned char c;
	if (infp != NULL) {
		fp = infp;
	} else if (strcmp(filename,"stdin") == 0) {
		fp = stdin; 
	} else {
		if ( (fp = fopen(filename,"r")) == NULL) {
			fprintf(stderr, "Can't open %s\n", filename);
			exit(1);
		}
	}
	selfile->fp = fp;
	if ((c = fgetc(fp)) == BIN_MAGIC) {
		selfile->read = read_seldata_bin;
	} else {
		selfile->read = read_seldata_ascii;
		ungetc(c, fp);
	}
	selfile->rsize = 0;
}
set_blastIn(SelFile *selfile) {
	selfile->read = read_blast_tabout;
}

static char buf[BUFSIZ];
read_seldata_ascii(FILE *fp, SelData *seldata)
{
	int scannum, readerr = 0;
	while (fgets(buf, BUFSIZ, fp) != NULL) {
		if (buf[0] == '#' || buf[0] == '\n') {
			fprintf(stderr, "skip: %s", buf);
			continue;
		}
#ifdef LARGE
		scannum = sscanf(buf, "%s %s %d %d %d %d %f %f %d",
#else
		scannum = sscanf(buf, "%s %s %hu %hu %hu %hu %f %f",
#endif
			seldata->name1, seldata->name2,
			&(seldata->from1), &(seldata->to1),
			&(seldata->from2), &(seldata->to2),
			&(seldata->dist), &(seldata->score)
#ifdef LARGE
			, &(seldata->dir)
#endif
			);
/*
printf("#####>>>>>>>%s,%s; %d,%d,%d,%d; dir=%d\n", seldata->name1, seldata->name2,
seldata->from1, seldata->to1, seldata->from2, seldata->to2, seldata->dir);
*/
		if (scannum < MIN_HOMDATA) {
			char *p = buf;
			while (*p && isspace(*p)) {
				p++;
			}
			if (*p==0) {
				fprintf(stderr, "skip: %s", buf);
			} else {
				fprintf(stderr, "read failed: %s", buf);
				readerr++;
				if (readerr > 100) {
					fprintf(stderr, "Too many read error: stop\n");
					exit(1);
				}
			}
			continue;
		} else if (scannum == MIN_HOMDATA) {
			if (Opt.sim) {
				seldata->score = seldata->dist;
			}
		}
#ifdef LARGE
#define LSIZ_HOMDATA 9
		if (scannum < LSIZ_HOMDATA) {
			seldata->dir = 1;
			if (seldata->from1 > seldata->to1) {
				int tmp = seldata->from1;
				seldata->from1 = seldata->to1;
				seldata->to1 = tmp;
				seldata->dir *= -1;
			}
			if (seldata->from2 > seldata->to2) {
				int tmp = seldata->from2;
				seldata->from2 = seldata->to2;
				seldata->to2 = tmp;
				seldata->dir *= -1;
			}
		}
#endif
		return scannum;
	}
	return 0;
}
read_seldata_bin(FILE *fp, SelData *seldata)
{
	int rsize;
	if ((rsize = fread(seldata, sizeof(SelData), 1, fp)) == 0) {
		if (feof(fp)) {
			return 0;
		} else {
			fprintf(stderr, "readerror\n");
			exit(1);
		}
	}
/*
	return rsize;
*/
	return MIN_HOMDATA;
}
read_blast_tabout(FILE *fp, SelData *seldata) {
	double ident, dist, evalue;
	int scannum;
	while (fgets(buf, BUFSIZ, fp) != NULL) {
		if (buf[0] == '#' || buf[0] == '\n') {
			continue;
		}
#ifdef LARGE
		scannum = sscanf(buf, "%s %s %f %*d %*d %*d %d %d %d %d %lg %f %d",
#else
		scannum = sscanf(buf, "%s %s %f %*d %*d %*d %hu %hu %hu %hu %lg %f",
#endif
			seldata->name1, seldata->name2, &ident,
			&(seldata->from1), &(seldata->to1), &(seldata->from2), &(seldata->to2),
			&evalue, &(seldata->score)
#ifdef LARGE
			, &(seldata->dir)
#endif

			);
		dist = 100.0 - ident;
		dist /= 100;
		dist = - log( 1 - dist - 0.2 * dist * dist);
		dist *= 100;
		seldata->dist = (float) dist;
		return scannum;
	}
	return 0;
}

read_genefile(char *filename, SimGraph *SimG, FILE **retfp)
{
	FILE *domfp;
/*
	static char buf[BUFSIZ];
*/
	char name[NAMELEN], spec[SPNAMELEN],
		spname[NAMELEN], *currname;
	int from, to, domnum, len, gpos, dir, prevdir, firstdir, truncFlag;
	int spnum = -1;
	int scannum, readerr = 0;
	DomID domid;
	SeqPos domlen[MAXDOM];
	int nodenum = 0;
	Node *node, *prevnode = NULL, *firstnode = NULL;
	specFlag spflag;
	char newseqflag = 0;
	int linenum = 0;
	Region tmpreg;

	if (! *filename || strcmp(filename,"stdin") == 0) {
		domfp = stdin;
	} else if ( (domfp = fopen(filename,"r")) == NULL) {
		fprintf(stderr, "Can't open %s\n", filename);
		exit(1);
	}
	domnum = 0;
	currname = NULL;
	while (fgets(buf, BUFSIZ, domfp) != NULL) {
		if (strncmp(buf, "//", 2) == 0) {
			/** end of the data **/
			break;
		} else if (buf[0] == '#') {
			/** end of the genome **/
			int i;
			for (i = 1; buf[i] && isspace(buf[i]); i++) {
			}
			if (buf[i] == '2') {
				/** circular genome **/
			    if (Opt.neighbor && ! Opt.adjInclRatio) {
				if (firstnode && prevnode) {
					if (gpos) {
						addNeighbor(prevnode, firstnode,
							prevdir, firstdir);
					}
				}
			    }
			} else {
				/** linear genome **/
			}
			continue;
		}
		gpos = dir = truncFlag = 0;
		scannum =
		  sscanf(buf, "%s %s %d %d %d %d", spec, name, &len, &gpos, &dir, &truncFlag);
		if (Opt.ignoreTrunc) {
			truncFlag = 0;
		}
		if (linenum++ == 0) {
			/* the first line */
			if (strlen(spec) > NAMELEN) {
				fprintf(stderr, "Too long species name\n");
				exit(1);
			}
		}
		if (scannum < 3) {
			fprintf(stderr, "read failed: %d: %s", linenum, buf);
			readerr++;
			if (readerr > 10) {
				fprintf(stderr, "Too many read error: stop\n");
				exit(1);
			}
			continue;
		}
		sprintf(spname, "%s:%s", spec, name);

		spnum = getSPid(spec);
		setSPflag(spflag, spnum);

		currname = addName(SimG->nhash, spname, nodenum);
		tmpreg.from = 1; tmpreg.to = len;
		node = addNode(SimG->nodes, currname, 1, len,
				&tmpreg, len, NULL, (NodeFlag) 0, spflag, truncFlag);
		if (dir == 0) dir = 1;
		node->dir = dir;
		if (newseqflag) {
			/* linear genome */
			firstnode = node;
			firstdir = dir;
		} else {
		    if (Opt.neighbor && ! Opt.adjInclRatio) {
			if (prevnode) {
				if (gpos) {
					addNeighbor(prevnode, node,
						prevdir, dir);
				}
			}
		    }
		}
		prevnode = node;
		prevdir = dir;
		nodenum++;
		newseqflag = 0;
	}
	if (retfp != NULL) {
		*retfp = domfp;
	} else if (domfp != stdin) {
		fclose(domfp);
	}
}
read_geneclustfile(char *filename, SimGraph *SimG)
{
	FILE *clstfp;
	char name[NAMELEN];
	char *spec;
	int nodeid;
	int scannum;
	int spnum;
	Node *node;
	if ( (clstfp = fopen(filename,"r")) == NULL) {
		fprintf(stderr, "Can't open %s\n", filename);
		exit(1);
	}
	while (fgets(buf, BUFSIZ, clstfp) != NULL) {
		if (buf[0] == '*') {
			scannum = sscanf(&buf[2], "%s", name);
			nodeid = getNameID(SimG->nhash, name);
			node = getNode(SimG->nodes, nodeid);
		} else if (buf[0] == ' ') {
			scannum = sscanf(&buf[2], "%s", name);
			spec = strtok(name, ":");
			if (spec == NULL) {
				spec = name;
			}
			spnum = getSPid(spec);
			addSPflag(node->spflag, spnum);
		}
	}
}

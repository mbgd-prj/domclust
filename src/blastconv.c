#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <dirent.h>
#include <ctype.h>

#define MAXMSP 40000
#define MAXCHR 128
#define AAnum 20
#define AAnum2 24

#define SEQBUFSIZ 1000000
#define MATRIXDIR "/bio/db/blast/matrix/aa/"

char *AA = "ARNDCQEGHILKMFPSTWYV";
char *AA2 = "ARNDCQEGHILKMFPSTWYVBZX*";
char AAidx[MAXCHR];
char AAidx2[MAXCHR];
int PSEUDOC = 5;
double EVAL_CUT = 0.01;
double BIAS_CUT = 0.8;
char seqbuf1[SEQBUFSIZ];
char seqbuf2[SEQBUFSIZ];
int Mat[AAnum2][AAnum2];
double Mexp;
int ignore_self = 1;
int ignore_second_for_allall = 1;
int do_not_output_alireg = 0;
int modified_eval = 0;
int verbose = 0;
int psi_maxiter = 0;
int flag = 0;
int recalc_score = 1;
double OvlpRatio = 0.8;

int GapOpen = 11;
int GapExt = 1;

int SCORE_CUT = -9999999;

double Prob[20] =
	{.07681,.05251,.04434,.05263,.01791,.04033,.06264,.07081,.02251,.05521,
	.09191,.05801,.02351,.03991,.05041,.07061,.05831,.01311,.03221,.06531};

typedef struct {
	double bitscore;
	int origscore;
	int score;
	double eval;
	int ident;
	char *qseq;
	char *sseq;
	int from1, to1, from2, to2;
	int delflag;
} MSP;
typedef struct {
	MSP Msp[MAXMSP];
	int MspNum;
	char qname[50];
	char sname[50];
	int qlen;
	int slen;
} MSPset;

MSPset MspSet;

char matfile[124] = MATRIXDIR "blosum62";
char *filenames[10000];
char filenamebuf[2000];
int filenum;

double calc_bias();
MSP *initmsp(), *nextmsp();

main(int argc, char **argv)
{
	int i;
	DIR *dir;
	struct dirent *e;
	char buf[250];

	if (argc <= 1) {
		die("Usage: blastconv [dirname|filenames ..]\n");
	}
	getargs(argc, argv);
	calc_prob();
	read_mat(matfile);
	if (filenum == 0) {
		read_blast(NULL);
	} else if (dir = opendir(filenames[0])) {
		i = 0;
		while (e = readdir(dir)) {
			if (e->d_name[0] == '.') continue;
			sprintf(filenamebuf,"%s/%s",filenames[0], e->d_name);
			read_blast(filenamebuf);
			if (++i % 10 == 0) {
				fprintf(stderr, ".");
			}
		}
	} else {
		for (i = 0; i < filenum; i++) {
			read_blast(filenames[i]);
		}
	}
}
getargs(int argc, char **argv)
{
	int i, j;
	char *p;
	filenum = 0;
	for (i = 1; i < argc; i++) {
		p = argv[i];
		if (*p == '-') {
			switch (*++p) {
			case 's':
				ignore_self = 0;
				break;
			case 'b':	/* output both */
				ignore_second_for_allall = 0;
				break;
			case 'o':
				do_not_output_alireg = 1;
				break;
			case 'e':
				EVAL_CUT = atof(++p);
				break;
			case 'E':
				modified_eval = 1;
				break;
			case 'B':
				BIAS_CUT = atof(++p);
				break;
			case 'v':
				verbose = 1;
				break;
			case 'j':
				psi_maxiter = atoi(++p);
				break;
			case 'c':
				SCORE_CUT = atoi(++p);
				break;
			case 'O':
				OvlpRatio = atof(++p);
				break;
			case 'M':
				strcpy(matfile, ++p);
				break;
			}
		} else {
			filenames[filenum++] = p;
		}
	}
}

char buf[BUFSIZ];
read_blast(char *filename)
{
	FILE *fp;
	char *p;
	MSP *msp = &(MspSet.Msp[0]);
	int pos1, pos2;
	double ratio;
	enum {BEGIN, NONE, QINFO, SINFO} status = BEGIN;
	int iternum;
	char msg[600];

	if (filename == NULL) {
		fp = stdin;
	} else if ((fp = fopen(filename, "r")) == NULL) {
		die("Can't open file\n");
	}
    while (fgets(buf, BUFSIZ, fp) != NULL) {
	if (status == BEGIN) {
		if (strncmp(buf, "BLAST", 5) != 0) {
			sprintf(msg, "%s: Not a blast file\n", filename);
			warn(msg);
			if (fp != stdin) {
				fclose(fp);
			}
			return 0;
		}
		status = NONE;
	} else if (strncmp(buf, "Query=", 6)==0) {
		for (p = &buf[6]; *p && *p != ' '; p++)
			;
		p = strtok(p, " \n");
		if (p) {
			strcpy(MspSet.qname, p);
		}
		status = QINFO;
		flag = 0;
	} else if (status == QINFO) {
		for (p = buf;  *p == ' '; p++) {
		}
		p++;
		MspSet.qlen = atoi(p);
		status = NONE;
	} else if (psi_maxiter && flag != 1) {
		if (strncmp(buf, "Results from round", 18)==0) {
			sscanf(&buf[19], "%d", &iternum);
			if (iternum >= psi_maxiter) {
				flag++;
			}
		} else if (strncmp(buf, "CONVERGED!", 10)==0) {
			flag++;
		}
		continue;
	} else if (strncmp(buf, ">", 1)==0) {
		scoreout(&MspSet);
		initmsp(&MspSet);
		msp = NULL;
		p = strtok(&buf[1], " \n");
		if (p) {
			strcpy(MspSet.sname, p);
		}
		status = SINFO;
	} else if (status == SINFO) {
		p = strstr(buf, "Length =");
		if (p) {
			MspSet.slen = atoi(p + 8);
			status = NONE;
		}
	} else if (strncmp(buf, " Score = ", 8) == 0) {
		if (msp) {
			ratio = calc_bias(msp);
			if (ratio <= BIAS_CUT) {
				if (verbose) {
					fprintf(stderr, "Ratio SKIP: %s:%s %lf\n",
					MspSet.qname,MspSet.sname,ratio);
				}
				/* cancel this msp */
				cancel_currentMSP(&MspSet);
			}
		}
		msp = nextmsp(&MspSet);
		p = strtok(&buf[8], " ");
		msp->bitscore = strtod(p, &p);
		if ((p = strtok(NULL, "(")) && (p = strtok(NULL, ")"))) {
			msp->score = msp->origscore = atoi(p);
			p = strtok(NULL, "=");
			p = strtok(NULL, " \n");
			if (*p == 'e') {
				static char tmpbuf[120];
				sprintf(tmpbuf, "1%s",p);
				msp->eval = strtod(tmpbuf, NULL);
			} else {
				msp->eval = strtod(p, NULL);
			}
		}
	} else if (strncmp(buf, " Identities = ", 14) == 0) {
		p = strtok(&buf[14], "(");
		p = strtok(NULL, ")");
		msp->ident = atoi(p);
	} else if (strncmp(buf, "Query:", 6) == 0) {
		p = strtok(&buf[7], " ");
		pos1 = atoi(p);
		p = strtok(NULL, " ");
		if (! p) {
			sprintf(msg, "parse error: %s: %s\n",filename,buf);
			warn(msg);
			continue;
		}
		strcat(msp->qseq, p);
		p = strtok(NULL, " \n");
		if (! p) {
			sprintf(msg, "parse error: %s: %s\n",filename,buf);
			warn(msg);
			continue;
		}
		pos2 = atoi(p);

		if (pos1 < pos2) {
			if (! msp->from1 || pos1 < msp->from1)
				msp->from1 = pos1;
			if (! msp->to1 || msp->to1 < pos2)
				msp->to1 = pos2;
		} else {
			if (! msp->from1 || pos2 < msp->from1)
				msp->from1 = pos2;
			if (! msp->to1 || msp->to1 < pos1)
				msp->to1 = pos1;
		}
	} else if (strncmp(buf, "Sbjct:", 6) == 0) {
		p = strtok(&buf[7], " ");
		pos1 = atoi(p);
		p = strtok(NULL, " ");
		if (! p) {
			sprintf(msg, "parse error: %s: %s\n",filename,buf);
			warn(msg);
			continue;
		}
		strcat(msp->sseq, p);
		p = strtok(NULL, " \n");
		if (! p) {
			sprintf(msg, "parse error: %s: %s\n",filename,buf);
			warn(msg);
			continue;
		}
		pos2 = atoi(p);

		if (pos1 < pos2) {
			if (! msp->from2 || pos1 < msp->from2)
				msp->from2 = pos1;
			if (! msp->to2 || msp->to2 < pos2)
				msp->to2 = pos2;
		} else {
			if (! msp->from2 || pos2 < msp->from2)
				msp->from2 = pos2;
			if (! msp->to2 || msp->to2 < pos1)
				msp->to2 = pos1;
		}
	} else if (strncmp(buf, "  Database:", 11)==0) {
		scoreout(&MspSet);
		initmsp(&MspSet);
		msp = NULL;
	}
    }
    fclose(fp);
}

MSP *initmsp(MSPset *MspSet)
{
	MspSet->MspNum = 0;
	MspSet->Msp[0].qseq = &seqbuf1[0];
	MspSet->Msp[0].sseq = &seqbuf2[0];
	MspSet->Msp[0].qseq[0] = '\0';
	MspSet->Msp[0].sseq[0] = '\0';
	MspSet->Msp[0].score = MspSet->Msp[0].origscore = 0;
	MspSet->Msp[0].delflag = 0;
	return &(MspSet->Msp[0]);
}
cancel_currentMSP(MSPset *MspSet)
{
	if (MspSet->MspNum > 0) {
		MspSet->MspNum--;
	}
}
MSP *nextmsp(MSPset *MspSet)
{
	MSP *msp;
	int qbufp, sbufp;
	msp = &(MspSet->Msp[MspSet->MspNum]);

	if (MspSet->MspNum > 0) {
		MSP *prevmsp = &(MspSet->Msp[MspSet->MspNum-1]);
		msp->qseq = prevmsp->qseq + strlen(prevmsp->qseq) + 1;
		msp->sseq = prevmsp->sseq + strlen(prevmsp->sseq) + 1;
		if (msp->qseq - seqbuf1 > SEQBUFSIZ ||
			msp->sseq - seqbuf2 > SEQBUFSIZ) {
			fprintf(stderr, "seqbuf overflows\n");
		}
	} else {
		msp->qseq = &seqbuf1[0];
		msp->sseq = &seqbuf2[0];
	}
	msp->qseq[0] = '\0';
	msp->sseq[0] = '\0';
	msp->score = msp->origscore = 0;
	msp->from1 = msp->to1 = msp->from2 = msp->to2 = 0;
	msp->delflag = 0;
	if (++MspSet->MspNum >= MAXMSP) {
		die(stderr, "Too many MSPs\n");
	}
	return msp;
}

scoreout(MSPset *MspSet)
{
	score_check(MspSet);
	check_ovlp(MspSet);
	print_mspset(MspSet);
}
score_check(MSPset *MspSet)
{
	int i;
	for (i = 0; i < MspSet->MspNum; i++) {
		if (modified_eval) {
			eval_mod(MspSet, &(MspSet->Msp[i]));
		}
		if (recalc_score) {
			calc_score(&(MspSet->Msp[i]));
		}
		score_check0(MspSet->qname, MspSet->sname, &(MspSet->Msp[i]));
	}
}
score_check0(char *qname, char *sname, MSP *msp)
{
	double ratio;
	if (msp->eval > EVAL_CUT) {
		if (verbose) {
			fprintf(stderr, "SKIP: low evalue: %s,%s\n", qname, sname);
		}
		msp->delflag = 1;
		return 0;
	}
	if (msp->score < SCORE_CUT) {
		if (verbose) {
			fprintf(stderr, "SKIP: low score: %s,%s\n", qname, sname);
		}
		msp->delflag = 1;
		return 0;
	}
	if (ignore_self && strcmp(qname, sname)==0
			&& msp->from1 == msp->from2
			&& msp->to1 == msp->to2) {
		if (verbose) {
			fprintf(stderr, "SKIP: identical seq: %s,%s\n",
				qname,sname);
		}
		msp->delflag = 1;
		return 0;
	}
	if (ignore_second_for_allall && strcmp(qname,sname) > 0) {
		msp->delflag = 1;
		return 0;
	}
	ratio = calc_bias(msp);
	if (ratio <= BIAS_CUT) {
		if (verbose) {
			fprintf(stderr, "Ratio SKIP: %s:%s %lf\n",
				qname,sname,ratio);
		}
		msp->delflag = 1;
	}
}
print_mspset(MSPset *MspSet)
{
	int i;
	for (i = 0; i < MspSet->MspNum; i++) {
		if (! MspSet->Msp[i].delflag) {
			print_msp(MspSet->qname, MspSet->sname, &(MspSet->Msp[i]));
			if (do_not_output_alireg) {
				/* output only the first hit */
				break;
			}
		}
	}
}
print_msp(char *qname, char *sname, MSP *msp)
{
	if (do_not_output_alireg) {
		printf("%s %s %d %d %lg\n",
			qname,sname,msp->ident,msp->score,msp->eval);
	} else {
		printf("%s %s %d %d %d %d %d %d %lg\n",
			qname,sname,msp->from1,msp->to1,
			msp->from2,msp->to2,
			msp->ident,msp->score,msp->eval);
	}
}

double calc_bias (MSP *msp)
{
	register int i, j;
	int qsiz, ssiz;
	double ProbQ[AAnum],ProbS[AAnum];
	int CntQ[AAnum], CntS[AAnum], totcntQ, totcntS;
	double totcnt;
	double Sexp, ratio;
	double sum;
	double pseudoc = PSEUDOC;


	for (i = 0; i < AAnum; i++) {
		CntQ[i] = CntS[i] = 0;
	}
	totcntQ = totcntS = totcnt = 0;

	qsiz = strlen(msp->qseq);
	ssiz = strlen(msp->sseq);
	for (i = 0; i < qsiz; i++) {
		if (AAidx[msp->qseq[i]]<0 || AAidx[msp->sseq[i]]<0){
			continue;
		}
		if (isalpha(msp->qseq[i])) {
			CntQ[AAidx[msp->qseq[i]]]++;
			totcntQ++;
		}
		if (isalpha(msp->sseq[i])) {
			CntS[AAidx[msp->sseq[i]]]++;
			totcntS++;
		}
	}
	for (i = 0; i < AAnum; i++) {
		ProbQ[i] = (Prob[i] * pseudoc + CntQ[i])
				/ (totcntQ + pseudoc);
		ProbS[i] = (Prob[i] * pseudoc + CntS[i])
				/ (totcntS + pseudoc);
	}
	sum = 0;
	for (i = 0; i < AAnum; i++) {
		for (j = 0; j < AAnum; j++) {
			sum += ProbQ[i] * ProbS[j] * Mat[i][j];
		}
	}
	totcnt = (double)(totcntS + totcntQ) / 2;
	Sexp = sum * totcnt;
	ratio = (msp->score - Sexp) / (msp->score - totcnt * Mexp);
	return ratio;
}
calc_score(MSP *msp)
{
	int i;
	int a1, a2, score;
	int gap_flag = 0;
	int qsiz;

	score = 0;
	qsiz = strlen(msp->qseq);
	for (i = 0; i < qsiz; i++) {
		if (msp->qseq[i] == '-') {
			if (gap_flag != 1) {
				score -= GapOpen;
			}
			score -= GapExt;
			gap_flag = 1;
		} else if (msp->sseq[i] == '-') {
			if (gap_flag != -1) {
				score -= GapOpen;
			}
			score -= GapExt;
			gap_flag = -1;
		} else {
			a1 = AAidx2[msp->qseq[i]];
			a2 = AAidx2[msp->sseq[i]];
			if (a1 < 0 || a2 < 0) {
				continue;
			}
			score += Mat[a1][a2];
			gap_flag = 0;
		}
	}
	msp->score = score;
}
cmpr_msp(const void *a, const void *b)
{
	const MSP *msp1 = a, *msp2 = b;
	return msp1->from1 - msp2->from1;
}
check_ovlp(MSPset *MspSet)
{
	int i, j, ov1, ov2;
	qsort(&(MspSet->Msp), MspSet->MspNum, sizeof(MSP), cmpr_msp);
	for (i = 0; i < MspSet->MspNum; i++) {
		if (MspSet->Msp[i].delflag) continue;
		for (j = 0; j < i; j++) {
			if (MspSet->Msp[j].delflag) continue;
			ovlp_msp(&(MspSet->Msp[i]), &(MspSet->Msp[j]));
		}
	}
}
ovlp_msp(MSP *msp1, MSP *msp2)
{
	int diff1, diff2, diag, gap;
	diff1 = msp2->from1 - msp1->to1;
	diff2 = msp2->from2 - msp1->to2;
	diag = abs(diff1 - diff2);
	gap = diff1 < diff2 ? diff1 : diff2;
	if (gap < 0) {
/*
		printf("%d,%d,%d,%d\n",
			msp1->from1,msp1->to1,msp1->from2,msp1->to2);
		printf("%d,%d,%d,%d\n",
			msp2->from1,msp2->to1,msp2->from2,msp2->to2);
*/
	}
}
ovlp_reg(int from1, int to1, int from2, int to2)
{
	int minov;
	int len1 = to1 - from1 + 1;
	int len2 = to2 - from2 + 2;
	minov = (len1 < len2 ? len1 : len2) * OvlpRatio;
	if (to1 < from2 + minov) {
		return 0;
	} else if (to2 < from1 + minov) {
		return 0;
	} else {
		return 1;
	}
}
eval_mod(MSPset *MspSet, MSP *msp)
{
	double orig_eval = msp->eval;
	msp->eval = MspSet->qlen * MspSet->slen * msp->eval / 100000;
	if (verbose) {
		fprintf(stderr, 
			"Modified Evalue: %lg -> %lg\n",orig_eval, msp->eval);
	}
}

calc_prob()
{
	int i;
	double sum = 0.0;
	for (i= 0; i < AAnum; i++){
		sum += Prob[i];
	}
	for (i= 0; i < AAnum; i++){
		Prob[i] /= sum;
	}
	for (i = 0; i < MAXCHR; i++) {
		AAidx[i] = -1;
		AAidx2[i] = -1;
	}
	for (i = 0; i < AAnum; i++) {
		AAidx[AA[i]] = i;
	}
	for (i = 0; i < AAnum2; i++) {
		AAidx2[AA2[i]] = i;
	}
}

read_mat(char *matfile)
{
	FILE *fp;
	int i, j;
	char *p;
	char *matdir;
	if ((fp = fopen(matfile, "r")) == NULL) {
		matdir = getenv("BLASTMAT");
		sprintf(matfile, "%s%s", matdir, matfile);
		if ((fp = fopen(matfile, "r")) == NULL) {
			fprintf(stderr, "Can't open matrix file %s\n", matfile);
			exit(1);
		}
	}
	i = 0;
	while (fgets(buf, BUFSIZ, fp) != NULL) {
		if (isalpha(buf[0])) {
			p = strtok(&buf[1], " ");
			j = 0;
			while (p && *p) {
				Mat[i][j] = atoi(p);
				p = strtok(NULL, " ");
				j++;
			}
			i++;
		}
	}
	Mexp = 0;
	for (i = 0; i < AAnum; i++) {
		for (j = 0; j < AAnum; j++) {
			Mexp += Mat[i][j] * Prob[i] * Prob[j];
		}
	}
}

die(char *msg)
{
	fprintf(stderr, msg);
	exit(1);
}
warn(char *msg)
{
	fprintf(stderr, msg);
}

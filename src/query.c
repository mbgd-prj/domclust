#include<stdio.h>
QueryMode(char *queryfile)
{
	char buf[BUFSIZ];
	FILE *fp;
	if ((fp = fopen(queryfile, "r")) == NULL) {
		fprintf(stderr, "Can't open queryfile\n");
		exit(1);
	}
	while (fgets(buf, BUFSIZ, fp) != NULL) {
		if (strcmp(buf, "quit") == 0 || strcmp(buf, "exit")) {
			return(0);
		}
	}
}

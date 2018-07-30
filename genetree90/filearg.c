#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define FARGS 50
#define FARGS_SIZE 20

/*  Open a file specified by @filename or @ filename in argv,
	and read the string tokens into fargv. Return fargc, the number
	of tokens. Call with in_file=NULL.
	Alternatively read from a file in_file.
*/

int get_file_arg(char *in_file, int argc, char **argv, char **fargv) {
	int i,j=0;
	char fname[100],buffer[100];
	FILE *in;
	fname[0]='\000';
	if(in_file!=NULL) strcpy(fname,in_file);
	else for(i=0;i<argc;i++) {
		if(argv[i][0]=='@'&&argv[i][1]=='\000') {
		strcpy(fname,argv[i+1]);
		break;
		}
		if(argv[i][0]=='@'&&argv[i][1]!='\000') {
		strcpy(fname,argv[i]+1);
		break;
		}
	}
	if(fname[0]=='\000') return 0;
	in=fopen(fname,"r");
	if(in==NULL) {
		printf("\n\tCan't open %s for fileargs\n",fname);
		exit(1);
	}
	while(feof(in)==0&&fscanf(in,"%s",buffer)==1) {
		if(j==FARGS) {
			printf("\n\tToo many entries in %s\n",fname);
			exit(1);
		}
		if(strlen(buffer)> FARGS_SIZE) {
			printf("\n\tfargv %d in %s exceeds length limit %d\n",
				j,fname,FARGS_SIZE);
			exit(1);
		}
		fargv[j][FARGS_SIZE]='\000'; /* Make sure! */
		strncpy(fargv[j],buffer,FARGS_SIZE);
		j++;
	}
	return j;
}


#ifdef TEST_FARGS
main(int argc, char **argv) {
	char **fargv;
	int i,fargc;
	fargv=(char**)malloc(FARGS*sizeof(char*));
	for(i=0;i<FARGS;i++)
		fargv[i]=(char*)malloc((FARGS_SIZE+1)*sizeof(char));
  fargc=get_file_arg(NULL,argc,argv,fargv);
	printf("\n\t");
	for(i=0;i<fargc;i++) printf("%s ",fargv[i]);
	printf("\n");
}
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
long mdum=(-1); /* random number seed */
#ifdef RAN2
float ran2(long *idum);  /* Numerical recipes long period generator */
#define RAND01 ((double)ran2(&mdum))
#else
double uni01(long *idum);
#define RAND01 (uni01(&mdum))
#endif
/* root_name is like myfile.# , where hash is to be replaced by
 0,1, ..., n_myfile - 1
 */

/* Get the number of files */
int n_files(char* root_name, char msym) {
	char *p,tmp_name[200];
	int count=0,flag=1;
	FILE *in;
	strcpy(tmp_name,root_name);
	p=strchr(tmp_name,msym);
	if(p==NULL) return -1; /* Error */
	while(flag==1) {
		*p='\000';
		sprintf(p,"%d",count);
		in=fopen(tmp_name,"r");
		flag=(in!=NULL);
		if(flag==1) { count++; fclose(in); }
	}
	return count;
}
		
/* Replace # successively by 0, 1, ..., n
	 first - string preceeding argv's (NULL if none)
	 last - string following argv's   (NULL if none)
	 start - start with argv[start]
	 flag - call with flag=0 to just print out commands
	 randarg - random number substitution for this arg
*/
void multiple_command(int argc, char **argv, char msym, int randarg,
	char *first, char *last, int n, int start, int flag) {
	char c_line[500],*p;
	int i,j;
	/* Random number substitution */
	if(randarg>=0&&randarg<argc)
	mdum=(-atol(argv[randarg]));
	for(i=0;i<n;i++) {
		for(j=0;j<500;j++) c_line[j]='\000';
		if(first!=NULL) strcpy(c_line,first);
		for(j=start;j<argc;j++) {
			if(j==randarg) {
				RAND01;
				sprintf(c_line+strlen(c_line),"%ld",mdum);
			}
			else strcat(c_line,argv[j]); 
			p=strchr(c_line,msym);
			if(p!=NULL) {
				*p='\000';
				sprintf(p,"%d",i);
			}
		if(j!=argc-1) strcat(c_line," ");
		}
		if(last!=NULL) strcat(c_line,last);
		if(flag==0) printf("%s\n",c_line); else system(c_line);
	}
}
			
#ifdef TESTM
main(int argc, char **argv) {
	multiple_command(argc,argv,'#',NULL,NULL,10,1);
	}
#endif

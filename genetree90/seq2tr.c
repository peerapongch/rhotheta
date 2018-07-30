/* Revised version 13/8/97 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <memory.h>


int **cmp_seq,n__seq,n0__seq=0,n__sites, max_str=5;
char *ancestor_seq[2];

/* Order columns of sequences, thinking of them as 
binary numbers.
*/

int b_cmp(const void *a, const void *b) {
	int i,flag=1,*pj,*pk;
	pj=(int*)a;pk=(int*)b;
	for(i=0;i<n0__seq;i++) {
		flag &= cmp_seq[i][*pj]==cmp_seq[i][*pk];
		if(flag==0) break;
	}
	if(flag==0&&cmp_seq[i][*pj]==1) return 1;
	else if(flag==0&&cmp_seq[i][*pj]==0) return -1;
	else return 1 - 2*(*pk < *pj);
}


int *site_sort(int *s) {
	int j;
	for(j=0;j<n__sites;j++) s[j]=j;
	qsort((void*)s,n__sites,sizeof(int),b_cmp);
return s;
}

void make_tree(int **path, int *s) {
	int i,j,*count;
	count=(int*)malloc(n0__seq*sizeof(int));
	for(j=0;j<n0__seq;j++) count[j]=0;
	for(j=0;j<n__sites;j++) for(i=0;i<n0__seq;i++)
		if(cmp_seq[i][s[j]]==1) path[i][count[i]++]=s[j]+1;
	for(i=0;i<n0__seq;i++) path[i][count[i]]=0;
	free(count);
}

/* return -1 not compatable in broard sense, 0 if not in narrow, 1 if ok */
int check_seq(int **q) {
	int u[2][2],i,j,k,l,x,ok=1;
	for(i=0;i<n__sites;i++) for(j=0;j<i;j++) {
		for(k=0;k<2;k++) for(l=0;l<2;l++) u[k][l]=0;
		for(x=0;x<n__seq;x++) {
			for(k=0;k<2;k++) for(l=0;l<2;l++)  
				u[k][l] += (q[x][i]==k && q[x][j]==l);
		}
		if(u[0][1] != 0&&u[1][0] != 0&&u[1][1] != 0) ok=0;
		if(u[0][1] != 0&&u[1][0] != 0&&u[1][1] != 0&&u[0][0] != 0) return (-1);
	}
	if(ok==0) return 0;
	return 1;
}

void compat_matrix(int **seq, int ** compat_matrix) {
	int u[2][2],i,j,k,l,x;
	for(i=0;i<n__sites;i++) for(j=0;j<i;j++) {
		compat_matrix[i][j]=0;
		for(k=0;k<2;k++) for(l=0;l<2;l++) u[k][l]=0;
		for(x=0;x<n__seq;x++) for(k=0;k<2;k++) for(l=0;l<2;l++) 
			u[k][l] += seq[x][i]==k && seq[x][j]==l;
		if(u[0][1] != 0&&u[1][0] != 0&&u[1][1] != 0) compat_matrix[i][j]=1;
		if(u[0][1] != 0&&u[1][0] != 0&&u[1][1] != 0&&u[0][0] != 0) 
			compat_matrix[j][i]=1;
	}
}

void tree_out(FILE *out, int **path, int row,int col, char *name, 
									int* seq_identity) {
	FILE *in;
	char buffer[100],*p,c;
	int i,j,k,count=0,check[3];
	if((in=fopen(name,"r"))==NULL) {
		printf("\n\tCannot open sequence file\n");
		exit(1);
	}
	for(i=0;i<row;i++) {
		c=' ';
		while(c!=':') {
			fscanf(in,"%c",&c);
			if(c!= '\n') fprintf(out,"%c",c);
		}
		for(j=0;path[seq_identity[i]][j];j++) 
			//edited here to handle H_n up to 9999
			if (path[seq_identity[i]][j]>=1000) fprintf(out,"%5d",path[seq_identity[i]][j]);
			else if (path[seq_identity[i]][j]>=100) fprintf(out,"%4d",path[seq_identity[i]][j]);
			else fprintf(out,"%3d",path[seq_identity[i]][j]);
			//fprintf(out,"%3d",path[seq_identity[i]][j]);
		fprintf(out,"%3d\n",0);
		for(j=0;j<col;j++) {
			c=' ';
			while(isspace(c)!=0) fscanf(in,"%c",&c);
		}
	}
	if(out!=stdout) fclose(out);
}

void get_size(int *prow,int *pcol,char *name) {
	FILE *in;
	char buffer[100],c;
	int r,x,count=0,bug=0;
	if((in=fopen(name,"r"))==NULL) {
		printf("\n\tCannot open sequence file\n");
		exit(1);
	}
	*prow=0;*pcol=0;
	while(fgets(buffer,100,in) != NULL) 
		if(strchr(buffer,':')!=NULL) (*prow)++;
	rewind(in);
	buffer[0]='\000';
	while(buffer[0]!=':' && count <max_str)  {
		fscanf(in,"%s",buffer);
		count++;
	}
	while(c!='\n' && bug <10000) {
		r=fscanf(in,"%c",&c);
		if(r==1&&isspace(c)==0) (*pcol)++;
		if(feof(in) != 0) break;
		bug++;
	}
	fclose(in);
	if(*prow <=0||*pcol<=0) {
		printf("\n\tProblem with data input file\n");
		exit(1);
	}
}
	
/*
void get_seq(int **q,int row,int col,char *name) {
	FILE *in;
	char buffer[100],c;
	int i,j,k,count=0,check[3],quit=0;
	if((in=fopen(name,"r"))==NULL) {
		printf("\n\tCannot open sequence file\n");
		exit(1);
	}
	for(i=0;i<row;i++) {
		buffer[0]='\000';
		count=0;
		while(buffer[0]!=':'&&count < max_str) {
			fscanf(in,"%s",buffer);
			count++;
		}
		for(j=0;j<col;j++) { 
			fscanf(in,"%d",&q[i][j]);
		}
	}
	for(j=0;j<col;j++) {
		check[0]=0;check[1]=0;
		for(i=0;i<row;i++) {
			check[0] += q[i][j]==0;
			check[1] += q[i][j]==1;
		}
		if(row>1 &&(check[0]==row||check[1]==row)) {
			printf("\n\tSite %d must be a segregating site",j+1);
			quit=1;
		}
		if((check[0]+check[1]) != row) {
			printf("\n\tSite %d must contain 0 or 1 data entries",j+1);
			quit=1;
		}
	}
	if(quit==1) {
		printf("\n\n");exit(1);
	}
	fclose(in);
}

*/

void get_seq(int **q,int row,int col,char *name, char *aname) {
	FILE *in;
	FILE* ain;
	char buffer[100],c;
	int i,j,k,count=0,check[3],quit=0;
	ancestor_seq[0]=(char*)malloc(col*sizeof(char));
	ancestor_seq[1]=(char*)malloc(col*sizeof(char));
	if(aname==NULL) for(j=0;j<col;j++) {
		ancestor_seq[0][j]='0';
		ancestor_seq[1][j]='$';
	}
	else {
		if((ain=fopen(aname,"r"))==NULL) {
			printf("\n\tCannot open ancestor file\n");
			exit(1);
		}
		for(j=0;j<col;j++) { 
			ancestor_seq[1][j]='$';
			c=' ';
			while(isspace(c)!=0) fscanf(ain,"%c",&c);
			ancestor_seq[0][j]=c;
		}
	}
	if((in=fopen(name,"r"))==NULL) {
		printf("\n\tCannot open sequence file\n");
		exit(1);
	}
	for(i=0;i<row;i++) {
		buffer[0]='\000';
		count=0;
		while(buffer[0]!=':'&&count < max_str) {
			fscanf(in,"%s",buffer);
			count++;
		}
		for(j=0;j<col;j++) { 
			c=' ';
			while(isspace(c)!=0) fscanf(in,"%c",&c);
			if(c==ancestor_seq[0][j]) q[i][j]=0;
			else if(ancestor_seq[1][j]=='$') { 
				ancestor_seq[1][j]=c; q[i][j]=1;
			}
			else if(c==ancestor_seq[1][j]) q[i][j]=1;
			else {
				printf("\n\tProblem with site types in sequence file");
				printf("\n\tsite %d, file %c, ancestors %c %c\n",
				j+1,c, ancestor_seq[0][j],ancestor_seq[1][j]);
				exit(1);
			}
		}
	}
	/* End read of file */
	for(j=0;j<col;j++) {
		check[0]=0;check[1]=0;
		for(i=0;i<row;i++) {
			check[0] += q[i][j]==0;
			check[1] += q[i][j]==1;
		}
		if(row>1 &&(check[0]==row||check[1]==row)) {
			printf("\n\tSite %d must be a segregating site",j+1);
			quit=1;
		}
	}
	if(quit==1) {
		printf("\n\n");exit(1);
	}
	fclose(in);
}

void recode(int **seq) {
	int i,j,u[2];
	for(j=0;j<n__sites;j++) {
		u[0]=0;u[1]=0;
		for(i=0;i<n__seq;i++) u[seq[i][j]]++;
		if(u[1]>u[0]) 
			for(i=0;i<n__seq;i++) seq[i][j]=1-seq[i][j];
	}
}

void do_tree(char *filename, char *outfile, char *aname) {
	int **q,*s,i,j,**path,**compat,ok,*seq_identity;
	FILE *out;
	get_size(&n__seq,&n__sites,filename);
	q=(int**)malloc(n__seq*sizeof(int*));
	for(i=0;i<n__seq;i++)
		q[i]=(int*)malloc(n__sites*sizeof(int));
	get_seq(q,n__seq,n__sites,filename,aname);
	ok=check_seq(q);
	if(ok==(-1)) {
		printf("\n\tSequences are incompatable in broard sense");
		compat=(int**)malloc(n__sites*sizeof(int*));
		for(i=0;i<n__sites;i++)
			compat[i]=(int*)malloc(n__sites*sizeof(int));
		compat_matrix(q,compat);
		printf("\n\tCompatability matrix (& broard, * narrow, . ok)\n");
		printf("\t  ");
		for(j=0;j<n__sites;j++) printf("%d",(j+1)%10);
			for(i=0;i<n__sites;i++) {
				printf("\n\t%d ",(i+1)%10);
				for(j=0;j<n__sites;j++) {
					if(compat[i][j]==1&&i<j) printf("&");
					else if(compat[i][j]==1&&i>j) printf("*");
					else printf(".");
				}
			}
			printf("\n");
			exit(1);
			}
			if(ok==0) {
				printf("\n\tSequences are incompatable in narrow sense,");
				printf("\n\tusing ancestral site (0) most common matrix.\n");
				recode(q);
			}
			s=(int*)malloc(n__sites*sizeof(int));
			path=(int**)malloc(n__seq*sizeof(int*));
			seq_identity=(int*)malloc(n__seq*sizeof(int));
			cmp_seq=(int**)malloc(n__seq*sizeof(int*));
			for(i=0;i<n__seq;i++) {
				path[i]=(int*)malloc((n__sites+1)*sizeof(int));
				seq_identity[i]=(-1);
				for(j=0;j<=n__sites;j++) path[i][j]=0;
			}
			/* Find duplicate sequences */
			for(i=0;i<n__seq;i++) {
				if(seq_identity[i]==(-1)) {
					cmp_seq[n0__seq]=q[i];
					seq_identity[i]=n0__seq;
					for(j=i+1;j<n__seq;j++) 
						if(seq_identity[j]==(-1)&&
							memcmp(q[i],q[j],n__sites*sizeof(int))==0) 
								seq_identity[j]=n0__seq;
				n0__seq++;
				}
			}
			site_sort(s);
			/* Binary sort order of columns of seqs, debug only
			printf("\t");
			for(i=0;i<n__sites;i++) printf("%d ",s[i]+1);
			printf("\n");
			*/
			make_tree(path,s);
			if(strcmp(outfile,"stdout")==0)
				tree_out(stdout,path,n__seq,n__sites,filename,seq_identity);
			else {
				out=fopen(outfile,"w");
				if(out==NULL) {
					printf("\n\tCan't open tree output file %s\n",outfile);
					exit(1);
				}
			tree_out(out,path,n__seq,n__sites,filename,seq_identity);
			}
			free(cmp_seq);
}

#ifdef SEQTOTR
char *example1=
"  0 6 : 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n\
  0 1 : 0 0 0 0 1 0 0 0 0 0 0 0 0 0\n\
  0 1 : 0 1 0 0 0 0 0 0 0 0 0 0 0 0\n\
  0 3 : 0 0 1 0 0 0 0 0 0 0 0 0 0 0\n\
  0 1 : 0 0 0 1 0 0 0 0 0 0 0 0 0 0\n\
  0 1 : 0 0 0 0 0 0 0 0 0 0 0 1 0 0\n\
  0 1 : 1 0 0 0 0 0 1 0 0 0 0 1 0 0\n\
  0 1 : 0 0 0 0 0 0 0 0 0 1 0 0 0 0\n\
  0 1 : 1 0 0 0 0 0 0 1 0 0 0 1 0 0\n\
  0 1 : 0 0 0 0 0 0 0 0 1 0 0 0 0 0\n\
  0 1 : 0 0 0 0 0 0 0 0 0 0 0 0 0 1\n\
  0 1 : 1 0 0 0 0 0 0 0 0 0 1 1 0 0\n\
  0 1 : 0 0 0 0 0 0 0 0 0 0 0 0 1 0\n\
  1 7 : 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n\
  1 3 : 0 0 0 0 1 0 0 0 0 0 0 0 0 0\n\
  1 4 : 1 0 0 0 0 0 0 0 0 0 0 1 0 0\n\
  1 2 : 0 0 0 1 0 0 0 0 0 0 0 0 0 0\n\
  1 1 : 0 0 0 1 0 1 0 0 0 0 0 0 0 0\n\
\n";

char *example2=
"  0 6 :  A . G Z\n\
  1 1 :  A . G X\n\
  1 1 :  T . G Z\n\
  2 3 :  A C G Z\n\
  2 1 :  A . C Z\n\
  2 1 :  A . G Z\n\
\n";

char *ancestor="        A . G Z\n";

void help() {
	printf("  Example 1\n");
	printf("  sequence file named example1.seq\n");
	printf(example1);
	printf("  subpopulations 0 and 1 on left are optional\n");
	printf("  0 ancestor base, 1 is mutant base\n");
	printf("  seq2tr example1.seq example1.tre\n\n");
	printf("  Example 2\n");
	printf("  sequence file named example2.seq\n");
	printf(example2);
	printf("  Spaces between bases are optional\n");
	printf("  ancestor base file named ancestor.bas\n");
	printf(ancestor);
	printf("  seq2tr example2.seq example2.tre ancestor.bas\n\n");
	exit(1);
}


char *usage="\n\tseq2tr seq_file tree_out_file [ancestor_base_file]\n";
char *version="\tVer 1.1 13/8/97, Bob Griffiths\n";
char *helpcmd="\tseq2tr -h for help\n\n";
main(int argc, char**argv) {
	if(argc>1&&(argv[1][1]=='h'||argv[1][1]=='h')) help();
	if(argc<3) {
		printf(usage);
		printf(version);
		printf(helpcmd);
		exit(1);
		}
	if(argc==3) do_tree(argv[1],argv[2],NULL);
	else if(argc==4) do_tree(argv[1],argv[2],argv[3]);
}
#endif

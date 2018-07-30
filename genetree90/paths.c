#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
	int n; /* total number of sequences */
	int types; /* number of distinct sequences */
	int *subpop; /* individuals in subpopulations */
	int **sequence; /* sequences as paths to the root */
	int *multiplicity; /* multiplicity of sequences */
	int singletons; /* number of singletons */
	int *singleton_list;  /* list of sequences which have singleton first
																coordinates */
	int *singleton_exterior; /* list of seq whose multiplicity is increased
		if a singleton is deleted (or -1 if no such sequence).
				*/
	int  *color; /* Subdivided populations that the genes are in */
} TREE;


int get_sites(TREE *t);

/* Return the sequence that site k is in. */
int find_path(int **p2root,int nseqs, int k) {
	int i,j;
	for(i=0;i<nseqs;i++) for(j=0;p2root[i][j] != 0;j++)
		if(p2root[i][j] == k) return i;
	return -1; /* Not found */
}

/* Is site k in sequence i */
int in_sequence(int **p2root,int i, int k) {
	int j;
	 for(j=0;p2root[i][j] != 0;j++) if(p2root[i][j] == k) return 1;
	 return 0;
}

char *int2str(char *buffer, int x) {
	sprintf(buffer,"%d",x);
	return buffer;
}

/* New tree with root to the leaf side of k */
/*
char mainf[30]="tree.";
*/
int newtree(TREE* t,int k, char *mainf) {
 int i,j,a,b,**p2root,nseqs;
 char ext[4],name[50];
 FILE *out;
 strcpy(name,mainf);
 int2str(ext,k);
 strcat(name,ext);
 p2root=t->sequence;
 if(k==0) { /* No change to tree */
 	out=fopen(name,"w");
 	for(i=0;i<t->types;i++) {
		fprintf(out,"\t%3d ",t->color[i]);
		fprintf(out,"%3d : ",t->multiplicity[i]);
			for(j=0;p2root[i][j]!=0;j++)
				fprintf(out,"%d ",p2root[i][j]);
		fprintf(out,"%d\n",0);
	}
		fclose(out);
		return 0;
 }
 nseqs=t->types;
 a=find_path(p2root,nseqs,k);
 /*
 if(a==(-1)) return -1;
 This will skip missing site labels in the tree, but still
 problems with multiple files without consecutive labeling
 */
 if(a==(-1)) {
 	printf("\n\tSite %d not found in tree, input problem ?\n",k);
	exit(1);
	}
 out=fopen(name,"w");
 for(i=0;i<nseqs;i++) {
	fprintf(out,"\t%3d ",t->color[i]);
	fprintf(out,"%3d : ",t->multiplicity[i]);
	if(in_sequence(p2root,i,k)==1) {
		for(j=0;p2root[i][j]!=k;j++) 
			fprintf(out,"%d ",p2root[i][j]);
		fprintf(out,"%d\n",0);
	}
	else {
		for(j=0;p2root[i][j]!=0;j++)
			if(in_sequence(p2root,a,p2root[i][j])==0) 
				fprintf(out,"%d ",p2root[i][j]);
		for(j=0;p2root[a][j]!=0;j++);	
		while(p2root[a][j]!=k) {
			j--;
			if(in_sequence(p2root,i,p2root[a][j])==0)
				fprintf(out,"%d ",p2root[a][j]);
		}
		fprintf(out,"%d\n",0);
	}
	}
	fclose(out);
	return 1;
}


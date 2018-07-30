#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
double **lu(double **a, int n, double **lum);
/* Solve LU vx[]=vb [] */
double *lu_solve(double **lum, int n, double *vb, double *vx);
/* g subpopulations, Migration matrix M */
#define c(i,j,k,l) c[i*g-i*(i-1)/2+j-i][k*g-k*(k-1)/2+l-k]
#define fs(i,j) fs[i*g-i*(i-1)/2+j-i]
#define fm1(i,j) fm1[i*g-i*(i-1)/2+j-i]
/* Triangular index scheme with i <= j */

void *gf_malloc(int bytes) {
	void *p;
	p=(void*)malloc(bytes);
	if(p==NULL) {
		printf("\n\tProblem allocating memory in gf_malloc()\n");
		exit(1);
	}
	return p;
}

double **get_matrix(int g, double *n, double theta, 
	double **m, double **c) {
	int i,j,k,l;
	double *q,scale;
	q=(double*)gf_malloc(g*sizeof(double));
	for(i=0;i<g;i++) {
		q[i]=0.0;
		for(j=0;j<g;j++) if(i!=j) q[i] += m[i][j];
	}
	for(i=0;i<g;i++) for(j=i;j<g;j++) {
		scale=q[i]+q[j]+theta;
		if(i==j) scale += 1/n[i];
		for(k=0;k<g;k++) for(l=k;l<g;l++) {
			if(i==k&&j==l) 
				c(i,j,k,l)= -scale;
			else if(i==k&&l>=k) 
				c(i,j,k,l)=m[j][l];
			else if(k<l&&l==i)
				c(i,j,k,l)=m[j][k];
			else if(j==k&&l>=k)
				c(i,j,k,l)=m[i][l];
			else if(k<l&&l==j)
				c(i,j,k,l)=m[i][k];
			else
				c(i,j,k,l)=0.0;
		}
		for(k=0;k<g;k++) for(l=k;l<g;l++) {
				c(i,j,k,l) /=scale;
				if(i==j&&!(k==i&&j==l)) c(i,j,k,l) *= 2.0;
		}
	}
	/* Debug 
	for(i=0;i<g;i++) for(j=i;j<g;j++)  
		for(k=0;k<g;k++) for(l=k;l<g;l++) 
		printf("\n %d %d %d %d %.4lf",i,j,k,l,c(i,j,k,l));
	printf("\n");
	*/
	free(q);
	return c;
}

double *get_r_matrix(int g, double *n, double theta, double **m,
	double *fm1, double *fs, int s) {
	int i,j;
	double *q,scale;
	q=(double*)gf_malloc(g*sizeof(double));
	for(i=0;i<g;i++) {
		q[i]=0.0;
		for(j=0;j<g;j++) if(i!=j) q[i] += m[i][j];
	}
	for(i=0;i<g;i++) for(j=i;j<g;j++) {
		scale=q[i]+q[j]+theta;
		if(i==j) scale += 1/n[i];
		if(s>0) fs(i,j)= -theta*fm1(i,j)/scale;
		if(s==0) {
			if(i==j) {
				fs(i,j)= -1.0/(n[i]*scale);
			}
			else fs(i,j)=0.0;
		}
		/* Debug
		printf("\n %d %d %.4lf",i,j,fs(i,j));
		*/
	}
	free(q);
	return fs;
}

void get_gf(int s, double **gf, int g, double *n, double theta, double **m) {
	double *fm1, *fs,**c, **clu,**a,**b; 
	int i,j,k,l;
	fs=(double*)gf_malloc(sizeof(double)*g*(g+1)/2);
	c=(double**)gf_malloc(sizeof(double*)*g*(g+1)/2);
	for(i=0;i<g*(g+1)/2;i++) 
		c[i]=(double*)gf_malloc(sizeof(double)*g*(g+1)/2);
	clu=(double**)gf_malloc(sizeof(double*)*g*(g+1)/2);
	for(i=0;i<g*(g+1)/2;i++) 
		clu[i]=(double*)gf_malloc(sizeof(double)*g*(g+1)/2);
	get_matrix(g,n,theta,m,c);
	lu(c,g*(g+1)/2,clu);
	/***  Debug
	a=(double**)gf_malloc(sizeof(double*)*g*(g+1)/2);
	b=clu;
	for(i=0;i<g*(g+1)/2;i++) 
		a[i]=(double*)gf_malloc(sizeof(double)*g*(g+1)/2);
	for(i=0;i<g*(g+1)/2;i++) for(j=0;j<g*(g+1)/2;j++) {
		a[i][j]=0.0;
		for(k=0;!(k>i||k>j);k++) {
			if(i!=k) a[i][j] += b[i][k]*b[k][j];
			else a[i][j] += b[i][j];
		}
	}
	printf("\nC matrix\n");
	for(i=0;i<g*(g+1)/2;i++) {
		printf("\n\t");
		for(j=0;j<g*(g+1)/2;j++) printf("%10.5lf ",c[i][j]);
	}
	printf("\nLU reconstructed\n");
	for(i=0;i<g*(g+1)/2;i++) {
		printf("\n\t");
		for(j=0;j<g*(g+1)/2;j++) printf("%10.5lf ",a[i][j]);
	}
	printf("\nLU residual\n");
	for(i=0;i<g*(g+1)/2;i++) {
		printf("\n\t");
		for(j=0;j<g*(g+1)/2;j++) printf("%10.5lf ",c[i][j]-a[i][j]);
	}
	printf("\n\n");
	****/
	fm1=NULL; /* k=0 */
	for(k=0;k<=s;k++) {
		get_r_matrix(g,n,theta,m,fm1,fs,k);
		lu_solve(clu,g*(g+1)/2,fs,gf[k]);
		fm1=gf[k];
	}
}
	
#define gen(s,i,j) gf[s][i*g-i*(i-1)/2+j-i]

#ifdef TWOGENES
char *usage="\n\ttwogenes n_subpop theta mig_file sites";
char *ver="\n\tVer 0.0 Bob Griffiths 11/10/96";
main(int argc, char **argv) {
	double *n, theta=2.0, **m, **gf, sum=0.0; 
	int i,j,k,l,g=4,s=5;
	FILE *in;
	if(argc<5) {
		printf(usage);
		printf(ver);
		exit(1);
	}
	g=atoi(argv[1]);
	theta=atof(argv[2]);
	in=fopen(argv[3],"r");
	s=atoi(argv[4]);
	m=(double**)gf_malloc(sizeof(double*)*g);
	n=(double*)gf_malloc(sizeof(double)*g);
	gf=(double**)gf_malloc(sizeof(double*)*(s+1));
	for(i=0;i<g;i++) 
		m[i]=(double*)gf_malloc(sizeof(double)*g);
	for(i=0;i<=s;i++) 
		gf[i]=(double*)gf_malloc(sizeof(double)*g*(g+1)/2);
	for(i=0;i<g;i++) fscanf(in,"%lf",n+i);
	for(i=0;i<g;i++) for(j=0;j<g;j++) fscanf(in,"%lf",&m[i][j]);
	/*
	n[0]=2.0;n[1]=2.0;n[2]=2.0;n[3]=2.0;
	for(i=0;i<g;i++) m[i][i]=0.0;
	for(i=0;i<g;i++) for(j=0;j<g;j++)  if(i!=j) m[i][j]=0.5;
	*/
	fclose(in);
	get_gf(s,gf,g,n,theta,m);
	for(i=0;i<g;i++) for(j=i;j<g;j++) {
		printf("\n\t%d %d\n\t",i,j);
		sum=0.0;
		for(k=0;k<=s;k++) {
			if(k!=0&&(k%10)==0) printf("\n\t");
			printf("%.4lf ",gen(k,i,j));
			sum += gen(k,i,j);
		}
		printf("%.4lf ",sum);
	}
}
#endif

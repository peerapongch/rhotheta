#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
/* LU matrix decomposition from Burden, Faires, p359,
implementation for non-singular matrices only.
Real indexing is from 0->n-1, but code is from 1->n.
Must have a(1,1) != 0.
L has diagonal elements 1.
*/
#define a(i,j) a[i-1][j-1]
#define l(i,j) lum[i-1][j-1]
#define u(i,j) lum[i-1][j-1]

double **lu(double **a, int n, double **lum) {
	int i,j,k;
	if(a(1,1)==0.0) {
		printf("\n\tCan't do LU\n");
		exit(1);
	}
	u(1,1)=a(1,1);
	for(j=2;j<=n;j++) {
		u(1,j)=a(1,j);
		l(j,1)=a(j,1)/u(1,1);
	}
	for(i=2;i<=n-1;i++) {
		u(i,i)=a(i,i);
		for(k=1;k<=i-1;k++) u(i,i) -= l(i,k)*u(k,i);
		if(u(i,i)==0.0) {
			printf("\n\tCan't do LU\n");
			exit(1);
		}
		for(j=i+1;j<=n;j++) {
			u(i,j)=a(i,j);l(j,i)=a(j,i);
			for(k=1;k<=i-1;k++) {
				u(i,j) -= l(i,k)*u(k,j);
				l(j,i) -= l(j,k)*u(k,i);
			}
			l(j,i) /= u(i,i);
		}
	}
	u(n,n) = a(n,n);
	for(k=1;k<=n-1;k++) u(n,n) -= l(n,k)*u(k,n);
	return lum;
}

#define vx(i) vx[i-1]
#define va(i) va[i-1]
#define vb(i) vb[i-1]

/* Solve LU vx[]=vb [] */

double *lu_solve(double **lum, int n, double *vb, double *vx) {
	int i,j,k;
	double x;
	vx(1)=vb(1);
	for(i=2;i<=n;i++) {
		vx(i)=vb(i);
		for(j=1;j<i;j++) vx(i) -= l(i,j)*vx(j);
	}
	vx(n) /= u(n,n);
	for(i=n-1;i>=1;i--) {
		x=vx(i);
		for(j=i+1;j<=n;j++) x -= u(i,j)*vx(j);
		x /= u(i,i);
		vx(i)=x;
	}
	return vx;
}

#ifdef TESTLU
main() {
	int i,j,k;
	double **a,**b,*va,*vb;
	a=(double**)malloc(5*sizeof(double*));
	b=(double**)malloc(5*sizeof(double*));
	va=(double*)malloc(5*sizeof(double));
	vb=(double*)malloc(5*sizeof(double));
	for(i=0;i<5;i++) {
		a[i]=(double*)malloc(5*sizeof(double));
		b[i]=(double*)malloc(5*sizeof(double));
	}
	for(i=0;i<5;i++) for(j=0;j<5;j++) {
		if(i==j) a[i][j]=5.0;
		if(i!=j) a[i][j]=2.0;
	}
	lu(a,5,b);
	for(i=0;i<5;i++) {
		printf("\n\t");
		for(j=0;j<5;j++) printf("%10.3lf ",b[i][j]);
	}

	for(i=0;i<5;i++) for(j=0;j<5;j++) {
		a[i][j]=0.0;
		for(k=0;!(k>i||k>j);k++) {
			if(i!=k) a[i][j] += b[i][k]*b[k][j];
			else a[i][j] += b[i][j];
		}
	}
	printf("\n\n");
	for(i=0;i<5;i++) {
		printf("\n\t");
		for(j=0;j<5;j++) printf("%10.3lf ",a[i][j]);
	}
	vb[0]=0.0;vb[1]=1.0;vb[2]=2.0;vb[3]=3.0;vb[4]=4.0;
	lu_solve(b,5,vb,va);
	printf("\n\n\t");
	for(j=0;j<5;j++) printf("%10.3lf ",va[j]);
}
#endif

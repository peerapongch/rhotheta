#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "factl.h"
extern int n_subpop;
extern int *subsize;
extern int growth;
extern double *growth_rate;
double integrate(double (*f)(double x, double *p),  double *p,
	double lower, double upper,
	double tolerance);
double lambda(int c, double tim);
double ilambda(int c, double tim0, double tim1);

/* Single island, constant population size */
double two_single(int a, int b, double theta) {
	if(a==0&&b==0) return 1.0/(1.0+theta);
	else return pow(theta/(1.0+theta),(double)(a+b))/(1.0+theta);
}

/* p[0] is number of seg sites, p[1] is tim, p[2] is theta */
double integrand(double x, double *p) {
	double z,y;
	int n;
	z=lambda(0,x+p[1])*exp(-p[2]*x)*exp(-ilambda(0,p[1],x+p[1]));
	n=(int)p[0];
	if(n>factdim) {
		printf("\n\tFactorial too big\n");
		exit(1);
	}
	if(n>0)
		y = z*pow(p[2]*x,p[0])/fact[n];
	else y=z;
	/*
		printf("\n\tx=%.4lf , z= %.4lf integ=%.4lf, time = %.4lf ",x,z,y,p[1]);
	*/
	return y;
}
		
/* p[0] is number of seg sites, p[1] is tim, p[2] is theta,
p[3] is magnification:  lambda(u) -> lambda(u*p[3]) */

double integrand_mag(double x, double *p) {
	double z,y;
	int n;
	z=lambda(0,p[3]*(x+p[1]))*exp(-p[2]*x)
		*exp(-ilambda(0,p[3]*p[1],p[3]*(x+p[1]))/p[3]);
	n=(int)p[0];
	if(n>factdim) {
		printf("\n\tFactorial too big\n");
		exit(1);
	}
	if(n>0)
		y = z*pow(p[2]*x,p[0])/fact[n];
	else y=z;
	/*
		printf("\n\tx=%.4lf , z= %.4lf integ=%.4lf, time = %.4lf ",x,z,y,p[1]);
	*/
	return y;
}

#define TEPS 1.0e-100
/* Single island, variable population size */
double two_single_var(int a, int b, double theta, double tim) {
	double p[3],z,tol,sum=0.0,x=0.0,dx=0.1,end_pt=5.0,inc=0.0,z1;
	int n=50,i;
	p[0]=(double)(a+b);p[1]=tim;p[2]=theta;
	/* Choose the end point, code for special schemes */
	if(growth==1)	 end_pt=log(5.0*growth_rate[0]+1.0)/growth_rate[0];
	dx=end_pt/n;
	/* Trapezoidal integration to get a rough idea to use for the
	relative tolerance */
	for(i=0;i<n;i++) {
		if(i==0||i==n-1) sum += 0.5*integrand(x,p);
		else sum += integrand(x,p);
		x += dx;
	}
	sum *= dx;
	/* End trapezoidal integration, sum is the integral approximation */
  tol=0.001*sum; /* Guess for relative accuracy */
	z=0.0;
	z1=integrate(integrand,p,0.0,0.1,tol);
	while(z1 > TEPS) {
		inc += 0.1;
		z += z1;
	 	z1=integrate(integrand,p,inc,inc+0.1,tol);
	}
	/*
	printf("\n\t sum= %.4e z = %.4e",sum,z);
	*/
	return z;
}

/* Single island, variable population size */
double two_single_var_mag(int a, int b, double theta,double tim,double mag) {
	double p[4],z,tol,sum=0.0,x=0.0,dx=0.01,end_pt=5.0,inc=0.0,z1;
	int n=50,i;
	p[0]=(double)(a+b);p[1]=tim;p[2]=theta;p[3]=mag;
	/* Choose the end point, code for special schemes */
	if(growth==1)	 end_pt=log(5.0*growth_rate[0]+1.0)/growth_rate[0];
	dx=end_pt/n;
	for(i=0;i<n;i++) {
		if(i==0||i==n-1) sum += 0.5*integrand(x,p);
		else sum += integrand_mag(x,p);
		x += dx;
	}
	sum *= dx;
  tol=0.001*sum; /* Guess for relative accuracy */
	z=0.0;
 	z1=integrate(integrand_mag,p,0.0,0.1,tol);
	while(z1 >TEPS) {
		inc += 0.1;
		z += z1;
 		z1=integrate(integrand_mag,p,inc,inc+0.1,tol);
	}
	/* Debug
	printf("\n\t XX mag = %.4lf, sum = %.4e, z = %.4e",mag,sum,z);
	*/
	return z;
}

	


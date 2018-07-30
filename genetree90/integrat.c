#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* Adapative quadrature integration
  function f(x,p), p is parameter list
	 N max number of steps
	 From Burden and Faires book p187
*/
#define N 100
#define f(z) (*f)(z,p)

double integrate(double (*f)(double x, double *p),  double *p,
	double lower, double upper,
	double tolerance) {
		double app,s1,s2,fd,fe,v[9],fa[N],fb[N],fc[N],h[N],s[N],a[N],tol[N],l[N];
		int i;
		app=0.0;
		i=1;
		tol[i]=10.0*tolerance;
		a[i]=lower;
		h[i]=(upper-lower)/2.0;
		fa[i]=f(lower);
		fc[i]=f(lower+h[i]);
		fb[i]=f(upper);
		s[i]=h[i]*(fa[i]+4.0*fc[i]+fb[i])/3.0;
		l[i]=1;
		while(i>0) {
			fd=f(a[i]+h[i]/2.0);
			fe=f(a[i]+3.0*h[i]/2.0);
			s1=h[i]*(fa[i]+4.0*fd+fc[i])/6.0;
			s2=h[i]*(fc[i]+4.0*fe+fb[i])/6.0;
			v[1]=a[i];
			v[2]=fa[i];
			v[3]=fc[i];
			v[4]=fb[i];
			v[5]=h[i];
			v[6]=tol[i];
			v[7]=s[i];
			v[8]=l[i];
			i--;
			if(fabs(s1+s2-v[7])<v[6]) 
				app=app+s1+s2;
				else {
					if(v[8]>=N) {
						printf("\n\tToo many steps (%d) in integrate\n",N);
						exit(1);
					}
					else {
						i++;
						a[i]=v[1]+v[5];
						fa[i]=v[3];
						fc[i]=fe;
						fb[i]=v[4];
						h[i]=v[5]/2.0;
						tol[i]=v[6]/2.0;
						s[i]=s2;
						l[i]=v[8]+1.0;
						i++;
						a[i]=v[1];
						fa[i]=v[2];
						fc[i]=fd;
						fb[i]=v[3];
						h[i]=h[i-1];
						tol[i]=tol[i-1];
						s[i]=s1;
						l[i]=l[i-1];
					}
				}
			}
			return app;
}


/*
double test(double x, double *p) {
	return 3.0*pow(3.0*x,4.0)*exp(- (*p)*x)/24.0;
}

double test(double x, double *p) {
	return exp(- (*p)*x);
}

main() {
	double z,a=3.0;
	z=integrate(test,&a,0.0,10.0,0.001);
	printf("\n\tIntegral is %.4lf\n",z);
}
*/
		

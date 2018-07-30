#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
extern unsigned _stklen=10480;
char* ver="9.0";
char *g_date="7/08/00";
int checktree=0;
#define MSYMBOL '#'
#define MFTEMP "genetree.tmp"
#define FARGS 50
#define FARGS_SIZE 20
#define THETAZERO 1.0e-30
FILE *mf_tmp;
/* Contains  xrep,sfnl,ssfnl,stmrca,sstmrca when multiple files */
long dum=(-1); /* random number seed */
#ifdef RAN2
float ran2(long *idum);  /* Numerical recipes long period generator */
#define RAND01 ((double)ran2(&dum))
#else
double uni01(long *idum);
#define RAND01 (uni01(&dum))
#endif
#define randint(n) (int)floor((n)*(RAND01))
int out_level=3;
int mf_out_level=100;
int n_files(char* root_name, char msym);
int multiple_files=0;
int n_m_files=0;
int check_m=0;
int two_stop=0;
int batch=0;
int notime=0;
int pairw=0;

int margin_flag=0;
int tmrca_distn_flag=0;
int mig_matrix_flag=0;
int age_distn_flag=0;
double *tmrca_distn;
double **age_distn;
int tmrca_cells=0;
double tmrca_top=5.0;
char distn_file[100];
double gene_pop=0.0,generation=0.0, c_scale=0.0;

void multiple_command(int argc, char **argv, char msym, int randarg,
	char *first, char *last, int n, int start, int flag);
int get_file_arg(char *in_file, int argc, char **argv, char **fargv);
double two_single(int a, int b, double theta);
double two_single_var_mag(int a, int b, double theta,double tim,double mag);
double two_single_var(int a, int b, double theta, double tim);
void get_gf(int s, double **gf, int g, double *n, double theta, double **m);
double **gf0; /* Indices are sites, subpops */
double ****s_gf; /* Indices are theta_s, m_s, sites, subpops */
#define pgf(s,i,j) gf0[s][i*n_subpop-i*(i-1)/2+j-i]
#define spgf(s,i,j,t_si,m_sj) s_gf[t_si][m_sj][s][i*n_subpop-i*(i-1)/2+j-i]
void do_tree(char *filename, char *tree_out, char *aname);

char tree_ps[100]="tree.ps";
char tree_gt[100]="tree.gt";
char tree_age[100]="tree.age";
/* Default is
 * treepic(name,age_name,pic_name,'p',12.0,14.0,17.0,25.0,0.8,2.4,1,1,c_scale);
 */
extern void treepic (char *name, char *age_name, char *pic_name, char g,
	double x_size,  double y_size, double x_paper,  double y_paper, 
	double l_width,double ps_r, int ax, int mu, double c_scale) ;

int debug=0;
int maxlength=30;
int maxsim=500;
int n_subpop=0; /* Number of subpopulations */
int set_n_subpop=0; /* Number of subpopulations */
int n_sites=0; /* Doesn't include 0 */
int n_sample=0; /* Total in the sample */
char tree_name[100];
char mainf[50];
double *subsize; /* Relative subsizes */
double *keep_subsize; /* Relative subsizes */
/* Age of sites */
double *m_age;
double *sage;
double *ssage;
double **mig_matrix;
double **csite; /* Color distribution of mutations */
int *m_csite;
int growth=0;
double *growth_rate; /* exponential growth rates */
double *t_rate; /* relative mutation rates */
int maxtype=30;
/* Surface variables */
int surface=0;
int tmrca_surface=0;
int t_surface=0;
int m_surface=0;
int g_surface=0;
double *surface_theta;
double *surface_m;
double *surface_g;
double t_s,m_s,g_s;
int theta_pts=0;
int m_pts=0;
int g_pts=0;
double ***surface_h; /* move */
double ***surface_fnl; /* 1st index is for theta, 2nd for m, for simulate */
double ***ssurface_fnl;  /* Replicate sum */
double ***sssurface_fnl; /* Replicate sumsq */
double ***tm_surface;  /* TMRCA surface mean */
double ***ttm_surface; /* TMRCA surface ss */
char tmrca_s_name[100];
int matrix_surface=0;
int use_new_mig=0;
int **move_migration;
/* migrations i->j in each simulation for a migration matrix surface */

char s_name[100]={'\000'};

double *move_prob;

int **anc_trace;
int *check_anc;
int *initial_color;
int *initial_subpop;
int *mrca_test;
double **mrca_c; /* Ancestor of sample colors i distribution in colors j*/
double *tmrca_c;
double *ss_tmrca_c;
int *m_mrca_c;
double *m_tmrca_c;

/* anc_trace[0][k]=i and anc_trace[1][k]=j
   means that individual k in the initial sample 
	 has current ancestor type i, and the ancestor is 
	 the jth individual in this group.
*/


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

int newtree(TREE* t,int k, char *mainname);

void *c_malloc(int bytes) {
	void *p;
	if(bytes ==0) {
		printf("\n\tTried to allocate zero bytes\n");
		exit(1);
	}
	if(bytes ==0) {
		printf("\n\tTried to allocate <0 bytes!\n");
		exit(1);
	}
	p=(void*)malloc(bytes);
	if(p==NULL) {
		printf("\n\tMemory allocation error, %d bytes",bytes);
		exit(1);
	}
	return p;
}

double sqrt_c(double x) {
	 if(x<0.0) {
		printf("\n\tTried to take sqrt(%2.4e)\n",x);
		exit(1);
	}
	if(fabs(x) < 1.0e-200) return 0.0;
	else return sqrt(x);
}

/* MAYBE KEEP THIS */
#define sqrt(x) sqrt_c((x))

double choose(int n, int j) {
	double x=1.0;
	int k,i;
	if(j>n) return 0.0;
	if(j<n-j) k=j; else k=n-j;
	if(k==0) return 1.0;
	if(k==1) return (double)n;
	if(k==2) return n*(n-1.0)/2.0;
	for(i=0;i<k;i++) 
		x *= (double)(n-i)/(i+1.0);
	return x;
}

/*  Surface for migration matrix */
char m_filein[100];
char m_fileout[100];
int mig_dim;
int *lexi;
double **new_mig;
int *m_index[2]; /* Surface indices */
double *m_end[2]; /* Surface endpts */
int *m_surface_pts;
int m_dim;

void read_migration(char *filename) {
	int count=0,scount,i,j,k;
	double x,y;
	FILE *m_in;
	m_in=fopen(filename,"r");
	if(m_in==NULL) {
		printf("\n\tCan't open migration parameter file %s\n",filename);
		exit(1);
	}
	while(feof(m_in)==0) {
		scount=fscanf(m_in,"%d %d %lf %lf %d",&i,&j,&x,&y,&k);
		if(scount!=5) {
			if(feof(m_in)!=0) break;
			else {
				printf("\n\tProblem reading migration parameter file %s\n",filename);
				exit(1);
			}
			break;
		}
		count++;
	}
	m_dim=count;
	rewind(m_in);
	m_index[0]=(int*)c_malloc(count*sizeof(int));
	m_index[1]=(int*)c_malloc(count*sizeof(int));
	m_surface_pts=(int*)c_malloc(count*sizeof(int));
	m_end[0]=(double*)c_malloc(count*sizeof(double));
	m_end[1]=(double*)c_malloc(count*sizeof(double));
	lexi=(int*)c_malloc(count*sizeof(int));
	new_mig=(double**)c_malloc(n_subpop*sizeof(double*));
	for(i=0;i<n_subpop;i++) 
		new_mig[i]=(double*)c_malloc(n_subpop*sizeof(double));
	k=0;
	mig_dim=count;
	m_pts=1;
	while(k<count) {
		scount=fscanf(m_in,"%d %d %lf %lf %d",
			m_index[0]+k,m_index[1]+k,
			m_end[0]+k,m_end[1]+k,m_surface_pts+k);
		if(scount!=5) {
			printf("\n\tProblem reading migration parameter file %s\n",filename);
			exit(1);
		}
		if(mig_matrix[m_index[0][k]][m_index[1][k]]==0.0) {
				printf("\n\tMigration matrix entry %d %d must be non-zero\n",
				m_index[0][k],m_index[1][k]);
				exit(1);
		}
		m_pts *= m_surface_pts[k];
		m_surface_pts[k]--;
		/* Set up for nextlex is 0<= index <= end */
		/* Debug only
		printf("%d %d %lf %lf %d\n",
			m_index[0][k],m_index[1][k],
			m_end[0][k],m_end[1][k],m_surface_pts[k]);
		*/
		k++;
	}
	fclose(m_in);
}

/* Run through in lexiographic order a string of length dim
    with elements lex[] such that 0 <= lex[i] <= range[i],
    i = 0,...,dim-1.
    Return the position which is incremented, or -1 when
    finished.
    Initially call with lex[] = {0,...,0}.
    Automatically resets when finished.
*/    

int nextlex(int *range,int *lex,int dim) {
    int last;   
    for(last = 0;last < dim & lex[last] == range[last];last++)
    lex[last] = 0;
    if(last == dim) return -1;
    lex[last]++;
    return last;
}

int new_migration_matrix(int flag) {
	int i,j,k,l=0;
	if(flag==0) {
		for(i=0;i<mig_dim;i++) lexi[i]=0;
		for(i=0;i<n_subpop;i++) for(j=0;j<n_subpop;j++) 
				new_mig[i][j]=mig_matrix[i][j];
		for(i=0;i<mig_dim;i++)
			new_mig[m_index[0][i]][m_index[1][i]]=m_end[0][i];
	}
	else {
		l=nextlex(m_surface_pts,lexi,mig_dim);
		for(i=0;i<mig_dim;i++)
			new_mig[m_index[0][i]][m_index[1][i]]
					= m_end[0][i] 
							+ lexi[i]*(m_end[1][i]-m_end[0][i])/m_surface_pts[i];
	}
	return l;
}




/*************  Example useage for nextlex
int lex[5]={0,0,0,0,0};
int range[5]={2,3,4,5,6};

main() {
	int dim=5,i,j;
	while(nextlex(range,lex,5)!=-1) {
		for(i=0;i<5;i++) printf("%d ",lex[i]);
		printf("\n");
	}
}

main(int argc, char **argv) {
	read_mutation(argv[1]);
}
***************/

void printtree(TREE *t, int flag) {
	int i,j;
	for(i=0;i<t->types;i++) {
		if(flag==2||debug==1) {
			if(n_subpop>1)
				printf("\n\t%3d; %3d %3d : ",i,t->color[i],t->multiplicity[i]);
			else
				printf("\n\t%3d; %3d : ",i,t->multiplicity[i]);
		}
		else {
			if(n_subpop>1)
				printf("\n\t%3d %3d : ",t->color[i],t->multiplicity[i]);
			else
				printf("\n\t%3d : ",t->multiplicity[i]);
		}
		j=0;
		while(t->sequence[i][j] != 0) {
			printf("%3d",t->sequence[i][j]);
			j++;
		}
		printf("%3d",0);
	}
	if(flag==1) { /* Verbose */
		printf("\n\tn = %d, types = %d, singletons= %d\n",
			t->n,t->types,t->singletons);
		printf("\n\n\tSingleton sequences (with interior parents below)\n\t");
		for(i=0;i<t->singletons;i++) printf("%d ",t->singleton_list[i]);
		printf("\n\t");
		for(i=0;i<t->singletons;i++) {
			if(t->singleton_exterior[i] == -1)
				printf("* ");
				else printf("%d ",t->singleton_exterior[i]);
		}
		printf("\n");
	printf("\n\tSubpopulation counts");
	printf("\n\t"); for(i=0;i<n_subpop;i++) printf("%4d",i);
	printf("\n\t"); for(i=0;i<n_subpop;i++) printf("%4d",t->subpop[i]);
	printf("\n");
	}
}

void copy_tree(TREE *t,TREE *s) {
	int i,flag=0;
	t->n=s->n;
	t->types=s->types;
	for(i=0;i<t->types;i++) 
		flag=(memcpy(t->sequence[i],s->sequence[i],maxlength*sizeof(int))==NULL);
	flag |=(memcpy(t->multiplicity,s->multiplicity,t->types*sizeof(int))==NULL);
	flag |=(memcpy(t->subpop,s->subpop,n_subpop*sizeof(int))==NULL);
	flag |=(memcpy(t->color,s->color,t->types*sizeof(int))==NULL);
	t->singletons=s->singletons;
	flag |=(memcpy(t->singleton_list,
								s->singleton_list,t->singletons*sizeof(int))==NULL);
	flag |=(memcpy(t->singleton_exterior,
								s->singleton_exterior,t->singletons*sizeof(int))==NULL);
	if(flag==1) {
		printf("\n\tCan't copy tree\n");
		printtree(s,1);
		exit(1);
	}
}

int get_sites(TREE* t) {
	int i,j,m=0;
	for(i=0;i<t->types;i++) for(j=0;t->sequence[i][j]!=0;j++)
		if(t->sequence[i][j]>m) m=t->sequence[i][j];
	return m;
}

int check_singleton(TREE *t,int z) {
	int x,j,k,flag;
	x=t->sequence[z][0];
	if(x == 0 || t->multiplicity[z] > 1) return 0;
	for(j=0;j<t->types;j++) {
		flag=1;
		for(k=0;flag==1;k++) {
			if(j != z && t->sequence[j][k]==x)
				return 0;
			flag = t->sequence[j][k] != 0;
		}	
	}
	return 1;
}

void get_connections(TREE *t) {
	int i,j,k;
	for(k=0;k<t->singletons;k++) {
		t->singleton_exterior[k]=(-1);
		j=t->singleton_list[k];
		for(i=0;i<t->types;i++)
			if(i !=j && t->sequence[i][0]==t->sequence[j][1]
				&&t->color[i]==t->color[j])
				t->singleton_exterior[k]=i;
	}
}

void get_singleton_list(TREE *t) {
	int i,j,k=0,l,flag;
	t->singletons=0;
	for(i=0;i<t->types;i++)  {
		if(check_singleton(t,i)==1) {
			t->singletons++;
			t->singleton_list[k++]=i;
		}
	}
	get_connections(t);
}

int delete_new[2];
int *delete_singleton(TREE* t,int j) {
	int i,k,flag=0;
	if(t->multiplicity[j]!=1) {
		printf("\n\tProblem in delete_singleton(), sequence %d\n",j);
		exit(1);
	}
	for(i=0;i<t->types;i++)
		if(i!=j && t->sequence[i][0] == t->sequence[j][1]
		 	&& t->color[i]==t->color[j]) {
			delete_new[0]=i;delete_new[1]=t->multiplicity[i];
			t->multiplicity[i]++;
			t->types--;
			for(k=j;k<t->types;k++) {
				memcpy(t->sequence[k],t->sequence[k+1],maxlength*sizeof(int));
				t->multiplicity[k]=t->multiplicity[k+1];
				t->color[k]=t->color[k+1];
			}
			/* t->subpop stays the same */
			flag=1;
			break;
	}
	if(flag==0)  { /* Delete first co-ordinate */
		k=0;
		while(t->sequence[j][k] != 0) {
			t->sequence[j][k]=t->sequence[j][k+1];
			k++;
		}
	}
	get_singleton_list(t);
	if(flag==0) return NULL; else return delete_new;
}

int c_printf(char *msg) {
	printf(msg);return 0;
}

char *c_err[] = {
"\n\tTree pointer is NULL\n",
"\n\tt->n is less than 1\n",
"\n\tt->types is <1 or > t->n\n",
"\n\tt->sequence is NULL\n",
"\n\tt->sequence[i] is NULL\n",
"\n\tt->multiplicity is NULL\n",
"\n\tSum of multiplicities is not t->n\n",
"\n\tt->singleton_list is NULL\n",
"\n\tt->singleton_exterior is NULL\n",
"\n\tt->singleton_list is incorrect\n",
"\n\tt->singleton_exterior is incorrect\n",
"\n\tt->subpopulation count is incorrect\n",
"\n\tt->multiplicity[i] is zero\n",
"\n\tt->types > maxtype\n",
};
int check_sub[20];
int check_tree(TREE *t) { 
	int sum=0,i,flag=1,j,k;
	if(t==NULL) flag=c_printf(c_err[0]);
	if(t->n < 1) flag=c_printf(c_err[1]);
	if(t->types < 1 || t->types > t->n) flag=c_printf(c_err[2]);
	if(t->types > maxtype) flag=c_printf(c_err[13]);
	if(t->sequence==NULL) flag=c_printf(c_err[3]);
	for(i=0;i<t->types;i++) if(t->sequence[i]==NULL) flag=c_printf(c_err[4]);
	if(t->multiplicity==NULL) flag=c_printf(c_err[5]);
	for(i=0;i<t->types;i++) {
		if(t->multiplicity[i]==0) flag=c_printf(c_err[12]);
		sum += t->multiplicity[i];
	}
	if(sum != t->n) flag=c_printf(c_err[5]);
	if(t->singleton_list==NULL) flag=c_printf(c_err[6]);
	if(t->singleton_exterior==NULL) flag=c_printf(c_err[7]);
	if(flag==0) printf("\n\tmaxlength=%d\n",maxlength);
	for(i=0;i<t->singletons;i++) 
		if(check_singleton(t,t->singleton_list[i])==0) 
			flag=c_printf(c_err[8]);
	for(k=0;k<t->singletons;k++) {
		j=t->singleton_list[k];
		for(i=0;i<t->types;i++)
			if(i !=j && t->sequence[i][0]==t->sequence[j][1]
			&& t->color[i]==t->color[j]
			)
				if(t->singleton_exterior[k] != i)
					flag=c_printf(c_err[9]);
	}
	for(i=0;i<n_subpop;i++) check_sub[i] = 0;
	for(i=0;i<t->types;i++) check_sub[t->color[i]] += t->multiplicity[i];
	for(i=0;i<n_subpop;i++) 
		if(check_sub[i] != t->subpop[i]) flag=c_printf(c_err[11]);
	return flag;
}

void get_memory(void) {
	int i,j;
	move_prob=(double*)c_malloc((2+maxtype*(2+n_subpop))*sizeof(double));
	subsize=(double*)c_malloc(n_subpop*sizeof(double));
	keep_subsize=(double*)c_malloc(n_subpop*sizeof(double));
	mig_matrix=(double**)c_malloc(n_subpop*sizeof(double*));
	m_age=(double*)c_malloc((n_sites+1)*sizeof(double));
	t_rate=(double*)c_malloc((n_sites+1)*sizeof(double));
	sage=(double*)c_malloc((n_sites+1)*sizeof(double));
	ssage=(double*)c_malloc((n_sites+1)*sizeof(double));
	anc_trace=(int**)c_malloc(2*sizeof(int*));
	anc_trace[0]=(int*)c_malloc(n_sample*sizeof(int));
	anc_trace[1]=(int*)c_malloc(n_sample*sizeof(int));
	check_anc=(int*)c_malloc(n_sample*sizeof(int));
	initial_color=(int*)c_malloc(n_sample*sizeof(int));
	mrca_test=(int*)c_malloc(n_sample*sizeof(int));
	csite=(double**)c_malloc((n_sites+1)*sizeof(double*));
	m_csite=(int*)c_malloc((n_sites+1)*sizeof(int));
	mrca_c=(double**)c_malloc(n_subpop*sizeof(double*));
	tmrca_c=(double*)c_malloc(n_subpop*sizeof(double));
	ss_tmrca_c=(double*)c_malloc(n_subpop*sizeof(double));
	initial_subpop=(int*)c_malloc(n_subpop*sizeof(int));
	m_tmrca_c=(double*)c_malloc(n_subpop*sizeof(double));
	m_mrca_c=(int*)c_malloc(n_subpop*sizeof(int));
		move_migration=(int**)malloc(n_subpop*sizeof(int*));
	for(i=0;i<n_subpop;i++) 
		move_migration[i]=(int*)malloc(n_subpop*sizeof(int));
	for(i=0;i<n_subpop;i++) {
		mig_matrix[i]=(double*)c_malloc(n_subpop*sizeof(double));
		mrca_c[i]=(double*)c_malloc(n_subpop*sizeof(double));
	}
	for(i=0;i<n_subpop;i++) {
		tmrca_c[i]=0.0;
		ss_tmrca_c[i]=0.0;
		for(j=0;j<n_subpop;j++) mrca_c[i][j]=0.0;
	}
	for(i=0;i<=n_sites;i++) 
		csite[i]=(double*)c_malloc(n_subpop*sizeof(double));
	for(i=0;i<n_sample;i++) {
		anc_trace[0][i]=anc_trace[1][i]=(-1);
	}
	for(i=0;i<=n_sites;i++) {
		t_rate[i]=1.0;
		m_age[i]=0.0;
		for(j=0;j<n_subpop;j++)  csite[i][j]=0.0;
	}
	if(n_subpop==1) subsize[0]=1.0;
	if(n_subpop >1) for(i=0;i<n_subpop;i++) {
		subsize[i] = 1.0/(double)n_subpop;
		mig_matrix[i][i]=0.0;
		for(j=0;j<n_subpop;j++) {
			if(i!=j) mig_matrix[i][j]=1.0/(n_subpop-1.0);
		}
	}
}

double ***get3array(int *d) {
	int i[3];
	double ***p=NULL;
		p=(double***)c_malloc(d[0]*sizeof(double**));
	for(i[0]=0;i[0]<d[0];i[0]++) {
		p[i[0]]=(double**)c_malloc(d[1]*sizeof(double*));
			for(i[1]=0;i[1]<d[1];i[1]++) {
				p[i[0]][i[1]]=(double*)c_malloc(d[2]*sizeof(double));
					for(i[2]=0;i[2]<d[2];i[2]++) p[i[0]][i[1]][i[2]]=0.0;
			}
	}
	return p;
}

void get_surface_mem(int *d) {
	int i,j,k;
	surface_h=get3array(d);
	surface_fnl=get3array(d);
	ssurface_fnl=get3array(d);
	sssurface_fnl=get3array(d);
	if(tmrca_surface==1) {
		tm_surface=get3array(d);
		ttm_surface=get3array(d);
	}
}

/* XXX Pairwise difference calculation */	

int common_site(int i, int j, TREE *t) {
	int k,l;
	for(k=0;t->sequence[i][k];k++) for(l=0;t->sequence[j][l];l++) 
		if(t->sequence[i][k]==t->sequence[j][l]) return t->sequence[i][k];
	return 0;
}

int lineage_sites(int i, int j, TREE *t) {
	int count=0,k,flag;
	flag=common_site(i,j,t);
	for(k=0;t->sequence[i][k]!=flag;k++) count++;
	for(k=0;t->sequence[j][k]!=flag;k++) count++;
	return count;
}

double harmonic(int n) {
	double x=0.0;
	int i;
	if(n==0) return 0.0;
	for(i=1;i<=n;i++) x += 1.0/(double)i;
	return x;
}

double invsq_series(int n) {
	double x=0.0,y;
	int i;
	if(n==0) return 0.0;
	for(i=1;i<=n;i++) {
		y = (double)i; x += 1.0/(y*y);
	}
	return x;
}

void pairwise(TREE *t) {
	int i,j,k,mn,mx,s,in,flag=0;
	double var=0.0,kp=0.0,n,m,varm,d,a1,a2,b1,b2,e1,e2;
	int **site_color;
	double **pair_diff; /* Average pairwise differences */
	site_color=(int**)c_malloc(n_subpop*sizeof(int*));
	pair_diff=(double**)c_malloc(n_subpop*sizeof(double*));
	for(i=0;i<n_subpop;i++) flag |= t->subpop[i]>1;
	for(i=0;i<n_subpop;i++) {
		site_color[i]=(int*)c_malloc((n_sites+1)*sizeof(int));
		pair_diff[i]=(double*)c_malloc(n_subpop*sizeof(double));
	  for(j=0;j<n_subpop;j++)	pair_diff[i][j]=0.0;
		for(j=0;j<n_sites+1;j++) site_color[i][j]=0;
	}
	for(i=0;i<t->types;i++) for(j=0;t->sequence[i][j];j++)
		site_color[t->color[i]][t->sequence[i][j]] += t->multiplicity[i];
	/* Count segregating sites */
	for(i=0;i<n_subpop;i++) for(j=1;j<n_sites+1;j++)  {
		if(site_color[i][j]>0&&site_color[i][j]!=t->subpop[i])
			site_color[i][0]++;
	}
	for(i=0;i<t->types;i++) for(j=i+1;j<t->types;j++) {
			mn=t->color[i]<t->color[j] ? t->color[i] : t->color[j];
			mx=t->color[i]>t->color[j] ? t->color[i] : t->color[j];
			pair_diff[mn][mx] += 
				(double)t->multiplicity[i]*(double)t->multiplicity[j]
													*(double)lineage_sites(i,j,t);
	}
	if(flag==1) {
		printf("\n\tPairwise difference averages, SD");
		printf("\n\tSample size, Segregating sites, Theta estimate M, SD, D\n");
	} else printf("\n\tPairwise difference averages");
	for(i=0;i<n_subpop;i++) if(t->subpop[i]>1) {
		n=t->subpop[i];
		in=t->subpop[i]-1;
		if(n>1) {
			pair_diff[i][i] /= n*(n-1.0)*0.5;
			kp=pair_diff[i][i];
			var=(3.0*n*(n+1)*kp+2.0*(n*n+n+3.0)*kp*kp)/
					(11.0*n*n-7.0*n+6.0);
		}
		else {
			kp=0.0;var=0.0;
		}
		s=site_color[i][0];
		a1=harmonic(in);
		a2=invsq_series(in);
		if(s>0&&a1>0) {
			m=s/a1;
			varm=(a1*a1*s+a2*s*s)/(a1*a1*(a1*a1+a2));
		}
		else { m=0.0; varm=0.0;}
		if(s>1&&a1>0) {
			b1=(n+1.0)/(3.0*(n-1.0));
			b2=2.0*(n*n+n+3)/(9.0*n*(n-1.0));
			e1=(b1-1.0/a1)/a1;
			e2=(b2-(n+2.0)/(a1*n) + a2/(a1*a1))/(a1*a1+a2);
			if(e1*s+e2*s*(s-1.0)>0.0) d=(kp-s/a1)/sqrt(e1*s+e2*s*(s-1.0)); 
				else d=0.0;
		} else d=0.0;
	  printf("\tPopulation    %d:  k=%.2lf, se=%.2lf\n",i,kp,sqrt(var));
		if( (d==0.0) && ((a1*kp-s)!=0.0) ) { /* Problems with var */
		printf("\tn=%4d            s=%2d, m=%.2lf, se=%.2lf\n",
											t->subpop[i],s,m,sqrt(varm));
		}
		else
		 printf("\tn=%4d            s=%2d, m=%.2lf, se=%.2lf, d=%12.4e\n",
											t->subpop[i],s,m,sqrt(varm),d);
	}
	if(flag==0) printf("\n");
	printf("\tBetween populations\n");
	for(i=0;i<n_subpop;i++) for(j=i;j<n_subpop;j++) if(i!=j) {
		pair_diff[i][j] /= t->subpop[i]*t->subpop[j];
		printf("\tPopulations %d %d:  k=%.2lf\n",i,j,pair_diff[i][j]);
	}
}

double mig(int i, int j, double tim, double m) {
	int k,l;
	if(n_subpop==1||i==j) return 0.0;
	else {
	if(use_new_mig==0) return m*mig_matrix[i][j];
	if(use_new_mig==1) return m*new_mig[i][j];
	}
}

int initial_gf(double theta, double m) {
	double **gf, **mg,*n;
	int g,s,i,j,k,si,sj;
	gf0=(double**)c_malloc((n_sites+1)*sizeof(double*));
	for(i=0;i<=n_sites;i++)
	gf0[i]=(double*)c_malloc((n_subpop*(n_subpop+1)/2)*sizeof(double));
	g=n_subpop;
	s=n_sites;
	n=subsize;
	gf=gf0;
	mg=(double**)c_malloc(g*sizeof(double*));
	for(i=0;i<g;i++) 
		mg[i]=(double*)c_malloc(g*sizeof(double));
	for(i=0;i<g;i++) for(j=0;j<g;j++) 
		if(i==j) mg[i][j]=0.0; else mg[i][j]=mig(i,j,1.0,m);
	get_gf(s,gf,g,n,theta,mg);
	if(surface==0) return 0;
	s_gf=(double****)c_malloc(theta_pts*sizeof(double***));
	for(i=0;i<theta_pts;i++) {
		s_gf[i]=(double***)c_malloc(m_pts*sizeof(double**));
		for(j=0;j<m_pts;j++) {
			s_gf[i][j]=(double**)c_malloc((n_sites+1)*sizeof(double*));
				for(k=0;k<=n_sites;k++)
					s_gf[i][j][k]=(double*)c_malloc((n_subpop*(n_subpop+1)/2)
																								*sizeof(double));
		} /* end j */
	} /* end i */
	if(matrix_surface==0) for(si=0;si<theta_pts;si++) for(sj=0;sj<m_pts;sj++)  {
		gf=s_gf[si][sj];
		for(i=0;i<g;i++) for(j=0;j<g;j++) if(i!=j)
			 mg[i][j]=mig(i,j,1.0,surface_m[sj]);
		get_gf(s,gf,g,n,surface_theta[si],mg);
	}
	if(matrix_surface==1) for(si=0;si<theta_pts;si++) for(sj=0;sj<m_pts;sj++)  {
		gf=s_gf[si][sj];
		new_migration_matrix(sj);
		for(i=0;i<g;i++) for(j=0;j<g;j++) mg[i][j]=new_mig[i][j];
		get_gf(s,gf,g,n,surface_theta[si],mg);
	}
	return 1;
}
		
void write_tree(TREE *t, char *total_name) {
	int i,j,k;
	int *s;
	FILE *combine;
	s=(int*)c_malloc((n_sites+1)*sizeof(int));
	for(i=0;i<=n_sites;i++) s[i]=0;
	for(i=0;i<t->types;i++)  s[t->sequence[i][0]] += t->multiplicity[i];
	combine=fopen(total_name,"w");
	if(combine==NULL) {
		printf("\n\tCan't open a file for the total tree\n");
		exit(1);
	}
	for(i=0;i<=n_sites;i++) if(s[i]>0) {
		for(j=0;j<t->types;j++) 
			if(t->sequence[j][0]==i) {
				fprintf(combine,"%3d :",s[i]);
				for(k=0;t->sequence[j][k];k++) 
					fprintf(combine," %d",t->sequence[j][k]);
			fprintf(combine," 0\n");
			break;
			}
	}
	fclose(combine);
	free(s);
}
			

int get_tree(char *filename,TREE *t) {
	FILE *in;
	int i,j,k,l,x,y,z=(-1),length=0,one_pop=0,flag;
	/* maxlength is global, so it can be used in copy_tree() */
	char separator[20];
	if((in=fopen(filename,"r"))==NULL) {
		printf("\n\tCannot open %s\n",filename);
		exit(1);
	}
	/* Check if colors given, if not assume a single population */
	fscanf(in,"%d",&x);
	fscanf(in,"%s",separator);
	if(separator[0]==':') one_pop=1; 
	rewind(in);
	/* Check done */
	t->types=0;
	t->n=0;
	while(feof(in) == 0) {
		if(one_pop==0) {
			if(fscanf(in,"%d %d",&x,&y)==2) {
				if(y>0) { /* Don't count multiplicity 0 sequences */
					t->types++;
					if(x>z) z=x;
				}
			}
		}
		else {
			if(fscanf(in,"%d",&y)==1) {
				if(y>0) {
					t->types++;
					z=x=0;
				}
			}
		}
		if(feof(in)!=0) break;
		fscanf(in,"%s",separator);
		x=(-1);
		while(x != 0 && feof(in) == 0) {
			fscanf(in,"%d",&x);
			length++;
		}
		if(length>maxlength) maxlength=length;
			length =0;
	}
	if(t->types==0) {
		printf("\n\tProblem with input file  %s\n",filename);
		exit(1);
	}
	if(t->types>maxtype) {
		printf("\n\tToo many types (> %d), set a  larger maxtype\n",
											maxtype);
		exit(1);
	}
	if(set_n_subpop<z+1) n_subpop=z+1;
	else n_subpop=set_n_subpop;
	t->sequence=(int**)c_malloc(maxtype*sizeof(int*));
	for(i=0;i<maxtype;i++)  
		t->sequence[i]=(int*)c_malloc(maxlength*sizeof(int));
	t->singleton_list=(int*)c_malloc(maxtype*sizeof(int));
	t->singleton_exterior=(int*)c_malloc(maxtype*sizeof(int));
	t->color=(int*)c_malloc(maxtype*sizeof(int));
	t->multiplicity=(int*)c_malloc(maxtype*sizeof(int));
	t->subpop=(int*)c_malloc(n_subpop*sizeof(int));
	for(i=0;i<n_subpop;i++) t->subpop[i]=0;
	rewind(in);
	z=(-1);
	i=0;
	while(i < t->types) {
		flag=1; /* Set to 0 id multiplicity read in is 0 */
		if(one_pop==0) {
			if(fscanf(in,"%d %d",&y,&x)!=2) {
				printf("\n\tProblem with reading data, sequence %d\n",i);
				exit(1);
			}
		}
		else {
			if(fscanf(in,"%d",&x)!=1) {
				printf("\n\tProblem with reading data, sequence %d\n",i);
				exit(1);
			}
			y=0;
		}
		if(x==0) flag=0;
		if(flag==1) {
			t->multiplicity[i]=x;
			t->color[i]=y;
			t->subpop[t->color[i]] += t->multiplicity[i];
			t->n += t->multiplicity[i];
		}
		fscanf(in,"%s",separator);
		if(separator[0] != ':') {
			printf("\n\tProblem with data format\n");
			exit(1);
		}
		j=0;
		x=(-1);
		while(x != 0 && feof(in) == 0) {
			fscanf(in,"%d",&x);
			if(flag==1) {
				t->sequence[i][j++]=x;
				if(x>z) z=x;
			}
		}
		if(flag==1) for(l=0;l<i;l++) 
			if(t->color[i]==t->color[l]&&t->sequence[i][0]==t->sequence[l][0]) {
				printf("\n\tError in tree file: %d,%d identical types?\n\n",i+1,l+1);
				exit(1);
		}
		if(flag==1) i++;
	}
	n_sample=t->n;
	if(n_sample < 2) {
		printf("\n\tSample size should be > 1\n");
		printf("\tproblem with input file  %s ?\n",filename);
		exit(1);
	}
	n_sites=z; /* Don't count 0 in site list */
	rewind(in);
	fclose(in);
  get_singleton_list(t);
	if(check_tree(t)==0) {
		printf("\n\tProblem with tree, inconsistent data\n");
		exit(1);
	}
	write_tree(t,tree_gt);
	if(debug==1) printf("\n\tTree read in and checked\n");
	return 1;
}

void free_tree(TREE* t) {
	int i;
	free(t->multiplicity);
	free(t->singleton_list);
	free(t->singleton_exterior);
	free(t->subpop);
	free(t->color);
	for(i=0;i<t->types;i++) free(t->sequence[i]);
	free(t->sequence);
	maxlength=0;
}

/* ancestor trace functions have types and multiplicities prior to changes*/

void print_trace(TREE *t) {
	int i,j,k,l;
	printf("\n\t          ");
	for(i=0;i<n_sample;i++) printf("%2d",initial_color[i]);
	printf("\n\t          ");
	for(i=0;i<n_sample;i++) printf("%2d",i%100);
	for(i=0,k=0;i<t->types;i++) {
		for(j=0;j<t->multiplicity[i];j++) {
			if(j==0) printf("\n\t%2d%4d",t->color[i],t->multiplicity[i]);
				else printf("\n\t      ");
			printf("%4d",k);
			for(l=0;l<n_sample;l++) 
				if(anc_trace[0][l]==i&&anc_trace[1][l]==j) printf(" +");
				else printf(" .");
			k++;
		}
	}
	printf("\n");
}

int anc_index(TREE *t, int i, int j) {
	int r,s,k=0;
	for(r=0,k=0;r<t->types;r++) for(s=0;s<t->multiplicity[r];s++) {
		if(r==i&&s==j) return k;
		k++;
	}
	printf("\n\tIndex outside range\n");
	exit(1);
	return (-1);
}

void check_trace(TREE *t) {
	int i,j,k;
	for(k=0;k<n_sample;k++) 
		if(
				anc_trace[0][k]<0
			||anc_trace[1][k]<0
			||anc_trace[0][k]>=t->types
		 	||anc_trace[1][k]>=t->multiplicity[anc_trace[0][k]]
		) {
			print_trace(t);
			printf("\n\tProblem with ancestor of %d\n",k);
			exit(1);
		}
		for(k=0;k<n_sample;k++) check_anc[k]=0;
		for(k=0;k<n_sample;k++) 
			check_anc[anc_index(t,anc_trace[0][k],anc_trace[1][k])]++;
		for(k=0;k<t->types;k++)  if(check_anc[k]==0) {
			print_trace(t);
			printf("\n\tAncestor %d has no decendents\n",k);
			exit(1);
		}
}

void initial_trace(TREE *t) {
	int i,j,k=0;
	for(i=0;i<t->types;i++) for(j=0;j<t->multiplicity[i];j++) {
		initial_color[k]=t->color[i];
		anc_trace[0][k]=i;anc_trace[1][k++]=j;
	}
}

int query_test[2];

int *anc_query(int c) {
 /* Return NULL if sample genes colour c have >1 ancestor in the
 ancestor tree , or return sequence number, number in multiplicity
 ordering query_test [0], [1].
 */
	int i,j,k,l,count=0;
	query_test[0]=(-1);
	for(k=0;k<n_sample;k++) {
 		if(initial_color[k]==c) {
			count++;
			if(query_test[0]==(-1)) {
				query_test[0]=anc_trace[0][k];
				query_test[1]=anc_trace[1][k];
			}
			else 
				if(query_test[0]!=anc_trace[0][k]||query_test[1]!=anc_trace[1][k])
		return NULL;
		}
		/* if(count==initial_subpop[c]) break; */
	}
	return query_test;
}

int *mrca_query(int c) {
	int *a;
	if(mrca_test[c]==1||(a=anc_query(c))==NULL) return NULL;
	mrca_test[c]=1;
	return a;
}
	

int coalesce_trace(int i,int mi) {
/* coalesce two from the i th type */
		int j,k,l,r,mini=(-1),maxi=(-1),count=0,flag=0,tmp;
		j=randint(mi); /* Random pair */
		k=randint((mi-1));
		if(k>=j) k++;
		else { l=j;j=k;k=l;} /* Now j is the minimum */
		if(out_level>=5) {
			for(r=0;r<n_sample;r++)  {
				if(anc_trace[0][r]==i) {
					if(anc_trace[1][r]==k) anc_trace[1][r]=j;
					if(mini==(-1)) mini=r; maxi=r;
				}
			}
	  	  for(r=mini;r<=maxi;r++) { 
				for(l=k+1;l<mi;l++)
				if(anc_trace[0][r]==i&&anc_trace[1][r]==l) anc_trace[1][r]--;
			}
		}
		return 1;
}

void delete_trace(int i) { /* Remove group i and close up indexing */
	int r;
	for(r=0;r<n_sample;r++) 
		if(anc_trace[0][r] > i) anc_trace[0][r]--;
}

void shift_trace(int i,int mi, int j, int mj) {
 /* Shift one line from i,mi group to j mj group */
  int k,r,count=0;
	if(out_level>=5) { /* Otherwise do nothing */
		k=randint(mi);
		for(r=0;r<n_sample;r++) {
				if(anc_trace[0][r]==i&&anc_trace[1][r]==k) {
					anc_trace[0][r]=j;
					anc_trace[1][r]=mj; /* Put it at the end of the group */
				}
				if(anc_trace[0][r]==i&&anc_trace[1][r]>k) anc_trace[1][r]--;
		}
		if(mi==1) delete_trace(i);
	}
}

void coalesce(TREE *t, int j) { /* Coalesce sequence j */
	if(t->multiplicity[j]<2) {
		printf("\n\tProblem in coalesce, sequence %d\n",j);
		printtree(t,1);
		exit(1);
	}
	t->n--;
	t->multiplicity[j]--;
	t->subpop[t->color[j]]--;
	if(check_singleton(t,j)==1) {
		get_singleton_list(t);
	}
}

/* Migrate 1 individual the jth sequence to color c,
	 return new sequence number */

int migrate(TREE *t, int j, int c) {
	int i,k,l,a;
	l=t->multiplicity[j]--;
	t->subpop[t->color[j]]--;
	t->subpop[c]++;
	for(i=0;i<t->types;i++)
		if(i!=j&&t->sequence[i][0]==t->sequence[j][0]&&t->color[i]==c) {
			shift_trace(j,t->multiplicity[j]+1,i,t->multiplicity[i]);
			t->multiplicity[i]++;
			if(l==1) {
				t->types--;
				for(k=j;k<t->types;k++) {
					memcpy(t->sequence[k],t->sequence[k+1],maxlength*sizeof(int));
					t->multiplicity[k]=t->multiplicity[k+1];
					t->color[k]=t->color[k+1];
				}
			get_singleton_list(t);
			if(j<i) return i-1; else return i;
			}
			get_singleton_list(t);
			return i;
		}

		if(t->multiplicity[j]==0) { /* Leave sequence in same place */
			t->multiplicity[j]=1;
			t->color[j]=c;
			get_singleton_list(t);
			return j;
		}

		if(t->multiplicity[j]>0) { /* Make a new type */
			i=t->types;
			t->types++;
			if(t->types>maxtype) {
				printf("\n\tToo many types (> %d), set a larger maxtype\n",
											maxtype);
				exit(1);
			}
			shift_trace(j,t->multiplicity[j]+1,i,0);
			t->multiplicity[i]=1;
			t->color[i]=c;
			k=0;
			while(t->sequence[j][k]!=0) {
				t->sequence[i][k]=t->sequence[j][k];
				k++;
			}
			t->sequence[i][k]=0;
			get_singleton_list(t);
			return i;
		}
		return -1;
}

double lambda(int c, double tim) {
	if(growth==1) {
		if(growth_rate[c]==0.0) return 1.0;
		else return exp(tim*growth_rate[c]);
	}
	return 1.0;
}

/* Integrated  (tim0,tim1) */
double ilambda(int c, double tim0, double tim1) {
	if(growth==1) {
		if(growth_rate[c]==0.0) return tim1-tim0;
		else if(tim0==0.0)
			return (exp(tim1*growth_rate[c])-1.0)/growth_rate[c];
		else return (exp(tim1*growth_rate[c])-exp(tim0*growth_rate[c]))
			/growth_rate[c];
	}
	return 1.0;
}

double rate(TREE *t, double tim, double theta, double m) {
	int i,j;
	double sum=0.0,margin;
	/* Mutation */
	sum = 0.5*t->n*theta;
	/* Coalescence and migration */
	if(n_subpop>1) for(i=0;i<n_subpop;i++) {
		if(t->subpop[i]>1)
			sum += 0.5*t->subpop[i]*(t->subpop[i]-1.0)
							*lambda(i,tim)/subsize[i];
		if(t->subpop[i]>0) {
			margin=0.0;
			for(j=0;j<n_subpop;j++) if(j!=i)  margin += mig(i,j,tim,m);
				sum += t->subpop[i]*margin;
		}
	}
	else sum += 0.5*t->n*(t->n-1.0)*lambda(0,tim)/subsize[0];
	if(sum<=0) {
		printf("\n\tBug with rate = %.4e",sum);
		printtree(t,1);
		exit(1);
	}
	return sum;
}

/* This function returning the integrated rate eventually needs rewriting 
for more clarity and generality */

double irate(TREE *t, double timb, double tima, double theta, double m) {
	int i,j;
	double sum=0.0,tum=0.0,margin;
	if(timb==tima) return 0.0;
	if(growth==0) return rate(t,tima,theta,m)*(timb-tima);
	sum = 0.5*t->n*theta;
	for(i=0;i<n_subpop;i++) {
		if(t->subpop[i]>1) {
			if(growth_rate[i]>0.0) tum += 0.5*t->subpop[i]*(t->subpop[i]-1.0)
		 *(lambda(i,timb)-lambda(i,tima))/(subsize[i]*growth_rate[i]);
		 else
			 tum += (timb-tima)*0.5*t->subpop[i]*(t->subpop[i]-1.0)/subsize[i];
		}
	/* Only works for constant migration rates,
	   modify so its ok for a mig() function
		 with exponential rates or more general
	*/
		if(t->subpop[i]>0) {
			margin=0.0;
			for(j=0;j<n_subpop;j++) if(j!=i)  margin += mig(i,j,tima,m);
			sum += t->subpop[i]*margin;
		}
	}
	if(sum<=0) {
		printf("\n\tBug with rate = %.4e",sum);
		printtree(t,1);
		exit(1);
	}
	return tum + (timb-tima)*sum;
}

double migration_time( TREE* t,int i, int j, double m, double tim) {
	double denom;
	denom=t->subpop[i]*mig(i,j,tim,m);
	if(denom == 0.0)  return -1.0; /* Used in next_time */
	return tim -log(RAND01)/denom;
/* Assumes a constant migration rate, modify otherwise */
}

double c_time( TREE* t,int c, double tim) {
	if(t->subpop[c] <= 1) {
		printf("\n\tError in c_time()\n");
		exit(1);
	}
	if(growth==0)
		return tim -2.0*subsize[c]*log(RAND01)/(t->subpop[c]*(t->subpop[c]-1.0));
	if(growth_rate[c]==0.0) 
		return tim -2.0*subsize[c]*log(RAND01)/(t->subpop[c]*(t->subpop[c]-1.0));
	else if(growth_rate[c]>0.0) 
		return log(-2.0*subsize[c]*growth_rate[c]*log(RAND01)/(t->subpop[c]*(t->subpop[c]-1.0))
		+ exp(tim*growth_rate[c]))/growth_rate[c];
	return -1.0; /* Add general schemes here */
}

double mutation_time ( TREE* t, double theta, double tim) {
	return tim - 2.0*log(RAND01)/(t->n*theta);
}
	
double next_time(TREE *t, double theta, double m,double tim) {
	double min_time,work_time;
	int i,j;
		if(notime==1) return tim + 0.05; /* Dummy return only */
		min_time=mutation_time(t,theta,tim);
		for(i=0;i<n_subpop;i++)  {
			if(t->subpop[i]>1) {
				work_time=c_time(t,i,tim);
				if(work_time<min_time) min_time=work_time;
			}
			for(j=0;j<n_subpop;j++) if(j!=i) {
				work_time=migration_time(t,i,j,m,tim);
				if(work_time>0.0 && work_time <min_time) min_time=work_time;
			}
		}
		return min_time;
}

/* If the jth sequence is also on island color c, return the multiplicity
on island c
*/
int in_sample(TREE *t, int j, int c) {
	int i;
	for(i=0;i<t->types ;i++)
		if(t->sequence[i][0]==t->sequence[j][0]&&t->color[i]==c)
			return t->multiplicity[i];
	return 0;
}

/* Single step functional */
double  f(TREE* t,double theta, double m, double tim, double *pa, double *pb) {
	int i,j,single,c,site,count=1,l=1;
	double x,cumulate =0.0,p,*q,*r;
	move_prob[0]=0.0;
	/*
	if(n_subpop==1) for(i=0;i<t->types;i++)  if(t->multiplicity[i]>1) {
					move_prob[l] = 0.5*t->n*(t->multiplicity[i]-1.0)
							*lambda(0,tim)/subsize[0] +move_prob[l-1];
					l++;
	}
	*/
	for(i=0;i<t->types;i++)  if(t->multiplicity[i]>1) {
					move_prob[l] = 0.5*t->subpop[t->color[i]]*(t->multiplicity[i]-1.0)
							*lambda(t->color[i],tim)/subsize[t->color[i]]+move_prob[l-1];
					l++;
	}
	for(i=0;i<t->singletons;i++) {
		site=t->sequence[ t->singleton_list[i] ][0];
		single=1;
		if(t->singleton_exterior[i] != -1)
			single += t->multiplicity[t->singleton_exterior[i]];
		 if(margin_flag==0)
		 	move_prob[l]= 0.5*single*theta*t_rate[site]+move_prob[l-1];
		 else
		 	move_prob[l]= 0.5*single*t_rate[site]+move_prob[l-1];
		 l++;
	}
	if(n_subpop>1) for(j=0;j<t->types;j++)  for(c=0;c<n_subpop;c++) 
		if(t->color[j]!=c) {
			move_prob[l] = (in_sample(t,j,c)+1.0)*mig(t->color[j],c,tim,m)
						*t->subpop[t->color[j]]/(t->subpop[c]+1.0)+move_prob[l-1];
		l++;
		}
	x=*pa= move_prob[l-1];
	for(i=1;i<l;i++) move_prob[i] /= x;
	*pb=rate(t,tim,theta,m);
	return x/(*pb);
}

#define MUTATE 1
#define COALESCE 2
#define MIGRATE 4

double move(TREE* t,double theta,double m,double tim,double tim0,int *pflag) {
	int i,j,k,l,si,sj,sk,single,c,site,*new,seq,c0,c1,flag,id,jd,count=1;
	double x,p,a,b,s_a,s_b,comb=1.0;
	double *q;
	if(checktree==1) {
  /*  DEBUG only  */
		if(check_tree(t)==0) {
			printtree(t,1);
			printf("\n\t********************** ERROR **********************\n");
			exit(1);
		}
		if(out_level>=5) check_trace(t);
	}
   /*********************************************************/
	/* Special cases of t->n=2 */
	if(two_stop==1&&n_subpop==1&&t->n==2) {
	  i=0;j=0;
		if(t->multiplicity[0]<2) {
			while(t->sequence[0][i]!=0) i++;
			while(t->sequence[1][j]!=0) j++;
			comb=pow(0.5,(double)(i+j))*choose(i+j,j)*2.0;
		}
		if(growth==0) x=comb*two_single(i,j,theta);
		if(growth==1) x=comb*two_single_var(i,j,theta,tim0);
		t->n=1;t->types=1;t->multiplicity[0]=1;
	for(si=0;si<theta_pts;si++) for(sj=0;sj<m_pts;sj++) 
		for(sk=0;sk<g_pts;sk++) {
		if(growth==0) 
				surface_h[si][sj][sk]=comb*two_single(i,j,surface_theta[si]);
		if(growth==1) 
			surface_h[si][sj][sk]=
				comb*two_single_var_mag(i,j,surface_theta[si],tim0,surface_g[sk]);
		}
		return x;
	}
	/****/
	if(two_stop==1&&n_subpop>1&&t->n==2&&growth==0) {
	  i=0;j=0;l=0;c0=t->color[0];c1=t->color[0];
		if(t->multiplicity[0]<2) {
			while(t->sequence[0][i]!=0) i++;
			while(t->sequence[1][j]!=0) {
				k=0;flag=0;
				while(t->sequence[0][k]!=0) {
					if(t->sequence[0][k]==t->sequence[1][j]) {
						flag=1;
						break;
					}
					k++;
				}
			if(flag==1) l++;
			j++;
			}
			c1=t->color[1];
		}
		/* i on sequence 0, j on sequence 1 and l in common,
		 distinct are id and jd */
		id=i-l;jd=j-l;
		if(id+jd==0) x=pgf(0,c0,c1);
		else {
			comb=pow(0.5,(double)(id+jd))*choose(id+jd,id);
			if(c0==c1) comb *= 2.0;
			x=pgf(id+jd,c0,c1)*comb;
		}
		t->n=1;t->types=1;t->multiplicity[0]=1;
	for(si=0;si<theta_pts;si++) for(sj=0;sj<m_pts;sj++) 
		for(sk=0;sk<g_pts;sk++) {
		if(id+jd==0) surface_h[si][sj][sk]= spgf(0,c0,c1,si,sj);
		else surface_h[si][sj][sk]= comb*spgf(id+jd,c0,c1,si,sj);
	}
	return x;
}
	/* end case of t->n = 2 */
   /*********************************************************/
	x=f(t,theta,m,tim,&a,&b);
	if(debug==1) 
		printf("\n\tf = %.4e, a = %.4e, b = %.4e, time %.4e",x,a,b,tim);
	for(sj=0;sj<m_pts;sj++) {
		if(matrix_surface==1) new_migration_matrix(sj);
		for(si=0;si<theta_pts;si++) for(sk=0;sk<g_pts;sk++) {
			/*
			f(t,surface_theta[si],surface_m[sj],surface_g[sk]*tim,&s_a,&s_b);
			in early versions.
			*/
			if(matrix_surface==0)
				s_b=rate(t,surface_g[sk]*tim,surface_theta[si],surface_m[sj]);
			if(matrix_surface==1) {
					use_new_mig=1;
					s_b=rate(t,surface_g[sk]*tim,surface_theta[si],1.0);
					use_new_mig=0;
			}
			/* A modified rate function is needed for a different surface
			 * migration matrix. use_new_mig is the flag for this.
			 */
			surface_h[si][sj][sk]=a/s_b;
		}
	}
	p=RAND01;
	*pflag=0;
	count=1;
	q=move_prob+1;
	for(i=0;i<t->types;i++)  
		if(t->multiplicity[i]>1&&p< *(q++)) {
			*pflag=COALESCE; 
			if(t->n==2) m_csite[0]=t->color[i];
			if(debug==1) printf("\n\tDecreased multiplicity %d\n",i);
  		coalesce_trace(i,t->multiplicity[i]);
			coalesce(t,i);
			for(si=0;si<theta_pts;si++) for(sj=0;sj<m_pts;sj++) 
				for(sk=0;sk<g_pts;sk++) 
					surface_h[si][sj][sk] *=
						lambda(t->color[i],surface_g[sk]*tim)/lambda(t->color[i],tim);
			return x;
	} /* end coalesce */
	for(i=0;i<t->singletons;i++) 
		if(p < *(q++)) {
			*pflag=MUTATE;
			/* The mutation that is being removed is the first on sequence
			number
				t->singleton_list[i]
			site number
				t->sequence[t->singleton_list[i]][0]
			of color
				t->color[t->singleton_list[i]]
			*/
			site=t->sequence[ t->singleton_list[i] ][0];
			seq=t->singleton_list[i];
			m_age[site] = tim;
			m_csite[site]= t->color[seq];
			if(debug==1) 
				printf("\n\tRemoved singleton %d, color %d",site,t->color[seq]);
			new=delete_singleton(t,seq);
			if(new != NULL)
				shift_trace(seq,1,new[0],new[1]);
			for(si=0;si<theta_pts;si++) for(sj=0;sj<m_pts;sj++) 
				for(sk=0;sk<g_pts;sk++) 
					if(margin_flag==0)
						surface_h[si][sj][sk] *=surface_theta[si]/theta;
			return x;
	} /* end mutate */
	if(n_subpop>1) for(j=0;j<t->types;j++)  for(c=0;c<n_subpop;c++) 
		if(t->color[j]!=c&&p< *(q++)) {
			*pflag=MIGRATE;
			if(debug==1) 
				printf("\n\tMigration, sequence %d, %d -> %d\n",j,t->color[j],c);
			if(matrix_surface==1) move_migration[t->color[j]][c]++;
			migrate(t,j,c);
			if(matrix_surface==0) for(si=0;si<theta_pts;si++) 
					for(sj=0;sj<m_pts;sj++) for(sk=0;sk<g_pts;sk++) 
					surface_h[si][sj][sk] *=surface_m[sj]/m;
			return x;
		} /* end migrate */
	printf("\n\tError in random move uniform=%lf, cumulate=%lf\n",
						p,*(q-1));
	exit(1);
	return -1;
}

void fix_tree(TREE *t,TREE *s,int *pflag,double frate, int a, int b, int c) {
	static int count=1;
	if(check_tree(t)==0) {
		printtree(t,1);
		printf("\n\t********************** ERROR **********************\n");
		exit(1);
	}
	printf("\n\t* Tree %d",count);
	if(*pflag&MUTATE) printf("\n\tMutate %d, site %d",a,b);
	if(*pflag&COALESCE) printf("\n\tCoalesce %d",a);
	if(*pflag&MIGRATE) printf("\n\tMigrate %d, %d -> %d",a,b,c);
	printf("\n\tRate= %.4e",frate);
	printtree(t,2);
	copy_tree(t,s);
	*pflag=0;
	count++;
}
	
void check_move(TREE* t,double theta, double m, char* fname) {
	TREE *s;
	int i,j,si,sj,sk,single,c,site,flag=0,col;
	double x,tim=0.0,p,a,b,s_a,s_b,fix_rate;
	tim=next_time(t,theta,m,0.0);
	x=f(t,theta,m,tim,&a,&b);
	printf("\n\tf = %.4e, a = %.4e, b = %.4e, time= %.4e",x,a,b,tim);
	flag=0;
	s=(TREE*)c_malloc(sizeof(TREE));
	get_tree(fname,s);
	for(i=0;i<t->singletons;i++) {
		site=t->sequence[ t->singleton_list[i] ][0];
		single=1;
		if(t->singleton_exterior[i] != -1)
			single += t->multiplicity[t->singleton_exterior[i]];
		fix_rate=0.5*single*theta*t_rate[site];
		flag=MUTATE;
		delete_singleton(t,t->singleton_list[i]);
		fix_tree(t,s,&flag,fix_rate,i,site,-1);
	}
	for(i=0;i<t->types;i++)  if(t->multiplicity[i]>1) {
		fix_rate=0.5*t->subpop[t->color[i]]*(t->multiplicity[i]-1.0)
							*lambda(t->color[i],tim)/subsize[t->color[i]];
		flag=COALESCE; 
		coalesce(t,i);
		fix_tree(t,s,&flag,fix_rate,i,-1,-1);
	}
	for(j=0;j<t->types;j++)  for(c=0;c<n_subpop;c++) 
		if((col=t->color[j]) != c) {
		fix_rate=(in_sample(t,j,c)+1.0)*mig(t->color[j],c,tim,m)
						*t->subpop[t->color[j]]/(t->subpop[c]+1.0);
		flag=MIGRATE;
		migrate(t,j,c);
		fix_tree(t,s,&flag,fix_rate,j,col,c);
		}
	printf("\n\n");
	exit(0);
}
				
void print_move_migration() {
	int i,j;
	printf("\n");
	for(i=0;i<n_subpop;i++) {
		for(j=0;j<n_subpop;j++) printf(" %d",move_migration[i][j]);
		printf("\n");
	}
}
				
double simulate(TREE* t,double theta,double m,double *ptmrca) {
	double current_time=0.0,last_time=0.0,fnl=1.0,temp;
	int flag=0,count=0,si,sj,sk,i,j,*a;
	for(si=0;si<theta_pts;si++) for(sj=0;sj<m_pts;sj++) 
		for(sk=0;sk<g_pts;sk++) surface_fnl[si][sj][sk]=1.0;
	if(matrix_surface==1) 
	 for(i=0;i<n_subpop;i++) 
		for(j=0;j<n_subpop;j++) move_migration[i][j]=0;
	if(out_level>=5) for(i=0;i<n_subpop;i++) {
		if(t->subpop[i]>0) { 
			m_tmrca_c[i]=(-1.0); m_mrca_c[i]=(-1); mrca_test[i]=0;
		}
		else  { 
			m_mrca_c[i]=i; m_tmrca_c[i]=0.0; mrca_test[i]=1;
		}
	}
	while(t->n>1 && count < maxsim) {
	  if(debug==1) {
			printf("\n\t###########################################");
			printtree(t,1);
			print_trace(t);
		}
		current_time=next_time(t,theta,m,last_time);
		fnl *= move(t,theta,m,current_time,last_time,&flag);
		for(si=0;si<theta_pts;si++) for(sj=0;sj<m_pts;sj++) 
			for(sk=0;sk<g_pts;sk++) 
				surface_fnl[si][sj][sk] *= surface_h[si][sj][sk];
		if(out_level>=5) for(i=0;i<n_subpop;i++) {
			a=mrca_query(i);
			if(a!=NULL) {
				m_mrca_c[i]=t->color[a[0]];
				m_tmrca_c[i]=current_time;
				if(debug==1)
					printf("\n\tMRCA of subpopulation %d at %.4lf",
						i,current_time);
			}
		}
		last_time=current_time;
		count++;
	}
	if(count >= maxsim) {
		printf("\n\tMaximum steps %d in simulation exceeded, use -x option.\n",
			maxsim);
		exit(1);
	}
	*ptmrca=current_time;
	if(debug==1) {
		printf("\n\tTotal functional is %.4e",fnl);
		printf("\n\t*******************************************");
	}
	if(matrix_surface==1) {
		for(sj=0;sj<m_pts;sj++)  {
			new_migration_matrix(sj);
			temp=1.0;
			for(i=0;i<n_subpop;i++) for(j=0;j<n_subpop;j++) {
				if(move_migration[i][j]>0 && new_mig[i][j]!=mig_matrix[i][j])
					temp *= 
						exp(move_migration[i][j]
									*log(new_mig[i][j]/mig_matrix[i][j]));
			}
			for(si=0;si<theta_pts;si++) for(sk=0;sk<g_pts;sk++) 
				surface_fnl[si][sj][sk] *= temp;
		}
	}
	return fnl;
}

FILE *p_s_out;
FILE *t_s_out;

char *psformat0="%.4e %.4e %4e %.4e %.4e";
char *psformat1="%.4e %3.0lf  %4e %.4e %.4e";

void print_surface(double xrep) {
	int si,sj,sk;
	double se=0.0;
	char *p;
	if(matrix_surface==0) p=psformat0;
	else p=psformat1;
	for(sk=0;sk<g_pts;sk++) {
		for(sj=0;sj<m_pts;sj++) {
			for(si=0;si<theta_pts;si++) {
				if(n_sample==2&&two_stop==1) se=0.0;
				else se=
					sqrt((sssurface_fnl[si][sj][sk]
					-ssurface_fnl[si][sj][sk]*ssurface_fnl[si][sj][sk]/xrep)
						/(xrep*(xrep-1.0)));
				fprintf(p_s_out,p,
					surface_g[sk],
					surface_m[sj],
					surface_theta[si],
					ssurface_fnl[si][sj][sk]/xrep,
					se);
		if(theta_pts > 1) fprintf(p_s_out,"\n");
			}
		if(m_pts > 1) fprintf(p_s_out,"\n");
		}
		if(g_pts > 1) fprintf(p_s_out,"\n");
	}
		fclose(p_s_out);
}

int print_tmrca_surface(void) {
	int si,sj,sk;
	double se=0.0,mean;
	if(two_stop==1) return 0;
	for(sk=0;sk<g_pts;sk++) {
		for(sj=0;sj<m_pts;sj++) {
			for(si=0;si<theta_pts;si++) {
				mean=tm_surface[si][sj][sk]/ssurface_fnl[si][sj][sk];
				if(n_sample==2&&two_stop==1) se=0.0;
				 se= sqrt(ttm_surface[si][sj][sk]/ssurface_fnl[si][sj][sk]
						-mean*mean);
				fprintf(t_s_out,"%.4e %.4e %4e %.4e %.4e",
					surface_g[sk],
					surface_m[sj],
					surface_theta[si],
					mean,
					se);
		if(theta_pts > 1) fprintf(t_s_out,"\n");
			}
		if(m_pts > 1) fprintf(t_s_out,"\n");
		}
		if(g_pts > 1) fprintf(t_s_out,"\n");
	}
		fclose(t_s_out);
		return 1;
}

/* Collate material from multiple file input */
#define MFENTRIES 5
#define SENTRIES 5
char *mf_str="%.4e %.4e %.4e %.4e %.4e";
char *mf1_str="%.4e %.4e %.4e";
char *unroot="Unroot";

int next_name(char *str) {
	char *p;
	int next;
	p=strchr(str,'.');
	if(p==NULL) {
		printf("\n\tFilename %s should contain a dot\n",str);
		exit(1);
	}
	p++;
	if(*p=='\000') {
		*p='0';
		*(p+1)='\000';
	}
	else {
		next=atoi(p);
		next++;
		sprintf(p,"%d%c",next,'\000');
	}
	return next;
}

void s_collate(double xrep) {
	double x[SENTRIES];
	char fname[100];
	FILE *in;
	int i,j,si,sj,sk,flag;
	for(sk=0;sk<g_pts;sk++) 
		for(sj=0;sj<m_pts;sj++) 
			for(si=0;si<theta_pts;si++) {
					ssurface_fnl[si][sj][sk]=0.0;
					sssurface_fnl[si][sj][sk]=0.0;
	}
	for(i=0;i<100;i++) fname[i]='\000';
	strcpy(fname,s_name);
	*(strchr(fname,'.')+1)='\000';
	for(i=0;i<n_m_files;i++) {
		next_name(fname);
		in=fopen(fname,"r");
		if(in==NULL) {
			printf("\n\tCan't open %s for reading\n",fname);
			exit(1);
		}
		for(sk=0;sk<g_pts;sk++) {
			for(sj=0;sj<m_pts;sj++) {
				for(si=0;si<theta_pts;si++) {
					flag=fscanf(in,"%lf %lf %lf %lf %lf",x,x+1,x+2,x+3,x+4);
					if(flag!=5) {
						printf("\n\tProblem collating surface files\n");
						exit(1);
					}
					ssurface_fnl[si][sj][sk] += x[3]*xrep;
					sssurface_fnl[si][sj][sk] += 
									x[4]*x[4]*(xrep-1.0)*xrep+ xrep*x[3]*x[3];
					/* Debug only
					printf("%.4e %.4e %4e %.4e %.4e\n",x[0],x[1],x[2],x[3],x[4]);
					*/
				}
			}
		}
		fclose(in);
	}
	*(strchr(fname,'.')+1)='\000';
	strcat(fname,"all");
	p_s_out=fopen(fname,"w");
	if(p_s_out==NULL) {
		printf("\n\tCan't open surface collate file %s\n",fname);
		exit(1);
	}
	print_surface(xrep);
}

					
void print_result(int i,double *a) {
/*  xrep,sfnl,ssfnl,stmrca,sstmrca for multiple files in a[],
    a[0] a[1] a[2]  a[3]   a[4]
		a[5] is total functional
*/
	if(i>=0) printf("\n\t%6d ",i);else printf("\n\tUnroot ");
	if(two_stop==0) {
		printf(mf_str,
			a[1]/a[0],
			sqrt( (a[2] -a[1]*a[1]/(a[0]*a[0]))/(a[0]-1.0) ),
			a[1]/a[5],
			a[3]/a[1],
			sqrt( a[4]/a[1] -a[3]*a[3]/(a[1]*a[1]))
		);
	}
	else if(two_stop==1) {
		printf(mf1_str,
			a[1]/a[0],
			sqrt( (a[2] -a[1]*a[1]/(a[0]*a[0]))/(a[0]-1.0) ),
			a[1]/a[5]
		);
	}
}

int collate(void) {
	double mf_result[MFENTRIES+1], u_result[MFENTRIES+1],total_fnl=0.0;
	double *mf_fnl; /* Multiple file functional values */
	FILE *in;
	int i,j,flag;
	if(mf_out_level<1) return 0;
	mf_fnl=(double*)c_malloc(n_m_files*sizeof(double));
	in=fopen(MFTEMP,"r");
	if(in==NULL) {
		printf("\n\tCan't open %s for reading\n,MFTEMP");
		exit(1);
	}
	/* 1st reading to get total functional */
	for(i=0;i<n_m_files;i++) {
		for(j=0;j<MFENTRIES;j++) {
			flag=fscanf(in,"%lf",mf_result+j);
			if(flag!=1) {
				printf("\n\tProblem reading %s\n",MFTEMP);
				exit(1);
			}
		}
		total_fnl += mf_result[1];
	}
	u_result[MFENTRIES]=total_fnl;
	mf_result[MFENTRIES]=total_fnl;
	rewind(in);
	/* Print headings */
	printf("\n\tMultiple tree summary");
	printf("\n\t  Tree ");
	printf("  Like     ");
	printf("  SD Like  ");
	printf("  Rel Like ");
	if(two_stop==0) {
		printf("  TMRCA    ");
		printf("  SD TMRCA ");
	}
/* in contains  xrep,sfnl,ssfnl,stmrca,sstmrca for multiple files */
	for(j=1;j<MFENTRIES;j++) u_result[j] = 0.0;
	for(i=0;i<n_m_files;i++) {
		for(j=0;j<MFENTRIES;j++) {
			flag=fscanf(in,"%lf",mf_result+j);
			if(flag!=1) {
				printf("\n\tProblem reading %s\n",MFTEMP);
				exit(1);
			}
		}
		for(j=1;j<MFENTRIES;j++) u_result[j] += mf_result[j];
		mf_fnl[i] = mf_result[1];
		print_result(i,mf_result);
		u_result[0] = mf_result[0];
	}
	print_result(-1,u_result);
	printf("\n\n");
	fclose(in);
	unlink(MFTEMP);
	return 1;
}

void replicate(TREE *t, double theta, double m, long rep) {
	TREE *s;
	FILE *t_age;
	FILE *distn_out;
	int i,j,si,sj,sk,unique=1,pos,min_pos=tmrca_cells,max_pos=0;
	long r;
	double probability,fnl,sfnl=0.0,ssfnl=0.0,xrep,fnl_check,
		tmrca, stmrca=0.0, sstmrca=0.0, zero=0.0;
	xrep=(double)rep;
	s=(TREE*)c_malloc(sizeof(TREE));
	get_tree(tree_name,s);
	for(i=0;i<=n_sites;i++) {
		sage[i]=0.0;ssage[i]=0.0;
		for(j=0;j<n_subpop;j++) csite[i][j]=0.0;
	}
	for(i=0;i<n_subpop;i++) initial_subpop[i]=t->subpop[i];
	for(r=0;r<rep;r++) {
		copy_tree(s,t);
		if(out_level>=5) initial_trace(s);
		fnl = simulate(s,theta,m,&tmrca);
		if(r==0) fnl_check=fnl;
		unique &= (fnl == fnl_check);
		sfnl += fnl;
		ssfnl += fnl*fnl;
		stmrca += tmrca*fnl;
		sstmrca += tmrca*tmrca*fnl;
		m_age[0]=tmrca;
		/*XX*/
		if(two_stop==0) {
			if(out_level>=5) for(i=0;i<n_subpop;i++) {
				mrca_c[i][m_mrca_c[i]] += fnl;
				tmrca_c[i] += m_tmrca_c[i]*fnl;
				ss_tmrca_c[i] += m_tmrca_c[i]*m_tmrca_c[i]*fnl;
			}
			for(i=0;i<=n_sites;i++) {
				sage[i] += fnl*m_age[i];
				ssage[i] += fnl*m_age[i]*m_age[i];
				csite[i][m_csite[i]] += fnl;
			}
		} /* end two_stop = 0 */
		/*Surface calculation */
		for(si=0;si<theta_pts;si++) for(sj=0;sj<m_pts;sj++) 
			for(sk=0;sk<g_pts;sk++) {
				ssurface_fnl[si][sj][sk] += surface_fnl[si][sj][sk];
				sssurface_fnl[si][sj][sk] += 
					surface_fnl[si][sj][sk]*surface_fnl[si][sj][sk];
		}
		if(tmrca_surface==1) 
			for(si=0;si<theta_pts;si++) for(sj=0;sj<m_pts;sj++) 
				for(sk=0;sk<g_pts;sk++) {
				tm_surface[si][sj][sk] += tmrca*surface_fnl[si][sj][sk];
				ttm_surface[si][sj][sk] += 
					tmrca*tmrca*surface_fnl[si][sj][sk];
		}
		/* TMRCA distribution collection */
		if(tmrca_distn_flag==1) {
			if(tmrca>tmrca_top) pos=tmrca_cells;
				else pos = (int)floor((tmrca/tmrca_top)*tmrca_cells);
				if(pos < min_pos) min_pos=pos;
				if(pos > max_pos) max_pos=pos;
				tmrca_distn[pos] += fnl;
		}
		if(age_distn_flag==1) {
			for(i=1;i<=n_sites;i++) {
			if(m_age[i]>tmrca_top) pos=tmrca_cells;
			else pos = (int)floor((m_age[i]/tmrca_top)*tmrca_cells);
			if(pos < min_pos) min_pos=pos;
			if(pos > max_pos) max_pos=pos;
			if(age_distn[pos]==NULL) {
				age_distn[pos]=c_malloc((n_sites+1)*sizeof(double));
				for(j=0;j<=n_sites;j++) age_distn[pos][j]=0.0;
			}
			age_distn[pos][i] += fnl;
			}
		}
	}  /* end replicate loop */
	probability=sfnl/xrep;

	if(tmrca_distn_flag==1) 
		for(i=0;i<=tmrca_cells;i++) tmrca_distn[i] /= sfnl;
	if(age_distn_flag==1) 
		for(j=0;j<=tmrca_cells;j++) 
			if(age_distn[j]!=NULL) for(i=1;i<=n_sites;i++) 
				 age_distn[j][i] /= sfnl;


	if(multiple_files==1) {
		fprintf(mf_tmp,"\t%.15e %.15e %.15e %.15e %.15e\n",
			xrep,sfnl,ssfnl,stmrca,sstmrca);
		fclose(mf_tmp);
	}
	if(out_level>=1) {
		if(t->n==2&&two_stop==1) 
		printf("\n\tExact likelihood of the tree %.4e",probability);
		else printf("\n\tEstimated likelihood of the tree %.4e",probability);
		if(rep>1&&unique==0)
			printf(", se %.4e",sqrt((ssfnl-sfnl*sfnl/xrep)/(xrep*(xrep-1.0))));
	}
	if(out_level>=3) {
		printf("\n\tMean TMRCA %.4e",stmrca/sfnl);
		if(rep>1)
			printf(", sd %.4e",
				sqrt((sstmrca/sfnl-stmrca*stmrca/(sfnl*sfnl))));
	}
	if(out_level>=4) {
		t_age=fopen(tree_age,"w");
		if(t_age==NULL) {
			printf("\n\tProblem writing age file\n");
			exit(1);
		}
		printf("\n\tAges of mutations, sd");
		fprintf(t_age,"%3d %.4lf\n",0,stmrca/sfnl);
		for(i=1;i<=n_sites;i++)  {
			if(sage[i]>0.0) {
				printf("\n\t%3d  %.2e",i,sage[i]/sfnl);
				if(rep>1)
					printf("  %.2e", 
						sqrt((ssage[i]/sfnl-sage[i]*sage[i]/(sfnl*sfnl))));
				fprintf(t_age,"%3d %8.4e\n",i,sage[i]/sfnl);
			}
			if(sage[i]<=0.0) 
				fprintf(t_age,"%3d %8.4e\n",i,-1.0);
		}
		fclose(t_age);
		c_scale=gene_pop*generation;
/*
  	treepic(tree_name,tree_age,tree_ps,'p',
										12.0,14.0,17.0,25.0,0.8,2.4,1,1,c_scale);
*/
	}
	if(out_level>=5) {
		if(n_subpop>1) {
			printf(
			"\n\tSample MRCA distribution in subpopulations");
			printf("\n\t     ");
			for(j=0;j<n_subpop;j++) printf("%7d   ",j);
				printf("\n\t     ");
				for(j=0;j<n_subpop;j++) {
					if(csite[0][j]==0.0) printf("  %.2e",zero);
					else printf("  %.2e",csite[0][j]/sfnl);
				}
		}
	}
	if(n_subpop>1) {
		if(out_level>=7) {
			printf(
				"\n\tSubpopulation TMRCA");
			printf("\n\t     ");
			for(j=0;j<n_subpop;j++) printf("%7d   ",j);
				printf("\n\tMean ");
				for(j=0;j<n_subpop;j++)
					if(tmrca_c[j]==0.0) printf("  %.2e",zero);
					else printf("  %.2e",tmrca_c[j]/sfnl);
				if(rep>1) {
					printf("\n\tSD   ");
					for(j=0;j<n_subpop;j++)
						if(tmrca_c[j]==0.0) printf("  %.2e",zero);
						else printf("  %.2e",
							sqrt(ss_tmrca_c[j]/sfnl - tmrca_c[j]*tmrca_c[j]/(sfnl*sfnl)));
				}
			}
			if(out_level>=6) {
				printf(
				"\n\tSubpopulation MRCA distribution in subpopulations");
				printf("\n\t     ");
				for(j=0;j<n_subpop;j++) printf("%7d   ",j);
				for(i=0;i<n_subpop;i++) {
					printf("\n\t%3d  ",i);
					for(j=0;j<n_subpop;j++) {
					if(mrca_c[i][j]==0.0) printf("  %.2e",zero);
					else printf("  %.2e",mrca_c[i][j]/sfnl);
				}
			}
			printf(
				"\n\tMutation distribution in subpopulations (site x subpopulation)");
			printf("\n\t     ");
			for(j=0;j<n_subpop;j++) printf("%7d   ",j);
			for(i=1;i<=n_sites;i++) if(sage[i] >0.0) {
				printf("\n\t%3d  ",i);
				for(j=0;j<n_subpop;j++) {
					if(csite[i][j]==0.0) printf("  %.2e",zero);
					else printf("  %.2e",csite[i][j]/sfnl);
				}
			}
		} /* end out_level 6 */
	} /* end if sub_pop > 1 */
	if(surface==1) print_surface(xrep);
	if(tmrca_surface==1) print_tmrca_surface();
	/* Collate multiple files if done */
	/*
	if(multiple_files==1 && (atoi(strchr(tree_name,'.')+1)==n_m_files-1) ) {
		collate();
		if(surface==1) s_collate(xrep);
	}
	*/
	if(tmrca_distn_flag) {
		distn_out=fopen(distn_file,"w");
		if(distn_out==NULL) {
				printf("\n\tCan't open distribution file %s\n\n",distn_file);
				distn_out=stdout;
		}
		fprintf(distn_out,"# TMRCA distribution");
		if(age_distn_flag==1) 
			fprintf(distn_out," and age of mutation distributions");
		fprintf(distn_out,"\n#\n");
		for(i=min_pos;i<=max_pos;i++) {
			fprintf(distn_out,"%7.3lf %.4e",
							i*tmrca_top/tmrca_cells,tmrca_distn[i]);
		if(age_distn_flag==1) {
			if(age_distn[i]!=NULL) 
				for(j=1;j<=n_sites;j++) fprintf(distn_out," %.4e",age_distn[i][j]);
			else
				for(j=1;j<=n_sites;j++) fprintf(distn_out," %.4e",0.0);
		}
		fprintf(distn_out,"\n");
	}
	fclose(distn_out);
	}
}



void get_migration_rates(char *fname) {
	FILE *in;
	int i,j;
	in=fopen(fname,"r");
	if(in==NULL) {
		printf("\n\tCan't open migration rate file %s\n",fname);
		exit(1);
	}
	mig_matrix_flag=1;
	for(i=0;i<n_subpop;i++) for(j=0;j<n_subpop;j++) {
		if(fscanf(in,"%lf",&mig_matrix[i][j])!=1) {
			printf("\n\tProblem with migration rate file %s\n",fname);
			exit(1);
		}
	}
	fclose(in);
	for(i=0;i<n_subpop;i++) if(mig_matrix[i][i]>0.0) {
		printf("\n\tDiagonal migration rates must be zero in %s\n",fname);
		exit(1);
	}
}

void get_subpop(char *fname) {
	FILE *in;
	int i,j;
	in=fopen(fname,"r");
	if(in==NULL) {
		printf("\n\tCan't open subpopulation size file %s\n",fname);
		exit(1);
	}
	for(i=0;i<n_subpop;i++) 
		if(fscanf(in,"%lf",&subsize[i])!=1) {
			printf("\n\tProblem with subpopulation size file %s\n",fname);
			exit(1);
		}
	fclose(in);
}

void get_growth(char *fname) {
	FILE *in;
	int i,j;
	in=fopen(fname,"r");
	if(in==NULL) {
		printf("\n\tCan't open growth rates file %s\n",fname);
		exit(1);
	}
	growth_rate=(double*)c_malloc(n_subpop*sizeof(double));
	for(i=0;i<n_subpop;i++)  {
		if(fscanf(in,"%lf",&growth_rate[i])!=1) {
			printf("\n\tProblem with growth rate file %s\n",fname);
			exit(1);
		}
	}
	fclose(in);
}

/* For relative rates think of genes of length k sites,
   divided into r blocks of k_1, ..., k_r whose sum is k.
	 Suppose the rates in these blocks are theta_1, ... , theta_k.
	 Then fix the overall rate at
	    theta = k_1*theta_1 + ... + k_r*theta_r.
	The value in t_rate[i] is the ratio
	     k_j*theta_j/theta
	if site i belongs to block j.
	The rates don't add to 1, since different sites can be in the
	same block. In the homogeneous case there is only 1 block, so 
	all entries are 1.
*/
  
   
void get_t_rate(char *fname) {
	FILE *in;
	int i,j;
	in=fopen(fname,"r");
	if(in==NULL) {
		printf("\n\tCan't open mutations rates file %s\n",fname);
		exit(1);
	}
	for(i=1;i<=n_sites;i++) 
		if(fscanf(in,"%lf",&t_rate[i])!=1) {
			printf("\n\tProblem with mutation rates file %s\n",fname);
			exit(1);
		}
	for(i=1;i<=n_sites;i++) if(t_rate[i]>1.0) {
		printf("\n\tRelative mutation rates in %s must be <= 1.0\n",fname);
		exit(1);
	}
	fclose(in);
}

char *usage =
"\n\tUsage: genetree tree_file theta runs seed {options}\n";

char *version=
"\n\tAncestral inference for gene trees\
\n\tVersion %s, Bob Griffiths, %s\n";

char *option0 =
"\tOptions on command line or in a file (or both)\
\n\t {file containing options: @file_name}\
\n\t-j input sequences [sequence file name],\
\n\t-J input sequences [sequence file name] [ancestor file name],\
\n\t   (tree will be output in tree_file)\
\n\t-z generate all possible trees with different roots\
\n\t-e exponential growth rate file (array of rates)\
\n";
char *option1 =
"\t-s number of subpopulations (default - number in tree_file)\
\n\t-m migration rate file (matrix of rates, diagonals zero)\
\n\t-p subpopulation relative size file (array of rates)\
\n";
char *option2 =
"\t-f surface outfile_name (g m theta likelihood sd_like)\
\n\t-g theta surface domain [theta0 theta1 surface_points]\
\n\t-i growth magnification surface domain [g0 g1 surface_points]\
\n\t-h migration surface domain [m0 m1 surface_points]\
\n\t-H migration matrix surface [infile] [outfile]\
\n\t-k tmrca distribution [n_cells top_value outfile]\
\n\t-K tmrca and age distributions [n_cells top_value outfile]\
\n";
char *option3 =
"\t-o output level [value 0-7] (default 3)\
\n\n\t      1 likelihood\
\n\t      2 show tree\
\n\t      3 TMRCA\
\n\t      4 age of mutations\
\n\t      5 MRCA distribution in subpopulations\
\n\t      6 mutation distribution in subpopulations\
\n\t      7 subpopulation TMRCA\n\
\n\t-c tree picture file name (default: tree.ps)\
\n\t-C [gene population size] [generation size] (for tree picture)\
\n\t-2 calculate exact likelihood from 2 ancestors\
\n\t-l multiple file summary [0 or 1] (default 1)\
\n\t-b generate batch commands for multiple trees\
\n\t-P calculate average pairwise differences\
\n\t-x maximum events in one simulation run (default - 500)\
\n\t-y maximum types in the ancestry of the sample (default - 30)\
\n\n";

/*
\n\t-r relative mutation rate file for sites (array of rates)\
\n\t-a migration rate magnification [value]\
\n\t-t tmrca surface outfile_name\
\n\t-v printout of immediate ancestral configurations\
\n\t-w check integrity of data at each move (debug only)\
\n\t-d debug\
*/

void get_args( int argc, char **argv, double *pm) {
	int i,j;
	double s_theta0,s_theta1,s_m0,s_m1,s_g0,s_g1;
	for(i=0;i<argc;i++) if(argv[i][0]=='-') {
		switch(argv[i][1]) {
			case 'H' : matrix_surface=1; 
								 strcpy(m_filein,argv[i+1]);
								 strcpy(m_fileout,argv[i+2]);
								 break;
			case 'b' : batch=1;break;
			case '2' : two_stop=1;break;
			case 'w' : checktree=1;break;
			case 'm' : get_migration_rates(argv[i+1]);break;
			case 'a' : *pm=atof(argv[i+1]);break;
			case 'c' : strcpy(tree_ps,argv[i+1]);break;
			case 'C' : gene_pop=atof(argv[i+1]);
								 generation=atof(argv[i+2]);
								 break;
			case 'o' : out_level=atoi(argv[i+1]);break;
			case 'p' : get_subpop(argv[i+1]);break;
			case 'e' : get_growth(argv[i+1]);growth=1;break;
			case 'r' : get_t_rate(argv[i+1]);break;
			case 'd' : debug=1;break;
			case 'Z' : multiple_files=1;
							  	n_m_files=atoi(argv[i+1]);break;
			case 'x' : maxsim=atoi(argv[i+1]);break;
			case 'y' : maxtype=atoi(argv[i+1]);break;
			case 'l' : mf_out_level = atoi(argv[i+1]);break;
			case 'f' : strcpy(s_name,argv[i+1]);
									break;
			case 't' : tmrca_surface=1; strcpy(tmrca_s_name,argv[i+1]);
									break;
			case 'g' : t_surface=1;
								 if(argc<i+4) {
								 	printf("\n\tError in -g option\n");
									exit(1);
								 }
								 s_theta0=atof(argv[i+1]);
								 s_theta1=atof(argv[i+2]);
								 theta_pts=atoi(argv[i+3]);
								 surface_theta=(double*)c_malloc(theta_pts*sizeof(double));
								 surface_theta[0]=s_theta0;
								 for(j=1;j<theta_pts;j++) 
										surface_theta[j]=surface_theta[j-1] 
											+ (s_theta1-s_theta0)/(theta_pts-1.0);
								 break;
			case 'h' : m_surface=1;
								 if(argc<i+4) {
								 	printf("\n\tError in -h option\n");
									exit(1);
								 }
								 s_m0=atof(argv[i+1]);
								 s_m1=atof(argv[i+2]);
								 m_pts=atoi(argv[i+3]);
								 surface_m=(double*)c_malloc(m_pts*sizeof(double));
								 surface_m[0]=s_m0;
								 for(j=1;j<m_pts;j++) 
										surface_m[j]=surface_m[j-1] 
											+ (s_m1-s_m0)/(m_pts-1.0);
								 break;
			case 'i' : g_surface=1;
								 if(argc<i+4) {
								 	printf("\n\tError in -i option\n");
									exit(1);
								 }
								 s_g0=atof(argv[i+1]);
								 s_g1=atof(argv[i+2]);
								 g_pts=atoi(argv[i+3]);
								 surface_g=(double*)c_malloc(g_pts*sizeof(double));
								 surface_g[0]=s_g0;
								 for(j=1;j<g_pts;j++) 
										surface_g[j]=surface_g[j-1] 
											+ (s_g1-s_g0)/(g_pts-1.0);
								 break;
			case 'k' : tmrca_cells=atoi(argv[i+1]);
								 tmrca_top=atof(argv[i+2]);
								 strcpy(distn_file,argv[i+3]);
								 tmrca_distn_flag=1;
								 break;
			case 'K' : tmrca_cells=atoi(argv[i+1]);
								 tmrca_top=atof(argv[i+2]);
								 strcpy(distn_file,argv[i+3]);
								 tmrca_distn_flag=1;
								 age_distn_flag=1;
								 break;
			case 'v' : check_m=1; break;
			default : break;
		}
	}
}

int two_stop_set(void) {
	if(two_stop==0) 
		return 0;
	else { /* Modify later */
		if(!(
					n_subpop==1
					||(n_subpop>1&&growth==0)
		)	) 
		{ two_stop=0; return 1; }
		if(two_stop==1&&out_level>2) {
			out_level=2; 
			printf("\n\tSetting output level to 2 since -2 option in used\n");
		}
	}
	return 1;
}

char *info=
	"\tDate %d-%d-19%02d, time %02d:%02d:%02d\n";


int more(void) {
	int c;
	c='\000';
	printf("\t   [cr for more; space-cr to quit]\n");
	c=getchar();
	if(c=='\n') return 1;
	else exit(0);
}

void do_options(void) {
	printf(option0);
	if(more())	printf(option1);
	if(more())	printf(option2);
	if(more())	printf(option3);
}


int main(int argc, char **argv) {
	TREE *t;
	char *p,str[20],**fargv,recall[500];
	int i,j,k,flag,d[3],fargc,done=0;
	long rep;
	FILE *matrixf;
	double theta,m=1.0,s_theta0,s_theta1,s_m0,s_m1,s_g0,s_g1,xrep;
	time_t start_time,end_time;
	struct tm *timeptr;
	time_t secsnow;
	if(argc>1) strcpy(tree_name,argv[1]); /* We at least need a tree */
	/* args that don't require all input */
	/* Set maximum number of types first */
	for(i=0;i<argc;i++) if(argv[i][0]=='-'&&argv[i][1]=='y') 
			maxtype=atoi(argv[i+1]);
	for(i=0;i<argc;i++) if(argv[i][0]=='-') {
		switch(argv[i][1]) {
			case 's' : set_n_subpop=atoi(argv[i+1]); break;
			case 'j' : do_tree(argv[i+1],argv[1],NULL); 
								 if(argc<9) done=1; break; 
			case 'J' : do_tree(argv[i+1],argv[1],argv[i+2]); 
								 if(argc<10) done=1; break; 
			case 'z' : 
				t=(TREE*)c_malloc(sizeof(TREE));
				get_tree(tree_name,t);
				strcpy(mainf,tree_name);
				p=strchr(mainf,'.');
				if(p==NULL) {
					p=mainf+strlen(mainf);
					p[0]='.';
					p[1]='\000';
				} else {
					p[0]='_';
					p = mainf + strlen(mainf);
					p[0]='.';
					p[1]='\000';
				}
				p=mainf;
				for(i=0;i<=n_sites;i++) newtree(t,i,p);
 		   	done=1;
			 	break;
			case 'P' : 
				t=(TREE*)c_malloc(sizeof(TREE));
				get_tree(tree_name,t);
				printtree(t,0);
				printf("\n");
			 	pairwise(t);
				exit(0);
			default : break;
		}
	}
	if(done==1) exit(0);
	if(argc<5) {
		printf(usage);
		printf(version,ver,g_date);
		do_options();
		exit(1);
	}
	/* File option input */
	fargv=(char**)c_malloc(FARGS*sizeof(char*));
	for(i=0;i<FARGS;i++)
		fargv[i]=(char*)malloc((FARGS_SIZE+1)*sizeof(char));
  fargc=get_file_arg(NULL,argc,argv,fargv);

	/* Possibly call with multiple files tree.0, tree.1, ... */
	if(strchr(tree_name,MSYMBOL)!=NULL)  {
		for(i=0;i<20;i++) str[i]='\000';
		n_m_files=n_files(tree_name,MSYMBOL);
		strcpy(str," -Z");
		sprintf(str+3,"%3d",n_m_files);
		for(i=5;i<argc;i++) if(argv[i][0]=='-' && argv[i][1]=='b')  {
			argv[i][0]=' ';argv[i][1]=' ';
  		multiple_command(argc,argv,MSYMBOL,4,NULL,str,n_m_files,0,0);
			exit(3);
		}
  	multiple_command(argc,argv,MSYMBOL,4,NULL,str,n_m_files,0,1);
		exit(2);
	}


 	theta=atof(argv[2]);
	margin_flag=(theta < THETAZERO);
	rep=atol(argv[3]);
	xrep=(double)rep;
	if(argv[4][0]=='-') {
		printf("\n\tCommand line error: Random seed forgotton ?\n");
		exit(1);
	}
	dum = -(atol(argv[4]));
	t=(TREE*)c_malloc(sizeof(TREE));
	/* Check memory order here -maybe need args 1st */
	get_tree(tree_name,t);
	get_memory();
	get_args(fargc,fargv,&m);
	get_args(argc,argv,&m);
	if(n_subpop>1&&mig_matrix_flag==0) {
		printf("\n\tA backward migration matrix is needed: -m filename\n\n");
		exit(1);
	}
	if(g_surface==1&&growth==0) {
		printf("\n\tA likelihood surface with -i requires -e filename\n\n");
		exit(1);
	}
	initial_trace(t);
	if(check_m==1)	check_move(t,theta,m,tree_name);
	if(debug==1) { out_level=1;mf_out_level=1; }
	if(out_level>=2) {
		printtree(t,0);
		printf("\n");
	}
	if(matrix_surface==1) {
		if(n_subpop==1) {
			printf("\n\tCan't use -H with a single population\n");
			exit(1);
		}
		m_surface=0;
	}
	if(margin_flag==1) two_stop=0;
	two_stop_set();
	if(matrix_surface==1) { 
		read_migration(m_filein);
	 surface_m=(double*)c_malloc(m_pts*sizeof(double));
	 for(j=0;j<m_pts;j++) surface_m[j]=j;
		k=0;
		matrixf=fopen(m_fileout,"w");
		if(matrixf==NULL) {
			printf("\n\tCan't open matrix file %s for writing\n\n",m_fileout);
			exit(1);
		}
		fprintf(matrixf,"Migration matrices in migration surface\n");
		while(new_migration_matrix(k)!=(-1)) {
			fprintf(matrixf,"\nMatrix %d\n",k);
			for(i=0;i<n_subpop;i++) {
				for(j=0;j<n_subpop;j++) fprintf(matrixf,"%.4lf ",new_mig[i][j]);
				fprintf(matrixf,"\n");
			}
		k++;
		}
		fclose(matrixf);
	} /* end matrix_surface==1 */
	/* Surface housekeeping */
	surface=t_surface==1||m_surface==1||g_surface==1||matrix_surface==1;
	if(surface==1) {
		if(s_name[0]=='\000') {
			printf("\n\tSurface options need  -f surface outfile_name\n");
			exit(1);
		}
		/* Check on surface outfile for multiple files */
	 if(multiple_files==1&&s_name[0]!='\000'&&strncmp(strchr(s_name,'.'),
 		strchr(tree_name,'.'),4)!=0) {
		printf("\n\tUse symbol pattern # in the -f file option\n");
		printf("\texample:  myoutfile.#\n");
		exit(1);
	}
	p_s_out=fopen(s_name,"w");
	if(p_s_out==NULL) {
		printf("\n\tCan't open  %s\n\n",s_name);
		exit(1);
	}
		t_s=theta;
		m_s=m;
		g_s=1.0;
		if(theta_pts==0) { theta_pts=1; surface_theta= &t_s;}
		if(m_pts==0) { m_pts=1; surface_m= &m_s;}
		if(g_pts==0) { g_pts=1; surface_g= &g_s;}
		d[0]=theta_pts;d[1]=m_pts;d[2]=g_pts;
		for(i=0;i<3;i++) if(d[i]==0) d[i]=1;
		get_surface_mem(d);
	}
	if(multiple_files==1) {
		if( *(strchr(tree_name,'.')+1) == '0')
		mf_tmp=fopen(MFTEMP,"w"); else mf_tmp=fopen(MFTEMP,"a");
		if(mf_tmp==NULL) {
			printf("\n\tProblem opening %s\n",MFTEMP);
			exit(1);
		}
	}
	if(two_stop==1&&n_subpop>1&&growth==0) initial_gf(theta,m);
	if(tmrca_surface==1) {
		t_s_out=fopen(tmrca_s_name,"w");
		if(t_s_out==NULL) {
			printf("\n\tCan't open  %s\n\n",tmrca_s_name);
			exit(1);
		}
	}
	if((out_level <= 2&&growth==0)||two_stop==1&&growth==0) notime=1;
	if(maxsim<t->n) maxsim=2*t->n;
	if(tmrca_distn_flag==1) {
		tmrca_distn=(double*)c_malloc((tmrca_cells+1)*sizeof(double));
		for(i=0;i<=tmrca_cells;i++) tmrca_distn[i]=0.0;
	}
	if(age_distn_flag==1) {
		age_distn=(double**)c_malloc((tmrca_cells+1)*sizeof(double*));
		for(i=0;i<=tmrca_cells;i++) {
			/* Do the memory allocation when values are put in cells */
			/*
			age_distn[i]=c_malloc((n_sites+1)*sizeof(double));
			for(j=0;j<=n_sites;j++) age_distn[i][j]=0.0;
			*/
			age_distn[i]=NULL;
		}
	}
#include "test.h" /* File for debug only */
	time(&start_time);
 	replicate(t,theta,m,rep);
	time(&end_time);
	if(out_level>=2) {
		printf("\n\tElapsed time = %ld sec\n",end_time-start_time);
		printf("\tgenetree");
		for(i=1;i<argc;i++) printf(" %s",argv[i]);
		printf("\n");
		if(fargc>0) {
			printf("\tfile: ");
			for(i=0;i<fargc;i++) printf(" %s",fargv[i]);
			printf("\n");
		}
	}
	if(n_subpop>1&&mig_matrix_flag==1) {
		printf("\n\tBackward migration rate matrix");
		for(i=0;i<n_subpop;i++) {
			printf("\n\t");
			for(j=0;j<n_subpop;j++) printf("%.3lf ",mig_matrix[i][j]);
		}
		printf("\n\n");
	} 
	printf("\tSoftware version %s, %s\n",ver,g_date);
	tzset();
	time(&secsnow);
	timeptr=localtime(&secsnow);
	if(out_level>=2) {
		printf(info,
			(timeptr->tm_mon)+1,timeptr->tm_mday,timeptr->tm_year,
					timeptr->tm_hour,timeptr->tm_min,timeptr->tm_sec);
	}
	/* Collate multiple files if done */

	if(multiple_files==1 && (atoi(strchr(tree_name,'.')+1)==n_m_files-1) ) {
		collate();
		if(surface==1) s_collate(xrep);
	}
	return 0;
}

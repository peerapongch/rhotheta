#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <memory.h>

#ifdef PSVIEW 
char psview[100]="ghostview";
#endif

char *treepic_ver="\n\tBob Griffiths Ver 6.00, 28/03/98\n";
double xpaper=17.0;
double ypaper=25.0;
double tmrca;
double c_tmrca;
double udy;
char f_str[100];
char f_check[100];
#define MAXNAME 50
int uaxis=0;
int aaxis=0;
int look=0;
int sub_offset=0;
int sub_matrix=0;
int axis_offset;
int lineage_label=0;
double fontsize=10.0;
double sfontsize=10.0;
double mfontsize=10.0;
double lwidth=0.8;
double ps_radius=2.4;
double bl[2],br[2];
double translate[2];
/* Image size all in cm */
double xsize=12.0;
double ysize=14.0;
char percent='%';
int axis=1;
int dot=1;
int label=1;
int t_landscape=0;
int mult=1;
int equal=0;
int example=0;
int age_flag=0;
int tex_flag=0;
int ps_flag=1;
int bare_flag=0;
int mxlength=50;
int types=0;
int nsites=0;
int *order_list;
int *multiplicity;
int *margin;
int **sequence=NULL;
int **descend; /* nsites+1 x nsites+1 size */
double *x_coordinate;
double *x_leaf;
double *agev;
int **pop_seq;
int npop=0;
int *color;
FILE *pic;
FILE *age;
char **site_name;
char **line_name;
char **subpop_name;
char site_fname[50];
char line_fname[50];
char subpop_fname[50];
int *label_name=NULL;
int m_label=0,l_label=0,s_label=0;

char *info_tree=
	"Date %d-%d-19%02d, time %02d:%02d:%02d\n";
char *create="Creator: Treepic; ";

#ifdef CROSS
void cross(double x, double y, double delta) {
	fprintf(pic,"%.4lf %.4lf moveto\n",x-delta,y);
	fprintf(pic,"%.4lf %.4lf lineto\n",x+delta,y);
	fprintf(pic,"%.4lf %.4lf moveto\n",x,y-delta);
	fprintf(pic,"%.4lf %.4lf lineto\n",x,y+delta);

}
#endif

void preamble() {
		double x0,dx,dy,delta=15.0,y,df;
		int i;
		struct tm *timeptr;
		time_t secsnow;
		tzset();
		time(&secsnow);
		timeptr=localtime(&secsnow);
	if(tex_flag==1) {
		x0=1.0+(1.0/(types-1.0));
		fprintf(pic,"%c%c %s ",percent,percent,create);
		fprintf(pic,info_tree,
			(timeptr->tm_mon)+1,timeptr->tm_mday,timeptr->tm_year,
					timeptr->tm_hour,timeptr->tm_min,timeptr->tm_sec);
		fprintf(pic,"\\magnification=\\magstep1\n");
		fprintf(pic,"\\hsize=%.4lf true cm\n",xpaper);
		fprintf(pic,"\\vsize=%.4lf true in\n",ypaper);
		fprintf(pic,"\\nopagenumbers\n");
		fprintf(pic,"\\input pictex\n");
		fprintf(pic,"\\font\\small=cmr9\n");
		fprintf(pic,"\\bigskip\n");
		fprintf(pic,"%c\\centerline{\\bf Heading1}\n%c\\bigskip\n",
										percent,percent);
		fprintf(pic,"%c\\centerline{\\bf Heading2}\n%c\\bigskip\n",
										percent,percent);
		fprintf(pic,"\\centerline\n");
		fprintf(pic,"{\n");
		fprintf(pic,"\\beginpicture\n");
		fprintf(pic,
			"\\setcoordinatesystem units <%.4lf true cm, %.4lf true cm>\n",
			xsize,ysize);
		fprintf(pic,"\\setplotarea x from 0 to %.2lf, y from 0.0 to 1.0\n",x0);
		fprintf(pic,"\\linethickness=%.2lf pt\n",lwidth);
		if(axis==1) {
			fprintf(pic,"\\axis right ticks\n");
			fprintf(pic,"  numbered from 0.0 to 1.0 by 0.1 /\n");
		}
	}
	if(ps_flag==1) {
	  translate[0]=(xpaper-xsize)*0.5*28.45+sub_offset*3.0*ps_radius;
		translate[1]=(ypaper-ysize)*0.5*28.45+npop*6.0*ps_radius;
		bl[0]= 0.0;bl[1]= 0.0;
		br[0]=xsize*28.45*(1.0+1.0/(types-1.0));
		br[1]=ysize*28.45;
		bl[0] += translate[0];
		br[0] += translate[0];
		bl[1] += translate[1];
		br[1] += translate[1];
		bl[0] -= 60+sub_offset*3.0*ps_radius;
		br[0] += 60;
		bl[1] -= 60+npop*6.0*ps_radius;
		br[1] += 30;
		fprintf(pic,"%c!PS-Adobe-3.0 EPSF-3.0\n",percent);
		fprintf(pic,"%c%c %s ",percent,percent,create);
		fprintf(pic,info_tree,
			(timeptr->tm_mon)+1,timeptr->tm_mday,timeptr->tm_year,
					timeptr->tm_hour,timeptr->tm_min,timeptr->tm_sec);
		fprintf(pic,"%c%cBoundingBox: %.2lf %.2lf %.2lf %.2lf\n",
				percent,percent,bl[0],bl[1],br[0],br[1]);
		fprintf(pic,"%c%cEndcomments\n",percent,percent);
		fprintf(pic,"%c%cEndProlog\n%c%c\n",
						percent,percent,percent,percent);
		fprintf(pic," /rightshow\n { \n dup stringwidth pop");
		fprintf(pic,"\n 4 -1 roll exch sub\n 3 -1 roll moveto\n show } def\n");
		fprintf(pic,"%c%c\n",percent,percent);
		/* 1mm = 2.845 pt */
#ifdef CROSS
		/*
		cross(0,0,delta);
		cross(0,ypaper*28.45,delta);
		cross(xpaper*28.45,0,delta);
		cross(xpaper*28.45,ypaper*28.45,delta);
		*/
		fprintf(pic,"%c%c\n",percent,percent);
		fprintf(pic,"%c%c Bounding Box crosses\n",percent,percent);
		cross(bl[0],bl[1],delta);
		cross(br[0],br[1],delta);
		fprintf(pic,"%c%c\n",percent,percent);
#endif
		fprintf(pic,"%c%c\n",percent,percent);
		fprintf(pic,"%.4lf %.4lf translate\n",translate[0],translate[1]);
		fprintf(pic,"%.4lf %.4lf scale\n",xsize*28.45,xsize*28.45);
		if(t_landscape==1)
			fprintf(pic,"90 rotate\n");
		fprintf(pic,"%.4lf setlinewidth\n",lwidth/(14.225*(xsize+xsize)));
		/*
		fprintf(pic,"%c%c/Times-Roman findfont %.4lf scalefont setfont\n",
						percent,percent,fontsize/(28.45*xsize));
		*/
		fprintf(pic,"/Helvetica findfont %.4lf scalefont setfont\n",
										fontsize/(28.45*xsize));
		df=fontsize/(28.45*xsize);
		/* Now do a right axis */
		if(axis==1) {
			fprintf(pic,"%c%c\n%c%c Axis\n%c%c\n",
								percent,percent,percent, percent, percent, percent);
			fprintf(pic,"%.4lf setlinewidth\n",0.25*lwidth/(14.225*(xsize+xsize)));
			dx=1.5*ps_radius/(28.45*xsize);
			dy=ysize/(5.0*xsize);
			fprintf(pic,"%.4lf %.4lf moveto %.4lf %.4lf lineto\n",
				1.0+1.0/(types-1.0),0.0,1.0+1.0/(types-1.0),5.0*dy);
			for(i=0;i<=5;i++) {
					fprintf(pic,"%.4lf %.4lf moveto %.4lf %.4lf lineto\n",
					1.0+1.0/(types-1.0),i*dy,1.0+1.0/(types-1.0)+dx,i*dy);
				fprintf(pic,"%.4lf %.4lf moveto\n",1.0+1.0/(types-1.0)+2.0*dx,i*dy);
				fprintf(pic,"(%.1lf) show\n",i/5.0);
			}
		}

		if(uaxis==1||aaxis==1) {
			sprintf(f_check,f_str,tmrca);
			axis_offset=strlen(f_check);
		}
		if(uaxis==1) {
			fprintf(pic,
			"%c%c\n%c%c Axis (user) -edit for different format numbers\n%c%c\n",
								percent,percent,percent, percent, percent, percent);
			fprintf(pic,"%.4lf setlinewidth\n",0.25*lwidth/(14.225*(xsize+xsize)));
			dx=1.5*ps_radius/(28.45*xsize);
			dy=ysize*udy/(tmrca*xsize);
			fprintf(pic,"%.4lf %.4lf moveto %.4lf %.4lf lineto\n",
				1.0+1.0/(types-1.0),0.0,1.0+1.0/(types-1.0),ysize/xsize);
			i=0;y=0.0;
			while(y<tmrca && i< 25) {
					fprintf(pic,"%.4lf %.4lf moveto %.4lf %.4lf lineto\n",
					1.0+1.0/(types-1.0),i*dy,1.0+1.0/(types-1.0)+dx,i*dy);
				/*
				fprintf(pic,"%.4lf %.4lf moveto\n",1.0+1.0/(types-1.0)+2.0*dx,i*dy);
				*/
				fprintf(pic,"%.4lf %.4lf ",
								1.0+1.0/(types-1.0)+0.5*dx*axis_offset,i*dy);
				fprintf(pic,f_str,y);
				y += udy;
				i++;
			}
			/* tmrca */
			fprintf(pic,"%.4lf %.4lf moveto %.4lf %.4lf lineto\n",
					1.0+1.0/(types-1.0),ysize/xsize,1.0+1.0/(types-1.0)+dx,ysize/xsize);
			/*
			fprintf(pic,"%.4lf %.4lf moveto\n",
				1.0+1.0/(types-1.0)+2.0*dx,ysize/xsize);
			*/
			fprintf(pic,"%.4lf %.4lf ",
				1.0+1.0/(types-1.0)+0.5*dx*axis_offset,ysize/xsize);
			fprintf(pic,f_str,tmrca);
		} /* end uaxis */
		if(aaxis==1) {
		sfontsize=0.75*fontsize;
		fprintf(pic,"%c%c Edit axis labels here if they clash\n",percent,percent);
		fprintf(pic,"%c%c Smaller font size\n",percent,percent);
		fprintf(pic,"/Helvetica findfont %.4lf scalefont setfont\n",
										sfontsize/(28.45*xsize));
			fprintf(pic,"%.4lf setlinewidth\n",0.25*lwidth/(14.225*(xsize+xsize)));
			dx=1.5*ps_radius/(28.45*xsize);
			dy=ysize/xsize;
			fprintf(pic,"%.4lf %.4lf moveto %.4lf %.4lf lineto\n",
				1.0+1.0/(types-1.0),0.0,1.0+1.0/(types-1.0),ysize/xsize);
			for(i=1;i<=nsites;i++) if(agev[i]>= 0.0) {
					fprintf(pic,"%.4lf %.4lf moveto %.4lf %.4lf lineto\n",
					1.0+1.0/(types-1.0),agev[i],1.0+1.0/(types-1.0)+dx,agev[i]);
				/*
				fprintf(pic,"%.4lf %.4lf moveto\n",
					1.0+1.0/(types-1.0)+2.0*dx,agev[i]);
				*/
				fprintf(pic,"%.4lf %.4lf ",
					1.0+1.0/(types-1.0)+0.5*dx*axis_offset,agev[i]);
				/* Modified 9/11/97
				fprintf(pic,f_str,agev[i]*tmrca);
				*/
				fprintf(pic,f_str,agev[i]*tmrca*xsize/ysize);
			}
			/* tmrca */
			fprintf(pic,"%.4lf %.4lf moveto %.4lf %.4lf lineto\n",
					1.0+1.0/(types-1.0),ysize/xsize,1.0+1.0/(types-1.0)+dx,ysize/xsize);
			/*
			fprintf(pic,"%.4lf %.4lf moveto\n",
				1.0+1.0/(types-1.0)+2.0*dx,ysize/xsize);
			*/
			fprintf(pic,"%.4lf %.4lf ",
				1.0+1.0/(types-1.0)+0.5*dx*axis_offset,ysize/xsize);
			fprintf(pic,f_str,tmrca);
		fprintf(pic,"/Helvetica findfont %.4lf scalefont setfont\n",
										fontsize/(28.45*xsize));
		} /* end aaxis */
		fprintf(pic,"stroke\n");
		fprintf(pic,"%c%c\n%c%c Draw the tree\n%c%c\n",
			percent, percent, percent, percent, percent, percent);
	  fprintf(pic,"newpath\n");
		fprintf(pic,"%.4lf setlinewidth\n",lwidth/(14.225*(xsize+xsize)));
	}
		
}

void post() {
	if(tex_flag==1) 
		fprintf(pic,"\\endpicture\n}\n\\bigskip\n\\end\n");
	if(ps_flag==1) 
		fprintf(pic,"stroke\nshowpage\n%c%cEOF\n",percent,percent);
}


void	treepic_vertical_leaf(double x,int j) {
	if(bare_flag==1)
		fprintf(pic,"%.4lf %.4lf %.4lf m%03d\n",x,0.0,x,j);
	if(tex_flag==1)
		fprintf(pic,"\\plot %.4lf %.4lf %.4lf %.4lf /\n",
				x,0.0,x,agev[j]);
	if(ps_flag==1) 
		fprintf(pic,"%.4lf %.4lf moveto %.4lf %.4lf lineto\n",
				x,0.0,x,agev[j]);
}

void treepic_vertical(int j,int k) {
	if(bare_flag==1)
		fprintf(pic,"%.4lf m%03d %.4lf m%03d\n",
					x_coordinate[j],j,x_coordinate[j],order_list[k]);
	if(tex_flag==1)
		fprintf(pic,"\\plot %.4lf %.4lf %.4lf %.4lf /\n",
					x_coordinate[j],agev[j],x_coordinate[j],agev[order_list[k]]);
	if(ps_flag==1) 
		fprintf(pic,"%.4lf %.4lf moveto %.4lf %.4lf lineto\n",
					x_coordinate[j],agev[j],x_coordinate[j],agev[order_list[k]]);
}

void	horizontal(double xmin,double xmax,int k) {
	if(bare_flag==1)
		fprintf(pic,"%.4lf m%03d %.4lf m%03d\n",
			xmin,order_list[k],xmax,order_list[k]);
	if(tex_flag==1)
		fprintf(pic,"\\plot %.4lf %.4lf %.4lf %.4lf /\n",
			xmin,agev[order_list[k]],xmax,agev[order_list[k]]);
	if(ps_flag==1)
		fprintf(pic,"%.4lf %.4lf moveto %.4lf %.4lf lineto\n",
			xmin,agev[order_list[k]],xmax,agev[order_list[k]]);
}

void point(int j) {
	if(bare_flag==1)
		fprintf(pic,"%.4lf m%03d\n",x_coordinate[j],j);
	if(tex_flag==1)  {
			if(dot==0) fprintf(pic,"%c%c ",percent,percent);
			fprintf(pic,"\\put{$\\bullet$} at %.4lf %.4lf\n",
											x_coordinate[j],agev[j]);
	}
	if(ps_flag==1) {
			fprintf(pic,"%c%c Point %d\n",percent,percent,j);
			fprintf(pic,"%.4lf %.4lf moveto\n",
											x_coordinate[j],agev[j]);
		/* 1mm = 2.845 pt */
			if(dot==0) fprintf(pic,"%c%c ",percent,percent);
			fprintf(pic,"%.4lf %.4lf %.4lf  0 360 arc closepath fill\n",
											x_coordinate[j],agev[j],ps_radius/(28.45*xsize));
	}
}

void point_text(int j) {
	int flag=0;
	if(tex_flag==1) {
		if(j==1) fprintf(pic,"{ \\small\n");
		if(x_coordinate[j]> -100.0)  {
			flag=(descend[j][j]<=(-2));
			if(label==0) fprintf(pic,"%c%c ",percent,percent);
			if(margin[j]+flag<2)
				fprintf(pic,"\\put{%d} [Bl] ",j);
			else
				fprintf(pic,"\\put{%d} [Bl] <0 true mm, 1 true mm> ",j);
			fprintf(pic,"at %.4lf %.4lf\n",x_coordinate[j],agev[j]);
		}
		if(j==nsites) fprintf(pic,"}\n");
	}
	if(ps_flag==1) {
		/* 9 pt Helvetica */
		if(j==1) {
			fprintf(pic,"%c%c\n%c%c Site labels - edit for different labels\n", 
												percent,percent,percent,percent);
			fprintf(pic,"%c%c\n",percent,percent);
			/*
			fprintf(pic,"%c%c/Times-Roman findfont %.4lf scalefont setfont\n",
						percent,percent,0.75*fontsize/(28.45*xsize));
			*/
			fprintf(pic,"/Helvetica findfont %.4lf scalefont setfont\n",
										mfontsize/(28.45*xsize));
		}
		if(x_coordinate[j]> -100.0)  {
			flag=(descend[j][j]<=(-2));
			if(margin[j]+flag<2) {
				fprintf(pic,"%.4lf %.4lf moveto\n",
									x_coordinate[j]+1.5*ps_radius/(28.45*xsize),agev[j]);
			}
			else
				fprintf(pic,"%.4lf %.4lf moveto\n",
					x_coordinate[j]+1.5*ps_radius/(28.45*xsize),
					agev[j]+1.0*ps_radius/(28.45*xsize));
			if(label==0) fprintf(pic,"%c%c ",percent,percent);
			if(m_label==0) fprintf(pic,"(%d) show\n",j);
			else fprintf(pic,"(%s) show\n",site_name[j-1]);
		}
	}
}

void text_line(int i) {
	if(mult==0) return;

	if(bare_flag==1)
		fprintf(pic,"%4d %.4lf yl\n",
				multiplicity[i],x_leaf[sequence[i][0]]);
	if(tex_flag==1)
		fprintf(pic,"\\put{%d}  <-2.8 true mm, 0 true mm> at %.4lf -0.05\n",
				multiplicity[i],x_leaf[sequence[i][0]]);
	if(ps_flag==1) {
	if(i==0) {
		fprintf(pic,"%c%c\n", percent,percent);
		fprintf(pic,"%c%c Multiplicities\n", percent,percent);
		fprintf(pic,"%c%c\n", percent,percent);
		/*
			fprintf(pic,"%c%c/Times-Roman findfont %.4lf scalefont setfont\n",
						percent,percent,fontsize/(28.45*xsize));
		*/
			fprintf(pic,"/Helvetica findfont %.4lf scalefont setfont\n",
										fontsize/(28.45*xsize));
	}
	if(lineage_label==0)
		fprintf(pic,"%.4lf %.4lf\n",
				x_leaf[sequence[i][0]] +1.5*ps_radius/(28.45*xsize)
				,-6.0*ps_radius/(28.45*xsize));
	else
		fprintf(pic,"%.4lf %.4lf\n",
				x_leaf[sequence[i][0]] +1.5*ps_radius/(28.45*xsize)
				,-12.0*ps_radius/(28.45*xsize));
		fprintf(pic,"(%d) rightshow\n",multiplicity[i]);
	}
}


void lineage_line(int i) {
	if(lineage_label==0)return;
	if(ps_flag==1) {
		if(lineage_label==1) {
		if(i==0) {
			fprintf(pic,"%c%c\n%c%c Lineage labels - edit for different labels\n",
											percent,percent,percent,percent);
		fprintf(pic,"%c%c\n",percent,percent);
		fprintf(pic,"/Helvetica findfont %.4lf scalefont setfont\n",
										fontsize/(28.45*xsize));
		}
		/*
		fprintf(pic,"%.4lf %.4lf moveto\n",
				x_leaf[sequence[i][0]] -1.0*ps_radius/(28.45*xsize)
				,-6.0*ps_radius/(28.45*xsize));
		if(l_label==0) fprintf(pic,"(%c) show\n",i+'A');
		else fprintf(pic,"(%s) show\n",line_name[i]);
		*/
		fprintf(pic,"%.4lf %.4lf ",
				x_leaf[sequence[i][0]] +1.5*ps_radius/(28.45*xsize)
				,-6.0*ps_radius/(28.45*xsize));
		if(l_label==0) fprintf(pic,"(%c) rightshow\n",i+'A');
		else fprintf(pic,"(%s) rightshow\n",line_name[i]);
		}
	}
}


void subpopulation_matrix() {
	int i,j,k;
	double yoffset=0.0;
	double xoffset=0.0;
	double coffset=0.0;
	if(sub_matrix==0) return;
	if(ps_flag==1) {
			fprintf(pic,"%c%c\n%c%c Subpopulation matrix\n",
											percent,percent,percent,percent);
		fprintf(pic,"%c%c\n",percent,percent);
		fprintf(pic,"/Helvetica findfont %.4lf scalefont setfont\n",
										fontsize/(28.45*xsize));
		coffset=0.5*sub_offset*fontsize/(28.45*xsize);
		if(mult==1&&lineage_label==1) yoffset=(-12.0);
		if(mult==1&&lineage_label==0) yoffset=(-6.0);
		if(mult==0&&lineage_label==1) yoffset=(-6.0);
		if(mult==0&&lineage_label==0) yoffset=(0.0);
		yoffset -= 2.0;   /* Make the matrix slightly distinct */
		fprintf(pic,"%c%c\n%c%c Population Labels, edit for real names\n",
											percent,percent,percent,percent);
		fprintf(pic,"%c%c\n",percent,percent);
		for(i=0;i<npop;i++) {
	 		fprintf(pic,"%.4lf %.4lf moveto\n",-coffset-10.0*ps_radius/(28.45*xsize)
					,(yoffset-6.0*(i+1.0))*ps_radius/(28.45*xsize));
			if(s_label==0) fprintf(pic,"(%d) show\n",i);
			else fprintf(pic,"(%s) show\n",subpop_name[i]);
		}
		fprintf(pic,"%c%c\n%c%c Subpopulation frequencies\n",
											percent,percent,percent,percent);
		for(i=0;i<npop;i++) {
			for(j=0;j<types;j++) {
		 		fprintf(pic,"%.4lf %.4lf\n",
					x_leaf[sequence[j][0]] +1.5*ps_radius/(28.45*xsize)
					,(yoffset-6.0*(i+1.0))*ps_radius/(28.45*xsize));
				if(pop_seq[i][j]>0)
					fprintf(pic,"(%d) rightshow\n",pop_seq[i][j]);
				else fprintf(pic,"(%c) rightshow\n",'.');
			}
		}
	}
}

/* Version that reads in subpopulation matrices too */


void write_tr(char *total_name) {
	int i,j,k,l,count=0;
	int *s,*site_o;
	FILE *combine;
	s=(int*)malloc((nsites+1)*sizeof(int));
	site_o=(int*)malloc((nsites+1)*sizeof(int));
	for(i=0;i<=nsites;i++) s[i]=0;
	for(i=0;i<types;i++)  s[sequence[i][0]] += multiplicity[i];
	combine=fopen(total_name,"w");
	if(combine==NULL) {
		printf("\n\tCan't open a file for the total tree\n");
		exit(1);
	}
	/* Get the best order */
	k=0; site_o[0]=(-1);
	for(i=0;i<types;i++) {
		for(j=0;j<i;j++) if(sequence[i][0]==sequence[j][0]) break;
		if(j==i) site_o[k++]=sequence[i][0];
	}
	site_o[k]=(-1);
	for(l=0;site_o[l]!=(-1);l++) {
		i=site_o[l];
		for(j=0;j<types;j++) 
			if(sequence[j][0]==i) {
				pop_seq[color[j]][count]=multiplicity[j];
			}
			count++;
	}
	for(l=0;site_o[l]!=(-1);l++) {
		i=site_o[l];
		for(j=0;j<types;j++) 
			if(sequence[j][0]==i) {
				fprintf(combine,"%3d :",s[i]);
				for(k=0;sequence[j][k];k++) 
					fprintf(combine," %d",sequence[j][k]);
				fprintf(combine," 0\n");
				break;
			}
			count++;
	}
	/* Debug only  
		for(j=0;j<types;j++)  printf("%d ",color[j]);
		printf("\n\n");
	for(i=0;i<npop;i++) {
		for(j=0;j<types;j++) printf("%d ",pop_seq[i][j]);
		printf("\n");
	}
	*/
	fclose(combine);
	free(s);
	free(site_o);
}

int read_tree(char *filename) {
	FILE *in;
	FILE *out;
	int i,j,x,y,z,length=0,sub=0;
	char separator[20];
	char a[20],b[20],c[20];
	if((in=fopen(filename,"r"))==NULL) {
		printf("\n\tCan't open %s\n\n",filename);
		exit(1);
	}
	/* Which file is it ? */
	fscanf(in,"%s %s %s",a,b,c);
	if(c[0]==':') sub=1;
	rewind(in);
	types=0;
	while(feof(in) == 0) {
		if(fscanf(in,"%d",&x)==1) types++;
		if(sub==1) {
			if(x>npop) npop=x;
		fscanf(in,"%d",&x);
		}
		fscanf(in,"%s",separator);
		while(x != 0 && feof(in) == 0) {
			fscanf(in,"%d",&x);
			length++;
		}
		if(length>mxlength) mxlength=length;
			length =0;
	}
	if(sub==1) npop++;
	if(sequence==NULL) {
		sequence=(int**)malloc(types*sizeof(int*));
		for(i=0;i<types;i++)  
			sequence[i]=(int*)malloc(mxlength*sizeof(int));
		multiplicity=(int*)malloc(types*sizeof(int));
		if(sub==1) {
			color=(int*)malloc(types*sizeof(int));
			pop_seq=(int**)malloc(npop*sizeof(int*));
			for(i=0;i<npop;i++) {
				pop_seq[i]=(int*)malloc(types*sizeof(int));
				for(j=0;j<types;j++) pop_seq[i][j]=0;
			}
		}
	}
	rewind(in);
	if(types==1) {
		fclose(in);
		return 0;
	}
	for(i=0;i<types;i++) {
		if(sub==1) fscanf(in,"%d",color+i);
		fscanf(in,"%d",&x);
		multiplicity[i]=x;
		fscanf(in,"%s",separator);
		if(separator[0] != ':') {
			return 0;
		}
		j=0;
		while(x != 0 && feof(in) == 0) {
			fscanf(in,"%d",&x);
			sequence[i][j++]=x;
			if(x>nsites) nsites=x;
		}
	}
	rewind(in);
	fclose(in);
	if(sub==1) {
		write_tr("treepic.all");
		read_tree("treepic.all");
	}
	xsize=0.75*(types+1.0);
	ysize=1.6*xsize;
	return 1;
}


/* descend[i][j]=k if site i has site j as a (direct) child, and
	sites i and j are in the kth leaf lineage path to the root.
	 descend[i][i]=(-(k+2)) if i is a leaf node of kth sequence.
	 descend[i][j]=(-1) otherwise
*/

void get_descend() {
	int i,j,k;
	/*
	descend=(int**)malloc((nsites+1)*(nsites+1)*sizeof(int*));
	*/
	descend=(int**)malloc((nsites+1)*sizeof(int*));
	order_list=(int*)malloc((nsites+1)*sizeof(int));
	x_coordinate=(double*)malloc((nsites+1)*sizeof(double));
	x_leaf=(double*)malloc((nsites+1)*sizeof(double));
	for(i=0;i<=nsites;i++) {
		descend[i]=(int*)malloc((nsites+1)*sizeof(int));
		for(j=0;j<=nsites;j++) descend[i][j]=(-1);
	}
	for(i=0;i<types;i++) {
			descend[ sequence[i][0] ][ sequence[i][0] ]=(-(2+i));
		for(j=0;sequence[i][j];j++) 
			descend[ sequence[i][j+1] ][ sequence[i][j] ]=i;
	}
	for(i=0;i<=nsites+1;i++) {
		order_list[i]=(-1);
		x_coordinate[i]=(-100.0);
	}
	margin=(int*)malloc((nsites+1)*sizeof(int));
	for(i=0;i<=nsites;i++) {
		margin[i]=0;
		for(j=0;j<=nsites;j++) if(i!=j&&descend[i][j]>=0) margin[i]++;
	}
}

void insert_sort(int *a, int *b, int n) {
	int i,j,k,l;
	for(j=1;j<n;j++) {
		k=a[j];
		l=b[j];
		i=j-1;
		while(i>=0&&a[i]>k) {
			a[i+1]=a[i];
			b[i+1]=b[i];
			i--;
		}
		a[i+1]=k;
		b[i+1]=l;
	}
}

int order_count=0;
int leaf_count=0;
int *leaf_order;

void order_site(int i) {
	int j,k,l=0;
	int**site_sort;
	/*
	printf("%d %d\n",i,order_count);
	*/
	order_list[order_count]=i;
	order_count++;
	if(order_count>nsites+1) exit(1);
	site_sort=(int**)malloc(2*sizeof(int*));
	site_sort[0]=(int*)malloc((nsites+1)*sizeof(int));
	site_sort[1]=(int*)malloc((nsites+1)*sizeof(int));
	for(j=0;j<=nsites;j++) if(descend[i][j]>=0) {
			site_sort[0][l]=descend[i][j];
			site_sort[1][l]=j;
		l++;
	}
	if(descend[i][i]<=(-2)) {
		site_sort[0][l]= -descend[i][i] - 2;
		site_sort[1][l]= i;
		l++;
	}
	if(l>0) {
	/* Sort sites */
		insert_sort(site_sort[0],site_sort[1],l);
		for(j=l-1;j>=0;j--) {
			if(site_sort[1][j]==i) 
				leaf_order[leaf_count++]=i;
			else order_site(site_sort[1][j]);
		}
	}
	free(site_sort[0]);
	free(site_sort[1]);
	free(site_sort);
}

void leaf_coordinates() {
	double dx=0.1,x=0.0;
	int i,j;
	dx=1.0/(types-1.0);
	for(j=types-1;j>=0;j--) {
			x_leaf[leaf_order[j]] = x;
			/* Vertical leaf lines */
			treepic_vertical_leaf(x,leaf_order[j]);
			x += dx;
	}
}

void get_coordinate(int k) {
	int j;
	double xmin=1.0e20,xmax=0.0,xc,temp;
		for(j=0;j<=nsites;j++) if(x_coordinate[j]>=0.0) {
			if(descend[order_list[k]][j]>=0) {
				xc=x_coordinate[j];
				xmax= xmax<xc ? xc:xmax;
				xmin= xmin>xc ? xc:xmin;
				treepic_vertical(j,k);
				/* Vertical line */
				/*
				fprintf(pic,"%.4lf m%03d %.4lf m%03d\n",
					x_coordinate[j],j,x_coordinate[j],order_list[k]);
				*/
			}
		}
		if(descend[order_list[k]][order_list[k]]<=(-2)) {
			if(xmin>x_leaf[order_list[k]]) xmin=x_leaf[order_list[k]];
			if(xmax<x_leaf[order_list[k]]) xmax=x_leaf[order_list[k]];
		}
		x_coordinate[order_list[k]]=0.5*(xmin+xmax);
		/* Horizontal line */

		if(xmax>xmin) horizontal(xmin,xmax,k);
		/* 
				fprintf(pic,"%.4lf m%03d %.4lf m%03d\n",
			xmin,order_list[k],xmax,order_list[k]);
    */
	if(k>0) get_coordinate(k-1);
	else { /* Points */
		if(ps_flag==1) {
			fprintf(pic,"stroke\n");
			fprintf(pic,"%c%c\n%c%c Draw the vertices\n%c%c\n",
				percent, percent, percent, percent, percent, percent);
	  	fprintf(pic,"newpath\n");
		}
		for(j=1;j<=nsites;j++) if(x_coordinate[j]> -100.0) point(j);
		for(j=1;j<=nsites;j++) point_text(j);
	}
			/*
			fprintf(pic,"%.4lf m%03d\n",x_coordinate[j],j);
			*/
}

void scale_array(double *a, double s, int n) {
	int i;
	for(i=0;i<n;i++) a[i] /= s;
}

void get_text() {
	int i;
	/* Text line under leaves */
	for(i=0;i<types;i++) text_line(i);
	for(i=0;i<types;i++) lineage_line(i);
	if(sub_matrix==1) subpopulation_matrix();
	/*
		fprintf(pic,"%4d %.4lf yl\n",
				multiplicity[i],x_leaf[sequence[i][0]]);
	*/
}
		
void equal_age() {
	int i,j,k;
	double m=0.0,x;
	agev=(double*)malloc((nsites+1)*sizeof(double));
	for(i=0;i<=nsites;i++) agev[i]=(-1.0);
	agev[0]=0.0;
	for(i=0;i<types;i++) {
		j=0;
		while(sequence[i][j]!=0) j++;
		for(k=j-1;k>=0;k--) agev[sequence[i][k]]=agev[sequence[i][k+1]]+1.0;
	}
	for(i=0;i<=nsites;i++) if(agev[i]>m) m=agev[i];
	for(i=0;i<=nsites;i++) agev[i]=m-agev[i]+1.0;;
	for(i=0;i<=nsites;i++) agev[i] /=m+1.0;
}

void get_age(char *filename) {			
	int i,j;
	char buffer[100];
	double a,sc;
	if((age=fopen(filename,"r"))==NULL) {
		printf("\n\tCan't open %s\n\n",filename);
		exit(1);
	}
	agev=(double*)malloc((nsites+1)*sizeof(double));
	for(i=0;i<=nsites;i++) agev[i]=(-1.0);
	while(fgets(buffer,100,age)!=NULL) {
		if(sscanf(buffer,"%d %lf",&j,&a)!=2) {
			printf("\n\tProblem reading age file\n");
			exit(1);
		}
	agev[j]=a;
	}
	if(agev[0]==0.0) {
		printf("\n\tProblem with TMRCA\n\n");
		exit(1);
	}
	c_tmrca=agev[0];
	sc=agev[0];
	if(ps_flag==1) sc=agev[0]*xsize/ysize;
	scale_array(agev,sc,nsites+1);
	fclose(age);
}

void free_mem() {
	int i;
	free(order_list);
	free(multiplicity);
	free(margin);
	free(x_coordinate);
	free(x_leaf);
	free(leaf_order);
	free(agev);
	for(i=0;i<types;i++) free(sequence[i]);
	free(sequence);
	for(i=0;i<=nsites;i++) free(descend[i]);
	free(descend);
}

int get_strings(int num,char **in_str, int *in_label,char *filename) {
	int i,test,mlen=0;
	char str[100];
	FILE *in;
	in=fopen(filename,"r");
	if(in==NULL) {
		printf("\n\tCan't open %s to get strings\n");
		exit(1);
	}
	for(i=0;i<num;i++) in_label[i]=(-1);
	for(i=0;i<num;i++) {
		test=fscanf(in,"%d %s",in_label+i,in_str[i]);
		if(test !=2) {
			printf("\n\tProblem reading string file %s\n",filename);
			exit(1);
		}
		if(mlen<strlen(in_str[i])) mlen=strlen(in_str[i]);
	}
	fclose(in);
	return mlen;
}


void example_files() {
	FILE* test;
	FILE* age;
	test=fopen("test.tre","w");
	age=fopen("age.dat","w");
	fprintf(test,"66 :  19  8  1  0\n");
	fprintf(test," 8 :  12 19  8  1  0\n");
	fprintf(test," 1 :   5 19  8  1  0\n");
	fprintf(test," 9 :   1  0\n");
	fprintf(test," 3 :  18  0\n");
	fprintf(test," 1 :  17 18  0\n");
	fprintf(test,"18 :   0\n");
	fprintf(test,"71 :  16  1  0\n");
	fprintf(test," 2 :  11 16  1  0\n");
	fprintf(test," 1 :  15  1  0\n");
	fprintf(test,"24 :  10  6  3 15  1  0\n");
	fprintf(test," 5 :   7 10  6  3 15  1  0\n");
	fprintf(test," 5 :   2  7 10  6  3 15  1  0\n");
	fprintf(test,"15 :  14  7 10  6  3 15  1  0\n");
	fprintf(test," 1 :   9  7 10  6  3 15  1  0\n");
	fprintf(test,"12 :  13  9  7 10  6  3 15  1  0\n");
	fprintf(age,"19 3.54\n");
	fprintf(age,"18 1.75\n");
	fprintf(age,"17 0.268\n");
	fprintf(age,"16 3.72\n");
	fprintf(age,"15 7.09\n");
	fprintf(age,"14 0.933\n");
	fprintf(age,"13 0.588\n");
	fprintf(age,"12 0.579\n");
	fprintf(age,"11 0.303\n");
	fprintf(age,"10 3.78\n");
	fprintf(age,"9 1.22\n");
	fprintf(age,"8 5.06\n");
	fprintf(age,"7 2.46\n");
	fprintf(age,"6 4.49\n");
	fprintf(age,"5 0.143\n");
	fprintf(age,"3 5.18\n");
	fprintf(age,"2 0.685\n");
	fprintf(age,"1 10.14\n");
	fprintf(age,"0 13.84\n");
	fprintf(age,"4  -1.0\n");
	fclose(age);
	fclose(test);
	printf("\tExample files test.tre and age.dat written\n");
}

/* Default is
 * treepic(name,age_name,pic_name,'p',12.0,14.0,19.0,28.0,0.8,2.4,1,1);
 *
 */

void treepic (char *name, char *age_name, char *pic_name, char g,
	double x_size,  double y_size, double x_paper,  double y_paper, 
	double l_width,double ps_r, int ax, int mu, double c_scale) {
	int i,j;
	if(read_tree(name)==0) return;
	sub_matrix=1;
	age_flag=1;
	lineage_label=1;
	bare_flag=0;
	if(g=='P')  {
			 tex_flag=1;ps_flag=0;
	}
	if(g=='p')  {
			 tex_flag=0;ps_flag=1;
	}
	mfontsize=10.0;
	sfontsize=10.0;
	xsize=x_size;
	ysize=y_size;
	xpaper=x_paper;
	ypaper=y_paper;
	lwidth=l_width;
	ps_radius=ps_r;
	axis=ax;
  f_str[0]='\000';
	if(c_scale>0.0)
					/*
		sprintf(f_str,"(%c.2e) rightshow\n",percent);
		*/
		sprintf(f_str,"(%s) rightshow\n","%.0lf");
	else
		sprintf(f_str,"(%s) rightshow\n","%.2lf");
  axis=0;uaxis=1;
	mult=mu;
	if(xsize*types/(types-1.0)>xpaper) xsize=xpaper-5.0;
	if(ysize>ypaper) ysize=xpaper-5.0;
	if((pic=fopen(pic_name,"w"))==NULL) {
		printf("\n\tCan't open %s\n\n",pic_name);
		exit(1);
	}
	get_descend();
	get_age(age_name);
  tmrca=c_tmrca;
	if(c_scale>0.0) tmrca*=c_scale;
  udy=tmrca/4.0;
	preamble();
	leaf_order=(int*)malloc(types*sizeof(int));
	order_site(0);
	leaf_coordinates();
	get_coordinate(order_count-1);
	get_text();
	post();
	fclose(pic);
	/* GCC hates free_mem(), at least on PC Linux, don't bother to free it 
	free_mem();
	*/
}

/* Recognised, not in menu
\n\t-P PicTex file\
\n\t-p Postscript file\
*/

#ifdef TREEPIC
char *usage="\n\ttreepic tree_file out_file [options]";
char *options =
"\tOptions\
\n\t @file_name : options are in a file\
\n\t-a age file\
\n\t-z label lineages\
\n\t-s subpopulation matrix\
\n\t-o [offset] in characters for subpopulation labels\
\n\t-nt No axis in picture\
\n\t-nm No multiplicities in picture\
\n\t-nd No dots in picture\
\n\t-nv No vertex labels in picture\
\n\t-u User axis [tmrca] [tick_increment] [format_string]\
\n\t-v Age of mutations axis [tmrca] [format_string]\
\n\t   format_string uses C syntax, examples: %c.2E, %c.0lf\
\n";

char *display =
"\t+  Display extra options (treepic + )";

char *extra =
"\tExtra Options\
\n\t+Nm filename with mutation labels\
\n\t+Nl filename with lineage labels\
\n\t+Ns filename with subpopulation labels\
\n\t+x xsize (cm)\
\n\t+y ysize\
\n\t+w page width\
\n\t+h page height\
\n\t+l line width (default 0.8 pt)\
\n\t+d dot width (default 4.8 pt)\
\n\t+f fontsize for frequencies, axis (default 10 pt)\
\n\t+m fontsize for mutation labels (default 10 pt)\
\n\t+r landscape\
\n";

char *examples=
"\n\t example: treepic test.tre test.ps -a age.dat\
\n\t treepic -e generates example files test.tre and age.dat\
\n";

#ifdef GENETREE
char *gt =
"\t (genetree generates files tree.gt and tree.age\
\n\t   with -o 4 (or > 4) option)\
\n";
#endif

#ifdef PSVIEW
char *view=
"\t-g view tree\
\n";
#endif

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


char age_name[100];

void get_options(char **argvo, int argco) {
	int i;
	for(i=0;i<argco;i++) if(argvo[i][0]=='-'||argvo[i][0]=='+') {
		switch(argvo[i][1]) {
			case 'a' : age_flag=1; strcpy(age_name,argvo[i+1]); break;
			case 'b' : bare_flag=1;ps_flag=0;age_flag=0;break;
			case 'z' : lineage_label=1;break;
			case 's' : sub_matrix=1;break;
			case 'o' : sub_offset=atoi(argvo[i+1]); break;
			case 'g' : look=1;break;
			case 'P' : tex_flag=1;ps_flag=0;break;
			case 'p' : ps_flag=1;tex_flag=0;break;
			case 'x' : xsize=atof(argvo[i+1]);break;
			case 'y' : ysize=atof(argvo[i+1]);break;
			case 'w' : xpaper=atof(argvo[i+1]);break;
			case 'h' : ypaper=atof(argvo[i+1]);break;
			case 'l' : lwidth=atof(argvo[i+1]);break;
			case 'd' : ps_radius=atof(argvo[i+1]);break;
			case 'f' : fontsize=atof(argvo[i+1]);break;
			case 'm' : mfontsize=atof(argvo[i+1]);break;
			case 'r' : t_landscape=1; break;
			case 'n' : 
								 if(argvo[i][2]=='a') axis=0;
					  		 if(argvo[i][2]=='m') mult=0;
					  		 if(argvo[i][2]=='d') dot=0;
					  		 if(argvo[i][2]=='v') label=0;
								 break;
			case 'e' : example_files();break;
			case 'u' : tmrca=atof(argvo[i+1]);
								 udy=atof(argvo[i+2]);
								 f_str[0]='\000';
								 if(argvo[i+3][0]=='%') 
									sprintf(f_str,"(%s) rightshow\n",argvo[i+3]);
								 else
									sprintf(f_str,"(%c.2e) rightshow\n",percent);
								 axis=0;uaxis=1;
								 break;

			case 'v' : tmrca=atof(argvo[i+1]);
								 f_str[0]='\000';
								 if(argvo[i+2][0]=='%') 
									sprintf(f_str,"(%s) rightshow\n",argvo[i+2]);
								 else
									sprintf(f_str,"(%c.2e) rightshow\n",percent);
								 axis=0;aaxis=1;uaxis=0;
								 break;
			case 'N' : 
								 if(argvo[i][2]=='m') {
									 m_label=1;
									 strcpy(site_fname,argvo[i+1]);
									}
								 if(argvo[i][2]=='l') {
									 l_label=1;
									 strcpy(line_fname,argvo[i+1]);
									}
								 if(argvo[i][2]=='s') {
									 s_label=1;
									 strcpy(subpop_fname,argvo[i+1]);
									}
									break;
			default  : break;
		}
	}
}


main(int argc, char** argv) {
	int i,j;
	char** fargv;
	int fargc;

	if(argc>1) if(argv[1][0]=='+') {
		printf(usage);
		printf(treepic_ver);
		printf(extra);
		exit(1);
	}

	if(argc>1) if(argv[1][0]=='-'&&argv[1][1]=='e') {
			example_files();
			exit(0);
	}

	if(argc<3) {
		printf(usage);
		printf(treepic_ver);
		printf(options,percent,percent,percent);
#ifdef PSVIEW
		printf(view);
#endif
		printf(display);
		printf(examples);
#ifdef GENETREE
		printf(gt);
#endif
		exit(1);
	}
	if(read_tree(argv[1])==0) {
			printf("\n\tProblem reading in tree file\n");
			exit(1);
	}
	if(argv[2][0]=='-') {
		printf("\n\tAn output file is needed as second argument\n");
		exit(1);
	}
	if((pic=fopen(argv[2],"w"))==NULL) {
		printf("\n\tCan't open %s\n\n",argv[2]);
		exit(1);
	}

	fargv=(char**)malloc(FARGS*sizeof(char*));
	for(i=0;i<FARGS;i++)
	fargv[i]=(char*)malloc((FARGS_SIZE+1)*sizeof(char));
  fargc=get_file_arg(NULL,argc,argv,fargv);
	get_options(argv+3,argc-3);
	get_options(fargv,fargc);
	/* XXXXXXXXX */
	if(m_label==1) {
		site_name=(char**)malloc(nsites*sizeof(char*));
		for(i=0;i<nsites;i++)
			site_name[i]=(char*)malloc(MAXNAME*sizeof(char));
		if(label_name!=NULL) free(label_name);
		label_name=(int*)malloc(nsites*sizeof(int));
		get_strings(nsites,site_name,label_name,site_fname);
	}
	if(l_label==1) {
		line_name=(char**)malloc(types*sizeof(char*));
		for(i=0;i<types;i++)
			line_name[i]=(char*)malloc(MAXNAME*sizeof(char));
		if(label_name!=NULL) free(label_name);
			label_name=(int*)malloc(types*sizeof(int));
		get_strings(types,line_name,label_name,line_fname);
	}
	if(s_label==1) {
		subpop_name=(char**)malloc(npop*sizeof(char*));
		for(i=0;i<npop;i++)
			subpop_name[i]=(char*)malloc(MAXNAME*sizeof(char));
		if(label_name!=NULL) free(label_name);
			label_name=(int*)malloc(npop*sizeof(int));
		sub_offset=get_strings(npop,subpop_name,label_name,subpop_fname);
	}
	if(xsize*types/(types-1.0)>xpaper) xsize=xpaper-5.0;
	if(ysize>ypaper) ysize=xpaper-5.0;
	if(tex_flag==1||ps_flag==1) {
		bare_flag=0;
		if(age_flag==0) { equal=1;axis=0;}
	}
	get_descend();
  if(age_flag!=0&&equal==0)	get_age(age_name);
	if(equal==1) equal_age();
	preamble();
	/*
	printf("\n");
	for(i=0;i<types;i++) {
		printf("\n");
	for(j=0;sequence[i][j];j++) printf("%d ",sequence[i][j]);
	}
	printf("\n");
	printf("\n**************************\n");
	for(j=0;j<=nsites;j++) printf("%2d ",j);
	printf("\n");
	for(i=0;i<=nsites;i++) {
		printf("\n");
	 	for(j=0;j<=nsites;j++) {
			if(descend[i][j]>=0) printf("%2d ",descend[i][j]);
			else printf("  .");
		}
	}
	printf("\n");
	printf("\n**************************\n");
	*/
	leaf_order=(int*)malloc(types*sizeof(int));
	order_site(0);
	/*
	printf("\n**************************\n");
	for(j=0;order_list[j]!=(-1);j++) printf("%2d ",order_list[j]);
	printf("\n");
	printf("\n**************************\n");
	*/
	leaf_coordinates();
	get_coordinate(order_count-1);
	get_text();
	post();
	fclose(pic);
#ifdef PSVIEW
  if(look==1) {
		strcat(psview," ");	
		strcat(psview,argv[2]);
		system(psview);
	}
#endif
	return 0;
}
#endif

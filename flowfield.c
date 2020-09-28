/*
 * Flow field calculation
 *
 */


//compile:
//cc flowfield.c -O2 -Wall -o flowfield -lm

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>	//getopt

extern  int      getopt(int, char *const*, const char *);
extern  int      optind;
extern  char    *optarg;


double EPSILON = 0.00001;
#define SQR(a) ((a)*(a))

int N=0;		//# of objects on the field
int nLx=0, nLy=0;	//# of lattice points
double Rx=0, Ry=0;	//size of sampling window on image



void usage (char **argv){
	fprintf(stdout,"\
usage:	%s <in.file >out.file [options]\n\
options:\n\
	-N 67	: number of objects in input file, eg. N=67\n\
	-nx 20	: the number of lattice points of the output (in the x \n\
		  direction, eg. 20).\n\
	-ny 30	: the number of lattice points of the output (in the y \n\
                  direction, eg. 30).\n\
	-rx 300 : sampling window size in the input field, in the direction\n\
		  parallel to the object's velocity vector. Typically \n\
		  measured in pixels (eg. 300)).\n\
	-ry 200 : sampling window size in the input field, in the direction\n\
                  perpendicular to the object's velocity vector. Typically \n\
		  measured in pixels (eg. 200)).\n\
	-h 	: display this help\n\
\n\
Options N, nx, ny, rx, ry can be fed through the stdin before the actual data, if no options are specified, for example:\n\
#67 \n\
#20 30 \n\
#300 200 \n\
<x> <y> <vx> <vy> \n\
 :   :   :    :\n\
\n\
\n\
",argv[0]);
	
	exit (0);
};


void get_args (int num, char** argv){
  int i;
  
  for (i=1;i<num;i++){
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
        
	case 'N':
	case 'n':
		if (argv[i][2] == 'x') 
			nLx=atoi(argv[++i]);	//lattice-geometry
		else if (argv[i][2] == 'y') 
			nLy=atoi(argv[++i]);
		else N=atoi(argv[++i]);		//# of objects
		break;
	case 'R':				//image geometry
	case 'r':
		if (argv[i][2] == 'x' || argv[i][2] == 'X') Rx=atof(argv[++i]);
		if (argv[i][2] == 'y' || argv[i][2] == 'Y') Ry=atof(argv[++i]);
                break;
	case 'h': 				//display help
		usage(argv);
		break;
        default:  ; 
      }
    }; 
  };
  
};

void check_args(char **argv){
  if (N==0) {
  	fprintf(stderr,"%s: invalid or no N parameter\n",argv[0]); 
	exit(0);};
  if (nLx==0) {
  	fprintf(stderr,"%s: invalid or no nx parameter\n",argv[0]); 
	exit(0);};
  if (nLy==0) {
  	fprintf(stderr,"%s: invalid or no ny parameter\n",argv[0]); 
	exit(0);};
  if (Rx==0) {
  	fprintf(stderr,"%s: invalid or no Rx parameter\n",argv[0]); 
	exit(0);};
  if (Ry==0) {
  	fprintf(stderr,"%s: invalid or no Ry parameter\n",argv[0]); 
	exit(0);};

}


int main (int argc, char** argv){
 int i,j,ix,iy;
 double *x,*y,*vx,*vy,dx,dy;
 double hRx=0,hRy=0,vopa[2],vope[2],a;
 
 get_args (argc, argv);

 if (argc <=2) {			//read parameters from input
	scanf("#%i\n",&N);
	scanf ("#%i %i\n",&nLx,&nLy);
	scanf ("#%lf %lf\n",&Rx,&Ry);
 };

 check_args(argv);
 
 hRx=0.5*Rx;
 hRy=0.5*Ry;

 x=(double*)calloc(N,sizeof(double));
 y=(double*)calloc(N,sizeof(double));
 vx=(double*)calloc(N,sizeof(double));
 vy=(double*)calloc(N,sizeof(double));
 

//read vectorfield:
 for (i=0;i<N;i++){
   scanf("%lf %lf %lf %lf\n",&x[i],&y[i],&vx[i],&vy[i]);
 };


//process:
 for (i=0;i<N;i++){
   a=sqrt(SQR(vx[i])+SQR(vy[i]));
   if ( a>EPSILON ) {	//define new coordinate system, attached to 'i'
     vopa[0]=vx[i]/a;	//unit vector of cell velocity
     vopa[1]=vy[i]/a;
     vope[0]=-vy[i]/a;	//unit vector perpendicular to velocity
     vope[1]=vx[i]/a;


//actual calculation:
     for (j=0;j<N;j++){
       if (j!=i){
         dx=(x[j]-x[i])*vopa[0]+(y[j]-y[i])*vopa[1];
         dy=(x[j]-x[i])*vope[0]+(y[j]-y[i])*vope[1];
         if ( (fabs(dx) < hRx) && (fabs(dy) < hRy) ){
//process only if j is within range hR of i
//the coordinates of point 'j' in the new system:
	 	ix=(int)((double)(nLx)*((double)(dx)+hRx)/(Rx));
		iy=(int)((double)(nLy)*((double)(dy)+hRy)/(Ry));
                      
// binned coordinate values (center of bin):
		dx=(double)(hRx*(2*ix+1))/((double)(nLx))-hRx;
		dy=(double)(hRy*(2*iy+1))/((double)(nLy))-hRy;

		
		fprintf(stdout, "%lf %lf %lf %lf %lf %lf\n", dx, dy,
					(vx[j] * vopa[0] + vy[j] * vopa[1]),
					(vx[j] * vope[0] + vy[j] * vope[1]),
					x[i], y[i]);
	 };
       };
     };
   };
 };



 free (x);
 free (y);
 free (vx);
 free (vy);

 return (0);
}

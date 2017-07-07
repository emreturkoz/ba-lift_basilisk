#include "axi.h"
#include "navier-stokes/centered.h"
#include "fractions.h"
#include "vof.h"
#include "tension.h"

scalar f12[], f23[];
scalar *interfaces = {f12, f23};

#define TAUABS 23.49E-9
#define H0ABS 8.27E-6
#define R0ABS 22.365E-6
#define VELBLISTER (H0ABS/TAUABS)

#define RHOLABS 1030.
#define MULABS 1.7E-3
#define MUSABS 1*MULABS
#define RHOAABS 1.1614
#define MUAABS 1.85E-5
#define RHOSABS 1*RHOLABS

#define R0 0.50
#define LSCALE (R0ABS/R0)
#define H0 (H0ABS/LSCALE)

#define SURFABS 0.04079 // Surface tension (N/m)
#define SURF 1.0

#define ReEXP (RHOLABS*VELBLISTER*R0ABS/MULABS)
#define WeEXP (RHOLABS*VELBLISTER*VELBLISTER*R0ABS/SURFABS)

// Dimensionless parameters
#define RHOL 1030
#define VELSCALE (sqrt(SURF*WeEXP/(RHOL*R0)))
#define TAU (H0/VELSCALE)
#define TSCALE (TAUABS/TAU)

#define MUL (RHOL*VELSCALE*R0/ReEXP)
#define MUA (MUL*MUAABS/MULABS)
#define RHOA (RHOL*RHOAABS/RHOLABS)
#define MUS (MUL*MUSABS/MULABS)
#define RHOS (RHOSABS*RHOL/RHOLABS)

  // Manually entered thickness values of the solid and fluid layers
#define HFABS 5.0E-6
#define HSABS 3.0E-6
double LS = HSABS/LSCALE;
double LF = HFABS/LSCALE; 
double CBLISTER = 1.25;
//#define CALCMU(f12,f23) ( (MUA*(1-f23)*(1-f12)) + (MUL*f23) + (MUS*f12)   )
//#define CALCRHO(f12,f23) ( (RHOA*(1-f23)*(1-f12)) + (RHOL*f23) + (RHOS*f12)  )
#define CALCRHO(f23) ( (RHOA*(1-f23)) + (RHOL*f23)  )
#define CALCMU(f23) (  (MUA*(1-f23) ) + (MUL*f23) )





int MAXLEVEL = 6;


vertex scalar psi12[];
vertex scalar psi23[];

face vector visc[];
face vector alphav[];



u.n[left] = neumann (0.0);
u.n[right] = dirichlet(0.0);
u.n[top] = neumann(0.0);



int main(){
  L0 = 2.0; // domain size
  DT = 0.025*TAU;
  origin(0.0, 0.0);
  mu = visc;
  alpha = alphav;
  //f12.sigma = 1000;
  f23.sigma = SURF;
  init_grid(1<<MAXLEVEL);
  run();
}

event init (t=0){
  
  mask( y > 1.0 ? top: none);

  
  foreach(){
    u.x[] = 0.0;
  }

  refine ( x < 0.4  && level < MAXLEVEL + 2);
  //refine ( y < R0/2.5  && level < MAXLEVEL + 1);
  
  // initialize the second fluid layer

  /*
  foreach_vertex(){
    if ( (x>LS) && x<(LS+LF)  ) 
      psi23[] = -x + LS + LF;
    else
      psi23[] = -x;
  }

  fractions(psi23, f23);
  */

  fraction(f23, (+LS + LF- x));
  
  fraction(f12, (+LS  - x)); 

  foreach(){
    psi12[] =  (+LS  - x);
  }

  static FILE *fp = fopen("configuration.txt", "w");
  fprintf(fp, "%g %g %g %g %g \n", TAU, TSCALE, LSCALE,  H0, LS+H0*(2.0/M_PI)*atan(2));
  fclose(fp);
    
}
 

event moving_blister (i++) {
  scalar f12[];

  /*
  if(t<=2*TAU)
     fraction (f12, ((LS  + (H0*pow( 1.0 - ((y)*(y)/(R0*R0)),CBLISTER )*(2.0/M_PI)*atan(t/TAU)))) ); 
  if(t>2*TAU)
    fraction (f12, ((LS  + (H0*pow( 1.0 - ((y)*(y)/(R0*R0)),CBLISTER )*(2.0/M_PI)*atan(2)))) );
  */

  /*
  if(t<=2.0*TAU){
  foreach(){
    foreach_dimension(){
      u.x[] = f12[]*H0*pow( 1- (y*y)/(R0*R0), CBLISTER)*(2.0/(M_PI*TAU*(1+ (t*t/TAU/TAU))) ) + (1-f12[])*u.x[];
    }
  }
  }
  */

  
  foreach(){
      // Velocity BC on the blister. 
    if( (t<=2.0*TAU) && (y<=R0) && ( x <= LS+ (H0*pow( 1.0 - ((y)*(y)/(R0*R0)),CBLISTER )*(2.0/M_PI)*atan(t/TAU)) )       ){      
      u.y[] = 0.0;
      u.x[] = H0*pow( 1- (y*y)/(R0*R0), CBLISTER)*( 2.0/( M_PI*TAU ) )*( 1 / ( 1+ (t*t/TAU/TAU) ) );
      }

    
    // No slip for the BC at the interface beyond the blister
    if( (y>R0) && (x<LS)  ){
      u.x[] = 0.0;
      u.y[] = 0.0;
    }
    

    // After the blister expansion stops - > no slip bc
    
    if( (t > 2.0*TAU) && (y<=R0) &&  (x < LS+ (H0*pow( 1.0 - ((y)*(y)/(R0*R0)),CBLISTER )*(2.0/M_PI)*atan(2)))  ){
      u.x[] = 0.0;
      u.y[] = 0.0;
      }
    

    
    if( (t>2.0*TAU) && (x<=LS) ){
      u.x[] = 0.0;
      u.y[] = 0.0;
    }
    
    
    
  }

  boundary ((scalar *){u});
  


// material properties

  foreach_face(){
    //double T12 = (f12[] + f12[-1,0])/2.;
    double T23 = (f23[] + f23[-1,0])/2.;

    visc.x[] = CALCMU(T23);
    alphav.x[] = 1./CALCRHO(T23);
  }


}


/*event properties(i++){

  
  if((t<=2*TAU)){
    foreach(){
      if(y<=R0)
      psi12[] = -LS - (H0*pow( 1.0 - ((y)*(y)/(R0*R0)),1.0 )*(2.0/M_PI)*atan(t/TAU)   )  + x;

    }
    
    fractions(psi12,f12);

      //fraction(f12, (-LS - (H0*pow( 1.0 - ((y)*(y)/(R0*R0)),CBLISTER )*(2.0/M_PI)*atan(t/TAU)   )  + x)); 

    foreach_vertex(){
      f12[] =  (1- f12[]);
    }
  } 
    
 
  if ( (t==2*TAU) || (t==TAU) || (t==0.5*TAU)  ){
    char name[99];  
    sprintf(name, "interface_%g.dat", t);
    FILE * fp = fopen(name, "w");
    output_facets(f12,fp);
    fclose(fp);
  }
  
  
*/
  
/*  foreach(){
    // Velocity BC on the blister. 
    if( (t<=2.0*TAU) && (y<=R0) && (f12[]>0.0) ){      
      u.y[] = 0.0;
      u.x[] = H0*pow( 1- (y*y)/(R0*R0), CBLISTER)*(2.0/(M_PI*TAU*(1+ (t*t/TAU/TAU))) );
      
    }

    // No slip for the BC at the interface beyond the blister
    if( (y>R0) && (f12[]>0.0)  ){
      u.x[] = 0.0;
      u.y[] = 0.0;
    }

    // After the blister expansion stops - > no slip bc
    if( (t > 2.0*TAU) && (f12[]>0.0) ){
      u.x[] = 0.0;
      u.y[] = 0.0;
    }
    
    
    
  }
  
  
}
*/


/*
event outputgfs (t += 0.2*TAU; t<=2.0*TAU){
  static int nf = 0;
  char name[80];
  sprintf(name, "output%d.gfs", nf);
  output_gfs (file = name);
  nf++;
}
*/

event gfsview (i += 10; t <= 50.0) {
  static FILE * fp = popen ("gfsview2D f23velocity_iso.gfv", "w");
  scalar omega[];
  vorticity (u, omega);  
  output_gfs (fp, t = t);
}


event images (t += 0.1*TAU; t<=10.0*TAU){
  static FILE * fp = fopen("f12.ppm", "w");
  output_ppm(f12, fp, min=0, max=1, linear=true);

  static FILE * fp2 = fopen("f23.ppm", "w");
  output_ppm(f23, fp2, min=0, max=1, linear=true);
  
}

// output of the blister profile
//event output (t+=0.001; t<=0.2){
event output (t+=0.025*TAU; t<=5*TAU){
  static int nf= 0;
  char name[100];
  sprintf(name,"blister_%g.dat",t/TAU);
  FILE *fp = fopen(name,"w");
  output_facets(f12,fp);
  fclose(fp);
  nf++;
  }

event output (t+=0.001; t<=0.2){
  static int nf= 0;
  char name[100];
  sprintf(name,"jet_%d.dat",nf);
  FILE *fp = fopen(name,"w");
  output_facets(f23,fp);
  fclose(fp);
  nf++;
  }
/*
event interface( t += TAU ) {
  char name[100];
  static int nf = 0;
  sprintf(name, "interface_%d.dat", nf);
  FILE * fp = fopen(name, "w");
  output_facets(f23,fp);
  fclose(fp);
  nf++;
}
*/


/*
event extractPosition (i += 10; t*TSCALE<100e-9){
  static FILE *fp = fopen("jetLength.txt", "a");

  vector h[];
  heights (f23,h);
  double yMax = -HUGE;;

  foreach(){
    if (h.x[] != nodata){
      double xi = x + height(h.x[])*Delta;
      if (xi > yMax)
	yMax = xi;
    }
  }
  fprintf(fp, "%g %g \n", t*TSCALE, LSCALE*(yMax-(LF+LS)));
  printf("%g %g \n", t*TSCALE, LSCALE*(yMax-(LF+LS)));
}
*/

//#include "grid/multigrid3D.h"
#include "navier-stokes/centered.h"
#include "fractions.h"
#include "vof.h"
#include "tension.h"


scalar f12[], f23[];
scalar psi23[];
scalar *interfaces = {f12, f23};

#define TAUABS 23.49E-9
#define H0ABS 4.5963E-6
#define R0ABS 17.9979E-6

#define RHOLABS 1030.
#define MULABS 1.7E-3
#define MUSABS (MULABS)
#define RHOAABS 1.1614
#define MUAABS 1.85E-5
#define RHOSABS (RHOLABS)

#define SURFABS 0.04079 // ABSOLUTE VALUE OF THE SURFACE TENSION (N/m)

// Calculating the inertia-capillary time scale and its ratio to the blister expansion time
#define TCABS (sqrt(RHOLABS*R0ABS*R0ABS*R0ABS/SURFABS))
#define TRATIOABS (TCABS/TAUABS)

// Experimental Ohnesorge number assuming that the forcing does not matter
#define OHEXP (MULABS/sqrt(RHOLABS*SURFABS*R0ABS))

#define R0 0.75
#define LSCALE (R0ABS/R0)
#define H0 (H0ABS/LSCALE)
#define TSCALE (LSCALE)

#define TAU (TAUABS/TSCALE)
#define RHOSIGMARATIO   (  (TAU*TRATIOABS/sqrt(R0*R0*R0)) * (TAU*TRATIOABS/sqrt(R0*R0*R0))   )

#define SURF 1
#define RHOL (SURF*RHOSIGMARATIO)


#define MUL ( OHEXP*sqrt(RHOL*SURF*R0) )


// PARAMETERS OF OTHER PHASES
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
//#define CALCMU(f12) ( MUA*(1-f12) + (MUL*f12)   )
#define CALCRHO(f12,f23) ( (RHOA*(1-f23)*(1-f12)) + (RHOL*f23) + (RHOS*f12)  )
//#define CALCRHO(f12) ( (RHOA*(1-f12)) + (RHOL*f12) )
#define CALCMU(f23) (  (MUA*(1-f23) ) + (MUL*f23) )





int MAXLEVEL = 7;


vertex scalar psi12[];

face vector visc[];
face vector alphav[];



u.n[left] = neumann(0.0);
//u.t[left] = dirichlet(0.0);
u.n[right] = neumann(0.0);
//u.n[right] = dirichlet(0.0);
u.n[top] = neumann(0.0);


u.n[bottom] = dirichlet(0.0);


//p[top] = dirichlet(0.0);
//p[left] = dirichlet(0.0);
//p[right] = neumann(0.0);


int main(){
  L0 = 4.0; // domain size
  //DT = 1e-9;
  DT = 0.05*TAU;
  //origin(0.0, 0.0);
  mu = visc;
  alpha = alphav;
  f23.sigma = SURF;
  init_grid(1<<MAXLEVEL);
  run();
}

event init (t=0){
  
  mask( y > 1.0 ? top: none);

  
  foreach(){
    u.x[] = 0.0;
  }
  

  //refine ( x < 0.7  && level < MAXLEVEL + 1);

 // refine ( x < 0.5  && level < MAXLEVEL + 1);

  fraction(f23, (+LS + LF- x));
 
  static FILE *fp = fopen("configuration.txt", "w");
  fprintf(fp, "%g %g %g %g %g \n", TAU, TSCALE, LSCALE,  H0, LS+H0*(2.0/M_PI)*atan(2));
  fclose(fp);
    
}
 

event moving_blister (i++) {

  // Need to update the interface as well


  vertex scalar f12_[];
  foreach_vertex(){
	if ( sqrt(y*y + z*z) <=R0 )
		f12_[] = (+LS + H0*pow( 1- ((y*y) + (z*z))/(R0*R0), CBLISTER)*( 2.0/( M_PI ) )*atan(t/TAU)  - x);
	else
		f12_[] = (+LS - x);
  }
  fractions(f12_, f12);
  

  
 
  
  foreach(){
      // Velocity BC on the blister. 
    if( (t<=2.0*TAU) && (sqrt(y*y + z*z) <= R0) &&  ( x <= LS+ (H0*pow( 1.0 - ((y*y + z*z)/(R0*R0)),CBLISTER )*(2.0/M_PI)*atan(t/TAU)) )    ){      
      u.y[] = 0.0;
      u.x[] = f12[]*H0*pow( 1- ((y*y + z*z)/(R0*R0)), CBLISTER)*( 2.0/( M_PI*TAU ) )*( 1 / ( 1+ (t*t/TAU/TAU) ) )  + (1. - f12[])*u.x[];
      u.z[] = 0.0;
      }

    
    // No slip for the BC at the interface beyond the blister
    if( sqrt(y*y + z*z) > R0 && (x<LS)  ){
      u.x[] = 0.0;
      u.y[] = 0.0;
      u.z[] = 0.0;
    }
    

    // After the blister expansion stops - > no slip bc
    
    if( (t > 2.0*TAU) && (sqrt(y*y + z*z) <= R0) && (x < LS+ (f12[])*(H0*pow( 1.0 - (( y*y + z*z )/(R0*R0)),CBLISTER )*(2.0/M_PI)*atan(t/TAU))) ){
      u.x[] = 0.0;
      u.y[] = 0.0;
      u.z[] = 0.0;
    }
   
    
  }

 
    
  boundary ((scalar *){u});


// material properties

  foreach_face(){
    double T12 = (f12[] + f12[-1,0])/2.;
    double T23 = (f23[] + f23[-1,0])/2.;
 
    visc.x[] = CALCMU(T23);
    alphav.x[] = 1./CALCRHO(T12,T23);
  }


}


event gfsview (i += 5; t <= 2000*TAU) {
  static FILE * fp = popen ("gfsview3D velocity3d.gfv", "w");
  scalar omega[];
  vorticity (u, omega);  
  

  scalar velocity[];
  foreach()
    velocity[] = (1-f12[])*sqrt(u.x[]*u.x[] + u.y[]*u.y[] + u.z[]*u.z[]);
  output_gfs (fp, t = t);



  

}


event images (t += 0.1*TAU; t<=5.0*TAU){
  static FILE * fp = fopen("f12.ppm", "w");
  output_ppm(f12, fp, 512, min=0, max=1, linear=true);

   
}


/*
event output (t += 0.1*TAU; t<=2.1*TAU){
  static int nf= 0;
  char name[100];
  sprintf(name,"blister_%g.dat",t/TAU);
  FILE *fp = fopen(name,"w");
  output_facets(f12,fp);
  fclose(fp);
  nf++;

  char name2[100];
  sprintf(name2,"pressure_%g.dat",t/TAU);
  FILE *fp2 = fopen(name2,"w");
  output_ppm(p,fp2, n=512);

  scalar velocity[];
  foreach()
    velocity[] = sqrt(u.x[]*u.x[] + u.y[]*u.y[]);
  boundary ({velocity});

  char name3[100];
  sprintf(name3,"velocity_%g.dat",t/TAU);
  FILE *fp3 = fopen(name3,"w");
  output_ppm(velocity,fp3, n=512);



  static FILE *fp4 = popen ("ppm2mpeg > velocity_all.mpg", "w");
  output_ppm(velocity, fp4, n=512);



  char name5[100];
  sprintf(name5,"f23_%g.dat",t/TAU);
  FILE *fp5 = fopen(name5,"w");
  output_ppm(f23,fp5, n=512);

  }


*/

/*
event output (t += 0.1*TAU; t<=2.1*TAU){
  static int nf= 0;
  char name[100];
  sprintf(name,"jet_%g.dat",t/TAU);
  FILE *fp = fopen(name,"w");
  output_facets(f23,fp);
  fclose(fp);
  nf++;
  }
*/



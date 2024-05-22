#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "mt19937ar.c"
unsigned long seed = 0;     //Seed for random number generator

//Maximum number of neighbors per particle
#define MAXNEIGH 35

#include "harddisks.h"
#include "mystring.c"

//Number of extra events (e.g. write, thermostat) to allocate space for
#define EXTRAEVENTS 12

//Pi (if not already defined)
#ifndef M_PI
#define M_PI 3.1415926535897932
#endif


double maxtime = 20000;           //Simulation stops at this time
double writeinterval = 1;     //Time between output to screen / data file

int makesnapshots = 1;          //Whether to make snapshots during the run (yes = 1, no = 0)
double snapshotinterval = 50;
int dumplogtime = 0 ;        // 1 if you want to dump configurations with cycles of exponentially spaced time steps
double dumplogbase = 1 ;
int dump_logcyclelength = 1 ;

int initialconfig = 2;    //= 0 load from file, 1 = SQUARE crystal, 2 = HEXAGONAL crystal, 3 = RANDOM, 4 = STAMPFLI (random wheels), 5 = ADDSMALL (random QC12 large in input file)
char inputfilename[100] = "init.sph"; //File to read as input snapshot (for initialconfig = 0)
int N = 5000;             //Number of particles (for FCC)
double areafrac = 0.83;               // nominal area fraction of the sedimented system (for bi-disperse mixture)
// To generate bi-disperse mixtures of hard disks having the same mass-density
double sizeratio = 0.54;              // ratio among the diameters of small and large particles
double large2totalfraction = 1.0;    // fraction of large particles, default = 1
// Non-additivity interactions can be considered, among particles with different radii
double nonadditivity = 1.0;

int initialvelocities = 0;   //=0 generate velocities according a Maxwell-Boltzmann distribution, 1 = loads them from file
char inputvelfile[100] = "init-vel.xyz";   //file to read input velocities (for initvelocities = 1)

//Variables related to the event queueing system. These can affect efficiency.
//The system schedules only events in the current block of time with length "eventlisttime" into a sorted binary search tree. 
//The rest are scheduled in unordered linked lists associated with the "numeventlists" next blocks.
//"numeventlists" is roughly equal to maxscheduletime / eventlisttime
//Any events occurring even later are put into an overflow list
//After every time block with length "eventlisttime", the set of events in the next linear list is moved into the binary search tree.
//All events in the overflow list are also rescheduled.

//After every "writeinterval", the code will output two listsizes to screen. 
//The first is the average number of events in the first list that gets moved into the event tree after each block.
//The second is the length of the overflow list at the last time it was looped over (0 if not processed in the last "writeinterval").
//Ideally, we set maxscheduletime large enough that the average overflow list size is negligible (i.e. <10 events)
//Also, there is some optimum value for the number of events per block (scales approximately linearly with "eventlisttime").
//I seem to get good results with an eventlisttime chosen such that there are a few hundred events per block, and dependence is pretty weak (similar performance in the range of e.g. 5 to 500 events per block...)
double maxscheduletime = 1.0;
int numeventlists;
double eventlisttimemultiplier = 0.5;  //event list time will be this / N
double eventlisttime;


//Neighbor lists
double shellsize = 1.5; //Shell size (equals 1+ \alpha)


//Internal variables
double simtime = 0 , simtimewindowlength = 1000 ;  //to avoid numerical precision on predicted event times from degrading with increasing simulation time
int timewindow = 0 ;                               // simtime is reset every simtimewindowlength, with timewindow taking into account how many times it happens
double reftime = 0 ;
int reftimewindow = 0 ;
int currentlist = 0;
int totalevents;

int listcounter1 = 0, listcounter2 = 0, mergecounter = 0;

particle** eventlists; //Last one is overflow list

particle* particles;
particle** celllist;
particle* root;
double xsize, ysize; //Box size
double hx, hy; //Half box size
double icxsize, icysize; //Inverse Cell size
int    cx, cy;  //Number of cells
double dptot[3] ;   //Momentum transfer projections (for calculating pressure tensor)
int XX = 0, YY = 1, XY = 2 ;
unsigned int colcounter = 0; //Collision counter (will probably overflow in a long run...)

int usethermostat = 0; //Whether to use a thermostat
double thermostatinterval = 0.01;


int main( int argc, char **argv )
{
    if (argc>1) setparametersfromfile(argv[1]) ;
    else setparametersfromfile("") ;
    init();
    if (argc>1) setattributesfromfile(argv[1]) ;
    else setattributesfromfile("input.dat") ;
    printf("Starting\n");

    checkoverlaps();
    int starting_time = time(0);
    while (simtime + simtimewindowlength*timewindow <= maxtime)
    {
      step();
    }
    timewindow = (int) (maxtime / simtimewindowlength) ;
    simtime = maxtime - simtimewindowlength * timewindow ;
    int elapsed_time = time(0) - starting_time ;
    int sec = elapsed_time % 60 ;
    int min = ( (elapsed_time-sec) / 60 ) % 60 ;
    int hour = (elapsed_time-sec-60*min) / 3600 ;
    printf( "\n   Elapsed time: %d:%d:%d\n", hour , min , sec ) ;

    printstuff();
    outputsnapshot();

    free(celllist);
    free(particles);
    free(eventlists);
    return 0;
}

/**************************************************
**                 COMPUTEKINENERGY
** Computes kinetic energy at simulation time
**************************************************/
void computeKinenergy(double *kinEn)
{
    int i;
    particle *p, up2datep;
    double v2tot = 0;
    *kinEn = 0.0 ;

    for (i = 0; i < N; i++)
    {
        p = particles + i;
	updatedparticle(p, &up2datep);
        v2tot += p->mass * (up2datep.vx * up2datep.vx + up2datep.vy * up2datep.vy);
    }

    *kinEn = 0.5 * v2tot ;
}



/**************************************************
**                 PRINTSTUFF
** Some data at the end of the simulation
**************************************************/
void printstuff()
{
    int i;
    particle* p;
    double vfilled = 0;
    double kinEn = 0.0 ;

    for (i = 0; i < N; i++)
    {
        p = particles + i;
        vfilled += p->radius * p->radius;
    }
    vfilled *= M_PI;
    computeKinenergy(&kinEn) ;
    printf("Average kinetic energy: %g\n", kinEn / N);
    double area = xsize * ysize ;
    double dens = N / area;
    double time = simtime + simtimewindowlength*timewindow ;
    double press = -(dptot[XX] + dptot[YY]) / (2.0 * area * time);
    double pressid = dens;
    double presstot = press + pressid;
    printf("Total time simulated  : %g\n", time);
    printf ("Density               : %g\n", dens);
    printf("Packing fraction      : %g\n", vfilled / area);
    printf("Measured pressure     : %g + %g = %g\n", press, pressid, presstot);

}


/**************************************************
**                    INIT
 **************************************************/
void init()
{
    int i;
    if( seed == 0 ) {
      FILE *fp=fopen("/dev/urandom","r");
      int tmp = fread(&seed,1,sizeof(unsigned long),fp);
      if (tmp != sizeof(unsigned long)) printf ("error with seed\n");
      fclose(fp);
    }
    printf("Seed: %u\n", (int)seed);
    init_genrand(seed);

    if (nonadditivity == -1) nonadditivity = sqrt(sizeratio) / (1 + sizeratio) * 2 ;

    if (initialconfig == 0)
    {
      loadparticles();
    } else if (initialconfig == 1)           squarelattice();
    else if (initialconfig == 2)             hexagonal();
    else if (initialconfig == 3)             randomconfiguration();
    else if (initialconfig == 4)             randomStampfli();
    else if (initialconfig == 5) {
      loadparticles();
      addsmallparticles() ;
      outputsnapshot();
      exit(0) ;
    }
    if (initialvelocities == 0) randommovement();
    else loadvelocities();
    hx = 0.5 * xsize ; hy = 0.5 * ysize ;	//Values used for periodic boundary conditions

    initevents();
    for (i = 0; i < N; i++)
    {
        particle* p = particles + i;
        p->boxestraveledx = 0;
        p->boxestraveledy = 0;
        p->nneigh = 0;
        p->counter = 0;
        p->t = simtime ;
        p->timewindow = timewindow ;
        p->xn = p->x;   //Set center of neighbor list to current position
        p->yn = p->y;
    }
    initcelllist();

    for (i = 0; i < N; i++)
    {
        makeneighborlist(particles + i);
    }
    printf("Done adding collisions\n");

    // initializing the pressure tensor
    dptot[XX] = 0.0 ;
    dptot[YY] = 0.0 ;
    dptot[XY] = 0.0 ;

}

/******************************************************
**               INITPARTICLES
** Allocate memory for particle array
** Note that this includes the extra events
******************************************************/
void initparticles(int n)
{
    N = n;
    totalevents = N + EXTRAEVENTS;
    particles = (particle*)calloc(totalevents, sizeof(particle));
    if (!particles)
    {
        printf("Failed to allocate memory for particles\n");
        exit(3);
    }
}


/******************************************************
**               MYGETLINE
** Reads a single line, skipping over lines
** commented out with #
******************************************************/
int mygetline(char* str, FILE* f)
{
    int comment = 1;
    while (comment)
    {
        if (!fgets(str, 255, f)) return -1;
        if (str[0] != '#') comment = 0;
    }
    return 0;
}


/**************************************************
**                    SQUARELATTICE
** Initialize system on a square lattice
**************************************************/
void squarelattice()
{
    int i, j;
    particle* p;

    int ncell = sqrt(N) + 0.0001;

    if (ncell * ncell != N)
    {
        printf("N should be a perfect square! (e.g. %d)\n", ncell * ncell);
        exit(3);
    } else if (areafrac > M_PI/4.)
    {
        printf("The area packing fraction of a square lattice of Hard Disks cannot be greater than ~0.785!\n");
        exit(3);
    }

    double hardcoreradius = 0.5 ;
    double area = N * M_PI * hardcoreradius * hardcoreradius / areafrac ;
    xsize = sqrt(area);
    ysize = xsize;

    double step = xsize / ncell;

    printf("step: %g\n", step);
    initparticles(N);
    printf("Placing particles\n");

    p = particles;
    for (i = 0; i < ncell; i++) for (j = 0; j < ncell; j++)
    {
        p->x = i * step;
        p->y = j * step;
        p++;
    }

    for (i = 0; i < N; i++)
    {
        p = particles + i;
        p->radius = hardcoreradius ;
        p->type = 0;
        p->mass = 1;
        backinbox(p);
    }

    printf("Area packing fraction: %g\n", M_PI * N * hardcoreradius * hardcoreradius / (xsize * ysize) ) ;
    printf("Starting configuration from square lattice crystal\n");
}




/**************************************************
**                    LOADPARTICLES
** Read particles from file
** First line: Number of particles
** Second line: box size along all three axes
** Rest: particle data
** Each particle line consists of :
** - A character indicating type (a = 0, b = 1, etc.)
** - 3 coordinates (x, y, z)
** - and the hard sphere radius
**************************************************/
void loadparticles()
{
    char tmp;
    int i, npart;
    particle* p;
    char buffer[255];
    double dummy = 0.0 ;
    double t0 = 0 ;

    FILE* file;
    file = fopen(inputfilename, "r");
    if (!file)
    {
        printf("File not found!\n");
        exit(3);
    }
    mygetline(buffer, file);
    int ftmp = sscanf(buffer, "%d %lf", &npart, &t0);
    if (ftmp != 2) ftmp = sscanf(buffer, "%d", &npart);
    if (ftmp != 1 && ftmp != 2) { printf("Read error (n or box)\n"); exit(3); }
    simtime = t0 ;
    timewindow = (int)(simtime / simtimewindowlength) ;
    simtime -= simtimewindowlength * timewindow ;
    mygetline(buffer, file);
    ftmp = sscanf(buffer, "%lf %lf\n", &xsize, &ysize);
    if (ftmp != 2) { printf("Read error (n or box)\n"); exit(3); }


    initparticles(npart);
    printf("Placing particles\n");
    double vfilled = 0;
    for (i = 0; i < N; i++)
    {
        p = &(particles[i]);
        mygetline(buffer, file);
        ftmp = sscanf(buffer, "%c %lf  %lf  %lf %lf\n", &tmp, &(p->x), &(p->y), &(dummy), &(p->radius));
        backinbox(p);
        if (ftmp != 5) { printf("Read error (particle) %d \n String: %s\n", ftmp, buffer); exit(3); }
        p->type = tmp - 'a';
        p->mass = 1;
        vfilled += p->radius * p->radius ;
    }
    fclose(file);

    printf("Area packing fraction: %lf\n", M_PI / (xsize * ysize) * vfilled);
    printf("Starting configuration read from %s\n", inputfilename);
}


/**************************************************
**                    HEXAGONAL
** Initialize system on an hexagonal lattice
**************************************************/
void hexagonal()
{
    int i, j;
    particle* p;

    double hardcoreradius = 0.5 ;
    double area = N * M_PI * hardcoreradius * hardcoreradius / areafrac ;
    double step = sqrt( area * 2.0 / sqrt(3) / N ) ;
    int ncelly = 2 * round( sqrt( N / sqrt(3) ) ) , ncellx = 0 ;
    while( N % ncelly != 0 && ncelly > 2 ) ncelly -= 2 ;
    ncellx = N / ncelly ;
    if( ncelly < 0.8 * 2 * ncellx / sqrt(3) || N % ncellx != 0 ) {
      printf( "*** ERROR: Please, choose a different particles number, like %d\n", (int)round(N*sqrt(3)/2.) ) ;
      exit(3) ;
    } else  printf( "\n  Generation of an hexagonal 2D lattice; cells per side : %d,%d ; lattice constant : %g\n", ncellx , ncelly , step ) ;

    if (areafrac > M_PI/2./sqrt(3))
    {
        printf("The area packing fraction of an hexagonal lattice of Hard Disks cannot be greater than ~0.9068!\n");
        exit(3);
    }

    xsize = step * ncellx ;
    ysize = step * ncelly * 0.5*sqrt(3.0) ;

    initparticles(N);
    printf("Placing particles\n");

    p = particles;
    for (i = 0; i < ncellx; i++) for (j = 0; j < ncelly; j++)
    {
        p->x = ( i + j * 0.5 ) * step ;
        p->y = j * 0.5*sqrt(3) * step ;
        p++;
    }

    for (i = 0; i < N; i++)
    {
        p = particles + i;
        p->radius = hardcoreradius ;
        p->type = 0;
        p->mass = 1;
        backinbox(p);
    }

    printf("Area packing fraction: %g\n", M_PI * N * hardcoreradius * hardcoreradius / (xsize * ysize) ) ;
    printf("Starting configuration from an hexagonal lattice crystal\n");
}



/**************************************************
**             RANDOMCONFIGURATION
** Randomly spreads out HD particles
**************************************************/
void randomconfiguration()
{
    particle* p ;
    double hardcoreradius = 0.5 ;
    int i, Nlarge = N * large2totalfraction ;
    large2totalfraction = (double)Nlarge / N ;

    initparticles(N);
    for (i = 0; i < Nlarge; i++) {
        p = particles + i ;
        p->radius = hardcoreradius ;
        p->type = 0 ;
        p->mass = 1.0 ;
    }
    for (i = Nlarge; i < N; i++) {
        p = particles + i ;
        p->radius = hardcoreradius * sizeratio ;
        p->type = 1 ;
        p->mass = sizeratio * sizeratio ;
    }

    //given the area fraction the x and y box sides are calculated
    xsize = sqrt( M_PI / 4 / areafrac * ((double)Nlarge + sizeratio*sizeratio*(N - Nlarge)) ) ;
    ysize = xsize ;
    hx = hy = 0.5 * xsize ;
    //the cell list is initialized to fastly check for overlaps when inserting new particles
    cx = (int)(xsize - 0.0001) / 1.1 ;
    cy = (int)(ysize - 0.0001) / 1.1 ;
    while (cx*cy > 8*N) {
        cx *= 0.9;
        cy *= 0.9;
    }
    celllist = (particle**) calloc(cx*cy, sizeof(particle*));
    icxsize = cx / xsize ;						//Set inverse cell size
    icysize = cy / ysize ;
    int cellx , celly , overlap ;
    int cdx, cdy ;
    double dx, dy, r2, rm ;
    particle *p2 ;
    printf("Placing particles\n");
    i = 0 ;
    while ( i < N ) {
        p = particles + i ;
	p->x = genrand_real2() * xsize ;
	p->y = genrand_real2() * ysize ;
	cellx = p->x * icxsize + cx ;
	celly = p->y * icysize + cy ;
	overlap = 0 ;
	for (cdx = cellx - 1; cdx < cellx + 2; cdx++)
	  for (cdy = celly - 1; cdy < celly + 2; cdy++) {
            p2 = celllist[celloffset(cdx % cx, cdy % cy)];
            while (p2) {
	      if( p2 != p ) {
		dx = p->x - p2->x ;
		dy = p->y - p2->y ;
		if (p2->nearboxedge) {
		  if (dx > hx) dx -= xsize; else if (dx < -hx) dx += xsize;  //periodic boundaries
		  if (dy > hy) dy -= ysize; else if (dy < -hy) dy += ysize;
		}
		r2 = dx * dx + dy * dy;
		rm = (p->radius + p2->radius) ;
		if(p->type != p2->type) rm *= nonadditivity ;
		if (r2 < rm * rm) {
		  overlap = 1 ;
		}
	      }
	      p2 = p2->next;
	    }
	  }
        if (overlap == 0) {
	  addtocelllist(p, p->x * icxsize, p->y * icysize) ;
	  i ++ ;
	}
    }

    //cell list is deleted
    for (i = 0; i < N; i++) {
        p = particles + i ;
	removefromcelllist(p) ;
        p->cell = 0 ;
    }
    free(celllist) ;
    celllist = NULL ;
    cx = cy = 0 ;
    icxsize = icysize = 0 ;

    double vfilled = 0 ;
    for (i = 0; i < N; i++) vfilled += particles[i].radius * particles[i].radius ;
    printf("Area packing fraction: %g\n", M_PI * vfilled / (xsize * ysize) ) ;
    printf("Starting configuration from random HD\n") ;
    if (N != Nlarge) {
      printf("Size ratio: %lf\n", sizeratio) ;
      printf("Fraction of large particles (R_L = 0.5): %lf\n", large2totalfraction) ;
    }
}




/**************************************************
**             RANDOMSTAMPFLI
** Stampfli QC12 configuration with randomly oriented wheels
**************************************************/
void randomStampfli()
{
  particle *p, **partarray ;
  double hardcoreradius = 0.5 ;
  int partarraylength, i, Nparts, Nlarge = 0 ;

  // starting seed
  xsize = 2 + sqrt(3) ;
  ysize = xsize ;
  partarraylength = 15 ;
  partarray = (particle **)calloc(partarraylength, sizeof(particle *)) ;
  for(i=0; i<partarraylength; i++) partarray[i] = (particle *)calloc(1, sizeof(particle)) ;
  partarray[0]->x = 1 + 0.5*sqrt(3) ;
  partarray[0]->y = 1 + 0.5*sqrt(3) ;
  for(i=0; i<6; i++) {
    partarray[i+1]->x = partarray[0]->x + cos(M_PI/3*i) ;
    partarray[i+1]->y = partarray[0]->y + sin(M_PI/3*i) ;
  }
  for(i=0; i<4; i++) {
    partarray[2*i+7]->x = partarray[0]->x + sqrt(2 + sqrt(3)) * cos(M_PI*((float)i/2+1./12)) ;
    partarray[2*i+7]->y = partarray[0]->y + sqrt(2 + sqrt(3)) * sin(M_PI*((float)i/2+1./12)) ;
    partarray[2*i+1+7]->x = partarray[0]->x + sqrt(2 + sqrt(3)) * cos(M_PI*((float)i/2+3./12)) ;
    partarray[2*i+1+7]->y = partarray[0]->y + sqrt(2 + sqrt(3)) * sin(M_PI*((float)i/2+3./12)) ;
  }
  Nparts = 15 ;

  while ( Nparts < N ) {
    xsize *= 2 + sqrt(3) ;
    ysize *= 2 + sqrt(3) ;
    hx = 0.5 * xsize ;
    hy = 0.5 * ysize ;
    for(i=0; i<partarraylength; i++) {
      partarray[i]->x *= ( 2 + sqrt(3) ) ;
      partarray[i]->y *= ( 2 + sqrt(3) ) ;
      backinbox( partarray[i] ) ;
    }
    partarray = (particle **)realloc(partarray, 19*partarraylength*sizeof(particle *)) ;
    for(i=partarraylength; i<19*partarraylength; i++) partarray[i] = (particle *)calloc(1, sizeof(particle)) ;
    for(i=0; i<partarraylength; i++) {
      int j ;
      float offset = (genrand_real2()<0.5)*1./6 ;
      for(j=0; j<6; j++) {
	partarray[i*18+partarraylength+j]->x = partarray[i]->x + cos(M_PI * ((float)j/3 + offset)) ;
	partarray[i*18+partarraylength+j]->y = partarray[i]->y + sin(M_PI * ((float)j/3 + offset)) ;
	backinbox( partarray[i*18+partarraylength+j] ) ;
      }
      for(j=0; j<12; j++) {
	partarray[i*18+partarraylength+j+6]->x = partarray[i]->x + sqrt(2 + sqrt(3)) * cos(M_PI*((float)j/6+1./12)) ;
	partarray[i*18+partarraylength+j+6]->y = partarray[i]->y + sqrt(2 + sqrt(3)) * sin(M_PI*((float)j/6+1./12)) ;
	backinbox( partarray[i*18+partarraylength+j+6] ) ;
      }
    }
    partarraylength *= 19 ;

    //the cell list is initialized to fastly check for overlaps when inserting new particles
    cx = (int)(xsize - 0.0001) / 1.1 ;
    cy = (int)(ysize - 0.0001) / 1.1 ;
    while (cx*cy > 8*partarraylength) {
        cx *= 0.9;
        cy *= 0.9;
    }
    celllist = (particle**) calloc(cx*cy, sizeof(particle*));
    icxsize = cx / xsize ;						//Set inverse cell size
    icysize = cy / ysize ;
    int cellx , celly , keepgoing ;
    int cdx, cdy , overlap = 0 ;
    double dx, dy, r2 ;
    particle *p2 ;
    for(i=0; i<partarraylength; i++) {
      backinbox( partarray[i] ) ;
      addtocelllist(partarray[i], partarray[i]->x * icxsize, partarray[i]->y * icysize) ;
    }
    for(i=0; i<partarraylength; i++) {
      if( partarray[i] != NULL ) {
	keepgoing = 1 ;
	cellx = partarray[i]->x * icxsize + cx ;
	celly = partarray[i]->y * icysize + cy ;
	for (cdx = cellx - 1; cdx < cellx + 2 && keepgoing; cdx++) {
	  for (cdy = celly - 1; cdy < celly + 2 && keepgoing; cdy++) {
	    p2 = celllist[celloffset(cdx % cx, cdy % cy)];
	    while (p2) {
	      if( p2 != partarray[i] ) {
		dx = partarray[i]->x - p2->x ;
		dy = partarray[i]->y - p2->y ;
		if (dx > hx) dx -= xsize; else if (dx < -hx) dx += xsize;  //periodic boundaries
		if (dy > hy) dy -= ysize; else if (dy < -hy) dy += ysize;
		r2 = dx * dx + dy * dy;
		if (r2 < 0.1) {
		  removefromcelllist(partarray[i]) ;
		  partarray[i]->cell = 0 ;
		  free( partarray[i] ) ;
		  overlap += 1 ;
		  partarray[i] = partarray[partarraylength-overlap] ;
		  partarray[partarraylength-overlap] = NULL ;
		  i -- ;
		  keepgoing = 0 ;
		  p2 = NULL ;
		} else p2 = p2->next;
	      } else p2 = p2->next;
	    }
	  }
	}
      }
    }
    Nparts = partarraylength - overlap ;
    partarray = (particle **)realloc(partarray, Nparts*sizeof(particle *)) ;
    partarraylength = Nparts ;
    Nlarge = Nparts ;

    //cell list is deleted
    for (i = 0; i < partarraylength; i++) {
      p = partarray[i] ;
      removefromcelllist(p) ;
      p->cell = 0 ;
    }
    free(celllist) ;
    celllist = NULL ;
    cx = cy = 0 ;
    icxsize = icysize = 0 ;
  }

  //the cell list is initialized to fastly check for overlaps when inserting new particles
  cx = (int)(xsize - 0.0001) / 1.1 ;
  cy = (int)(ysize - 0.0001) / 1.1 ;
  while (cx*cy > 8*partarraylength) {
    cx *= 0.9;
    cy *= 0.9;
  }
  celllist = (particle**) calloc(cx*cy, sizeof(particle*));
  icxsize = cx / xsize ;						//Set inverse cell size
  icysize = cy / ysize ;
  int cellx , celly ;
  int cdx, cdy ;
  double dx, dy, r2 ;
  particle *p2 ;
  int **bondlist = (int **)calloc(partarraylength, sizeof(int *)) ;
  for(i=0; i<partarraylength; i++) {
    backinbox( partarray[i] ) ;
    addtocelllist(partarray[i], partarray[i]->x * icxsize, partarray[i]->y * icysize) ;
    bondlist[i] = (int *)calloc(6, sizeof(int)) ;
    int j;
    for( j=0; j<6; j++) bondlist[i][j] = -1 ;
    partarray[i]->idx = i ;
  }
  for(i=0; i<partarraylength; i++) {
    cellx = partarray[i]->x * icxsize + cx ;
    celly = partarray[i]->y * icysize + cy ;
    for (cdx = cellx - 1; cdx < cellx + 2; cdx++) {
      for (cdy = celly - 1; cdy < celly + 2; cdy++) {
	p2 = celllist[celloffset(cdx % cx, cdy % cy)];
	while (p2) {
	  if( p2 != partarray[i] ) {
	    dx = partarray[i]->x - p2->x ;
	    dy = partarray[i]->y - p2->y ;
	    if (dx > hx) dx -= xsize; else if (dx < -hx) dx += xsize;  //periodic boundaries
	    if (dy > hy) dy -= ysize; else if (dy < -hy) dy += ysize;
	    r2 = dx * dx + dy * dy;
	    if (r2 < 1.1) {
	      int j = 0 ;
	      while(bondlist[i][j] != -1) j++ ;
	      bondlist[i][j] = p2->idx ;
	    }
	  }
	  p2 = p2->next;
	}
      }
    }
  }
  partarray = (particle **)realloc(partarray, 2*partarraylength*sizeof(particle *)) ;
  partarraylength *= 2 ;
  for(i=0; i<Nlarge; i++) {
    int j1, j2, i1, i2;
    for(j1=0; j1<6; j1++) {
      if( bondlist[i][j1] != -1 ) {
	i1 = bondlist[i][j1] ;
	double dx1 = partarray[i1]->x - partarray[i]->x ;
	double dy1 = partarray[i1]->y - partarray[i]->y ;
	if (dx1 > hx) dx1 -= xsize; else if (dx1 < -hx) dx1 += xsize;  //periodic boundaries
	if (dy1 > hy) dy1 -= ysize; else if (dy1 < -hy) dy1 += ysize;
	for(j2=0; j2<6; j2++) {
	  if( bondlist[i1][j2] != -1 && bondlist[i1][j2] != i ) {
	    i2 = bondlist[i1][j2] ;
	    double dx2 = partarray[i2]->x - partarray[i1]->x ;
	    double dy2 = partarray[i2]->y - partarray[i1]->y ;
	    if (dx2 > hx) dx2 -= xsize; else if (dx2 < -hx) dx2 += xsize;  //periodic boundaries
	    if (dy2 > hy) dy2 -= ysize; else if (dy2 < -hy) dy2 += ysize;
	    double orientedarea = dx1 * dy2 - dx2 * dy1 ;
	    if( orientedarea < 1.01 && orientedarea > 0.99 && i < i1 && i < i2 ) {
	      partarray[Nparts] = (particle *)calloc(1, sizeof(particle)) ;
	      partarray[Nparts]->x = partarray[i]->x + 0.5 * ( dx1 + dx2 ) ;
	      partarray[Nparts]->y = partarray[i]->y + 0.5 * ( dy1 + dy2 ) ;
	      backinbox( partarray[Nparts] ) ;
	      addtocelllist(partarray[Nparts], partarray[Nparts]->x * icxsize, partarray[Nparts]->y * icysize) ;
	      Nparts ++ ;
	    }
	  }
	}
      }
    }
  }
  // deleting overlapping small particles
  int overlap = 0 , keepgoing ;
  for(i=Nlarge; i<Nparts; i++) {
    if( partarray[i] != NULL ) {
      keepgoing = 1 ;
      cellx = partarray[i]->x * icxsize + cx ;
      celly = partarray[i]->y * icysize + cy ;
      for (cdx = cellx - 1; cdx < cellx + 2 && keepgoing; cdx++) {
	for (cdy = celly - 1; cdy < celly + 2 && keepgoing; cdy++) {
	  p2 = celllist[celloffset(cdx % cx, cdy % cy)];
	  while (p2 && keepgoing) {
	    if( p2 != partarray[i] ) {
	      dx = partarray[i]->x - p2->x ;
	      dy = partarray[i]->y - p2->y ;
	      if (dx > hx) dx -= xsize; else if (dx < -hx) dx += xsize;  //periodic boundaries
	      if (dy > hy) dy -= ysize; else if (dy < -hy) dy += ysize;
	      r2 = dx * dx + dy * dy;
	      if (r2 < 0.1) {
		removefromcelllist(partarray[i]) ;
		partarray[i]->cell = 0 ;
		free( partarray[i] ) ;
		overlap += 1 ;
		partarray[i] = partarray[Nparts-overlap] ;
		partarray[Nparts-overlap] = NULL ;
		i -- ;
		keepgoing = 0 ;
		p2 = NULL ;
	      } else p2 = p2->next;
	    } else p2 = p2->next;
	  }
	}
      }
    }
  }
  Nparts = Nparts - overlap ;
  //cell list is deleted
  for (i = 0; i < Nparts; i++) {
    p = partarray[i] ;
    removefromcelllist(p) ;
    p->cell = 0 ;
  }
  free(celllist) ;
  celllist = NULL ;
  cx = cy = 0 ;
  icxsize = icysize = 0 ;
  // bondlist is deleted
  for(i=0; i<Nlarge; i++) free( bondlist[i] ) ;
  free( bondlist ) ;
  bondlist = NULL ;

  //given the area fraction the x and y box sides are calculated
  double vfilled = (float)Nlarge * hardcoreradius * hardcoreradius + (float)(Nparts-Nlarge) * hardcoreradius * hardcoreradius * sizeratio * sizeratio ;
  double infareafrac = M_PI * vfilled / (xsize * ysize) ;
  double scale = sqrt( infareafrac / areafrac ) ;

  // initializing particles
  N = Nparts ;
  large2totalfraction = (double)Nlarge / N ;
  initparticles(N);
  for (i = 0; i < Nlarge; i++) {
    p = particles + i ;
    p->x = partarray[i]->x * scale ;
    p->y = partarray[i]->y * scale ;
    p->radius = hardcoreradius ;
    p->type = 0 ;
    p->mass = 1.0 ;
    p->idx = i ;
  }
  for (i = Nlarge; i < N; i++) {
    p = particles + i ;
    p->x = partarray[i]->x * scale ;
    p->y = partarray[i]->y * scale ;
    p->radius = hardcoreradius * sizeratio ;
    p->type = 1 ;
    p->mass = sizeratio * sizeratio ;
    p->idx = i ;
  }
  xsize *= scale ;
  ysize *= scale ;
  // deleting partarray
  for (i = 0; i < N; i++) {
    backinbox( particles + i ) ;
    free( partarray[i] ) ;
  }
  free( partarray ) ;
  partarray = NULL ;
  partarraylength = 0 ;

  printf("Area packing fraction: %g\n", M_PI * vfilled / (xsize * ysize) ) ;
  printf("Starting from a quasi-random Stampfli QC12 configuration\n") ;
  printf("Size ratio: %lf\n", sizeratio) ;
  printf("Fraction of large particles (R_L = 0.5): %lf\n", large2totalfraction) ;
}




/**************************************************
**             ADDSMALLPARTICLES
** Adds to a QC12 configuration of large particles the small ones in the center of squares
**************************************************/
void addsmallparticles()
{
  particle *p = NULL, **partarray ;
  double hardcoreradius = 0.5 ;
  int partarraylength, i, Nparts = N, Nlarge = N ;

  partarraylength = N ;
  partarray = (particle **)calloc(partarraylength, sizeof(particle *)) ;
  for(i=0; i<partarraylength; i++) {
    partarray[i] = (particle *)calloc(1, sizeof(particle)) ;
    partarray[i]->x = particles[i].x ;
    partarray[i]->y = particles[i].y ;
    partarray[i]->idx = i ;
  }

  //the cell list is initialized to fastly check for overlaps when inserting new particles
  hx = 0.5 * xsize ;
  hy = 0.5 * ysize ;
  cx = (int)(xsize - 0.0001) / 2.0 ;
  cy = (int)(ysize - 0.0001) / 2.0 ;
  while (cx*cy > 8*partarraylength) {
    cx *= 0.9;
    cy *= 0.9;
  }
  celllist = (particle**) calloc(cx*cy, sizeof(particle*));
  icxsize = cx / xsize ;						//Set inverse cell size
  icysize = cy / ysize ;
  int cellx , celly ;
  int cdx, cdy ;
  double dx, dy, r2 ;
  particle *p2 ;

  // building of the list of bonds among large particles
  int **bondlist = (int **)calloc(partarraylength, sizeof(int *)) ;
  for(i=0; i<partarraylength; i++) {
    backinbox( partarray[i] ) ;
    addtocelllist(partarray[i], partarray[i]->x * icxsize, partarray[i]->y * icysize) ;
    bondlist[i] = (int *)calloc(6, sizeof(int)) ;
    int j;
    for( j=0; j<6; j++) bondlist[i][j] = -1 ;
  }

  // calculation of the maximum reachable packing fraction
  for (i = 0; i < partarraylength; i++) {
    partarray[i]->radius = hardcoreradius ;
    partarray[i]->type = 0 ;
  }
  double scale = 1 , md ;
  for(i=0; i<partarraylength; i++) {
    cellx = partarray[i]->x * icxsize + cx ;
    celly = partarray[i]->y * icysize + cy ;
    for (cdx = cellx - 1; cdx < cellx + 2; cdx++) {
      for (cdy = celly - 1; cdy < celly + 2; cdy++) {
	p2 = celllist[celloffset(cdx % cx, cdy % cy)];
	while (p2) {
	  if( p2 != partarray[i] ) {
	    dx = partarray[i]->x - p2->x ;
	    dy = partarray[i]->y - p2->y ;
	    if (dx > hx) dx -= xsize; else if (dx < -hx) dx += xsize;  //periodic boundaries
	    if (dy > hy) dy -= ysize; else if (dy < -hy) dy += ysize;
	    r2 = dx * dx + dy * dy;
	    md = partarray[i]->radius + p2->radius;
	    if (partarray[i]->type != p2->type) md *= nonadditivity ;
	    if (r2 / md / md < scale) scale = r2 / md ;
	  }
	  p2 = p2->next;
	}
      }
    }
  }
  scale = sqrt( scale ) ;
  if( scale < 1 ) scale = 1. / scale ;
  xsize *= scale ;
  ysize *= scale ;
  hx = 0.5 * xsize ;
  hy = 0.5 * ysize ;
  for (i = 0; i < partarraylength; i++) {
    partarray[i]->x *= scale ;
    partarray[i]->y *= scale ;
    backinbox( partarray[i] ) ;
  }
  // celllist rebuilding
  for (i = 0; i < Nparts; i++) {
    removefromcelllist( partarray[i] ) ;
    partarray[i]->cell = 0 ;
  }
  free(celllist) ;
  celllist = NULL ;
  cx = cy = 0 ;
  icxsize = icysize = 0 ;
  cx = (int)(xsize - 0.0001) / 2.0 ;
  cy = (int)(ysize - 0.0001) / 2.0 ;
  while (cx*cy > 8*partarraylength) {
    cx *= 0.9;
    cy *= 0.9;
  }
  celllist = (particle**) calloc(cx*cy, sizeof(particle*));
  icxsize = cx / xsize ;						//Set inverse cell size
  icysize = cy / ysize ;
  for(i=0; i<partarraylength; i++) addtocelllist(partarray[i], partarray[i]->x * icxsize, partarray[i]->y * icysize) ;
  // celllist rebuilt

  // building of the neighbour list
  for(i=0; i<partarraylength; i++) {
    cellx = partarray[i]->x * icxsize + cx ;
    celly = partarray[i]->y * icysize + cy ;
    for (cdx = cellx - 1; cdx < cellx + 2; cdx++) {
      for (cdy = celly - 1; cdy < celly + 2; cdy++) {
	p2 = celllist[celloffset(cdx % cx, cdy % cy)];
	while (p2) {
	  if( p2 != partarray[i] ) {
	    dx = partarray[i]->x - p2->x ;
	    dy = partarray[i]->y - p2->y ;
	    if (dx > hx) dx -= xsize; else if (dx < -hx) dx += xsize;  //periodic boundaries
	    if (dy > hy) dy -= ysize; else if (dy < -hy) dy += ysize;
	    r2 = dx * dx + dy * dy;
	    if (sqrt(r2) < 1.1) {
	      int j = 0 ;
	      while(bondlist[i][j] != -1) j++ ;
	      if(j>5) printf("ERROR in small particles creation");
	      bondlist[i][j] = p2->idx ;
	    }
	  }
	  p2 = p2->next;
	}
      }
    }
  }

  // addition of the small particles in the center of squares
  partarray = (particle **)realloc(partarray, 2*partarraylength*sizeof(particle *)) ;
  partarraylength *= 2 ;
  for(i=0; i<Nlarge; i++) {
    int j1, j2, i1, i2;
    for(j1=0; j1<6; j1++) {
      if( bondlist[i][j1] != -1 ) {
	i1 = bondlist[i][j1] ;
	double dx1 = partarray[i1]->x - partarray[i]->x ;
	double dy1 = partarray[i1]->y - partarray[i]->y ;
	if (dx1 > hx) dx1 -= xsize; else if (dx1 < -hx) dx1 += xsize;  //periodic boundaries
	if (dy1 > hy) dy1 -= ysize; else if (dy1 < -hy) dy1 += ysize;
	for(j2=0; j2<6; j2++) {
	  if( bondlist[i1][j2] != -1 && bondlist[i1][j2] != i ) {
	    i2 = bondlist[i1][j2] ;
	    double dx2 = partarray[i2]->x - partarray[i1]->x ;
	    double dy2 = partarray[i2]->y - partarray[i1]->y ;
	    if (dx2 > hx) dx2 -= xsize; else if (dx2 < -hx) dx2 += xsize;  //periodic boundaries
	    if (dy2 > hy) dy2 -= ysize; else if (dy2 < -hy) dy2 += ysize;
	    double orientedarea = dx1 * dy2 - dx2 * dy1 ;
	    if( orientedarea < 1.05 && orientedarea > 0.95 && i < i1 && i < i2 ) {
	      partarray[Nparts] = (particle *)calloc(1, sizeof(particle)) ;
	      partarray[Nparts]->x = partarray[i]->x + 0.5 * ( dx1 + dx2 ) ;
	      partarray[Nparts]->y = partarray[i]->y + 0.5 * ( dy1 + dy2 ) ;
	      backinbox( partarray[Nparts] ) ;
	      partarray[Nparts]->type = 1 ;
	      addtocelllist(partarray[Nparts], partarray[Nparts]->x * icxsize, partarray[Nparts]->y * icysize) ;
	      Nparts ++ ;
	    }
	  }
	}
      }
    }
  }
  // deleting overlapping small particles
  int overlap = 0 , keepgoing ;
  for(i=Nlarge; i<Nparts; i++) {
    if( partarray[i] != NULL ) {
      keepgoing = 1 ;
      cellx = partarray[i]->x * icxsize + cx ;
      celly = partarray[i]->y * icysize + cy ;
      for (cdx = cellx - 1; cdx < cellx + 2 && keepgoing; cdx++) {
	for (cdy = celly - 1; cdy < celly + 2 && keepgoing; cdy++) {
	  p2 = celllist[celloffset(cdx % cx, cdy % cy)];
	  while (p2 && keepgoing) {
	    if( p2 != partarray[i] ) {
	      dx = partarray[i]->x - p2->x ;
	      dy = partarray[i]->y - p2->y ;
	      if (dx > hx) dx -= xsize; else if (dx < -hx) dx += xsize;  //periodic boundaries
	      if (dy > hy) dy -= ysize; else if (dy < -hy) dy += ysize;
	      r2 = dx * dx + dy * dy;
	      if (r2 < 0.1) {
		removefromcelllist(partarray[i]) ;
		partarray[i]->cell = 0 ;
		free( partarray[i] ) ;
		overlap += 1 ;
		partarray[i] = partarray[Nparts-overlap] ;
		partarray[Nparts-overlap] = NULL ;
		i -- ;
		keepgoing = 0 ;
		p2 = NULL ;
	      } else p2 = p2->next;
	    } else p2 = p2->next;
	  }
	}
      }
    }
  }
  Nparts = Nparts - overlap ;

  // calculation of the maximum reachable packing fraction WITH SMALL PARTICLES, too
  scale = 1 ;
  for(i=0; i<Nparts; i++) {
    cellx = partarray[i]->x * icxsize + cx ;
    celly = partarray[i]->y * icysize + cy ;
    for (cdx = cellx - 1; cdx < cellx + 2; cdx++) {
      for (cdy = celly - 1; cdy < celly + 2; cdy++) {
	p2 = celllist[celloffset(cdx % cx, cdy % cy)];
	while (p2) {
	  if( p2 != partarray[i] ) {
	    dx = partarray[i]->x - p2->x ;
	    dy = partarray[i]->y - p2->y ;
	    if (dx > hx) dx -= xsize; else if (dx < -hx) dx += xsize;  //periodic boundaries
	    if (dy > hy) dy -= ysize; else if (dy < -hy) dy += ysize;
	    r2 = dx * dx + dy * dy;
	    md = partarray[i]->radius + p2->radius;
	    if (partarray[i]->type != p2->type) md *= nonadditivity ;
	    if (r2 / md / md < scale) scale = r2 / md ;
	  }
	  p2 = p2->next;
	}
      }
    }
  }
  scale = sqrt( scale ) ;
  if( scale < 1 ) scale = 1. / scale ;
  for (i = 0; i < Nparts; i++) {
    partarray[i]->x *= scale ;
    partarray[i]->y *= scale ;
  }
  xsize *= scale ;
  ysize *= scale ;
  hx = 0.5 * xsize ;
  hy = 0.5 * ysize ;
  //cell list is deleted
  for (i = 0; i < Nparts; i++) {
    p = partarray[i] ;
    removefromcelllist(p) ;
    p->cell = 0 ;
  }
  free(celllist) ;
  celllist = NULL ;
  cx = cy = 0 ;
  icxsize = icysize = 0 ;
  // bondlist is deleted
  for(i=0; i<Nlarge; i++) free( bondlist[i] ) ;
  free( bondlist ) ;
  bondlist = NULL ;

  // initializing particles
  N = Nparts ;
  large2totalfraction = (double)Nlarge / N ;
  free( particles ) ;
  particles = NULL ;
  initparticles(N) ;
  for (i = 0; i < Nlarge; i++) {
    p = particles + i ;
    p->x = partarray[i]->x ;
    p->y = partarray[i]->y ;
    p->radius = hardcoreradius ;
    p->type = 0 ;
    p->mass = 1.0 ;
    p->idx = i ;
  }
  for (i = Nlarge; i < N; i++) {
    p = particles + i ;
    p->x = partarray[i]->x ;
    p->y = partarray[i]->y ;
    p->radius = hardcoreradius * sizeratio ;
    p->type = 1 ;
    p->mass = sizeratio * sizeratio ;
    p->idx = i ;
  }
  // deleting partarray
  for (i = 0; i < N; i++) {
    backinbox( particles + i ) ;
    free( partarray[i] ) ;
  }
  free( partarray ) ;
  partarray = NULL ;
  partarraylength = 0 ;

  double vfilled = (float)Nlarge * hardcoreradius * hardcoreradius + (float)(Nparts-Nlarge) * hardcoreradius * hardcoreradius * sizeratio * sizeratio ;
  double infareafrac = M_PI * vfilled / (xsize * ysize) ;
  printf("Maximum area packing fraction: %.16g\n", infareafrac ) ;
  if( areafrac < infareafrac ) {
    scale = sqrt( infareafrac / areafrac ) ;
    for (i = 0; i < N; i++) {
      p = particles + i ;
      p->x *= scale ;
      p->y *= scale ;
    }
    xsize *= scale ;
    ysize *= scale ;
    hx = 0.5 * xsize ;
    hy = 0.5 * ysize ;
  }
  printf("Area packing fraction: %.16g\n", M_PI * vfilled / (xsize * ysize) ) ;
  printf("Generatin random QC12 configuration\n") ;
  printf("Size ratio: %lf\n", sizeratio) ;
  printf("Fraction of large particles (R_L = 0.5): %lf\n", large2totalfraction) ;
}




/**************************************************
**                RANDOMMOVEMENT
** Assign random velocities to all particles
** Center-of-mass velocity = 0
** Kinetic energy per particle = kT
**************************************************/
void randommovement()
{
    particle* p;
    double v2tot = 0, vxtot = 0, vytot = 0;
    double mtot = 0;
    int i;

    for (i = 0; i < N; i++)
    {
        p = particles + i;
        double imsq = 1.0 / sqrt(p->mass);

        p->vx = imsq * random_gaussian();
        p->vy = imsq * random_gaussian();
        vxtot += p->mass * p->vx;					//Keep track of total v
        vytot += p->mass * p->vy;
        mtot += p->mass;
    }


    vxtot /= mtot; vytot /= mtot;
    for (i = 0; i < N; i++)
    {
        p = &(particles[i]);
        p->vx -= vxtot;					//Make sure v_cm = 0
        p->vy -= vytot;
        v2tot += p->mass * (p->vx * p->vx + p->vy * p->vy);
    }
    double fac = sqrt(2.0 / (v2tot / N));
    for (i = 0; i < N; i++)
    {
        p = &(particles[i]);
        p->vx *= fac;					//Fix energy
        p->vy *= fac;
    }
}




/**************************************************
**                LOADVELOCITIES
** Read in velocities from a file
**************************************************/
void loadvelocities(void)
{
    particle* p;
    char *buffer , c , dummystring[100] ;
    int len = 1 , i , nparts = 0 ;
    FILE *file ;
    buffer = (char*)calloc( len , sizeof(char) ) ;
    buffer[0] = '\0' ;

    file = fopen( inputvelfile , "r" ) ;
    if ( !file ) {
      printf("File %s not found, starting with random velocities\n" , inputvelfile);
      randommovement();

    } else {
      myreadline( &buffer , &len , file ) ;  // line containing the number of particles in the system
      sscanf( buffer , "%c %c %d" , &c , &c , &nparts ) ;
      if( nparts != N ) { printf( "Error in the file format!" ) ; exit(3); }
      myreadline( &buffer , &len , file ) ;  // line with timestep information
      double t0 = 0 ;
      sscanf( buffer , "%c %s %lf" , &c , dummystring , &t0 ) ;
      if( t0 > 0 ) {
	simtime = t0 ;
	timewindow = (int)(simtime / simtimewindowlength) ;
	simtime -= simtimewindowlength * timewindow ;
	printf( "Read time from velocity file: %lf\n" , t0 ) ;
      }

      for (i = 0; i < N; i++)
	{
	  myreadline( &buffer , &len , file ) ;
	  p = particles + i;
	  sscanf( buffer , "%lf %lf" , &(p->vx) , &(p->vy) ) ;
	}
    }
}

/**************************************************
**                UPDATE
** Update particle to the current time
**************************************************/
void update(particle* p1)
{
    double dt = simtime - p1->t + simtimewindowlength * (timewindow - p1->timewindow) ;
    p1->t = simtime;
    p1->timewindow = timewindow;
    p1->x += dt * p1->vx;
    p1->y += dt * p1->vy;
}

/**************************************************
**                UPDATEDPARTICLE
** Update particle p1 to the current time and puts
**   it into p2, without modifying p1
**************************************************/
inline void updatedparticle(particle* p1, particle* p2)
{
    p2->x = p1->x;
    p2->y = p1->y;
    p2->vx = p1->vx;
    p2->vy = p1->vy;
    double dt = simtime - p1->t + simtimewindowlength * (timewindow - p1->timewindow) ;
    if (dt > 0) {
      p2->x += dt * p1->vx;
      p2->y += dt * p1->vy;
    }
    p2->t = simtime;
    p2->timewindow = timewindow ;
}

/**************************************************
**                 INITCELLLIST
** Initialize the cell list
** Cell size should be at least:
**    shellsize * [diameter of largest sphere]
**    or sqrt(2./pi*4/(N/A)+1) * [diameter of largest sphere] , if the density is low enough
** Can use larger cells at low densities
**************************************************/
void initcelllist()
{
    int i;
    double maxdiameter = 0.0 ;
    for (i = 0; i < N; i++) 
    {
        particle* p = particles + i;
        if (p->radius > maxdiameter) maxdiameter = p->radius;
    }
    maxdiameter *= 2. ;
    if( shellsize < sqrt(2./M_PI*4/N*xsize*ysize + 1.0) / maxdiameter )  shellsize = sqrt(2./M_PI*4/N*xsize*ysize + 1.0) / maxdiameter ;
    cx = (int)(xsize - 0.0001) / shellsize / maxdiameter;				//Set number of cells
    cy = (int)(ysize - 0.0001) / shellsize / maxdiameter;

    while (cx*cy > 4*N)
    {
        cx *= 0.9;
        cy *= 0.9;
    }

    printf("Cells: %d, %d\n", cx, cy);
    celllist = (particle**) calloc(cx*cy, sizeof(particle*));

    icxsize = cx / xsize;						//Set inverse cell size
    icysize = cy / ysize;
    for (i = 0; i < N; i++) 
    {
        particle* p = particles + i;
	int xcell = p->x * icxsize ;
	int ycell = p->y * icysize ;
	if( xcell == cx ) xcell = 0 ;
	if( ycell == cy ) ycell = 0 ;
        addtocelllist(p, xcell, ycell);
    }
}

/**************************************************
**               CELLOFFSET
** Convert x y-coordinates of cell to
** a single-integer cell index
**************************************************/
int celloffset(int a, int b)
{
    return a + b*cx;
}

/**************************************************
**                    ADDTOCELLLIST
** Add particle to the cell list at cell index
** (cellx, celly)
** Note that each cell is implemented as a doubly
** linked list.
**************************************************/
void addtocelllist(particle* p, int cellx, int celly)
{
    p->cell = celloffset(cellx, celly);
    p->next = celllist[p->cell];	//Add particle to celllist
    if (p->next) p->next->prev = p;			//Link up list
    celllist[p->cell] = p;
    p->prev = NULL;
    //Check if particle is near the box edge, where nearest image convention should be checked
    p->nearboxedge = (cellx == 0 || celly == 0 || cellx == cx - 1 || celly == cy - 1);

}

/**************************************************
**               REMOVEFROMCELLLIST
**************************************************/
void removefromcelllist(particle* p1)
{
    if (p1->prev) p1->prev->next = p1->next;    //Remove particle from celllist
    else          celllist[p1->cell] = p1->next;
    if (p1->next) p1->next->prev = p1->prev;
}

/**************************************************
**                     STEP
** Process a single event
**************************************************/
void step()
{
    particle* ev;
    ev = root->right;
    while (ev == NULL)                  //Need to include next event list?
    {
        addnexteventlist();
        ev = root->right;
    }

    while (ev->left) ev = ev->left;		//Find first event

    if ( ev->eventtime + simtimewindowlength * ev->eventtimewindow < simtime + simtimewindowlength * timewindow ) {
      printf("\n *** Negative time step at event (simtime,simtime-time,type) :\t%g ,\t%g ,\t%d\n\n" , simtime , simtime - ev->eventtime + simtimewindowlength * (timewindow - ev->eventtimewindow) , ev->eventtype) ;
      if( ev->eventtype == 0 ) printf("%d %lf %lf, %lf %lf, %lf %d %lf; %d %lf %lf, %lf %lf, %lf\n\n", ev->type, ev->x, ev->y, ev->vx, ev->vy, ev->eventtime, ev->eventtimewindow, simtimewindowlength, ev->p2->type, ev->p2->x, ev->p2->y, ev->p2->vx, ev->p2->vy, ev->p2->eventtime);
    }
    simtime = ev->eventtime;
    timewindow = ev->eventtimewindow ;
    removeevent(ev);
    switch(ev->eventtype)
    {
        case 0:
            corescollision(ev);
            break;
        case 8:
            makeneighborlist(ev);
            break;
        case 100:
            write(ev);
            break;
        case 101:
            thermostat(ev);
            break;
        case 102:
            dumpsnapshot(ev);
            break;
    }
}



/**************************************************
**                MAKENEIGHBORLIST
** Update neighbor list for a single particle
** In this version cylindrical cells are used to 
**   make the neighbor list
**************************************************/
void makeneighborlist(particle* p1)
{

    int cdx, cdy;
    particle* p2;
    double dx, dy, r2, rm;

    update(p1);

    //Put particle back in the box
    if (p1->x >= xsize) { p1->x -= xsize; p1->boxestraveledx++; }
    else if (p1->x < 0) { p1->x += xsize; p1->boxestraveledx--; }
    if (p1->y >= ysize) { p1->y -= ysize; p1->boxestraveledy++; }
    else if (p1->y < 0) { p1->y += ysize; p1->boxestraveledy--; }

    p1->xn = p1->x;     //Set Neighbor list position to current position
    p1->yn = p1->y;

    removefromcelllist(p1);

    //Remove particle i from the neighborlist of its old neighbors
    int i, j;
    for (i = 0; i < p1->nneigh; i++)
    {
        p2 = p1->neighbors[i];
        for (j = 0; j < p2->nneigh; j++)
        {
            if (p2->neighbors[j] == p1)
            {
                p2->nneigh--;
                p2->neighbors[j] = p2->neighbors[p2->nneigh];
                break;
            }
        }
    }

    int cellx = p1->x * icxsize, celly = p1->y * icysize;
    addtocelllist(p1, cellx, celly);

    cellx += cx;
    celly += cy;
    p1->nneigh = 0;

    for (cdx = cellx - 1; cdx < cellx + 2; cdx++)
        for (cdy = celly - 1; cdy < celly + 2; cdy++)
	{
            p2 = celllist[celloffset(cdx % cx, cdy % cy)];
            while (p2)
	    {
	      if (p2 != p1)
		{
		  update(p2);
		  dx = p1->xn - p2->xn;
		  dy = p1->yn - p2->yn;
		  if (p1->nearboxedge)
		  {
		    if (dx > hx) dx -= xsize; else if (dx < -hx) dx += xsize;  //periodic boundaries
		    if (dy > hy) dy -= ysize; else if (dy < -hy) dy += ysize;
		  }
		  r2 = dx * dx + dy * dy;
		  rm = (p1->radius + p2->radius) * shellsize;
		  // if (p1->type != p2->type) rm *= nonadditivity ; // wrong
		  if (r2 < rm * rm)       //infinite vertical cylinders overlapping condition
		  {
		    if (p1->nneigh >= MAXNEIGH || p2->nneigh >= MAXNEIGH)
		    {
		      printf("Too many neighbors. Increase MAXNEIGH.\n");
		      exit(3);
		    }
		    p1->neighbors[p1->nneigh++] = p2;
		    p2->neighbors[p2->nneigh++] = p1;
		  }
		}
	      p2 = p2->next;
	    }
	}

    findcollisions(p1);


}


/**************************************************
**                FINDNEIGHBORLISTUPDATE
** Assumes p1 is up to date
** Note that the particle is always in the same
** box as its neighborlist position (p->xn)
** This version considers cylindrical neigh lists
**************************************************/
double findneighborlistupdate(particle* p1)
{
    double dx = p1->x - p1->xn ;
    double dy = p1->y - p1->yn ;

    double dvx = p1->vx, dvy = p1->vy ;

    double dv2 = dvx * dvx + dvy * dvy;
    if (dv2 == 0) return maxtime ;        //q2f: is it worth it?

    double dr2 = dx * dx + dy * dy;
    double md = (shellsize - 1) * p1->radius;
    double b = dx * dvx + dy * dvy ;                  //drh.dvh
    // double disc = b * b - dv2 * (dr2 - md * md) ;
    if (b > 0) {
      double q = - b - sqrt( b * b - dv2 * (dr2 - md * md) ) ;     //should be more numerically robust
      return (dr2 - md * md) / q ;     //time to go out from the lateral surface
    } else {
      double q = - b + sqrt( b * b - dv2 * (dr2 - md * md) ) ;     //should be more numerically robust
      return q / dv2 ;     //time to go out from the lateral surface
    }
    // return (-b + sqrt(disc)) / dv2 ;     //time to go out from the lateral surface
}

/**************************************************
**                FINDCORESCOLLISION
** Detect the next cores collision for two particles
** Note that p1 is always up to date in
** findcorescollision
**************************************************/
int findcorescollision(particle* p1, particle* p2, double* tmin, int* type)
{
    particle p2updated;
    updatedparticle(p2, &p2updated);
    double dx = p1->x - p2updated.x;    //relative distance at current time
    double dy = p1->y - p2updated.y;
    if (p1->nearboxedge)
    {
        if (dx > hx) dx -= xsize; else if (dx < -hx) dx += xsize;  //periodic boundaries
        if (dy > hy) dy -= ysize; else if (dy < -hy) dy += ysize;
    }

    double dvx = p1->vx - p2updated.vx;                               //relative velocity
    double dvy = p1->vy - p2updated.vy;

    double b = dx * dvx + dy * dvy ;                  //dr.dv
    if (b > 0) return 0;                                      //Particles flying apart
    double dr2 = dx * dx + dy * dy ;

    double md = p1->radius + p2->radius;
    if (p1->type != p2->type) md *= nonadditivity ;
    double A = md * md - dr2;
    if (2 * b * *tmin > A) return 0;                        //Approximate check to discard hits at times > tmin

    double dv2 = dvx * dvx + dvy * dvy;
    double disc = b * b + dv2 * A;
    if (disc < 0) return 0;
    double t = A / (b - sqrt(disc)) ;
    if (t < *tmin) {
      *tmin = t;
      *type = 0;
      return 1;
    }

    return 0;
}


/**************************************************
**                FINDCOLLISIONS
** Find all collisions for particle p1.
** The particle 'not' isn't checked.
**************************************************/
void findcollisions(particle* p1)    //All collisions of particle p1
{
    int i;
    double tmin = findneighborlistupdate(p1);
    int type = 8;
    particle* partner = p1;
    particle* p2;
    for (i = 0; i < p1->nneigh; i++)
    {
        p2 = p1->neighbors[i];
        if( findcorescollision(p1, p2, &tmin, &type) )  partner = p2;
    }
    createevent(tmin, p1, partner, type);
    p1->counter2 = partner->counter;
}



/**************************************************
**                FINDALLCOLLISION
** All collisions of all particle pairs
**************************************************/
void findallcollisions()       //All collisions of all particle pairs
{
    int i, j;
    for (i = 0; i < N; i++)
    {
        particle* p1 = particles + i;
        particle* partner = p1;
	particle* p2;
        double tmin = findneighborlistupdate(p1);
        int type = 8;
        for (j = 0; j < p1->nneigh; j++)
        {
            p2 = p1->neighbors[j];
            if (p2 > p1)
            {
	      if(findcorescollision(p1, p2, &tmin, &type))  partner = p2;
            }
        }
        createevent(tmin, p1, partner, type);
        p1->counter2 = partner->counter;
    }
}







/**************************************************
**                  CORESCOLLISION
** Process a single hard core collision event
**************************************************/
void corescollision(particle* p1)
{
    particle* p2 = p1->p2;
    update(p1);
    if (p1->counter2 != p2->counter)
    {
      findcollisions(p1);
        return;
    }

    update(p2);
    p1->counter++;
    p2->counter++;

    double m1 = p1->mass, r1 = p1->radius;
    double m2 = p2->mass, r2 = p2->radius;

    double r = r1 + r2;
    if (p1->type != p2->type) r *= nonadditivity ;
    double rinv = 1.0 / r;
    double dx = (p1->x - p2->x);			//Normalized distance vector
    double dy = (p1->y - p2->y);
    if (p1->nearboxedge)
    {
        if (dx > hx) dx -= xsize; else if (dx < -hx) dx += xsize;  //periodic boundaries
        if (dy > hy) dy -= ysize; else if (dy < -hy) dy += ysize;
    }
    dx *= rinv;  dy *= rinv;

    double dvx = p1->vx - p2->vx;                               //relative velocity
    double dvy = p1->vy - p2->vy;

    double b = dx * dvx + dy * dvy;                  //dr.dv/|dr|
    b *= 2.0 / (m1 + m2);
    double dv1 = b * m2, dv2 = b * m1;

    p1->vx -= dv1 * dx;         //Change velocities after collision
    p1->vy -= dv1 * dy;         //delta v = (-) dx2.dv2
    p2->vx += dv2 * dx;
    p2->vy += dv2 * dy;

    // dvtot += dv1*m1 * r;
    dptot[XX] += dx*r * dv1*m1*dx;
    dptot[YY] += dy*r * dv1*m1*dy;
    dptot[XY] += dx*r * dv1*m1*dy;
    colcounter++;

    removeevent(p2);

    findcollisions(p1);
    findcollisions(p2);
}










/**************************************************
**                 INITEVENTS
**************************************************/
void initevents()
{
    eventlisttime = eventlisttimemultiplier / N;
    numeventlists = ceil(maxscheduletime / eventlisttime);
    maxscheduletime = numeventlists * eventlisttime;
    printf("number of event lists: %d\n", numeventlists);

    eventlists = (particle**)calloc(numeventlists + 1, sizeof(particle*));
    if (!eventlists)
    {
        printf("Failed to allocate memory for eventlists\n");
        exit(3);
    }


    root = particles + N;				//Create root event
    root->eventtime = -99999999999.99;				//Root event is empty, but makes sure every other event has a parent
    root->eventtimewindow = -8;
    root->eventtype = 127;					//This makes sure we don't have to keep checking this when adding/removing events
    root->parent = NULL;

    double time = simtime + simtimewindowlength * timewindow ;
    particle* writeevent = particles + N + 1;		//Set up write event
    if (time == 0) createevent(0, writeevent, NULL, 100) ;
    else createevent(writeinterval, writeevent, NULL, 100) ;

    particle* dumpevent = particles + N + 3;		//Set up write event
    if (time == 0) createevent(0, dumpevent, NULL, 102) ;
    else createevent(snapshotinterval, dumpevent, NULL, 102) ;


    printf("Event tree initialized.\n");

    if (usethermostat)
    {
        particle* thermostatevent = particles + N + 2;
	createevent(thermostatinterval, thermostatevent, NULL, 101) ;
        printf("Started thermostat\n");
    }

}

/**************************************************
**                  ADDEVENTTOTREE
** Add event into binary tree
**************************************************/
void addeventtotree(particle* newevent)
{
    double time = newevent->eventtime;
    int Tw = newevent->eventtimewindow;
    particle* loc = root;
    int busy = 1;
    while (busy)						//Find location to add event into tree (loc)
    {
        if (time - loc->eventtime < simtimewindowlength * (loc->eventtimewindow - Tw))  //Go left  // should keep numerical precision even at long times
        {
            if (loc->left) loc = loc->left;
            else
            {
                loc->left = newevent;
                busy = 0;
            }
        }
        else						//Go right
        {
            if (loc->right) loc = loc->right;
            else
            {
                loc->right = newevent;
                busy = 0;
            }
        }
    }
    newevent->parent = loc;
    newevent->left = NULL;
    newevent->right = NULL;

}

/**************************************************
**                  ADDEVENT
**************************************************/
void addevent(particle* newevent)
{
    double dt = newevent->eventtime - reftime + simtimewindowlength * (newevent->eventtimewindow - reftimewindow) ;

    if (dt < eventlisttime) //Event in near future: Put it in the tree
    {
        newevent->queue = currentlist;
        addeventtotree(newevent);
    }
    else                    //Put it in one of the event lists
    {
        int list_id;
        if (dt >= numeventlists * eventlisttime) list_id = numeventlists;    //This also handles int overflow when calculating list_id
        else
        {
            list_id = currentlist + dt / eventlisttime;
            if (list_id >= numeventlists)
            {
                list_id -= numeventlists;
            }
        }

        newevent->queue = list_id;
        newevent->right = eventlists[list_id]; //Add to linear list
        newevent->left = NULL;
        if (newevent->right) newevent->right->left = newevent;
        eventlists[list_id] = newevent;
    }
}
/**************************************************
**                  CREATEEVENT
**  Here the absolute time of the event is simtime + timeinterval
**************************************************/
void createevent(double timeinterval, particle* p1, particle* p2, int type)
{
    double dt = simtime + timeinterval ;
    int dw = (int) (dt / simtimewindowlength) ;
    p1->eventtime = dt - simtimewindowlength * dw ;
    p1->eventtimewindow = timewindow + dw ;
    p1->eventtype = type;
    p1->p2 = p2;
    addevent(p1);
}

/**************************************************
**                     ADDNEXTEVENTLIST
** Sort all events from the first event list 
** into the binary tree.
**************************************************/
void addnexteventlist()
{
    currentlist++;
    reftime += eventlisttime;
    int dw = (int) (reftime / simtimewindowlength) ;
    reftimewindow += dw ;
    reftime -= simtimewindowlength * dw ;
    if (currentlist == numeventlists) //End of array of event lists?
    {
        currentlist = 0;    //Start at the beginning again

        //Process overflow queue
        particle* ev = eventlists[numeventlists];
        eventlists[numeventlists] = NULL;
        listcounter2 = 0;
        while (ev)
        {
            particle* nextev = ev->right;
            addevent(ev);
            ev = nextev;
            listcounter2 += 1;      //Count how many events there were in overflow
        }
    }

    particle* ev = eventlists[currentlist];
    while (ev)
    {
        particle* nextev = ev->right;
        addeventtotree(ev);
        ev = nextev;
        listcounter1++;         //Count how many events there were in event list
    }
    eventlists[currentlist] = NULL;
    mergecounter++;
}

/**************************************************
**                  REMOVEEVENT
** Remove an event from the event calendar
**************************************************/
void removeevent(particle* oldevent)
{
    if (oldevent->queue != currentlist) //Not in the binary tree
    {
        if (oldevent->right) oldevent->right->left = oldevent->left;
        if (oldevent->left) oldevent->left->right = oldevent->right;
        else
        {
            eventlists[oldevent->queue] = oldevent->right;
        }
        return;
    }

    particle* parent = oldevent->parent;
    particle* node;					//This node will be attached to parent in the end


    if (oldevent->left == NULL)			//Only one child: easy to delete
    {
        node = oldevent->right;			//Child2 is attached to parent
        if (node)
        {
            node->parent = parent;
        }
    }
    else if (oldevent->right == NULL)		//Only one child again
    {
        node = oldevent->left;			//Child1 is attached to parent
        node->parent = parent;
    }
    else		  //Node to delete has 2 children
    {               //In this case: a) Find first node after oldevent     (This node will have no left)
                    //              b) Remove this node from the tree     (Attach node->right to node->parent)
                    //              c) Put this node in place of oldevent (Oldevent's children are adopted by node)
        node = oldevent->right;
        while (node->left) node = node->left;	//Find first node of right tree of descendants of oldevent
        particle* pnode = node->parent;
        if (pnode != oldevent)			//node is not a child of oldevent
        {						//Both of oldevent's children should be connected to node
            pnode->left = node->right;		//Remove node from right tree
            if (node->right) node->right->parent = pnode;
            oldevent->left->parent = node;
            node->left = oldevent->left;
            oldevent->right->parent = node;
            node->right = oldevent->right;
        }
        else					//This means node == oldevent->right
        {						//Only left has to be attached to node
            oldevent->left->parent = node;
            node->left = oldevent->left;
        }
        node->parent = parent;
    }
    if (parent->left == oldevent) parent->left = node;
    else                          parent->right = node;
}


/**************************************************
**                  OUTPUTSNAPSHOT
** Write a final snapshot to disk
**************************************************/
void outputsnapshot()
{
    char filename[200];
    sprintf(filename, "snapshot_end.sph");
    FILE* file = fopen(filename, "w");
    int i;
    particle *p, up2datep;
    double time = simtime + simtimewindowlength * timewindow ;

    fprintf(file, "%d %lf\n%.12g %.12g 0.0\n", (int)N, time, xsize, ysize);
    for (i = 0; i < N; i++)
    {
        p = particles + i;
	updatedparticle(p, &up2datep);

        fprintf(file, "%c %.12g %.12g %.12g %g\n", 'a' + p->type, up2datep.x + xsize * p->boxestraveledx, up2datep.y + ysize * p->boxestraveledy, 0.0, p->radius);
    }
    fclose(file);

}


/**************************************************
**                    DUMPSNAPSHOT
** Writes a movie
**************************************************/
void dumpsnapshot(particle* dumpevent)
{
  // static int first = 1;
    int i;
    particle *p, up2datep;
    FILE *file , *vel_file ;
    static int logstep = 1 ;
    double time = simtime + simtimewindowlength * timewindow ;

    char filename[200];
    sprintf(filename, "mov.n%d.v%.4g.sph", N, xsize * ysize);
    // if (first) { first = 0;  file = fopen(filename, "w"); }
    // else                     file = fopen(filename, "a");
    file = fopen(filename, "a");
    fprintf(file, "%d %g\n%.12g %.12g 0.0\n", (int)N, time, xsize, ysize);
    sprintf(filename, "vel.last.xyz") ;
    vel_file = fopen(filename, "w") ;
    fprintf(vel_file , "# N %d\n# timestep %.12g\n", (int)N, time) ;
    for (i = 0; i < N; i++)
      {
	p = &(particles[i]);
	updatedparticle(p, &up2datep);   //maybe not so efficient to compute 2 times the same quantities...

	fprintf(file, "%c %.12g %.12g %.12g %g\n", 
                'a' + p->type, 
                up2datep.x + xsize * p->boxestraveledx, 
                up2datep.y + ysize * p->boxestraveledy, 
                0.0, p->radius);
	fprintf(vel_file, "%.16g  %.16g\n", 
                up2datep.vx, 
                up2datep.vy);
      }
    fclose(file);
    fclose(vel_file);

    //Schedule next write event
    if(dumplogtime == 1) {
      createevent(pow(dumplogbase, logstep-1.0) * (dumplogbase - 1.0), dumpevent, NULL, 102) ;
      if( logstep == dump_logcyclelength ) logstep = 1 ;
      else logstep ++ ;
    } else createevent(snapshotinterval, dumpevent, NULL, 102) ;
}



/**************************************************
**                    WRITE
** Writes a movie
**************************************************/
void write(particle* writeevent)
{
    static int counter = 0 ;
    static double dptotlast[3] ;
    static double timelast = 0 ;   
    int i ;
    particle *p ;
    FILE *file ;

    double kinEn = 0;
    int maxneigh = 0, minneigh = 100;
    for (i = 0; i < N; i++)
    {
        p = particles + i;
        if (p->nneigh > maxneigh) maxneigh = p->nneigh;
        if (p->nneigh < minneigh) minneigh = p->nneigh;
    }
    computeKinenergy(&kinEn);
    double temperature = kinEn / (float)N ;

    double area = xsize * ysize;
    double pressid = (double)N / area;
    double time = simtime + simtimewindowlength * timewindow ;
    double expressnow[3] ;
    if(time!=timelast) {
      expressnow[XX] = -(dptot[XX] - dptotlast[XX]) / (area * (time - timelast));
      expressnow[YY] = -(dptot[YY] - dptotlast[YY]) / (area * (time - timelast));
      expressnow[XY] = -(dptot[XY] - dptotlast[XY]) / (area * (time - timelast));
    } else {
      expressnow[XX] = 0.0, expressnow[YY] = 0.0, expressnow[XY] = 0.0 ;
    }
    double pressnow = 0.5*(expressnow[XX] + expressnow[YY]) + pressid;
    dptotlast[XX] = dptot[XX] ;
    dptotlast[YY] = dptot[YY] ;
    dptotlast[XY] = dptot[XY] ;
    timelast = time;
    if (colcounter == 0) pressnow = 0;

    double listsize1 = (double)listcounter1 / mergecounter;     //Average number of events in the first event list
    int listsize2 = listcounter2;                               //Number of events in overflow list during last rescheduling (0 if not recently rescheduled)
    if (mergecounter == 0) listsize1 = 0;
    listcounter1 = listcounter2 = mergecounter = 0;

    printf("Simtime: %g, Collisions: %u, Press: %g, T: %g, TotEn: %g, Listsizes: (%lf, %d), Neigh: %d - %d\n", 
	   time, colcounter, pressnow, temperature, kinEn/N, listsize1, listsize2, minneigh, maxneigh);

    //Print some data to a file
    char filename[200];
    sprintf(filename, "press.n%d.v%.4g.sph", N, xsize * ysize);
    // if (counter == 0) file = fopen(filename, "w");
    // else              file = fopen(filename, "a");
    file = fopen(filename, "a");
    fprintf(file, "%g %.16g %.16g %.16g %.16g %.16g\n", time, pressnow, expressnow[XX], expressnow[YY], expressnow[XY], kinEn/N);
    fclose(file);

    counter++;

    //Schedule next write event
    createevent(writeinterval, writeevent, NULL, 100) ;
}



/**************************************************
**                    BACKINBOX
** Apply periodic boundaries
** Just for initialization !
**************************************************/
void backinbox(particle* p)
{
    p->x -= xsize * floor(p->x / xsize);
    p->y -= ysize * floor(p->y / ysize);
}


/**************************************************
**                    THERMOSTAT
** Simple thermostat
** Only used when "usethermostat" is non-zero
** Kicks a random selection of particles 
** periodically, giving them a new velocity from
** the Maxwell-Boltzmann distribution
**************************************************/
void thermostat(particle* thermostatevent)
{
    int i, num;
    particle* p;
    int freq = N / 100;
    if (freq == 0) freq = 1;
    for (i = 0; i < freq; i++)
    {
        num = genrand_real2() * N;			//Random particle
        p = particles + num;
        double imsq = 1.0 / sqrt(p->mass);
        update(p);
        p->vx = random_gaussian() * imsq;			//Kick it
        p->vy = random_gaussian() * imsq;
        p->counter ++;
        removeevent(p);
        findcollisions(p);
    }
    //Schedule next thermostat event
    createevent(thermostatinterval, thermostatevent, NULL, 101) ;
}


/**************************************************
**                    RANDOM_GAUSSIAN
** Generates a random number from a
** Gaussian distribution (Box–Muller)
**************************************************/
double random_gaussian()
{
    static int have_deviate = 0;
    static double u1, u2;
    double  x1, x2, w;

    if (have_deviate)
    {
        have_deviate = 0;
        return u2;
    }
    else
    {
        do
        {
            x1 = 2.0 * genrand_real2() - 1.0;
            x2 = 2.0 * genrand_real2() - 1.0;
            w = x1 * x1 + x2 * x2;
        } while (w >= 1.0);
        w = sqrt((-2.0 * log(w)) / w);
        u1 = x1 * w;
        u2 = x2 * w;
        have_deviate = 1;
        return u1;
    }
}




/******************************************************
**               SETPARAMETERSFROMFILE
** Sets simulation parameters read in from a file
******************************************************/
void setparametersfromfile( char * filename )
{
  char *buffer , **words ;
  int len = 1 , nwords = 0 , i ;
  FILE *paramfile ;
  buffer = (char*)calloc( len , sizeof(char) ) ;
  buffer[0] = '\0' ;

  if ( strcmp( filename , "" ) ) paramfile = fopen( filename , "r" ) ;
  else paramfile = fopen( "input.dat" , "r" ) ;
  if ( !paramfile ) printf("File not found, starting with default parameters\n");

  else {
    while ( myreadline( &buffer , &len , paramfile ) ) {
      words = mysplitline( &nwords , buffer , len ) ;
      if ( nwords > 0 ) {
	if( ! strcmp( words[0] , "area_fraction" ) ) sscanf( words[1] , "%lf" , &areafrac ) ;

	else if( ! strcmp( words[0] , "N_particles" ) ) sscanf( words[1] , "%d" , &N ) ;

	else if( ! strcmp( words[0] , "initial_configuration" ) ) {
	  if( ! strcmp( words[1] , "file" ) ) {
	    initialconfig = 0 ;
	    sprintf( inputfilename , "%s" , words[2] ) ;
	  } else if( ! strcmp( words[1] , "square" ) ) initialconfig = 1 ;
	  else if( ! strcmp( words[1] , "hexagonal" ) ) initialconfig = 2 ;
	  else if( ! strcmp( words[1] , "random" ) ) initialconfig = 3 ;
	  else if( ! strcmp( words[1] , "stampflirnd" ) ) initialconfig = 4 ;
	  else if( ! strcmp( words[1] , "addsmallparticles" ) ) {
	    initialconfig = 5 ;
	    sprintf( inputfilename , "%s" , words[2] ) ;
	  }
	}

	else if( ! strcmp( words[0] , "initial_velocities" ) ) {
	  if( ! strcmp( words[1] , "file" ) ) {
	    initialvelocities = 1 ;
	    sprintf( inputvelfile , "%s" , words[2] ) ;
	  }

	} else if( ! strcmp( words[0] , "size_ratio" ) ) sscanf( words[1] , "%lf" , &sizeratio ) ;

	else if( ! strcmp( words[0] , "large_spheres_fraction" ) ) sscanf( words[1] , "%lf" , &large2totalfraction ) ;

	else if( ! strcmp( words[0] , "nonadditivity" ) ) nonadditivity = -1 ;

	else if( ! strcmp( words[0] , "time" ) ) sscanf( words[1] , "%lf" , &maxtime ) ;

	else if( ! strcmp( words[0] , "snapshots" ) ) {
	  makesnapshots = 1 ;
	  if( ! strcmp( words[1] , "log" ) ) {
	    dumplogtime = 1 ;
	    sscanf( words[2] , "%lf" , &dumplogbase ) ;
	    sscanf( words[3] , "%d" , &dump_logcyclelength ) ;
	    if(dumplogbase <= 1) {
	      printf("*** dumplogbase and dump_logcyclelength must be strictly greater than 1 !\n");
	      exit(3);
	    }
	  } else sscanf( words[1] , "%lf" , &snapshotinterval ) ;

	} else if( ! strcmp( words[0] , "thermostat" ) ) {
	  usethermostat = 1 ;
	  sscanf( words[1] , "%lf" , &thermostatinterval ) ;

	} else if( ! strcmp( words[0] , "write_interval" ) ) sscanf( words[1] , "%lf" , &writeinterval ) ;

        else if( ! strcmp( words[0] , "seed" ) ) sscanf( words[1] , "%ld" , &seed ) ;

	for ( i=0 ; i<nwords ; i++ ) free(words[i]) ;
	nwords = 0 ;
	free(words) ;
      }
    }
  }
}




/******************************************************
**               SETATTRIBUTESFROMFILE
** Sets particles attributes reading options from a file
******************************************************/
void setattributesfromfile( char * filename )
{
  char *buffer , **words ;
  int len = 1 , nwords = 0 , i ;
  FILE *paramfile ;
  buffer = (char*)calloc( len , sizeof(char) ) ;
  buffer[0] = '\0' ;

  paramfile = fopen( filename , "r" ) ;
  if ( !paramfile ) printf("\n");

  else {
    while ( myreadline( &buffer , &len , paramfile ) ) {
      words = mysplitline( &nwords , buffer , len ) ;
      if ( nwords > 0 ) {
	if( ! strcmp( words[0] , "set" ) ) {
	  if( ! strcmp( words[1] , "type" ) ) {
	    char tmpc ;
	    sscanf( words[2] , "%c" , &tmpc ) ;
	    int type2change = tmpc - 'a' ;
	    if( ! strcmp( words[3] , "mass" ) ) {
	      double newmass = 1.0 ;
	      sscanf( words[4] , "%lf" , &newmass ) ;
	      printf( "Setting mass %lf for particles of type %c\n" , newmass , tmpc ) ;
	      for (i = 0; i < N; i++)  if( particles[i].type == type2change ) particles[i].mass = newmass ;
	    } else {
	      printf( "*** ERROR: set: unrecognized option %s\n" , words[3] ) ;
	      exit(3) ;
	    }
	  } else {
	    printf( "*** ERROR: set: unrecognized option %s\n" , words[1] ) ;
	    exit(3) ;
	  }
	}

	for ( i=0 ; i<nwords ; i++ ) free(words[i]) ;
	nwords = 0 ;
	free(words) ;
      }
    }
  }
}



/******************************************************
**               CHECKOVERLAP
** It looks for overlaps among a particle with the neighbors
** JUST BEFORE STARTING THE SIMULATION
******************************************************/
int checkoverlap(particle *p)
{
  double TOLERANCE = 1e-9 ;
  int overlap_found = 0 , j, i_ov = -1, j_ov = -1 ;
  double dx, dy, rmin, r ;
  particle *p2, up1, up2;

  updatedparticle(p, &up1);
  for (j = 0; j < p->nneigh; j++) {
    p2 = p->neighbors[j];
    updatedparticle(p2, &up2);
    dx = up1.x - up2.x;
    dy = up1.y - up2.y;
    if (dx > hx) dx -= xsize; else if (dx < -hx) dx += xsize;  //periodic boundaries
    if (dy > hy) dy -= ysize; else if (dy < -hy) dy += ysize;
    r = sqrt(dx*dx + dy*dy);
    rmin = p->radius + p2->radius;
    if (p->type != p2->type) rmin *= nonadditivity ;
    if (rmin - TOLERANCE > r) {
      i_ov = (int)(p-particles) ;
      j_ov = (int)(p2-particles) ;
      printf("*** OVERLAP FOUND: %g\t%lf ;\t%d\t%d\n", (rmin-r)/rmin, rmin, i_ov, j_ov);
      overlap_found = 1 ;
    }
  }
  return overlap_found ;
}


/******************************************************
**               CHECKOVERLAPS
** It looks for overlapping particles
** JUST BEFORE STARTING THE SIMULATION
******************************************************/
void checkoverlaps()
{
  int overlap_found = 0 , i = 0 ;
  particle *p1;

  while( i < N ) {
    p1 = particles + i;
    overlap_found += checkoverlap(p1);
    i ++ ;
  }

  if (overlap_found > 0) {
    printf("*** Found %d overlaps in the current configuration\n", overlap_found/2) ;
    exit(3) ;
  }
}

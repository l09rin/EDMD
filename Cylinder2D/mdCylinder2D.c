#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "mt19937ar.c"
unsigned long seed = 1;     //Seed for random number generator

//Maximum number of neighbors per particle
#define MAXNEIGH 20

#include "mdCylinder2D.h"
#include "mystring.c"

//Number of extra events (e.g. write, thermostat) to allocate space for
#define EXTRAEVENTS 12

//Pi (if not already defined)
#ifndef M_PI
#define M_PI 3.1415926535897932
#endif


double maxtime = 20000;           //Simulation stops at this time
int makesnapshots = 1;          //Whether to make snapshots during the run (yes = 1, no = 0)
double writeinterval = 1;     //Time between output to screen / data file
double snapshotinterval = 50;  //Time between snapshots (should be a multiple of writeinterval)

int initialconfig = 2;    //= 0 load from file, 1 = FCC crystal, 2 = random bi-disperse quasi-2D layer of equal mass density particles
char inputfilename[100] = "init.sph"; //File to read as input snapshot (for initialconfig = 0)
double packfrac = 0.2;                     //Packing fraction (for initialconfig = 1)
int N = 5000;             //Number of particles (for FCC)
double sizeratio = 0.54;  // ratio among the diameters of small and large particles (for bi-disperse mixture)
double large2totalfraction = 0.65;    // fraction of large particles (for bi-disperse mixture)
double areafrac = 0.83;               // naive area fraction of the sedimented system (for bi-disperse mixture)

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
double maxscheduletime = 0.1;
int numeventlists;
double eventlisttimemultiplier = 0.05;  //event list time will be this / N
double eventlisttime;


//Neighbor lists
const double shellsize = 1.5; //Shell size (equals 1+ \alpha)


//Internal variables
double simtime = 0;
double reftime = 0;
int currentlist = 0;
int totalevents;

int listcounter1 = 0, listcounter2 = 0, mergecounter = 0;

particle** eventlists; //Last one is overflow list

particle* particles;
particle** celllist;
particle* root;
double xsize, ysize, zsize; //Box size
double hx, hy, hz; //Half box size
double icxsize, icysize; //Inverse Cell size
int    cx, cy;  //Number of cells
double dvtot = 0;   //Momentum transfer (for calculating pressure)
unsigned int colcounter = 0; //Collision counter (will probably overflow in a long run...)


const int usethermostat = 1; //Whether to use a thermostat
double thermostatinterval = 0.01;

// with infinite vertical cells walls have to be present, to correctly predict collisions (with the correct image)
// hard walls on top and bottom of the simulation box along the z-direction
#define WALLCOLLISION_TYPE 16
unsigned int topZwallcolcounter = 0 , btmZwallcolcounter = 0 ;
double topZwalldvtot = 0 , btmZwalldvtot = 0 ;   //Momentum transfer to the wall

double g = 264.0 ; //constant acceleration along the z-axis, only positive values, oriented towards z negative


int main( int argc, char **argv )
{
    if (argc>1) setparametersfromfile(argv[1]) ;
    else setparametersfromfile("") ;
    init();
    printf("Starting\n");

    while (simtime <= maxtime)
    {
      step();
    }
    simtime = maxtime;

    printstuff();
    outputsnapshot();

    free(celllist);
    free(particles);
    free(eventlists);
    return 0;
}

/**************************************************
**                 COMPUTEENERGY
** Computes kinetic and potential energies at simulation time
**************************************************/
void computeenergy(double *kinEn, double *potEn)
{
    int i;
    particle *p, up2datep;
    double v2tot = 0;
    *kinEn = 0.0 ;
    *potEn = 0.0 ;

    for (i = 0; i < N; i++)
    {
        p = particles + i;
	updatedparticle(p, &up2datep);
        v2tot += p->mass * (up2datep.vx * up2datep.vx + up2datep.vy * up2datep.vy + up2datep.vz * up2datep.vz);
        *potEn += p->mass * g * up2datep.z ;
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
    double potEn = 0.0 , kinEn = 0.0 ;

    for (i = 0; i < N; i++)
    {
        p = particles + i;
        vfilled += p->radius * p->radius * p->radius * 8;
    }
    vfilled *= M_PI / 6.0;
    computeenergy(&kinEn, &potEn) ;
    printf("Average kinetic energy: %lf\n", kinEn / N);
    printf("Average potential energy: %lf ,\tgravity: %lf\n", potEn / N, g);
    printf("Total energy: %lf\n", (potEn + kinEn) / N);
    double volume = xsize * ysize * zsize;
    double dens = N / volume;
    double press = -dvtot / (3.0 * volume * simtime);
    double pressid = dens;
    double presstot = press + pressid;
    printf("Total time simulated  : %lf\n", simtime);
    //  printf ("Density               : %lf\n", (double) N / volume);
    printf("Packing fraction      : %lf\n", vfilled / volume);
    printf("Measured pressure     : %lf + %lf = %lf\n", press, pressid, presstot);
    double topZwallpress = -topZwalldvtot / (xsize * ysize * simtime) , btmZwallpress = btmZwalldvtot / (xsize * ysize * simtime) ;   //pressure acting on the wall
    printf("    pressure on the top wall     : %lf\n", topZwallpress);
    printf("    pressure on the bottom wall     : %lf\n", btmZwallpress);

}


/**************************************************
**                    INIT
 **************************************************/
void init()
{
    int i;
    //   FILE *fp=fopen("/dev/urandom","r");
    //   int tmp = fread(&seed,1,sizeof(unsigned long),fp);
    //   if (tmp != sizeof(unsigned long)) printf ("error with seed\n");
    //   fclose(fp);
    printf("Seed: %u\n", (int)seed);
    init_genrand(seed);
    
    if (g < 0) {
      printf("** The gravity acceleration is considered as negative, downward\n");
      g *= -1.0 ;
    }

    if (initialconfig == 0)
    {
      loadparticles();
      checkifinsideZwalls();
    } else if (initialconfig == 1)           fcc();
    else if (initialconfig == 2)             binary2Dlayer();
    randommovement();
    hx = 0.5 * xsize; hy = 0.5 * ysize; hz = 0.5 * zsize;	//Values used for periodic boundary conditions

    initevents();
    for (i = 0; i < N; i++)
    {
        particle* p = particles + i;
        p->boxestraveledx = 0;
        p->boxestraveledy = 0;
        p->boxestraveledz = 0;
        p->nneigh = 0;
        p->counter = 0;
        p->t = 0;
        p->xn = p->x;   //Set center of neighbor list to current position
        p->yn = p->y;
        p->zn = p->z;
    }
    initcelllist();

    for (i = 0; i < N; i++)
    {
        makeneighborlist(particles + i);
    }
    printf("Done adding collisions\n");



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
**                    FCC
** Initialize system on a face-centered cubic
** lattice
**************************************************/
void fcc()
{
    int i, j, k;
    particle* p;

    int ncell = cbrt(N / 4) + 0.0001;

    if (ncell * ncell * ncell * 4 != N)
    {
        printf("N should be 4 * a perfect cube! (e.g. %d)\n", ncell * ncell * ncell * 4);
        exit(3);
    } else if (packfrac > M_PI/(3*sqrt(2.)))
    {
        printf("The packing fraction of an FCC lattice of HS cannot be greater than ~0.74!\n");
        exit(3);
    }

    double volume = N / (6.0 / M_PI * packfrac);
    xsize = cbrt(volume);
    ysize = xsize;
    zsize = xsize;

    double step = xsize / ncell;

    double zoffset = 0.0 ;
    if ( step < 2.0 )  //max packing fraction for a cubic ensemble of FCC cells
                                      //with no particles overlapping the walls phi_max = M_PI / 12 < 0.2618 
    {                                 //particles are intended to have radius 0.5
      printf("*** At such a packing fraction an fcc() crystal enclosed between walls cannot be arranged in a cubic shape !\n") ;
      // To safely enclose the crystal within z-walls we roughly rescale zsize,xzize,ysize to keep the same packing fraction
      //   and a tiny non-zero distance of the first and last layers of particles from the wall
      //   (used the cardano formula to solve the cubic eqn for lambda: rescaled i-sizes as follows, such that the distance
      //    between the wall and the particles is dz = 0.01*part->radius -- here = (0.5) )
      double a = 2*1.01*(0.5)/ncell/step , b = 1-0.5/ncell ;
      double delta = b*b/4-a*a*a/27 ;
      double lambda = cbrt( 0.5*b+sqrt(delta) ) + cbrt( 0.5*b-sqrt(delta) ) ;
      zsize = lambda*lambda*zsize ;
      xsize = xsize/lambda ;
      ysize = ysize/lambda ;
      step = step/lambda ;
      zoffset = 1.01*(0.5)-0.25*step ;
      printf("Box size (x,y,z): %.3lf, %.3lf, %.3lf\n", xsize, ysize, zsize) ;
    }

    printf("step: %lf\n", step);
    initparticles(N);
    printf("Placing particles\n");

    p = particles;
    for (i = 0; i < ncell; i++) for (j = 0; j < ncell; j++) for (k = 0; k < ncell; k++)
    {
        p->x = (i + 0.25) * step;
        p->y = (j + 0.25) * step;
        p->z = (k + 0.25) * step + zoffset;
        p++;
        p->x = (i + 0.75) * step;
        p->y = (j + 0.75) * step;
        p->z = (k + 0.25) * step + zoffset;
        p++;
        p->x = (i + 0.75) * step;
        p->y = (j + 0.25) * step;
        p->z = (k + 0.75) * step + zoffset;
        p++;
        p->x = (i + 0.25) * step;
        p->y = (j + 0.75) * step;
        p->z = (k + 0.75) * step + zoffset;
        p++;
    }

    for (i = 0; i < N; i++)
    {
        p = particles + i;
        p->radius = 0.5;
        p->type = 0;
        p->mass = 1;
    }

    printf("Packing fraction: %lf\n", M_PI / (6.0 * xsize * ysize * zsize) * N);
    printf("Starting configuration from fcc crystal\n");
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
** - The radius
**************************************************/
void loadparticles()
{
    char tmp;
    int i, npart;
    particle* p;
    char buffer[255];

    FILE* file;
    file = fopen(inputfilename, "r");
    if (!file)
    {
        printf("File not found!\n");
        exit(3);
    }
    mygetline(buffer, file);
    int ftmp = sscanf(buffer, "%d", &npart);
    if (ftmp != 1) { printf("Read error (n or box)\n"); exit(3); }
    mygetline(buffer, file);
    ftmp = sscanf(buffer, "%lf %lf %lf\n", &xsize, &ysize, &zsize);
    if (ftmp != 3) { printf("Read error (n or box)\n"); exit(3); }


    initparticles(npart);
    printf("Placing particles\n");
    double vfilled = 0;
    for (i = 0; i < N; i++)
    {
        p = &(particles[i]);
        mygetline(buffer, file);
        ftmp = sscanf(buffer, "%c %lf  %lf  %lf %lf\n", &tmp, &(p->x), &(p->y), &(p->z), &(p->radius));
        backinbox(p);
        if (ftmp != 5) { printf("Read error (particle) %d \n String: %s\n", ftmp, buffer); exit(3); }
        p->type = tmp - 'a';
        p->mass = 1;
        vfilled += 8*p->radius*p->radius*p->radius;
    }
    fclose(file);

    printf("Packing fraction: %lf\n", M_PI / (6.0 * xsize * ysize * zsize) * vfilled);
    printf("Starting configuration read from %s\n", inputfilename);
}


/**************************************************
**                    BINARY2DLAYER
** Randomly spreads out particles of 2 sizes 
** within a thin xy-layer
**************************************************/
void binary2Dlayer()
{
    particle* p ;
    int i, Nlarge = N * large2totalfraction ;
    large2totalfraction = (double)Nlarge / N ;

    initparticles(N);
    for (i = 0; i < Nlarge; i++) {
        p = particles + i ;
        p->radius = 0.5 ;
        p->type = 0 ;
        p->mass = 1.0 ;
    }
    for (i = Nlarge; i < N; i++) {
        p = particles + i ;
        p->radius = 0.5 * sizeratio ;
        p->type = 1 ;
        p->mass = 8. * p->radius * p->radius * p->radius ;
    }

    //given the area fraction the x and y box sides are calculated
    xsize = sqrt( M_PI / 4 * N / areafrac * (large2totalfraction + sizeratio*sizeratio*(1. - large2totalfraction)) ) ;
    ysize = xsize ;
    hx = hy = 0.5 * xsize ;
    //zsize is initially given a value to conveniently randomly spread all the particles near to the bottom of the box
    packfrac = 0.25 ;
    zsize = M_PI / 6 * N / packfrac / xsize / ysize * (large2totalfraction + sizeratio*sizeratio*sizeratio*(1. - large2totalfraction)) ;
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
    double dx, dy, dz, r2, rm ;
    particle *p2 ;
    printf("Placing particles\n");
    i = 0 ;
    while ( i < N ) {
        p = particles + i ;
	p->x = genrand_real2() * xsize ;
	p->y = genrand_real2() * ysize ;
	p->z = genrand_real2() * (zsize - 2.*p->radius) + p->radius ;
	cellx = p->x * icxsize + cx ;
	celly = p->y * icysize + cy ;
	overlap = 0 ;
	for (cdx = cellx - 1; cdx < cellx + 2; cdx++)
	  for (cdy = celly - 1; cdy < celly + 2; cdy++) {
            p2 = celllist[celloffset(cdx % cx, cdy % cy)];
            while (p2) {
	      dx = p->x - p2->x ;
	      dy = p->y - p2->y ;
	      dz = p->z - p2->z ;
	      if (p2->nearboxedge) {
		if (dx > hx) dx -= xsize; else if (dx < -hx) dx += xsize;  //periodic boundaries
		if (dy > hy) dy -= ysize; else if (dy < -hy) dy += ysize;
	      }
	      r2 = dx * dx + dy * dy + dz * dz;
	      rm = (p->radius + p2->radius) ;
	      if (r2 < rm * rm) {
		overlap = 1 ;
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

    zsize = xsize ;
    printf("Starting configuration from randomly bidisperse HS in a quasi-2D layer\n") ;
    printf("Area fraction: %lf\n", areafrac) ;
    printf("Size ratio: %lf\n", sizeratio) ;
    printf("Fraction of large particles (R_L = 0.5): %lf\n", large2totalfraction) ;
}




/**************************************************
**                RANDOMMOVEMENT
** Assign random velocities to all particles
** Center-of-mass velocity = 0
** Kinetic energy per particle = 3kT/2
**************************************************/
void randommovement()
{
    particle* p;
    double v2tot = 0, vxtot = 0, vytot = 0, vztot = 0;
    double mtot = 0;
    int i;

    for (i = 0; i < N; i++)
    {
        p = particles + i;
        double imsq = 1.0 / sqrt(p->mass);

        p->vx = imsq * random_gaussian();
        p->vy = imsq * random_gaussian();
        p->vz = imsq * random_gaussian();
        vxtot += p->mass * p->vx;					//Keep track of total v
        vytot += p->mass * p->vy;
        vztot += p->mass * p->vz;
        mtot += p->mass;
    }


    vxtot /= mtot; vytot /= mtot; vztot /= mtot;
    for (i = 0; i < N; i++)
    {
        p = &(particles[i]);
        p->vx -= vxtot;					//Make sure v_cm = 0
        p->vy -= vytot;
        p->vz -= vztot;
        v2tot += p->mass * (p->vx * p->vx + p->vy * p->vy + p->vz * p->vz);
    }
    double fac = sqrt(3.0 / (v2tot / N));
    for (i = 0; i < N; i++)
    {
        p = &(particles[i]);
        p->vx *= fac;					//Fix energy
        p->vy *= fac;
        p->vz *= fac;
    }
}

/**************************************************
**                UPDATE
** Update particle to the current time
**************************************************/
void update(particle* p1)
{
    double dt = simtime - p1->t;
    p1->t = simtime;
    p1->x += dt * p1->vx;
    p1->y += dt * p1->vy;
    p1->z += dt * p1->vz - 0.5 * g * dt*dt;
    p1->vz -= g * dt;
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
    p2->z = p1->z;
    p2->vx = p1->vx;
    p2->vy = p1->vy;
    p2->vz = p1->vz;
    if (p1->t < simtime) {
      double dt = simtime - p1->t;
      p2->x += dt * p1->vx;
      p2->y += dt * p1->vy;
      p2->z += dt * p1->vz - 0.5 * g * dt*dt;
      p2->vz -= g * dt;
    }
    p2->t = simtime;
}

/**************************************************
**                 INITCELLLIST
** Initialize the cell list
** Cell size should be at least:
**    shellsize * [diameter of largest sphere]
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
    cx = (int)(xsize - 0.0001) / shellsize / maxdiameter;				//Set number of cells
    cy = (int)(ysize - 0.0001) / shellsize / maxdiameter;

    while (cx*cy > 8*N)
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
        addtocelllist(p, p->x * icxsize, p->y * icysize);
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

    simtime = ev->eventtime;
    removeevent(ev);
    switch(ev->eventtype)
    {
        case 0:
            collision(ev);
            break;
        case 8:
            makeneighborlist(ev);
            break;
        case WALLCOLLISION_TYPE:
            zwallcollision(ev);
            break;
        case 100:
            write(ev);
            break;
        case 101:
            thermostat(ev);
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
    if (p1->z >= zsize) { p1->z -= zsize; p1->boxestraveledz++; }
    else if (p1->z < 0) { p1->z += zsize; p1->boxestraveledz--; }

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
    double dr2 = dx * dx + dy * dy;
    double md = (shellsize - 1) * p1->radius;
    double b = dx * dvx + dy * dvy ;                  //drh.dvh
    double disc = b * b - dv2 * (dr2 - md * md) ;

    return  (-b + sqrt(disc)) / dv2;     //time to go out from the lateral surface
}

/**************************************************
**                FINDCOLLISION
** Detect the next collision for two particles
** Note that p1 is always up to date in
** findcollision
**************************************************/
int findcollision(particle* p1, particle* p2, double* tmin)
{
    particle p2updated;
    updatedparticle(p2, &p2updated);
    double dx = p1->x - p2updated.x;    //relative distance at current time
    double dy = p1->y - p2updated.y;
    double dz = p1->z - p2updated.z;
    if (p1->nearboxedge)
    {
        if (dx > hx) dx -= xsize; else if (dx < -hx) dx += xsize;  //periodic boundaries
        if (dy > hy) dy -= ysize; else if (dy < -hy) dy += ysize;
	//        if (dz > hz) dz -= zsize; else if (dz < -hz) dz += zsize;
    }

    double dvx = p1->vx - p2updated.vx;                               //relative velocity
    double dvy = p1->vy - p2updated.vy;
    double dvz = p1->vz - p2updated.vz;

    double b = dx * dvx + dy * dvy + dz * dvz;                  //dr.dv
    if (b > 0) return 0;                                      //Particles flying apart
    double dr2 = dx * dx + dy * dy + dz * dz;
    double md = p1->radius + p2->radius;
    double A = md * md - dr2;
    if (2 * b * *tmin > A) return 0;                        //Approximate check to discard hits at times > tmin

    double dv2 = dvx * dvx + dvy * dvy + dvz * dvz;

    double disc = b * b + dv2 * A;
    if (disc < 0) return 0;
    double t = (-b - sqrt(disc)) / dv2;
    if (t < *tmin) 
    {
        *tmin = t;
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
    if( findZwallscollision(p1, &tmin) ) type = WALLCOLLISION_TYPE ;
    for (i = 0; i < p1->nneigh; i++)
    {
        p2 = p1->neighbors[i];
        if(findcollision(p1, p2, &tmin))
        {
            partner = p2;
            type = 0;
        }
    }
    createevent(tmin + simtime, p1, partner, type);
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
        double tmin = findneighborlistupdate(p1);
        int type = 8;
	if( findZwallscollision(p1, &tmin) ) type = WALLCOLLISION_TYPE ;
        for (j = 0; j < p1->nneigh; j++)
        {
            particle* p2 = p1->neighbors[j];
            if (p2 > p1)
            {
                if(findcollision(p1, p2, &tmin))
                {
                    partner = p2;
                    type = 0;
                }
            }
        }
        createevent(tmin, p1, partner, type);
        p1->counter2 = partner->counter;
    }
}







/**************************************************
**                  COLLISION
** Process a single collision event
**************************************************/
void collision(particle* p1)
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
    double rinv = 1.0 / r;
    double dx = (p1->x - p2->x);			//Normalized distance vector
    double dy = (p1->y - p2->y);
    double dz = (p1->z - p2->z);
    if (p1->nearboxedge)
    {
        if (dx > hx) dx -= xsize; else if (dx < -hx) dx += xsize;  //periodic boundaries
        if (dy > hy) dy -= ysize; else if (dy < -hy) dy += ysize;
	//        if (dz > hz) dz -= zsize; else if (dz < -hz) dz += zsize;
    }
    dx *= rinv;  dy *= rinv;  dz *= rinv;

    double dvx = p1->vx - p2->vx;                               //relative velocity
    double dvy = p1->vy - p2->vy;
    double dvz = p1->vz - p2->vz;

    double b = dx * dvx + dy * dvy + dz * dvz;                  //dr.dv
    b *= 2.0 / (m1 + m2);
    double dv1 = b * m2, dv2 = b * m1;

    p1->vx -= dv1 * dx;         //Change velocities after collision
    p1->vy -= dv1 * dy;         //delta v = (-) dx2.dv2
    p1->vz -= dv1 * dz;
    p2->vx += dv2 * dx;
    p2->vy += dv2 * dy;
    p2->vz += dv2 * dz;

    dvtot += dv1*m1 * r;
    colcounter++;

    removeevent(p2);


    findcollisions(p1);
    findcollisions(p2);
}




/**************************************************
**                  ZWALLCOLLISION
** Process a single collision event with the wall
**************************************************/
void zwallcollision(particle* p)
{
    update(p);
    p->counter++;

    p->vz *= -1.0;         //Change velocities after collision

    if (p->vz > 0)
    {
      btmZwalldvtot += 2*p->vz*p->mass;
      btmZwallcolcounter++;

    } else {
      topZwalldvtot += 2*p->vz*p->mass;
      topZwallcolcounter++;
    }

    findcollisions(p);
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
    root->eventtype = 127;					//This makes sure we don't have to keep checking this when adding/removing events
    root->parent = NULL;

    particle* writeevent = particles + N + 1;		//Set up write event
    writeevent->eventtime = 0;
    writeevent->eventtype = 100;
    writeevent->p2 = NULL;
    addevent(writeevent);


    printf("Event tree initialized.\n");

    if (usethermostat)
    {
        particle* thermostatevent = particles + N + 2;
        thermostatevent->eventtime = thermostatinterval;
        thermostatevent->eventtype = 101;
        thermostatevent->p2 = NULL;
        addevent(thermostatevent);

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
    particle* loc = root;
    int busy = 1;
    while (busy)						//Find location to add event into tree (loc)
    {
        if (time < loc->eventtime)				//Go left
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
    double dt = newevent->eventtime - reftime;

    if (dt < eventlisttime) //Event in near future: Put it in the tree
    {
        newevent->queue = currentlist;
        addeventtotree(newevent);
    }
    else                    //Put it in one of the event lists
    {
        int list_id = currentlist + dt / eventlisttime;
        if (list_id >= numeventlists)
        {
            list_id -= numeventlists;
            if (list_id > currentlist - 1) list_id = numeventlists; //Overflow
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
**************************************************/
void createevent(double time, particle* p1, particle* p2, int type)
{
    p1->eventtime = time;
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

    fprintf(file, "%d\n%.12lf %.12lf %.12lf\n", (int)N, xsize, ysize, zsize);
    for (i = 0; i < N; i++)
    {
        p = particles + i;
	updatedparticle(p, &up2datep);

        fprintf(file, "%c %.12lf  %.12lf  %.12lf  %lf\n", 'a' + p->type, up2datep.x + xsize * p->boxestraveledx, up2datep.y + ysize * p->boxestraveledy, up2datep.z + zsize * p->boxestraveledz, p->radius);
    }
    fclose(file);

}


/**************************************************
**                    WRITE
** Writes a movie
**************************************************/
void write(particle* writeevent)
{
    static int counter = 0;
    static int first = 1;
    static double lastsnapshottime = -999999999.9;
    static double dvtotlast = 0;
    static double btmZwalldvtotlast = 0 , topZwalldvtotlast = 0 ;
    static double timelast = 0;   
    int i;
    particle *p, up2datep;
    FILE* file;

    double kinEn = 0, potEn = 0, totEn = 0;
    int maxneigh = 0, minneigh = 100;
    for (i = 0; i < N; i++)
    {
        p = particles + i;
        if (p->nneigh > maxneigh) maxneigh = p->nneigh;
        if (p->nneigh < minneigh) minneigh = p->nneigh;
    }
    computeenergy(&kinEn, &potEn);
    totEn = kinEn + potEn ;
    double temperature = kinEn / (float)N / 1.5;

    double volume = xsize * ysize * zsize;
    double pressid = (double)N / volume;
    double pressnow = -(dvtot - dvtotlast) / (3.0 * volume * (simtime - timelast));
    pressnow += pressid;
    dvtotlast = dvtot;
    double area = xsize * ysize ;
    double btmpressnow = (btmZwalldvtot - btmZwalldvtotlast) / (area * (simtime - timelast));
    double toppressnow = -(topZwalldvtot - topZwalldvtotlast) / (area * (simtime - timelast));
    btmZwalldvtotlast = btmZwalldvtot ;
    topZwalldvtotlast = topZwalldvtot ;
    timelast = simtime;
    if (colcounter == 0) pressnow = 0;
    if (btmZwallcolcounter == 0) btmpressnow = 0 ;
    if (topZwallcolcounter == 0) toppressnow = 0 ;

    double listsize1 = (double)listcounter1 / mergecounter;     //Average number of events in the first event list
    int listsize2 = listcounter2;                               //Number of events in overflow list during last rescheduling (0 if not recently rescheduled)
    if (mergecounter == 0) listsize1 = 0;
    listcounter1 = listcounter2 = mergecounter = 0;

    printf("Simtime: %lf, Collisions: %u, Press: %lf, T: %lf, PotEn: %lf, TotEn: %lf, Listsizes: (%lf, %d), Neigh: %d - %d\n", 
	   simtime, colcounter, pressnow, temperature, potEn/N, totEn/N, listsize1, listsize2, minneigh, maxneigh);

    char filename[200];
    if (makesnapshots && simtime - lastsnapshottime > snapshotinterval - 0.001)
    {
        sprintf(filename, "mov.n%d.v%.4lf.sph", N, xsize * ysize * zsize);
        if (first) { first = 0;  file = fopen(filename, "w"); }
        else                     file = fopen(filename, "a");
        fprintf(file, "%d\n%.12lf %.12lf %.12lf\n", (int)N, xsize, ysize, zsize);
        for (i = 0; i < N; i++)
        {
            p = &(particles[i]);
	    updatedparticle(p, &up2datep);   //maybe not so efficient to compute 2 times the same quantities...

            fprintf(file, "%c %.12lf  %.12lf  %.12lf  %lf\n", 
                'a' + p->type, 
                up2datep.x + xsize * p->boxestraveledx, 
                up2datep.y + ysize * p->boxestraveledy, 
                up2datep.z + zsize * p->boxestraveledz, p->radius);
        }
        fclose(file);
        lastsnapshottime = simtime;
    }

    //Print some data to a file
    sprintf(filename, "press.n%d.v%.4lf.sph", N, xsize * ysize * zsize);
    if (counter == 0) file = fopen(filename, "w");
    else              file = fopen(filename, "a");
    fprintf(file, "%lf %lf %lf %lf %lf %lf\n", simtime, pressnow, btmpressnow, toppressnow, potEn/N, totEn/N);
    fclose(file);

    counter++;

    //Schedule next write event    
    writeevent->eventtime = simtime + writeinterval;
    addevent(writeevent);
}



/**************************************************
**                    BACKINBOX
** Apply periodic boundaries
** Just for initialization
**************************************************/
void backinbox(particle* p)
{
    p->x -= xsize * floor(p->x / xsize);
    p->y -= ysize * floor(p->y / ysize);
    p->z -= zsize * floor(p->z / zsize);
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
        p->vz = random_gaussian() * imsq;
        p->counter ++;
        removeevent(p);
        findcollisions(p);
    }
    //Schedule next thermostat event
    thermostatevent->eventtime = simtime + thermostatinterval;
    addevent(thermostatevent);
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
**               CHECKIFINSIDEZWALLS
** Assumes that the system is enclosed between hard
**   walls along the z-direction
** It checks if superposition between particles and the
**   walls can be found
******************************************************/
void checkifinsideZwalls()
{
  int overlap_found = 0 , i = 0 , i_ov = -1 ;
  double z_ov = 0 ;
  particle* p ;

  while( !overlap_found && i < N )
  {
    p = particles + i ;
    backinbox(p) ;
    if ( p->z < p->radius || (p->z+p->radius) > zsize ) {
      overlap_found = 1 ;
      i_ov = i ;
      z_ov = p->z ;
    }
    i++ ;
  }
  if (overlap_found)
  {
    printf("At least one particle in the configuration overlaps the walls (idx, z): %d , %lf", i_ov, z_ov) ;
    exit(3) ;
  }
}


/******************************************************
**               FINDZWALLCOLLISION
** Detect the next collision between a particle and the
**   walls in the z-direction
** Assumes that p1 is up to date and within the box
******************************************************/
int findZwallscollision(particle* p, double* tmin)
{
  double t = -1 ;

  if (g==0)
  {
    if (p->vz<0)  t = (p->radius - p->z) / p->vz ;
    else if (p->vz>0)  t = (zsize - p->z - p->radius) / p->vz ;
    if (t < *tmin) {
      *tmin = t;
      return 1;
    }

  } else {
    double disc ;
    if (p->vz > 0) {
      disc = p->vz * p->vz + 2*g*(p->z-zsize+p->radius) ;
      if (disc>=0)  t = (p->vz - sqrt(disc) ) / g ;
    }
    if (t < 0) {
      disc = p->vz * p->vz + 2*g*(p->z-p->radius) ;
      t = ( p->vz + sqrt(disc) ) / g ;
    }
    if (t < *tmin) {
      *tmin = t;
      return 1;
    }
  }

  return 0 ;
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

	else if( ! strcmp( words[0] , "packing_fraction" ) ) sscanf( words[1] , "%lf" , &packfrac ) ;

	else if( ! strcmp( words[0] , "N_particles" ) ) sscanf( words[1] , "%d" , &N ) ;

	else if( ! strcmp( words[0] , "initial_configuration" ) ) {
	  if( ! strcmp( words[1] , "file" ) ) {
	    initialconfig = 0 ;
	    sprintf( inputfilename , "%s" , words[2] ) ;
	  } else if( ! strcmp( words[1] , "fcc" ) || ! strcmp( words[1] , "FCC" ) ) initialconfig = 1 ;
	  else if( ! strcmp( words[1] , "random_bidisperse" ) ) initialconfig = 2 ;
	}

	else if( ! strcmp( words[0] , "size_ratio" ) ) sscanf( words[1] , "%lf" , &sizeratio ) ;

	else if( ! strcmp( words[0] , "large_spheres_fraction" ) ) sscanf( words[1] , "%lf" , &large2totalfraction ) ;

	else if( ! strcmp( words[0] , "gravity" ) ) sscanf( words[1] , "%lf" , &g ) ;

	else if( ! strcmp( words[0] , "time" ) ) sscanf( words[1] , "%lf" , &maxtime ) ;

	else if( ! strcmp( words[0] , "snapshots" ) ) {
	  makesnapshots = 1 ;
	  sscanf( words[1] , "%lf" , &snapshotinterval ) ;

	} else if( ! strcmp( words[0] , "write_interval" ) ) sscanf( words[1] , "%lf" , &writeinterval ) ;

	else if( ! strcmp( words[0] , "seed" ) ) sscanf( words[1] , "%ld" , &seed ) ;

	for ( i=0 ; i<nwords ; i++ ) free(words[i]) ;
	nwords = 0 ;
	free(words) ;
      }
    }
  }
}

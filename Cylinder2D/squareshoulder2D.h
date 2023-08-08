
//Particle structure
typedef struct sparticle
{
	double x, y;				//Position
	double vx, vy;			//Velocity
	double xn, yn;			//Neighbor list center
	struct sparticle* neighbors[MAXNEIGH];
	uint8_t nneigh;				//Number of neighbors
        double t;					//Last update time
        int timewindow ;
	double radius;
	double mass;
        uint8_t nearboxedge;		//Is this particle in a cell near the box edge?
	int cell;					//Current cell
	int boxestraveledx, boxestraveledy;	//Keeps track of dynamics across periodic boundaries
	unsigned int counter;		//Number of collision events experienced
	unsigned int counter2;		//Number of collision events of collision partner at time of scheduling
	struct sparticle* prev, *next;	//Doubly linked cell list
	uint8_t type;				//Particle type
	double eventtime;
	int eventtimewindow ;
	struct sparticle* left;		//Left child in tree or previous event in event list
	struct sparticle* right;	//Right child in tree or next event in event list
	struct sparticle* parent;	//Parent in tree
	struct sparticle* p2;		//Collision partner
	int queue;					//Index of the event queue the event is in
	unsigned char eventtype;
} particle;

int main();
void printstuff();
void init();
void setparametersfromfile( char * filename );


void initevents();
void squarelattice();
void loadparticles();
void hexagonal();
void randommovement();
void loadvelocities();
void initcelllist();
void addtocelllist(particle* p, int cellx, int celly);
void removefromcelllist(particle* p);
int celloffset(int a, int b);

void step();
int findcollision(particle*, particle*, double*);
void findallcollisions();
void findcollisions(particle*);
void collision(particle*);
void update(particle *);
inline void updatedparticle(particle *, particle *);

void addevent(particle*);
void removeevent(particle*);
void createevent(double time, particle* p1, particle* p2, int type);
void addnexteventlist();
double findneighborlistupdate(particle* p1);
void makeneighborlist(particle* p1);

void outputsnapshot();
void dumpsnapshot(particle* ev);
void write(particle* ev);
void thermostat(particle* ev);
double random_gaussian();
void backinbox(particle* p);

// 2 eliminate
int findZwallscollision(particle* p, double* tmin);
void zwallcollision(particle* p);
void ramp(particle* ev);

void computeenergy(double *, double *);

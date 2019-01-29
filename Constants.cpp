double cellRadius = 0.34;		// radius of each cell cap (constant for now, microns)
double L_divide = 5.75;		// 3 length when cells divide (microns)
double k_cc = 100000.0;			// elastic constant for cell-cell interactions (atm?)
double k_wc = 100000.0;			// elastic constant between cells and wall
double varL = 0.0;			// variation in the length of daughter cells
double varAngle = 0.0;
double var_pos = 0.0;
double viscosity = 1.0;		// dimensionless viscosity of the medium: vis*growth_rate/p_thresh
double wall_rough = 0.0;
double gamma_t = 10000.0;
double wall_mu = 0.15;
double cell_mu = 0.15;
double density_threshold = 0.6;   // Mya set this value to 0.6 previously.
double tension = 10.0;
double DH = 0.0;			// determines how tightly the surface tension holds the cells (+ value is low agar concentration, - value is high agar concentration)

// time
double t_max = 20.0; // total simulation time
double initial_dt = 0.00001;	// initial time step
double OutputTime = 100*initial_dt;		// how often to output
double t0 = 0.0;	// time at start of simulation
double UpdateTime = 100*initial_dt;

// non-physical constants
int BoxX = 150;	// number of boxes per side of the grid
int BoxY = 150;
int BoxZ = 50;
int BoxZAgar = 30;
int maxLevels = 1;
double BoxLength = 2;
int FilterLen = 5;

// directory name
char DirName[500] = "";

// maximum number of cells in the simulation
int maxCells = 100000;

// nutrient constants
double Tortuosity = 2.0;
double KC = 0.001;
double C_rate = 3.7e-5*Tortuosity*5.0;
double maxGrowthRate = 1.0;
double Maintenance_rate = 0.0;
double Rc = 3.0;
double DiffColony = 1.0;
double DiffAgar = 3.0;
double maxCarbon = 0.01;

// colony constants


double Cdt = 0.1;
double ConvCrit = 1.0e-4;
int minIter = 500;
int maxIter = 20000;
int InterfaceCondition = 1; // 1: continuous $\partial C/partial n$; 2: flux continuity with qC; 3: continuous $\partial C/\partial t$;
bool NutrientGSI = 0;
int refinementGridHeight = 4;

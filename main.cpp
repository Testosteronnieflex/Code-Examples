#include <stdio.h>
#include <time.h>
#include <math.h>
#include <vector>

#include "mt19937.h"
#include "LinkedList.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NDIM 3

#pragma warning(disable : 4996)

const double cuberoot_2 = pow(2.0, 1.0 / 3.0);
const double hexaroot_2 = pow(2.0, 1.0 / 6.0);


/*Calculates a random double between min and max*/
double randd(double min, double max)
{
	return (dsfmt_genrand() * (max - min)) + min;
}

/*Generate random number according to gaussian distribution*/
double random_BoxMuller(double sigma, double mu)
{
	double u, v;
	double s = 2;
	double z0;

	do {
		u = dsfmt_genrand() * 2.0 - 1.0;
		v = dsfmt_genrand() * 2.0 - 1.0;
		s = u * u + v * v;
	} while (s >= 1.0 || s == 0);

	z0 = u * sqrt(-2 * log(s) / s);
	return z0 * sigma + mu;
}

/*Calculates a mod n*/
int mod(int a, int n)
{
	int r = a % n;
	if (r < 0)
		r += n;

	return r;
}

struct Particle {
	double coords[NDIM];		//Coordinates of particle
	//double coords_old[NDIM];//Coordinates of the previous timestep
	double dr[NDIM];			//r(t) - r(t - dt)
	double v[NDIM];				//Velocity of particle
	double v_init[NDIM];		//Velocity at t=0
	double f[NDIM];				//Force on particle
	double f_old[NDIM];			//Force on particle in previous step
	double sigma;				//Diameter of particle
	double m;					//Mass
	double gamma;				//Friction coefficient (3 pi sigma mu)
	double kappa;				//Spring constant that models interaction strenght with interface
	int cell_coords[NDIM];		//Coordinates of cell list where particle is stored
	std::vector<int> verlet;	//List of neighbours (so verlet list)
	double dr_verlet[NDIM];		//Displacement since last verlet list update
	double rcut_max;			//Maximum cutoff distance for particle-particle interactions
};

class Config {
public:
	Particle* config;		//Array of particle objects
	double box[NDIM];
	Cells cell_lists;
	double r_cell;			//Minimum size of the cells
	double rv;				//Cut off radius on top of interaction cutoff for verlet lists (so thickness of skin)
	int n;					//Current number of timesteps in between verlet list updates
	int n_last;				//Number of timesteps in between the last two verlet list updates

	int Ntot;			//Number of particles in conf
	int Ns;				//N of small particles
	int Nl;				//N of large particles
	double Nr;			//Ratio Ns/Nl
	double sigmar;		//Ratio sigma_l / sigma_s
	double T;			//Temperature
	double Epot;		//Potential energy
	double Ekin;		//Kinetic energy
	double nu;			//Coupling strength to heat bath
	double v0_vt;		//Average velocity autocorrelation
	double H;			//Height of ceiling
	double mu;			//Viscosity of the solvent
	double std_s;  //The standard deviaton in the radius of the small particles
	double std_l;  //The standard deviation in the radius of the big particles
	double max_r;  //The maximum radius of the particles for the cell lists


	/*------------------Methods-----------------*/
public:
	/*Generate a cubic lattice of Nx x Ny x Nz particles with volume fraction eta_0*/
	void genCubic(int Nx, int Ny, int Nz, double eta_0, double gamma_s);

	/*Generate a cubic lattice, but now with a distrubtion of radii*/
	void genCubicGaussian(int Nx, int Ny, int Nz, double eta_0, double gamma_s, double std_s, double std_l);

	/*Calculates the maximum radius of the particles*/
	void getMaxRadius();

	/*Initialize the cell lists*/
	void initCells();

	/*Check if particle needs to be moved to different cell and do move if necessary*/
	void moveCell(int ID);

	/*Initialize velocities at system temperature T*/
	void initVelocity(double Tinit);

	/*Store configuration in file*/
	void write_data(const char filename[]);

	/*Store configuration in single cell with coords cx, cy, cz in file*/
	void write_cell(const char filename[], int cx, int cy, int cz);

	/*Store the particles in the verlet list of particle with index ID*/
	void write_verlet(const char filename[], int ID);

	/*Calculate force between particle and interface*/
	void Fint(int partID);

	/*Calculate force between particles due to WCA potential*/
	void FWCA();

	void FWCA_cell();

	void FWCA_verlet();

	/*Calculate the forces on the particles. Also calculates Epot*/
	void forces(double dt);

	/*Calculates particle interactions and drag force*/
	void forcesLangevin();

	/*Integrate one timestep dt*/
	void integrateVelocityVerlet(double dt);

	/*Integrate one timestep dt with brownian dynamics*/
	void integrateBrownianDynamics(double dt);

	/*Integrate one timestep dt with langevin dynamics*/
	void integrateLangevinDynamics(double dt);

	/*Calculate the mean sqared displacement at current time*/
	double MSD();

	/*Calculate all the distances between particles and generate a verlet list, rv is the extended radius on top of cutoff radius*/
	void VerletList();

	/*Constructs verlet list using cell lists*/
	void VerletList_cell();

	/*calculates the average velocity*/
	double averageVelocity();

	/*Stores a density profile in z direction of both species of particle*/
	void densityProfile();
};

void Config::getMaxRadius()
{
	max_r=config[0].sigma;
	for (int i=1;i<Ntot;i++)
	{
		if (max_r<config[i].sigma)
		{
			max_r=config[i].sigma;
		}
	}
}

void Config::VerletList()
{	//r_v is the extra distance on top of the radii of the particles that are kept in the list
	n_last = n;
	n = 0;

	//If the verlet is calculated again, empty the verlet list
	for (int i = 0; i < Ntot; i++)
	{
		config[i].verlet.clear();
		config[i].dr_verlet[0] = 0.0;
		config[i].dr_verlet[1] = 0.0;
		config[i].dr_verlet[2] = 0.0;
	}
	//loop over all the particles to check there distances
	double distance[NDIM];	//Contains the distances in the x,y and z direction
	double r2;				//Squared the distance
	double r;				//The distance
	for (int i = 0; i < Ntot - 1; i++){
		for (int j = i + 1; j < Ntot; j++){
			//calculate the distance
			distance[0] = config[i].coords[0] - config[j].coords[0];
			distance[1] = config[i].coords[1] - config[j].coords[1];
			distance[2] = config[i].coords[2] - config[j].coords[2];
			distance[0] -= (int)(2.0 * distance[0] / box[0]) * box[0];		//Nearest image
			distance[1] -= (int)(2.0 * distance[1] / box[1]) * box[1];		//Nearest image

			r2 = distance[0] * distance[0] + distance[1] * distance[1] + distance[2] * distance[2];
			r = sqrt(r2);
			//Correct the distance for the radii of the particles
			//r -= hexaroot_2 * (config[i].sigma + config[j].sigma);
			//If the distance is less than r_v append the the neighbouring atom to the verlet list
			if (r < (rv + (config[i].sigma + config[j].sigma) / 2.0)){
				config[i].verlet.push_back(j);
				config[j].verlet.push_back(i);
			}
		}
	}
}

void Config::VerletList_cell()
{	//r_v is the extra distance on top of the radii of the particles that are kept in the list
	n_last = n;
	n = 0;

	//If the verlet is calculated again, empty the verlet list
	for (int i = 0; i < Ntot; i++)
	{
		config[i].verlet.clear();
		config[i].dr_verlet[0] = 0.0;
		config[i].dr_verlet[1] = 0.0;
		config[i].dr_verlet[2] = 0.0;
	}

	double dist[NDIM];
	double r2;
	double r;

	/*Go over all cell lists*/
	int cx, cy, cz;
	int ix, iy;
	int czmax, czmin;
	int ID;
	LinkedList* cell_ptr = NULL;
	for (int ip = 0; ip < Ntot; ip++) {
		/*Particle - Particle interaction*/

		/*Cell in which particle ip is located*/
		cx = config[ip].cell_coords[0];
		cy = config[ip].cell_coords[1];
		cz = config[ip].cell_coords[2];

		if (cz == 0) {
			czmin = 0;
			czmax = 1;
		}
		else if (cz == cell_lists.nz - 1) {
			czmin = cz - 1;
			czmax = cz;
		}
		else {
			czmin = cz - 1;
			czmax = cz + 1;
		}

		for (int i = cx - 1; i <= cx + 1; i++) {
			for (int j = cy - 1; j <= cy + 1; j++) {
				for (int iz = czmin; iz <= czmax; iz++) {
					ix = mod(i, cell_lists.nx);
					iy = mod(j, cell_lists.ny);
					cell_ptr = cell_lists.access(ix, iy, iz);

					if (ix == cx && iy == cy && iz == cz) {//Particle i is in this cell
						for (int j = 0; j < cell_ptr->len(); j++) {
							ID = cell_ptr->iterate();
							if (ID == ip)//Don't calculate force of particle with itself
								continue;

							dist[0] = config[ip].coords[0] - config[ID].coords[0];
							dist[1] = config[ip].coords[1] - config[ID].coords[1];
							dist[2] = config[ip].coords[2] - config[ID].coords[2];
							dist[0] -= (int)(2.0 * dist[0] / box[0]) * box[0];		//Nearest image
							dist[1] -= (int)(2.0 * dist[1] / box[1]) * box[1];		//Nearest image

							r2 = dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2];

							r = sqrt(r2);
							//Correct the distance for the radii of the particles
							//r -= hexaroot_2 * (config[ip].sigma + config[ID].sigma);
							//If the distance is less than r_v append the the neighbouring atom to the verlet list
							if (r < (rv + (config[ip].sigma + config[ID].sigma) / 2.0)) {
								config[ip].verlet.push_back(ID);
							}
						}
					}
					else {
						for (int j = 0; j < cell_ptr->len(); j++) {
							ID = cell_ptr->iterate();

							dist[0] = config[ip].coords[0] - config[ID].coords[0];
							dist[1] = config[ip].coords[1] - config[ID].coords[1];
							dist[2] = config[ip].coords[2] - config[ID].coords[2];
							dist[0] -= (int)(2.0 * dist[0] / box[0]) * box[0];		//Nearest image
							dist[1] -= (int)(2.0 * dist[1] / box[1]) * box[1];		//Nearest image

							r2 = dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2];

							r = sqrt(r2);
							//Correct the distance for the radii of the particles
							//r -= hexaroot_2 * (config[ip].sigma + config[ID].sigma);
							//If the distance is less than r_v append the the neighbouring atom to the verlet list
							if (r < (rv + (config[ip].sigma + config[ID].sigma) / 2.0)) {
								config[ip].verlet.push_back(ID);
							}
						}
					}
					cell_ptr->reset_iterate();
				}
			}
		}
	}
}

double Config::averageVelocity()
{
	double vx, vy, vz, vnet;
	//The sum of all velocities
	double vtot = 0.0;
	//Loop over all particles
	for (int i = 0; i < Ntot; i++)
	{
		vx = config[i].v[0];
		vy = config[i].v[1];
		vz = config[i].v[2];
		vnet = sqrt(vx * vx + vy * vy + vz * vz);
		vtot += vnet;
	}
	//Return the average velocity
	return vtot / Ntot;
}

void Config::genCubic(int Nx, int Ny, int Nz, double eta_0, double gamma_s)
{
	double x, y, z;
	double PNs;		//Probability of small particle
	double sigma_l;
	int index;
	int Ns_actual = 0;

	Ntot = Nx * Ny * Nz;
	Ns = (int)((Nr / (Nr + 1.0)) * Ntot);
	PNs = (double)Ns / (double)Ntot;
	sigma_l = sigmar;
	n = 0;

	/*Generate the particle array of the right size*/
	//config = (Particle*)malloc(sizeof(Particle) * Ntot);
	config = new Particle[Ntot];
	if (config == NULL) {
		Ntot = 0;
		box[0] = 0;
		box[1] = 0;
		box[2] = 0;
		return;
	}

	/*Place each particle at the right coordinates*/
	index = 0;
	z = 0.0;
	for (int iz = 0; iz < Nz; iz++) {
		y = 0.0;
		for (int iy = 0; iy < Ny; iy++) {
			x = 0.0;
			for (int ix = 0; ix < Nx; ix++) {
				if (index >= Ntot)
					return;
				config[index].coords[0] = x;
				config[index].coords[1] = y;
				config[index].coords[2] = z;
				config[index].dr[0] = 0.0;
				config[index].dr[1] = 0.0;
				config[index].dr[2] = 0.0;
				if (dsfmt_genrand() < PNs) {
					config[index].sigma = 1.0;
					config[index].m = 1.0;
					config[index].kappa = 50;
					config[index].gamma = gamma_s;
					config[index].rcut_max = hexaroot_2 * (1.0 + sigma_l) / 2.0;
					Ns_actual++;
				}
				else {
					config[index].sigma = sigma_l;
					config[index].m = pow(sigma_l, 3);
					config[index].kappa = 50 * sigmar * sigmar;
					config[index].gamma = sigmar * gamma_s;
					config[index].rcut_max = hexaroot_2 * sigma_l;
				}
				//config[index].gamma = 3 * M_PI * config[index].sigma * viscosity;
				config[index].verlet.reserve(50);
				index++;
				x += 1.0;
			}
			y += 1.0;
		}
		z += 1.0;
	}

	Ns = Ns_actual;
	Nl = Ntot - Ns;
	Nr = (double)Ns / (double)Nl;

	/*Store the dimensions of the box*/
	box[0] = (double)Nx;
	box[1] = (double)Ny;
	box[2] = (double)Nz;

	/*Set the right volume fraction*/
	double volume = box[0] * box[1] * box[2];
	double Vpart = (1.0 / 6.0) * M_PI * ((double)Ns + Nl * pow(sigma_l, 3.0));
	double Vtot = Vpart / eta_0;
	double scale_factor = pow(Vtot / volume, (1.0 / 3.0));

	box[0] *= scale_factor;
	box[1] *= scale_factor;
	box[2] *= scale_factor;
	volume = box[0] * box[1] * box[2];
	printf("eta_0 = %lf\n", Vpart / volume);

	for (int i = 0; i < Ntot; i++) {
		config[i].coords[0] *= scale_factor;
		config[i].coords[1] *= scale_factor;
		config[i].coords[2] = config[i].coords[2] * scale_factor + sigmar / 2.0;
	}

	H = box[2];	//Set height of interface
}

void Config::genCubicGaussian(int Nx, int Ny, int Nz, double eta_0, double gamma_s, double std_s, double std_l)
{
	double x, y, z;
	double PNs;		//Probability of small particle
	double sigma_l;
	double radius; //The radius of each particle
	int index;
	int Ns_actual = 0;

	Ntot = Nx * Ny * Nz;
	Ns = (int)((Nr / (Nr + 1.0)) * Ntot);
	PNs = (double)Ns / (double)Ntot;
	sigma_l = sigmar;
	n = 0;

	/*Generate the particle array of the right size*/
	//config = (Particle*)malloc(sizeof(Particle) * Ntot);
	config = new Particle[Ntot];
	if (config == NULL) {
		Ntot = 0;
		box[0] = 0;
		box[1] = 0;
		box[2] = 0;
		return;
	}

	/*Place each particle at the right coordinates*/
	index = 0;
	z = 0.0;
	for (int iz = 0; iz < Nz; iz++) {
		y = 0.0;
		for (int iy = 0; iy < Ny; iy++) {
			x = 0.0;
			for (int ix = 0; ix < Nx; ix++) {
				if (index >= Ntot)
					return;
				config[index].coords[0] = x;
				config[index].coords[1] = y;
				config[index].coords[2] = z;
				config[index].dr[0] = 0.0;
				config[index].dr[1] = 0.0;
				config[index].dr[2] = 0.0;
				if (dsfmt_genrand() < PNs) {
					radius=random_BoxMuller(std_s,1.0);
					//config[index].rcut_max = hexaroot_2 * (radius + sigma_l) / 2.0;
					Ns_actual++;
				}
				else {
					radius=random_BoxMuller(std_l,sigma_l);
					//config[index].rcut_max = hexaroot_2 * sigma_l;
				}
				config[index].sigma = radius;
				config[index].m = pow(radius, 3);
				config[index].kappa = 50*radius*radius;
				config[index].gamma = radius * gamma_s;
				//config[index].gamma = 3 * M_PI * config[index].sigma * viscosity;
				config[index].verlet.reserve(50);
				index++;
				x += 1.0;
			}
			y += 1.0;
		}
		z += 1.0;
	}
	//Imidiatly do the calculation for the max radius here so it doesnt have to be done anymore
	getMaxRadius();
	Ns = Ns_actual;
	Nl = Ntot - Ns;
	Nr = (double)Ns / (double)Nl;

	/*Store the dimensions of the box*/
	box[0] = (double)Nx;
	box[1] = (double)Ny;
	box[2] = (double)Nz;

	/*Set the right volume fraction*/
	double volume = box[0] * box[1] * box[2];
	double Vpart = (1.0 / 6.0) * M_PI * ((double)Ns + Nl * pow(sigma_l, 3.0));
	double Vtot = Vpart / eta_0;
	double scale_factor = pow(Vtot / volume, (1.0 / 3.0));

	box[0] *= scale_factor;
	box[1] *= scale_factor;
	box[2] *= scale_factor;
	volume = box[0] * box[1] * box[2];
	printf("eta_0 = %lf\n", Vpart / volume);

	for (int i = 0; i < Ntot; i++) {
		config[i].coords[0] *= scale_factor;
		config[i].coords[1] *= scale_factor;
		config[i].coords[2] = config[i].coords[2] * scale_factor + sigmar / 2.0;
	}

	H = box[2];	//Set height of interface
}

void Config::initCells()
{
	int nx, ny, nz;
	r_cell = max_r * hexaroot_2;

	nx = (int)floor(box[0] / r_cell);
	ny = (int)floor(box[1] / r_cell);
	nz = (int)floor((box[2] + r_cell) / r_cell);
	cell_lists.create(nx, ny, nz);
	cell_lists.sx = box[0] / (double)nx;
	cell_lists.sy = box[1] / (double)ny;
	cell_lists.sz = box[2] / (double)nz;

	printf("Generated %d x %d x %d cells\n", nx, ny, nz);

	/*Place all particles*/
	int ix, iy, iz;
	for (int i = 0; i < Ntot; i++) {
		ix = (int)floor(config[i].coords[0] / cell_lists.sx);
		iy = (int)floor(config[i].coords[1] / cell_lists.sx);
		iz = (int)floor(config[i].coords[2] / cell_lists.sy);
		cell_lists.access(ix, iy, iz)->push(i);
		config[i].cell_coords[0] = ix;
		config[i].cell_coords[1] = iy;
		config[i].cell_coords[2] = iz;
	}
}

void Config::write_data(const char filename[])
{
	double lbox;
	FILE* f = fopen(filename, "w");

	fprintf(f, "%d\n", Ntot);

	for (int i = 0; i < NDIM; i++) {
		lbox = box[i] / 2.0;
		fprintf(f, "%lf %lf\n", -lbox, lbox);
	}

	for (int i = 0; i < Ntot; i++) {
		for (int d = 0; d < NDIM; d++)
			fprintf(f, "%lf\t", config[i].coords[d] - box[d] / 2.0);
		if (config[i].sigma > 1.0)
			fprintf(f, "%lf\t%d\n", config[i].sigma, 1);
		else
			fprintf(f, "%lf\t%d\n", config[i].sigma, 2);
	}

	fclose(f);
}

void Config::write_cell(const char filename[], int cx, int cy, int cz)
{
	double lbox;
	LinkedList* cell_ptr = cell_lists.access(cx, cy, cz);
	if (cell_ptr == NULL) {
		return;
	}

	FILE* f = fopen(filename, "w");

	fprintf(f, "%d\n", cell_ptr->len());

	for (int i = 0; i < NDIM; i++) {
		lbox = box[i] / 2.0;
		fprintf(f, "%lf %lf\n", -lbox, lbox);
	}

	int index;
	for (int i = 0; i < cell_ptr->len(); i++) {
		index = cell_ptr->iterate();
		for (int d = 0; d < NDIM; d++)
			fprintf(f, "%lf\t", config[index].coords[d] - box[d] / 2.0);
		if (config[i].sigma > 1.0)
			fprintf(f, "%lf\t%d\n", config[index].sigma, 1);
		else
			fprintf(f, "%lf\t%d\n", config[index].sigma, 2);
	}

	cell_ptr->reset_iterate();

	fclose(f);
}

void Config::write_verlet(const char filename[], int ID)
{
	double lbox;
	FILE* f = fopen(filename, "w");

	fprintf(f, "%d\n", (int)config[ID].verlet.size() + 1);

	for (int i = 0; i < NDIM; i++) {
		lbox = box[i] / 2.0;
		fprintf(f, "%lf %lf\n", -lbox, lbox);
	}

	for (int d = 0; d < NDIM; d++) {
		fprintf(f, "%lf\t", config[ID].coords[d] - box[d] / 2.0);
	}
	if (config[ID].sigma > 1.0)
		fprintf(f, "%lf\t%d\n", config[ID].sigma, 1);
	else
		fprintf(f, "%lf\t%d\n", config[ID].sigma, 2);

	for (int i = 0; i < config[ID].verlet.size(); i++) {
		for (int d = 0; d < NDIM; d++)
			fprintf(f, "%lf\t", config[config[ID].verlet[i]].coords[d] - box[d] / 2.0);
		if (config[config[ID].verlet[i]].sigma > 1.0)
			fprintf(f, "%lf\t%d\n", config[config[ID].verlet[i]].sigma, 1);
		else
			fprintf(f, "%lf\t%d\n", config[config[ID].verlet[i]].sigma, 2);
	}

	fclose(f);
}

//void Config::initVelocity(double Tinit)
//{
//	double v2sum = 0.0;
//	double vsum[NDIM] = { 0.0, 0.0, 0.0 };
//	double vi;
//	double fs;
//
//	T = Tinit;
//
//	for (int i = 0; i < Ntot; i++) {
//		for (int d = 0; d < NDIM; d++) {
//			vi = randd(-1.0, 1.0);
//			config[i].v[d] = vi;
//			v2sum += config[i].m * vi * vi;
//			vsum[d] += vi;
//		}
//	}
//
//	v2sum /= (double)Ntot;
//	for (int d = 0; d < NDIM; d++) {
//		vsum[d] /= (double)Ntot;
//	}
//	fs = sqrt(3.0 * Tinit / v2sum);
//	for (int i = 0; i < Ntot; i++) {
//		for (int d = 0; d < NDIM; d++) {
//			config[i].v[d] = (config[i].v[d] - vsum[d]) * fs;
//			config[i].v_init[d] = config[i].v[d];
//			config[i].f[d] = 0.0;
//		}
//	}
//}

void Config::initVelocity(double Tinit)
{
	T = Tinit;

	for (int i = 0; i < Ntot; i++) {
		for (int d = 0; d < NDIM; d++) {
			config[i].v[d] = 0.0;
			config[i].v_init[d] = config[i].v[d];
			config[i].f[d] = 0.0;
		}
	}
}

void Config::FWCA()
{
	double dist[NDIM];
	double r2;
	double r2i;
	double r6i;
	double fradial;
	double sigma12, s2, s6;
	double r2_cut;
	double r2_floor;

	Epot = 0.0;

	for (int i = 0; i < Ntot; i++) {
		for (int d = 0; d < NDIM; d++) {
			config[i].f_old[d] = config[i].f[d];
			config[i].f[d] = 0.0;
		}
	}

	for (int i = 0; i < Ntot; i++) {
		/*Particle - Particle forces*/
		for (int j = i + 1; j < Ntot; j++) {
			dist[0] = config[i].coords[0] - config[j].coords[0];
			dist[1] = config[i].coords[1] - config[j].coords[1];
			dist[2] = config[i].coords[2] - config[j].coords[2];
			dist[0] -= (int)(2.0 * dist[0] / box[0]) * box[0];		//Nearest image
			dist[1] -= (int)(2.0 * dist[1] / box[1]) * box[1];		//Nearest image

			r2 = dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2];

			sigma12 = (config[i].sigma + config[j].sigma) / 2.0;
			s2 = sigma12 * sigma12;
			r2_cut = cuberoot_2 * s2;

			if (r2 < r2_cut) {
				r2i = 1.0 / r2;
				r6i = pow(r2i, 3.0);
				s6 = pow(s2, 3.0);
				fradial = 48.0 * r2i * r6i * s6 * (s6 * r6i - 0.5);
				Epot += 4.0 * r6i * s6 * (s6 * r6i - 1.0) + 1.0;
				for (int d = 0; d < NDIM; d++) {
					config[i].f[d] += fradial * dist[d];
					config[j].f[d] -= fradial * dist[d];
				}
			}
		}

		/*Particle - floor forces and test particle - ceiling forces*/
		r2_floor = config[i].coords[2] * config[i].coords[2];
		s2 = config[i].sigma * config[i].sigma / 4.0;
		r2_cut = cuberoot_2 * s2;


		if (r2_floor < r2_cut) {
			r2i = 1.0 / r2_floor;
			r6i = pow(r2i, 3.0);
			s6 = pow(s2, 3.0);
			fradial = 48.0 * r2i * r6i * s6 * (s6 * r6i - 0.5);
			Epot += 4.0 * r6i * s6 * (s6 * r6i - 1.0) + 1.0;
			config[i].f[2] += fradial * config[i].coords[2];
		}

		Fint(i); //Interaction with interface

	}
}

void Config::FWCA_cell()
{
	double dist[NDIM];
	double r2;
	double r2i;
	double r6i;
	double fradial;
	double sigma12, s2, s6;
	double r2_cut;
	double r2_floor;

	Epot = 0.0;

	for (int i = 0; i < Ntot; i++) {
		for (int d = 0; d < NDIM; d++) {
			config[i].f_old[d] = config[i].f[d];
			config[i].f[d] = 0.0;
		}
	}

	/*Go over all cell lists*/
	int cx, cy, cz;
	int ix, iy;
	int czmax, czmin;
	int ID;
	LinkedList* cell_ptr = NULL;
	for (int ip = 0; ip < Ntot; ip++) {
		/*Particle - Particle interaction*/

		/*Cell in which particle ip is located*/
		cx = config[ip].cell_coords[0];
		cy = config[ip].cell_coords[1];
		cz = config[ip].cell_coords[2];

		if (cz == 0) {
			czmin = 0;
			czmax = 1;
		}
		else if (cz == cell_lists.nz - 1) {
			czmin = cz - 1;
			czmax = cz;
		}
		else {
			czmin = cz - 1;
			czmax = cz + 1;
		}

		for (int i = cx - 1; i <= cx + 1; i++) {
			for (int j = cy - 1; j <= cy + 1; j++) {
				for (int iz = czmin; iz <= czmax; iz++) {
					ix = mod(i, cell_lists.nx);
					iy = mod(j, cell_lists.ny);
					cell_ptr = cell_lists.access(ix, iy, iz);

					if (ix == cx && iy == cy && iz == cz) {//Particle i is in this cell
						for (int j = 0; j < cell_ptr->len(); j++) {
							ID = cell_ptr->iterate();
							if (ID == ip)//Don't calculate force of particle with itself
								continue;

							r2 = 0.0;
							dist[0] = config[ip].coords[0] - config[ID].coords[0];
							dist[1] = config[ip].coords[1] - config[ID].coords[1];
							dist[2] = config[ip].coords[2] - config[ID].coords[2];
							dist[0] -= (int)(2.0 * dist[0] / box[0]) * box[0];		//Nearest image
							dist[1] -= (int)(2.0 * dist[1] / box[1]) * box[1];		//Nearest image

							r2 += dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2];

							sigma12 = (config[ip].sigma + config[ID].sigma) / 2.0;
							s2 = sigma12 * sigma12;
							r2_cut = cuberoot_2 * s2;

							if (r2 < r2_cut) {
								r2i = 1.0 / r2;
								r6i = pow(r2i, 3.0);
								s6 = pow(s2, 3.0);
								fradial = 48.0 * r2i * r6i * s6 * (s6 * r6i - 0.5);
								Epot += 4.0 * r6i * s6 * (s6 * r6i - 1.0) + 1.0;
								for (int d = 0; d < NDIM; d++) {
									config[ip].f[d] += fradial * dist[d];
								}
							}
						}
					}
					else {
						for (int j = 0; j < cell_ptr->len(); j++) {
							ID = cell_ptr->iterate();

							r2 = 0.0;
							dist[0] = config[ip].coords[0] - config[ID].coords[0];
							dist[1] = config[ip].coords[1] - config[ID].coords[1];
							dist[2] = config[ip].coords[2] - config[ID].coords[2];
							dist[0] -= (int)(2.0 * dist[0] / box[0]) * box[0];		//Nearest image
							dist[1] -= (int)(2.0 * dist[1] / box[1]) * box[1];		//Nearest image

							r2 += dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2];

							sigma12 = (config[ip].sigma + config[ID].sigma) / 2.0;
							s2 = sigma12 * sigma12;
							r2_cut = cuberoot_2 * s2;

							if (r2 < r2_cut) {
								r2i = 1.0 / r2;
								r6i = pow(r2i, 3.0);
								s6 = pow(s2, 3.0);
								fradial = 48.0 * r2i * r6i * s6 * (s6 * r6i - 0.5);
								Epot += 4.0 * r6i * s6 * (s6 * r6i - 1.0) + 1.0;
								for (int d = 0; d < NDIM; d++) {
									config[ip].f[d] += fradial * dist[d];
								}
							}
						}
					}
					cell_ptr->reset_iterate();
				}
			}
		}

		/*Particle - floor forces*/
		r2_floor = config[ip].coords[2] * config[ip].coords[2];
		s2 = (config[ip].sigma + 1.0) * (config[ip].sigma + 1.0) / 4.0;
		r2_cut = cuberoot_2 * s2;


		if (r2_floor < r2_cut) {
			r2i = 1.0 / r2_floor;
			r6i = pow(r2i, 3.0);
			s6 = pow(s2, 3.0);
			fradial = 48.0 * r2i * r6i * s6 * (s6 * r6i - 0.5);
			Epot += 4.0 * r6i * s6 * (s6 * r6i - 1.0) + 1.0;
			config[ip].f[2] += fradial * config[ip].coords[2];
		}

		Fint(ip); //Interaction with interface
	}
}

void Config::FWCA_verlet()
{
	double dist[NDIM];
	double r2;
	double r2i;
	double r6i;
	double fradial;
	double sigma12, s2, s6;
	double r2_cut;
	double r2_floor;

	Epot = 0.0;
	//reset the forces
	for (int i = 0; i < Ntot; i++) {
		for (int d = 0; d < NDIM; d++) {
			config[i].f_old[d] = config[i].f[d];
			config[i].f[d] = 0.0;
		}
	}

	//Loop over all particles and calculate all the forces from the neighbours in the Verlet list
	int neighbours;
	int j;
	for (int i = 0; i < Ntot; i++) {
		neighbours = config[i].verlet.size();
		/*Particle - Particle forces*/
		for (int k = 0; k < neighbours; k++) {
			j = config[i].verlet[k];
			r2 = 0.0;
			dist[0] = config[i].coords[0] - config[j].coords[0];
			dist[1] = config[i].coords[1] - config[j].coords[1];
			dist[2] = config[i].coords[2] - config[j].coords[2];
			dist[0] -= (int)(2.0 * dist[0] / box[0]) * box[0];		//Nearest image
			dist[1] -= (int)(2.0 * dist[1] / box[1]) * box[1];		//Nearest image

			r2 += dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2];

			sigma12 = (config[i].sigma + config[j].sigma) / 2.0;
			s2 = sigma12 * sigma12;
			r2_cut = cuberoot_2 * s2;

			if (r2 < r2_cut) {
				r2i = 1.0 / r2;
				r6i = pow(r2i, 3.0);
				s6 = pow(s2, 3.0);
				fradial = 48.0 * r2i * r6i * s6 * (s6 * r6i - 0.5);
				Epot += 4.0 * r6i * s6 * (s6 * r6i - 1.0) + 1.0;
				for (int d = 0; d < NDIM; d++) {
					config[i].f[d] += fradial * dist[d];
				}
			}
		}

		/*Particle - floor forces and test particle - ceiling forces*/
		r2_floor = config[i].coords[2] * config[i].coords[2];
		s2 = config[i].sigma * config[i].sigma / 4.0;
		r2_cut = cuberoot_2 * s2;

		//config[i].f[2] -= config[i].m * 1.0; //gravity for test

		if (r2_floor < r2_cut) {
			r2i = 1.0 / r2_floor;
			r6i = pow(r2i, 3.0);
			s6 = pow(s2, 3.0);
			fradial = 48.0 * r2i * r6i * s6 * (s6 * r6i - 0.5);
			Epot += 4.0 * r6i * s6 * (s6 * r6i - 1.0) + 1.0;
			config[i].f[2] += fradial * config[i].coords[2];
		}

		Fint(i); //Interaction with interface
	}
}

void Config::forces(double dt)
{
	double dist[NDIM];
	double r2;
	double r2i;
	double r6i;
	double fradial;
	double sigma12, s2, s6;
	double r2_cut;
	double r2_floor;

	Epot = 0.0;

	for (int i = 0; i < Ntot; i++) {
		for (int d = 0; d < NDIM; d++) {
			config[i].f_old[d] = config[i].f[d];
			config[i].f[d] = 0.0;
		}
	}

	/*Go over all cell lists*/
	int cx, cy, cz;
	int ix, iy;
	int czmax;
	int czmin;
	int ID;
	LinkedList* cell_ptr = NULL;
	for (int ip = 0; ip < Ntot; ip++) {
		/*Particle - Particle interaction*/

		/*Cell in which particle ip is located*/
		cx = config[ip].cell_coords[0];
		cy = config[ip].cell_coords[1];
		cz = config[ip].cell_coords[2];

		if (cz == 0) {
			czmin = 0;
			czmax = 1;
		}
		else if (cz == cell_lists.nz - 1) {
			czmin = cz - 1;
			czmax = cz;
		}
		else {
			czmin = cz - 1;
			czmax = cz + 1;
		}

		for (int i = cx - 1; i <= cx + 1; i++) {
			for (int j = cy - 1; j <= cy + 1; j++) {
				for (int iz = czmin; iz <= czmax; iz++) {
					ix = mod(i, cell_lists.nx);
					iy = mod(j, cell_lists.ny);
					cell_ptr = cell_lists.access(ix, iy, iz);

					if (ix == cx && iy == cy && iz == cz){//Particle i is in this cell
						for (int j = 0; j < cell_ptr->len(); j++) {
							ID = cell_ptr->iterate();
							if (ID == ip)//Don't calculate force of particle with itself
								continue;

							r2 = 0.0;
							dist[0] = config[ip].coords[0] - config[ID].coords[0];
							dist[1] = config[ip].coords[1] - config[ID].coords[1];
							dist[2] = config[ip].coords[2] - config[ID].coords[2];
							dist[0] -= (int)(2.0 * dist[0] / box[0]) * box[0];		//Nearest image
							dist[1] -= (int)(2.0 * dist[1] / box[1]) * box[1];		//Nearest image

							r2 += dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2];

							sigma12 = (config[ip].sigma + config[ID].sigma) / 2.0;
							s2 = sigma12 * sigma12;
							r2_cut = cuberoot_2 * s2;

							if (r2 < r2_cut) {
								r2i = 1.0 / r2;
								r6i = pow(r2i, 3.0);
								s6 = pow(s2, 3.0);
								fradial = 48.0 * r2i * r6i * s6 * (s6 * r6i - 0.5);
								Epot += 4.0 * r6i * s6 * (s6 * r6i - 1.0) + 1.0;
								for (int d = 0; d < NDIM; d++) {
									config[ip].f[d] += fradial * dist[d];
								}
							}
						}
					}
					else {
						for (int j = 0; j < cell_ptr->len(); j++) {
							ID = cell_ptr->iterate();

							r2 = 0.0;
							dist[0] = config[ip].coords[0] - config[ID].coords[0];
							dist[1] = config[ip].coords[1] - config[ID].coords[1];
							dist[2] = config[ip].coords[2] - config[ID].coords[2];
							dist[0] -= (int)(2.0 * dist[0] / box[0]) * box[0];		//Nearest image
							dist[1] -= (int)(2.0 * dist[1] / box[1]) * box[1];		//Nearest image

							r2 += dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2];

							sigma12 = (config[ip].sigma + config[ID].sigma) / 2.0;
							s2 = sigma12 * sigma12;
							r2_cut = cuberoot_2 * s2;

							if (r2 < r2_cut) {
								r2i = 1.0 / r2;
								r6i = pow(r2i, 3.0);
								s6 = pow(s2, 3.0);
								fradial = 48.0 * r2i * r6i * s6 * (s6 * r6i - 0.5);
								Epot += 4.0 * r6i * s6 * (s6 * r6i - 1.0) + 1.0;
								for (int d = 0; d < NDIM; d++) {
									config[ip].f[d] += fradial * dist[d];
								}
							}
						}
					}
					cell_ptr->reset_iterate();
				}
			}
		}

		/*Particle - floor forces*/
		r2_floor = config[ip].coords[2] * config[ip].coords[2];
		s2 = (config[ip].sigma + 1.0) * (config[ip].sigma + 1.0) / 4.0;
		r2_cut = cuberoot_2 * s2;


		if (r2_floor < r2_cut) {
			r2i = 1.0 / r2_floor;
			r6i = pow(r2i, 3.0);
			s6 = pow(s2, 3.0);
			fradial = 48.0 * r2i * r6i * s6 * (s6 * r6i - 0.5);
			Epot += 4.0 * r6i * s6 * (s6 * r6i - 1.0) + 1.0;
			config[ip].f[2] += fradial * config[ip].coords[2];
		}

		Fint(ip); //Interaction with interface

		/*Random brownian forces and friction*/
		for (int d = 0; d < NDIM; d++) {
			config[ip].f[d] += random_BoxMuller(sqrt(2 * config[ip].gamma * T), 0.0);
			config[ip].f[d] -= config[ip].gamma * config[ip].v[d];
		}
	}

	//for (int i = 0; i < Ntot; i++) {
	//	/*Particle - Particle forces*/
	//	for (int j = i + 1; j < Ntot; j++) {
	//		r2 = 0.0;
	//		dist[0] = config[i].coords[0] - config[j].coords[0];
	//		dist[1] = config[i].coords[1] - config[j].coords[1];
	//		dist[2] = config[i].coords[2] - config[j].coords[2];
	//		dist[0] -= (int)(2.0 * dist[0] / box[0]) * box[0];		//Nearest image
	//		dist[1] -= (int)(2.0 * dist[1] / box[1]) * box[1];		//Nearest image

	//		r2 += dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2];

	//		sigma12 = (config[i].sigma + config[j].sigma) / 2.0;
	//		s2 = sigma12 * sigma12;
	//		r2_cut = cuberoot_2 * s2;

	//		if (r2 < r2_cut) {
	//			r2i = 1.0 / r2;
	//			r6i = pow(r2i, 3.0);
	//			s6 = pow(s2, 3.0);
	//			fradial = 48.0 * r2i * r6i * s6 * (s6 * r6i - 0.5);
	//			Epot += 4.0 * r6i * s6 * (s6 * r6i - 1.0) + 1.0;
	//			for (int d = 0; d < NDIM; d++) {
	//				config[i].f[d] += fradial * dist[d];
	//				config[j].f[d] -= fradial * dist[d];
	//			}
	//		}
	//	}

	//	/*Particle - floor forces*/
	//	r2_floor = config[i].coords[2] * config[i].coords[2];
	//	s2 = (config[i].sigma + 1.0) * (config[i].sigma + 1.0) / 4.0;
	//	r2_cut = cuberoot_2 * s2;


	//	if (r2_floor < r2_cut) {
	//		r2i = 1.0 / r2_floor;
	//		r6i = pow(r2i, 3.0);
	//		s6 = pow(s2, 3.0);
	//		fradial = 48.0 * r2i * r6i * s6 * (s6 * r6i - 0.5);
	//		Epot += 4.0 * r6i * s6 * (s6 * r6i - 1.0) + 1.0;
	//		config[i].f[2] += fradial * config[i].coords[2];
	//	}

	//	Fint(i); //Interaction with interface

	//	/*Random brownian forces and friction*/
	//	for (int d = 0; d < NDIM; d++) {
	//		config[i].f[d] += random_BoxMuller(sqrt(2 * config[i].gamma * T / dt), 0.0);
	//		config[i].f[d] -= config[i].gamma * config[i].v[d];
	//	}
	//}
}

void Config::forcesLangevin()
{
	double dist[NDIM];
	double r2;
	double r2i;
	double r6i;
	double fradial;
	double sigma12, s2, s6;
	double r2_cut;
	double r2_floor;

	Epot = 0.0;

	for (int i = 0; i < Ntot; i++) {
		for (int d = 0; d < NDIM; d++) {
			config[i].f_old[d] = config[i].f[d];
			config[i].f[d] = 0.0;
		}
	}

	/*Go over all cell lists*/
	int cx, cy, cz;
	int ix, iy;
	int czmax;
	int czmin;
	int ID;
	LinkedList* cell_ptr = NULL;
	for (int ip = 0; ip < Ntot; ip++) {
		/*Particle - Particle interaction*/

		/*Cell in which particle ip is located*/
		cx = config[ip].cell_coords[0];
		cy = config[ip].cell_coords[1];
		cz = config[ip].cell_coords[2];

		if (cz == 0) {
			czmin = 0;
			czmax = 1;
		}
		else if (cz == cell_lists.nz) {
			czmin = cell_lists.nz - 1;
			czmax = cell_lists.nz;
		}
		else {
			czmin = cz - 1;
			czmax = cz + 1;
		}

		for (int i = cx - 1; i <= cx + 1; i++) {
			for (int j = cy - 1; j <= cy + 1; j++) {
				for (int iz = czmin; iz <= czmax; iz++) {
					ix = mod(i, cell_lists.nx);
					iy = mod(j, cell_lists.ny);
					cell_ptr = cell_lists.access(ix, iy, iz);

					if (ix == cx && iy == cy && iz == cz) {//Particle i is in this cell
						for (int j = 0; j < cell_ptr->len(); j++) {
							ID = cell_ptr->iterate();
							if (ID == ip)//Don't calculate force of particle with itself
								continue;

							r2 = 0.0;
							dist[0] = config[ip].coords[0] - config[ID].coords[0];
							dist[1] = config[ip].coords[1] - config[ID].coords[1];
							dist[2] = config[ip].coords[2] - config[ID].coords[2];
							dist[0] -= (int)(2.0 * dist[0] / box[0]) * box[0];		//Nearest image
							dist[1] -= (int)(2.0 * dist[1] / box[1]) * box[1];		//Nearest image

							r2 += dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2];

							sigma12 = (config[ip].sigma + config[ID].sigma) / 2.0;
							s2 = sigma12 * sigma12;
							r2_cut = cuberoot_2 * s2;

							if (r2 < r2_cut) {
								r2i = 1.0 / r2;
								r6i = pow(r2i, 3.0);
								s6 = pow(s2, 3.0);
								fradial = 48.0 * r2i * r6i * s6 * (s6 * r6i - 0.5);
								Epot += 4.0 * r6i * s6 * (s6 * r6i - 1.0) + 1.0;
								for (int d = 0; d < NDIM; d++) {
									config[ip].f[d] += fradial * dist[d];
								}
							}
						}
					}
					else {
						for (int j = 0; j < cell_ptr->len(); j++) {
							ID = cell_ptr->iterate();

							r2 = 0.0;
							dist[0] = config[ip].coords[0] - config[ID].coords[0];
							dist[1] = config[ip].coords[1] - config[ID].coords[1];
							dist[2] = config[ip].coords[2] - config[ID].coords[2];
							dist[0] -= (int)(2.0 * dist[0] / box[0]) * box[0];		//Nearest image
							dist[1] -= (int)(2.0 * dist[1] / box[1]) * box[1];		//Nearest image

							r2 += dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2];

							sigma12 = (config[ip].sigma + config[ID].sigma) / 2.0;
							s2 = sigma12 * sigma12;
							r2_cut = cuberoot_2 * s2;

							if (r2 < r2_cut) {
								r2i = 1.0 / r2;
								r6i = pow(r2i, 3.0);
								s6 = pow(s2, 3.0);
								fradial = 48.0 * r2i * r6i * s6 * (s6 * r6i - 0.5);
								Epot += 4.0 * r6i * s6 * (s6 * r6i - 1.0) + 1.0;
								for (int d = 0; d < NDIM; d++) {
									config[ip].f[d] += fradial * dist[d];
								}
							}
						}
					}
					cell_ptr->reset_iterate();
				}
			}
		}

		/*Particle - floor forces*/
		r2_floor = config[ip].coords[2] * config[ip].coords[2];
		s2 = (config[ip].sigma + 1.0) * (config[ip].sigma + 1.0) / 4.0;
		r2_cut = cuberoot_2 * s2;


		if (r2_floor < r2_cut) {
			r2i = 1.0 / r2_floor;
			r6i = pow(r2i, 3.0);
			s6 = pow(s2, 3.0);
			fradial = 48.0 * r2i * r6i * s6 * (s6 * r6i - 0.5);
			Epot += 4.0 * r6i * s6 * (s6 * r6i - 1.0) + 1.0;
			config[ip].f[2] += fradial * config[ip].coords[2];
		}

		Fint(ip); //Interaction with interface

		/*Random brownian forces and friction*/
		for (int d = 0; d < NDIM; d++) {
			config[ip].f[d] -= config[ip].gamma * config[ip].v[d];
		}
	}
}

void Config::Fint(int partID)
{
	if (config[partID].coords[2] > H)
		config[partID].f[2] -= config[partID].kappa * (config[partID].coords[2] - H);
}

void Config::moveCell(int ID)
{
	/*Determine coords of required cell*/
	int ix, iy, iz;
	ix = (int)floor(config[ID].coords[0] / cell_lists.sx);
	iy = (int)floor(config[ID].coords[1] / cell_lists.sx);
	iz = (int)floor(config[ID].coords[2] / cell_lists.sy);

	/*Check if cell list needs updating*/
	if (ix == config[ID].cell_coords[0] && iy == config[ID].cell_coords[1] && iz == config[ID].cell_coords[2])
		return;

	//if (iz >= cell_lists.nz || iz < 0);
	//	write_data("ESCAPE.dat");

	cell_lists.access(config[ID].cell_coords[0], config[ID].cell_coords[1], config[ID].cell_coords[2])->delete_ID(ID);
	cell_lists.access(ix, iy, iz)->push(ID);
	config[ID].cell_coords[0] = ix;
	config[ID].cell_coords[1] = iy;
	config[ID].cell_coords[2] = iz;
}

double Config::MSD()
{
	double msd = 0.0;

	for (int i = 0; i < Ntot; i++) {
		for (int d = 0; d < NDIM; d++) {
			if (config[i].sigma == 1.0)
				msd += config[i].dr[d] * config[i].dr[d] / (double)Ns;
		}
	}

	return msd;
}

void Config::integrateVelocityVerlet(double dt)
{
	double rnew;
	int nbox;

	Ekin = 0.0;

	for (int i = 0; i < Ntot; i++) {
		for (int d = 0; d < NDIM - 1; d++) {
			rnew = config[i].coords[d] + config[i].v[d] * dt + (config[i].f[d] / (config[i].m * 2.0)) * dt * dt;
			config[i].dr[d] += rnew - config[i].coords[d];

			/*Periodic boundary conditions*/
			if (rnew >= box[d]) {
				nbox = (int)(rnew / box[d]);
				rnew -= nbox * box[d];
			}
			else if (rnew < 0) {
				nbox = (int)(rnew / box[d]) - 1;
				rnew -= nbox * box[d];
			}
			config[i].coords[d] = rnew;
		}

		rnew = config[i].coords[2] + config[i].v[2] * dt + (config[i].f[2] / (config[i].m * 2.0)) * dt * dt;
		config[i].dr[2] = rnew - config[i].coords[2];
		config[i].coords[2] = rnew;
		//config[i].coords[2] += config[i].v[2] * dt + (config[i].f[2] / (config[i].m * 2.0)) * dt * dt;
		moveCell(i);
	}

	forces(dt);

	for (int i = 0; i < Ntot; i++) {
		for (int d = 0; d < NDIM; d++) {
			config[i].v[d] += (config[i].f[d] + config[i].f_old[d]) * dt / (2.0 * config[i].m);
			Ekin += config[i].m * config[i].v[d] * config[i].v[d];
		}
	}

	//T = Ekin / (3.0 * (double)Ntot);
	Ekin *= 0.5;
}

void Config::integrateBrownianDynamics(double dt)
{
	/*Calculate particle-particle interaction forces and forces with wall/ceiling*/
	FWCA_verlet();

	Ekin = 0.0;

	/*Euler integration*/
	double rnew, dr, dr2_verlet;
	double rv_over2_squared = rv * rv / 4.0;
	int nbox;
	bool update_verlet = false;
	double random;
	for (int i = 0; i < Ntot; i++) {
		for (int d = 0; d < NDIM - 1; d++) {
			random = random_BoxMuller(1.0, 0.0);
			config[i].v[d] =  sqrt(2.0 * T / (config[i].gamma)) * random + config[i].f[d] / config[i].gamma;
			rnew = config[i].coords[d] + sqrt(2.0 * T * dt / config[i].gamma) * random + dt * config[i].f[d] / config[i].gamma;
			dr = rnew - config[i].coords[d];
			config[i].dr[d] += dr;
			config[i].dr_verlet[d] += dr;

			Ekin += config[i].m * config[i].v[d] * config[i].v[d];

			/*Periodic boundary conditions*/
			if (rnew >= box[d]) {
				nbox = (int)(rnew / box[d]);
				rnew -= nbox * box[d];
			}
			else if (rnew < 0) {
				nbox = (int)(rnew / box[d]) - 1;
				rnew -= nbox * box[d];
			}
			config[i].coords[d] = rnew;
		}

		random = random_BoxMuller(1.0, 0.0);
		config[i].v[2] = sqrt(2.0 * T / (config[i].gamma)) * random + config[i].f[2] / config[i].gamma;
		rnew = config[i].coords[2] + sqrt(2.0 * T * dt / config[i].gamma) * random + dt * config[i].f[2] / config[i].gamma;
		dr = rnew - config[i].coords[2];
		config[i].dr[2] += dr;
		config[i].dr_verlet[2] += dr;
		config[i].coords[2] = rnew;
		//config[i].coords[2] += sqrt(2.0 * T * dt / config[i].gamma) * random + dt * config[i].f[2] / config[i].gamma;
		/*config[i].v[2] = sqrt(2.0 * T / (config[i].gamma * dt)) * random_BoxMuller(1.0, 0.0) + config[i].f[2] / config[i].gamma;
		config[i].coords[2] += config[i].v[2] * dt;*/
		moveCell(i);

		Ekin += config[i].m * config[i].v[2] * config[i].v[2];
		dr2_verlet = config[i].dr_verlet[0] * config[i].dr_verlet[0] + config[i].dr_verlet[1] * config[i].dr_verlet[1] + config[i].dr_verlet[2] * config[i].dr_verlet[2];
		if (dr2_verlet > rv_over2_squared) {
			update_verlet = true;
		}
	}

	if (update_verlet) {
		//printf("Need to update verlet list!\n");
		VerletList_cell();
	}

	Ekin *= 0.5;
}

void Config::integrateLangevinDynamics(double dt)
{
	forcesLangevin();

	Ekin = 0.0;

	double rnew;
	int nbox;
	for (int i = 0; i < Ntot; i++) {
		for (int d = 0; d < NDIM - 1; d++) {
			rnew = config[i].coords[d] + config[i].v[d] * dt + config[i].f[d] * dt * dt / (config[i].m * 2);
			config[i].dr[d] += rnew - config[i].coords[d];
			Ekin += config[i].m * config[i].v[d] * config[i].v[d];


			/*Periodic boundary conditions*/
			if (rnew >= box[d]) {
				nbox = (int)(rnew / box[d]);
				rnew -= nbox * box[d];
			}
			else if (rnew < 0) {
				nbox = (int)(rnew / box[d]) - 1;
				rnew -= nbox * box[d];
			}
			config[i].coords[d] = rnew;

			config[i].v[d] += random_BoxMuller(sqrt(2.0 * dt * T *config[i].gamma) / config[i].m, 0.0) + config[i].f[d] * dt / config[i].m;
		}

		rnew = config[i].coords[2] + config[i].v[2] * dt + config[i].f[2] * dt * dt / (config[i].m * 2);
		config[i].dr[2] += rnew - config[i].coords[2];
		Ekin += config[i].m * config[i].v[2] * config[i].v[2];
		config[i].v[2] += random_BoxMuller(sqrt(2.0 * dt * T * config[i].gamma) / config[i].m, 0.0) + config[i].f[2] * dt / config[i].m;
		config[i].coords[2] = rnew;
		moveCell(i);
	}
	Ekin *= 0.5;
}

void Config::densityProfile()
{
	int profile_l[100];
	int profile_s[100];
	double dbin = (H + sigmar) / (double)100;

	for (int i = 0; i < 100; i++) {
		profile_l[i] = 0;
		profile_s[i] = 0;
	}

	int n;
	for (int i = 0; i < Ntot; i++) {
		n = (int)(config[i].coords[2] / dbin);
		if (n > 99) {
			printf("Density profile out of range!");
			return;
		}
		if (config[i].sigma == 1.0) {
			profile_s[n]++;
		}
		else
			profile_l[n]++;
	}

	char buffer[50];
	sprintf(buffer, "profile_H%.2lf.dat", H);
	FILE* fprofile = fopen(buffer, "w");
	for (int i = 0; i < 100; i++) {
		fprintf(fprofile, "%lf\t%lf\t%lf\n", i * dbin / H, (double)profile_s[i] / (double)Ns, (double)profile_l[i] / (double)Nl);
	}
	fclose(fprofile);
}

int main()
{
	dsfmt_seed(time(NULL));

	double dt = 0.0005;
	double gamma_s = 2.5;

	Config conf;
	conf.sigmar = 7.0;
	conf.Nr = 150.0;
	conf.genCubic(8, 8, 70, 0.02, gamma_s);
	char buffer[100];

	//clock_t start = clock();

	FILE* fmsd;

	conf.rv = conf.sigmar * 0.2;
	conf.initVelocity(1.0);
	conf.initCells();
	conf.write_data("init.dat");
	//conf.forces(dt);
	printf("Nr: %.2lf\tNl: %d\tNs: %d\tNtot: %d\n", conf.Nr, conf.Nl, conf.Ns, conf.Ntot);
	printf("Epot: %lf\n", conf.Epot);

	conf.VerletList();
	for (int i = 0; i < 10000; i++) { //Equilibrate
		conf.integrateBrownianDynamics(dt);
	}

	for (int i = 0; i < conf.Ntot; i++) {
		conf.config[i].dr[0] = 0.0;
		conf.config[i].dr[1] = 0.0;
		conf.config[i].dr[2] = 0.0;
		//conf.config[i].gamma = gamma_s * conf.config[i].sigma;
	}

	//conf.VerletList();
	//conf.write_verlet("verlet_init.dat", 1500);
	int i = 0;
	while (conf.H > 100){
		conf.n++;
		conf.integrateBrownianDynamics(dt);
		if (i % 50000 == 0) {
			sprintf(buffer, "snapshot_%d.dat", i);
			printf("snapshot: %d\tT: %lf\tH: %.2lf\tn: %d\n", i, conf.Ekin * 2.0 / (3.0 * conf.Ntot), conf.H, conf.n_last);
			conf.write_data(buffer);
			//fmsd = fopen("msd.dat", "a");
			//fprintf(fmsd, "%lf\t%lf\t%lf\n", i * dt, conf.MSD(), 6.0 * conf.T * i * dt / gamma_s);
			//fclose(fmsd);
			conf.densityProfile();
		}
		conf.H -= 0.02 * dt;
		i++;
	}

	conf.densityProfile();
	//conf.write_verlet("verlet_final.dat", 1500);


	//printf("%d\n", conf.cell_lists.access(1, 1, 1)->len());
	//conf.write_cell("cell_final.dat", 1, 1, 1);
	//conf.write_cell("celltest2.dat", 5, 5, 5);

	//clock_t end = clock();
	//printf("Time taken: %.2lf seconds\n", (double)(end - start) / (double)CLOCKS_PER_SEC);

	return 0;
}

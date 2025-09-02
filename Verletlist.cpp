#include <stdio.h>
#include <time.h>
#include <math.h>


#include "mt19937.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NDIM 3

#pragma warning(disable : 4996)

const double cuberoot_2 = pow(2.0, 1.0 / 3.0);
const double hexaroot_2=pow(2.0,1.0/6.0);

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

struct Particle {
	double coords[NDIM];	//Coordinates of particle
	//double coords_old[NDIM];//Coordinates of the previous timestep
	double dr[NDIM];		//r(t) - r(t - dt)
	double v[NDIM];			//Velocity of particle
	double v_init[NDIM];	//Velocity at t=0
	double f[NDIM];			//Force on particle
	double f_old[NDIM];		//Force on particle in previous step
	double sigma;			//Diameter of particle
	double m;				//Mass
	double gamma;			//Friction coefficient (3 pi sigma mu)
	double kappa;			//Spring constant that models interaction strenght with interface
	double verlet[200];    //List of all the neighbours TODO make the array shorterand adaptive in size
	int nver; 						//The number of particles in the Verlet list
};

class Config {
public:
	Particle* config;		//Array of particle objects
	double box[NDIM];

public:
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


	/*------------------Methods-----------------*/
public:
	/*Generate a cubic lattice of Nx x Ny x Nz particles with volume fraction eta_0*/
	void genCubic(int Nx, int Ny, int Nz, double eta_0, double viscosity);

	/*Initialize velocities at system temperature T*/
	void initVelocity(double Tinit);

	/*Store configuration in file*/
	void write_data(const char filename[]);

	/*Calculate force between particle and interface*/
	void Fint(int partID);

	/*Calculate force between particles due to WCA potential*/
	void FWCA();
	//same as FWCA but with verlet lists
	void FWCAverlet();

	///*Calculate the forces on the particles. Also calculates Epot*/
	void forces();

	/*Integrate one timestep dt*/
	void integrateVelocityVerlet(double dt);

	/*Integrate one timestep dt with brownian dynamics*/
	void integrateBrownianDynamics(double dt);

	/*Calculate all the distances between particles and generate a verlet list*/
	void VerletList(double rv);

	//calculates the average velocity
	double averageVelocity();
};


void Config::VerletList(double rv)
{	//r_v is the extra distance on top of the radii of the particles that are kept in the list

	//If the verlet is calculated again, initialise the number of particles in al the verletlists at 0
	for (int i=0;i<Ntot;i++)
	{
		config[i].nver=0;
	}
	//loop over all the particles to check there distances
	for (int i=0;i<Ntot-1;i++)
	{
		for (int j=i+1;j<Ntot;j++)
		{
				double distance[NDIM];	//Contains the distances in the x,y and z direction
				double r2;							//Squared the distance
				double r;								//The distance
				int n1,n2;              //The indices for the elements in the verlet list
				//calculate the distance
				distance[0] = config[i].coords[0] - config[j].coords[0];
				distance[1] = config[i].coords[1] - config[j].coords[1];
				distance[2] = config[i].coords[2] - config[j].coords[2];
				distance[0] -= (int)(2.0 * distance[0] / box[0]) * box[0];		//Nearest image
				distance[1] -= (int)(2.0 * distance[1] / box[1]) * box[1];		//Nearest image

				r2 = distance[0] * distance[0] + distance[1] * distance[1] + distance[2] * distance[2];
				r=sqrt(r2);
				//Correct the distance for the radii of the particles
				r-=hexaroot_2*(config[i].sigma+config[j].sigma);
				//If the distance is less than r_v append the the neighbouring atom to the verlet list
				if (r<rv)
				{
					n1=config[i].nver;
					n2=config[j].nver;
					//add the neighbour to both
					config[i].verlet[n1]=j;
					config[j].verlet[n2]=i;
					//Update the toal number of verlet neighbours for both
					config[i].nver=n1+1;
					config[j].nver=n2+1;
				}
			}
	}
}


void Config::genCubic(int Nx, int Ny, int Nz, double eta_0, double viscosity)
{
	double x, y, z;
	double PNs;		//Probability of small particle
	double sigma_l;
	int index;
	int Ns_actual = 0;

	mu = viscosity;

	Ntot = Nx * Ny * Nz;
	Ns = (int)((Nr / (Nr + 1.0)) * Ntot);
	PNs = (double)Ns / (double)Ntot;
	sigma_l = sigmar;

	/*Generate the particle array of the right size*/
	config = (Particle*)malloc(sizeof(Particle) * Ntot);
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
				config[index].coords[0] = x;
				config[index].coords[1] = y;
				config[index].coords[2] = z;
				if (dsfmt_genrand() < PNs) {
					config[index].sigma = 1.0;
					config[index].m = 1.0;
					config[index].kappa = 50;
					Ns_actual++;
				}
				else {
					config[index].sigma = sigma_l;
					config[index].m = pow(sigma_l, 3);
					config[index].kappa = 50 * sigmar * sigmar;
				}
				config[index].gamma = 3 * M_PI * config[index].sigma * mu;
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

void Config::write_data(const char filename[])
{
	double lbox;
	FILE* f = fopen(filename, "w");

	fprintf(f, "%d\n", Ntot);

	for (int i = 0; i < NDIM; i++) {
		lbox = box[i] / 2.0;
		fprintf(f, "%lf %lf\n", -lbox, lbox);
	}

	for (unsigned int i = 0; i < Ntot; i++) {
		for (int d = 0; d < NDIM; d++)
			fprintf(f, "%lf\t", config[i].coords[d] - box[d] / 2.0);
		if (config[i].sigma > 1.0)
			fprintf(f, "%lf\t%d\n", config[i].sigma, 1);
		else
			fprintf(f, "%lf\t%d\n", config[i].sigma, 2);
	}

	fclose(f);
}

void Config::initVelocity(double Tinit)
{
	double v2sum = 0.0;
	double vsum[NDIM] = { 0.0, 0.0, 0.0 };
	double vi;
	double fs;

	T = Tinit;

	for (int i = 0; i < Ntot; i++) {
		for (int d = 0; d < NDIM; d++) {
			vi = randd(-1.0, 1.0);
			config[i].v[d] = vi;
			v2sum += config[i].m * vi * vi;
			vsum[d] += vi;
		}
	}

	v2sum /= (double)Ntot;
	for (int d = 0; d < NDIM; d++) {
		vsum[d] /= (double)Ntot;
	}
	fs = sqrt(3.0 * Tinit / v2sum);
	for (int i = 0; i < Ntot; i++) {
		for (int d = 0; d < NDIM; d++) {
			config[i].v[d] = (config[i].v[d] - vsum[d]) * fs;
			config[i].v_init[d] = config[i].v[d];
			config[i].f[d] = 0.0;
		}
	}
}


double Config::averageVelocity()
{
	//The sum of all velocities
	double vtot=0.0;
	//Loop over all particles
	for (int i=0;i<Ntot;i++)
	{
		double vx=config[i].v[0];
		double vy=config[i].v[1];
		double vz=config[i].v[2];
		double vnet=sqrt(vx*vx+vy*vy+vz*vz);
		vtot+=vnet;
	}
	//Return the average velocity
	return vtot/Ntot;
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

	for (int i = 0; i < Ntot - 1; i++) {
		/*Particle - Particle forces*/
		for (int j = i + 1; j < Ntot; j++) {
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
					config[j].f[d] -= fradial * dist[d];
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
	/*Particle - floor forces and test particle - ceiling forces*/
	r2_floor = config[Ntot - 1].coords[2] * config[Ntot - 1].coords[2];
	s2 = config[Ntot - 1].sigma * config[Ntot - 1].sigma / 4.0;
	r2_cut = cuberoot_2 * s2;

	//config[Ntot - 1].f[2] -= config[Ntot - 1].m * 1.0; //gravity for test

	if (r2_floor < r2_cut) {
		r2i = 1.0 / r2_floor;
		r6i = pow(r2i, 3.0);
		s6 = pow(s2, 3.0);
		fradial = 48.0 * r2i * r6i * s6 * (s6 * r6i - 0.5);
		Epot += 4.0 * r6i * s6 * (s6 * r6i - 1.0) + 1.0;
		config[Ntot - 1].f[2] += fradial * config[Ntot - 1].coords[2];
	}

	Fint(Ntot - 1);
}


void Config::FWCAverlet()
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

	for (int i = 0; i < Ntot; i++) {
		int neighbours;
		neighbours=config[i].nver;
		/*Particle - Particle forces*/
		for (int k=0; k < neighbours; k++) {
			int j=config[i].verlet[k];
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
					//config[j].f[d] -= fradial * dist[d];
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
	/*Particle - floor forces and test particle - ceiling forces*/
	r2_floor = config[Ntot - 1].coords[2] * config[Ntot - 1].coords[2];
	s2 = config[Ntot - 1].sigma * config[Ntot - 1].sigma / 4.0;
	r2_cut = cuberoot_2 * s2;

	//config[Ntot - 1].f[2] -= config[Ntot - 1].m * 1.0; //gravity for test

	if (r2_floor < r2_cut) {
		r2i = 1.0 / r2_floor;
		r6i = pow(r2i, 3.0);
		s6 = pow(s2, 3.0);
		fradial = 48.0 * r2i * r6i * s6 * (s6 * r6i - 0.5);
		Epot += 4.0 * r6i * s6 * (s6 * r6i - 1.0) + 1.0;
		config[Ntot - 1].f[2] += fradial * config[Ntot - 1].coords[2];
	}

	Fint(Ntot - 1);
}


void Config::forces()
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
					config[j].f[d] -= fradial * dist[d];
				}
			}
		}

		/*Particle - floor forces and test particle - ceiling forces*/
		r2_floor = config[i].coords[2] * config[i].coords[2];
		s2 = (config[i].sigma + 1.0) * (config[i].sigma + 1.0) / 4.0;
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

		/*Random brownian forces and friction*/
		for (int d = 0; d < NDIM; d++) {
			config[i].f[d] += random_BoxMuller(sqrt(2 * config[i].gamma * T), 0.0);
			config[i].f[d] -= config[i].gamma * config[i].v[d];
		}
	}
}

void Config::Fint(int partID)
{
	if (config[partID].coords[2] > H)
		config[partID].f[2] -= config[partID].kappa * (config[partID].coords[2] - H);
}

void Config::integrateVelocityVerlet(double dt)
{
	double rnew;
	int nbox;

	Ekin = 0.0;

	for (int i = 0; i < Ntot; i++) {
		for (int d = 0; d < NDIM - 1; d++) {
			rnew = config[i].coords[d] + config[i].v[d] * dt + (config[i].f[d] / (config[i].m * 2.0)) * dt * dt;

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

		config[i].coords[2] += config[i].v[2] * dt + (config[i].f[2] / (config[i].m * 2.0)) * dt * dt;
	}

	forces();

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
	FWCAverlet();

	/*Euler integration*/
	double rnew;
	int nbox;
	for (int i = 0; i < Ntot; i++) {
		for (int d = 0; d < NDIM - 1; d++) {
			rnew = config[i].coords[d] + sqrt(2.0 * T * dt / config[i].gamma) * random_BoxMuller(sqrt(T * config[i].gamma), 0.0) + (dt / config[i].gamma) * config[i].f[d];

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

		config[i].coords[2] += sqrt(2.0 * T * dt / config[i].gamma) * random_BoxMuller(sqrt(T * config[i].gamma), 0.0) + (dt / config[i].gamma) * config[i].f[2];

	}
}


int main()
{
	dsfmt_seed(time(NULL));

	Config conf;
	conf.sigmar = 7.0;
	conf.Nr = 20.0;
	conf.genCubic(8, 8, 50, 0.04, 0.25);
	char buffer[100];

	/*conf.Nr = 1.0;
	conf.sigmar = 7.0;
	conf.Nl = 1;
	conf.Ns = 1;
	conf.Ntot = 2;
	conf.config = new Particle[2];
	conf.box[0] = 6;
	conf.box[1] = 6;
	conf.box[2] = 6;
	conf.config[0].coords[0] = 1;
	conf.config[0].coords[1] = 3;
	conf.config[0].coords[2] = 3;
	conf.config[1].coords[0] = 2;
	conf.config[1].coords[1] = 3;
	conf.config[1].coords[2] = 3;
	conf.config[0].m = 1;
	conf.config[1].m = 1;
	conf.config[0].sigma = 1;
	conf.config[1].sigma = 2;*/

	conf.initVelocity(1.0);
	/*
	conf.forces();
	printf("Nr: %.2lf\tNl: %d\tNs: %d\tNtot: %d\n", conf.Nr, conf.Nl, conf.Ns, conf.Ntot);
	printf("Epot: %lf\n", conf.Epot);
	int i = 0;
	while (conf.H > conf.box[2] / 2.0) {
		conf.integrateVelocityVerlet(0.01);
		if (i % 200 == 0) {
			sprintf(buffer, "snapshot_%d.dat", i);
			printf("snapshot: %d\tT: %lf\n", i, conf.Ekin * 2.0 / (3.0 * conf.Ntot));
			conf.write_data(buffer);
		}
		conf.H -= 0.1 * 0.01;
		i++;
	}
	*/

	conf.write_data("Initial.dat");
	int N_ver=10;//We update the verlet list every 5 steps
	double v_ave;//the average velocity
	double dt=0.001;
	v_ave=conf.averageVelocity();
	time_t t0=time(NULL);
	time_t t1;
	for (int i = 0; i < 1000000; i++) {
		if (i%N_ver==0)
		{
			conf.VerletList(v_ave*N_ver*3*dt);
		}
		conf.integrateBrownianDynamics(dt);
		if (i % 1000 == 0) {
			t1=time(NULL);
			printf("%d\n",t1-t0);
			sprintf(buffer, "snapshot_%d.dat", i);
			printf("snapshot: %d\n", i);
			conf.write_data(buffer);
		}
		conf.H -= 0.05 * 0.005;
	}

	conf.write_data("Final.dat");
	printf("Epot: %lf\n", conf.Epot);

	return 1;
}

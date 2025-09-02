#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include "mt19937.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif // M_PI

#define MIN(X,Y) (((X)<(Y)) ? (X):(Y))
//LATTICESIZE ONLY WORKS FOR THIS ANGLE
#define Latticesize 42
#define mc_steps 1000000
#define init_steps 0
#define out_steps 10000
#define NDIMspin 3
#define NDIM 3

/* Initialization variables */
//The array containing the spins
//L is the number of unit cells for the honeycomb lattice in both directions
int L=Latticesize;
int TweeL=2*Latticesize;
int n_part=2*Latticesize*Latticesize;
double Lattice[2*Latticesize][Latticesize][NDIMspin]={0.};
double Coordinates[2*Latticesize][Latticesize][NDIM]={0.};
//For each particles do the neighbours have both x and y coordinates
int Neighbourtable[2*Latticesize][Latticesize][3][2]={0};
//Rotated coordinates not needed here as it is already put into the calculation of the field
//double Rotatedcoordinates[(2*Latticesize*Latticesize-Latticesize)][NDIMspin]={0};
double Bfield[2*Latticesize][Latticesize]={0};
//Contains an array of wether the thirdneibhbour is up or down
bool Updownarray[2*Latticesize][Latticesize];

//The derivative of the x and y spin
double DSx[2*Latticesize][Latticesize][NDIM]={0.};
double DSy[2*Latticesize][Latticesize][NDIM]={0.};
//Array with the vorticity
double Vortarr[2*Latticesize][Latticesize][NDIM]={0.};
double Topcurr[2*Latticesize][Latticesize]={0.};
//Array with the correlation in the topological charge
double Chi_c[2*Latticesize][Latticesize]={0.};
//Arrays with the correlation for the spin in the z direction, and of those of the spins in the x-y direciton
double Chi_z[2*Latticesize][Latticesize]={0.};
double Chi_p[2*Latticesize][Latticesize]={0.};

//The parameters for our simulation. These give Ld=8.0, but these are the targetvalues after annealing
//ONLY CHANGE JEX AND A
double J_target=1.0;
double D_target=0.5;//0.205034;
double Jex_target=0.03;
//Ac is doorgaans ongeveer gelijk aan D^2/2J
double Ac_target=0.0;//0.125/2.;
//As variert qua sterkte
double As_target=0.0;
double kbT_target=0.01;

//The size of the steps
double alpha=1.0;

double theta=0.0280998; //This corresponds to 1.61 degrees

double acceptance=0.0;

//The reciprical lattice vectors for the moire super lattice
double G1[2]={0.};
double G2[2]={0.};
double G3[2]={0.};
//The reciprical lattice vectors for regular lattice
double recipvectorarray[6][2];

//vector between + and - in AFM lattice
double xi[2]={1./sqrt(3.),0.};

//Contains the topological charge, and currents
double Jmu[3]={0.};

//The vectors from one side to another
/*
double e1up[3]={sqrt(3.)/2.,0.5,0.};
double e1down[3]={sqrt(3.)/2.,-0.5,0.};
double e2up[3]={-sqrt(3.)/2.,0.5,0.};
double e2down[3]={-sqrt(3.)/2.,-0.5,0.};
double e3up[3]={0.,1.,0.};
double e3down[3]={0.,-1.,0.0};
*/

//double evec[2][3][3]={{{sqrt(3.)/2.,-0.5,0.},{-sqrt(3.)/2.,-0.5,0.},{0.,1.,0.}},{{sqrt(3.)/2.,0.5,0.},{-sqrt(3.)/2.,0.5,0.},{0.,-1.,0.0}}};
double evecup[3][3]={{sqrt(3.)/2.,-0.5,0.},{-sqrt(3.)/2.,-0.5,0.},{0.,1.,0.}};
double evecdown[3][3]={{sqrt(3.)/2.,0.5,0.},{-sqrt(3.)/2.,0.5,0.},{0.,-1.,0.0}};
//The unitvector in the z direction
double zvec[3]={0.,0.,1.};

//cross product
void Cross(double vec1[],double vec2[],double crossvec[])
{
    crossvec[0]=vec1[1]*vec2[2]-vec1[2]*vec2[1];
    crossvec[1]=vec1[2]*vec2[0]-vec1[0]*vec2[2];
    crossvec[2]=vec1[0]*vec2[1]-vec1[1]*vec2[0];
}

//inner product
double Dot(double vec1[],double vec2[])
{
    double inprod=0.;
    inprod=vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2];
    return inprod;
}

//matmul of a 2d matrix
void Matvec2d(double mat[2][2],double vec[2],double prod[2])
{
    prod[0]=mat[0][0]*vec[0]+mat[0][1]*vec[1];
    prod[1]=mat[1][0]*vec[0]+mat[1][1]*vec[1];
}

//The rotation matrix
void Rmat(double theta,double res[2][2])
{
    res[0][0]=cos(theta);
    res[1][1]=cos(theta);
    res[0][1]=-sin(theta);
    res[1][0]=sin(theta);
}

void Buildrandomlattice()
{   //Builds a lattice with a random configuration

    for(int i=0;i<(TweeL);++i)
    {
        for(int j=0;j<L;++j)
            {
                double theta=0.;
                double phi=0.;

                theta=2.0*M_PI*dsfmt_genrand();
                phi=M_PI*dsfmt_genrand();

                Lattice[i][j][0]=sin(phi)*cos(theta);
                Lattice[i][j][1]=sin(phi)*sin(theta);
                Lattice[i][j][2]=cos(phi);
                Lattice[i][j][3]=theta;
                Lattice[i][j][4]=phi;
            }
    }
}

void Buildfmlattice()
{
    //Builds a lattice with all spins pointing upwards
    for(int i=0;i<(TweeL);++i)
    {
        for(int j=0;j<L;++j)
        {
            double theta=0.;
            double phi=0.;

            theta=0.0;
            phi=0.0;
            Lattice[i][j][0]=0.0;
            Lattice[i][j][1]=0.0;
            Lattice[i][j][2]=1.0;
            Lattice[i][j][3]=theta;
            Lattice[i][j][4]=phi;
        }
    }
}

void Buildfmylattice()
{
    //Builds a lattice with all spins pointing in the +y direction
    for(int i=0;i<(TweeL);++i)
    {
        for(int j=0;j<L;++j)
        {
            double theta;
            double phi;

            theta=M_PI/2,0;
            phi=M_PI/2.0;

            Lattice[i][j][0]=0.0;
            Lattice[i][j][1]=1.0;
            Lattice[i][j][2]=0.0;
            Lattice[i][j][3]=theta;
            Lattice[i][j][4]=phi;
        }
    }

}


void Loadcoordinates(char filename[128])
{
  //printf("(opening) filename = %s\n", filename);

  FILE* initfile = fopen (filename, "r");				//Open the file

  for (int i = 0; i < TweeL; i++)
    for(int j=0;j<L;++j)
  {
      {
        int test = fscanf (initfile, "%lf\t%lf\t%lf\t%lf\t%lf\n", &(Lattice[i][j][0]), &(Lattice[i][j][1]), &(Lattice[i][j][2]), &(Lattice[i][j][3]), &(Lattice[i][j][4]));
        if (test < 5) { printf ("Problem reading particles!\n"); exit(3);}

      }
  }

  fclose (initfile);

}


void Buildlattice(int latticetype,char namelowest[128])
{
    //builds an initial lattice, the given integer tells us which type of lattice we want
    //0 gives a random lattice, 1 a ferromagnetic lattice, and 2 a lattice with all the spins in the y direction
    //all other numbers will also give an FM lattice
    char nameskyr[128]="Skyrmion_config.dat";
    char namespiral[128]="Spiral_config.dat";
    //!TODO DEZE BESTCONFIG NOG FF FIXEN
    //char namelowest[128]="Best_config.dat";
    switch(latticetype)
    {
    case 0:
        Loadcoordinates(namelowest);
    case 1:
        Buildrandomlattice();
        break;
    case 2:
        Buildfmlattice();
        break;
    case 3:
        Loadcoordinates(nameskyr);
        break;
    case 4:
        Loadcoordinates(namespiral);
        break;

    default:
        Buildfmlattice();
        break;
    }
    //printf("Initial lattice is build! \n");
}


void Calculate_neighbourtable()
{
    //Construct the table with all the neighbours
    for(int i=0;i<(TweeL);++i)
    {
        for(int j=0;j<L;++j)
        {
            //The neighbour to the right
            Neighbourtable[i][j][0][0]=(i+TweeL+1)%(TweeL);
            Neighbourtable[i][j][0][1]=j;
            //The neighbour to the left
            Neighbourtable[i][j][1][0]=(i+TweeL-1)%(TweeL);
            Neighbourtable[i][j][1][1]=j;

            //wanneer i even is en j ook is het boven
            if(j%2==0)
            {
                if(i%2==0)
                {
                    //The neighbour is above
                    Neighbourtable[i][j][2][0]=i;
                    Neighbourtable[i][j][2][1]=(j+1+L)%L;
                    Updownarray[i][j]=true;
                }
                else
                {
                    //The neighbour is below
                    Neighbourtable[i][j][2][0]=i;
                    Neighbourtable[i][j][2][1]=(j-1+L)%L;
                    Updownarray[i][j]=false;
                }

            }
            else
            {
                if(i%2!=0)
                {
                    //The neighbour is above
                    Neighbourtable[i][j][2][0]=i;
                    Neighbourtable[i][j][2][1]=(j+1+L)%L;
                    Updownarray[i][j]=true;
                }
                else
                {
                    //The neighbour is below
                    Neighbourtable[i][j][2][0]=i;
                    Neighbourtable[i][j][2][1]=(j-1+L)%L;
                    Updownarray[i][j]=false;
                }

            }

        }
    }

}


void Generate_positions()
{
    for(int i=0;i<(TweeL);++i)
    {
        for(int j=0;j<L;++j)
        {
            double x;
            double y;
            if(j%2==0)
            //!BIETJE OMSLAGTIG
            {
                x=(double)i*sqrt(3)/2.0;
                y=(double)j*1.5-0.5*((double)(i%2));
            }
            else
            {
                x=(double)i*sqrt(3)/2.0;
                y=((double)j-1.0)*1.5+1.0+0.5*((double)(i%2));
            }
            Coordinates[i][j][0]=x;
            Coordinates[i][j][1]=y;
        }
    }
}


void Generate_recip_vectors()
{
    //!TODO, DEZE FACTOR MISSCHIEN PRECIEZER
    G1[1]=4.*M_PI/(sqrt(3.)*36.373067);
    double xrot=4.*M_PI/(sqrt(3.)*36.373067)*-sin(M_PI/3.); //FACTOR GEDEELD DOOR 35 WEER TOEVOEGEN
    double yrot=4.*M_PI/(sqrt(3.)*36.373067)*cos(M_PI/3.); //FACTOR GEDEELD DOOR 35 WEER TOEVOEGEN
    G2[0]=xrot;
    G3[0]=xrot;
    G2[1]=yrot;
    G3[1]=-yrot;

}

void Calculate_B_field()
{
    for(int i=0;i<(TweeL);++i)
    {
        for(int j=0;j<L;++j)
        {
            double B=0.;
            //The coordinates
            double x=Coordinates[i][j][0];
            double y=Coordinates[i][j][1];
            //EQUATION A1 FROM THE HIERARCY PAPER
            //WEER COS VAN MAKEN
            //!WEER OPLETTEN, DE FACTOR JEX ZIT NU DUS IN CALCULATEMAILTONIAN ZODAT HET MAKKELIJKER IS MET DE ANNEALING
            B=(sin(x*G1[0]+y*G1[1])+sin(x*G2[0]+y*G2[1])+sin(x*G3[0]+y*G3[1]));

            //The reciprical lattice vectors for a nonrotated system
            //B+=2.0*Jex*(cos(x*(G1[0]-xi[0])+y*(G1[1]-xi[1]))+cos(x*(G2[0]-xi[0])+y*(G2[1]-xi[1]))+cos(x*(G3[0]-xi[0])+y*(G3[1]-xi[1])));
            //B+=2.0*Jex*(cos(x*(G1[0]+xi[0])+y*(G1[1]+xi[1]))+cos(x*(G2[0]+xi[0])+y*(G2[1]+xi[1]))+cos(x*(G3[0]+xi[0])+y*(G3[1]+xi[1])));
            Bfield[i][j]=B;
        }
    }

}


double Calculatehamiltonian(int Latti,int Lattj,double J, double Jex, double Ac,double As, double D)
{
    //Calculates the Hamiltonian, it furthermore has the parameters J, Jex, A and D as optional arguments
    double Energy=0.0;
    //Contains both the x and y coordinates
    int Nb[3][2]={0};

    //The spin of our lattice site
    double Spin[3];
    for(int i=0;i<3;++i)
    {
        Spin[i]=Lattice[Latti][Lattj][i];
    }

    Nb[0][0]=Neighbourtable[Latti][Lattj][0][0]; //Right neighbour
    Nb[0][1]=Neighbourtable[Latti][Lattj][0][1];
    Nb[1][0]=Neighbourtable[Latti][Lattj][1][0]; //left neigbhour
    Nb[1][1]=Neighbourtable[Latti][Lattj][1][1];
    Nb[2][0]=Neighbourtable[Latti][Lattj][2][0];//up or down neighbour
    Nb[2][1]=Neighbourtable[Latti][Lattj][2][1];

    //Calculate the energy from the different therms
    //Energy from Heisenberg coupling (divided by two because it is calculated for two particles)
    for(int i=0;i<3;++i)
    {
        Energy-=2*J*Dot(Spin,Lattice[Nb[i][0]][Nb[i][1]]);
    }

    //ALL THE OTHER STUFF WITH CROSS VECTORS ETC IS NOT NEEDED HERE
    //but we still need the anisitropy therm and the coupling with the other field
    Energy-=As*Spin[2]*Spin[2];
    //!ONTHOUDEN, HIER STAAT DE FACTOR JEX BUITEN HET B-VELD, BIJ DE VORIGE CODE ERIN
    Energy-=Jex*Bfield[Latti][Lattj]*Spin[2];
    double evecupdown[3][3];
    //Calculate the DMI coupling, first we see which of the latticevectors we need by checking if it has a neighbour below or above
    //!TODODODODODOD DEZE DIT IS ZO DOMMMMM
    //Chech which vectors are needed to describe the orientatins
    //Also calculate the energy due to the compass anisitropy with it
    if(Updownarray[Latti][Lattj]==true)
    {
        //now calculate the real DMI coupling
        double dvec[3]={0.};
        double spincross[3]={0.};
        for(int i=0;i<3;++i)
        {
            //calculate d
            Cross(zvec,evecup[i],dvec);
            //The neighbour spin is:
            int Nbi=Nb[i][0];
            double Nbj=Nb[i][1];
            //Get the spin of the neighbour
            Cross(Spin,Lattice[Nb[i][0]][Nb[i][1]],spincross);
            double dotspin=0.;
            dotspin=Dot(dvec,spincross);
            //This is the DMI therm
            Energy-=2.0*D*dotspin;
            //This is the Compass therm
            Energy-=2.0*Ac*Dot(Spin,dvec)*Dot(Lattice[Nb[i][0]][Nb[i][1]],dvec);
        }

    }
    else
    {
       //now calculate the real DMI coupling
        double dvec[3]={0.};
        double spincross[3]={0.};
        for(int i=0;i<3;++i)
        {
            //calculate d
            Cross(zvec,evecdown[i],dvec);
            //The neighbour spin is:
            int Nbi=Nb[i][0];
            double Nbj=Nb[i][1];
            //Get the spin of the neighbour

            Cross(Spin,Lattice[Nb[i][0]][Nb[i][1]],spincross);
            double dotspin=0.;
            dotspin=Dot(dvec,spincross);
            //This is the DMI therm
            Energy-=2.0*D*dotspin;
            //This is the Compass therm
            Energy-=2.0*Ac*Dot(Spin,dvec)*Dot(Lattice[Nb[i][0]][Nb[i][1]],dvec);
        }
    }

    return Energy;

}



void Flipspins(double kbT,double J,double Jex, double Ac,double As, double D)
{
    //BECAUSE IT IS NOW SIMPLER, WE ONLY NEED THE HAMILTONIAN OF ONE PARTICLE

    //loop over all the particles
    for(int i=0;i<n_part;++i)
    {
        int Ri;
        int Rj;
        Ri=(int)(TweeL*dsfmt_genrand());
        Rj=(int)(L*dsfmt_genrand());

        //de vectoren voor de willekeurige verandering
        double rx;
        double ry;
        double rz;

        double xikwad=2.;

        while(xikwad>=1.)
        {
            //willekeurige getallen op (0,1) zie pagina 529 van de pdf vna simulaitons fo fluids
            double ran1=dsfmt_genrand();
            double ran2=dsfmt_genrand();
            double xi1=2.*ran1-1.;
            double xi2=2.*ran2-1.;
            xikwad=xi1*xi1+xi2*xi2;
            rx=2.*xi1*sqrt(1.-xikwad);
            ry=2.*xi2*sqrt(1.-xikwad);
            rz=1-2.*xikwad;
        }
        /*
        if(i==1)
        {
            printf("%lf\t%lf\t%lf\n",rx,ry,rz);
        }
        */

        /*
        double Thetanew;
        double Phinew;

        //Get the new angles of the spins
        Thetanew=fmod(Lattice[id][Ri][Rj][3]+alpha*M_PI*(2.0*dsfmt_genrand()-1.0),2.0*M_PI);
        Phinew=Lattice[id][Ri][Rj][4]+alpha*M_PI*(2.0*dsfmt_genrand()-1.0);


        //Making sure that the angles are between 0,pi and 0,2pi
        if(Thetanew<0.)
        {
            Thetanew=2.0*M_PI-Thetanew;
        }
        if(Phinew>M_PI)
        {
            Phinew=M_PI-fmod(Phinew,M_PI);
            Thetanew=fmod(Thetanew+M_PI,2.0*M_PI);
        }
        else if(Phinew<0.)
        {
            Phinew*=-1.;
            Thetanew=fmod(Thetanew+M_PI,2.0*M_PI);
        }
        */

        double Oldspin[3]={0.};
        Oldspin[0]=Lattice[Ri][Rj][0];
        Oldspin[1]=Lattice[Ri][Rj][1];
        Oldspin[2]=Lattice[Ri][Rj][2];

        double nx=Oldspin[0]+alpha*rx;
        double ny=Oldspin[1]+alpha*ry;
        double nz=Oldspin[2]+alpha*rz;

        double n=nx*nx+ny*ny+nz*nz;

        /*
        if(i==1)
        {
            printf("%lf\t%lf\t%lf\t%lf\n",nx,ny,nz,n);
        }
        */

        //maak er weer eenheidsvectoren van
        nx=nx/sqrt(n);
        ny=ny/sqrt(n);
        nz=nz/sqrt(n);

        /*
        if(i==1)
        {
            printf("%lf\t%lf\t%lf\t%lf\n",nx,ny,nz,n);
        }
        */

        //NU DE SPINS UPDATEN NAAR DE NIEUWE VARIANT
        Lattice[Ri][Rj][0]=nx;
        Lattice[Ri][Rj][1]=ny;
        Lattice[Ri][Rj][2]=nz;


        //calculate the energy after the flip

        Enew=Calculatehamiltonian(Ri,Rj,J,Jex,Ac,As,D);

        double Boltzmanweight=0.;
        double randomnumb=0.;

        //check if the change is accepted or not

        Boltzmanweight=MIN(1.0,exp(-(Enew-Eold)/kbT));
        randomnumb=dsfmt_genrand();

        if(randomnumb<Boltzmanweight)
        {
            acceptance+=1.0/n_part;
        }
        else
        {
            Lattice[Ri][Rj][0]=Oldspin[0];
            Lattice[Ri][Rj][1]=Oldspin[1];
            Lattice[Ri][Rj][2]=Oldspin[2];
        }

    }
}

//!TODO KLOPT NIET MEEEEEEEEERRRRR
void Write_data(int Latticetype,double J,double Jex, double Ac,double As, double D,double Etot)
{
    //Latticetypes work the same as before. This code is meant to only return the final state
    //The filename will consists out of the type, J, Jex, A and D and the first line in the results will contain
    //also all these, but also the total Energy


    //contains the filename
    char buffer[128];
    //makes the filename depend on the paramers
    switch(Latticetype)
    {
        case 0:
            sprintf(buffer,"Periodic_Random_J_%.2f_Jex_%.2f_Ac_%.2f_As_%.2f_D_%.2f.dat",J,Jex,Ac,As,D);
            break;
        case 1:
            sprintf(buffer,"Periodic_Up_J_%.2f_Jex_%.2f_Ac_%.2f_As_%.2f_D_%.2f.dat",J,Jex,Ac,As,D);
            break;
        case 2:
            sprintf(buffer,"Periodic_Skyr_J_%.2f_Jex_%.2f_Ac_%.2f_As_%.2f_D_%.2f.dat",J,Jex,Ac,As,D);
            break;
        case 3:
            sprintf(buffer,"Periodic_SP_J_%.2f_Jex_%.2f_Ac_%.2f_As_%.2f_D_%.2f.dat",J,Jex,Ac,As,D);
            break;
        default:
            sprintf(buffer,"HIERISIETSMISGEGAAN_J_%.2f_Jex_%.2f_Ac_%.2f_As_%.2f_D_%.2f.dat",J,Jex,Ac,As,D);
            break;
    }

    FILE* fp=fopen(buffer,"w");
    //First also print lattice type and all the rest on the first row
    fprintf(fp,"%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",Latticetype,J,Jex,Ac,As,D,Etot);


    //Now print all the coordinates
    for(int i=0;i<(TweeL);++i)
    {
        for(int j=0;j<L;++j)
            {
            double x=Coordinates[i][j][0];
            double y=Coordinates[i][j][1];
            double z=Coordinates[i][j][2];
            fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",x,y,z,Lattice[i][j][0],Lattice[i][j][1],Lattice[i][j][2]);
            }

    }
    fclose(fp);
}

void Lattice_Saver(char name[128])
{
    FILE* fp=fopen(name,"w");
    for(int i=0;i<(TweeL);++i)
    {
        for(int j=0;j<L;++j)
        {
            double Sx=Lattice[i][j][0];
            double Sy=Lattice[i][j][1];
            double Sz=Lattice[i][j][2];
            double Theta=Lattice[i][j][3];
            double Phi=Lattice[i][j][4];
            fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\n",Sx,Sy,Sz,Theta,Phi);
        }
    }
    fclose(fp);
}

void Temperature_annealing(int steps, int mcsteps,double kbT_init,double J, double Jex,double Ac,double As, double D)
{
    /*Perform temperature annealing, by starting at a high temperature and then cooling down
    //steps is the number of steps during which the kbt needs to be lowered, and mc_steps is the number of sweeps for each step
    //kbT the initial kbt and the target kbt is already defined above in the script
    J etc also needs to be specified as this leaves the option open for both temperature and parameter annealing
    */
    //! mc_steps is a constant, but mcsteps is what goes into this function OP BLIJVEN LETTEN
    double kbT;
    for(int i=0;i<steps;++i)
    {
        double acceptratio=0.;
        kbT=kbT_init-(kbT_init-kbT_target)*(double)i/((double)steps);
        for(int j=0;j<mcsteps;++j)
        {
            Flipspins(kbT,J,Jex,Ac,As,D);
        }
        acceptratio=acceptance/((double)mcsteps);

        if(acceptratio<0.4)
            {
                if(alpha>0.11)alpha-=0.1;
                else if(alpha>0.02)alpha-=0.01;
            }
        if((acceptratio>0.6)&&(alpha<1.0))alpha+=0.1;

        acceptance=0.;
    }
    //reset alpha again
    alpha=1.0;

}


void Parameter_annealing(int steps,int mcsteps,double kbT_init, double kbT_tar,double J_init,double J_tar,double Jex_init, double Jex_tar,double Ac_init, double Ac_tar,double As_init, double As_tar,double D_init, double D_tar)
{
    /*
    This function performs parameter annealing on all the parameters. For each parameter you must supply the target value and the initial value and the target value
    If you use tar=init there will ofcouce be no annealing for this parater
    */
    double kbT,J,Jex,Ac,As,D;
    for(int i=0;i<steps;++i)
    {
        double acceptratio=0.;
        kbT=kbT_init-(kbT_init-kbT_tar)*(double)i/((double)steps);
        J=J_init-(J_init-J_tar)*(double)i/((double)steps);
        Jex=Jex_init-(Jex_init-Jex_tar)*(double)i/((double)steps);
        Ac=Ac_init-(Ac_init-Ac_tar)*(double)i/((double)steps);
        As=As_init-(As_init-As_tar)*(double)i/((double)steps);
        D=D_init-(D_init-D_tar)*(double)i/((double)steps);
        for(int j=0;j<mcsteps;j++)
        {
            Flipspins(kbT,J,Jex,Ac,As,D);
        }
        acceptratio=acceptance/((double)mcsteps);

        if(acceptratio<0.4)
            {
                if(alpha>0.11)alpha-=0.1;
                else if(alpha>0.02)alpha-=0.01;
            }
        if((acceptratio>0.6)&&(alpha<1.0))alpha+=0.1;

        acceptance=0.;
    }
    //reset alpha again
    alpha=1.0;
}

//HIER VEEL CODE VOOR INTERPERTATIE VAN DE DATA
//Ik heb nodig: magnetisme, afgeleides van de spin, vorticiteit, topologische lading,Chi_C, Chi_z and Chiplus

double Calculate_magnetism()
{
    //calculates the total magnetism of the lattice
    double magnetism=0.;
    for(int i=0;i<(2*Latticesize);++i)
    {
        for(int j=0;j<Latticesize;++j)
        {
            magnetism+=Lattice[i][j][2]/(double)n_part;
        }
    }
    return magnetism;
}

//!TODODODDO FACTOR 2 NIET VERGETEN
double Calculate_E_tot(double J, double Jex, double Ac,double As, double D)
{
    double Etot=0.0;
    for(int i=0;i<(TweeL);++i)
    {
        for(int j=0;j<L;++j)
        {
            Etot+=Calculatehamiltonian(i,j,J,Jex,Ac,As,D)/2.;
        }
    }
    return Etot;
}


void Calculate_Derivatives()
{
    //Loop over all the particles and calculate the derivates of the spins
    double dx=sqrt(3.)/2.;
    double dy1=1.0;
    double dy2=0.5;
    for(int i=0;i<TweeL;++i)
    {
        for(int j=0;j<L;++j)
        {
            int Nb[3][2]={0};
            //Get the neighbours
            Nb[0][0]=Neighbourtable[i][j][0][0]; //Right neighbour
            Nb[0][1]=Neighbourtable[i][j][0][1];
            Nb[1][0]=Neighbourtable[i][j][1][0]; //left neigbhour
            Nb[1][1]=Neighbourtable[i][j][1][1];
            Nb[2][0]=Neighbourtable[i][j][2][0];//up or down neighbour
            Nb[2][1]=Neighbourtable[i][j][2][1];
            //The spins of the neighbours,the first one says which neighbour it is, and the second what the value of that spincomponent is
            double Sn[3][3];
            for (int k=0;k<3;++k)
            {
                for(int l=0;l<3;++l)
                {
                    Sn[k][l]=Lattice[Nb[k][0]][Nb[k][1]][l];
                }
            }

            for(int k=0;k<3;++k)
            {
                //Calculate the x derivative
                DSx[i][j][k]=(Sn[1][k]-Sn[0][k])/(2.0*dx);
            }
            //Calculate the y derivative
            if (Updownarray[i][j]==true)
            {
                for(int k=0;k<3;++k)
                {
                    //calculate the y derivative
                    DSy[i][j][k]=(Sn[2][k]-(Sn[1][k]+Sn[0][k])/2.0)/(dy1+dy2);
                }
            }
            else
            {
                for(int k=0;k<3;++k)
                {
                    //Here is an min in front as in this case the neighbour is below
                    DSy[i][j][k]=-(Sn[2][k]-(Sn[1][k]+Sn[0][k])/2.0)/(dy1+dy2);
                }
            }


        }
    }
}

void Calculate_Vorticity()
{
    for (int i=0;i<TweeL;++i)
    {
        for(int j=0;j<L;++j)
        {
            double vortvec[3];
            Cross(DSx[i][j],DSy[i][j],vortvec);
            for (int k=0;k<3;++k)
            {
                Vortarr[i][j][k]=vortvec[k];
            }
        }
    }
}

void Calculate_Topcurr()
{
    for(int i=0;i<TweeL;++i)
    {
        for(int j=0;j<L;++j)
        {
            //!NIET HELEMAAL NETJES MAAR WERKT ALS HET GOED IS
            Topcurr[i][j]=Dot(Lattice[i][j],Vortarr[i][j]);
        }
    }
}


void Calculate_correlations()
{
    //This function calculates the various calculation functions
    for(int i=0;i<TweeL;++i)
    {
        for(int j=0;j<L;++j)
        {
            double rho;
            double Sx;
            double Sy;
            double Sz;
            double Chiz=0.;
            double Chic=0.;
            double Chip=0.;
            rho=Topcurr[i][j];
            Sx=Lattice[i][j][0];
            Sy=Lattice[i][j][1];
            Sz=Lattice[i][j][2];

            for(int k=0;k<TweeL;++k)
            {
                for(int l=0;l<L;++l)
                {
                    Chic+=rho*Topcurr[k][l];
                    Chiz+=Sz*Lattice[k][l][2];
                    Chip+=(Sx*Lattice[k][l][0]-Sy*Lattice[k][l][1]);
                }
            }
            Chi_c[i][j]=Chic;
            Chi_z[i][j]=Chiz;
            Chi_p[i][j]=Chip;
        }
    }
}

double *Calculate_Properties(double *results)
{
    //This function calculates the properties of the lattice, which includes the magnetization, avarege topological charge and the averages of the correlations

    //calculate the average magnetization, and takes the absolute value because of invariance
    Calculate_Derivatives();
    Calculate_Vorticity();
    Calculate_Topcurr();
    Calculate_correlations();
    results[0]=fabs(Calculate_magnetism());
    double rho=0.0;
    double Chic=0.0;
    double Chiz=0.0;
    double Chip=0.0;
    double dnpart=(double)n_part;
    for(int i=0;i<TweeL;++i)
    {
        for(int j=0;j<L;++j)
        {
            rho+=Topcurr[i][j]/dnpart;
            Chic+=Chi_c[i][j]/dnpart;
            Chiz+=Chi_z[i][j]/dnpart;
            Chip+=Chi_p[i][j]/dnpart;
        }
    }
    results[1]=rho;
    results[2]=Chic;
    results[3]=Chiz;
    results[4]=Chip;
    return results;
}

int main()
{
    //!A SCRIPT FOR BASIS TEMPERATURE ANNEALING

    //!DEZE FF DEFINIEREN
    int steps=100;
    int mcsteps=500;

    dsfmt_seed(time(NULL));
    //loop over all the possible startlattices
    //These can be recycled every time, which is very nice
    Generate_positions();
    Generate_recip_vectors();
    Calculate_B_field();
    Calculate_neighbourtable();

    //Do the temperature annealing from kbt=0.2 to kbt=0.01
    double kbT_init=0.2;

    double J=J_target;
    double Ac=Ac_target;
    double As=As_target;

    double J_ex=Jex_target;
    //!TODO HIER WEER resultsarray[26][6]={0.}; van maken
    double resultsarray[51][6]={0.};
    double resultsarray2[51][6]={0.};
    char namelowest[128];
    sprintf(namelowest,"Best_config.dat",J_ex);
    //!todo hier ook weer 26
    for(int i=0;i<51;++i)
    {
        //Calculate the value of D
        double D;
        D=(double)i*D_target/50.;


        //The lowest energy found thus far for this configuration
        double E_lowest=0.;
        for(int startlattice=0;startlattice<5;++startlattice)
        {
            //Build the startlattice
            Buildlattice(startlattice,namelowest);

            //Perform the temperature annealing
            Temperature_annealing(steps,mcsteps,kbT_init,J,J_ex,Ac,As,D);

            double Etot;

            //Calculate the energy and shit of the lattice
            Etot=Calculate_E_tot(J,J_ex,Ac,As,D);
            double results[5]={0.};

            Calculate_Properties(results);

            if(Etot<E_lowest)
            {
                //If this is the best lattice thus far, save it so it can be used for furhter calculations
                //!TODO DE BEST CONFIG FF FIXEN DAT DAT JUUST GOAT
                Lattice_Saver(namelowest);
                //also safe all the caratheristics of the lattice
                //!TODO HIER WEER resultsarray[i][0]=Etot; van maken
                resultsarray[i][0]=Etot;
                for(int j=1;j<6;++j)
                {
                    //!HIER OOK resultsarray[i][j]
                    resultsarray[i][j]=results[j-1];
                }
            }


        }
        char gavenaam2[128];
        sprintf(gavenaam2,"Struc_D=%.2f_Jex=%.2f.dat",D,J_ex);
        Lattice_Saver(gavenaam2);
    //!TODO DIT HAAKJE NIETVERGETEN
    }

    for(int i=0;i<51;++i)
    {
        //Calculate the value of D
        double D;
        D=D_target-(double)i*D_target/50.;


        //The lowest energy found thus far for this configuration
        double E_lowest=0.;
        for(int startlattice=0;startlattice<5;++startlattice)
        {
            //Build the startlattice
            Buildlattice(startlattice,namelowest);

            //Perform the temperature annealing
            Temperature_annealing(steps,mcsteps,kbT_init,J,J_ex,Ac,As,D);

            double Etot;

            //Calculate the energy and shit of the lattice
            Etot=Calculate_E_tot(J,J_ex,Ac,As,D);
            double results[5]={0.};

            Calculate_Properties(results);

            if(Etot<E_lowest)
            {
                //If this is the best lattice thus far, save it so it can be used for furhter calculations
                //!TODO DE BEST CONFIG FF FIXEN DAT DAT JUUST GOAT
                Lattice_Saver(namelowest);
                //also safe all the caratheristics of the lattice
                //!TODO HIER WEER resultsarray[i][0]=Etot; van maken
                resultsarray2[i][0]=Etot;
                for(int j=1;j<6;++j)
                {
                    //!HIER OOK resultsarray[i][j]
                    resultsarray2[i][j]=results[j-1];
                }
            }


        }
        char gavenaam[128];
        sprintf(gavenaam,"Struc_Jex=%.2f_D=%.2f.dat",J_ex,D);
        Lattice_Saver(gavenaam);
    //!TODO DIT HAAKJE NIETVERGETEN
    }

    //Safe all the results to a file
    char resultsname[128];
    char resultsname2[128];

    sprintf(resultsname,"Poging3_Jex_%.2f_annealing_results.dat",J_ex);
    sprintf(resultsname2,"Poging4_Jex_%.2f_annealing_results.dat",J_ex);
    //!TODO, OOK WEER GEWOON TERUGVOEGE
    //!todo hier ook weer 26
    FILE* fp=fopen(resultsname,"w");
    for(int i=0;i<51;i++)
    {
        fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",resultsarray[i][0],resultsarray[i][1],resultsarray[i][2],resultsarray[i][3],resultsarray[i][4],resultsarray[i][5]);
    }

    fclose(fp);

    FILE* fp2=fopen(resultsname2,"w");
    for(int i=0;i<51;i++)
    {
        fprintf(fp2,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",resultsarray2[i][0],resultsarray2[i][1],resultsarray2[i][2],resultsarray2[i][3],resultsarray2[i][4],resultsarray2[i][5]);
    }

    fclose(fp2);


    return 0;
}

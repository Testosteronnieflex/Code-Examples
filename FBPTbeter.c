#include <stdio.h>
#include <time.h>
//#include <omp.h>
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include "mt19937.h"

#define NTHREADS 45 //The number of cores the code is ran on 90

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif // M_PI


#define MIN(X,Y) (((X)<(Y)) ? (X):(Y))
//LATTICESIZE ONLY WORKS FOR THIS ANGLE
#define Latticesize 42
#define TweeL 84
//!EVENTUEEL WEER GOED ZETTEN
#define n_part 1764
#define NDIMspin 3
#define NDIM 3

//De hoeveelheid MC stappen
#define MC_steps 1000000 ///factor 10 drbij
#define init_steps 500000 ///factor 10 erbij

/* Initialization variables */
//First the elements that are constants
double J=1.0;
//double D=0.20;
//double Jex=0.20;
//Ac is doorgaans ongeveer gelijk aan D^2/2J
//double Ac=0.0;
//double As=0.;

//double theta=0.0280998;
//double alpha=0.2;//initial value of alpha for all simulations

char nameskyr[128]="Skyrmion_config.dat";

int L=Latticesize;

//Now the arrays that are constant for all grids
double Coordinates[Latticesize][Latticesize][NDIM]={0.}; //contains the xy coordinates for all points
int Neighbourtable[Latticesize][Latticesize][3][2]={0}; //contains all the neigbhrous of all the particles
double Bfield[Latticesize][Latticesize]={0}; //The Bfield for all the lattice points
//Contains an array of wether the thirdneibhbour is up or down
bool Updownarray[Latticesize][Latticesize];

//De moiréroosterconstante
double a_moir=36.373067;

//The reciprical lattice vectors for the moire super lattice
double G1[2]={0.};
double G2[2]={0.};
double G3[2]={0.};
//The reciprical lattice vectors for regular lattice
double recipvectorarray[6][2];

//vector between + and - in AFM lattice
double xi[2]={1./sqrt(3.),0.};
//contains the vectors in the case there is a particle in either the +y direction or -y direction
double evecup[3][3]={{sqrt(3.)/2.,-0.5,0.},{-sqrt(3.)/2.,-0.5,0.},{0.,1.,0.}};
double evecdown[3][3]={{sqrt(3.)/2.,0.5,0.},{-sqrt(3.)/2.,0.5,0.},{0.,-1.,0.0}};
//The unitvector in the z direction
double zvec[3]={0.,0.,1.};

//This array contains which core does which temperature
int Lattice_index[NTHREADS]={0};
//An array with the temperatues for each run
double Temperature_array[NTHREADS]={0.};
//Array with the temperature of all grids
double E_tot_array[NTHREADS]={0.};
//Array with modified alpha's
double Alpha_array[NTHREADS]={0.3};
//Array with all the values of the magnetisation
//double Mag_res[NTHREADS][MC_steps];

//Contains all the data
double Lattice[NTHREADS][Latticesize][Latticesize][NDIMspin]={0.};

//Een ding met de histogrammen van de stroming van de deeltjes
double nup[NTHREADS]={0.};
double ndown[NTHREADS]={0.};

bool Nup_arr[NTHREADS];
bool Ndown_arr[NTHREADS];

//Dit is de function ft en etat
double ft[NTHREADS];
double etat[NTHREADS-1];

double DT_arr[NTHREADS-1];

double Temperature_array_mod[NTHREADS]={0.};

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

//Generate the Lattice index array
void Generate_Lattice_index()
{
    for(int i=0;i<NTHREADS;++i)
    {
        Lattice_index[i]=i;
    }
}

//generate the temperatures for each thread, the input are the minum temperature and the factor that is the difference between the two
void Generate_Temperatures_Geometric(double T_min, double factor)
{
    double T=T_min;
    for(int i=0;i<NTHREADS;++i)
    {
        Temperature_array[i]=T;
        T*=factor;

    }
}

//!JA HALLOOOOOOO
void Generate_Temperatures_Linear(double T_min,double dT)
{
    double T=T_min;
    for(int i;i<NTHREADS;++i)
    {
        Temperature_array[i]=T;
        T+=dT;
    }
}

//Berekend wat van elke positie de buren zijn
void Calculate_neighbourtable()
{
    //Construct the table with all the neighbours
    for(int i=0;i<(Latticesize);++i)
    {
        for(int j=0;j<Latticesize;++j)
        {
            //The neighbour to the right
            Neighbourtable[i][j][0][0]=(i+Latticesize+1)%(Latticesize);
            Neighbourtable[i][j][0][1]=j;
            //The neighbour to the left
            Neighbourtable[i][j][1][0]=(i+Latticesize-1)%(Latticesize);
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

//Bereken de xy-coordinaten van alle roosterpunten
void Generate_positions()
{
    for(int i=0;i<(Latticesize);++i)
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

//Reken de reciproceroostervecotren uit
void Generate_recip_vectors()
{
    //!TODO, DEZE FACTOR MISSCHIEN PRECIEZER
    G1[1]=4.*M_PI/(sqrt(3.)*a_moir);
    double xrot=4.*M_PI/(sqrt(3.)*a_moir)*-sin(M_PI/3.);
    double yrot=4.*M_PI/(sqrt(3.)*a_moir)*cos(M_PI/3.);
    G2[0]=xrot;
    G3[0]=xrot;
    G2[1]=yrot;
    G3[1]=-yrot;

}

//Reken het magneetveld uit
void Calculate_B_field()
{
    for(int i=0;i<(Latticesize);++i)
    {
        for(int j=0;j<L;++j)
        {
            double B=0.;
            //The coordinates
            double x=Coordinates[i][j][0];
            double y=Coordinates[i][j][1];
            //EQUATION A1 FROM THE HIERARCY PAPER

            //!WEER OPLETTEN, DE FACTOR JEX ZIT NU DUS IN CALCULATEMAILTONIAN ZODAT HET MAKKELIJKER IS MET DE ANNEALING
            B=(cos(x*G1[0]+y*G1[1])+cos(x*G2[0]+y*G2[1])+cos(x*G3[0]+y*G3[1]));

            Bfield[i][j]=B;
        }
    }

}

//Load the coordinates for the nth thread
void Load_Coordinates_Parallel(char filename[128],int id)
{
    FILE* initfile=fopen(filename,"r");
    for (int i = 0; i < Latticesize; i++)
    for(int j=0;j<L;++j)
  {
      {
          /*
        double cum;
        double mals;
        int test = fscanf (initfile, "%lf\t%lf\t%lf\t%lf\t%lf\n", &(Lattice[id][i][j][0]), &(Lattice[id][i][j][1]), &(Lattice[id][i][j][2]), &(cum), &(mals));
        if (test < 5) { printf ("Problem reading particles!\n"); exit(3);}
        */
        Lattice[id][i][j][0]=0.;
        Lattice[id][i][j][1]=0.;
        Lattice[id][i][j][2]=1.;

      }
  }

  fclose (initfile);

}


//!HIER NOG ECHT FF GOED NAAR KIJKEN
double Calculate_Flip_Hamiltonian(int Latti,int Lattj,double J, double Jex, double Ac,double As, double D, int id)
{
    //Berekent de hamiltoniaan van de punten i,j met een roosterid van Li
    double Energy=0.0;
    //reeks met de buren
    int Nb[3][2]={0};

    //de spin van onze spin
    double Spin[3]={0.};
    for(int i=0;i<3;++i)
    {
        Spin[i]=Lattice[id][Latti][Lattj][i];
    }

    Nb[0][0]=Neighbourtable[Latti][Lattj][0][0]; //Right neighbour
    Nb[0][1]=Neighbourtable[Latti][Lattj][0][1];
    Nb[1][0]=Neighbourtable[Latti][Lattj][1][0]; //left neigbhour
    Nb[1][1]=Neighbourtable[Latti][Lattj][1][1];
    Nb[2][0]=Neighbourtable[Latti][Lattj][2][0];//up or down neighbour
    Nb[2][1]=Neighbourtable[Latti][Lattj][2][1];

    //Nu is het tied om de energie uit te rekenen
    //E van de Heisenbergkoppeling
    for(int i=0;i<3;++i)
    {
        Energy-=2.*J*Dot(Spin,Lattice[id][Nb[i][0]][Nb[i][1]]);
    }
    //Bereken de enkelionanisitropy en de koppeling met het externe magneetveld
    Energy-=As*Spin[2]*Spin[2];
    //!ONTHOUDEN, HIER STAAT DE FACTOR JEX BUITEN HET B-VELD, BIJ DE VORIGE CODE ERIN
    Energy-=Jex*Bfield[Latti][Lattj]*Spin[2];
    double evecupdown[3][3];

    //Kijk welke vectoren nodig zijn om de DMI en Ac uit te rekenen
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
            Cross(Spin,Lattice[id][Nb[i][0]][Nb[i][1]],spincross);
            double dotspin=0.;
            dotspin=Dot(dvec,spincross);
            //This is the DMI therm
            Energy-=2.0*D*dotspin;
            //This is the Compass therm
            Energy-=2.0*Ac*Dot(Spin,dvec)*Dot(Lattice[id][Nb[i][0]][Nb[i][1]],dvec);
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

            Cross(Spin,Lattice[id][Nb[i][0]][Nb[i][1]],spincross);
            double dotspin=0.;
            dotspin=Dot(dvec,spincross);
            //This is the DMI therm
            Energy-=2.0*D*dotspin;
            //This is the Compass therm
            Energy-=2.0*Ac*Dot(Spin,dvec)*Dot(Lattice[id][Nb[i][0]][Nb[i][1]],dvec);
        }
    }

    return Energy;

}

//!HIER NOG FF GOED NAAR KIJKEN
double Calculate_Single_Hamiltonian(int Latti,int Lattj,double J, double Jex, double Ac,double As, double D, int id)
{
    //Berekent de hamiltoniaan van de punten i,j met een roosterid van Li
    double Energy=0.0;
    //reeks met de buren
    int Nb[3][2]={0};

    //de spin van onze spin
    double Spin[3]={0.};
    for(int i=0;i<3;++i)
    {
        Spin[i]=Lattice[id][Latti][Lattj][i];
    }

    Nb[0][0]=Neighbourtable[Latti][Lattj][0][0]; //Right neighbour
    Nb[0][1]=Neighbourtable[Latti][Lattj][0][1];
    Nb[1][0]=Neighbourtable[Latti][Lattj][1][0]; //left neigbhour
    Nb[1][1]=Neighbourtable[Latti][Lattj][1][1];
    Nb[2][0]=Neighbourtable[Latti][Lattj][2][0];//up or down neighbour
    Nb[2][1]=Neighbourtable[Latti][Lattj][2][1];

    //Nu is het tied om de energie uit te rekenen
    //E van de Heisenbergkoppeling
    for(int i=0;i<3;++i)
    {
        Energy-=J*Dot(Spin,Lattice[id][Nb[i][0]][Nb[i][1]]);
    }
    //Bereken de enkelionanisitropy en de koppeling met het externe magneetveld
    Energy-=As*Spin[2]*Spin[2];
    //!ONTHOUDEN, HIER STAAT DE FACTOR JEX BUITEN HET B-VELD, BIJ DE VORIGE CODE ERIN
    Energy-=Jex*Bfield[Latti][Lattj]*Spin[2];
    double evecupdown[3][3];

    //Kijk welke vectoren nodig zijn om de DMI en Ac uit te rekenen
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
            Cross(Spin,Lattice[id][Nb[i][0]][Nb[i][1]],spincross);
            double dotspin=0.;
            dotspin=Dot(dvec,spincross);
            //This is the DMI therm
            Energy-=D*dotspin;
            //This is the Compass therm
            Energy-=Ac*Dot(Spin,dvec)*Dot(Lattice[id][Nb[i][0]][Nb[i][1]],dvec);
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

            Cross(Spin,Lattice[id][Nb[i][0]][Nb[i][1]],spincross);
            double dotspin=0.;
            dotspin=Dot(dvec,spincross);
            //This is the DMI therm
            Energy-=D*dotspin;
            //This is the Compass therm
            Energy-=Ac*Dot(Spin,dvec)*Dot(Lattice[id][Nb[i][0]][Nb[i][1]],dvec);
        }
    }

    return Energy;

}

double Calculate_E_tot(double J, double Jex, double Ac,double As, double D, int id)
{
    double Energy=0.;
    for(int i=0;i<L;++i)
    {
        for(int j=0;j<L;++j)
        {
            Energy+=Calculate_Single_Hamiltonian(i,j,J,Jex,Ac,As,D,id);
        }
    }

    return Energy;
}

void Flipspins(double kbT,double J,double Jex, double Ac,double As, double D,double alpha, int id)
{

    double acceptance=0.0;
    //Door de simplicifatie hebben we maar een hamiltoniaan nodig, de flip_hamiltonian
    //De extra energie die er nu bijkomt en afgaat
    double delta_E=0.0;

    //loop over all the partilces
    for(int i=0;i<n_part;++i)
    {
        int Ri;
        int Rj;
        Ri=(int)(L*dsfmt_genrand());
        Rj=(int)(L*dsfmt_genrand());

        double Eold=0.;
        double Enew=0.;
        //!DIT FIXEN DIT MOET ECHT SNELLER oude code werkte niet
        //Eold=Calculate_E_tot(J,Jex,Ac,As,D,id);
        Eold=Calculate_Flip_Hamiltonian(Ri,Rj,J,Jex,Ac,As,D,id);

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
        Oldspin[0]=Lattice[id][Ri][Rj][0];
        Oldspin[1]=Lattice[id][Ri][Rj][1];
        Oldspin[2]=Lattice[id][Ri][Rj][2];

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
        Lattice[id][Ri][Rj][0]=nx;
        Lattice[id][Ri][Rj][1]=ny;
        Lattice[id][Ri][Rj][2]=nz;


        //calculate the energy after the flip
        //!DIT FIXEN DIT MOET ECHT SNELLER oude code werkte niet
        //Enew=Calculate_E_tot(J,Jex,Ac,As,D,id);
        Enew=Calculate_Flip_Hamiltonian(Ri,Rj,J,Jex,Ac,As,D,id);

        double Boltzmanweight=0.;
        double randomnumb=0.;

        //check if the change is accepted or not

        Boltzmanweight=MIN(1.0,exp(-(Enew-Eold)/kbT));
        randomnumb=dsfmt_genrand();

        if(randomnumb>Boltzmanweight)
        {

            Lattice[id][Ri][Rj][0]=Oldspin[0];
            Lattice[id][Ri][Rj][1]=Oldspin[1];
            Lattice[id][Ri][Rj][2]=Oldspin[2];
        }
        else
        {
            //bereken de acceptatieratio en de extra energie
            acceptance+=1.0/n_part;
            delta_E+=(Enew-Eold);
        }

    }
    //Update de E en alfa's
    E_tot_array[id]+=delta_E;
    /*
    double alpha_old=Alpha_array[id];
    if(acceptance<0.4)
            {
                if(alpha_old>0.2)Alpha_array[id]=alpha_old-0.1;
                else if(alpha>0.05)Alpha_array[id]=alpha_old-0.04;
            }
    if((acceptance>0.6)&&(alpha_old<1.0))Alpha_array[id]=alpha_old+0.1;
    */
    //return acceptance;
    /*
    printf("ACCPETATIEMEYK\n");
    printf("%lf\t%lf\n",kbT,acceptance);
    */
}





double *Calculate_Properties_Parallel(double *results,int id)
{
    //This function calculates the properties of the lattice, which includes the magnetization, avarege topological charge and the averages of the correlations

    //calculate the average magnetization, and takes the absolute value because of invariance
    //Loop over all the particles and calculate the derivates of the spins
    double dx=sqrt(3.)/2.;
    double dy1=1.0;
    double dy2=0.5;
    double DSx[Latticesize][Latticesize][NDIM]={0.};
    double DSy[Latticesize][Latticesize][NDIM]={0.};
    for(int i=0;i<L;++i)
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
                    Sn[k][l]=Lattice[id][Nb[k][0]][Nb[k][1]][l];
                }
            }

            for(int k=0;k<3;++k)
            {
                //Calculate the x derivative
                //!STAAT DEZE WEL GOED???????, we doen hier links, min rechts nu wel goed denk
                DSx[i][j][k]=(Sn[0][k]-Sn[1][k])/(2.0*dx);
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

    double Vortarr[Latticesize][Latticesize][NDIM]={0.};
    for(int i=0;i<L;++i)
    {
        for(int j=0;j<L;++j)
        {
            //rekent de vorticiteit voor elk punt uut
            double vortvec[3]={0.};
            Cross(DSx[i][j],DSy[i][j],vortvec);
            for(int k=0;k<3;++k)
            {
                Vortarr[i][j][k]=vortvec[k];
            }
        }
    }

    //Reken de topologische lading uut
    double Topcurr[Latticesize][Latticesize]={0.};

    for(int i=0;i<L;++i)
    {
        for(int j=0;j<L;++j)
        {
            Topcurr[i][j]=Dot(Lattice[id][i][j],Vortarr[i][j]);
        }
    }

    //We gaan nu de correlaties uitrekenen, eigenlijk alleen chi_c, zie vergelijking C3 in hierchy paper
    double Chi_c=0.;
    //!TODO NORMALISATIEFACTOR BIJ DE RHO
    double rho_tot=0.;
    //Dit is de correlatie van de spin in de z-richting met zn naaste buren
    double Chi_z=0.;
    //Dit is de topologische lading in het kwadraat
    double rho2=0.;
    for(int i=0;i<L;++i)
    {
        for(int j=0;j<L;++j)
        {
            double rho=0.;
            double sz=0.;
            //de spin in de z-richting
            sz=Lattice[id][i][j][2];
            rho=Topcurr[i][j];
            //bereken de totale topologische lading
            rho_tot+=rho;
            //Nu aangepast dat voor Chic alleen de nearest neighbours worden gepakt
            int Nb[3][2]={0};
            //Get the neighbours
            Nb[0][0]=Neighbourtable[i][j][0][0]; //Right neighbour
            Nb[0][1]=Neighbourtable[i][j][0][1];
            Nb[1][0]=Neighbourtable[i][j][1][0]; //left neigbhour
            Nb[1][1]=Neighbourtable[i][j][1][1];
            Nb[2][0]=Neighbourtable[i][j][2][0];//up or down neighbour
            Nb[2][1]=Neighbourtable[i][j][2][1];
            for(int k=0;k<3;++k)
            {
                int i_i=Nb[k][0];
                int j_i=Nb[k][1];
                Chi_c+=rho*(Topcurr[i_i][j_i])/n_part;
                Chi_z+=sz*Lattice[id][i_i][j_i][2];
            }
            //bereken rho^2
            rho2+=rho*rho;
        }
    }

    results[0]=rho_tot;
    results[1]=Chi_c;
    results[2]=Chi_z;
    results[3]=rho2;


    return results;
}

double Calculate_magnetism_Parallel(int id)
{
    //calculates the total magnetism of the lattice
    double magnetism=0.;
    for(int i=0;i<Latticesize;++i)
    {
        for(int j=0;j<Latticesize;++j)
        {
            magnetism+=Lattice[id][i][j][2]/(double)n_part;
        }
    }
    return magnetism;
}

void Lattice_Saver_Parallel(char name[128],int id)
{
    FILE* fp=fopen(name,"w");
    for(int i=0;i<(L);++i)
    {
        for(int j=0;j<L;++j)
        {
            double Sx=Lattice[id][i][j][0];
            double Sy=Lattice[id][i][j][1];
            double Sz=Lattice[id][i][j][2];
            fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\n",Sx,Sy,Sz);
        }
    }
    fclose(fp);
}


int main(int argc, char *argv[])
{
    //De standaardwaarden voor D,Jex en As
    double D=0.18;
    double Jex=0.05;//0.1
    double Asfactor=0.;
    double Acfactor=0.;
    double As=0.0;
    double Ac=0.0;
    //ding voor saus
    char *eptr;
    if(argc==1)
    {
        printf("Geen argumenten, gebruikt standaardwaarden\n");
        printf("D=%lf\tJex=%lf\n",D,Jex);
    }
    else if(argc==2)
    {
        printf("Alleen een waarde voor D ingevoerd\n");
        D=strtod(argv[1],&eptr);
        printf("D=%lf\tJex=%lf\n",D,Jex);
    }
    else if(argc==3)
    {
        printf("Zowel voor D als J een waarde ingevoerd\n");
        D=strtod(argv[1],&eptr);
        Jex=strtod(argv[2],&eptr);
        printf("D=%lf\tJex=%lf\n",D,Jex);
    }
    else if(argc==4)
    {
        printf("Voor D, Jex en de Ac-factor een waarde ingevoerd\n");
        D=strtod(argv[1],&eptr);
        Jex=strtod(argv[2],&eptr);
        Acfactor=strtod(argv[3],&eptr);
        Ac=Acfactor*D*D;
        printf("D=%lf\tJex=%lf\tAs/D^2=%lf\n",D,Jex,Acfactor);
        printf("As=%lf\n",Ac);
    }
    else if(argc==5)
    {
        printf("Voor D, Jex, Acfactor en As facot een waarde ingevoerd\n");
        D=strtod(argv[1],&eptr);
        Jex=strtod(argv[2],&eptr);
        Acfactor=strtod(argv[3],&eptr);
        Ac=Acfactor*D*D;
        Asfactor=strtod(argv[4],&eptr);
        As=Asfactor*D*D;
        printf("D=%lf\tJex=%lf\tAc/D^2=%lf\n",D,Jex,Acfactor);
        printf("Ac=%lf\n",Ac);
        printf("D=%lf\tJex=%lf\tAs/D^2=%lf\n",D,Jex,Asfactor);
        printf("As=%lf\n",As);

    }
    else
    {
        printf("HIER GAAT IETS FOUT\n");
    }


    //De eq step is de hoeveelheid stappen voordat er een temperatuurwissel wordt gepbropeerd
    int eq_step=1;
    //Data_step zijn de hoeveelheid temperatuurwissels voordat er een datapunt wordt opgeslagen
    int data_step=100; ///KAN MISSCHIEN NOG WEL EEN STUK HOGER
    int data_points=(MC_steps-init_steps)/data_step;
    double ddp=(double)data_points;

    double Mag_res[NTHREADS]={0.};
    double E_res[NTHREADS]={0.};
    double E2_res[NTHREADS]={0.};
    double Cv_res[NTHREADS]={0.};
    double Rho_res[NTHREADS]={0.};
    double Chic_res[NTHREADS]={0.};
    double Chiz_res[NTHREADS]={0.};
    double Rho2_res[NTHREADS]={0.};

    dsfmt_seed(time(NULL));

    //Het aantal kernen dat ik wil gebruiken
    //omp_set_num_threads(NTHREADS);
    //De temperatuurfactoren, dit geeft voor 12 kernen een tmax van 2.98
    double T_min=0.1; //Dit zodat er 10 extra bufferplekken tussen dit punt en T=0.1 zitten
    //zorgt voor een maximum T van ongeveer 1.0 enzo, lekker smullen
    double T_factor=1.05;

    //Calculate the positons and stuff
    Generate_positions();
    Generate_recip_vectors();
    Calculate_B_field();
    Calculate_neighbourtable();
    //initalize the alttice indexes
    Generate_Lattice_index();
    //Calculate the temperatures
    //Hier initialiseren we eerst de temperaturen  op een geometrische manier
    Generate_Temperatures_Geometric(T_min,T_factor);
    double T_max=Temperature_array[NTHREADS-1];

    //!NIET ZO SNEL ALS IE ZOU KUNNEN, MAAR GAAT ZO TE ZIEN IGG WEL GOED maar dus nog wel false sharing
    //#pragma omp parallel
    for(int i=0;i<NTHREADS;++i)
    {
        //int id=omp_get_thread_num();
        Load_Coordinates_Parallel(nameskyr,i);
    }

    //Laten equilibreren
    //for(int i=0;i<MC_eq;++i)
    //{
        //!PARALLEL
        //Hier is ook een hele hoop false sharing, maar dit zou in principe moeten werken, maar doet het dus niet we verliezen altijd een of meerdere kernen
    for(int id=0;id<NTHREADS;++id)
    {
        //double kbT=Temperature_array[id];
        //Flipspins(kbT,J,Jex,Ac,As,D,alpha,id);
        E_tot_array[id]=Calculate_E_tot(J,Jex,Ac,As,D,id);
    }

    //}
    /*
    printf("ENERGIËN NAAT EQUILIBRUIEREN\n");
    for(int i=0;i<NTHREADS;++i)
    {
        printf("T=%lf\t E=%lf\n",Temperature_array[i],E_tot_array[i]);
    }
    */

    printf("De temperaturen aan het begin:\n");
    for(int i=0;i<NTHREADS;++i)
    {
        printf("%.3f\t",Temperature_array[i]);
    }
    printf("\n");


    //initialiseer de array die bijhoudt of een temperatuur een n_up of een n_down is
    for(int i=0;i<NTHREADS;++i)
    {
        //!TODO KAN MISSCHIEN EFFICIENTER ZODAT ER NIET CONTINU OVER TWEE ARRAY'S GELOOPT HOEFT TE WORDEN
        Nup_arr[i]=false;
        Ndown_arr[i]=false;
    }

    //!We hebben nu startgrids en de energieën, we gaan nu de temperaturen optimaliseren
    int N_sw=100000; //De hoeveelheid sweeps per temepratuur equilbriatiestap
    for(int i=0;i<5;++i)
    {
        //doe 5 stappen om de temperaturen de configureren
        //we voeren nu ook alle sweepstappen uit om temperaturen te configueren
        for(int j=0;j<N_sw;++j)
        {
            //Voor elke sweep registeren of iets een n_up of n_down moet krijgen en dan deze dingen ook bij die histogrammen erbij optellen

            //eerst die MC stappen
            for(int id=0;id<NTHREADS;++id)
            {
                //!TODO kijken of het zo wel helemaal goed gaat met die alpha  NEE, DAT GAAT HET NIET
                double kbT=Temperature_array[id];
                double alpha=0.85*sqrt(kbT);
                Flipspins(kbT,J,Jex,Ac,As,D,alpha,id);
            }
            //nu zijn alle spins geflipt, nu de array's updaten
            //verander de waarden van de array's bij de hoogste en laagste temperatuur naar true

            //nu zijn ze geupdate en kunnen we pogen de temperaturen om te wisselen

            /// ff voor zorgen dat dit alleen met buren gaat
            ///Nu lekker de temperatuurtjes gaan proberen te wisselen.
            if(j%2==0)
            {
                for(int k=0;k<(NTHREADS-1);k=k+2)
                {
                    //Dit izjn de indexen van de gewenste temperaturen
                    int L_i=Lattice_index[k];
                    int L_j=Lattice_index[k+1];
                    //En haal daaruit de energieën enzo
                    double kbT_i=Temperature_array[L_i];
                    double kbT_j=Temperature_array[L_j];
                    double E_i=E_tot_array[L_i];
                    double E_j=E_tot_array[L_j];


                    double Boltzmanweight=0.;
                    double randomnumb=0.;

                    Boltzmanweight=MIN(1.0,exp((1./kbT_j-1./kbT_i)*(E_j-E_i)));
                    randomnumb=dsfmt_genrand();

                    if(randomnumb<Boltzmanweight)
                    {
                        //De nieuwe temperaturen
                        Temperature_array[L_i]=kbT_j;
                        Temperature_array[L_j]=kbT_i;
                        //en de nieuwe indexen van de temperaturen
                        Lattice_index[k]=L_j;
                        Lattice_index[k+1]=L_i;
                        //acceptance+=1/127.;
                    }
                }
            }
            else
            {
                ///probeer hem hier voor de oneven temperaturen eerst
                for(int k=1;k<(NTHREADS-1);k=k+2)
                {
                    //Dit izjn de indexen van de gewenste temperaturen
                    int L_i=Lattice_index[k];
                    int L_j=Lattice_index[k+1];
                    //En haal daaruit de energieën enzo
                    double kbT_i=Temperature_array[L_i];
                    double kbT_j=Temperature_array[L_j];
                    double E_i=E_tot_array[L_i];
                    double E_j=E_tot_array[L_j];


                    double Boltzmanweight=0.;
                    double randomnumb=0.;

                    Boltzmanweight=MIN(1.0,exp((1./kbT_j-1./kbT_i)*(E_j-E_i)));
                    randomnumb=dsfmt_genrand();

                    if(randomnumb<Boltzmanweight)
                    {
                        //De nieuwe temperaturen
                        Temperature_array[L_i]=kbT_j;
                        Temperature_array[L_j]=kbT_i;
                        //en de nieuwe indexen van de temperaturen
                        Lattice_index[k]=L_j;
                        Lattice_index[k+1]=L_i;
                        //acceptance+=1/127.;
                    }
                }
            }

            ///DEZE MISSCHIEN WEER TERUG EENTJE NAAR BOVEN
            int L_laagste=Lattice_index[0];
            int L_hoogste=Lattice_index[NTHREADS-1];
            Nup_arr[L_laagste]=true;
            Ndown_arr[L_laagste]=false;
            Nup_arr[L_hoogste]=false;
            Ndown_arr[L_hoogste]=true;

            //temperaturen zijn geupdate, nu de histogrammen updaten
            //doe bij nup+1 als er een array met up inzit, en vice versa, als een ding nog geen markering heeft, gewoon skippen jwz
            for(int id=0;id<NTHREADS;++id)
            {
                //kijken welke temperatuur t ding heeft
                int L_i=Lattice_index[id];
                bool up=Nup_arr[L_i];
                if(up==true)
                {
                    nup[id]+=1;
                }
                else
                {
                    bool down=Ndown_arr[L_i];
                    if (down==true)
                    {
                        ndown[id]+=1;
                    }
                }
            }
            ///klaar met de for loop voor de sweeps nu, de waardens voor ndown enzo zijn geupdate
        }
        //Die meuk is nu geupdate, dus we gaan nu ft en etat uitrekenen
        ///DIT KAN MET EEN PRECIEZERE MANIER, WSS LATER FF IMPLEMENTEREN
        printf("De waardes voor enzo zijn:\n");

        for(int id=0;id<NTHREADS;++id)
        {
            //pak weer de index voor de juiste temperatuur
            double dnup=(double)nup[id];
            double dndown=(double)ndown[id];
            double sum=dnup+dndown;
            if(sum!=0.)
            {
                ft[id]=dnup/sum;
            }
            else
            {
                ft[id]=0.5;
            }
            printf("dnup: %lf, dndown: %lf, ft: %lf\n",dnup,dndown,ft[id]);
        }

        //en nu eta uitrekenen
        //dit zijn de waarden van eta en de waarde van de constante
        double inteta=0.0;
        double C=0.;


        ///HIER WEL WEER EEN GEDEELTE MET DIE F MODIFICEREN
        printf("DE NIEUW F DINGEN:");
        for(int id=0;id<NTHREADS;++id)
        {
            double ftid=ft[id];
            ft[id]=ftid;
            printf("%lf\n",ft[id]);
        }



        for(int id=0;id<NTHREADS-1;++id)
        {
            //pak weer de index voor de juiste temperatuur
            int L_i=Lattice_index[id];
            int L_j=Lattice_index[id+1];
            double DT;
            DT=Temperature_array[L_j]-Temperature_array[L_i];
            double df;
            df=ft[id+1]-ft[id];
            double etaa;
            ///HIER GAAT HET FOUT misschien ook niet die DT in het kwadraat
            if(df<0.)etaa=sqrt(-df/(DT*DT));
            else etaa=sqrt(df/(DT*DT));/// HIER DIT ZO GEDAAN DAT NIET T HELE GEBIED MET 0 WORDT GESKIPT-sqrt(df/(DT*DT));
            etat[id]=etaa;
            DT_arr[id]=DT;
            //reken dat integraaltje uut
            inteta+=etaa*DT;
            printf("id::\t DT=%lf\t df=%lf\t etaa=%lf\t inteta=%lf\n",DT,df,etaa,inteta);
        }
        printf("De waarde van inteta is::: %lf\n",inteta);
        C=1.0/inteta;
        printf("DE WAARDE VAN C IS::: %lf\n",C);
        //we hebben nu waardes en het is nu dus tijd om de optimale waardes van de temperatuur uit te rekenen
        //loop over alle waarden heen die tussen de twee extrama's liggen
        //!DIT INTEGREREN KAN SWS WEL BETERRRRRRR WANT DIT IS WEL ECHT MATIGE SHIT AH SAHBI
        for(int id=1;id<(NTHREADS-1);++id)
        {
            //Dit is de waarde die het integraal moet hebben
            double val=((double)id)/((double)NTHREADS-1.);
            //int L_i=Lattice_index[id]; //ja gewoon weer de lattice index zoals altijd
            double etaint=0.;//de waarde van het integraal tot dan toe
            int id_t=0; //dit is de index van de temperatuur totaan waar we tot nu toe hebben geintergreerd
            int id_t1=1;
            ///allebei op id_t geïnitialiseerd want valt staks toch weer weg
            ///NOG FF NETJES DOEN MET DIE ID_T ENZO
            int L_t=Lattice_index[id_t];
            int L_t1=Lattice_index[id_t1];
            while(etaint<val)
            {
                //verhoog de waarde van het integraal totdat t te hoog is
                etaint+=C*etat[id_t]*DT_arr[id_t];
                //update the waardes voor id_t enzo
                id_t++;
                id_t1++;
                L_t=L_t1;
                L_t1=Lattice_index[id_t1];
            }
            //nu hebben we het integraal, maar de waarde is dus nu te hoog, dus rekenen we t verschil uit en kijken hoeveel lager de temperatuur moet
            double diff=etaint-val;
            double minT;
            ///Dit is hoeveel de temperatuur lager moet
            ///NOG KIJKEN NAAR DIE -1 DINGEN HIERO
            ///FACTOR C DRBIE HIER
            minT=-diff/(C*etat[id_t-1]);
            //Dit is de gemodificeerde temperatuur voor deze index
            ///ALS HET GOED IS DIT NU WEL GOED
            double Tmod=Temperature_array[L_t]+minT;
            Temperature_array_mod[id]=Tmod;
        }
        //Nu zijn alle gemodificeerde temperaturen berekent, dus we gaan de T_array updaten
        for(int id=1;id<NTHREADS-1;++id)
        {
            ///ALS HET GOED IS, ZORGT DIT ERVOOR DAT ER GOED MET DE VERSCHILLENDE GRIDS WEER WORDT DOORGEREKEND EN DAT t_MIN EN MAX HETZELFDE BLIJVEN
            int L_i=Lattice_index[id];
            Temperature_array[L_i]=Temperature_array_mod[id];
        }


        //De temperaturen zijn geupdate, Nsweeps updaten en de histogrammen resetten en dan kunnen we nog een keer
        ///N_sw*=2;
        for(int j=0;j<NTHREADS;++j)
        {
            //!TODO KAN MISSCHIEN EFFICIENTER ZODAT ER NIET CONTINU OVER TWEE ARRAY'S GELOOPT HOEFT TE WORDEN
            Nup_arr[j]=false;
            Ndown_arr[j]=false;
            nup[j]=0;
            ndown[j]=0;
        }
        char naam[128];
        sprintf(naam,"ft_run%d_D=%.3f_Jex=%.3f.txt",i,D,Jex);
        FILE* fp=fopen(naam,"w");
        printf("De temperaturen zijn aangepast naar:\n");
        for(int id=0;id<NTHREADS;++id)
        {
            int L_i=Lattice_index[id];
            printf("%.3f\t",Temperature_array[L_i]);
            fprintf(fp,"%lf\t%lf\n",Temperature_array[L_i],ft[id]);
        }
        fclose(fp);
        printf("\n");


        for(int i=0;i<NTHREADS;++i)
        {
            //int id=omp_get_thread_num();
            Load_Coordinates_Parallel(nameskyr,i);
        }

        //Laten equilibreren
        //for(int i=0;i<MC_eq;++i)
        //{
            //!PARALLEL
            //Hier is ook een hele hoop false sharing, maar dit zou in principe moeten werken, maar doet het dus niet we verliezen altijd een of meerdere kernen
        for(int id=0;id<NTHREADS;++id)
        {
            //double kbT=Temperature_array[id];
            //Flipspins(kbT,J,Jex,Ac,As,D,alpha,id);
            E_tot_array[id]=Calculate_E_tot(J,Jex,Ac,As,D,id);
        }

    }


    ///TEMPERATUREN ZIJN HOPELIJK GOED NU, TIJD OM HELEMAAL LOESOE TE GAAN
    ///Nu de daadwerkelijke MC uitvoeren met PT
    /*
    for(int i=0;i<MC_steps;++i)
    {
        //!PARALLEL
        for(int id=0;id<NTHREADS;++id)
        {
            ///ALFA MOET NOG GEUPDATE WORDEN
            //double alpha=0.3;//Alpha_array[id];
            double kbT=Temperature_array[id];
            //Doe ff vijf sweeps voordat we de T gaan wisselen
            for(int j=0;j<eq_step;++j)
            {
                //!HIER GAAT IETS GIGANTISCH MIS
                Flipspins(kbT,J,Jex,Ac,As,D,alpha,id);
                //printf("%lf\n",E_tot_array[j]);
            }
            //Na het doen van de MC-steppen de energieën uitrekenen
            //!TODO SNELLER, doen met +- enzo
            //E_tot_array[id]=Calculate_E_tot(J,Jex,Ac,As,D,id);
        }
        //De magnetisisatie uitrekenen, dit per temperatuur, maak op de datapuntstappen een datapunt

        double acceptance=0.0;
        //Nu lekker de temperatuurtjes gaan proberen te wisselen. Als i even is pakken we gewoon alle temperaturen, anders de laagste en hoogste nie
        for(int j=0;j<(NTHREADS-1);++j)
        {

            //Dit izjn de indexen van de gewenste temperaturen
            int L_i=Lattice_index[j];
            int L_j=Lattice_index[j+1];
            //En haal daaruit de energieën enzo
            double kbT_i=Temperature_array[L_i];
            double kbT_j=Temperature_array[L_j];
            double E_i=E_tot_array[L_i];
            double E_j=E_tot_array[L_j];
            //printf("%d\t%d",L_i,L_j);
            //printf("%lf\t%lf\n",E_i,E_j);

            double Boltzmanweight=0.;
            double randomnumb=0.;
            /*
            printf("DE KBTS\n");
            printf("De lage:%lf\n",kbT_i);
            printf("DE HOGE %lf\n",kbT_j);
            /

            //check if the change is accepted or not

            Boltzmanweight=MIN(1.0,exp((1./kbT_j-1./kbT_i)*(E_j-E_i)));
            randomnumb=dsfmt_genrand();
            //Als ie wordt geaccepteerd de dingen updaten, anders niks doen
            //printf("B=%lf,random=%lf\n",Boltzmanweight,randomnumb);
            if(randomnumb<Boltzmanweight)
            {
                //De nieuwe temperaturen
                Temperature_array[L_i]=kbT_j;
                Temperature_array[L_j]=kbT_i;
                //en de nieuwe indexen van de temperaturen
                Lattice_index[j]=L_j;
                Lattice_index[j+1]=L_i;
                //printf("HUUUTTTS\n");
                acceptance+=1/((double)NTHREADS);
            }


        }

        //printf("De acceptatie:::%lf\n",acceptance);


        if(((i%data_step)==0)&&(i>init_steps))
        {
            for(int j=0;j<NTHREADS;++j)
            {
                //We kiekn eerst ff welke temperatuur welke index heeft
                int L_i=Lattice_index[j];
                //Reken alle eigenschappen enzo uit, alleen Cv nog niet, die komt aant einde
                double m=fabs(Calculate_magnetism_Parallel(L_i)); //! ABSOLUTE WAARDE
                Mag_res[j]+=m/ddp;
                double E=E_tot_array[L_i];
                E_res[j]+=E/ddp;
                ///KIJKEN WAT ER MISGAAT BIJ DE E2
                E2_res[j]+=(E*E)/(ddp);
                //reken al die dingen zoals topoligische lading enzo uit
                //!CHIZ ENZO ZIJN WEGGEHAALD
                double results[4]={0.};
                Calculate_Properties_Parallel(results,L_i);
                Rho_res[j]+=(results[0]/ddp); //!OOk hier de absolute waarde //maar weghalen lijkt me denk beter?????
                Chic_res[j]+=results[1]/ddp;
                Chiz_res[j]+=results[2]/ddp;
                Rho2_res[j]+=results[3]/ddp;
            }
        }
        /*
        for(int j=0;j<NTHREADS;++j)
        {
            printf("%lf\n",E_tot_array[j]);
        }
        /


    }
    */

    return 0;
}

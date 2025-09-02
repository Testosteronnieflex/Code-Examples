#include <stdio.h>
#include <time.h>
//#include <omp.h>
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <complex.h>
//#include <windows.h>
#include "mt19937.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif // M_PI

#define MIN(X,Y) (((X)<(Y)) ? (X):(Y))

//The size of the lattice we will concider
#define L 1024//4096/
#define N (L*L)

#define N_core 32

#define N_core2 (N_core/2)

//The number of MC steps in our simulation (probably not needed)
#define MC_steps 1000

//The maximum number of allowed iterations
#define MAX_I 10000

//The number of sweeps for the autocorrelation funciton
#define correlation_sweeps 20000


//The coupling constant
float J=1.;
//The uniformly applied magnetic field_sum
float H=-0.02;
//The variance in the values of h_i
float std_h_i=1.0;
//The initial water fraction
float F_in=0.5;

double T=0.5;

double P_add;

//The array with the external field
float H_field[L][L]={0.};

//The array with all the spins, +1 represents water, -1 represents ice
//!TODO 1D MAKEN als dat misschien wat uitmaakt
short S[N_core][L][L]={0};

unsigned short core_arr[N_core];

bool active_array[N_core2][L][L]={false};

bool active_cluster[N_core2][L][L]={false};

//This contains the average magnetization for all the samples
double Correlation_array[correlation_sweeps];


double Correlation_time[correlation_sweeps];


//!DEFINING ALL THE FUNCTIONS
double gaussrand(double mean, double std);
double r2();
double random_power(double xmin, double xmax,double degree);
int biased_binary();
short random_binary();
void Init_H();
void Init_S_parallel();
double calculate_E_parallel(int core,int i,int j);
double calculate_E_v2_parallel(int core,int i,int j);
void sweep_T_parallel(int core);
void finite_T_sim_parallel(int core,int it);
void Init_core_index();
void Shuffle_core_index();
void Active_sites(int core_1,int core_2,int order);
void Reset_cluster(int order);
void Calculate_P_add();
void Cluster_sweep(int core_1, int core_2, int order);
void Perform_Cluster(int steps);
int Find(int x, int labels[]);
void Union(int x, int y, int labels[]);
int maxValue(int myArray[], size_t size);
//void Cluster_sizes_parallel(char filename[],char filename2[],int core);
//void print_S_parallel(char filename[],int core);
double Mean_parallel(int core);

double gaussrand(double mean, double std)
{
    //Generatates a random number that is gaussian distributed
	static double U, V;
	static int phase = 0;
	double Z;

	if(phase == 0) {
		U = (rand() + 1.) / (RAND_MAX + 2.);
		V = rand() / (RAND_MAX + 1.);
		Z = sqrt(-2 * log(U)) * sin(2 * M_PI * V);
	} else
		Z = sqrt(-2 * log(U)) * cos(2 * M_PI * V);

	phase = 1 - phase;

	return mean+std*Z;
}

double r2()
{
    //Generates a random nubmer between 0 and 1
    return (double)rand() / (double)RAND_MAX ;
}


double random_power(double xmin,double xmax,double degree)
{
    //Generates a random number with power law distribution

    //!DIT WORDT ONMEUIG VAAK HERHAALD, KAN VEUL SNELLUR
    double dp1=degree+1.;
    double pmin = pow(xmin, dp1);
    double pmax = pow(xmax, dp1);
    double v=pmin+(pmax-pmin)*r2();
    return pow(v,1.0/(dp1));
}

int biased_binary()
{
    //Reteruns either -1 or 1 depending with a bias determined by F_in
    double ran=r2();
    int res;

    if(ran<F_in)
    {
        res=1;
    }
    else
    {
        res=-1;
    }
    return res;
}

short random_binary()
{
    double ran=r2();
    short res;

    if(ran<0.5)
    {
        res=1;
    }
    else{
        res=-1;
    }
    return res;
}

void Init_H()
{
    //Initializes the H-field
    for (int i=0;i<L;i++)
    {
        for (int j=0;j<L;j++)
        {
            //Here we no get a random number with the specified variance
            H_field[i][j]=gaussrand(H,std_h_i)+H;
        }
    }
}

void Init_S_parallel()
{
    //generates the initial grids for all the systems
    for(int k=0;k<N_core;k++)
    {
        for (int i=0;i<L;i++)
        {
            for (int j=0;j<L;j++)
            {
                //This creates a random number that is either -1 or 1
                S[k][i][j]=biased_binary();
            }
        }
    }
}

double calculate_E_parallel(int core,int i,int j)
{
    //This function calculates the energy of a gridpoint i,j
    //we use periodic BC's
    //First calculate the energy from the interaction with the neighbours
    double sum=0.;
    int spin=S[core][i][j];
    //Callculate the interaction with the four neighbours
    sum-=J*(double)(spin*S[core][((i+1+L)%L)][j]);
    sum-=J*(double)(spin*S[core][((i-1+L)%L)][j]);
    sum-=J*(double)(spin*S[core][i][((j+1+L)%L)]);
    sum-=J*(double)(spin*S[core][i][((j-1+L)%L)]);
    //Add the external field
    sum-=H_field[i][j]*(double)spin;
    return sum;
}

double calculate_E_v2_parallel(int core,int i,int j)
{
    //Uses the weird update rules from the paper
    int sum=0;
    int spin=S[core][i][j];
    //Check if there is a majority in the allignment
    sum-=(spin*S[core][((i+1+L)%L)][j]);
    sum-=(spin*S[core][((i-1+L)%L)][j]);
    sum-=(spin*S[core][i][((j+1+L)%L)]);
    sum-=(spin*S[core][i][((j-1+L)%L)]);
    //if it tied update according to the external field
    if (sum==0)
    {
        double field=H_field[i][j]*(double)spin;
        if (field<0.0)
        {
            sum=-1;
        }
        else{
            sum=1;
        }
    }
    return (double)sum;
}

/*
void sweep_T_parallel(double T,int core)
{
    //!DIT KAN SNELLER DOOR HET EXTERNE VELD EN DE INTERACTIES APART TE ZIEN
    //Performs sweeps at finite temperature
    for(int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++)
        {
            double E_old=calculate_E_parallel(core,i,j);
            //Now perform a swap to see what the energy will become
            S[core][i][j]=-S[core][i][j];
            double E_new=calculate_E_parallel(core,i,j);
            //Now accept the move according to the Metropolis rules
            double Boltzmanweight=MIN(1.0,exp(-(E_new-E_old)/T));
            double randomnumb=dsfmt_genrand();
            if(randomnumb>Boltzmanweight)
            {
                //undo the move
                S[core][i][j]=-S[core][i][j];
            }
        }
    }
}
*/

void finite_T_sim_parallel(int core,int it)
{
    //Performs a simulation for a specified temperature and specified number of iterations
    for(int i=0;i<it;i++)
    {
        sweep_T_parallel(core);
    }
}

void Init_core_index()
{
    //initializes an array with the core array indices
    for(int i=0;i<N_core;i++)
    {
        core_arr[i]=i;
    }
}


void Shuffle_core_index()
{
    //This creates a random permutation of the core_arr
    for(int i=N_core-1;i>0;i--)
    {
        //here we switch arround the values
        int j=rand()%(i+1);
        unsigned short temp=core_arr[i];
        core_arr[i]=core_arr[j];
        core_arr[j]=temp;
    }
}

void Active_sites(int core_1,int core_2,int order)
{
    //This generates an array with the active sites, is true for points that are active
    for(int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++)
        {
            if(S[core_1][i][j]==S[core_2][i][j])
            {
                active_array[order][i][j]=false;
            }
            else
            {
                active_array[order][i][j]=true;
            }
        }
    }

}


void Reset_cluster(int order)
{
    for(int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++)
        {
            active_cluster[order][i][j]=false;
        }
    }
}


void Calculate_P_add()
{
    P_add=1-exp(-4.0/T);
}


void Cluster_sweep(int core_1, int core_2, int order)
{
    //Before we can do a cluster sweep, we first need to reset the active clustersb
    //!WSS NIE NODIG WANT WE LOOPEN ALLUS WEG
    //Reset_cluster();
    //We first calculate all the active sites



    Active_sites(core_1,core_2,order);
    //We now have all the active sites and will thus build a cluster
    //For this we will loop over all elements and try to build a cluster from there
    for(int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++)
        {
            if(active_array[order][i][j])
            {
                unsigned short cluster_index_array[N][2];
                int cluster_elements=1;
                int tries=0;
                //We have found an active site so we will now try to build a cluster
                cluster_index_array[0][0]=i;
                cluster_index_array[0][1]=j;
                //Because it is now already part of a cluster, we can delete it from the active sites
                active_array[order][i][j]=false;

                short s_cluster=S[core_1][i][j];

                //we will try to add elements until we have tried to add as many elements in the cluster as there have been tries
                while(tries<cluster_elements)
                {
                    int i_try=cluster_index_array[cluster_elements-1][0];
                    int j_try=cluster_index_array[cluster_elements-1][1];

                    //Now we check for each neighbour if it is an active site, and if it is, we try to add it, if added we delete it from the active sites


                    if(active_array[order][((i+1+L)%L)][j])
                    {
                        //This is an active site, so we will now look if they actually have the same sign
                        if(s_cluster=S[core_1][((i+1+L)%L)][j])
                        {
                            //generate a random number to see if it can be added, it it can be added modify all this shit
                            if(r2()<P_add)
                            {
                                cluster_index_array[cluster_elements][0]=((i+1+L)%L);
                                cluster_index_array[cluster_elements][1]=j;
                                cluster_elements++;
                                active_array[order][((i+1+L)%L)][j]=false;

                            }

                        }

                    }
                    if(active_array[order][((i-1+L)%L)][j])
                    {
                        if(s_cluster=S[core_1][((i-1+L)%L)][j])
                        {
                            if(r2()<P_add)
                            {
                                cluster_index_array[cluster_elements][0]=((i-1+L)%L);
                                cluster_index_array[cluster_elements][1]=j;
                                cluster_elements++;
                                active_array[order][((i-1+L)%L)][j]=false;
                            }

                        }

                    }
                    if(active_array[order][i][((j+1+L)%L)])
                    {
                        if(s_cluster=S[core_1][i][((j+1+L)%L)])
                        {
                            if(r2()<P_add)
                            {
                                cluster_index_array[cluster_elements][0]=i;
                                cluster_index_array[cluster_elements][1]=((j+1+L)%L);
                                cluster_elements++;
                                active_array[order][i][((j+1+L)%L)]=false;
                            }
                        }

                    }
                    if(active_array[order][i][((j-1+L)%L)])
                    {
                        if(s_cluster=S[core_1][i][((j-1+L)%L)])
                        {
                            if(r2()<P_add)
                            {
                                cluster_index_array[cluster_elements][0]=i;
                                cluster_index_array[cluster_elements][1]=((j-1+L)%L);
                                cluster_elements++;
                                active_array[order][i][((j-1+L)%L)]=false;
                            }

                        }

                    }
                    //We have now done one more try, so it needs to be increased by one
                    tries++;

                }

                //We now have one cluster, we will now give the elements in the clusters a random spin in both clusters as described in the paper


                //These here are the new new spins for the two clusters
                short s_new_1=random_binary();
                short s_new_2=-s_new_1;

                for(int k=0;k<cluster_elements;k++)
                {
                    //We will now update all the spins in both the clusters
                    //first we get the indices of both and after that we just update them for both
                    int i_update=cluster_index_array[k][0];
                    int j_update=cluster_index_array[k][1];

                    //update step
                    S[core_1][i_update][j_update]=s_new_1;
                    S[core_2][i_update][j_update]=s_new_2;

                }


            }
        }
    }



    //!AT THE MOMENT RUN IN SERIAL, MAY NOT BE THE MOST EFFICIENT GODVERDOMME

    sweep_T_parallel(core_1);
    sweep_T_parallel(core_2);





}

void sweep_T_parallel(int core)
{
    //!DIT KAN SNELLER DOOR HET EXTERNE VELD EN DE INTERACTIES APART TE ZIEN
    //Performs sweeps at finite temperature
    for(int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++)
        {
            double E_old=calculate_E_parallel(core,i,j);
            //Now perform a swap to see what the energy will become
            S[core][i][j]=-S[core][i][j];
            double E_new=calculate_E_parallel(core,i,j);
            //Now accept the move according to the Metropolis rules
            double Boltzmanweight=MIN(1.0,exp(-(E_new-E_old)/T));
            double randomnumb=dsfmt_genrand();
            if(randomnumb>Boltzmanweight)
            {
                //undo the move
                S[core][i][j]=-S[core][i][j];
            }
        }
    }
}


void Perform_Cluster(int steps)
{

    //We start by calculating the P_add
    Calculate_P_add();
    //reset the order of the core-indices
    Init_core_index();


    //loop over the number of steps
    //for each step we will do the following:
    //create a random permutation of which grids are coupled with which
    //divide the different clusters over the different cores
    //perform a cluster sweep over them (both the clustermove and regular sweeps are done at the same time
    //


    for(int i=0;i<steps;i++)
    {
        //perform the random swapping of the indices
        Shuffle_core_index();
        //we now have the shuffeld core and will now just do the main simulation
            for(int i=0;i<N_core;i+=2)
            {
                //printf("DEZE RUND OP DEZE CORE:::  %d\n",omp_get_thread_num());
                int core_1=core_arr[i];
                int core_2=core_arr[i+1];

                Cluster_sweep(core_1,core_2,i/2);
            }

    }
}

int Find(int x, int labels[])
{
    int y=x;
    while(labels[y]!=y)
    {
        y=labels[y];
    }
    while (labels[x]!=x)
    {
        int z = labels[x];
        labels[x] = y;
        x = z;
    }
    return y;
}

void Union(int x, int y,int labels[])
{
    labels[Find(x,labels)] = Find(y,labels);
}


int maxValue(int myArray[], size_t size)
{
    /* enforce the contract */
    assert(myArray && size);
    size_t i;
    int maxValue = myArray[0];

    for (i = 1; i < size; ++i) {
        if ( myArray[i] > maxValue ) {
            maxValue = myArray[i];
        }
    }
    return maxValue;
}


void Cluster_Sizes_Parallel(int core, char filename[128],char filename2[128])
{
    printf("test\n");
    //Generates the cluster characheristics of all the samples
    //label_arr contains the label for each point
    int largest_label=0;

    int label_arr[L][L]={0};

    int labels[N];
    //contains the circumferences
    int circumferences[L][L]={0};

    //we are now going to calculate the total number of clusters while also do all the other stuff at the same time
    int total_clusters=0;
    //initialize the labels and also the circumferences
    for(int i=0;i<N;i++)
    {
        labels[i]=i;
    }



    //Now loop over the entire grid to label it
    for(int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++)
        {
            //check if there is actually a melt pond at the site
            if(S[core][i][j]==1.)
            {
                int left=label_arr[((i-1+L)%L)][j];
                int up=label_arr[i][((j-1+L)%L)];
                if (left==0 && up==0)
                {
                    //neither a label above or to the bottom
                    largest_label=largest_label+1;
                    label_arr[i][j]=largest_label;
                    //No neighbors to the left or above
                    circumferences[i][j]=4;

                }
                else if(left!=0 && up==0)
                {
                    //only a neighbour to the left
                    label_arr[i][j]=Find(left,labels);
                    circumferences[i][j]=2;
                }
                else if(left==0 && up!=0)
                {
                    //only a neighbour above
                    label_arr[i][j]=Find(up,labels);
                    circumferences[i][j]=2;
                }
                else
                {
                    //neighbours both to the left and above
                    Union(left,up,labels); //link the clusters
                    label_arr[i][j]=Find(left,labels);
                    //As it has both neighbours above and below the circumference remains zero
                }
            }
        }
    }


    int n_labels=labels[N-1];

    //Apply periodic BC's
    int *new_labels_pbc = calloc(sizeof(int), n_labels); // allocate array, initialized to zero
    for(int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++)
        {
            if(S[core][i][j]==1.)
            {
                int x=Find(label_arr[i][j],labels);
                if (new_labels_pbc[x] == 0)
                {
                    new_labels_pbc[0]++;
                    new_labels_pbc[x] = new_labels_pbc[0];
                }
                label_arr[i][j]=new_labels_pbc[x];
            }
        }
    }

    //Nu hen we totaal aantal clusters en kunnen we dust geheugen vrij maken
    total_clusters = new_labels_pbc[0];
    free(new_labels_pbc);


    int sizes_arr[total_clusters];
    int circum_arr[total_clusters];
    //Set all the elements to zero for the sizes
    for(int i=0;i<total_clusters;i++)
    {
        sizes_arr[i]=0;
        circum_arr[i]=0;
    }

    //now add the cluster
    for(int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++)
        {
            if(S[core][i][j]==1.)
            {
                //add one extra element to the label
                int label=label_arr[i][j];
                sizes_arr[label]+=1;
                circum_arr[label]+=circumferences[i][j];
            }
        }
    }

    //now we print the data
    FILE* fr=fopen(filename,"w");
    for (int i=1;i<total_clusters;i++)
    {
        fprintf(fr,"%d\t",sizes_arr[i]);
    }
    fprintf(fr,"\n");
    fclose(fr);

    //now we print the circumferences
    FILE* fp=fopen(filename2,"w");
    for (int i=1;i<total_clusters;i++)
    {
        fprintf(fp,"%d\t",circum_arr[i]);
    }
    fprintf(fp,"\n");
    fclose(fp);


}

void print_S_parallel(char filename[128],int core)
{
    FILE* fr=fopen(filename,"w");
    for (int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++)
        {
            fprintf(fr,"%d\t",S[core][i][j]);
        }
        fprintf(fr,"\n");
    }
    fclose(fr);
}


double Mean_parallel(int core)
{
    double res=0.;
    for (int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++)
        {
            res+=(double)S[core][i][j]/((double)N);
        }
    }
    return res;
}


double Calculate_chi()
{
    double mean=0.;
    //calculate the mean of the correlation array
    for(int i=0;i<correlation_sweeps;i++)
    {
        mean+=Correlation_array[i]/correlation_sweeps;
    }
    //We now calculate the actualy chi(t)
    for(int i=0;i<correlation_sweeps;i++)
    {
        double chi=0.;
        //this makes sure that we dont measure a longer time than possible
        for(int j=0;j<(correlation_sweeps-i);j++)
        {
            //This is thus the correlation times for i
            chi+=(Correlation_array[j]*Correlation_array[j+i]-mean*mean)/(correlation_sweeps-i);
        }
        Correlation_time[i]=chi;
    }
}



int main(int argc, char *argv[])
{
    char *eptr;
    if(argc==1)
    {
        printf("Geen argumenten, gebruikt standaardwaarden\n");
        printf("T=%lf\tstd_h_iex=%lf\n",T);
    }
    else if(argc==2)
    {
        printf("Alleen een waarde voor T ingevoerd\n");
        T=strtod(argv[1],&eptr);
        printf("T=%lf\n",T);
    }


    dsfmt_seed(time(NULL));
    time_t t;
    /* Intializes random number generator */
    srand((unsigned) time(&t));
    //Generate the initial spins and external random field
    Init_S_parallel();
    Init_H();

    //FILE* fr=fopen("magnetisaties.txt","w");
    //set the number of the cores we will run parallel on
    //omp_set_num_threads(N_core2);
    for(int i=0;i<10;i++)
    {
        /*
        for(int j=0;j<N_core;j++)
        {
            //printf("Core %d heeeft een magnetisatie van: %lf\n",j,Mean_parallel(j));
            //fprintf(fr,"%lf\n",Mean_parallel(j));
        }
        //printf("\n\n");
        */
        Perform_Cluster(100);
        printf("Iteratie/100: %d\n",i);

    }

    for(int i=0;i<10;i++)
    {
        /*
        for(int j=0;j<N_core;j++)
        {
            //printf("Core %d heeeft een magnetisatie van: %lf\n",j,Mean_parallel(j));
            //fprintf(fr,"%lf\n",Mean_parallel(j));
        }
        //printf("\n\n");
        */
        Perform_Cluster(1000);
        printf("Iteratie/1000: %d\n",i);
        char file_name1[128];
        sprintf(file_name1,"spincluster_T=%.2f_%d.txt",T,i);
        print_S_parallel(file_name1,0);

    }

    //fclose(fr);

    //3.5 UUURRR VEEELS TE LANG

    /*
    FILE* fr=fopen("fase_test4.txt","w");

    for(int j=0;j<1;j++)
    {
        //Generate a new random field and spins


        for(int p=0;p<1;p++)
        {

            Init_S_parallel();
            Init_H();
            //First perform sweeps to calibrate the system
            Perform_Cluster(correlation_sweeps);

            for(int i=0;i<correlation_sweeps;i++)
            {
                //Correlation_array[i]=Mean_parallel(0);
                Perform_Cluster(1);
                //calculate the fields
                if(i%10==0)
                {
                    for(int k=0;k<N_core;k++)
                    {
                        double mean_grid=Mean_parallel(k);
                        fprintf(fr,"%lf\t%lf\n",H,Mean_parallel(k));
                    }
                }
            }

        }
        H-=0.025;

    }

    fclose(fr);

    */

    //Calculate the actual correlation times
    //Calculate_chi();


    //Save the correlation times to a file

    /*
    FILE* fr=fopen("correlationtime_test8.txt","w");

    for(int i=0;i<correlation_sweeps;i++)
    {
        fprintf(fr,"%lf\t%lf\n",Correlation_time[i],Correlation_array[i]);
    }
    fclose(fr);
    */



    //omp_set_num_threads(N_core);
    /*
        for(int j=0;j<N_core;j++)
            {

                //printf("DEZE RUND OP DEZE CORE:::  %d\n",omp_get_thread_num());

                //we slaan hier de grid op, de omtrekken en de grootes
                //eerst hen we wat naampies nodig

                char file_name1[128];
                char file_name2[128];
                char file_name3[128];

                sprintf(file_name1,"spincluster%d.txt",j);
                sprintf(file_name2,"sizescluster%d.txt",j);
                sprintf(file_name3,"circluster%d.txt",j);
                //This prints out the clustersizes and circumferences
                //printf("cum\n");
                //Cluster_Sizes_Parallel(j,file_name2,file_name3);
                //This prints the actual grids
                print_S_parallel(file_name1,j);



                //int core_1=core_arr[i];
                //int core_2=core_arr[i+1];
                //Cluster_sweep(core_1,core_2,i/2);

            }


    */


    printf("het werrukts\n");

    return 0;
}



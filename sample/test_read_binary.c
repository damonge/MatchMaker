#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fitsio.h>

//Header of the MatchMaker binary file
typedef struct {
  long n_halos_total; //Total number of halos for this run
  long n_halos_here;  //Number of halos in this file
  int n_files;        //Total number of files for this run
  float boxsize;      //Box size in Mpc/h
  float redshift;     //Redshift
  float omega_m;      //Matter parameter (Omega_M)
  float omega_l;      //Cosmological constant parameter (Omega_Lambda)
  float hubble;       //Hubble constant / (100 km/s/Mpc) (h)
  char fill[256-40];  //Meaningless buffer space
} FoFHeader;

//FoF group structure
typedef struct {
  int np;         //Number of particles
  float m_halo;   //Halo mass (in M_sun/h)
  float x_avg[3]; //CM position
  float x_rms[3]; //Position variance
  float v_avg[3]; //CM velocity
  float v_rms[3]; //Velocity dispersion
  float lam[3];   //Angular momentum
  float b;        //second largest to largest eigenvalue of the inertia tensor
  float c;        //smallest to largest eigenvalue of the inertia tensor
  float ea[3];    //Eigenvector for the largest eigenvalue of the inertia tensor
  float eb[3];    //Eigenvector for the second largest eigenvalue of the inertia tensor
  float ec[3];    //Eigenvector for the smallest eigenvalue of the inertia tensor
} FoFHalo;

//This reads a list of FoF groups from a single MatchMaker
//binary file
FoFHalo *read_binary_single(char *fname_in,long *nhalos)
{
  int bsize;
  FoFHeader head;
  FoFHalo *halos;
  FILE *fi=fopen(fname_in,"rb");

  //Read header
  fread(&bsize,sizeof(int),1,fi);
  fread(&head,sizeof(FoFHeader),1,fi);
  fread(&bsize,sizeof(int),1,fi);

  //Allocate memory
  *nhalos=head.n_halos_total;
  halos=malloc(head.n_halos_total*sizeof(FoFHalo));

  //Read halos
  fread(&bsize,sizeof(int),1,fi);
  fread(halos,sizeof(FoFHalo),head.n_halos_total,fi);
  fread(&bsize,sizeof(int),1,fi);
  fclose(fi);

  return halos;
}

void msg_abort(int a,char *msg)
{
  fprintf(stderr,"%s",msg);
  exit(a);
}
 
int main(int argc,char **argv)
{ 
  FoFHalo *halos;
  long nhalos;
  char fname_in[256];

  if(argc!=2) {
    printf("Usage: test_read_binary input_binary\n");
    exit(0);
  }
  sprintf(fname_in,"%s",argv[1]);

  printf("Reading\n");
  halos=read_binary_single(fname_in,&nhalos);
  printf("There are %ld halos\n",nhalos);

  //Do whatever

  free(halos);

  return 0;
}

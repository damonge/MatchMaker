///////////////////////////////////////////////////////////////////////
//                                                                   //
//   Copyright 2012 David Alonso                                     //
//                                                                   //
//                                                                   //
// This file is part of MatchMaker.                                  //
//                                                                   //
// MatchMaker is free software: you can redistribute it and/or       //
// modify it under the terms of the GNU General Public License as    //
// published by the Free Software Foundation, either version 3 of    //
// the License, or (at your option) any later version.               //
//                                                                   //
// MatchMaker is distributed in the hope that it will be useful, but //
// WITHOUT ANY WARRANTY; without even the implied warranty of        //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU //
// General Public License for more details.                          //
//                                                                   //
// You should have received a copy of the GNU General Public License //
// along with MatchMaker.  If not, see                               //
//        <http://www.gnu.org/licenses/>.                            //
//                                                                   //
///////////////////////////////////////////////////////////////////////
#ifndef _COMMON_H_
#define _COMMON_H_

#include <mpi.h>

//////
// Defined in parameters.h
typedef struct {
  //Read from param file
  double np_alloc_factor;  //Alloc factor for particles
  char init_prefix[256];   //Snapshot directory
  char init_dir[256];      //Snapshot prefix
  char output_prefix[256]; //Output prefix
  int output_pernode;
  double dx_extra;         //Buffer for parallelization (in Mpc/h)
  double b_fof;            //Percolation fraction
  int np_min;              //Minimum number of particles

  //Read from IC
  unsigned long long n_part; //Total number of particles
  double boxsize; //Box size
  double omega_m; //Omega_M
  double h;       //Hubble parameter
  double mp;      //Particle mass
  int n_files;    //Number of snapshot files
  double a_init;  //Snapshot time

  //MPI parameters
  int n_nodes;
  int i_node;
  int i_node_left;
  int i_node_right;
  float x_offset;
} Parameters;
Parameters Param;

int read_parameters(char *fname);


//////
// Defined in particle.h
typedef struct {
  float x[3];
  float v[3];   // velocity
#ifdef _LONGIDS
  unsigned long long id;
#else //_LONGIDS
  unsigned int id;
#endif //_LONGIDS
  int fof_id;
  int cll_id;
} Particle;

typedef struct {
  Particle* p;
  Particle* p_back;

  int np_allocated; //Total number of particles allocated
  int np_local; //Total number of particles saved
  int np_indomain; //Total number of particles found in domain
  int np_toleft; //Number of particles that will be passed on to the node on the left
  int np_centre; //Number of particles that won't be passed to nor from any other node
  int np_fromright; //Number of particles passed from the right node
  unsigned long long np_total;
  float np_average;
} Particles;

MPI_Datatype ParticleMPI;
MPI_Datatype HaloMPI;

//////
// Defined in fof.c
typedef struct {
  int np;
  double m_halo;
  double x_avg[3];
  double x_rms[3];
  double v_avg[3];
  double v_rms[3];
} FoFHalo;

FoFHalo *fof_get_halos(int *n_halos_out,Particles *particles);


//////
// Defined in snap_io.c
typedef struct {
  int    np[6];
  double mass[6];
  double time;
  double redshift;
  int    flag_sfr;
  int    flag_feedback;
  unsigned int np_total[6];
  int    flag_cooling;
  int    num_files;
  double boxsize;
  double omega0;
  double omega_lambda;
  double hubble_param; 
  int flag_stellarage;
  int flag_metals;
  unsigned int np_total_highword[6];
  int  flag_entropy_instead_u;
  char fill[60];
} GadgetHeader;

void print_header(GadgetHeader header);
GadgetHeader read_header_multiple(char *dirname,char *prefix);
Particles *read_input_snapshot();
void write_halos(int n_fof,FoFHalo *fof);


//////
// Defined in msg.c
void msg_init();
void msg_printf(const char *fmt, ...);
void msg_abort(const int errret,const char *fmt, ...);

#endif //_COMMON_H_

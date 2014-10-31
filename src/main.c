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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <math.h>

#include "common.h"

int mpi_init(int* p_argc,char ***p_argv);

int main(int argc,char* argv[])
{
  mpi_init(&argc,&argv);
  msg_init();

  if(argc<2)
    msg_abort(1,"Error: Parameter file not specified. cola_halo param.init\n");
  read_parameters(argv[1]);

  msg_printf("Reading particles\n");
  Particles* particles=read_input_snapshot();

  msg_printf("Getting halos\n");
  int n_halos;
  FoFHalo *fh=fof_get_halos(&n_halos,particles);

  msg_printf("Writing output\n");
  write_halos(n_halos,fh);

  free(particles->p);
  free(particles->p_back);
  free(particles);

  MPI_Finalize();
  return 0;
}

int mpi_init(int* p_argc,char*** p_argv)
{
  MPI_Init(p_argc,p_argv);

  //Define MPI Particle and Halo structures
  MPI_Aint p_off[5];
  MPI_Datatype p_ot[5];
  int p_bc[5];
  p_bc[0]=3; //3*x
  p_bc[1]=3; //3*v
  p_bc[2]=1; //id
  p_bc[3]=1; //fof_id
  p_bc[4]=1; //cll_id
  p_off[0]=offsetof(Particle,x);
  p_off[1]=offsetof(Particle,v);
  p_off[2]=offsetof(Particle,id);
  p_off[3]=offsetof(Particle,fof_id);
  p_off[4]=offsetof(Particle,cll_id);
  p_ot[0]=MPI_FLOAT;
  p_ot[1]=MPI_FLOAT;
  p_ot[2]=MPI_UNSIGNED_LONG_LONG;
  p_ot[3]=MPI_INT;
  p_ot[4]=MPI_INT;
  MPI_Type_struct(5,p_bc,p_off,p_ot,&ParticleMPI);
  MPI_Type_commit(&ParticleMPI);

  MPI_Aint h_off[6];
  MPI_Datatype h_ot[6];
  int h_bc[6];
  h_bc[0]=1; // np
  h_bc[1]=1; // m
  h_bc[2]=3; // 3*x_avg
  h_bc[3]=3; // 3*x_rms
  h_bc[4]=3; // 3*v_avg
  h_bc[5]=3; // 3*v_rms
  h_off[0]=offsetof(FoFHalo,np);
  h_off[1]=offsetof(FoFHalo,m_halo);
  h_off[2]=offsetof(FoFHalo,x_avg);
  h_off[3]=offsetof(FoFHalo,x_rms);
  h_off[4]=offsetof(FoFHalo,v_avg);
  h_off[5]=offsetof(FoFHalo,v_rms);
  h_ot[0]=MPI_INT;
  h_ot[1]=MPI_DOUBLE;
  h_ot[2]=MPI_DOUBLE;
  h_ot[3]=MPI_DOUBLE;
  h_ot[4]=MPI_DOUBLE;
  h_ot[5]=MPI_DOUBLE;
  MPI_Type_struct(6,h_bc,h_off,h_ot,&HaloMPI);
  MPI_Type_commit(&HaloMPI);

  /*
  MPI_Aint ext_f,ext_llu,ext_int,ext_d;
  MPI_Type_extent(MPI_FLOAT,&ext_f);
  MPI_Type_extent(MPI_DOUBLE,&ext_d);
  MPI_Type_extent(MPI_INT,&ext_int);
  MPI_Type_extent(MPI_UNSIGNED_LONG_LONG,&ext_llu);
  printf("%d %d\n",p_off[0],0);
  printf("%d %d\n",p_off[1],3*ext_f);
  printf("%d %d\n",p_off[2],6*ext_f);
  printf("%d %d\n",p_off[3],6*ext_f+ext_llu);
  printf("%d %d\n",p_off[4],6*ext_f+ext_llu+ext_int);
  printf("%d %d\n",h_off[0],0);
  printf("%d %d\n",h_off[1],ext_int);
  printf("%d %d\n",h_off[2],ext_int+ext_d);
  printf("%d %d\n",h_off[3],ext_int+4*ext_d);
  printf("%d %d\n",h_off[4],ext_int+7*ext_d);
  printf("%d %d\n",h_off[5],ext_int+10*ext_d);
  */
  return 0;
}

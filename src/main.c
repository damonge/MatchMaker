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

static int mpi_init(int* p_argc,char*** p_argv)
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
#ifdef _LONGIDS
  p_ot[2]=MPI_UNSIGNED_LONG_LONG;
#else //_LONGIDS
  p_ot[2]=MPI_UNSIGNED;
#endif //_LONGIDS
#ifdef _LONG_INT
  p_ot[3]=MPI_LONG;
  p_ot[4]=MPI_LONG;
#else //_LONG_INT
  p_ot[3]=MPI_INT;
  p_ot[4]=MPI_INT;
#endif //_LONG_INT
  MPI_Type_create_struct(5,p_bc,p_off,p_ot,&ParticleMPI);
  MPI_Type_commit(&ParticleMPI);

  MPI_Aint h_off[12];
  MPI_Datatype h_ot[12];
  int h_bc[12];
  h_bc[0]=1; // np
  h_bc[1]=1; // m
  h_bc[2]=3; // 3*x_avg
  h_bc[3]=3; // 3*x_rms
  h_bc[4]=3; // 3*v_avg
  h_bc[5]=3; // 3*v_rms
  h_bc[6]=3; // 3*lam
  h_bc[7]=1; // b
  h_bc[8]=1; // c
  h_bc[9]=3; // 3*ea
  h_bc[10]=3; // 3*ea
  h_bc[11]=3; // 3*ea
  h_off[0]=offsetof(FoFHalo,np);
  h_off[1]=offsetof(FoFHalo,m_halo);
  h_off[2]=offsetof(FoFHalo,x_avg);
  h_off[3]=offsetof(FoFHalo,x_rms);
  h_off[4]=offsetof(FoFHalo,v_avg);
  h_off[5]=offsetof(FoFHalo,v_rms);
  h_off[6]=offsetof(FoFHalo,lam);
  h_off[7]=offsetof(FoFHalo,b);
  h_off[8]=offsetof(FoFHalo,c);
  h_off[9]=offsetof(FoFHalo,ea);
  h_off[10]=offsetof(FoFHalo,eb);
  h_off[11]=offsetof(FoFHalo,ec);
  h_ot[0]=MPI_INT;
  h_ot[1]=MPI_FLOAT;
  h_ot[2]=MPI_FLOAT;
  h_ot[3]=MPI_FLOAT;
  h_ot[4]=MPI_FLOAT;
  h_ot[5]=MPI_FLOAT;
  h_ot[6]=MPI_FLOAT;
  h_ot[7]=MPI_FLOAT;
  h_ot[8]=MPI_FLOAT;
  h_ot[9]=MPI_FLOAT;
  h_ot[10]=MPI_FLOAT;
  h_ot[11]=MPI_FLOAT;
  MPI_Type_create_struct(12,h_bc,h_off,h_ot,&HaloMPI);
  MPI_Type_commit(&HaloMPI);

  return 0;
}

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
  lint n_halos;
  FoFHalo *fh=fof_get_halos(&n_halos,particles);

  msg_printf("Writing output\n");
  write_halos(n_halos,fh);

  free(particles->p);
  free(particles->p_back);
  free(particles);

  MPI_Finalize();
  return 0;
}

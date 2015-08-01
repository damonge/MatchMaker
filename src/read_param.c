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
#include <math.h>
#include <assert.h>

#include "common.h"

static int read_parameter_file(char const *fname,Parameters* const param);
static void bcast_string(char* string,int len);

int read_parameters(char *filename)
{
  //Only process 0 reads params then communicates
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  if(myrank==0) {
    int ret=read_parameter_file(filename,&Param);
    if(ret!=0)
      msg_abort(1001,"Error: Unable to read parameter file: %s\n",filename);
  }

  // Share parameters with other nodes
  MPI_Bcast(&Param,sizeof(Parameters),MPI_BYTE,0,MPI_COMM_WORLD);

  bcast_string(Param.init_prefix,256);
  bcast_string(Param.init_dir,256);
  bcast_string(Param.output_prefix,256);

  //Get MPI parameters
  MPI_Comm_rank(MPI_COMM_WORLD,&(Param.i_node));
  MPI_Comm_size(MPI_COMM_WORLD,&(Param.n_nodes));
  Param.x_offset=Param.i_node*((float)(Param.boxsize)/Param.n_nodes);
  Param.i_node_left=(Param.i_node-1+Param.n_nodes)%Param.n_nodes;
  Param.i_node_right=(Param.i_node+1)%Param.n_nodes;

  return 0;
}

static int linecount(FILE *f)
{
  //////
  // Counts #lines from file
  int i0=0;
  char ch[1000];
  while((fgets(ch,sizeof(ch),f))!=NULL) {
    i0++;
  }
  return i0;
}

static int read_parameter_file(char const *fname,Parameters *param)
{
  FILE *fi;
  int n_lin,ii;
  
  //Read parameters from file
  fi=fopen(fname,"r");
  if(fi==NULL)
    msg_abort(1000, "Error: couldn't open file %s\n",fname);

  param->output_format=0;

  n_lin=linecount(fi);
  rewind(fi);
  for(ii=0;ii<n_lin;ii++) {
    char s0[512],s1[64],s2[256];
    if(fgets(s0,sizeof(s0),fi)==NULL)
      msg_abort(1001,"Error reading line %d, file %s\n",ii+1,fname);
    if((s0[0]=='#')||(s0[0]=='\n')) continue;
    int sr=sscanf(s0,"%s %s",s1,s2);
    if(sr!=2)
      msg_abort(1002,"Error reading line %d, file %s\n",ii+1,fname);
    
    if(!strcmp(s1,"np_alloc_factor="))
      param->np_alloc_factor=atof(s2);
    else if(!strcmp(s1,"init_dir="))
      sprintf(param->init_dir,"%s",s2);
    else if(!strcmp(s1,"init_prefix="))
      sprintf(param->init_prefix,"%s",s2);
    else if(!strcmp(s1,"output_prefix="))
      sprintf(param->output_prefix,"%s",s2);
    else if(!strcmp(s1,"output_pernode="))
      param->output_pernode=atoi(s2);
    else if(!strcmp(s1,"output_format=")) {
      if(!strncmp(s2,"ASCII",4))
	param->output_format=0;
      else if(!strncmp(s2,"FITS",4))
	param->output_format=1;
      else if(!strncmp(s2,"BINARY",4))
	param->output_format=2;
      else
	msg_abort(1002,"Unrecognized format %s\n",s2);
    }
    else if(!strcmp(s1,"dx_extra="))
      param->dx_extra=atof(s2);
    else if(!strcmp(s1,"b_fof="))
      param->b_fof=atof(s2);
    else if(!strcmp(s1,"np_min="))
      param->np_min=atoi(s2);
    else
      msg_printf("Unknown parameter %s\n",s1);
  }
  fclose(fi);

  //Quote parameters
  msg_printf("Read parameters\n");
  msg_printf("  np_alloc_factor= %.3lf\n",param->np_alloc_factor);
  msg_printf("  init_dir= %s\n",param->init_dir);
  msg_printf("  init_prefix= %s\n",param->init_prefix);
  msg_printf("  output_prefix= %s\n",param->output_prefix);
  msg_printf("  output_pernode= %d\n",param->output_pernode);
  msg_printf("  output_format= %d\n",param->output_format);
  msg_printf("  dx_extra= %.3lf\n",param->dx_extra);
  msg_printf("  b_fof= %.3lf\n",param->b_fof);
  msg_printf("  np_min= %d\n",param->np_min);
  msg_printf("\n");
  
  //Read IC header and check
  GadgetHeader head=read_header_multiple(param->init_dir,param->init_prefix);
  print_header(head);
  param->n_part=(unsigned long long)(head.np_total[1])+
    ((unsigned long long)(head.np_total_highword[1]) << 32);
  param->boxsize=head.boxsize;
  param->omega_m=head.omega0;
  param->h=head.hubble_param;
  param->n_files=head.num_files;
  param->a_init=1./(1+head.redshift);
  param->mp=head.mass[1];

  //Sanity checks:
  //Signs and values
  if(param->np_alloc_factor<0 || param->np_alloc_factor>3)
    msg_abort(123,"Allocation factor seems weird %lf\n",param->np_alloc_factor);
  if(param->b_fof<0 || param->b_fof>1)
    msg_abort(123,"Percolation fraction seems weird %lf\n",param->b_fof);

  //Buffer size
  int nnode; MPI_Comm_size(MPI_COMM_WORLD,&nnode);
  if(param->dx_extra>=0.5*param->boxsize/nnode)
    msg_abort(1001,"Buffer size might be too big!\n"
	      "Reduce it or the number of cores\n");
  if(param->dx_extra<4.0)
    msg_printf("WARNING: buffer size might be too small!\n");

  //Check flatness
  if(fabs(1-head.omega0-head.omega_lambda)>0.001)
    msg_printf("WARNING: non-flat cosmology!!\n");

  //Check particle mass
  double mp_here=27.7459457*head.omega0*pow(head.boxsize,3)/param->n_part;
  if(fabs(mp_here/param->mp-1)>=0.01)
    msg_abort(1002,"Particle mass doesn't fit... %lE!=%lE\n",
	      mp_here,param->mp);

  msg_printf("\n");

  return 0;
}

static void bcast_string(char* pstring,int len)
{
  const int n=len;

  const int ret2=MPI_Bcast(pstring,n,MPI_CHAR,0,MPI_COMM_WORLD);
  assert(ret2==MPI_SUCCESS);
}

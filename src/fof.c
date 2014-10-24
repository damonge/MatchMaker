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

#include <dirent.h>

#include "common.h"

static int Ngrid[3];
static int Ngrid_tot;
static float Ac[3];
static float I_ac[3];
static float Dfof;
static float D2fof;
static float Lbox;
static float Lbox_half;

typedef struct {
  int np_in;
  int *ids_in;
  int np_out;
  int *ids_out;
} Cell;

typedef struct {
  int np;
  int *ids;
} FoFGroup;

static int compare_fof(const void *p1,const void *p2)
{
  const FoFGroup *f1=p1;
  const FoFGroup *f2=p2;
  
  if(f1->np>f2->np)
    return -1;
  else if(f1->np<f2->np)
    return 1;
  else
    return 0;
}

static void init_fof(void)
{
  int ax;
  float ipd=Param.boxsize/pow((double)(Param.n_part),1./3.);
  float dx[3];
  
  dx[0]=Param.boxsize/Param.n_nodes+Param.dx_extra;
  dx[1]=Param.boxsize;
  dx[2]=Param.boxsize;

  for(ax=0;ax<3;ax++) {
    Ngrid[ax]=(int)(dx[ax]/ipd)+1;
    Ac[ax]=dx[ax]/Ngrid[ax];
    I_ac[ax]=1./Ac[ax];
  }

  Ngrid_tot=Ngrid[0]*Ngrid[1]*Ngrid[2];
  Dfof=ipd*Param.b_fof;
  D2fof=Dfof*Dfof;
  Lbox=(float)(Param.boxsize);
  Lbox_half=(float)(Param.boxsize/2);
#ifdef _DEBUG
  msg_printf(" Grid size = %d (%d x %d x %d)\n",
	 Ngrid_tot,Ngrid[0],Ngrid[1],Ngrid[2]);
  msg_printf(" Cell sizes : %lf x %lf x %lf\n",
	 Ac[0],Ac[1],Ac[2]);
  msg_printf(" FoF separation : %lf\n",Dfof);
#endif //_DEBUG
}

static void get_neighbors(int ip0,Cell *cll,Particle *p,int fof_id)
{
  int ii;
  Particle *p0=&(p[ip0]);
  Cell c0=cll[p0->cll_id];
  float *x0=p0->x;

  //Couple with particles in the same cell
  for(ii=0;ii<c0.np_in;ii++) {
    int id=c0.ids_in[ii];

    if(id!=ip0) {
      if(p[id].fof_id==0) {
	float *x1=p[id].x;
	float dx;
	
	dx=fabs(x0[0]-x1[0]);
	if(dx<=Dfof) {
	  float dy,dz,d2;
	  dy=fabs(x0[1]-x1[1]);
	  dz=fabs(x0[2]-x1[2]);
	  if(dy>Lbox_half) dy=Lbox-dy;
	  if(dz>Lbox_half) dz=Lbox-dz;
	  d2=dx*dx+dy*dy+dz*dz;
	  if(d2<=D2fof) {
	    if(p0->fof_id==0) //New halo!
	      p0->fof_id=fof_id;
	    
	    p[id].fof_id=fof_id;
	    get_neighbors(id,cll,p,fof_id);
	  }
	}
      }
    }
  }

  //Couple with particles in the same cell
  for(ii=0;ii<c0.np_out;ii++) {
    int id=c0.ids_out[ii];

    if(id!=ip0) {
      if(p[id].fof_id==0) {
	float dx;
	float *x1=p[id].x;
	
	dx=fabs(x0[0]-x1[0]);
	if(dx<=Dfof) {
	  float dy,dz,d2;
	  dy=fabs(x0[1]-x1[1]);
	  dz=fabs(x0[2]-x1[2]);
	  if(dy>Lbox_half) dy=Lbox-dy;
	  if(dz>Lbox_half) dz=Lbox-dz;
	  d2=dx*dx+dy*dy+dz*dz;
	  if(d2<=D2fof) {
	    if(p0->fof_id==0) //New halo!
	      p0->fof_id=fof_id;
	    
	    p[id].fof_id=fof_id;
	    get_neighbors(id,cll,p,fof_id);
	  }
	}
      }
    }
    else
      msg_abort(123,"This shouldn't have happened\n");
  }
}

static Cell *gather_particles_in_cells(Particles *particles) 
{
  int i;
  Cell *cll=malloc(Ngrid_tot*sizeof(Cell));
  if(cll==NULL)
    msg_abort(123,"Failed to allocate memory for NGBS cells\n");
  
  for(i=0;i<Ngrid_tot;i++) {
    cll[i].np_in=0;
    cll[i].np_out=0;
  }

  Particle *p=particles->p;
  for(i=0;i<particles->np_local;i++) {
    int ax;
    int ic[3],ic_close[3];
    int icell;

    for(ax=0;ax<3;ax++) {
      float x=p[i].x[ax];
      float ix=x*I_ac[ax];
      int ic_in=(int)ix;
      int ic_cl=(int)(ix+0.5);
      float dist_close=fabs(ic_cl*Ac[ax]-x);

      if(ic_in>=Ngrid[ax]) {
	if(ax==0)
	  ic[ax]=Ngrid[ax]-1;
	else
	  ic[ax]=ic_in-Ngrid[ax];
      }
      else if(ic_in<0) { 
	if(ax==0)
	  ic[ax]=0;
	else
	  ic[ax]=ic_in+Ngrid[ax];
      }
      else ic[ax]=ic_in;

      if(dist_close<=Dfof) {
	ic_close[ax]=ic[ax]-1+2*(ic_cl-ic_in);
	if(ic_close[ax]>=Ngrid[ax]) {
	  if(ax==0)
	    ic_close[ax]=-1;
	  else
	    ic_close[ax]-=Ngrid[ax];
	}
	else if(ic_close[ax]<0) {
	  if(ax==0)
	    ic_close[ax]=-1;
	  else
	    ic_close[ax]+=Ngrid[ax];
	}
      }
      else
	ic_close[ax]=-1;
    }
    icell=ic[0]+Ngrid[0]*(ic[1]+Ngrid[1]*ic[2]);
    cll[icell].np_in++;
    p[i].cll_id=icell;
    p[i].fof_id=0;

    if(ic_close[0]>=0) {
      icell=ic_close[0]+Ngrid[0]*(ic[1]+Ngrid[1]*ic[2]);
      cll[icell].np_out++;
      if(ic_close[1]>=0) {
	icell=ic[0]+Ngrid[0]*(ic_close[1]+Ngrid[1]*ic[2]);
	cll[icell].np_out++;
	icell=ic_close[0]+Ngrid[0]*(ic_close[1]+Ngrid[1]*ic[2]);
	cll[icell].np_out++;
	if(ic_close[2]>=0) {
	  icell=ic[0]+Ngrid[0]*(ic[1]+Ngrid[1]*ic_close[2]);
	  cll[icell].np_out++;
	  icell=ic_close[0]+Ngrid[0]*(ic[1]+Ngrid[1]*ic_close[2]);
	  cll[icell].np_out++;
	  icell=ic[0]+Ngrid[0]*(ic_close[1]+Ngrid[1]*ic_close[2]);
	  cll[icell].np_out++;
	  icell=ic_close[0]+Ngrid[0]*(ic_close[1]+Ngrid[1]*ic_close[2]);
	  cll[icell].np_out++;
	}
      }
      else if(ic_close[2]>=0) {
	icell=ic[0]+Ngrid[0]*(ic[1]+Ngrid[1]*ic_close[2]);
	cll[icell].np_out++;
	icell=ic_close[0]+Ngrid[0]*(ic[1]+Ngrid[1]*ic_close[2]);
	cll[icell].np_out++;
      }
    }
    else if(ic_close[1]>=0) {
      icell=ic[0]+Ngrid[0]*(ic_close[1]+Ngrid[1]*ic[2]);
      cll[icell].np_out++;
      if(ic_close[2]>=0) {
	icell=ic[0]+Ngrid[0]*(ic[1]+Ngrid[1]*ic_close[2]);
	cll[icell].np_out++;
	icell=ic[0]+Ngrid[0]*(ic_close[1]+Ngrid[1]*ic_close[2]);
	cll[icell].np_out++;
      }
    }
    else if(ic_close[2]>=0) {
      icell=ic[0]+Ngrid[0]*(ic[1]+Ngrid[1]*ic_close[2]);
      cll[icell].np_out++;
    }
  }

  for(i=0;i<Ngrid_tot;i++) {
    int np_in=cll[i].np_in;
    int np_out=cll[i].np_out;
    if(np_in>0) {
      cll[i].ids_in=malloc(np_in*sizeof(int));
      if(cll[i].ids_in==NULL)
	msg_abort(123,"Failed to allocate ids in cells\n");
    }
    if(np_out>0) {
      cll[i].ids_out=malloc(np_out*sizeof(int));
      if(cll[i].ids_out==NULL)
	msg_abort(123,"Failed to allocate ids in cells\n");
    }
    cll[i].np_in=0;
    cll[i].np_out=0;
  }

  for(i=0;i<particles->np_local;i++) {
    int ax;
    int ic[3],ic_close[3];
    int icell;

    for(ax=0;ax<3;ax++) {
      float x=p[i].x[ax];
      float ix=x*I_ac[ax];
      int ic_in=(int)ix;
      int ic_cl=(int)(ix+0.5);
      float dist_close=fabs(ic_cl*Ac[ax]-x);

      if(ic_in>=Ngrid[ax]) {
	if(ax==0)
	  ic[ax]=Ngrid[ax]-1;
	else
	  ic[ax]=ic_in-Ngrid[ax];
      }
      else if(ic_in<0) { 
	if(ax==0)
	  ic[ax]=0;
	else
	  ic[ax]=ic_in+Ngrid[ax];
      }
      else ic[ax]=ic_in;

      if(dist_close<=Dfof) {
	ic_close[ax]=ic[ax]-1+2*(ic_cl-ic_in);
	if(ic_close[ax]>=Ngrid[ax]) {
	  if(ax==0)
	    ic_close[ax]=-1;
	  else
	    ic_close[ax]-=Ngrid[ax];
	}
	else if(ic_close[ax]<0) {
	  if(ax==0)
	    ic_close[ax]=-1;
	  else
	    ic_close[ax]+=Ngrid[ax];
	}
      }
      else
	ic_close[ax]=-1;
    }
    icell=ic[0]+Ngrid[0]*(ic[1]+Ngrid[1]*ic[2]);
    if(icell!=p[i].cll_id)
      msg_abort(123,"Error in cell assignment, node %d\n",Param.i_node);
    cll[icell].ids_in[cll[icell].np_in++]=i;

    if(ic_close[0]>=0) {
      icell=ic_close[0]+Ngrid[0]*(ic[1]+Ngrid[1]*ic[2]);
      cll[icell].ids_out[cll[icell].np_out++]=i;
      if(ic_close[1]>=0) {
	icell=ic[0]+Ngrid[0]*(ic_close[1]+Ngrid[1]*ic[2]);
	cll[icell].ids_out[cll[icell].np_out++]=i;
	icell=ic_close[0]+Ngrid[0]*(ic_close[1]+Ngrid[1]*ic[2]);
	cll[icell].ids_out[cll[icell].np_out++]=i;
	if(ic_close[2]>=0) {
	  icell=ic[0]+Ngrid[0]*(ic[1]+Ngrid[1]*ic_close[2]);
	  cll[icell].ids_out[cll[icell].np_out++]=i;
	  icell=ic_close[0]+Ngrid[0]*(ic[1]+Ngrid[1]*ic_close[2]);
	  cll[icell].ids_out[cll[icell].np_out++]=i;
	  icell=ic[0]+Ngrid[0]*(ic_close[1]+Ngrid[1]*ic_close[2]);
	  cll[icell].ids_out[cll[icell].np_out++]=i;
	  icell=ic_close[0]+Ngrid[0]*(ic_close[1]+Ngrid[1]*ic_close[2]);
	  cll[icell].ids_out[cll[icell].np_out++]=i;
	}
      }
      else if(ic_close[2]>=0) {
	icell=ic[0]+Ngrid[0]*(ic[1]+Ngrid[1]*ic_close[2]);
	cll[icell].ids_out[cll[icell].np_out++]=i;
	icell=ic_close[0]+Ngrid[0]*(ic[1]+Ngrid[1]*ic_close[2]);
	cll[icell].ids_out[cll[icell].np_out++]=i;
      }
    }
    else if(ic_close[1]>=0) {
      icell=ic[0]+Ngrid[0]*(ic_close[1]+Ngrid[1]*ic[2]);
      cll[icell].ids_out[cll[icell].np_out++]=i;
      if(ic_close[2]>=0) {
	icell=ic[0]+Ngrid[0]*(ic[1]+Ngrid[1]*ic_close[2]);
	cll[icell].ids_out[cll[icell].np_out++]=i;
	icell=ic[0]+Ngrid[0]*(ic_close[1]+Ngrid[1]*ic_close[2]);
	cll[icell].ids_out[cll[icell].np_out++]=i;
      }
    }
    else if(ic_close[2]>=0) {
      icell=ic[0]+Ngrid[0]*(ic[1]+Ngrid[1]*ic_close[2]);
      cll[icell].ids_out[cll[icell].np_out++]=i;
    }
  }

  return cll;
}

static FoFGroup *assign_particles_to_fof(Particles *particles,Cell *cll,int *n_fof_out)
{
  int i;
  int n_fof=1;
  Particle *p=particles->p;
  
  //Do FoF
  for(i=0;i<particles->np_indomain;i++) {
    if(p[i].fof_id==0) {
      get_neighbors(i,cll,p,n_fof);
      if(p[i].fof_id==n_fof) //Halo found
	n_fof++;
    }
  }
  n_fof--;

  //Free memory from cells
  for(i=0;i<Ngrid_tot;i++) {
    if(cll[i].np_in>0)
      free(cll[i].ids_in);
    if(cll[i].np_out>0)
      free(cll[i].ids_out);
  }
  free(cll);

  //Allocate FoF groups
  FoFGroup *fof=malloc(n_fof*sizeof(FoFGroup));
  if(fof==NULL)
    msg_abort(123,"Unable to allocate memory for FoF groups\n");
  for(i=0;i<n_fof;i++)
    fof[i].np=0;

  //Compare particles
  //Receive buffer particles from left and send to right
  int tag=200;
  MPI_Status stat;
  MPI_Sendrecv(&(particles->p[particles->np_indomain]),particles->np_fromright,
			ParticleMPI,Param.i_node_right,tag,
			particles->p_back,particles->np_toleft,
			ParticleMPI,Param.i_node_left,tag,
			MPI_COMM_WORLD,&stat);

  //Resolve redundancies
  Particle *p_fromleft=particles->p_back;
  for(i=0;i<particles->np_toleft;i++) {
    if(p_fromleft[i].fof_id!=0)
      p[i].fof_id=0;
  }

  //Allocate id arrays in halos
  for(i=0;i<particles->np_local;i++) {
    if(p[i].fof_id!=0)
      fof[p[i].fof_id-1].np++;
  }

  int n_fof_true=0;
  for(i=0;i<n_fof;i++) {
    int np=fof[i].np;
    if(np>0) {
      n_fof_true++;
      fof[i].ids=malloc(np*sizeof(int));
      if(fof[i].ids==NULL)
	msg_abort(123,"Failed to allocate memory for particle ids in FoF groups\n");
      fof[i].np=0;
    }
  }
#ifdef _DEBUG
  printf("In node %d there were initially %d FoF groups "
	 "but only %d after talking to node %d. Solved %d redundancies\n",
	 Param.i_node,n_fof,n_fof_true,Param.i_node_left,n_fof-n_fof_true);
#endif //_DEBUG  

  //Assign particle ids to FoF groups
  for(i=0;i<particles->np_local;i++) {
    if(p[i].fof_id!=0) {
      int id=i;
      int i_fof=p[i].fof_id-1;
      fof[i_fof].ids[fof[i_fof].np]=id;
      fof[i_fof].np++;
    }
  }

  //Sort by number of particles
  qsort(fof,n_fof,sizeof(FoFGroup),compare_fof);

  *n_fof_out=n_fof_true;
  
  return fof;
}

static FoFHalo *get_halos(int *n_halos_out,int n_fof,FoFGroup *fg,Particles *particles)
{
  int i,n_halos;
  Particle *p=particles->p;
  FoFHalo *fh;

  int np_current=fg[0].np;
  n_halos=0;
  for(i=0;i<n_fof;i++) {
    int np=fg[i].np;
    if(np>np_current)
      msg_abort(123,"Ordering seems to be wrong!\n");
    if(np>=Param.np_min)
      n_halos++;
    np_current=np;
  }
  if(n_halos==0)
    msg_abort(123,"Couldn't find any halos above np_min=%d\n",Param.np_min);
#ifdef _DEBUG
  printf("In node %d there were initially %d FoF groups "
	 "but only %d of them have >%d particles\n",
	 Param.i_node,n_fof,n_halos,Param.np_min);
#endif //_DEBUG  
  
  fh=malloc(n_halos*sizeof(FoFHalo));
  if(fh==NULL)
    msg_abort(123,"Unable to allocate memory for halos\n");  

  for(i=0;i<n_halos;i++) {
    int j;
    int np=fg[i].np;

    fh[i].np=np;
    fh[i].m_halo=np*Param.mp;
    for(j=0;j<3;j++) {
      fh[i].x_avg[j]=0;
      fh[i].x_rms[j]=0;
      fh[i].v_avg[j]=0;
      fh[i].v_rms[j]=0;
    }
    
    for(j=0;j<np;j++) {
      int ax;
      int ip=fg[i].ids[j];
      for(ax=0;ax<3;ax++) {
	float x=p[ip].x[ax];
	float v=p[ip].v[ax];

	fh[i].x_avg[ax]+=x;
	fh[i].x_rms[ax]+=x*x;
	fh[i].v_avg[ax]+=v;
	fh[i].v_rms[ax]+=v*v;
      }
    }

    for(j=0;j<3;j++) {
      fh[i].x_avg[j]/=np;
      fh[i].x_rms[j]=sqrt(fh[i].x_rms[j]/np-fh[i].x_avg[j]*fh[i].x_avg[j]);
      fh[i].v_avg[j]/=np;
      fh[i].v_rms[j]=sqrt(fh[i].v_rms[j]/np-fh[i].v_avg[j]*fh[i].v_avg[j]);
    }
    fh[i].x_avg[0]+=Param.x_offset;

    if(fg[i].np>0)
      free(fg[i].ids);
  }

  free(fg);

  *n_halos_out=n_halos;
  return fh;
}

FoFHalo *fof_get_halos(int *n_halos_out,Particles *particles)
{
  FoFGroup *fof_g;
  FoFHalo *fof_h;
  Cell *cll;
  int n_fof,n_halos;

  init_fof();
  msg_printf(" Gathering particles in cells for NB searching\n");
  cll=gather_particles_in_cells(particles);
  msg_printf(" Matchmaking\n");
  fof_g=assign_particles_to_fof(particles,cll,&n_fof);
  msg_printf(" Computing halo properties\n");
  fof_h=get_halos(&n_halos,n_fof,fof_g,particles);
  msg_printf(" Done\n\n");

  *n_halos_out=n_halos;

  return fof_h;
}

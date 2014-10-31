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
#include <dirent.h>

#include "common.h"

static void bad_block(void)
{
  msg_abort(123,"Bad block!\n");
}

static float wrap_float(float x,float lb)
{
  if(x<0) return x+lb;
  else if(x>=lb) return x-lb;
  else return x;
}

static void my_fread(void *ptr,size_t size,size_t count,FILE *stream)
{
  size_t stat=fread(ptr,size,count,stream);
  if(stat!=count) {
    msg_abort(123,"Error freading\n");
  }
}

static FILE *my_fopen(char *fname,char *mode)
{
  FILE *f=fopen(fname,mode);
  if(f==NULL) {
    msg_abort(123,"Couldn't open file %s\n",fname);
  }

  return f;
}

static int compare_halo(const void *p1,const void *p2)
{
  const FoFHalo *f1=p1;
  const FoFHalo *f2=p2;
  
  if(f1->np>f2->np)
    return -1;
  else if(f1->np<f2->np)
    return 1;
  else
    return 0;
}

static int compare_parts(const void *e1,const void *e2)
{
  const Particle *p1=e1;
  const Particle *p2=e2;

  if(p1->x[0]<p2->x[0])
    return -1;
  else if(p1->x[0]>p2->x[0])
    return 1;
  else
    return 0;
}

static Particles *allocate_particles(void)
{
  Particles *particles=malloc(sizeof(Particles));
  const int nnode=Param.n_nodes;

  //Check alloc factor
  double alloc_factor_buffer=(1+nnode*Param.dx_extra/Param.boxsize);
  double alloc_factor_true=Param.np_alloc_factor*alloc_factor_buffer;

  const int np_alloc=(int)((double)(alloc_factor_true*Param.n_part)/nnode);
  particles->p=malloc(sizeof(Particle)*np_alloc);
  if(particles->p==0)
    msg_abort(0010,"Error: Failed to allocate memory for particles\n");
#ifdef _DEBUG
  printf("Node %d: %d Mbytes allocated for %d particles (alloc_factor= %.2lf)\n",
	 Param.i_node,(int)(sizeof(Particle)*np_alloc/(1024*1024.)),
	 np_alloc,Param.np_alloc_factor);
#endif //_DEBUG  

  particles->np_allocated=np_alloc;

  particles->np_total=Param.n_part;
  particles->np_average=(float)((double)(Param.n_part)/nnode);

  return particles;
}

void print_header(GadgetHeader header)
{
  int i;

  msg_printf("Snapshot header\n");
  msg_printf("  #part. in this file: ");
  for(i=0;i<6;i++) msg_printf("%d, ",header.np[i]);
  msg_printf("\n");

  msg_printf("  #part. in total:     ");
  for(i=0;i<6;i++) {
    unsigned long long npa=(unsigned long long)(header.np_total[i]);
    unsigned long long npb=((unsigned long long)(header.np_total_highword[i]) << 32);
    unsigned long long np=npa+npb;
    msg_printf("%llu, ",np);
  }
  msg_printf("\n");

  msg_printf("  #part. mass (10^10 M_sun/h) :");
  for(i=0;i<6;i++) msg_printf("%.2lE, ",header.mass[i]);
  msg_printf("\n");

  msg_printf("  Redshift : %.2lE, Scale factor: %.2lE\n",
	 header.redshift,header.time);

  msg_printf("  O_M = %.2lf, O_L = %.2lf, h = %.2lf\n",
	 header.omega0,header.omega_lambda,header.hubble_param);

  msg_printf("  L_box = %.2lE Mpc/h\n",header.boxsize);

  msg_printf("  Total number of files: %d\n",header.num_files);
}

GadgetHeader read_header_multiple(char *dirname,char *prefix)
{
  int i;
  struct dirent *ent;
  DIR *dir;
  char **fname_array;
  int n_files=0;
  int prefix_len=strlen(prefix);

  dir=opendir(dirname);
  //Find number of files
  if(dir==NULL) {
    msg_abort(123,"wrong directory %s\n",dirname);
  }
  while((ent=readdir(dir))!=NULL) {
    if(!strncmp(ent->d_name,prefix,prefix_len))
      n_files++;
  }
  closedir(dir);
  if(n_files!=0)
    msg_printf("Found %d files\n",n_files);
  else
    msg_abort(123,"Found no input snapshots %s %s\n",dirname,prefix);

  //Collect filenames
  fname_array=(char **)malloc(n_files*sizeof(char *));
  for(i=0;i<n_files;i++)
    fname_array[i]=(char *)malloc(256*sizeof(char));

  dir=opendir(dirname);
  n_files=0;
  while((ent=readdir(dir))!=NULL) {
    if(!strncmp(ent->d_name,prefix,prefix_len)) {
      sprintf(fname_array[n_files],"%s/%s",dirname,ent->d_name);
      n_files++;
    }
  }
  closedir(dir);

  //Check first file's header
  int blklen1;
  GadgetHeader header;
  FILE *fp=my_fopen(fname_array[0],"r");
  my_fread(&blklen1,sizeof(int),1,fp);
  if(blklen1!=sizeof(GadgetHeader)) bad_block();
  my_fread(&header,sizeof(GadgetHeader),1,fp);
  my_fread(&blklen1,sizeof(int),1,fp);
  if(blklen1!=sizeof(GadgetHeader)) bad_block();
  fclose(fp);

  if(header.num_files!=n_files) {
    fprintf(stderr,"Wrong number of files! %d != %d\n",
	    header.num_files,n_files);
  }

  for(i=0;i<n_files;i++)
    free(fname_array[i]);
  free(fname_array);

  return header;
}

static void read_snapshot(Particles *particles)
{
  int i;
  struct dirent *ent;
  DIR *dir;
  char **fname_array;
  int n_files=Param.n_files;
  int prefix_len=strlen(Param.init_prefix);
  int this_node=Param.i_node;

  //Collect filenames
  fname_array=(char **)malloc(n_files*sizeof(char *));
  for(i=0;i<n_files;i++)
    fname_array[i]=(char *)malloc(256*sizeof(char));

  dir=opendir(Param.init_dir);
  int i_file=0;
  while((ent=readdir(dir))!=NULL) {
    if(!strncmp(ent->d_name,Param.init_prefix,prefix_len)) {
      sprintf(fname_array[i_file],"%s/%s",Param.init_dir,ent->d_name);
      i_file++;
    }
  }
  closedir(dir);

  //Check first file's header
  int blklen1,blklen2;
  GadgetHeader header;
  FILE *fp=my_fopen(fname_array[0],"r");
  my_fread(&blklen1,sizeof(int),1,fp);
  if(blklen1!=sizeof(GadgetHeader)) bad_block();
  my_fread(&header,sizeof(GadgetHeader),1,fp);
  my_fread(&blklen1,sizeof(int),1,fp);
  if(blklen1!=sizeof(GadgetHeader)) bad_block();
  fclose(fp);
  if(header.num_files!=Param.n_files) {
    msg_abort(123,"Wrong number of files! %d != %d\n",
	      header.num_files,n_files);
  }
  
  //Compute total number of particles
  unsigned long long np_tot=(unsigned long long)(header.np_total[1])+
    ((unsigned long long)(header.np_total_highword[1]) << 32);
  Particle *p=particles->p;
  
  //Iterate through files
  //Read all particles within domain
  int np_saved=0;
  long long np_read=0;
  float lbox_f=(float)(Param.boxsize);
  float dx_domain=Param.boxsize/Param.n_nodes;
  float edge_total_left=Param.x_offset;
  float edge_domain_right=edge_total_left+dx_domain;

  for(i=0;i<n_files;i++) {
    int np_here;
    int j;
    int file_index=(Param.i_node+i)%n_files;

#ifdef _DEBUG
    printf("%d reading %d-th file\n",Param.i_node,file_index);
#endif //_DEBUG

    fp=my_fopen(fname_array[file_index],"r");
    my_fread(&blklen1,sizeof(int),1,fp);
    if(blklen1!=sizeof(GadgetHeader)) bad_block();
    my_fread(&header,sizeof(GadgetHeader),1,fp);
    my_fread(&blklen1,sizeof(int),1,fp);
    if(blklen1!=sizeof(GadgetHeader)) bad_block();
    
    np_here=header.np[1];
    if(np_read+np_here>np_tot) {
      fprintf(stderr,"There's something wrong with the number of particles!\n");
      exit(1);
    }
    else {
      np_read+=np_here;
    }

    //Read positions
    my_fread(&blklen1,sizeof(int),1,fp);
    int np_got_here=0;
    for(j=0;j<np_here;j++) {
      int ax;
      float x[3];

      my_fread(x,sizeof(float),3,fp);
      for(ax=0;ax<3;ax++)
	x[ax]=wrap_float(x[ax],lbox_f);

      if(x[0]>=edge_total_left && x[0]<edge_domain_right) {
	if(np_saved>=particles->np_allocated)
	  msg_abort(123,"Too many particles in file, enlarge alloc_factor %d %d %d\n",
		    np_saved,particles->np_allocated,particles->np_total);

	x[0]-=edge_total_left;
	for(ax=0;ax<3;ax++)
	  p[np_saved].x[ax]=x[ax];
	p[np_saved].id=j;
	p[np_saved].fof_id=0;
	p[np_saved].cll_id=-1;

	np_saved++;
	np_got_here++;
      }
    }
    my_fread(&blklen2,sizeof(int),1,fp);
    if(blklen1!=blklen2) bad_block();

    if(np_got_here==0) continue;

    //Read velocities
    my_fread(&blklen1,sizeof(int),1,fp);
    int np_new=0;
    for(j=0;j<np_here;j++) {
      float v[3];
      int new_p_index=np_saved-np_got_here+np_new;
      my_fread(v,sizeof(float),3,fp);
      if(p[new_p_index].id==j) {
	int ax;
	for(ax=0;ax<3;ax++)
	  p[new_p_index].v[ax]=v[ax];
	
	np_new++;
      }
    }
    if(np_new!=np_got_here)
      msg_abort(1001,"Error reading file\n");
    my_fread(&blklen2,sizeof(int),1,fp);
    if(blklen1!=blklen2) bad_block();
    
    //Read ids
    my_fread(&blklen1,sizeof(int),1,fp);
    np_new=0;
    for(j=0;j<np_here;j++) {
      int new_p_index=np_saved-np_got_here+np_new;
#ifdef _LONGIDS
      unsigned long long id;
      my_fread(&id,sizeof(unsigned long long),1,fp);
#else //_LONGIDS
      unsigned int id;
      my_fread(&id,sizeof(unsigned int),1,fp);
#endif //_LONGIDS
      if(p[new_p_index].id==j) {
	p[new_p_index].id=id;
	np_new++;
      }
    }
    if(np_new!=np_got_here)
      msg_abort(1001,"Error reading file\n");
    my_fread(&blklen2,sizeof(int),1,fp);
    if(blklen1!=blklen2) bad_block();

    fclose(fp);
  }
  if(np_read!=np_tot) {
    fprintf(stderr,"There's something wrong with the number of particles!\n");
    exit(1);
  }
#ifdef _DEBUG
  printf(" Node %d read %d particles\n",this_node,np_saved);
#endif //_DEBUG
  for(i=0;i<n_files;i++)
    free(fname_array[i]);
  free(fname_array);
  msg_printf(" Done reading\n");

  //Check total number of particles
  unsigned long long np_saved_here=np_saved;
  unsigned long long np_global=0;
  MPI_Reduce(&np_saved_here,&np_global,1,MPI_UNSIGNED_LONG_LONG,
	     MPI_SUM,0,MPI_COMM_WORLD);
  if(this_node == 0 && np_global!=np_tot)
    msg_abort(123,"Not all particles were saved %llu != %llu\n",np_global,np_tot);

  particles->np_indomain=np_saved;

  msg_printf(" Sorting particles and sharing buffer with neighbor nodes\n");
  //Sort particles by x[0]
  qsort(particles->p,particles->np_indomain,sizeof(Particle),compare_parts);

  //Find particles that will be sent to left
  int np_toleft=0;
  int found_last=0;
  float np_toleft_expected=particles->np_indomain*(Param.dx_extra/dx_domain);
  while(found_last==0) {
    float x=particles->p[np_toleft].x[0];
    if(x>Param.dx_extra)
      found_last=1;
    else {
      np_toleft++;
      if(np_toleft>=particles->np_indomain)
	msg_abort(123,"There seem to be too many particles on the "
		  "left side of the domain for node %d\n",Param.i_node);
    }
  }
  if(np_toleft>3*np_toleft_expected || np_toleft<0.3*np_toleft_expected)
    msg_printf("WARNING: there seem to be too many particles on the "
	       "left side of the domain for node %d. %d found, about %.1f expected\n",
	       Param.i_node,np_toleft,np_toleft_expected);
  particles->np_toleft=np_toleft;
  particles->np_centre=particles->np_indomain-np_toleft;

  //Allocate the buffer particles that will be sent back.
  particles->p_back=malloc(particles->np_toleft*sizeof(Particle));
  if(particles->p_back==NULL)
    msg_abort(123,"Failed to allocate memory for buffer particles\n");

  //Send particles to left and receive from right
  int tag=100;
  int np_fromright;
  MPI_Status stat;
  MPI_Sendrecv(&np_toleft,   1,MPI_INT,Param.i_node_left ,tag,
	       &np_fromright,1,MPI_INT,Param.i_node_right,tag,
	       MPI_COMM_WORLD,&stat);
#ifdef _DEBUG
  printf(" Node %d will send %d particles to node %d and will receive %d from node %d\n",
	 Param.i_node,np_toleft,Param.i_node_left,np_fromright,Param.i_node_right);
#endif //_DEBUG

  if(particles->np_indomain+np_fromright>particles->np_allocated)
    msg_abort(123,"Not enough space for received particles in node %d\n",Param.i_node);
  particles->np_fromright=np_fromright;
  particles->np_local=particles->np_indomain+np_fromright;

  MPI_Sendrecv(particles->p,particles->np_toleft,
	       ParticleMPI,Param.i_node_left,tag,
	       &(particles->p[particles->np_indomain]),particles->np_fromright,
	       ParticleMPI,Param.i_node_right,tag,
	       MPI_COMM_WORLD,&stat);

  //Add offset to received particles
  int i_off=particles->np_indomain;
  for(i=0;i<particles->np_fromright;i++)
    particles->p[i_off+i].x[0]+=dx_domain;
  msg_printf(" Done\n\n");
}

Particles *read_input_snapshot(void)
{
  Particles *particles=allocate_particles();
  read_snapshot(particles);

  return particles;
}

static void write_halos_pernode(int n_halos,FoFHalo *fhal)
{
  int ii;
  char fname[256];
  FILE *fo;
  
  //Sort by mass
  qsort(fhal,n_halos,sizeof(FoFHalo),compare_halo);

  sprintf(fname,"%s.%d.dat",Param.output_prefix,Param.i_node);
  msg_printf(" Writing to file %s\n",fname);
  fo=my_fopen(fname,"w");
  fprintf(fo,"#ID NP Mass x_avg[3] x_rms[3] v_avg[3] v_rms[3]\n");
  for(ii=0;ii<n_halos;ii++) {
    fprintf(fo,"%d %d %lE ",ii,fhal[ii].np,fhal[ii].m_halo);
    fprintf(fo,"%lE %lE %lE ",
	    fhal[ii].x_avg[0],fhal[ii].x_avg[1],fhal[ii].x_avg[2]);
    fprintf(fo,"%lE %lE %lE ",
	    fhal[ii].x_rms[0],fhal[ii].x_rms[1],fhal[ii].x_rms[2]);
    fprintf(fo,"%lE %lE %lE ",
	    fhal[ii].v_avg[0],fhal[ii].v_avg[1],fhal[ii].v_avg[2]);
    fprintf(fo,"%lE %lE %lE ",
	    fhal[ii].v_rms[0],fhal[ii].v_rms[1],fhal[ii].v_rms[2]);
    fprintf(fo,"\n");
  }
  fclose(fo);

  free(fhal);
}

static void write_halos_root(int n_halos,FoFHalo *fhal)
{
  //Gather halos
  int ii;
  int n_halos_total=0;
  int n_halos_here=n_halos;
  int *n_halos_array=NULL;
  int *n_disp_array=NULL;
  FoFHalo *fh_tot=NULL;

  if(Param.i_node==0) {
    n_halos_array=malloc(Param.n_nodes*sizeof(int));
    if(n_halos_array==NULL)
      msg_abort(123,"Out of memory\n");
    n_disp_array=malloc(Param.n_nodes*sizeof(int));
    if(n_disp_array==NULL)
      msg_abort(123,"Out of memory\n");
  }

  MPI_Gather(&n_halos_here,1,MPI_INT,
	     n_halos_array,1,MPI_INT,
	     0,MPI_COMM_WORLD);

  if(Param.i_node==0) {
    for(ii=0;ii<Param.n_nodes;ii++)
      n_halos_total+=n_halos_array[ii];
    if(n_halos_total<=0)
      msg_abort(123,"Something went wrong\n");
    msg_printf(" %d halos found in total with %d particles of more\n",
	       n_halos_total,Param.np_min);

    //Calculate array of displacements
    n_disp_array[0]=0;
    for(ii=1;ii<Param.n_nodes;ii++)
      n_disp_array[ii]=n_halos_array[ii-1];

    fh_tot=malloc(n_halos_total*sizeof(FoFHalo));
    if(fh_tot==NULL)
      msg_abort(123,"Failed to allocate memory for all halos\n");
  }
  MPI_Gatherv(fhal,n_halos_here,HaloMPI,
	      fh_tot,n_halos_array,n_disp_array,HaloMPI,
	      0,MPI_COMM_WORLD);
  free(fhal);

  if(Param.i_node==0) {
    char fname[256];
    FILE *fo;

    free(n_halos_array);
    free(n_disp_array);

    //Sort by mass
    qsort(fh_tot,n_halos_total,sizeof(FoFHalo),compare_halo);

    sprintf(fname,"%s.dat",Param.output_prefix);
    msg_printf(" Writing to file %s\n",fname);
    fo=my_fopen(fname,"w");
    fprintf(fo,"#ID NP Mass x_avg[3] x_rms[3] v_avg[3] v_rms[3]\n");
    for(ii=0;ii<n_halos_total;ii++) {
      fprintf(fo,"%d %d %lE ",ii,fh_tot[ii].np,fh_tot[ii].m_halo);
      fprintf(fo,"%lE %lE %lE ",
	      fh_tot[ii].x_avg[0],fh_tot[ii].x_avg[1],fh_tot[ii].x_avg[2]);
      fprintf(fo,"%lE %lE %lE ",
	      fh_tot[ii].x_rms[0],fh_tot[ii].x_rms[1],fh_tot[ii].x_rms[2]);
      fprintf(fo,"%lE %lE %lE ",
	      fh_tot[ii].v_avg[0],fh_tot[ii].v_avg[1],fh_tot[ii].v_avg[2]);
      fprintf(fo,"%lE %lE %lE ",
	      fh_tot[ii].v_rms[0],fh_tot[ii].v_rms[1],fh_tot[ii].v_rms[2]);
      fprintf(fo,"\n");
    }
    fclose(fo);

    free(fh_tot);
  }

  msg_printf(" Done\n\n");
}

void write_halos(int n_halos,FoFHalo *fhal)
{
  if(Param.output_pernode)
    write_halos_pernode(n_halos,fhal);
  else
    write_halos_root(n_halos,fhal);
}

#include "allvars.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#define DL for (int d=0;d<3;d++) /* Dimension loop */
#define BF 1e30 /* Big factor */
#define MAXNAMELENGTH 10
/* #define BGC_MAGIC ((uint64_t)0x1234567801010101ll) */
#define BGC2_HEADER_SIZE 1024

int read_unfbi77(char *fname, float ***p, int decfact, int ThisTask);
int read_gadget(char *fname0, float ***p, int decfact, int p_type, int ThisTask);
int read_bgtreebi(char *fname, float ***p, int decfact, int nskip, int ThisTask);
int read_bgc2(char *fname0, float ***p, int decfact, int haloid, int ThisTask);
int read_tipsy(char *fname, float ***p, int decfact, int p_type, int ThisTask);

void namefile(char *filename, char *name, char *suffix, char *ext) {
  char dumstr[100];
  sprintf(filename,"%s%s/%s",All.OutputDir,name,All.FilePrefix);
  /* Decimation factor */
  if (All.DecimateInputBy > 1) {
    sprintf(dumstr,"_dec%d",All.DecimateInputBy);
    strcat(filename,dumstr);
  }
  /* Particle type */
  if (All.ParticleType >= 0) {
    if ((All.PositionFileFormat == 1) || (All.PositionFileFormat == 4)) {
      sprintf(dumstr,"_%d",All.ParticleType);
      strcat(filename,dumstr);
    }
  }
  /* Particle type */
  if ((All.PositionFileFormat == 3) && (All.Bgc2HaloId >= 0)) {
    sprintf(dumstr,"_%d",All.Bgc2HaloId);
    strcat(filename,dumstr);
  }
  /* Squish factors */
  if ((All.SquishY > 0) && (All.SquishY < 1)) {
    /* sprintf(dumstr,"_sqY%d",(int)(1.0/All.SquishY)); */
    sprintf(dumstr,"_sqY%4.2f",All.SquishY);
    dumstr[5] = 'p';
    strcat(filename,dumstr);
  }
  if ((All.SquishZ > 0) && (All.SquishZ < 1)) {
    /* sprintf(dumstr,"_sqZ%d",(int)(1.0/All.SquishZ)); */
    sprintf(dumstr,"_sqZ%4.2f",All.SquishZ);
    dumstr[5] = 'p';
    strcat(filename,dumstr);
  }
  /* File extension */
  sprintf(dumstr,"%s.%s%s",suffix,name,ext);
  strcat(filename,dumstr);
}

void namefile_blk(char *filename, char *suffix, int xsplit, int ysplit, int zsplit) {
  char ext[MAXNAMELENGTH];
  sprintf(ext,".%02d.%02d.%02d",xsplit,ysplit,zsplit);
  namefile(filename,"part",suffix,ext);
}

void namefile_adj(char *filename, char *suffix) { namefile(filename,"adjs",suffix,""); }

void namefile_vol(char *filename, char *suffix) { namefile(filename,"vols",suffix,""); }

int read_positions_file(char *fname, float ***p, int ThisTask) {
  int i,np = 0;
  float rmin[3],rmax[3],rmin_tot,rmax_tot;

  /* Read position file based on file format on process 0 */
  printf("Thread %03d: -----------------------------\n",ThisTask);
  printf("Thread %03d: Reading in position file...\n",ThisTask);
  if (All.PositionFileFormat == 0) {
    np = read_unfbi77(fname,p,All.DecimateInputBy,ThisTask);
  } else if (All.PositionFileFormat == 1) {
    np = read_gadget(fname,p,All.DecimateInputBy,All.ParticleType,ThisTask);
  } else if (All.PositionFileFormat == 2) {
    np = read_bgtreebi(fname,p,All.DecimateInputBy,All.BgTreebiNskip,ThisTask);
  } else if (All.PositionFileFormat == 3) {
    np = read_bgc2(fname,p,All.DecimateInputBy,All.Bgc2HaloId,ThisTask);
  } else if (All.PositionFileFormat == 4) {
    np = read_tipsy(fname,p,All.DecimateInputBy,All.ParticleType,ThisTask);
  } else {
    printf("Thread %03d: Invalid PositionFileFormat = %d.\n",ThisTask,All.PositionFileFormat);
    np = -2;
  }
  if (np > 0) {
    printf("Thread %03d: %d particles read from file\n",ThisTask,np);
    /* for (i=0; i<3; i++) { */
    /*   printf("Thread %03d: test1 %d %+4.2f %+4.2f %+4.2f\n",ThisTask,i,(*p)[i][0],(*p)[i][1],(*p)[i][2]); */
    /*   printf("Thread %03d: test2 %d %+4.2f %+4.2f %+4.2f\n",ThisTask,i,*(*((*p)+i)+0),*(*((*p)+i)+1),*(*((*p)+i)+2)); */
    /*   printf("\n"); */
    /* } */
    printf("Thread %03d: -----------------------------\n\n",ThisTask);
    fflush(stdout);
    
    /* Get minimum & maximum */
    printf("Thread %03d: -----------------------------\n",ThisTask);
    printf("Thread %03d: Scaling positions to unit box...\n",ThisTask);
    DL { rmin[d] = BF; rmax[d] = -BF; }
    for (i=0; i<np; i++) {
      DL {
	if ((*p)[i][d] < rmin[d]) rmin[d] = (*p)[i][d];
	if ((*p)[i][d] > rmax[d]) rmax[d] = (*p)[i][d];
      }
    }
    printf("Thread %03d: Raw positions:\n",ThisTask);
    printf("Thread %03d:     xmin=%+4.2f ymin=%+4.2f zmin=%+4.2f\n",
	   ThisTask,rmin[0],rmin[1],rmin[2]);
    printf("Thread %03d:     xmax=%+4.2f ymax=%+4.2f zmax=%+4.2f\n",
	   ThisTask,rmax[0],rmax[1],rmax[2]);

    /* Squish in some direction */
    if ((All.SquishY != 1) || (All.SquishZ != 1)) {
      DL { rmin[d] = BF; rmax[d] = -BF; }
      for (i=0; i<np; i++) {
	printf("SquishY: %f\n",All.SquishY);
	printf("SquishZ: %f\n",All.SquishZ);
	(*p)[i][1] *= All.SquishY;
	(*p)[i][2] *= All.SquishZ;
	DL {
	  (*p)[i][d] *= pow(All.SquishY*All.SquishZ,-1./3.);
	  if ((*p)[i][d] < rmin[d]) rmin[d] = (*p)[i][d];
	  if ((*p)[i][d] > rmax[d]) rmax[d] = (*p)[i][d];
	}
      }
      printf("Thread %03d: Squished positions:\n",ThisTask);
      printf("Thread %03d:     xmin=%+4.2f ymin=%+4.2f zmin=%+4.2f\n",
	     ThisTask,rmin[0],rmin[1],rmin[2]);
      printf("Thread %03d:     xmax=%+4.2f ymax=%+4.2f zmax=%+4.2f\n",
	     ThisTask,rmax[0],rmax[1],rmax[2]);
    }
    
    /* Total min/max */
    rmin_tot = BF;
    rmax_tot = -BF;
    DL {
      if (rmin[d] < rmin_tot) rmin_tot = rmin[d];
      if (rmax[d] > rmax_tot) rmax_tot = rmax[d];
    }
    if (!All.PeriodicBoundariesOn) {
      All.BoxSize = rmax_tot-rmin_tot;
      rmin_tot -= All.Border*All.BoxSize;
      rmax_tot += All.Border*All.BoxSize;
      All.BoxSize*=(1.+2.*All.Border);
    }

    /* Scale positions to be on 0 to 1 */
    DL { rmin[d] = BF; rmax[d] = -BF; }
    for (i=0; i<np; i++) {
      DL {
	(*p)[i][d] -= rmin_tot;
	(*p)[i][d] *= 1./All.BoxSize;
	  if ((*p)[i][d] < rmin[d]) rmin[d] = (*p)[i][d];
	  if ((*p)[i][d] > rmax[d]) rmax[d] = (*p)[i][d];
      }
    }
    printf("Thread %03d: Scaled positions:\n",ThisTask);
    printf("Thread %03d:     xmin=%+4.2f ymin=%+4.2f zmin=%+4.2f\n",
	   ThisTask,rmin[0],rmin[1],rmin[2]);
    printf("Thread %03d:     xmax=%+4.2f ymax=%+4.2f zmax=%+4.2f\n",
	   ThisTask,rmax[0],rmax[1],rmax[2]);
    printf("Thread %03d: -----------------------------\n\n",ThisTask);
    fflush(stdout);
  }

  return(np);
}

int read_unfbi77(char *fname, float ***p, int decfact, int ThisTask)
{
  FILE *fd;
  int np,npdec,pc,dum,i;
  float *ptemp;
  /* Fortran77 4-byte headers and footers */
  /* Delete "dum" statements if you don't need them */

  /* Open file */
  if (!(fd=fopen(fname,"r"))) {
    printf("Thread %03d: can't open file `%s`\n",ThisTask,fname);
    return(-10);
  }

  /* Read number of particles */
  fread(&dum,1,4,fd); 
  fread(&np,1, sizeof(int),fd); 
  fread(&dum,1,4,fd);

  /* Decimate */
  npdec = np/decfact;
  printf("Thread %03d: Decimating %d particle by %d results in %d particles\n",
	 ThisTask,np,decfact,npdec);
  if (npdec == 0) {
    fclose(fd);
    return(-50);
  }

  /* Allocate */
  ptemp = (float *)malloc(np*sizeof(float));
  if (ptemp == NULL) {
    printf("Thread %03d: Unable to allocate temporary particle arry in read_unfbi77!\n",
	   ThisTask);
    fclose(fd);
    return(-20);
  }
  (*p) = (float **)malloc(npdec*sizeof(float *));
  if ((*p) == NULL) {
    printf("Thread %03d: Unable to allocate particle array in read_unfbi77!\n",
	   ThisTask);
    fclose(fd);
    free(ptemp);
    return(-20);
  }

  /* Read masses */
  fread(&dum,1,4,fd);
  fread(ptemp,np,4,fd);
  fread(&dum,1,4,fd);

  /* Read x */
  fread(&dum,1,4,fd);
  fread(ptemp,np,4,fd);
  for (i=0, pc=0; i<np; i++) {
    if ((i % decfact) == 0) {
      (*p)[pc] = (float *)malloc(3*sizeof(float));
      if ((*p)[pc] == NULL) {
	printf("Thread %03d: Unable to allocate particle array in read_unfbi77!\n",ThisTask);
	fclose(fd); free(ptemp); free((*p));
	return(-20);
      }
      (*p)[pc][0] = ptemp[i];
      pc++;
    }
  }
  fread(&dum,1,4,fd);

  /* Read y */
  fread(&dum,1,4,fd);
  fread(ptemp,np,4,fd);
  for (i=0, pc=0; i<np; i++) {
    if ((i % decfact) == 0) {
      (*p)[pc][1] = ptemp[i];
      pc++;
    }
  }
  fread(&dum,1,4,fd);

  /* Read z */
  fread(&dum,1,4,fd);
  fread(ptemp,np,4,fd);
  for (i=0, pc=0; i<np; i++) {
    if ((i % decfact) == 0) {
      (*p)[pc][2] = ptemp[i];
      pc++;
    }
  }
  fread(&dum,1,4,fd);

  /* Close, free, and return */
  fclose(fd);
  free(ptemp);
  return(npdec);
}

int read_bgtreebi(char *fname, float ***p, int decfact, int nskip, int ThisTask)
{
  int ns,np,i,pc,npdec;
  FILE *fd;
  int dim,idum1,idum2;
  float fdum1,fdum2,fdum3;

  /* Prevent nskip from being negative */
  if (nskip < 0) {
    ns = 0;
  } else {
    ns = nskip;
  }
  
  /* Open file */
  if (!(fd=fopen(fname,"r"))) {
    printf("Thread %03d: can't open file `%s`\n",ThisTask,fname);
    return(-10);
  }

  /* Read header */
  fscanf(fd, " %d     %d     %d", &np, &idum1, &idum2);
  fscanf(fd, "      %d", &dim);
  fscanf(fd, "  %g",&fdum1);
  printf("Thread %03d: %d particles are in this file\n",ThisTask,np);
  np-= ns;
  printf("Thread %03d: %d are halo particles\n",ThisTask,np);

  /* Decimate */
  npdec = np/decfact;
  printf("Thread %03d: Decimating %d particle by %d results in %d particles\n",
	 ThisTask,np,decfact,npdec);
  if (npdec == 0) {
    fclose(fd);
    return(-50);
  }

  /* Read masses */
  for (i=0; i<(np+ns); i++)
    fscanf(fd, "  %g",&fdum1);

  /* Allocate */
  (*p) = (float **)malloc(npdec*sizeof(float *));
  if ((*p) == NULL) {
    printf("Thread %03d: Unable to allocate particle array in read_bgtreebi!\n",
	   ThisTask);
    fclose(fd);
    return(-20);
  }
  (**p)--; 

  /* Read positions */
  for (i=0, pc=0; i<(np+ns); i++) {
    if ((i < ns) || (((i-ns) % decfact)!=0)) {
      fscanf(fd, "  %g %g %g",&fdum1,&fdum2,&fdum3);
    } else {
      (*p)[pc] = (float *)malloc(dim*sizeof(float));
      if ((*p)[pc] == NULL) {
	printf("Thread %03d: Unable to allocate %d particle in read_bgtreebi!\n",
	       ThisTask,pc);
	fclose(fd);
	return(-21);
      }
      fscanf(fd, "  %g %g %g", &(*p)[pc][0], &(*p)[pc][1], &(*p)[pc][2]);
      pc++;
    }
  }
  fclose(fd);
  return(npdec);
}


/* this routine loads particle positions from Rockstar's Bgc2 files */
int read_bgc2(char *fname0, float ***p, int decfact, int haloid, int ThisTask)
{
  FILE   *fd;
  char   fname[MAXLEN_FILENAME];
  char   buf[200];
  char   *extstr;
  int    i,k,files;
  uint32_t dummy;
  int    n,ntot,pc_new;
  int    np,npdec;
  int    ng,gtot;
  int    ngroups;

  struct bgc2_header {
    uint64_t magic;             /* A magic number to identify this as a BGC file. */
    int64_t version;            /* File version number. */

    int64_t num_files;          /* number of files output is distributed into */
    int64_t file_id;            /* this files ID (number) if multiple files are output */
    int64_t snapshot;           /* Snapshot ID */

    int64_t format_group_data;  /* output group data format identifier, see enum gdata_format below */
    int64_t format_part_data;   /* output particle data format identifier, see enum pdata_format below */

    int64_t group_type;         /* FOF, SO, etc */
    int64_t ngroups;            /* number of groups stored LOCALLY in this file */
    int64_t ngroups_total;      /* number of groups stored GLOBALLY over all output BGC files */

    int64_t npart;              /* number of particles bound in groups LOCALLY in this file */
    int64_t npart_total;        /* number of particles bound in groups GLOBALLY over all output BGC files */
    int64_t npart_orig;         /* number of particles from original simulation input */

    int64_t max_npart;          /* maximum number of particles in one group LOCALLY in this file */
    int64_t max_npart_total;    /* maximum number of particles in one group GLOBALLY over all output BGC files */

    int64_t min_group_part;     /* minimum number of particles in a group */
    int64_t valid_part_ids;     /* valid particle IDs mean they match input snapshot */

    double linkinglength;       /* for FOF halos, what linking length is used */
    double overdensity;         /* mostly SO: overdensity with respect to mean */

    double time;                /* time of the input snapshot */
    double redshift;            /* redshift of the input snapshot */

    double box_size;            /* input BoxSize */
    double box_min[3];          /* alternative to center, Gadget assumes (0,0,0) */
    double bounds[6];           /* Spatial bounds of the halo centers contained in this file */

    double part_mass;           /* mass of particles if only one */

    double Omega0;              /* Matter density at z=0 */
    double OmegaLambda;         /* Dark energy density at z=0 */

    /* these define the units, but are NOT ALWAYS SET */
    double Hubble0;             /* in whatever units are used */
    double GravConst;           /* in whatever units are used */

    uint8_t padding[BGC2_HEADER_SIZE - ( 36 * 8 )];
  } header1;

  typedef struct {
    int64_t id;
    int64_t parent_id;
    uint64_t npart;
    uint64_t npart_self; /* Excluding substructure */
    float radius;
    float mass;
    float pos[3];
    float vel[3];
    float vmax, rvmax;
  } BGC_GROUP_DATA;

  typedef struct {
    int64_t part_id;
    float pos[3];
    float vel[3];
  } BGC_PART_DATA;
  
  BGC_GROUP_DATA dummy_group;
  BGC_PART_DATA dummy_part;
  BGC_GROUP_DATA **groups;
  
#define SKIP fread(&dummy, sizeof(dummy), 1, fd);
#define SKIP_GROUP fread(&dummy_group, sizeof(dummy_group), 1, fd); 
#define SKIP_PART fread(&dummy_part, sizeof(dummy_part), 1, fd);
  
  /* Remove trailing split if present */
  extstr = strrchr(fname0,'.')-5;
  dummy = sscanf(extstr,".%04d.bgc2",&i);
  if (dummy == 1) {
    k = extstr-fname0;
    strncpy(fname,fname0,k);
    fname[k] = '\0';
  } else {
    strcpy(fname,fname0);
  }

  /* Read header from first file to get number of files */
  sprintf(buf,"%s",fname);
  if(!(fd=fopen(buf,"rb")))
    {
      sprintf(buf,"%s.%04d.bgc2",fname,0);
      if(!(fd=fopen(buf,"r")))
	{
	  printf("Thread %03d: can't open file `%s`\n",ThisTask,fname);
	  return(-10);
	}
    }
  SKIP;
  fread(&header1, sizeof(header1), 1, fd);
  SKIP;
  if (header1.format_group_data!=40) {
    printf("Thread %03d: group    format = %d\n",ThisTask,(int)(header1.format_group_data));
    return(-10);
  }
  if (header1.format_part_data!=30) {
    printf("Thread %03d: particle format = %d\n",ThisTask,(int)(header1.format_part_data));
    return(-10);
  }
  files = (int)(header1.num_files);
  ngroups = (int)(header1.ngroups_total);
  fclose(fd);
  printf("Thread %03d: %d groups split between %d files.\n",ThisTask,ngroups,files);

  /* Allocate groups */
  groups = (BGC_GROUP_DATA **)malloc(ngroups*sizeof(BGC_GROUP_DATA *));
  if (groups == NULL) {
    printf("Thread %03d: Unable to allocate groups array in read_bgc2!\n",ThisTask);
    return(-20);
  }

  /* Read groups from each file, counting particles */
  for(i=0,ng=0,gtot=0,np=0;i<files;i++) {

    /* File name */
    if (files > 1) {
      sprintf(buf,"%s.%04d.bgc2",fname,i);
    } else {
      sprintf(buf,"%s",fname);
    }
    if(!(fd=fopen(buf,"r")))
      {
	printf("Thread %03d: can't open file `%s`\n",ThisTask,buf);
	return(-10);
      }

    /* Header */
    SKIP;
    fread(&header1, sizeof(header1), 1, fd);
    SKIP;
    /* printf("Thread %03d: Reading %d groups from file %s\n",ThisTask,(int)header1.ngroups,buf); */

    /* Group info */
    SKIP;
    for(n=0;n<header1.ngroups;n++,gtot++) {
      groups[gtot] = (BGC_GROUP_DATA *)malloc(sizeof(BGC_GROUP_DATA));
      if (groups[gtot] == NULL) {
	printf("Thread %03d: Unable to allocate group %d in read_bgc2!\n",ThisTask,gtot);
	return(-20);
      }
      /*SKIP;*/
      fread(groups[gtot],sizeof(BGC_GROUP_DATA),1,fd);
      /* SKIP; */
      /* Negative particle number means something is wrong */
      if ((int)groups[gtot]->npart < 0) {
	printf("Thread %03d: Negative number of particles in group %d.\n",ThisTask,gtot);
	return(-25);
      }
      /* Set haloid to first group if not provided */
      if ((haloid < 0) && (gtot==0)) haloid = groups[gtot]->id;
      if (groups[gtot]->id == haloid) { /* Use parent_id if substructure not included in parent */
	np+=groups[gtot]->npart;
	ng++;
      }
    }
    fclose(fd);
  }
  printf("Thread %03d: Found %d particles across %d groups in halo %d.\n",ThisTask,np,ng,haloid);
  if (gtot != ngroups) {
    printf("Thread %03d: Found %d groups, but %d were expected.\n",ThisTask,gtot,ngroups);
    return(-30);
  }
  
  /* Decimate */
  npdec = np/decfact;
  if ((np % decfact) != 0) npdec++;
  printf("Thread %03d: Decimating %d particle by %d results in %d particles\n",
	 ThisTask,np,decfact,npdec);
  if (npdec < 4) { 
    printf("Thread %03d: Cannot create a volume with <4 particles\n",ThisTask);
    return(-50);
  }

  /* Allocate */
  (*p) = (float **)malloc(npdec*sizeof(float *));
  if ((*p) == NULL) {
    printf("Thread %03d: Unable to allocate particle array in read_bgc2!\n",
	   ThisTask);
    return(-20);
  }
  (**p)--; 

  /* Read particle data */
  for(i=0,gtot=0,ntot=0,pc_new=0;i<files;i++) {
    /* File name */
    if (files > 1) {
      sprintf(buf,"%s.%04d.bgc2",fname,i);
    } else {
      sprintf(buf,"%s",fname);
    }
    if(!(fd=fopen(buf,"r"))) {
      printf("Thread %03d: can't open file `%s`\n",ThisTask,buf);
      return(-10);
    }

    /* Header */
    SKIP;
    fread(&header1, sizeof(header1), 1, fd);
    SKIP;

    /* Group info */
    SKIP;
    for(k=0;k<header1.ngroups;k++) {
      SKIP_GROUP;
    }
    SKIP;

    /* Particle info */
    for(k=0;k<header1.ngroups;k++,gtot++) {
      SKIP;
      if (groups[gtot]->id == haloid) {
	for(n=0;n<(int)(groups[gtot]->npart);n++,ntot++) {
	  if ((ntot % decfact) == 0) {
	    (*p)[pc_new] = (float *)malloc(3*sizeof(float));
	    if ((*p)[pc_new] == NULL) {
	      printf("Thread %03d: Unable to allocate %d particle in read_gadget!\n",
		     ThisTask,pc_new);
	      return(-21);
	    }
	    SKIP_PART;
	    memmove((*p)[pc_new], dummy_part.pos, 3*sizeof(float));
	    pc_new++;
	  } else {
	    SKIP_PART;
	  }
	}
      } else {
	for(n=0;n<(int)(groups[gtot]->npart);n++) {
	  SKIP_PART;
	}
      }
      SKIP;
    }

    /* Close this file */
    fclose(fd);
  }

  /* Check that decimate was handled correctly */
  if (pc_new != npdec) {
    printf("Thread %03d: Read in %d particles, but %d were expected from the header.\n",
	   ThisTask,pc_new,npdec);
    return(-51);
  }

  return(npdec);
}

/* this routine loads particle positions from Gadget's default
 * binary file format. (A snapshot may be distributed into
 * multiple files so long as the header reflects this)
 */
int read_gadget(char *fname0, float ***p, int decfact, int p_type, int ThisTask)
{
  char   fname[MAXLEN_FILENAME];
  FILE   *fd;
  char   buf[200];
  char   *extstr;
  int    i,k,dummy,maxtype;
  float  dummy_pos[3];
  int    n,ntot,pc,pc_new;
  int    files,np,npall,npdec;

  struct io_header_1
  {
    int      npart[6];
    double   mass[6];
    double   time;
    double   redshift;
    int      flag_sfr;
    int      flag_feedback;
    int      npartTotal[6];
    int      flag_cooling;
    int      num_files;
    double   BoxSize;
    double   Omega0;
    double   OmegaLambda;
    double   HubbleParam; 
    char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
  } header1;

#define SKIP fread(&dummy, sizeof(dummy), 1, fd);
#define SKIP_POS fread(&dummy_pos, sizeof(dummy_pos), 1, fd);

  np = 0;
  npdec = 0;
  pc_new = 0;
  ntot = 0;

  /* Set maximum type number that should be read */
  if (p_type > 0)
    maxtype = p_type + 1;
  else
    maxtype = 6;

  /* Remove trailing split if present */
  extstr = strrchr(fname0,'.');
  dummy = sscanf(extstr,".%d",&i);
  if (dummy == 1) {
    k = extstr-fname0;
    strncpy(fname,fname0,k);
    fname[k] = '\0';
  } else {
    strcpy(fname,fname0);
  }

  /* Read header from first file to get number of files */
  sprintf(buf,"%s",fname);
  if(!(fd=fopen(buf,"r")))
    {
      sprintf(buf,"%s.%d",fname,0);
      if(!(fd=fopen(buf,"r")))
	{
	  printf("Thread %03d: can't open file `%s`\n",ThisTask,fname);
	  return(-10);
	}
    }
  fread(&dummy, sizeof(dummy), 1, fd);
  fread(&header1, sizeof(header1), 1, fd);
  fread(&dummy, sizeof(dummy), 1, fd);
  files = header1.num_files;
  fclose(fd);

  /* Loop over files */
  for(i=0, pc=0; i<files; i++, pc=pc_new)
    {
      /* Create file name */
      if(files>1)
	sprintf(buf,"%s.%d",fname,i);
      else
	sprintf(buf,"%s",fname);

      /* Open fle */
      if(!(fd=fopen(buf,"r")))
	{
	  printf("Thread %03d: can't open file `%s`\n",ThisTask,buf);
	  return(-10);
	}
      printf("Thread %03d: Reading Gadget file `%s' ...\n",ThisTask,buf);

      /* Read header */
      fread(&dummy, sizeof(dummy), 1, fd);
      fread(&header1, sizeof(header1), 1, fd);
      fread(&dummy, sizeof(dummy), 1, fd);

      /* Get numbers from header */
      if(files==1)
	{
	  for(k=0, npall=0; k<5; k++) 
	    npall+= header1.npart[k]; 
	  if (p_type > 0) {
	    np = header1.npart[p_type];
	  } else {
	    np = npall;
	  }
	}
      else
	{
	  for(k=0, npall=0; k<5; k++)
	    npall+= header1.npartTotal[k];
	  if (p_type > 0) {
	    np = header1.npartTotal[p_type];
	  } else {
	    np = npall;
	  }
	}

      /* For first file only */
      if (i==0) {
	/* Decimate */
	npdec = np/decfact;
	if ((np % decfact) != 0) {
	  npdec++;
	}
	printf("Thread %03d: Decimating %d particle by %d results in %d particles\n",
	       ThisTask,np,decfact,npdec);
	if (npdec < 4) { 
	  printf("Thread %03d: Cannot create a volume with <4 particles\n",ThisTask);
	  return(-50);
	}

	/* Allocate */
	(*p) = (float **)malloc(npdec*sizeof(float *));
	if ((*p) == NULL) {
	  printf("Thread %03d: Unable to allocate particle array in read_gadget!\n",
		 ThisTask);
	  return(-20);
	}
	(**p)--; 
      }

      /* Read positions */
      SKIP;
      for(k=0,pc_new=pc;k<maxtype;k++)
	{
	  /* Read positions if no type specified or type matches */
	  if ((p_type < 0) || (p_type == k)) {
	    for(n=0;n<header1.npart[k];n++,ntot++) {
	      if ((ntot % decfact) == 0) {
		(*p)[pc_new] = (float *)malloc(3*sizeof(float));
		if ((*p)[pc_new] == NULL) {
		  printf("Thread %03d: Unable to allocate %d particle in read_gadget!\n",
			 ThisTask,pc_new);
		  return(-21);
		}
		fread((*p)[pc_new], sizeof(float), 3, fd);
		pc_new++;
	      } else {
		SKIP_POS;
	      }
	    }
	  }
	  /* Otherwise skip over other types advancing the file */
	  else {
	    for(n=0;n<header1.npart[k];n++)
	      SKIP_POS;
	  }
	}
      fclose(fd);
    }

  if (pc_new != npdec) {
    printf("Thread %03d: Read in %d particles, but %d were expected from the header.\n",
	   ThisTask,pc_new,npdec);
    return(-51);
  }

  return(npdec);
}


int read_tipsy(char *fname, float ***p, int decfact, int p_type, int ThisTask)
{
  FILE   *fd;
  int    k,n,ntot,pc,dummy;
  int    np, npdec;

  struct io_header 
  {
    double   time;
    int      ntot;
    int      ndim;
    int      npart[3];
  } header1;

  typedef struct particle_gas
  {
    float   mass;
    float   pos[3];
    float   vel[3];
    float   rho;
    float   eps;
    float   metals;
    float   phi;
  } TIPSY_GAS_PART;

  typedef struct particle_dm
  {
    float   mass;
    float   pos[3];
    float   vel[3];
    float   eps;
    float   phi;
  } TIPSY_DM_PART;

  typedef struct particle_star
  {
    float   mass;
    float   pos[3];
    float   vel[3];
    float   metals;
    float   tform;
    float   eps;
    float   phi;
  } TIPSY_STAR_PART;

  TIPSY_GAS_PART dummy_gas;
  TIPSY_DM_PART dummy_dm;
  TIPSY_STAR_PART dummy_star;

  /* Open file */
  if(!(fd=fopen(fname,"r")))
    {
      printf("Thread %03d: can't open file `%s`\n",ThisTask,fname);
      return(-10);
    }

  /* Read header and weird trailing 4 bytes */
  fread(&header1, sizeof(header1), 1, fd);
  fread(&dummy, sizeof(dummy), 1, fd);

  /* Get numbers from header */
  if (p_type > 0) {
    np = header1.npart[p_type];
  } else {
    np = header1.ntot;
  }

  /* Decimate */
  npdec = np/decfact;
  if ((np % decfact) != 0) {
    npdec++;
  }
  printf("Thread %03d: Decimating %d particle by %d results in %d particles\n",
	 ThisTask,np,decfact,npdec);
  if (npdec < 4) { 
    printf("Thread %03d: Cannot create a volume with <4 particles\n",ThisTask);
    return(-50);
  }

  /* Allocate */
  (*p) = (float **)malloc(npdec*sizeof(float *));
  if ((*p) == NULL) {
    printf("Thread %03d: Unable to allocate particle array in read_tipsy!\n",
	   ThisTask);
    return(-20);
  }
  (**p)--; 

  /* Initialize things */
  ntot = 0;
  pc = 0;

  /* Gas particles */
  k = 0;
  if ((p_type < 0) || (p_type == k)) {
    for(n=0;n<header1.npart[k];n++,ntot++) {
      fread(&dummy_gas, sizeof(dummy_gas), 1, fd);
      if ((ntot % decfact) == 0) {
	(*p)[pc] = (float *)malloc(3*sizeof(float));
	if ((*p)[pc] == NULL) {
	  printf("Thread %03d: Unable to allocate %d particle in read_tipsy!\n",
		 ThisTask,pc);
	  return(-21);
	}
	(*p)[pc] = dummy_gas.pos;
	pc++;
      }
    }
  } else {
    for(n=0;n<header1.npart[k];n++) {
      fread(&dummy_gas, sizeof(dummy_gas), 1, fd);
    }
  }

  /* Dark matter particles */
  k++;
  if ((p_type < 0) || (p_type == k)) {
    for(n=0;n<header1.npart[k];n++,ntot++) {
      fread(&dummy_dm, sizeof(dummy_dm), 1, fd);
      if ((ntot % decfact) == 0) {
	(*p)[pc] = (float *)malloc(3*sizeof(float));
	if ((*p)[pc] == NULL) {
	  printf("Thread %03d: Unable to allocate %d particle in read_tipsy!\n",
		 ThisTask,pc);
	  return(-21);
	}
	(*p)[pc] = dummy_dm.pos;
	pc++;
      }
    }
  } else {
    for(n=0;n<header1.npart[k];n++) {
      fread(&dummy_dm, sizeof(dummy_dm), 1, fd);
    }
  }

  /* Star particles */
  k++;
  if ((p_type < 0) || (p_type == k)) {
    for(n=0;n<header1.npart[k];n++,ntot++) {
      fread(&dummy_star, sizeof(dummy_star), 1, fd);
      if ((ntot % decfact) == 0) {
	(*p)[pc] = (float *)malloc(3*sizeof(float));
	if ((*p)[pc] == NULL) {
	  printf("Thread %03d: Unable to allocate %d particle in read_tipsy!\n",
		 ThisTask,pc);
	  return(-21);
	}
	(*p)[pc] = dummy_star.pos;
	pc++;
      }
    }
  } else {
    for(n=0;n<header1.npart[k];n++) {
      fread(&dummy_star, sizeof(dummy_star), 1, fd);
    }
  }
  
  /* Close the file */
  fclose(fd);

  /* Ensure that we have the right number of particles */
  if (pc != npdec) {
    printf("Thread %03d: Read in %d particles, but %d were expected from the header.\n",
	   ThisTask,pc,npdec);
    return(-51);
  }

  return(npdec);
}


int read_parameter_file(char *fname) {
#define FLOAT 1
#define STRING 2
#define INT 3
#define MAXTAGS 300
  
  const int maxlen=1000;
  FILE *fd;
  char buf[maxlen], buf1[maxlen], buf2[maxlen], buf3[maxlen];
  int i, j, nt;
  int id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  int  errorFlag = 0;

  printf("-----------------------------------------\n");
  printf("Reading parameter file %s...\n",fname);

  /* Check that sizes are correct */
  if(sizeof(long long) != 8)
    {
      printf("Type `long long' is not 64 bit on this platform. Stopping.\n\n");
      errorFlag = 5;
    }
  if(sizeof(int) != 4)
    {
      printf("\nType `int' is not 32 bit on this platform. Stopping.\n\n");
      errorFlag = 5;
    }
  if(sizeof(float) != 4)
    {
      printf("\nType `float' is not 32 bit on this platform. Stopping.\n\n");
      errorFlag = 5;
    }
  if(sizeof(double) != 8)
    {
      printf("\nType `double' is not 64 bit on this platform. Stopping.\n\n");
      errorFlag = 5;
    }

  /* Read parameter file */
  if (!errorFlag) 
    {
      nt = 0;
  
      /* Identify tags and set type */
      strcpy(tag[nt], "FilePrefix");
      addr[nt] = All.FilePrefix;
      id[nt++] = STRING;

      strcpy(tag[nt], "FileSuffix");
      addr[nt] = All.FileSuffix;
      id[nt++] = STRING;

      strcpy(tag[nt], "NumDivide");
      addr[nt] = &All.NumDivide;
      id[nt++] = INT;
  
      strcpy(tag[nt], "PeriodicBoundariesOn");
      addr[nt] = &All.PeriodicBoundariesOn;
      id[nt++] = INT;

      strcpy(tag[nt], "BoxSize");
      addr[nt] = &All.BoxSize;
      id[nt++] = FLOAT;

      strcpy(tag[nt], "Border");
      addr[nt] = &All.Border;
      id[nt++] = FLOAT;

      strcpy(tag[nt], "PositionFile");
      addr[nt] = All.PositionFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "PositionFileFormat");
      addr[nt] = &All.PositionFileFormat;
      id[nt++] = INT;
  
      strcpy(tag[nt], "ParticleType");
      addr[nt] = &All.ParticleType;
      id[nt++] = INT;

      strcpy(tag[nt], "BgTreebiNskip");
      addr[nt] = &All.BgTreebiNskip;
      id[nt++] = INT;

      strcpy(tag[nt], "Bgc2HaloId");
      addr[nt] = &All.Bgc2HaloId;
      id[nt++] = INT;

      strcpy(tag[nt], "OutputDir");
      addr[nt] = All.OutputDir;
      id[nt++] = STRING;

      strcpy(tag[nt], "OutputAdjacenciesOn");
      addr[nt] = &All.OutputAdjacenciesOn;
      id[nt++] = INT;

      strcpy(tag[nt], "DecimateInputBy");
      addr[nt] = &All.DecimateInputBy;
      id[nt++] = INT;

      strcpy(tag[nt], "MaxNumSnapshot");
      addr[nt] = &All.MaxNumSnapshot;
      id[nt++] = INT;

      strcpy(tag[nt], "SquishY");
      addr[nt] = &All.SquishY;
      id[nt++] = FLOAT;

      strcpy(tag[nt], "SquishZ");
      addr[nt] = &All.SquishZ;
      id[nt++] = FLOAT;


      /* Open file and begin reading lines */
      if((fd = fopen(fname, "r")))
        {
          while(!feof(fd))
            {
              *buf = 0;
              fgets(buf, 200, fd);
              if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
                continue;

              if(buf1[0] == '%')
                continue;

              for(i = 0, j = -1; i < nt; i++)
                if(strcmp(buf1, tag[i]) == 0)
                  {
                    j = i;
                    tag[i][0] = 0;
                    break;
                  }
              if(j >= 0)
                {
                  switch (id[j])
                    {
                    case FLOAT:
                      *((float *) addr[j]) = atof(buf2);
                      break;
                    case STRING:
                      strcpy(addr[j], buf2);
                      break;
                    case INT:
                      *((int *) addr[j]) = atoi(buf2);
                      break;
                    }
                }
              else
                {
                  fprintf(stdout, "Error in file %s:\n   Tag '%s' not allowed or multiple defined.\n",
                          fname, buf1);
                  errorFlag = 1;
                }
            } 
          fclose(fd);

        }
      else
        {
          printf("\nParameter file %s not found.\n\n",fname);
          errorFlag = 2;
        }

      /* Handle parameters that are not required */
      if(errorFlag != 2) {
        for(i = 0; i < nt; i++) {
	  if(*tag[i]) {
	    if (strcmp(tag[i],"FileSuffix")==0) {
	      /* if (id[i]==STRING) { */
	      printf("String value for tag '%s' is missing. Setting it to an empty string.\n",tag[i]);
	      strcpy(addr[i],"");
	    } else if ((strcmp(tag[i],"PositionFile")==0) && (All.PositionFileFormat==-1)) {
	      strcpy(addr[i],"");
	    } else if ((strcmp(tag[i],"ParticleType")==0) && (All.PositionFileFormat!=1) && (All.PositionFileFormat!=4)) {
	      *((int *) addr[i]) = -1;
	    } else if ((strcmp(tag[i],"BgTreebiNskip")==0) && (All.PositionFileFormat!=2)) {
	      *((int *) addr[i]) = -1;
	    } else if ((strcmp(tag[i],"Bgc2HaloId")==0) && (All.PositionFileFormat!=3)) {
	      *((int *) addr[i]) = -1;
	    
	    } else {
	      printf("Error. I miss a value for tag '%s' in parameter file '%s'.\n", 
		     tag[i], fname);
	      errorFlag = 1;
	    }
	  }
	}
      }
    }

  /* Check parameters */
  if (!errorFlag) 
    {
      /* Number of divides */
      if (All.NumDivide < 1) {
	All.NumDivide = 1;
      }

      /* Periodic boundary conditions */
      if (All.PeriodicBoundariesOn <= 0) {
	All.PeriodicBoundariesOn = 0;
	All.BoxSize = 0.;
      } else {
	All.PeriodicBoundariesOn = 1;
	if (All.BoxSize <= 0) {
	  fprintf(stdout,"BoxSize must be creater than 0.\n");
	  errorFlag = 1;
	}
      }
      
      /* Border */
      if (All.Border <= 0) {
	fprintf(stdout,"Border must be greater than 0.\n");
	errorFlag = 1;
      }

      /* Adjacencies */
      if (All.OutputAdjacenciesOn <= 0) {
	All.OutputAdjacenciesOn = 0;
      } else {
	All.OutputAdjacenciesOn = 1;
      }

      /* Decimate factor */
      if (All.DecimateInputBy <= 1) {
	All.DecimateInputBy = 1;
      }

      /* Add trailing slash to output dir if not present */
      i = strlen(All.OutputDir);
      if(i > 0)
	if(All.OutputDir[i - 1] != '/')
	  strcat(All.OutputDir, "/");
      
      /* Make sure that squish factors not greater than 1 */
      if ((All.SquishY > 1.0) || (All.SquishZ > 1.0)) {
	fprintf(stdout,"SquishY and SquishZ must be less than 1.\n");
	errorFlag = 1;
      }
      if (All.SquishY <= 0) All.SquishY = 1.0;
      if (All.SquishZ <= 0) All.SquishZ = 1.0;

    }

#undef FLOAT
#undef STRING
#undef INT
#undef MAXTAGS

  printf("-----------------------------------------\n");
  return(errorFlag);

}

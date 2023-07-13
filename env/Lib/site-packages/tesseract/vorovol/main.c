#include "qhull_a.h"
#include "allvars.h"
#include "proto.h"

#include "glob.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define DL for (int d=0;d<3;d++) /* Dimension loop */
#define BF 1e30 /* Big factor */

char *strdup(const char *s);

int main(int argc, char *argv[]) {
  int *exitcode;
  char **snaplist;
  char *fmtstr1, *fmtstr2, ftemp[MAXLEN_FILENAME], **fsnap;
  int nsect,errorFlag,i,nsnap,isnap,maxsnap;
  int nerror,nskipped,nullcode=-999;
  glob_t globbuf;
  int ThisTask = 0;
  int FileExtension,FormatFlag=0;

  /* Pars arguments */
  if(argc == 2) {
    strcpy(ParameterFile, argv[1]);
    FileExtension = -1;
  } else if (argc == 3) {
    strcpy(ParameterFile, argv[1]);
    FileExtension = atoi(argv[2]);
  } else {
    printf("Parameter file is missing missing.\n");
    printf("Call with <ParameterFile> [<FileExtension>]\n");
    endrun(1, ThisTask);
  }

  /* Some initialization stuff */
  if (ThisTask == 0) {
    
    /* Read parameter file */
    errorFlag = read_parameter_file(ParameterFile);
    if (errorFlag) endrun(errorFlag, ThisTask);
    nsect = All.NumDivide*All.NumDivide*All.NumDivide;
    
    /* Get file template */
    printf("-----------------------------------------\n");
    printf("Creating format string from %s...\n",All.PositionFile);
    fmtstr1 = strchr(All.PositionFile,'%');
    if (fmtstr1) {
      fmtstr2 = strchr(fmtstr1,'d');
      if (!fmtstr2) {
	printf("Format string must be integer (fmtstr = %s)\n",fmtstr1);
	endrun(1, ThisTask);
      }
      FormatFlag = 1;
      if (FileExtension >= 0) {
	sprintf(ftemp,All.PositionFile,FileExtension);
      } else {
	strncpy(ftemp,All.PositionFile,fmtstr1-All.PositionFile);
	ftemp[fmtstr1-All.PositionFile] = '\0';
	strcat(ftemp,"*");
	strcat(ftemp,fmtstr2+1);
      }
    } else {
      if (FileExtension >= 0) {
	printf("FileExtension provided, but PositionFile in parameters does not contain a format.\n");
	endrun(1, ThisTask);
      }
      strcpy(ftemp,All.PositionFile);
    }

    /* Search for files */
    printf("Searching for files matching %s...\n",ftemp);
    errorFlag = glob(ftemp,GLOB_ERR|GLOB_TILDE,NULL,&globbuf);
    if (errorFlag) {
      printf("glob exited with error code %d\n",errorFlag);
      endrun(errorFlag, ThisTask);
    }
    if ((int)globbuf.gl_pathc == 0) {
      printf("No files found.\n");
      endrun(1, ThisTask);
    } else if (((int)globbuf.gl_pathc > 1) || (FileExtension >= 0)) {
      fmtstr1 = strchr(All.FileSuffix,'%');
      if (!fmtstr1) {
	printf("There are %d files, but the suffix %s does not contain a format.\n",
	       (int)globbuf.gl_pathc,All.FileSuffix);
	globfree( &globbuf );
	endrun(1, ThisTask);
      }
      fmtstr2 = strchr(fmtstr1,'d');
      if (!fmtstr2) {
	printf("File suffix %s does not contain an integer format specifier.\n",
	       All.FileSuffix);
	globfree( &globbuf );
	endrun(1, ThisTask);
      }
    }
    nsnap = (int)globbuf.gl_pathc;
    printf("Found %d files matching %s.\n",nsnap,ftemp);
    printf("-----------------------------------------\n\n");
    fflush(stdout);
    
    /* Copy over files */
    exitcode = (int *)malloc(sizeof(int)*nsnap);
    snaplist = (char **)malloc(sizeof(char *)*nsnap);
    for (fsnap = globbuf.gl_pathv, isnap = 0; *fsnap != NULL; ++fsnap, isnap++ ){
      exitcode[isnap] = nullcode;
      snaplist[isnap] = strdup(*fsnap);
    }
    if( nsnap > 0 ) globfree( &globbuf );

    /* Get maximum number of snapshots */
    if (All.MaxNumSnapshot > 0) {
      maxsnap = All.MaxNumSnapshot;
    } else {
      maxsnap = nsnap;
    }
  }

  /* Loop over files */
  nerror = 0;
  for (isnap = 0; isnap < maxsnap; isnap++) {
    /* Declare private things initialized here */
    float **p;
    int np,iext,isect,sector[3],ThisTask=0;
    char fext[100],blockfile[MAXLEN_FILENAME];
    exitcode[isnap] = 0;

    /* Get file extension */
    if ((nsnap == 1) && (FormatFlag != 1)) {
      iext = 0;
      strcpy(fext,All.FileSuffix);
    } else {
      sscanf(snaplist[isnap],All.PositionFile,&iext);
      sprintf(fext,All.FileSuffix,iext);
    }
    printf("Thread %03d: File %d is %s, iext = %d, fext = %s\n",
	   ThisTask,isnap,snaplist[isnap],iext,fext);
    
    /* Skip reading things if using parts file */
    if (All.PositionFileFormat != -1) {

      /* Read position file */
      np = read_positions_file(snaplist[isnap],&p,ThisTask);
      if (np < 0) {
	exitcode[isnap] = -np;
	np = 0;
      }
      
      if (!exitcode[isnap]) {
	/* Loop over sectors */
	for (isect = 0; isect < nsect ; isect++) {
	  if (exitcode[isnap]) break;
	  sector[0] = isect/(All.NumDivide*All.NumDivide);
	  sector[1] = (isect - (All.NumDivide*All.NumDivide)*sector[0])/All.NumDivide;
	  sector[2] = isect % All.NumDivide;
	  namefile_blk(blockfile,fext,sector[0],sector[1],sector[2]);
	  exitcode[isnap] = voz1b1(sector,blockfile,p,np,ThisTask);
	}
      }
    
      /* Free positions (change this if voztie will use positions directly) */
      if (p != NULL) {
	for (i = 0; i < np ; i++) 
	  if (p[i] != NULL) free(p[i]);
	free(p);
      }
    }

    /* Tie sectors together */
    if (!exitcode[isnap]) exitcode[isnap] = voztie(fext,ThisTask);

    /* Count errors */
    if (exitcode[isnap]) nerror++;

    printf("Thread %03d: exitcode = %d\n",ThisTask,exitcode[isnap]);
    fflush(stdout);
  }

  /* Free file list */
  nskipped = 0;
  if (nerror > 0) 
    printf("\nThe following %d snapshot(s) generated an error:\n",nerror);
  for ( isnap = 0; isnap < nsnap; isnap++) {
    if (exitcode[isnap] == nullcode)
      nskipped++;
    else if (exitcode[isnap] != 0) 
      printf("    (exitcode = %4d) %s\n",exitcode[isnap],snaplist[isnap]);
    free(snaplist[isnap]);
  }
  printf("%d snapshots were not processed.\n",nskipped);
  free(exitcode);
  free(snaplist);

  /* Gracefull exit */
  endrun(0, 0);
}

/* int make_block(int sector[3], float **p, int np, coordT **pblock, int **orig, int ThisTask) { */
/*   float width, width2, totwidth2; */
/*   float guard_space, guard_gap; */
/*   PARTADJ *adjs; */
/*   float center[3], rmin[3], rmax[3]; */
/*   coordT rtemp[3]; */
/*   coordT deladjs[3*MAXVERVER]; */
/*   int nvp, nvpbuf, nvpall, nvpgrd; */
/*   int i,j,inbuff, inmain; */
/*   float predict; */
/*   int exitcode = 0; */
/*   int d; */

/*   /\* Print info to start *\/ */
/*   printf("Thread %03d: Block %02d.%02d.%02d\n", */
/* 	 ThisTask,sector[0],sector[1],sector[2]); */
/*   printf("Thread %03d:     -------------------------\n",ThisTask); */

/*   /\* Determine the size of the box *\/ */
/*   width = 1./(float)All.NumDivide; */
/*   width2 = 0.5*width; */
/*   /\* totwidth = width + 2.*All.Border; *\/ */
/*   totwidth2 = width2 + All.Border; */

/*   /\* Determine spacing of guard points *\/ */
/*   DL center[d] = ((float)sector[d]+0.5)*width; */
/*   guard_space = width/(float)NGUARD; */
/*   if ((All.Border*All.Border - 2.*guard_space*guard_space) < 0.) { */
/*     printf("Thread %03d: border = %f, guard_space = %f.\n", */
/* 	   ThisTask,All.Border,guard_space); */
/*     printf("Thread %03d: Not enough guard points for given border.\nIncrease guards to >= %f.\n", */
/*            ThisTask,sqrt(2.)*width/All.Border); */
/*     exitcode = 1021; */
/*   } */
/*   if (!exitcode) { */
/*     guard_gap = (All.Border / 2.)*(1. + sqrt(1 - 2.*guard_space*guard_space/(All.Border*All.Border))); */
/*     printf("Thread %03d:     guard_space = %f, border = %f, guard_gap = %f.\n", */
/* 	   ThisTask,guard_space,All.Border,guard_gap); */
/*     printf("Thread %03d:     center = (%4.2f, %4.2f, %4.2f)\n", */
/* 	   ThisTask,center[0],center[1],center[2]); */
/*     printf("\n"); */

/*     /\* Count guard particles *\/ */
/*     nvpgrd = 0; */
/*     DL { */
/*       if (All.PeriodicBoundariesOn || (sector[d]!=0)) nvpgrd += (NGUARD+1)*(NGUARD+1); */
/*       if (All.PeriodicBoundariesOn || (sector[d]!=(All.NumDivide-1))) nvpgrd+= (NGUARD+1)*(NGUARD+1); */
/*     } */
/*     printf("Thread %03d:     %d guard points will be added.\n",ThisTask,nvpgrd); */
  
/*     /\* Count particles in this block *\/ */
/*     nvp = 0;    /\* number of particles in sector without buffer *\/ */
/*     nvpbuf = 0; /\* number of particles in sector with buffer *\/ */
    
/*     DL { rmin[d] = BF; rmax[d] = -BF; } */
/*     for (i=0; i<np; i++) { */
      
/*       inbuff = 1; */
/*       inmain = 1; */
/*       DL { */
/* 	/\* printf("%d %d\n",i,d); *\/ */
/* 	/\* printf("%f\n",p[1][d]); *\/ */
/* 	rtemp[d] = (double)p[i][d] - (double)center[d]; */
/* 	if (rtemp[d] < rmin[d]) rmin[d] = rtemp[d]; */
/* 	if (rtemp[d] > rmax[d]) rmax[d] = rtemp[d]; */
/* 	/\* Wrap periodic boundary *\/ */
/* 	if (All.PeriodicBoundariesOn) { */
/* 	  if (rtemp[d] > 0.5) rtemp[d] --; */
/* 	  if (rtemp[d] < -0.5) rtemp[d] ++; */
/* 	} */
/* 	inbuff = inbuff && (fabs(rtemp[d]) < totwidth2); */
/* 	inmain = inmain && (fabs(rtemp[d]) <= width2); */
/*       } */
/*       if (inbuff) nvpbuf++; */
/*       if (inmain) nvp++; */
/*     }   */
/*     nvpbuf += nvpgrd; */
/*     printf("Thread %03d:     %d particles in this block\n",ThisTask,nvp); */
/*     printf("Thread %03d:     xmin=%+4.2f ymin=%+4.2f zmin=%+4.2f\n", */
/* 	   ThisTask,rmin[0],rmin[1],rmin[2]); */
/*     printf("Thread %03d:     xmax=%+4.2f ymax=%+4.2f zmax=%+4.2f\n", */
/* 	   ThisTask,rmax[0],rmax[1],rmax[2]); */
/*     printf("Thread %03d:     -------------------------\n\n",ThisTask); */
/*     fflush(stdout); */
/*   } */

/*   /\* Allocate particles in this sector *\/ */
/*   if (!exitcode) { */
/*     (*pblock) = (coordT **)malloc(3*nvpbuf*sizeof(coordT)); */
/*     if ((*pblock) == NULL) { */
/*       printf("Thread %03d: Unable to allocate pblock\n",ThisTask); */
/*       exitcode = 1100; */
/*     } */
/*   } */
/*   if (!exitcode) { */
/*     orig = (int *)malloc(nvpbuf*sizeof(int)); */
/*     if (orig == NULL) { */
/*       printf("Thread %03d: Unable to allocate orig\n",ThisTask); */
/*       exitcode = 1200; */
/*     } */
/*   } */

/*   /\* Particles in this sector *\/ */
/*   if (!exitcode) { */
/*     printf("Thread %03d:     -------------------------\n",ThisTask); */
/*     printf("Thread %03d:     Selecting particles in this block...\n",ThisTask); */
/*     /\* Select particles in this sector *\/ */
/*     nvp = 0;    /\* number of particles without buffer *\/ */
/*     nvpall = 0; /\* number of particles including buffer and guard *\/ */
/*     DL { rmin[d] = BF; rmax[d] = -BF; } */
/*     for (i=0; i<np; i++) { */
/*       inmain = 1; */
/*       DL { */
/* 	rtemp[d] = p[i][d] - center[d]; */
/* 	/\* Wrap periodic boundary *\/ */
/* 	if (All.PeriodicBoundariesOn) { */
/* 	  if (rtemp[d] > 0.5) rtemp[d] --; */
/* 	  if (rtemp[d] < -0.5) rtemp[d] ++; */
/* 	} */
/* 	inmain = inmain && (fabs(rtemp[d]) <= width2); */
/*       } */
/*       if (inmain) { */
/* 	(*pblock)[3*nvp] = rtemp[0]; */
/* 	(*pblock)[3*nvp+1] = rtemp[1]; */
/* 	(*pblock)[3*nvp+2] = rtemp[2]; */
/* 	orig[nvp] = i; */
/* 	nvp++; */
/* 	DL { */
/* 	  if (rtemp[d] < rmin[d]) rmin[d] = rtemp[d]; */
/* 	  if (rtemp[d] > rmax[d]) rmax[d] = rtemp[d]; */
/* 	} */
/*       } */
/*       else { */
/* 	/\* printf("%d: rtemp=(%f,%f,%f) width2=%f\n",i,rtemp[0],rtemp[1],rtemp[2],width2); *\/ */
/*       } */
/*     } */
/*     printf("Thread %03d:     %d particles in this block (without buffer)\n",ThisTask,nvp); */
/*     printf("Thread %03d:     xmin=%+4.2f ymin=%+4.2f zmin=%+4.2f\n", */
/* 	   ThisTask,rmin[0],rmin[1],rmin[2]); */
/*     printf("Thread %03d:     xmax=%+4.2f ymax=%+4.2f zmax=%+4.2f\n", */
/* 	   ThisTask,rmax[0],rmax[1],rmax[2]); */
/*     printf("Thread %03d:     -------------------------\n\n",ThisTask); */
/*     fflush(stdout); */
  
/*     /\* Buffer particles *\/ */
/*     printf("Thread %03d:     -------------------------\n",ThisTask); */
/*     printf("Thread %03d:     Creating buffer...\n",ThisTask); */
/*     nvpbuf = nvp; */
/*     for (i=0; i<np; i++) { */
/*       inbuff = 1; */
/*       DL { */
/* 	rtemp[d] = p[i][d] - center[d]; */
/* 	/\* Wrap periodic boundary *\/ */
/* 	if (All.PeriodicBoundariesOn) { */
/* 	  if (rtemp[d] > 0.5) rtemp[d] --; */
/* 	  if (rtemp[d] < -0.5) rtemp[d] ++; */
/* 	} */
/* 	inbuff = inbuff && (fabs(rtemp[d])<totwidth2); */
/*       } */
/*       if ((inbuff > 0) && */
/* 	  ((fabs(rtemp[0])>width2)||(fabs(rtemp[1])>width2)||(fabs(rtemp[2])>width2))) { */
/* 	DL (*pblock)[3*nvpbuf+d] = rtemp[d]; */
/* 	orig[nvpbuf] = i; */
/* 	nvpbuf++; */
/* 	DL { */
/* 	  if (rtemp[d] < rmin[d]) rmin[d] = rtemp[d]; */
/* 	  if (rtemp[d] > rmax[d]) rmax[d] = rtemp[d]; */
/* 	} */
/*       } */
/*     } */
/*     nvpall = nvpbuf; */
/*     printf("Thread %03d:     %d particles in this block (with buffer)\n",ThisTask,nvpbuf); */
/*     printf("Thread %03d:     xmin=%+4.2f ymin=%+4.2f zmin=%+4.2f\n", */
/* 	   ThisTask,rmin[0],rmin[1],rmin[2]); */
/*     printf("Thread %03d:     xmax=%+4.2f ymax=%+4.2f zmax=%+4.2f\n", */
/* 	   ThisTask,rmax[0],rmax[1],rmax[2]); */
/*     printf("\n"); */

/*     /\* Predict if the box were uniform density *\/ */
/*     predict = pow(width+2.*All.Border,3)*(float)np; */
/*     printf("Thread %03d:     There should be ~ %f points if the box were uniform; there are %d\n", */
/* 	   ThisTask,predict,nvpbuf); */
/*     printf("Thread %03d:     -------------------------\n\n",ThisTask); */
/*     fflush(stdout); */

/*     /\* Guard points *\/ */
/*     if (nvpgrd > 0) { */
/*       printf("Thread %03d:     -------------------------\n",ThisTask); */
/*       printf("Thread %03d:     Adding %d guard points...\n",ThisTask,nvpgrd); */
/*       /\* Z faces *\/ */
/*       for (i=0; i<NGUARD+1; i++) { */
/* 	for (j=0; j<NGUARD+1; j++) { */
/* 	  /\* Bottom *\/ */
/* 	  if (All.PeriodicBoundariesOn || (sector[2]!=0)) { */
/* 	    (*pblock)[3*nvpall]   = -width2 + (float)i * guard_space; */
/* 	    (*pblock)[3*nvpall+1] = -width2 + (float)j * guard_space; */
/* 	    (*pblock)[3*nvpall+2] = -width2 - guard_gap; */
/* 	    nvpall++; */
/* 	  } */
/* 	  /\* Top *\/ */
/* 	  if (All.PeriodicBoundariesOn || (sector[2]!=(All.NumDivide-1))) { */
/* 	    (*pblock)[3*nvpall]   = -width2 + (float)i * guard_space; */
/* 	    (*pblock)[3*nvpall+1] = -width2 + (float)j * guard_space; */
/* 	    (*pblock)[3*nvpall+2] = width2 + guard_gap; */
/* 	    nvpall++; */
/* 	  } */
/* 	} */
/*       } */
/*       /\* Y faces *\/ */
/*       for (i=0; i<NGUARD+1; i++) { /\* Don't want to overdo the corners*\/ */
/* 	for (j=0; j<NGUARD+1; j++) { */
/* 	  /\* Left *\/ */
/* 	  if (All.PeriodicBoundariesOn || (sector[1]!=0)) { */
/* 	    (*pblock)[3*nvpall]   = -width2 + (float)i * guard_space; */
/* 	    (*pblock)[3*nvpall+1] = -width2 - guard_gap; */
/* 	    (*pblock)[3*nvpall+2] = -width2 + (float)j * guard_space; */
/* 	    nvpall++; */
/* 	  } */
/* 	  /\* Right *\/ */
/* 	  if (All.PeriodicBoundariesOn || (sector[1]!=(All.NumDivide-1))) { */
/* 	    (*pblock)[3*nvpall]   = -width2 + (float)i * guard_space; */
/* 	    (*pblock)[3*nvpall+1] = width2 + guard_gap; */
/* 	    (*pblock)[3*nvpall+2] = -width2 + (float)j * guard_space; */
/* 	    nvpall++; */
/* 	  } */
/* 	} */
/*       } */
/*       /\* X faces *\/ */
/*       for (i=0; i<NGUARD+1; i++) { */
/* 	for (j=0; j<NGUARD+1; j++) { */
/* 	  /\* Front *\/ */
/* 	  if (All.PeriodicBoundariesOn || (sector[0]!=0)) { */
/* 	    (*pblock)[3*nvpall]   = -width2 - guard_gap; */
/* 	    (*pblock)[3*nvpall+1] = -width2 + (float)i * guard_space; */
/* 	    (*pblock)[3*nvpall+2] = -width2 + (float)j * guard_space; */
/* 	    nvpall++; */
/* 	  } */
/* 	  /\* Back *\/ */
/* 	  if (All.PeriodicBoundariesOn || (sector[0]!=(All.NumDivide-1))) { */
/* 	    (*pblock)[3*nvpall]   = width2 + guard_gap; */
/* 	    (*pblock)[3*nvpall+1] = -width2 + (float)i * guard_space; */
/* 	    (*pblock)[3*nvpall+2] = -width2 + (float)j * guard_space; */
/* 	    nvpall++; */
/* 	  } */
/* 	} */
/*       } */
/*       DL { rmin[d] = BF; rmax[d] = -BF; } */
/*       for (i=nvpbuf;i<nvpall;i++) { */
/* 	DL { */
/* 	  if ((*pblock)[3*i+d] < rmin[d]) rmin[d] = (*pblock)[3*i+d]; */
/* 	  if ((*pblock)[3*i+d] > rmax[d]) rmax[d] = (*pblock)[3*i+d]; */
/* 	} */
/*       } */
/*       printf("Thread %03d:     %d particles in this block (with guard points)\n", */
/* 	     ThisTask,nvpall); */
/*       printf("Thread %03d:     There should be %d.\n", */
/* 	     ThisTask,nvpbuf + nvpgrd); */
/*       /\* if (All.PeriodicBoundariesOn) { *\/ */
/*       /\* } *\/ */
/*       printf("Thread %03d:     xmin=%+4.2f ymin=%+4.2f zmin=%+4.2f\n", */
/* 	     ThisTask,rmin[0],rmin[1],rmin[2]); */
/*       printf("Thread %03d:     xmax=%+4.2f ymax=%+4.2f zmax=%+4.2f\n", */
/* 	     ThisTask,rmax[0],rmax[1],rmax[2]); */
/*       printf("Thread %03d:     -------------------------\n\n",ThisTask); */
/*       fflush(stdout); */
/*     } */
/*   } */

/*   endrun(0, 0); */
/* } */

/* int vozblock(int sector[3], char* outfile, float **p, int np, int ThisTask) { */
/*   float width, width2, totwidth2; /\*totwidth2; *\/ */
/*   float guard_space, guard_gap; */
/*   PARTADJ *adjs; */
/*   float center[3], rmin[3], rmax[3]; */
/*   coordT rtemp[3], *parts = 0; */
/*   coordT deladjs[3*MAXVERVER]; */
/*   float *vols = 0; */
/*   int nvp, nvpbuf, nvpall, nvpgrd; */
/*   int i,j,inbuff, inmain; */
/*   int *orig = 0; */
/*   FILE *out; */
/*   double totalvol; */
/*   float predict; */
/*   int exitcode = 0; */
/*   int ninfinite; */
/*   int d; */
/*   int NTask = 1; */

/*   if (vols != NULL) { */
/*     printf("Pointer is not null on initialization\n"); */
/*   } */

/*   /\* Select particles in this block *\/ */
/*   exitcode = make_block(sector, p, np, &parts, &orig, ThisTask); */

/*   /\* Allocate adjacencies *\/ */
/*   if (!exitcode) { */
/*     adjs = (PARTADJ *)malloc(np*sizeof(PARTADJ)); */
/*     if (adjs == NULL) { */
/*       printf("Thread %03d: Unable to allocate adjs\n",ThisTask); */
/*       exitcode = 1300; */
/*     } */
/*   } */

/*   /\* Compute adjacencies *\/ */
/*   if (!exitcode) { */
/*     printf("\n"); */
/*     printf("Thread %03d:     File read.  Tessellating ...\n",ThisTask); */
/*     if (NTask == 1) fflush(stdout); */
/*     exitcode = delaunadj(parts, nvp, nvpbuf, nvpall, All.PeriodicBoundariesOn, &adjs, ThisTask); */
/*     printf("Thread %03d:     Done Tessellating. exitcode=%d\n",ThisTask,exitcode); */
/*     fflush(stdout); */
/*   } */

/*   /\* Allocate volumes *\/ */
/*   if (!exitcode) { */
/*     vols = (float *)malloc(nvp*sizeof(float)); */
/*     if (vols == NULL) { */
/*       printf("Thread %03d: Unable to allocate vols\n",ThisTask); */
/*       exitcode = 1400; */
/*     } */
/*   } */

/*   /\* Calculate volumes*\/ */
/*   if (!exitcode) { */
/*     printf("Thread %03d:     -------------------------\n",ThisTask); */
/*     printf("Thread %03d:     Finding volumes ...\n",ThisTask);  */
/*     if (NTask == 1) fflush(stdout); */
/*     for (i=0; i<nvp; i++) { /\* Just the original particles *\/ */
/*       if (exitcode) break; */
/*       /\* Particles that have adjacencies *\/ */
/*       if (adjs[i].nadj > 0) { */
/* 	for (j = 0; j < adjs[i].nadj; j++) { */
/* 	  for (d = 0; d < 3; d++) { */
/* 	    /\* Distance between particle and adjacent particle *\/ */
/* 	    deladjs[3*j + d] = parts[3*adjs[i].adj[j]+d] - parts[3*i+d]; */
/* 	    /\* Wrap periodic *\/ */
/* 	    if (All.PeriodicBoundariesOn) { */
/* 	      if (deladjs[3*j+d] < -0.5) deladjs[3*j+d]++; */
/* 	      if (deladjs[3*j+d] >  0.5) deladjs[3*j+d]--; */
/* 	    } */
/* 	  } */
/* 	} */
/* 	exitcode = adj2vol(deladjs, adjs[i].nadj, &(vols[i]), ThisTask); */
/* 	if (exitcode) break; */
/* 	vols[i] *= (float)np; /\* Why? *\/ */
/*       } */
/*       /\* Particles bordering upper delaunay facet *\/ */
/*       else if (adjs[i].nadj == -1) { */
/* 	vols[i] = (float)(-1); /\* INFINITY; *\/ */
/*       } */
/*       /\* Particles that have a weird number of adjacencies *\/ */
/*       else { */
/* 	printf("Thread %03d: Particle %d has nadj=%d\n",ThisTask,i,adjs[i].nadj); */
/* 	exitcode = 666; */
/*       } */
/*     } */
/*   } */

/*   /\* Continue if no exit code *\/ */
/*   if (!exitcode) { */
/*     /\* Get the adjacencies back to their original values *\/ */
/*     for (i=0; i<nvp; i++) { */
/*       for (j = 0; j < adjs[i].nadj; j++) { */
/* 	adjs[i].adj[j] = orig[adjs[i].adj[j]]; */
/*       } */
/*     } */
  
/*     /\* Calculate total and average volumne *\/ */
/*     totalvol = 0.; */
/*     ninfinite = 0; */
/*     for (i=0;i<nvp; i++) { */
/*       if (isfinite(vols[i]) && (vols[i]>=0)) { */
/* 	totalvol += (double)vols[i]; */
/*       } */
/*       else ninfinite++; */
/*     } */
/*     printf("Thread %03d:     Total volume = %g\n",ThisTask,totalvol); */
/*     printf("Thread %03d:     Average volume = %g\n",ThisTask,totalvol/(float)nvp); */
/*     printf("Thread %03d:     Number of infinite volumes = %d\n",ThisTask,ninfinite); */
/*     printf("Thread %03d:     -------------------------\n\n",ThisTask); */
/*     fflush(stdout); */
/*   } */

/*   /\* Output *\/ */
/*   if (!exitcode) { */
/*     printf("Thread %03d:     -------------------------\n",ThisTask); */
/*     printf("Thread %03d:     Writing output to %s...\n",ThisTask,outfile); */
/*     out = fopen(outfile,"w"); */
/*     if (out == NULL) { */
/*       printf("Thread %03d: Unable to open %s\n",ThisTask,outfile); */
/*       exitcode = 1001; */
/*     } */
/*   } */
   
/*   /\* Write if no exit code so far *\/ */
/*   if (!exitcode) { */
/*     /\* Info on particles in this file *\/ */
/*     fwrite(&np,1, sizeof(int),out); */
/*     fwrite(&nvp,1, sizeof(int),out); */
/*     /\* Tell us where the original particles were *\/ */
/*     fwrite(orig,sizeof(int),nvp,out); */
/*     /\* Volumes*\/ */
/*     fwrite(vols,sizeof(float),nvp,out); */
/*     /\* Adjacencies *\/ */
/*     for (i=0;i<nvp;i++) { */
/*       fwrite(&(adjs[i].nadj),1,sizeof(int),out); */
/*       if (adjs[i].nadj > 0) */
/* 	fwrite(adjs[i].adj,adjs[i].nadj,sizeof(int),out); */
/*       /\* else printf("0"); *\/ */
/*     } */
/*     fclose(out); */
/*     printf("Thread %03d:     -------------------------\n\n",ThisTask); */
/*     fflush(stdout); */
/*   } */

/*   /\* Free things *\/ */
/*   if (adjs != NULL) { */
/*     for (i = 0; i < nvp; i++) { */
/*       if (adjs[i].nadj > 0) free(adjs[i].adj); */
/*     } */
/*     free(adjs); */
/*   } */
/*   if (parts != NULL) free(parts); */
/*   if (orig  != NULL) free(orig); */
/*   if (vols  != NULL) free(vols); */

/*   return(exitcode); */
/* } */

int voz1b1(int sector[3], char* outfile, float **p, int np, int ThisTask) {
  float width, width2, totwidth2; /*totwidth2; */
  float guard_space, guard_gap;
  PARTADJ *adjs;
  float center[3], rmin[3], rmax[3];
  coordT rtemp[3], *parts = 0;
  coordT deladjs[3*MAXVERVER];
  float *vols = 0;
  int nvp, nvpbuf, nvpall, nvpgrd;
  int i,j,inbuff, inmain;
  int *orig = 0;
  FILE *out;
  double totalvol;
  float predict;
  int exitcode = 0;
  int ninfinite;
  int d;
  int NTask = 1;
  int countnadj2 = 0;

  if (vols != NULL) {
    printf("Pointer is not null on initialization\n");
  }

  /* Print info to start */
  printf("Thread %03d: Block %02d.%02d.%02d\n",
	 ThisTask,sector[0],sector[1],sector[2]);
  printf("Thread %03d:     -------------------------\n",ThisTask);

  /* Determine the size of the box */
  width = 1./(float)All.NumDivide;
  width2 = 0.5*width;
  /* totwidth = width + 2.*All.Border; */
  totwidth2 = width2 + All.Border;

  /* Determine spacing of guard points */
  DL center[d] = ((float)sector[d]+0.5)*width;
  guard_space = width/(float)NGUARD;
  if ((All.Border*All.Border - 2.*guard_space*guard_space) < 0.) {
    printf("Thread %03d: border = %f, guard_space = %f.\n",
	   ThisTask,All.Border,guard_space);
    printf("Thread %03d: Not enough guard points for given border.\nIncrease guards to >= %f.\n",
           ThisTask,sqrt(2.)*width/All.Border);
    exitcode = 1021;
  }
  if (!exitcode) {
    guard_gap = (All.Border / 2.)*(1. + sqrt(1 - 2.*guard_space*guard_space/(All.Border*All.Border)));
    printf("Thread %03d:     guard_space = %f, border = %f, guard_gap = %f.\n",
	   ThisTask,guard_space,All.Border,guard_gap);
    printf("Thread %03d:     center = (%4.2f, %4.2f, %4.2f)\n",
	   ThisTask,center[0],center[1],center[2]);
    printf("\n");

    /* Count guard particles */
    nvpgrd = 0;
    DL {
      if (All.PeriodicBoundariesOn || (sector[d]!=0)) nvpgrd += (NGUARD+1)*(NGUARD+1);
      if (All.PeriodicBoundariesOn || (sector[d]!=(All.NumDivide-1))) nvpgrd+= (NGUARD+1)*(NGUARD+1);
    }
    printf("Thread %03d:     %d guard points will be added.\n",ThisTask,nvpgrd);
  
    /* Count particles in this block */
    nvp = 0;    /* number of particles in sector without buffer */
    nvpbuf = 0; /* number of particles in sector with buffer */
    
    DL { rmin[d] = BF; rmax[d] = -BF; }
    for (i=0; i<np; i++) {
      
      inbuff = 1;
      inmain = 1;
      DL {
	/* printf("%d %d\n",i,d); */
	/* printf("%f\n",p[1][d]); */
	rtemp[d] = (double)p[i][d] - (double)center[d];
	if (rtemp[d] < rmin[d]) rmin[d] = rtemp[d];
	if (rtemp[d] > rmax[d]) rmax[d] = rtemp[d];
	/* Wrap periodic boundary */
	if (All.PeriodicBoundariesOn) {
	  if (rtemp[d] > 0.5) rtemp[d] --;
	  if (rtemp[d] < -0.5) rtemp[d] ++;
	}
	inbuff = inbuff && (fabs(rtemp[d]) < totwidth2);
	inmain = inmain && (fabs(rtemp[d]) <= width2);
      }
      if (inbuff) nvpbuf++;
      if (inmain) nvp++;
    }  
    nvpbuf += nvpgrd;
    printf("Thread %03d:     %d particles in this block\n",ThisTask,nvp);
    printf("Thread %03d:     xmin=%+4.2f ymin=%+4.2f zmin=%+4.2f\n",
	   ThisTask,rmin[0],rmin[1],rmin[2]);
    printf("Thread %03d:     xmax=%+4.2f ymax=%+4.2f zmax=%+4.2f\n",
	   ThisTask,rmax[0],rmax[1],rmax[2]);
    printf("Thread %03d:     -------------------------\n\n",ThisTask);
    fflush(stdout);
  }

  /* Allocate particles in this sector */
  if (!exitcode) {
    parts = (coordT *)malloc(3*nvpbuf*sizeof(coordT));
    if (parts == NULL) {
      printf("Thread %03d: Unable to allocate parts\n",ThisTask);
      exitcode = 1100;
    }
  }
  if (!exitcode) {
    orig = (int *)malloc(nvpbuf*sizeof(int));
    if (orig == NULL) {
      printf("Thread %03d: Unable to allocate orig\n",ThisTask);
      exitcode = 1200;
    }
  }

  if (!exitcode) {
    printf("Thread %03d:     -------------------------\n",ThisTask);
    printf("Thread %03d:     Selecting particles in this block...\n",ThisTask);
    /* Select particles in this sector */
    nvp = 0;    /* number of particles without buffer */
    nvpall = 0; /* number of particles including buffer and guard */
    DL { rmin[d] = BF; rmax[d] = -BF; }
    for (i=0; i<np; i++) {
      inmain = 1;
      DL {
	rtemp[d] = p[i][d] - center[d];
	/* Wrap periodic boundary */
	if (All.PeriodicBoundariesOn) {
	  if (rtemp[d] > 0.5) rtemp[d] --;
	  if (rtemp[d] < -0.5) rtemp[d] ++;
	}
	inmain = inmain && (fabs(rtemp[d]) <= width2);
      }
      if (inmain) {
	parts[3*nvp] = rtemp[0];
	parts[3*nvp+1] = rtemp[1];
	parts[3*nvp+2] = rtemp[2];
	orig[nvp] = i;
	nvp++;
	DL {
	  if (rtemp[d] < rmin[d]) rmin[d] = rtemp[d];
	  if (rtemp[d] > rmax[d]) rmax[d] = rtemp[d];
	}
      }
      else {
	/* printf("%d: rtemp=(%f,%f,%f) width2=%f\n",i,rtemp[0],rtemp[1],rtemp[2],width2); */
      }
    }
    printf("Thread %03d:     %d particles in this block (without buffer)\n",ThisTask,nvp);
    printf("Thread %03d:     xmin=%+4.2f ymin=%+4.2f zmin=%+4.2f\n",
	   ThisTask,rmin[0],rmin[1],rmin[2]);
    printf("Thread %03d:     xmax=%+4.2f ymax=%+4.2f zmax=%+4.2f\n",
	   ThisTask,rmax[0],rmax[1],rmax[2]);
    printf("Thread %03d:     -------------------------\n\n",ThisTask);
    fflush(stdout);
  
    /* Buffer particles */
    printf("Thread %03d:     -------------------------\n",ThisTask);
    printf("Thread %03d:     Creating buffer...\n",ThisTask);
    nvpbuf = nvp;
    for (i=0; i<np; i++) {
      inbuff = 1;
      DL {
	rtemp[d] = p[i][d] - center[d];
	/* Wrap periodic boundary */
	if (All.PeriodicBoundariesOn) {
	  if (rtemp[d] > 0.5) rtemp[d] --;
	  if (rtemp[d] < -0.5) rtemp[d] ++;
	}
	inbuff = inbuff && (fabs(rtemp[d])<totwidth2);
      }
      if ((inbuff > 0) &&
	  ((fabs(rtemp[0])>width2)||(fabs(rtemp[1])>width2)||(fabs(rtemp[2])>width2))) {
	DL parts[3*nvpbuf+d] = rtemp[d];
	orig[nvpbuf] = i;
	nvpbuf++;
	DL {
	  if (rtemp[d] < rmin[d]) rmin[d] = rtemp[d];
	  if (rtemp[d] > rmax[d]) rmax[d] = rtemp[d];
	}
      }
    }
    nvpall = nvpbuf;
    printf("Thread %03d:     %d particles in this block (with buffer)\n",ThisTask,nvpbuf);
    printf("Thread %03d:     xmin=%+4.2f ymin=%+4.2f zmin=%+4.2f\n",
	   ThisTask,rmin[0],rmin[1],rmin[2]);
    printf("Thread %03d:     xmax=%+4.2f ymax=%+4.2f zmax=%+4.2f\n",
	   ThisTask,rmax[0],rmax[1],rmax[2]);
    printf("\n");

    /* Predict if the box were uniform density */
    predict = pow(width+2.*All.Border,3)*(float)np;
    printf("Thread %03d:     There should be ~ %f points if the box were uniform; there are %d\n",
	   ThisTask,predict,nvpbuf);
    printf("Thread %03d:     -------------------------\n\n",ThisTask);
    fflush(stdout);

    /* Guard points */
    if (nvpgrd > 0) {
      printf("Thread %03d:     -------------------------\n",ThisTask);
      printf("Thread %03d:     Adding %d guard points...\n",ThisTask,nvpgrd);
      /* Z faces */
      for (i=0; i<NGUARD+1; i++) {
	for (j=0; j<NGUARD+1; j++) {
	  /* Bottom */
	  if (All.PeriodicBoundariesOn || (sector[2]!=0)) {
	    parts[3*nvpall]   = -width2 + (float)i * guard_space;
	    parts[3*nvpall+1] = -width2 + (float)j * guard_space;
	    parts[3*nvpall+2] = -width2 - guard_gap;
	    nvpall++;
	  }
	  /* Top */
	  if (All.PeriodicBoundariesOn || (sector[2]!=(All.NumDivide-1))) {
	    parts[3*nvpall]   = -width2 + (float)i * guard_space;
	    parts[3*nvpall+1] = -width2 + (float)j * guard_space;
	    parts[3*nvpall+2] = width2 + guard_gap;
	    nvpall++;
	  }
	}
      }
      /* Y faces */
      for (i=0; i<NGUARD+1; i++) { /* Don't want to overdo the corners*/
	for (j=0; j<NGUARD+1; j++) {
	  /* Left */
	  if (All.PeriodicBoundariesOn || (sector[1]!=0)) {
	    parts[3*nvpall]   = -width2 + (float)i * guard_space;
	    parts[3*nvpall+1] = -width2 - guard_gap;
	    parts[3*nvpall+2] = -width2 + (float)j * guard_space;
	    nvpall++;
	  }
	  /* Right */
	  if (All.PeriodicBoundariesOn || (sector[1]!=(All.NumDivide-1))) {
	    parts[3*nvpall]   = -width2 + (float)i * guard_space;
	    parts[3*nvpall+1] = width2 + guard_gap;
	    parts[3*nvpall+2] = -width2 + (float)j * guard_space;
	    nvpall++;
	  }
	}
      }
      /* X faces */
      for (i=0; i<NGUARD+1; i++) {
	for (j=0; j<NGUARD+1; j++) {
	  /* Front */
	  if (All.PeriodicBoundariesOn || (sector[0]!=0)) {
	    parts[3*nvpall]   = -width2 - guard_gap;
	    parts[3*nvpall+1] = -width2 + (float)i * guard_space;
	    parts[3*nvpall+2] = -width2 + (float)j * guard_space;
	    nvpall++;
	  }
	  /* Back */
	  if (All.PeriodicBoundariesOn || (sector[0]!=(All.NumDivide-1))) {
	    parts[3*nvpall]   = width2 + guard_gap;
	    parts[3*nvpall+1] = -width2 + (float)i * guard_space;
	    parts[3*nvpall+2] = -width2 + (float)j * guard_space;
	    nvpall++;
	  }
	}
      }
      DL { rmin[d] = BF; rmax[d] = -BF; }
      for (i=nvpbuf;i<nvpall;i++) {
	DL {
	  if (parts[3*i+d] < rmin[d]) rmin[d] = parts[3*i+d];
	  if (parts[3*i+d] > rmax[d]) rmax[d] = parts[3*i+d];
	}
      }
      printf("Thread %03d:     %d particles in this block (with guard points)\n",
	     ThisTask,nvpall);
      printf("Thread %03d:     There should be %d.\n",
	     ThisTask,nvpbuf + nvpgrd);
      /* if (All.PeriodicBoundariesOn) { */
      /* } */
      printf("Thread %03d:     xmin=%+4.2f ymin=%+4.2f zmin=%+4.2f\n",
	     ThisTask,rmin[0],rmin[1],rmin[2]);
      printf("Thread %03d:     xmax=%+4.2f ymax=%+4.2f zmax=%+4.2f\n",
	     ThisTask,rmax[0],rmax[1],rmax[2]);
      printf("Thread %03d:     -------------------------\n\n",ThisTask);
      fflush(stdout);
    }
  }

  /* Allocate adjacencies */
  if (!exitcode) {
    adjs = (PARTADJ *)malloc(np*sizeof(PARTADJ));
    if (adjs == NULL) {
      printf("Thread %03d: Unable to allocate adjs\n",ThisTask);
      exitcode = 1300;
    }
  }

  /* Tesselate */
  if (!exitcode) {
    printf("\n");
    printf("Thread %03d:     File read.  Tessellating ...\n",ThisTask);
    if (NTask == 1) fflush(stdout);
    exitcode = delaunadj(parts, nvp, nvpbuf, nvpall, All.PeriodicBoundariesOn, &adjs, ThisTask);
    printf("Thread %03d:     Done Tessellating. exitcode=%d\n",ThisTask,exitcode);
    fflush(stdout);
  }

  /* Allocate volumes */
  if (!exitcode) {
    vols = (float *)malloc(nvp*sizeof(float));
    if (vols == NULL) {
      printf("Thread %03d: Unable to allocate vols\n",ThisTask);
      exitcode = 1400;
    }
  }

  /* Calculate volumes*/
  if (!exitcode) {
    printf("Thread %03d:     -------------------------\n",ThisTask);
    printf("Thread %03d:     Finding volumes ...\n",ThisTask); 
    if (NTask == 1) fflush(stdout);
    for (i=0; i<nvp; i++) { /* Just the original particles */
      if (exitcode) break;
      /* Particles that have adjacencies */
      if (adjs[i].nadj > 0) {
	for (j = 0; j < adjs[i].nadj; j++) {
	  for (d = 0; d < 3; d++) {
	    /* Distance between particle and adjacent particle */
	    deladjs[3*j + d] = parts[3*adjs[i].adj[j]+d] - parts[3*i+d];
	    /* Wrap periodic */
	    if (All.PeriodicBoundariesOn) {
	      if (deladjs[3*j+d] < -0.5) deladjs[3*j+d]++;
	      if (deladjs[3*j+d] >  0.5) deladjs[3*j+d]--;
	    }
	  }
	}
	exitcode = adj2vol(deladjs, adjs[i].nadj, &(vols[i]), ThisTask);
	if (exitcode) break;
	vols[i] *= (float)np; /* Why? */
      }
      /* Particles bordering upper delaunay facet */
      else if (adjs[i].nadj == -1) {
	vols[i] = (float)(-1); /* INFINITY; */
      }
      /* Particles with nadj = -2 (invalid vertex?) */
      else if (adjs[i].nadj == -2) {
	countnadj2++;
	printf("Thread %03d: Particle %d has nadj=%d (%d total)\n",
	       ThisTask,i,adjs[i].nadj,countnadj2);
      }
      /* Particles that have a weird number of adjacencies */
      else {
	printf("Thread %03d: Particle %d has nadj=%d\n",ThisTask,i,adjs[i].nadj);
	exitcode = 666;
      }
    }
  }

  /* Continue if no exit code */
  if (!exitcode) {
    /* Get the adjacencies back to their original values */
    for (i=0; i<nvp; i++) {
      for (j = 0; j < adjs[i].nadj; j++) {
	adjs[i].adj[j] = orig[adjs[i].adj[j]];
      }
    }
  
    /* Calculate total and average volumne */
    totalvol = 0.;
    ninfinite = 0;
    for (i=0;i<nvp; i++) {
      if (isfinite(vols[i]) && (vols[i]>=0)) {
	totalvol += (double)vols[i];
      }
      else ninfinite++;
    }
    printf("Thread %03d:     Total volume = %g\n",ThisTask,totalvol);
    printf("Thread %03d:     Average volume = %g\n",ThisTask,totalvol/(float)nvp);
    printf("Thread %03d:     Number of infinite volumes = %d\n",ThisTask,ninfinite);
    printf("Thread %03d:     -------------------------\n\n",ThisTask);
    fflush(stdout);
  }

  /* Output */
  if (!exitcode) {
    printf("Thread %03d:     -------------------------\n",ThisTask);
    printf("Thread %03d:     Writing output to %s...\n",ThisTask,outfile);
    out = fopen(outfile,"w");
    if (out == NULL) {
      printf("Thread %03d: Unable to open %s\n",ThisTask,outfile);
      exitcode = 1001;
    }
  }
   
  /* Write if no exit code so far */
  if (!exitcode) {
    /* Info on particles in this file */
    fwrite(&np,1, sizeof(int),out);
    fwrite(&nvp,1, sizeof(int),out);
    /* Tell us where the original particles were */
    fwrite(orig,sizeof(int),nvp,out);
    /* Volumes*/
    fwrite(vols,sizeof(float),nvp,out);
    /* Adjacencies */
    for (i=0;i<nvp;i++) {
      fwrite(&(adjs[i].nadj),1,sizeof(int),out);
      if (adjs[i].nadj > 0)
	fwrite(adjs[i].adj,adjs[i].nadj,sizeof(int),out);
      /* else printf("0"); */
    }
    fclose(out);
    printf("Thread %03d:     -------------------------\n\n",ThisTask);
    fflush(stdout);
  }

  /* Free things */
  if (adjs != NULL) {
    for (i = 0; i < nvp; i++) {
      if (adjs[i].nadj > 0) free(adjs[i].adj);
    }
    free(adjs);
  }
  if (parts != NULL) free(parts);
  if (orig  != NULL) free(orig);
  if (vols  != NULL) free(vols);

  return(exitcode);
}


/* Tie together blocks */
int voztie(char *suffix, int ThisTask) {
  FILE *part,*vol,*adj;
  char partfile[80], volfile[80], adjfile[80];
  float *vols, volstemp;
  PARTADJ *adjs;
  int i,j,k,p,nout,nsect,isect;
  int np,np2,nvp,nvpmax,nvpsum,npnotdone,*orig,na;
  double avgnadj, avgvol;
  int exitcode = 0;
  int NTask = 1;

  /* Loop over sectors counting particles in each block */
  printf("Thread %03d: -----------------------------\n",ThisTask);
  printf("Thread %03d: Counting particles in each block...\n",ThisTask);
  np = -1; nvpmax = -1; nvpsum = 0;
  nsect = All.NumDivide*All.NumDivide*All.NumDivide;
  for (isect = 0; isect < nsect ; isect++) {
    if (exitcode) break;
    i = isect/(All.NumDivide*All.NumDivide);
    j = (isect - (All.NumDivide*All.NumDivide)*i)/All.NumDivide;
    k = isect % All.NumDivide;
    namefile_blk(partfile,suffix,i,j,k);
    part = fopen(partfile,"r");
    if (part == NULL) {
      printf("Thread %03d: Unable to open file %s.\n\n",ThisTask,partfile);
      exitcode = 2500+isect;
    }
    if (!exitcode) {
      fread(&np2,1,sizeof(int),part);
      fread(&nvp,1,sizeof(int),part);
      if (np == -1)
	np = np2;
      else
	if (np2 != np) {
	  printf("Thread %03d: Incompatible total particle numbers: %d,%d\n\n",
		 ThisTask,np,np2);
	  exitcode = 2600+isect;
	}  
      if (nvp > nvpmax) nvpmax = nvp;
      fclose(part);
    }
  }
  if (!exitcode) {
    printf("Thread %03d: We have %d particles to tie together.\n",ThisTask,np); 
    if (NTask == 1) fflush(stdout);
    printf("Thread %03d: The largest file contains %d particles.\n",ThisTask,nvpmax);
    printf("Thread %03d: -----------------------------\n\n",ThisTask);
    fflush(stdout);
  }

  /* Allocate arrays to be read in */
  if (!exitcode) {
    adjs = (PARTADJ *)malloc(np*sizeof(PARTADJ));
    if (adjs == NULL) {
      printf("Thread %03d: Couldn't allocate adjs.\n",ThisTask);
      exitcode = 2100;
    }
  }
  if (!exitcode) {
    vols = (float *)malloc(np*sizeof(float));
    if (vols == NULL) {
      printf("Thread %03d: Couldn't allocate vols.\n",ThisTask);
      exitcode = 2200;
    }
  }
  if (!exitcode) {
    orig = (int *)malloc(nvpmax*sizeof(int));
    if (orig == NULL) {
      printf("Thread %03d: Couldn't allocate orig.\n",ThisTask);
      exitcode = 2300;
    }
  }

  /* Initialize volumes to -1 */
  if (!exitcode) {
    for (p=0;p<np;p++)
      vols[p] = -1.;
  }

  /* Loop over sectors, reading in particle info */
  printf("Thread %03d: -----------------------------\n",ThisTask);
  printf("Thread %03d: Reading in particle data...\n",ThisTask);
  for (isect = 0; isect < nsect ; isect++) {
    if (exitcode) break;
    i = isect/(All.NumDivide*All.NumDivide);
    j = (isect - (All.NumDivide*All.NumDivide)*i)/All.NumDivide;
    k = isect % All.NumDivide;
    namefile_blk(partfile,suffix,i,j,k);
    printf("\n");
    printf("Thread %03d: Reading sector file %s ...\n",ThisTask,partfile);
    part = fopen(partfile,"r");
    if (part == NULL) {
      printf("Thread %03d: Unable to open file %s.\n\n",ThisTask,partfile);
      exitcode = 2700+isect;
    }
    if (!exitcode) {
      /* Read header info */
      fread(&np2,1,sizeof(int),part);
      fread(&nvp,1,sizeof(int),part);
      nvpsum += nvp;
      
      /* Read original indicies */
      fread(orig,nvp,sizeof(int),part);
      int count_mismatch = 0;
      for (p=0;p<nvp;p++) {
	if (exitcode) break;
	fread(&volstemp,1,sizeof(float),part);
	if (isfinite(volstemp) && (volstemp > -1)) {
	  if (vols[orig[p]] > -1.)
	    if (vols[orig[p]] != volstemp) {
	      printf("Thread %03d: Inconsistent volumes for p. %d: (%g,%g)!\n",
		     ThisTask,orig[p],vols[orig[p]],volstemp);
	      count_mismatch+=1;
	      /* exitcode = 2800+isect; */
	    }
	  /* Copy over and convert units */
	  if (!exitcode)
	    vols[orig[p]] = volstemp*All.BoxSize*All.BoxSize*All.BoxSize/np;
	}
      }
      if (count_mismatch > 0) {
	printf("Thread %03d: %d inconsistent volumes detected.\n",ThisTask,count_mismatch);
	/* exitcode = 2800+isect; */
      }	
    }

    /* Read adjacencies */
    if (!exitcode) {
      for (p=0;p<nvp;p++) {
	if (exitcode) break;
	fread(&na,1,sizeof(int),part);
	adjs[orig[p]].nadj = na;
	if (na > 0) {
	  adjs[orig[p]].adj = (int *)malloc(na*sizeof(int));
	  if (adjs[orig[p]].adj == NULL) {
	    printf("Thread %03d: Couldn't allocate adjs[orig[%d]].adj.\n",ThisTask,p);
	    exitcode = 2400;
	  }
	  if (!exitcode)
	    fread(adjs[orig[p]].adj,na,sizeof(int),part);
	}
      }
    }
    
    /* Close file */
    if (part != NULL) fclose(part);
  }

  /* Print some info on the particles */
  if (!exitcode) {
    printf("\n");
    npnotdone = 0; avgnadj = 0.; avgvol = 0.;
    for (p=0;p<np;p++) {
      if (vols[p] == -1.) npnotdone++;
      else {
	avgnadj += (double)(adjs[p].nadj);
	avgvol += (double)(vols[p]);
      }
    }
    if (npnotdone > 0)
      printf("Thread %03d: %d particles not done or infinite (should be 0 if periodic).\n",
	     ThisTask,npnotdone);
    printf("Thread %03d: %d particles done more than once.\n",
	   ThisTask,nvpsum-np);
    printf("Thread %03d: Average # adjacencies = %lf (%f for Poisson)\n",
	   ThisTask,avgnadj/(double)np,
	   48.*3.141593*3.141593/35.+2.);
    printf("Thread %03d: Total volume = %lf\n",ThisTask,avgvol);
    printf("Thread %03d: Average volume = %lf\n",ThisTask,avgvol/(double)np);
    printf("Thread %03d: -----------------------------\n\n",ThisTask);
    fflush(stdout);
  }

  /* Output adjacencies */
  if (All.OutputAdjacenciesOn && !exitcode) {
    namefile_adj(adjfile,suffix);
    printf("Thread %03d: -----------------------------\n",ThisTask);
    printf("Thread %03d: Outputting adjacencies to %s... \n",ThisTask,adjfile);
    adj = fopen(adjfile,"w");
    if (adj == NULL) {
      printf("Thread %03d: Unable to open %s\n",ThisTask,adjfile);
      exitcode = 2003;
    }
    if (!exitcode) {
      /* Write total number and number of adjacencies for each particle to file */
      fwrite(&np,1, sizeof(int),adj);
      for (i=0;i<np;i++)
	fwrite(&adjs[i].nadj,1,sizeof(int),adj);
      /* Write adjacency arrays to file, only including adjacencies that are further down the list */
      for (i=0;i<np;i++) {
	/* if (i%100000 == 0) printf("Thread %03d: P=%d: Starting to write %d adj\n",ThisTask,i,adjs[i].nadj); */
	if (adjs[i].nadj > 0) {
	  nout = 0;
	  for (j=0;j<adjs[i].nadj; j++) if (adjs[i].adj[j] > i) nout++;
	  /* if (i%100000 == 0) printf("Thread %03d: P=%d: determined nout=%d\n",ThisTask,i,nout); */
	  fwrite(&nout,1,sizeof(int),adj);
	  /* if (i%100000 == 0) printf("Thread %03d: P=%d: wrote nout to file\n",ThisTask,i); */
	  for (j=0;j<adjs[i].nadj; j++)
	    if (adjs[i].adj[j] > i)
	      fwrite(&(adjs[i].adj[j]),1,sizeof(int),adj);
	  /* if (i%100000 == 0) printf("Thread %03d: P=%d: wrote adjs to file\n",ThisTask,i); */
	}
      }
      printf("Thread %03d: Wrote adjs to file\n",ThisTask);
    }

    /* Close file */
    if (adj != NULL) fclose(adj);
    printf("Thread %03d: -----------------------------\n\n",ThisTask);
    fflush(stdout);
  }

  /* Volumes (can write if error only in writing adjacencies) */
  if ((exitcode == 0) || (exitcode == 2003)) {
    namefile_vol(volfile,suffix);
    printf("Thread %03d: -----------------------------\n",ThisTask);
    printf("Thread %03d: Outputting volumes to %s... \n",ThisTask,volfile);
    vol = fopen(volfile,"w");
    if (vol == NULL) {
      printf("Thread %03d: Unable to open %s\n",ThisTask,volfile);
      exitcode = 2004;
    }
    if (exitcode != 2004) {
      fwrite(&np,1, sizeof(int),vol);
      fwrite(vols,sizeof(float),np,vol);
      fclose(vol);
    }
    printf("Thread %03d: -----------------------------\n\n",ThisTask);
    fflush(stdout);
  }

  /* Free things */
  printf("Freeing adjs\n"); fflush(stdout);
  if (adjs != NULL) {
    for (i = 0; i < np; i++) {
      if (adjs[i].nadj > 0 && (adjs[i].adj != NULL)) 
	free(adjs[i].adj);
    }
    free(adjs);
  }
  printf("Freeing vols\n"); fflush(stdout);
  if (vols != NULL) free(vols);
  printf("Freeing orig\n"); fflush(stdout);
  if (orig != NULL) free(orig);

  return(exitcode);
}

void endrun(int ierr, int ThisTask)
{
  if(ierr)
    {
      printf("Thread %03d: endrun called with an error level of %d\n\n\n", ThisTask, ierr);
      fflush(stdout);
      /* MPI_Abort(MPI_COMM_WORLD, ierr); */
      exit(ierr);
    }
  /* MPI_Finalize(); */
  printf("Exiting \n");
  exit(ierr);
}

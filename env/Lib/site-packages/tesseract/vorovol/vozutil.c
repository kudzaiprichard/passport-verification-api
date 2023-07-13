#include "qhull_a.h"
#include "allvars.h"
#include "proto.h"

#define FOREACHvertex2_(vertices) FOREACHsetelement_(vertexT, vertices2,vertex2)

int compar(const void * n1, const void * n2) {
  int i1,i2;
  i1 = *(int *)n1;
  i2 = *(int *)n2;
  return 2*(i1 > i2) - 1 + (i1 == i2);
}

int tessellate(coordT *points, int nvp, int periodic, float **vols, int ThisTask) {
  PARTADJ *adjs;
  coordT deladjs[3*MAXVERVER];
  int i,j;
  double totalvol;
  int exitcode = 0;
  int ninfinite;
  int d;
  int countnadj2 = 0;

  /* Allocate adjacencies */
  if (!exitcode) {
    adjs = (PARTADJ *)malloc(nvp*sizeof(PARTADJ));
    if (adjs == NULL) {
      printf("Thread %03d: Unable to allocate adjs\n",ThisTask);
      exitcode = 1300;
    }
  }
  
  /* Tesselate */
  if (!exitcode) {
    printf("\n");
    printf("Thread %03d:     File read.  Tessellating ...\n",ThisTask);
    fflush(stdout);
    exitcode = delaunadj(points, nvp, nvp, nvp, periodic, &adjs, ThisTask);
    printf("Thread %03d:     Done Tessellating. exitcode=%d\n",ThisTask,exitcode);
    fflush(stdout);
  }

  /* Allocate volumes */
  if (!exitcode) {
    (*vols) = (float *)malloc(nvp*sizeof(float));
    if ((*vols) == NULL) {
      printf("Thread %03d: Unable to allocate vols\n",ThisTask);
      exitcode = 1400;
    }
  }

  /* Calculate volumes*/
  if (!exitcode) {
    printf("Thread %03d:     -------------------------\n",ThisTask);
    printf("Thread %03d:     Finding volumes ...\n",ThisTask);
    for (i=0; i<nvp; i++) { /* Just the original particles */
      if (exitcode) break;
      /* Particles that have adjacencies */
      if (adjs[i].nadj > 0) {
        for (j = 0; j < adjs[i].nadj; j++) {
          for (d = 0; d < 3; d++) {
            /* Distance between particle and adjacent particle */
            deladjs[3*j + d] = points[3*adjs[i].adj[j]+d] - points[3*i+d];
            /* Wrap periodic */
            if (periodic) {
              if (deladjs[3*j+d] < -0.5) deladjs[3*j+d]++;
              if (deladjs[3*j+d] >  0.5) deladjs[3*j+d]--;
            }
          }
        }
        exitcode = adj2vol(deladjs, adjs[i].nadj, &((*vols)[i]), ThisTask);
        if (exitcode) break;
        (*vols)[i] *= (float)nvp; /* Why? */
      /* Particles bordering upper delaunay facet */
      } else if (adjs[i].nadj == -1) {
        (*vols)[i] = (float)(-1); /* INFINITY; */
      /* Particles with nadj = -2 (invalid vertex?) */
      } else if (adjs[i].nadj == -2) {
	countnadj2++;
        printf("Thread %03d: Particle %d has nadj=%d (%d total)\n", 
	       ThisTask,i,adjs[i].nadj,countnadj2);
      /* Particles that have a weird number of adjacencies */
      } else {
        printf("Thread %03d: Particle %d has nadj=%d\n",ThisTask,i,adjs[i].nadj);
        exitcode = 666;
      }
    }
  }

  /* Continue if no exit code */
  if (!exitcode) {
    /* Calculate total and average volumne */
    totalvol = 0.;
    ninfinite = 0;
    for (i=0;i<nvp; i++) {
      if (isfinite((*vols)[i]) && ((*vols)[i]>=0)) {
        totalvol += (double)(*vols)[i];
      }
      else ninfinite++;
    }
    printf("Thread %03d:     Total volume = %g\n",ThisTask,totalvol);
    printf("Thread %03d:     Average volume = %g\n",ThisTask,totalvol/(float)nvp);
    printf("Thread %03d:     Number of infinite volumes = %d\n",ThisTask,ninfinite);
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

  return(exitcode);
}

/* Finds Delaunay adjacencies of a set of points */
int delaunadj (coordT *points, int nvp, int nvpbuf, int nvpall, int periodic, PARTADJ **adjs, int ThisTask) {
  int dim= 3;               /* dimension of points */
  boolT ismalloc= False;    /* True if qhull should free points in qh_freeqhull() or reallocation */
  char flags[250];          /* option flags for qhull, see qh_opt.htm */
  FILE *outfile= stdout;    /* output from qh_produce_output()
                               use NULL to skip qh_produce_output() */
  FILE *errfile= stderr;    /* error messages from qhull code */
  int exitcode;             /* 0 if no error from qhull */
  int curlong, totlong;     /* memory remaining after qh_memfreeshort */

  int i,ver,count;
  int numfacets, numsimplicial, numridges, totneighbors, 
    numcoplanars, numtricoplanars;

  PARTADJ adjst;

  setT *vertices, *vertices2, *vertex_points, *coplanar_points;
  vertexT *vertex, **vertexp;
  vertexT *vertex2, **vertex2p;
  int vertex_i, vertex_n;
  facetT *facet, *neighbor, **neighborp;
  pointT *point, **pointp;

  int nguard_count;

  /* Allocate temporary adjacencies */
  adjst.adj = (int *)malloc(MAXVERVER*sizeof(int));
  if (adjst.adj == NULL) {
    printf("Thread %03d: Unable to allocate adjst.adj\n",ThisTask);
    return(3100);
  }

  /* Delaunay triangulation*/
  /* 'qh facet_list' will contain the convex hull */
  sprintf(flags, "qhull s d");
  printf("Thread %03d: -----------------------------\n",ThisTask);
  exitcode = qh_new_qhull(dim, nvpall, points, ismalloc,
			  flags, outfile, errfile);
  printf("Thread %03d: -----------------------------\n",ThisTask);

  /* if no error */
  if (!exitcode) {
    /* Add vertices */
    qh_countfacets(qh facet_list, NULL, 0, &numfacets, &numsimplicial,
                   &totneighbors, &numridges, &numcoplanars, &numtricoplanars);
    qh_vertexneighbors();
    vertices = qh_facetvertices(qh facet_list, NULL, 0);
    vertex_points = qh_settemp(nvpall);
    coplanar_points= qh_settemp (nvpall);
    qh_setzero(vertex_points, 0, nvpall);
    qh_setzero (coplanar_points, 0, nvpall);
    FOREACHvertex_(vertices) {
      qh_point_add (vertex_points, vertex->point, vertex);
    }
    FORALLfacet_(qh facet_list) {
      FOREACHpoint_(facet->coplanarset) {
        qh_point_add (coplanar_points, point, facet);
      }
    }
    /* Loop over vertices */
    ver = 0;
    FOREACHvertex_i_(vertex_points) {
      (*adjs)[ver].nadj = 0;
      /* Count neighboring vertices and check that they are real */
      adjst.nadj = 0;
      if (vertex) {
	/* Loop over neighboring facets */
        FOREACHneighbor_(vertex) {
          if ((*adjs)[ver].nadj > -1) {
            if (neighbor->visitid) { /* Change this to check for uppder delaunay? */
              vertices2 = neighbor->vertices;
	      /* Loop over facet vertices, adding them if it is not the original vertex */
              FOREACHvertex2_(vertices2) {
                if (ver != qh_pointid(vertex2->point)) {
                  adjst.adj[adjst.nadj] = (int)qh_pointid(vertex2->point);
                  adjst.nadj ++;
                  if (adjst.nadj > MAXVERVER) {
                    printf("Thread %03d: Increase MAXVERVER to at least %d!\n",
			   ThisTask,adjst.nadj);
		    free(adjst.adj);
		    return(3023);
                  }
                }
              }
            } else {
	      if (periodic || !neighbor->upperdelaunay) {
		printf("Thread %03d:     Unreal verticies for ver=%d: visitid=%d upperdelaunay=%d\n",
		       ThisTask,ver,neighbor->visitid,neighbor->upperdelaunay);
	      }
              (*adjs)[ver].nadj = -1; /* There are unreal vertices here */
            }
          }
        }
      } else {
	/* */
	(*adjs)[ver].nadj = -2;
      }

      /* Enumerate the unique adjacencies*/
      if (adjst.nadj >= 4) {
        qsort((void *)adjst.adj, adjst.nadj, sizeof(int), &compar);
        count = 1;
        nguard_count = 0;
        for (i=1; i<adjst.nadj; i++)
          if (adjst.adj[i] != adjst.adj[i-1]) {
            if (adjst.adj[i] >= nvpbuf) {
              nguard_count++;
              /* printf("Thread %03d:     Guard point encountered.  Increase border and/or nguard.\n",ThisTask); */
              /* printf("Thread %03d:         P:(%f,%f,%f), G: (%f,%f,%f)\n",ThisTask, */
              /*        points[3*ver],points[3*ver+1],points[3*ver+2], */
              /*        points[3*adjst.adj[i]],points[3*adjst.adj[i]+1],points[3*adjst.adj[i]+2]); */
            }
            count++;
          }
        (*adjs)[ver].adj = (int *)malloc(count*sizeof(int));
        if ((*adjs)[ver].adj == NULL) {
          printf("Thread %03d: Unable to allocate (*adjs)[ver].adj\n",ThisTask);
	  free(adjst.adj);
	  return(3110);
        }
        (*adjs)[ver].adj[0] = adjst.adj[0];
        count = 1;
        for (i=1; i<adjst.nadj; i++)
          if (adjst.adj[i] != adjst.adj[i-1]) {
            (*adjs)[ver].adj[count] = adjst.adj[i];
            count++;
          }
        (*adjs)[ver].nadj = count;
        if (nguard_count > 0) {
          printf("Thread %03d:     %d Guard points encountered.  Increase border and/or nguard.\n",
		 ThisTask,nguard_count);
          printf("Thread %03d:         P:(%f,%f,%f)\n",
		 ThisTask,points[3*ver],points[3*ver+1],points[3*ver+2]);
        }
      } else {
	if (periodic) {
	  printf("Thread %03d: Number of adjacencies %d < 4, particle %d -> %d\n",
		 ThisTask,adjst.nadj,ver,ver);
	  free(adjst.adj);
	  return(3024);
	}
      }
      ver++;
      if (ver == nvp) break;
    }
    qh_settempfree (&coplanar_points);
    qh_settempfree (&vertex_points);
    qh_settempfree (&vertices);
  }
  qh_freeqhull(!qh_ALL);                 /* free long memory */
  qh_memfreeshort (&curlong, &totlong);  /* free short memory and memory allocator */
  if (curlong || totlong)
    fprintf (errfile, "Thread %03d: qhull internal warning (delaunadj): did not free %d bytes of long memory (%d pieces)\n", 
	     ThisTask, totlong, curlong);
  free(adjst.adj);
  return(exitcode);
}


/* Calculates the Voronoi volume from a set of Delaunay adjacencies */
int adj2vol (coordT *deladjs, int numpoints, float *vol, int ThisTask) {
  coordT points[3*MAXVERVER];
  pointT intpoints[3*MAXVERVER];
  int dim= 3;               /* dimension of points */
  boolT ismalloc= False;    /* True if qhull should free points in qh_freeqhull() or reallocation */
  char flags[250];          /* option flags for qhull, see qh_opt.htm */
  FILE *outfile= NULL;      /* output from qh_produce_output()
			     use NULL to skip qh_produce_output() */
  FILE *errfile= stderr;    /* error messages from qhull code */
  int exitcode;             /* 0 if no error from qhull */
  facetT *facet;            /* set by FORALLfacets */
  int curlong, totlong;     /* memory remaining after qh_memfreeshort */

  coordT *point, *normp, *coordp, *feasiblep, *deladj;
  int i, j, k;
  boolT zerodiv;
  float runsum;

  /* Create point array (not sure why) */
  for (i=0; i<numpoints; i++) {
    runsum = 0.;
    deladj = deladjs + dim*i;
    point = points + (dim+1)*i;
    for (j=0;j<dim;j++) {
      runsum += deladj[j]*deladj[j];
      point[j] = deladj[j];
    }
    point[dim] = -0.5*runsum;
  }

  /* Create convex hull for just these adjacencies */
  sprintf(flags, "qhull H0");
  exitcode = qh_new_qhull(dim+1, numpoints, points, ismalloc,
			  flags, outfile, errfile);

  numpoints = 0;
  /* if no error */
  if (!exitcode) {
    /* Loop over facets counting them*/
    FORALLfacets {
      numpoints++;
    }

    /* Loop over facets creating a list of half space points */
    j = 0;
    FORALLfacets {
      /* Check that there's a point inside all of the half spaces */
      if (!qh feasible_point) {
        fprintf(stdout, "Thread %03d: qhull input error (qh_printafacet): option 'Fp' needs qh feasible_point\n",
		ThisTask);
        qh_errexit(qh_ERRinput, NULL, NULL);
      }
      /* Select pointer for the point associated with this half space */
      point = coordp = intpoints + j*dim;
      j++;
      normp = facet->normal;
      feasiblep = qh feasible_point;
      /* Compute the coordinate of the halfspace point */
      if (facet->offset < -qh MINdenom) {
        for (k= qh hull_dim; k--; )
          *(coordp++) = (*(normp++) / - facet->offset) + *(feasiblep++);
      }
      else {
        for (k= qh hull_dim; k--; ) {
          *(coordp++) = qh_divzero(*(normp++), facet->offset, qh MINdenom_1,
                                   &zerodiv) + *(feasiblep++);
          if (zerodiv) {
            qh_memfree (point, qh normal_size);
            printf("Thread %03d: LABELprintinfinite\n",ThisTask);
	    return(3022);
          }
        }
      }
    }
  }
  qh_freeqhull (!qh_ALL);
  qh_memfreeshort (&curlong, &totlong);

  /* Now we calculate the volume from the half space points */
  sprintf (flags, "qhull FA");
  exitcode = qh_new_qhull(dim, numpoints, intpoints, ismalloc,
			  flags, outfile, errfile);
  qh_getarea(qh facet_list);
  *vol = qh totvol;

  /* Free memory */
  qh_freeqhull (!qh_ALL);
  qh_memfreeshort (&curlong, &totlong);
  if (curlong || totlong)
    fprintf (errfile, "Thread %03d: qhull internal warning (adj2vol): did not free %d bytes of long memory (%d pieces)\n", 
	     ThisTask,totlong, curlong);
  /*free(points); free(intpoints);*/

  return(exitcode);
}

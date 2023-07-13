/*! \file proto.h
 *  \brief this file contains all function prototypes of the code
 */

int voz1b1(int sector[3], char *outfile, float **p, int np, int ThisTask);
int voztie(char *suffix, int ThisTask);

int tessellate(coordT *points, int nvp, int periodic, float **vols, int ThisTask);
int delaunadj (coordT *points, int nvp, int nvpbuf, int nvpall, int periodic, PARTADJ **adjs, int ThisTask);
int adj2vol (coordT *deladjs, int numpoints, float *vol, int ThisTask);

int read_parameter_file(char *fname);
int read_positions_file(char *fname, float ***p, int ThisTask);

void namefile_blk(char *filename, char *suffix, int xsplit, int ysplit, int zsplit);
void namefile_adj(char *filename, char *suffix);
void namefile_vol(char *filename, char *suffix);

void endrun(int ierr, int ThisTask);

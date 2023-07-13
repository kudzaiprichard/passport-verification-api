#define MAXVERVER 2000
#define NGUARD 42 /*Actually, the number of SPACES between guard points
		    in each dim */
#define MAXLEN_FILENAME  1000

extern char ParameterFile[MAXLEN_FILENAME];

/* extern float **P; */

typedef struct Partadj 
{
  int nadj;
  int *adj;
} 
  PARTADJ;

extern struct parameters
{
  char FilePrefix[MAXLEN_FILENAME];
  char FileSuffix[MAXLEN_FILENAME];
  int NumDivide;
  int PeriodicBoundariesOn;
  float Border;
  float BoxSize;
  char PositionFile[MAXLEN_FILENAME];
  int PositionFileFormat;
  int ParticleType;
  int BgTreebiNskip;
  int Bgc2HaloId;
  char OutputDir[MAXLEN_FILENAME];
  int OutputAdjacenciesOn;
  int DecimateInputBy;
  int MaxNumSnapshot;
  float SquishY;
  float SquishZ;
}
  All;

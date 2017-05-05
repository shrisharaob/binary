//#include <string>
//#include <cstring>

#define NE 100
#define NI 100
#define NFF 100
#define K 10
#define SQRT_K (sqrt((double)K))
#define THRESHOLD_E 1.0
#define THRESHOLD_I 1.0
#define JEE 1.0
#define JIE 1.0
#define JEI -1.5
#define JII -1.0
#define JE0 2.0
#define JI0 1.0
#define TAU_FF 1.0
#define TAU_E 1.0
#define TAU_I 1.0
#define JE0_K (JE0 / sqrt((double)K))
#define JI0_K (JI0 / sqrt((double)K))
#define JEE_K (JEE / sqrt((double)K))
#define JEI_K (JEI / sqrt((double)K))
#define JIE_K (JIE / sqrt((double)K))
#define JII_K (JII / sqrt((double)K))
#define N_NEURONS (NE + NI)

#define cFF 1.0 //0.20 // FOR LATER USE !

#define T_STOP 100
#define T_TRANSIENT (T_STOP / 2)
#define IF_GEN_MAT_DEFAULT 0

double m0, tStop, recModulation, ffModulation;
double *conMat, *conMatFF;
int trialNumber;
unsigned int *nPostNeurons, *sparseConVec, *idxVec;
unsigned int *nPostNeuronsFF, *sparseConVecFF, *idxVecFF;

int IF_GEN_MAT;
unsigned long long int nSteps;

double phi_ext, m0_ext, m1_ext;

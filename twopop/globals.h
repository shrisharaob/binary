//#include <string>
//#include <cstring>

#define NE 10000
#define NI 10000
#define NFF 10000
#define K 1000
#define SQRT_K (sqrt((double)K))
#define THRESHOLD_E 1.0
#define THRESHOLD_I 1.0
#define JEE 1.0
#define JIE 1.0
#define JEI -1.5
#define JII -1.0
#define JE0 2.0
#define JI0 1.0
#define TAU_FF 4.0
#define TAU_E 4.0
#define TAU_I 2.0
#define JE0_K (JE0 / (cFF * sqrt((double)K)))
#define JI0_K (JI0 / (cFF * sqrt((double)K)))
#define JEE_K (JEE / sqrt((double)K))
#define JEI_K (JEI / sqrt((double)K))
#define JIE_K (JIE / sqrt((double)K))
#define JII_K (JII / sqrt((double)K))
#define N_NEURONS (NE + NI)

#define cFF 0.2
#define SQRT_KFF (sqrt((double)K * cFF))

#define T_STOP 10000
#define T_TRANSIENT (T_STOP * 0.05)
#define IF_GEN_MAT_DEFAULT 0
#define IF_SAVE_SPKS 0
#define IF_STEP_PHI0 0

#define N_SEGMENTS 20

double m0, tStop, recModulation, ffModulation;
double *conMat, *conMatFF;
int trialNumber;
unsigned int *nPostNeurons, *sparseConVec, *idxVec;
unsigned int *nPostNeuronsFF, *sparseConVecFF, *idxVecFF;

int IF_GEN_MAT;
unsigned long long int nSteps;

double phi_ext, m0_ext, m1_ext;

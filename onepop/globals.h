//#include <string>
//#include <cstring>

#define NE 0U
#define NI 1000U
#define K 100
#define SQRT_K (sqrt((double)K))
#define THRESHOLD_E 1.0
#define THRESHOLD_I 1.0
#define JEE 0 //1.0
#define JIE 0 //1.0
#define JEI 0 //-1.5
#define JII -1.0
#define JE0 0 //2.0
#define JI0 1.0
#define TAU_E 1e-10
#define TAU_I 1.0
#define JEE_K (JEE / sqrt(K))
#define JEI_K (JEI / sqrt(K))
#define JIE_K (JIE / sqrt(K))
#define JII_K (JII / sqrt(K))
#define N_NEURONS (NE + NI)
#define cFF 0.0 // FOR LATER USE !
#define T_STOP 10000
#define T_TRANSIENT (T_STOP / 2)
#define IF_GEN_MAT_DEFAULT 0

double m0, tStop, recModulation, ffModulation;
double *conMat;
int trialNumber;
unsigned int *nPostNeurons, *sparseConVec, *idxVec;
int IF_GEN_MAT, IF_SAVE_MAT;
unsigned long long int nSteps;

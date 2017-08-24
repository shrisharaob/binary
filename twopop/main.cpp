#include <stdio.h>
#include <math.h>
// #include <boost/random/uniform_01.hpp>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <random>
#include <ctime>
#include <string>
#include <cstring>
// #include <boost/filesystem>
#include "globals.h"

using namespace::std ;
std::string folderName;

void ProgressBar(float progress, float me, float mi) {
    int barWidth = 31;
    std::cout << "Progress: [";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
      if (i < pos)  std::cout << "\u25a0"; //std::cout << "=";
      else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << "% done | mE = " << me << " mI = " << mi << "\r";
    std::cout.flush();
    if(progress == 1.) std::cout << std::endl;
}

void ProgressBar(float progress, float me, float mi, float m0) {
    int barWidth = 31;
    std::cout << "Progress: [";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
      if (i < pos)  std::cout << "\u25a0"; //std::cout << "=";
      else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << "% done | m0 = " << m0 << " mE = " << me << " mI = " << mi << "\r";
    std::cout.flush();
    if(progress == 1.) std::cout << std::endl;
}


void M1ComponentI(vector<double> &x, unsigned int n, double* m1, double* phase) {
  double dPhi = M_PI / (double)n;
  double xCord = 0, yCord = 0;
  for(unsigned int i = NE; i < NE + n; i++) {
    xCord += x[i] * cos(2.0 * i * dPhi);
    yCord += x[i] * sin(2.0 * i * dPhi);
  }
  *m1 = (2.0 / (double)n) * sqrt(xCord * xCord + yCord * yCord);
  *phase = 0.5 * atan2(yCord, xCord);
  if(*phase < 0) {
    *phase = *phase + M_PI;
  }
}


void M1ComponentE(vector<double> &x, unsigned int n, double* m1, double* phase) {
  double dPhi = M_PI / (double)n;
  double xCord = 0, yCord = 0;
  for(unsigned int i = 0; i < n; i++) {
    xCord += x[i] * cos(2.0 * i * dPhi);
    yCord += x[i] * sin(2.0 * i * dPhi);
  }
  *m1 = (2.0 / (double)n) * sqrt(xCord * xCord + yCord * yCord);
  *phase = 0.5 * atan2(yCord, xCord);
  if(*phase < 0) {
    *phase = *phase + M_PI;
  }
}


double UniformRand() {
    // printf("\n hello from space 0\n");
    return (double)rand() / (double)RAND_MAX ;
}

void Dummy() {
  printf("me too\n");
}

unsigned int RandFFNeuron() { //unsigned int min, unsigned int max) {
  unsigned int result = NFF;
  while(result == NFF) {
    result = (unsigned int) (rand() / (double) (RAND_MAX + 1UL) * (NFF + 1));
  }
  return result;
}

unsigned int RandENeuron() { //unsigned int min, unsigned int max) {
  unsigned int result = NE;
  while(result == NE) {
    result = (unsigned int) (rand() / (double) (RAND_MAX + 1UL) * (NE + 1));
  }
  return result;
}

unsigned int RandINeuron() { //unsigned int min, unsigned int max) {
  unsigned int result = NI;
  while(result == NI) {
    result = NE + (unsigned int) (rand() / (double) (RAND_MAX + 1UL) * (NI + 1));
  }
  return result;  
}


// double ConProb(double phiI, double phiJ, unsigned int N) {
//   double out = ((double)K / (double)N) * (1 + 2.0 * (recModulation / SQRT_K) * cos(2.0 * (phiI - phiJ)));
//   if((out < 0) | (out > 1)) {
//     cout << "connection probability not in range [0, 1]!" << endl;
//     exit(1);
//   }
//   return out;
// }


double ConProb(double phiI, double phiJ, unsigned int N, double recModulationAB) {
  double out = ((double)K / (double)N) * (1 + 2.0 * (recModulationAB / SQRT_K) * cos(2.0 * (phiI - phiJ)));
  if((out < 0) | (out > 1)) {
    cout << "connection probability not in range [0, 1]!" << endl;
    exit(1);
  }
  return out;
}


double ConProbFF(double phiI, double phiJ, unsigned int N) {
  double out = ((double)K * cFF / (double)N) * (1 + 2.0 * (ffModulation / SQRT_KFF) * cos(2.0 * (phiI - phiJ)));
  if((out < 0) | (out > 1)) {
    cout << "connection probability not in range [0, 1]!" << endl;
    exit(1);
  }
  return out;
}

double FFTuningCurve(unsigned int i, double phiOft) {
  double out = m0_ext + m1_ext * cos(2 * (phiOft - (double)i * M_PI / (double)NFF));
  if((out < 0) | (out > 1)) {
    cout << "external rate [0, 1]!" << endl;
    exit(1);
  }
  //  out = m0_ext;
  return out;
}

void GenSparseMat(unsigned int *conVec,  unsigned int rows, unsigned int clms, unsigned int* sparseVec, unsigned int* idxVec, unsigned int* nPostNeurons ) {
  /* generate sparse representation
     conVec       : input vector / flattened matrix 
     sparseVec    : sparse vector
     idxVec       : every element is the starting index in sparseVec for ith row in matrix conVec
     nPostNeurons : number of non-zero elements in ith row 
  */
  printf("\n MAX Idx of conmat allowed = %u \n", rows * clms);
  unsigned long long int i, j, counter = 0, nPost;
  for(i = 0; i < rows; ++i) {
    nPost = 0;
    for(j = 0; j < clms; ++j) {
      // printf("%llu %llu %llu %llu\n", i, j, i + clms * j, i + rows * j);
      if(conVec[i + rows * j]) { /* i --> j  */
        sparseVec[counter] = j;
        counter += 1;
        nPost += 1;
      }
    }
    nPostNeurons[i] = nPost;
  }
  
  idxVec[0] = 0;
  for(i = 1; i < rows; ++i) {
    idxVec[i] = idxVec[i-1] + nPostNeurons[i-1];
  }
  
  // FILE *fp = NULL;
  // fp = fopen("check_prob.csv", "w");
  // unsigned int ccount = 0;
  // vector <unsigned int> shiftedVec(NE);
  // vector <unsigned int> shiftedVecBuff(NE);
  // for(i = 0; i < N_NEURONS; ++i) {
  //   ccount = 0;
  //   for(j =0; j < NE; ++j) {
  //     shiftedVecBuff[j] = conVec[j + i * N_NEURONS];
  //   }
  //   std::rotate(shiftedVecBuff.begin(), shiftedVecBuff.begin() + i, shiftedVecBuff.end());
  //   for(unsigned int l = 0; l < NE; l++) {
  //     shiftedVec[l] += shiftedVecBuff[l];
  //   }
  // }
  
  // FILE *fp2 = NULL;
  // fp2 = fopen("prob.csv", "w");
  // for(i = 0; i < NE; i++) {
  //   // printf("%llu:\n", i); //    %f\n", i, shiftedVec[i] / (double)NE);    
  //   //fprintf(fp, "%f\n", shiftedVec[i] / (double)NE);
  //   fprintf(fp2, "%f\n", ConProb(0.0, (double)i * M_PI / (double)NE, NE, recModulationEE));
  // }
  // fclose(fp2);
  // shiftedVecBuff.clear();
  // shiftedVec.clear();
}

void GenFFConMat() {
  std::random_device rd;
  std::default_random_engine gen(rd());
  std::uniform_real_distribution<double> UniformRand(0.0, 1.0);
  unsigned long long int nConnections = 0;
  cout << "generating FF conmat" << endl;
  unsigned int *conMatFF = new unsigned int [(unsigned long int)NFF * N_NEURONS];
  for (unsigned long int i = 0; i < NFF; i++)  {
    for (unsigned long int j = 0; j < N_NEURONS; j++)  {
      // i --> j
      if(j < NE) { //E-to-E
	if(UniformRand(gen) <= ConProbFF(i * M_PI / (double)NFF, j * M_PI / (double)NE, NFF)) {
	  conMatFF[i + NFF * j] = 1;
  	  // conMatFF[i + NFF * j] = 1;
	  nConnections += 1;
	}
      }
      else {
	if(UniformRand(gen) <= (double)K * cFF / (double)NFF) {
	  conMatFF[i + NFF * j] = 1;
	  nConnections += 1;
	}
      }
    }
  }
  cout << "done" << endl;
  cout << "computing sparse rep" << endl;    
  sparseConVecFF = new unsigned int[nConnections];
  GenSparseMat(conMatFF, NFF, N_NEURONS, sparseConVecFF, idxVecFF, nPostNeuronsFF);
  cout << "done" << endl;


  
  FILE *ffp;
  ffp = fopen("kcount_ff.csv", "w");
  for(unsigned int lll = 0; lll < NFF; lll++) {
    fprintf(ffp, "%u\n", nPostNeuronsFF[lll]);
  }
  fclose(ffp);
  delete [] conMatFF;

 
  
  FILE *fpSparseConVec, *fpIdxVec, *fpNpostNeurons;
  unsigned long int nElementsWritten;
  printf("done\n#connections = %llu\n", nConnections);
  printf("writing to file ... "); fflush(stdout);
  fpSparseConVec = fopen("sparseConVecFF.dat", "wb");
  nElementsWritten = fwrite(sparseConVecFF, sizeof(*sparseConVecFF), nConnections, fpSparseConVec);
  fclose(fpSparseConVec);
  if(nElementsWritten != nConnections) {
    printf("\n Error: All elements not written \n");
  }
  fpIdxVec = fopen("idxVecFF.dat", "wb");
  fwrite(idxVecFF, sizeof(*idxVecFF), NFF,  fpIdxVec);
  fclose(fpIdxVec);
  fpNpostNeurons = fopen("nPostNeuronsFF.dat", "wb");
  fwrite(nPostNeuronsFF, sizeof(*nPostNeuronsFF), NFF, fpNpostNeurons);
  fclose(fpNpostNeurons);
  printf("done\n");



  //-----------------------------------------
  // printf("testing sparsevecff read\n");
  // unsigned long int dummy;
  // FILE *fpSparseConVecFF = fopen("sparseConVecFF.dat", "rb");
  // dummy = fread(sparseConVecFF, sizeof(*sparseConVecFF), nConnections, fpSparseConVecFF);
  //   printf("sparseConvec read: %lu %llu \n", dummy, nConnections);
  // if(dummy != nConnections) {
  //   printf("sparseConvec read error ? %lu %llu \n", dummy, nConnections);
  // fclose(fpSparseConVecFF);
  // }
  
}
    
void GenConMat() {
  std::random_device rd;
  std::default_random_engine gen(rd());
  std::uniform_real_distribution<double> UniformRand(0.0, 1.0);
  //  int dummy;
  unsigned long long int nConnections = 0;
  cout << "generating conmat" << endl;
  unsigned int *conMat = new unsigned int [(unsigned long int)N_NEURONS * N_NEURONS];

  for (unsigned long int i = 0; i < NE; i++)  {
    for (unsigned long int j = 0; j < NE; j++)  {
      conMat[i + N_NEURONS * j] = 0;
    }
  }
  
  for (unsigned long int i = 0; i < N_NEURONS; i++)  {
    for (unsigned long int j = 0; j < N_NEURONS; j++)  {
      // i --> j
      if(i < NE && j < NE) { //E-to-E
	if(UniformRand(gen) <= ConProb(i * M_PI / (double)NE, j * M_PI / (double)NE, NE, recModulationEE)) {
	  conMat[i + N_NEURONS * j] = 1;
	  nConnections += 1;
	}
      }
      if(i < NE && j >= NE) { //E-to-I
	if(UniformRand(gen) <= ConProb(i * M_PI / (double)NE, j * M_PI / (double)NI, NE, recModulationIE)) {
	  conMat[i + N_NEURONS * j] = 1;
	  nConnections += 1;
	}
      }
      if(i >= NE && j < NE) { //I-to-E
	if(UniformRand(gen) <= ConProb(i * M_PI / (double)NI, j * M_PI / (double)NE, NI, recModulationEI)) {
	  conMat[i + N_NEURONS * j] = 1;
	  nConnections += 1;
	}
      }
      if(i >= NE && j >= NE) { //I-to-I
	if(UniformRand(gen) <= (double)K / (double)NI) {
	  conMat[i + N_NEURONS * j] = 1;
	  nConnections += 1;
	}
      }
    }
  }


  
  cout << "done" << endl;


  for (unsigned long int i = 0; i < NE; i++)  {
    for (unsigned long int j = 0; j < NE; j++)  {
      if((conMat[i + N_NEURONS * j] < 0) || (conMat[i + N_NEURONS * j] > 1)) {
  	  printf("%lu %lu !!!!!!! oh la la !!!!\n", i, j);
  	}
    }
  }

  
  unsigned long int nElementsWritten;
  // FILE *fpcmat;
  // fpcmat = fopen("cmat.dat", "wb");
  // nElementsWritten = fwrite(conMat, sizeof(unsigned int), N_NEURONS * N_NEURONS, fpcmat);
  // fclose(fpcmat);

  // for (unsigned long int i = 0; i < NE; i++)  {
  //   printf("\n");    
  //   for (unsigned long int j = 0; j < NE; j++)  {
  //     printf("%u ", conMat[i + N_NEURONS * j]);
  //     // printf("%lu: %u\n", i + N_NEURONS * j, conMat[i + N_NEURONS * j]);
  //   }
  // }
  // printf("\n\n");


  
  cout << "computing sparse rep" << endl;    
  sparseConVec = new unsigned int[nConnections];
  GenSparseMat(conMat, N_NEURONS, N_NEURONS, sparseConVec, idxVec, nPostNeurons);
  cout << "done" << endl;
  FILE *ffp;
  ffp = fopen("kcount.csv", "w");
  for(unsigned int lll = 0; lll < N_NEURONS; lll++) {
    fprintf(ffp, "%u\n", nPostNeurons[lll]);
  }
  fclose(ffp);
  delete [] conMat;
  // write to file
  FILE *fpSparseConVec, *fpIdxVec, *fpNpostNeurons;

  printf("done\n#connections = %llu\n", nConnections);
  printf("writing to file ... "); fflush(stdout);
  fpSparseConVec = fopen("sparseConVec.dat", "wb");
  nElementsWritten = fwrite(sparseConVec, sizeof(*sparseConVec), nConnections, fpSparseConVec);
  fclose(fpSparseConVec);
  if(nElementsWritten != nConnections) {
    printf("\n Error: All elements not written \n");
  }
  fpIdxVec = fopen("idxVec.dat", "wb");
  fwrite(idxVec, sizeof(*idxVec), N_NEURONS,  fpIdxVec);
  fclose(fpIdxVec);
  fpNpostNeurons = fopen("nPostNeurons.dat", "wb");
  fwrite(nPostNeurons, sizeof(*nPostNeurons), N_NEURONS, fpNpostNeurons);
  fclose(fpNpostNeurons);
  printf("done\n");


  // for (unsigned long int i = 0; i < 10; i++)  {
  //     printf("%u ", sparseConVec[i]);
  // }
  // exit(1);  
}

void VectorSum(vector<double> &a, vector<double> &b) { 
  for(unsigned int i = 0; i < N_NEURONS; ++i) {
    a[i] += b[i];
  }
}

void VectorSumFF(vector<double> &a, vector<double> &b) { 
  for(unsigned int i = 0; i < NFF; ++i) {
    a[i] += b[i];
  }
}

void VectorDivide(vector<double> &a, double z) { 
   for(unsigned int i = 0; i < N_NEURONS; ++i) {
     a[i] /= z;
   }
}

void VectorDivideFF(vector<double> &a, double z) { 
   for(unsigned int i = 0; i < NFF; ++i) {
     a[i] /= z;
   }
}


void VectorCopy(vector<double> &a, vector<double> &b) {
  // copy elements of a to b
   for(unsigned int i = 0; i < N_NEURONS; ++i) {
     b[i] = a[i];
   }
}


double PopAvg(vector<double> &a, unsigned int start, unsigned int end) {
  double result = 0.0;
  for(unsigned int i = start; i < end; ++i) {
    result += a[i];
  }
  return result / (double)(end - start);
}

double AddElements(vector<double> &a, unsigned int start, unsigned int end) {
  // add elements of a[start] to a[end]
  double result = 0.0;
  for(unsigned int i = start; i < end; ++i) {
    result += a[i];
  }
  return result;
}

double Heavside(double input) {
  if(input > 0) {
    return 1.0;
  }
  else {
    return 0.0;
  }
}

void Create_Dir(string dir) {
  string mkdirp = "mkdir -p " ;
  mkdirp += dir;
  const char * cmd = mkdirp.c_str();
  const int dir_err = system(cmd);
  if(-1 == dir_err) {
    cout << "error creating directories" << endl ;
  }
  cout << "Created directory : " ;
  cout << mkdirp << endl ;
}

void LoadSparseConMat() {
  FILE *fpSparseConVec, *fpIdxVec, *fpNpostNeurons;
  fpSparseConVec = fopen("sparseConVec.dat", "rb");
  fpIdxVec = fopen("idxVec.dat", "rb");
  fpNpostNeurons = fopen("nPostNeurons.dat", "rb");
  //  int ;
  unsigned long int nConnections = 0, dummy = 0;
  printf("%p %p %p\n", fpIdxVec, fpNpostNeurons, fpSparseConVec);
  dummy = fread(nPostNeurons, sizeof(*nPostNeurons), N_NEURONS, fpNpostNeurons);
  fclose(fpNpostNeurons);
  for(unsigned int i = 0; i < N_NEURONS; ++i) {
    nConnections += nPostNeurons[i];
  }
  printf("#Post read\n #rec cons = %lu \n", nConnections);  
  sparseConVec = new unsigned int[nConnections] ;
  dummy = fread(sparseConVec, sizeof(*sparseConVec), nConnections, fpSparseConVec);
  if(dummy != nConnections) {
    printf("sparseConvec read error ? \n");
  }
  printf("#sparse cons = %lu \n\n", dummy);    
  printf("sparse vector read\n");  
  dummy = fread(idxVec, sizeof(*idxVec), N_NEURONS, fpIdxVec);
  printf("#idx vector read\n");    
  fclose(fpSparseConVec);
  fclose(fpIdxVec);
}

void LoadRewiredCon(unsigned long nElements) {
  unsigned long nElementsRead = 0;
  FILE *fpRewired;
  fpRewired = fopen("newPostNeurons.dat", "rb");
  nElementsRead = fread(IS_REWIRED_LINK, sizeof(*IS_REWIRED_LINK), nElements, fpRewired);
  if(nElements != nElementsRead) {
    printf("rewired read error ? \n");
  }
}

void LoadFFSparseConMat() {
  FILE *fpSparseConVecFF, *fpIdxVecFF, *fpNpostNeuronsFF;
  //  fpSparseConVecFF = fopen("sparseConVecFF.dat", "rb");
  perror("ohlala");
  fpIdxVecFF = fopen("idxVecFF.dat", "rb");
  fpNpostNeuronsFF = fopen("nPostNeuronsFF.dat", "rb");
  unsigned long int long nConnections = 0, dummy = 0;
  // printf("%p %p %p\n", fpIdxVecFF, fpNpostNeuronsFF, fpSparseConVecFF);
  dummy = fread(nPostNeuronsFF, sizeof(*nPostNeuronsFF), NFF, fpNpostNeuronsFF);
  fclose(fpNpostNeuronsFF);
  printf("#Post read\n");
  for(unsigned int i = 0; i < NFF; ++i) {
    nConnections += nPostNeuronsFF[i];
  }
  printf("#connections = %llu\n", nConnections);  
  sparseConVecFF = new unsigned int[nConnections];
  fpSparseConVecFF = fopen("sparseConVecFF.dat", "rb");
  perror("lalala");
  dummy = fread(sparseConVecFF, sizeof(*sparseConVecFF), nConnections, fpSparseConVecFF);
  perror("ooplala");
  //    dummy = fread(sparseConVecFF, sizeof(*sparseConVecFF), 1, fpSparseConVecFF);
  printf("sparseConvec read: %llu %llu \n", dummy, nConnections);
  if(dummy != nConnections) {
    printf("sparseConvec read error ? %llu %llu \n", dummy, nConnections);
  }
  printf("sparse vector read\n");  
  dummy = fread(idxVecFF, sizeof(*idxVecFF), NFF, fpIdxVecFF);
  printf("#idx vector read\n");    
  fclose(fpSparseConVecFF);
  fclose(fpIdxVecFF);
}

void RunSim() {
  double dt, probUpdateFF, probUpdateE, probUpdateI, uNet, spinOld = 0, spinOldFF = 0;
  // uExternalE, uExternalI;
  unsigned long int nSteps, i, nLastSteps, nInitialSteps, intervalLen = (NE + NI) / 100;
  unsigned int updateNeuronIdx = 0, chunkIdx = 0;
  int updatePop = 0; // 0: FF, 1: E, 2: I
  double runningFre = 0, runningFri = 0, runningFrFF = 0;
  // FILE *fpMeanRates;
  vector<double> spins(N_NEURONS);
  vector<double> spinsFF(NFF);  
  vector<double> spkTimes;
  vector<unsigned long int> spkNeuronIdx;
  vector<double> firingRates(N_NEURONS);
  vector<double> firingRatesFF(NFF);  
  vector<double> frLast(N_NEURONS);
  vector<double> totalInput(N_NEURONS);
  vector<double> FFInput(N_NEURONS);  
  vector<double> netInputVec(N_NEURONS);
  vector<double> firingRatesChk(N_NEURONS);
  vector<double> firingRatesChkTMP(N_NEURONS);
  vector<double> firingRatesAtT(N_NEURONS);  
  vector<double> ratesAtInterval(N_NEURONS);
  double popME1, popME1Phase;  // popMI1, popMI1Phase,
  double phiExtOld = phi_ext;

  std::string m1FileName;     
  m1FileName = "MI1_inst_theta" + std::to_string(phi_ext * 180 / M_PI) + "_tr" + std::to_string(trialNumber) + ".txt";
  FILE *fpInstM1 = fopen(m1FileName.c_str(), "w");  
  
  dt = (TAU_FF * TAU_E * TAU_I) / (NE * TAU_FF * TAU_I + NI * TAU_FF * TAU_E + NFF * TAU_E * TAU_I);
  nSteps = (unsigned long)((tStop / dt) + (T_TRANSIENT / dt));
  nInitialSteps = (T_TRANSIENT / dt);
  probUpdateFF = dt * (double)NFF / TAU_FF;  
  probUpdateE = dt * (double)NE / TAU_E;
  probUpdateI = dt * (double)NI / TAU_I;
  
  // uExternalE = sqrt((double) K) * JE0 * m0;
  // uExternalI = sqrt((double) K) * JI0 * m0;
  // printf("%f\n", uExternalE);
  // printf("%f\n", uExternalI);


  if(IF_STEP_PHI0) {
    // used for computing m_E^1(t)
    intervalLen = (unsigned long)floor((nSteps - nInitialSteps) / 80.0);
  }
  
  printf("dt = %f, #steps = %lu, T_STOP = %d\n", dt, nSteps, T_STOP);
  printf("prob updateFF = %f, prob update E = %f, prob update I = %f\n", probUpdateFF, probUpdateE, probUpdateI);
  nLastSteps = nSteps - (unsigned long int )((float)T_TRANSIENT / dt);
  printf("\n");
  ProgressBar(0.0, runningFre, runningFri, runningFrFF);
  FILE *fpInstRates = fopen("fr_inst.txt", "w");

  std::random_device rdFF;  
  std::default_random_engine generatorFF(rdFF());
  std::discrete_distribution<int> MultinomialDistr {probUpdateFF, probUpdateE, probUpdateI};

  int nPhisInSteps = 2;
  for(i = 0; i < nSteps; i++) {
    if(IF_STEP_PHI0) {
      
      if((i > 0) && i % (unsigned long)floor((nSteps - nInitialSteps) / (double)nPhisInSteps) == 0) {
	// change stimulus after half time
	phi_ext = phi_ext +  0.5 * M_PI / (double)nPhisInSteps;
	// printf("hello %f\n", phi_ext * 180  / M_PI);      
      }
    }
    if(i > 0 && i % (nSteps / 100) == 0) {
      runningFrFF = PopAvg(firingRatesFF, 0, NFF) / (double)(i + 1);      
      runningFre = PopAvg(firingRates, 0, NE) / (double)(i + 1);
      runningFri = PopAvg(firingRates, NE, N_NEURONS) / (double)(i + 1);
      fprintf(fpInstRates, "%f;%f\n", runningFre, runningFri);
      fflush(fpInstRates);
      ProgressBar((float)i / nSteps, runningFre, runningFri, runningFrFF);    
    }

    uNet = 0.0;
    updatePop = MultinomialDistr(generatorFF);

    if(updatePop == 0) {
      updateNeuronIdx = RandFFNeuron();
      spinOldFF = spinsFF[updateNeuronIdx];
      spinsFF[updateNeuronIdx] = 0;
      if(UniformRand() <= FFTuningCurve(updateNeuronIdx, phi_ext)) {
        // FFTuning returons prob of state_i = 1	
	spinsFF[updateNeuronIdx] = 1;
      }

//////////////////////////////////////////////////
      // FF input
      // if(i >= nLastSteps) { 
      // 	unsigned int tmpIdx1, cntr1, iidxx;
      // 	cntr = 0;
      // 	for(iidxx = 0; iidxx < NFF; iidxx++) {
      // 	  if(spinsFF[iidxx]) {
      // 	    tmpIdx1 = idxVecFF[iidxx];
      // 	    while(cntr < nPostNeuronsFF[iidxx]) {
      // 	      unsigned int kk = sparseConVecFF[tmpIdx1 + cntr];
      // 	      cntr1 += 1;
      // 	      if(kk < NE) {
      // 		FFInput[kk] += JE0_K;
      // 	      }
      // 	      else {
      // 		FFInput[kk] += JI0_K;
      // 	      }
      // 	    }
      // 	  }
      // 	}
      // }
//////////////////////////////////////////////////      

      
      
      if(spinOldFF == 0 && spinsFF[updateNeuronIdx] == 1) {
	unsigned int tmpIdx, cntr;
	cntr = 0;
	tmpIdx = idxVecFF[updateNeuronIdx];
	while(cntr < nPostNeuronsFF[updateNeuronIdx]) {
	  unsigned int kk = sparseConVecFF[tmpIdx + cntr];
	  cntr += 1;
	  if(kk < NE) {
	    netInputVec[kk] += JE0_K;
	    // if(i >= nLastSteps) { FFInput[kk] += JE0_K;}
	  }
	  else {
	    netInputVec[kk] += JI0_K;
	  }
	}
      }
      else if(spinOldFF == 1 && spinsFF[updateNeuronIdx] == 0) {
	unsigned int tmpIdx, cntr = 0;
	cntr = 0;      
	tmpIdx = idxVecFF[updateNeuronIdx];
	while(cntr < nPostNeuronsFF[updateNeuronIdx]) {
	  unsigned int kk = sparseConVecFF[tmpIdx + cntr];
	  cntr += 1;
	  if(kk < NE) {
	    netInputVec[kk] -= JE0_K;
	    //if(i >= nLastSteps) { FFInput[kk] -= JE0_K;}
	  }
	  else {
	    netInputVec[kk] -= JI0_K;
	  }
	}
      }
    }
    else {
      if(updatePop == 1) {
	updateNeuronIdx = RandENeuron();
	uNet = netInputVec[updateNeuronIdx] - THRESHOLD_E;
      }
      else if(updatePop == 2)  {
	updateNeuronIdx = RandINeuron();
	uNet = netInputVec[updateNeuronIdx] - THRESHOLD_I;     
      }
      spinOld = spins[updateNeuronIdx];    
      spins[updateNeuronIdx] = Heavside(uNet);
      if(spinOld == 0 && spins[updateNeuronIdx] == 1) {
	unsigned int tmpIdx, cntr;
	cntr = 0;
	tmpIdx = idxVec[updateNeuronIdx];
	while(cntr < nPostNeurons[updateNeuronIdx]) {
	  unsigned int kk = sparseConVec[tmpIdx + cntr];
	  cntr += 1;
	  if(updateNeuronIdx < NE)  {
	    unsigned int IS_STRENGTHENED = 0;
	    if(IF_LOADREWIREDCON) {
	      IS_STRENGTHENED = IS_REWIRED_LINK[tmpIdx + cntr - 1];
	    }
	    if(kk < NE) {
	      if(IS_STRENGTHENED) {
		netInputVec[kk] += (rewiredEEWeight * JEE_K);
		printf("rewired synapse!\n");
	      }
	      else {
		netInputVec[kk] += JEE_K;
	      }
	    }
	    else {
	      netInputVec[kk] += JIE_K;
	    }
	  }
	  else {
	    if(kk < NE) {
	      netInputVec[kk] += JEI_K;
	    }
	    else {
	      netInputVec[kk] += JII_K;
	    }
	  }
	}
      }
      else if(spinOld == 1 && spins[updateNeuronIdx] == 0) {
	// unsigned int tmpIdx;
	// tmpIdx = idxVec[updateNeuronIdx];
	// for(unsigned int kk = sparseConVec[tmpIdx]; kk < sparseConVec[tmpIdx + nPostNeurons[updateNeuronIdx]]; kk++) {
	unsigned int tmpIdx, cntr = 0;
	cntr = 0;      
	tmpIdx = idxVec[updateNeuronIdx];
	while(cntr < nPostNeurons[updateNeuronIdx]) {
	  unsigned int kk = sparseConVec[tmpIdx + cntr];
	  cntr += 1;
	  if(updateNeuronIdx < NE)  {
	    unsigned int IS_STRENGTHENED = 0;
	    if(IF_LOADREWIREDCON) {
	      IS_STRENGTHENED = IS_REWIRED_LINK[tmpIdx + cntr - 1];
	    }
	    if(kk < NE) {
	      if(IS_STRENGTHENED) {
		netInputVec[kk] -= (rewiredEEWeight * JEE_K);
	      }
	      else {
		netInputVec[kk] -= JEE_K;
	      }
	    }
	    else {
	      netInputVec[kk] -= JIE_K;
	    }
	  }
	  else {
	    if(kk < NE) {
	      netInputVec[kk] -= JEI_K;
	    }
	    else {
	      netInputVec[kk] -= JII_K;
	    }
	  }
	}
      }
    }

    VectorSum(ratesAtInterval, spins);           
    if(i > 0 && (i % intervalLen == 0)){
      VectorDivide(ratesAtInterval, (double)intervalLen);
      M1ComponentE(ratesAtInterval, NE, &popME1, &popME1Phase);
      	 
      if(IF_STEP_PHI0) {
	fprintf(fpInstM1, "%f;%f;%f\n", popME1, popME1Phase, phi_ext);
        fflush(fpInstM1);
      }
      else {
	fprintf(fpInstM1, "%f;%f\n", popME1, popME1Phase);
	fflush(fpInstM1);
      }
      for(unsigned int ijk = 0; ijk < N_NEURONS; ijk++) {
	ratesAtInterval[ijk] = 0;
      }
    }

    
    if(spinOld == 0 && spins[updateNeuronIdx] == 1) {
      spkTimes.push_back(i * dt);
      spkNeuronIdx.push_back(updateNeuronIdx);
    }
    
   VectorSum(firingRates, spins);
   VectorSumFF(firingRatesFF, spinsFF);   
   if(i >= nLastSteps) {
     VectorSum(frLast, spins);
     VectorSum(totalInput, netInputVec);
   }

   //   Store Chunks 
   if(i >= nInitialSteps) {
     VectorSum(firingRatesChk, spins);
     if(i > nInitialSteps) {
       if(!((i - nInitialSteps) %  (unsigned int)((nSteps - nInitialSteps - 1) / (double)N_SEGMENTS))) {
	 double chunckLength = (double)((nSteps - nInitialSteps - 1) / (double)N_SEGMENTS);
	 // if(!((i - nInitialSteps) %  (unsigned int)((nSteps - nInitialSteps - 1) / 2))) {	 
	 printf("\n chk i = %lu \n", i - nInitialSteps);
	 VectorCopy(firingRatesChk, firingRatesChkTMP);
	 // VectorDivide(firingRatesChkTMP, (i - nInitialSteps));
	 VectorDivide(firingRatesChkTMP, chunckLength);
	 std::string txtFileName = "meanrates_theta" + std::to_string(phi_ext * 180 / M_PI) + "_tr" + std::to_string(trialNumber) + "_chnk" + std::to_string(chunkIdx) + ".txt";
	 chunkIdx += 1;       
	 FILE *fpRates = fopen(txtFileName.c_str(), "w");
	 for(unsigned int ii = 0; ii < N_NEURONS; ii++) {
	   fprintf(fpRates, "%f\n", firingRatesChkTMP[ii]);
	 }
	 fclose(fpRates);
	 for(unsigned int ijk = 0; ijk < N_NEURONS; ijk++) {
	   firingRatesChk[ijk] = 0;
	 }
	 
       }
     }
   }
  }
    
  VectorDivide(firingRates, nSteps);
  VectorDivideFF(firingRatesFF, nSteps);  
  VectorDivide(frLast, (nSteps - nLastSteps));
  VectorDivide(totalInput, (nSteps - nLastSteps));
  // VectorDivide(FFInput, (nSteps - nLastSteps));    
  fclose(fpInstRates);
  fclose(fpInstM1);

  unsigned int tmpIdx1, cntr1, iidxx;

  for(iidxx = 0; iidxx < NFF; iidxx++) {
      tmpIdx1 = idxVecFF[iidxx];
      cntr1 = 0;      
      while(cntr1 < nPostNeuronsFF[iidxx]) {
	unsigned int kk = sparseConVecFF[tmpIdx1 + cntr1];
	cntr1 += 1;
	if(kk < NE) {
	  FFInput[kk] += JE0_K * firingRatesFF[iidxx];
	}
	else {
	  FFInput[kk] += JI0_K * firingRatesFF[iidxx];
	}
      }
  }

  
  std::string txtFileName = "meanrates_theta" + std::to_string(phiExtOld * 180 / M_PI) + "_tr" + std::to_string(trialNumber) + ".txt";  
  FILE *fpRates = fopen(txtFileName.c_str(), "w");
  for(unsigned int ii = 0; ii < N_NEURONS; ii++) {
    fprintf(fpRates, "%f\n", firingRates[ii]);
  }
  fclose(fpRates);

  txtFileName = "meanrates_theta" + std::to_string(phiExtOld * 180 / M_PI) + "_tr" + std::to_string(trialNumber) + "_last.txt";  
  fpRates = fopen(txtFileName.c_str(), "w");
  for(unsigned int ii = 0; ii < N_NEURONS; ii++) {
    fprintf(fpRates, "%f\n", frLast[ii]);
  }
  fclose(fpRates);

  if(IF_SAVE_INPUT) {
    txtFileName = "meaninput_theta" + std::to_string(phiExtOld * 180 / M_PI) + "_tr" + std::to_string(trialNumber) + "_last.txt";  
    FILE *fpInputs = fopen(txtFileName.c_str(), "w");
    for(unsigned int ii = 0; ii < N_NEURONS; ii++) {
      fprintf(fpInputs, "%f\n", totalInput[ii]);
    }
    fclose(fpInputs);
    //
    txtFileName = "meanFFinput_theta" + std::to_string(phiExtOld * 180 / M_PI) + "_tr" + std::to_string(trialNumber) + "_last.txt";  
    FILE *fpFFInputs = fopen(txtFileName.c_str(), "w");
    for(unsigned int ii = 0; ii < NE+NI; ii++) {
      // fprintf(fpFFInputs, "%f\n", firingRatesFF[ii]);
      fprintf(fpFFInputs, "%f\n", FFInput[ii]);      
    }
    fclose(fpFFInputs);
  }

  
  if(IF_SAVE_SPKS) {
      txtFileName = "spktimes_theta" + std::to_string(phiExtOld * 180 / M_PI) + "_tr" + std::to_string(trialNumber) + ".txt";
      FILE *fpSpks = fopen(txtFileName.c_str(), "w");
      for(unsigned long long int ii = 0; ii < spkNeuronIdx.size(); ii++) {
	fprintf(fpSpks, "%lu;%f\n", spkNeuronIdx[ii], spkTimes[ii]);
      }
      fclose(fpSpks);
  }

  spins.clear();
  spkTimes.clear();
  spkNeuronIdx.clear();
  firingRates.clear();
  totalInput.clear();
  FFInput.clear();
  netInputVec.clear();
  firingRatesChk.clear();
  firingRatesAtT.clear();  
  ratesAtInterval.clear();
  printf("\nsee you later, alligator\n");
}
  
int main(int argc, char *argv[]) {
  printf("#args = %d\n", argc);
  if(argc > 1) {
    m0_ext = atof(argv[1]);
  }
  if(argc > 2) {
    printf("hello mExtOn1\n");
    m1_ext = atof(argv[2]);
  }
  if(argc > 3) {
    recModulationEE = atof(argv[3]); // parameter p
    recModulationIE = 0; //-1.0 * recModulationEE;
    recModulationEI = 0; //-1.0 * recModulationEE;
  }
  if(argc > 4) {
    ffModulation = atof(argv[4]); // parameter gamma
  }
  if(argc > 5) {
    phi_ext = M_PI * atof(argv[5]) / 180.0; // parameter external orientation
  }
  if(argc > 6) {
    trialNumber = atof(argv[6]); // parameter gamma
  }

  if(argc > 7) {
    rewiredEEWeight = atof(argv[7]); // parameter strengthened weight prefactor
    IF_LOADREWIREDCON = 1;
  }

  
  tStop = T_STOP;

  cout << "NE = " << NE << " NI = " << NI << " NFF = " << NFF << " K = " << K << " p = " << recModulationEE << " m0 = " << m0_ext << " m0_One = " << m1_ext << endl;
  cout << "gamma = " << ffModulation << " KFF = " << cFF * K << " Phi_Ext = " << phi_ext * 180.0 / M_PI << endl;
  cout << "TAU_E = " << TAU_E << " Tau_I = " << TAU_I << endl;
  cout << "Trial# = " << trialNumber << endl;
  //  sprintf(folderName, "N%uK%um0%dpgamma%dT%d", N_NEURONS, (int)(m0 * 1e3), K, recModulation, ffModulation, (int)(tStop * 1e-3));
  // folderName = "./data/N" + std::to_string(N_NEURONS) + "K" + std::to_string(K) + "m0" + std::to_string((int)(m0 * 1e3)) + "p" + std::to_string((unsigned int)(10 * recModulation)) + "gamma" + std::to_string((unsigned int)(ffModulation)) + std::to_string((int)(tStop * 1e-3));

  folderName = "./data/N" + std::to_string(N_NEURONS) + "K" + std::to_string(K) + "m0" + std::to_string((int)(m0 * 1e3)) + "p" + std::to_string((unsigned int)(10 * recModulationEE)) + "gamma0" + std::to_string((int)(tStop * 1e-3));  

  // string mkdirp = "mkdir -p " ;
  // Create_Dir(folderName);
  
  nPostNeurons = new unsigned int[N_NEURONS];
  idxVec = new unsigned int[N_NEURONS];

  nPostNeuronsFF = new unsigned int[NFF];
  idxVecFF = new unsigned int[NFF];



  if(trialNumber == 0) {
    clock_t timeStartCM = clock(); 
    GenFFConMat();    
    GenConMat();
    clock_t timeStopCM = clock();
    double elapsedTimeCM = (double)(timeStopCM - timeStartCM) / CLOCKS_PER_SEC;  
    cout << "\n connection gen, elapsed time= " << elapsedTimeCM << "s, or " << elapsedTimeCM / 60.0 << "min" << endl;
  }
  else {
    printf("loading FF Sparse matrix\n");      
    LoadFFSparseConMat();
    printf("loading Sparse matrix\n");  
    LoadSparseConMat();

    unsigned long nEE_IEConnections = 0;
    for(unsigned i = 0; i < NE; i++) {
      nEE_IEConnections += nPostNeurons[i];
    }
    IS_REWIRED_LINK = new unsigned int[nEE_IEConnections];
    if(IF_LOADREWIREDCON) {
      LoadRewiredCon(nEE_IEConnections);
    }
    else {
      for(unsigned long i = 0; i < nEE_IEConnections; i++) {
	IS_REWIRED_LINK[i] = 0;
      }
    }
  }

  
  clock_t timeStart = clock(); 
  RunSim();
  clock_t timeStop = clock();
  double elapsedTime = (double)(timeStop - timeStart) / CLOCKS_PER_SEC;
  cout << "\n elapsed time = " << elapsedTime << "s, or " << elapsedTime / 60.0 << "min" << endl; 
  FILE *fpet = fopen("elapsedTime.txt", "a");
  fprintf(fpet, "%llu %f\n", nSteps, elapsedTime);
  fclose(fpet);

  // delete [] conMat;
  delete [] nPostNeurons;
  delete [] idxVec;
  delete [] sparseConVec;
  delete [] nPostNeuronsFF;
  delete [] idxVecFF;
  delete [] sparseConVecFF;
  delete [] IS_REWIRED_LINK;
  return 0; //EXIT_SUCCESS;
}

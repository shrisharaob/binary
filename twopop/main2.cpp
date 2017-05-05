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


double ConProb(double phiI, double phiJ, unsigned int N) {
  double out = ((double)K / (double)N) * (1 + 2.0 * (recModulation / SQRT_K) * cos(2.0 * (phiI - phiJ)));
  if((out < 0) | (out > 1)) {
    cout << "connection probability not in range [0, 1]!" << endl;
    exit(1);
  }
  return out;
}

double ConProbFF(double phiI, double phiJ, unsigned int N) {
  double out = ((double)K * cFF / (double)N) * (1 + 2.0 * (ffModulation / SQRT_K) * cos(2.0 * (phiI - phiJ)));
  if((out < 0) | (out > 1)) {
    cout << "connection probability not in range [0, 1]!" << endl;
    exit(1);
  }
  return out;
}

double FFTuningCurve(unsigned int i) {
  double out = m0_ext + m1_ext * cos(2 * (phi_ext - (double)i * M_PI / (double)NFF));
  if((out < 0) | (out > 1)) {
    cout << "external rate [0, 1]!" << endl;
    exit(1);
  }
  out = m0_ext;
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
  
  FILE *fp2 = NULL;
  fp2 = fopen("prob.csv", "w");
  for(i = 0; i < NE; i++) {
    // printf("%llu:\n", i); //    %f\n", i, shiftedVec[i] / (double)NE);    
    //fprintf(fp, "%f\n", shiftedVec[i] / (double)NE);
    fprintf(fp2, "%f\n", ConProb(0.0, (double)i * M_PI / (double)NE, NE));
  }
  fclose(fp2);
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
  unsigned long int nElementsWritten, dummy;
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
  printf("testing sparsevecff read\n");
  FILE *fpSparseConVecFF = fopen("sparseConVecFF.dat", "rb");
  dummy = fread(sparseConVecFF, sizeof(*sparseConVecFF), nConnections, fpSparseConVecFF);
    printf("sparseConvec read: %lu %llu \n", dummy, nConnections);
  if(dummy != nConnections) {
    printf("sparseConvec read error ? %lu %llu \n", dummy, nConnections);
  fclose(fpSparseConVecFF);
  }
}
    
void GenConMat() {
  std::random_device rd;
  std::default_random_engine gen(rd());
  std::uniform_real_distribution<double> UniformRand(0.0, 1.0);
  //  int dummy;
  unsigned long long int nConnections = 0;
  cout << "generating conmat" << endl;
  unsigned int *conMat = new unsigned int [(unsigned long int)N_NEURONS * N_NEURONS];
  for (unsigned long int i = 0; i < N_NEURONS; i++)  {
    for (unsigned long int j = 0; j < N_NEURONS; j++)  {
      // i --> j
      if(i < NE && j < NE) { //E-to-E
	if(UniformRand(gen) <= ConProb(i * M_PI / (double)NE, j * M_PI / (double)NE, NE)) {
	  conMat[i + N_NEURONS * j] = 1;
	  nConnections += 1;
	}
      }
      if(i < NE && j >= NE) { //E-to-I
	if(UniformRand(gen) <= (double)K / (double)NE) {
	  conMat[i + N_NEURONS * j] = 1;
	  nConnections += 1;
	}
      }
      if(i >= NE) { //I-to-E and i-to-i
	if(UniformRand(gen) <= (double)K / (double)NI) {
	  conMat[i + N_NEURONS * j] = 1;
	  nConnections += 1;
	}
      }
    }
  }

  cout << "done" << endl;
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

  FILE *fpSparseConVec, *fpIdxVec, *fpNpostNeurons;
  unsigned long int nElementsWritten;
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
  printf("#Post read\n");
  for(unsigned int i = 0; i < N_NEURONS; ++i) {
    nConnections += nPostNeurons[i];
  }
  sparseConVec = new unsigned int[nConnections] ;
  dummy = fread(sparseConVec, sizeof(*sparseConVec), nConnections, fpSparseConVec);
  if(dummy != nConnections) {
    printf("sparseConvec read error ? \n");

  }
  printf("sparse vector read\n");  
  dummy = fread(idxVec, sizeof(*idxVec), N_NEURONS, fpIdxVec);
  printf("#idx vector read\n");    
  fclose(fpSparseConVec);
  fclose(fpIdxVec);
}


void LoadFFSparseConMat() {
  FILE *fpSparseConVecFF, *fpIdxVecFF, *fpNpostNeuronsFF;
  //  fpSparseConVecFF = fopen("sparseConVecFF.dat", "rb");
  perror("ohlala");
  fpIdxVecFF = fopen("idxVecFF.dat", "rb");
  fpNpostNeuronsFF = fopen("nPostNeuronsFF.dat", "rb");
  unsigned long int long nConnections = 0, dummy = 0;
  printf("%p %p %p\n", fpIdxVecFF, fpNpostNeuronsFF, fpSparseConVecFF);
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
  double dt, probUpdateFF, probUpdateE, probUpdateI, uNet, spinOld, spinOldFF;
  // uExternalE, uExternalI;
  unsigned long int nSteps, i, nLastSteps;
  unsigned int updateNeuronIdx;
  int updatePop = 0; // 0: FF, 1: E, 2: I
  double runningFre = 0, runningFri = 0, runningFrFF = 0;
  // FILE *fpMeanRates;
  vector<double> spins(N_NEURONS);
  vector<double> spinsFF(NFF);  
  vector<double> spkTimes;
  vector<double> spkNeuronIdx;
  vector<double> firingRates(N_NEURONS);
  vector<double> firingRatesFF(NFF);  
  vector<double> frLast(N_NEURONS);  
  vector<double> netInputVec(N_NEURONS);
  dt = (TAU_FF * TAU_E * TAU_I) / (NE * TAU_FF * TAU_I + NI * TAU_FF * TAU_E + NFF * TAU_E * TAU_I);
  nSteps = (unsigned long)(tStop / dt);
  probUpdateFF = dt * (double)NFF / TAU_FF;  
  probUpdateE = dt * (double)NE / TAU_E;
  probUpdateI = dt * (double)NI / TAU_I;
  // uExternalE = sqrt((double) K) * JE0 * m0;
  // uExternalI = sqrt((double) K) * JI0 * m0;
  // printf("%f\n", uExternalE);
  // printf("%f\n", uExternalI);    
  printf("dt = %f, #steps = %lu, T_STOP = %d\n", dt, nSteps, T_STOP);
  printf("prob updateFF = %f, prob update E = %f, prob update I = %f\n", probUpdateFF, probUpdateE, probUpdateI);
  nLastSteps = nSteps - (unsigned long int )((float)T_TRANSIENT / dt);
  printf("\n");
  ProgressBar(0.0, runningFre, runningFri, runningFrFF);
  FILE *fpInstRates = fopen("fr_inst.txt", "w");

  std::random_device rdFF;  
  std::default_random_engine generatorFF(rdFF());
  std::discrete_distribution<int> MultinomialDistr {probUpdateFF, probUpdateE, probUpdateI};
  
  for(i = 0; i < nSteps; i++) {
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
      spinsFF[updateNeuronIdx] = UniformRand() <= FFTuningCurve(updateNeuronIdx); // FFTuning returons prob of state_i = 1
      if(spinOldFF == 0 && spinsFF[updateNeuronIdx] == 1) {
	unsigned int tmpIdx, cntr;
	cntr = 0;
	tmpIdx = idxVecFF[updateNeuronIdx];
	while(cntr < nPostNeuronsFF[updateNeuronIdx]) {
	  unsigned int kk = sparseConVecFF[tmpIdx + cntr];
	  cntr += 1;
	  if(kk < NE) {
	    netInputVec[kk] += JE0_K;
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
	    if(kk < NE) {
	      netInputVec[kk] += JEE_K;
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
	    if(kk < NE) {
	      netInputVec[kk] -= JEE_K;
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
    if(spinOld == 0 && spins[updateNeuronIdx] == 1) {
      spkTimes.push_back(i * dt);
      spkNeuronIdx.push_back(updateNeuronIdx);
    }
   VectorSum(firingRates, spins);
   VectorSumFF(firingRatesFF, spinsFF);   
   if(i >= nLastSteps) {
     VectorSum(frLast, spins);
   }
  }
  VectorDivide(firingRates, nSteps);
  VectorDivide(firingRatesFF, nSteps);  
  VectorDivide(frLast, (nSteps - nLastSteps));  
  fclose(fpInstRates);
    
  std::string txtFileName = "meanrates_theta" + std::to_string(phi_ext * 180 / M_PI) + "_tr" + std::to_string(trialNumber) + ".txt";  
  FILE *fpRates = fopen(txtFileName.c_str(), "w");
  for(unsigned int ii = 0; ii < N_NEURONS; ii++) {
    fprintf(fpRates, "%f\n", firingRates[ii]);
  }
  fclose(fpRates);
  txtFileName = "meanrates_theta" + std::to_string(phi_ext * 180 / M_PI) + "_last.txt";  
  fpRates = fopen(txtFileName.c_str(), "w");
  for(unsigned int ii = 0; ii < N_NEURONS; ii++) {
    fprintf(fpRates, "%f\n", frLast[ii]);
  }
  fclose(fpRates); 
  
  spins.clear();
  spkTimes.clear();
  spkNeuronIdx.clear();
  firingRates.clear();
  netInputVec.clear();
  printf("\nsee you later, alligator\n");
}
  
int main(int argc, char *argv[]) {
  if(argc > 1) {
    m0_ext = atof(argv[1]);
  }
  m1_ext = 0.1 * m0_ext;
  if(argc > 2) {
    recModulation = atof(argv[2]); // parameter p
  }
  if(argc > 3) {
    ffModulation = atof(argv[3]); // parameter gamma
  }
  if(argc > 4) {
    phi_ext = M_PI * atof(argv[4]) / 180.0; // parameter external orientation
  }
  if(argc > 5) {
    trialNumber = atof(argv[5]); // parameter gamma
  }

  tStop = 1000; //T_STOP;

  cout << "NE = " << NE << " NI = " << NI << " NFF = " << NFF << " K = " << K << " p = " << recModulation << " m0 = " << m0 << endl ;
  cout << "TAU_E = " << TAU_E << " Tau_I = " << TAU_I << endl;
  //  sprintf(folderName, "N%uK%um0%dpgamma%dT%d", N_NEURONS, (int)(m0 * 1e3), K, recModulation, ffModulation, (int)(tStop * 1e-3));

  // folderName = "./data/N" + std::to_string(N_NEURONS) + "K" + std::to_string(K) + "m0" + std::to_string((int)(m0 * 1e3)) + "p" + std::to_string((unsigned int)(10 * recModulation)) + "gamma" + std::to_string((unsigned int)(ffModulation)) + std::to_string((int)(tStop * 1e-3));

  folderName = "./data/N" + std::to_string(N_NEURONS) + "K" + std::to_string(K) + "m0" + std::to_string((int)(m0 * 1e3)) + "p" + std::to_string((unsigned int)(10 * recModulation)) + "gamma0" + std::to_string((int)(tStop * 1e-3));  

  // string mkdirp = "mkdir -p " ;
  // Create_Dir(folderName);
  
  nPostNeurons = new unsigned int[N_NEURONS];
  idxVec = new unsigned int[N_NEURONS];

  nPostNeuronsFF = new unsigned int[NFF];
  idxVecFF = new unsigned int[NFF];



  if(phi_ext == 0) {
    GenFFConMat();    
    clock_t timeStartCM = clock(); 
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
  return 0; //EXIT_SUCCESS;
}

#include <stdio.h>
#include <math.h>
// #include <boost/random/uniform_01.hpp>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <random>
#include <time.h>
#include <string>
#include <cstring>
#include <pthread.h>
//#include <complex>
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

// double VectorDotProd(double a, vector<double> &b, unsigned n) { //, unsigned int n) {
//   for(unsigned int i = 0; i < n; ++i) {
//     c[i] = a[i] * b[i];
//   }
// }

void M1Component(vector<double> &x, unsigned int n, double* m1, double* phase) {
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


// double UniformRand() {
//   std::random_device rd2;
//   std::default_random_engine gen2(rd());
//   std::uniform_real_distribution<double> UniformRand2(0.0, 1.0);
//   return UniformRand()
//   // return (double)rand() / (double)RAND_MAX ;
// }

double UniformRand() {
  std::random_device rd;
  std::default_random_engine gen(rd());
  std::uniform_real_distribution<double> UniformRand_Aux(0.0, 1.0);
  return UniformRand_Aux(gen);
  // return (double)rand() / (double)RAND_MAX ;
}


void Dummy() {
  printf("me too\n");
}

std::default_random_engine generator;
std::uniform_int_distribution<int> IntDistr(0, NI);
unsigned int RandINeuron() { //unsigned int min, unsigned int max) {
  unsigned int result = NI;
  while(result == NI) {
    result = IntDistr(generator);
    
    // result = NE + (unsigned int) (rand() / (double) (RAND_MAX + 1UL) * (NI + 1));
  }
  return result;
}

void VectorCopy(vector<double> &a, vector<double> &b) { //, unsigned int n) {
  // copy elements of a to b
  for(unsigned int i = 0; i < N_NEURONS; ++i) {
    b[i] = a[i];
  }
}

void VectorSum(vector<double> &a, vector<double> &b) { //, unsigned int n) {
  // add elements of b to a
  for(unsigned int i = 0; i < N_NEURONS; ++i) {
    a[i] += b[i];
  }
}

void VectorDivide(vector<double> &a, double z) { //, unsigned int n) {
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

double ConProb(double phiI, double phiJ, unsigned int N) {
  double out = ((double)K / (double)N) * (1 - 2.0 * (recModulation / SQRT_K) * cos(2.0 * (phiI - phiJ)));
  if((out < 0) | (out > 1)) {
    cout << "connection probability not in range [0, 1]!" << endl;
    exit(1);
  }
  return out;
}



void GenSparseMat(unsigned int *conVec,  unsigned int rows, unsigned int clms, unsigned int* sparseVec, unsigned int* idxVec, unsigned int* nPostNeurons ) {
  /* generate sparse representation
     conVec       : input vector / flattened matrix 
     sparseVec    : sparse vector
     idxVec       : every element is the starting index in sparseVec for ith row in matrix conVec
     nPostNeurons : number of non-zero elements in ith row 
  */
  unsigned long long int i, j, counter = 0, nPost;
  for(i = 0; i < rows; ++i) {
    nPost = 0;
    for(j = 0; j < clms; ++j) {
      if(conVec[i + clms * j]) { /* i --> j  */
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
  FILE *fp = NULL;
  fp = fopen("check_prob.csv", "w");
  vector <unsigned int> shiftedVec(N_NEURONS);
  vector <unsigned int> shiftedVecBuff(N_NEURONS);  
  //  unsigned int ccount = 0;
  for(i = 0; i < N_NEURONS; ++i) {
    //    ccount = 0;
    for(j =0; j < N_NEURONS; ++j) {
      shiftedVecBuff[j] = conVec[j + i * N_NEURONS];
    }
    std::rotate(shiftedVecBuff.begin(), shiftedVecBuff.begin() + i, shiftedVecBuff.end());
    for(unsigned int l = 0; l < N_NEURONS; l++) {
      shiftedVec[l] += shiftedVecBuff[l];
    }
    // // for(i = 0; i < N_NEURONS; ++i) {
    // if(conVec[i + j * N_NEURONS]) {
    //   ccount += 1;
    // }
    // // }
    // fprintf(fp, "%u\n", ccount);
  }
  FILE *fp2 = NULL;
  fp2 = fopen("prob.csv", "w");
  for(i = 0; i < N_NEURONS; ++i) {
    fprintf(fp, "%f\n", shiftedVec[i] / (double)N_NEURONS);
    fprintf(fp2, "%f\n", ConProb(0.0, i * M_PI / (double)NI, NI));
  }
  fclose(fp); 
  fclose(fp2); 
  shiftedVecBuff.clear();
  shiftedVec.clear();
}

void GenConMat() {
  std::random_device rd;
  std::default_random_engine gen(rd());
  std::uniform_real_distribution<double> UniformRand(0.0, 1.0);
  //  int dummy;
  unsigned long long int nConnections = 0;
  unsigned int N;
  cout << "generating conmat" << endl;
  unsigned int *conMat = new unsigned int [(unsigned long int)N_NEURONS * N_NEURONS];
  for (unsigned long int i = 0; i < N_NEURONS; i++)  {
    for (unsigned long int j = 0; j < N_NEURONS; j++)  {
      N = NI;
      if(j < NE) {
	N = NE;
      }
      if(UniformRand(gen) <= ConProb(i * M_PI / (double)N, j * M_PI / (double)N, N)) {
	conMat[i + N_NEURONS * j] = 1;
	//conMat[j + N_NEURONS * i] = 1;	
	nConnections += 1;
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

  // Write sparse vector to file
  if(IF_SAVE_MAT) {
    FILE *fpSparseConVec, *fpIdxVec, *fpNpostNeurons;
    fpSparseConVec = fopen("sparseConVec.dat", "wb");
    unsigned long long int nElementsWritten;
    printf("done\n#connections = %llu\n", nConnections);
    printf("writing to file ... "); fflush(stdout);
    fpSparseConVec = fopen("sparseConVec.dat", "wb");
    nElementsWritten = fwrite(sparseConVec, sizeof(*sparseConVec), nConnections, fpSparseConVec);
    if(nElementsWritten != nConnections) {
      cout << "incorrect #of elemens written to sparvec file!" << endl;
      exit(1);
    }
    fclose(fpSparseConVec);
    fpIdxVec = fopen("idxVec.dat", "wb");
    fwrite(idxVec, sizeof(*idxVec), N_NEURONS,  fpIdxVec);
    fclose(fpIdxVec);
    fpNpostNeurons = fopen("nPostNeurons.dat", "wb");
    fwrite(nPostNeurons, sizeof(*nPostNeurons), N_NEURONS, fpNpostNeurons);
    fclose(fpNpostNeurons);
    printf("done\n");
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
  unsigned long long int dummy;
  printf("%p %p %p\n", fpIdxVec, fpNpostNeurons, fpSparseConVec);
  dummy = fread(nPostNeurons, sizeof(*nPostNeurons), N_NEURONS, fpNpostNeurons);
  fclose(fpNpostNeurons);
  printf("#Post read\n");
  unsigned long long int nConnections = 0;
  for(unsigned int i = 0; i < N_NEURONS; ++i) {
    // printf("%u", i);
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

void RunSim() {
  double uNet, spinOld, uExternalI;
  unsigned long long int i, nLastSteps = 2000 * NI, intervalLen = NI / 1000;
  unsigned int updateNeuronIdx;
  //  int updatePop;
  double runningFre = 0;
  double runningFri = 0;
  // FILE *fpMeanRates;
  vector<double> spins(N_NEURONS);
  vector<double> spkTimes;
  vector<double> spkNeuronIdx;
  vector<double> firingRates(N_NEURONS);
  vector<double> firingRatesAtT(N_NEURONS);
  vector<double> ratesAtInterval(N_NEURONS);  
  vector<double> frLast(N_NEURONS);  
  vector<double> netInputVec(N_NEURONS);
  nSteps = NI * 1000;
  uExternalI = sqrt((double) K) * JI0 * m0;
  printf("external input = %f\n", uExternalI);    
  printf("#Steps = %llu, T_STOP = %d, trial# = %d \n", nSteps, T_STOP, trialNumber);
  nLastSteps = nSteps - NI*2; 
  printf("\n");
  cout << "interval len: " << intervalLen << endl;
  std::string fileName; cout << "J_sqrtK = " << JII_K << endl;
  cout << ' ' << endl;
  ProgressBar(0.0, runningFre, runningFri);
  FILE *fpInstRates = fopen("fr_inst.txt", "w");
  std::string m1FileName;     
  m1FileName = "MI1_inst_tr" + std::to_string(trialNumber) + ".txt";
  FILE *fpInstM1 = fopen(m1FileName.c_str(), "w");  
  FILE *fpRatesAtInter2 = NULL; // *fpRatesAtT = NULL,
  double popMI1, popMI1Phase;

  std::string popAvgFn;     
  popAvgFn = "popAvg_tr" + std::to_string(trialNumber) + ".txt";
  FILE *fpSpkCnt = fopen(popAvgFn.c_str(), "w");
    
  // Random initialization
  // sleep(10);
  
  srand (time(NULL));
  for(unsigned tstl = 0; tstl < 10; tstl++) {
    cout << UniformRand() << endl;
  }
  //  exit(1);
  // double initialProb = UniformRand();
  // cout << "init prob = " << initialProb << endl;
  for(unsigned l = 0; l < N_NEURONS; l++) {
    frLast[l] = 0;
    spins[l] = 0;
    if(UniformRand() < 0.5) {
      spins[l] = 1;
    }
    unsigned int tmpIdxInit, cntrInit;
    if(spins[l]) {
      cntrInit = 0;
      tmpIdxInit = idxVec[l];
      while(cntrInit < nPostNeurons[l]) {
	unsigned int kkInit = sparseConVec[tmpIdxInit + cntrInit];
	cntrInit += 1;
	netInputVec[kkInit] += JII_K;	
      }
    }
  }

  double spkCntAvg = 0;
  spkCntAvg = PopAvg(spins, NE, N_NEURONS);
  fprintf(fpSpkCnt, "%f\n", spkCntAvg);
  cout << "Initial spk cnt / N = "  << spkCntAvg << endl;
  
  // Run 
  for(i = 0; i < nSteps; i++) { 
    if(i > NI*100 && i % (nSteps / 100) == 0) {
      runningFri = PopAvg(spins, NE, N_NEURONS);
      fprintf(fpInstRates, "%f;%f\n\n", runningFre, runningFri);
      fflush(fpInstRates);
      ProgressBar((float)i / nSteps, runningFre, runningFri);    
    }

    // if(i > 0 && i % (nSteps / 1000) == 0) {
    //   runningFri = PopAvg(spins, NE, N_NEURONS);
    //   M1Component(spins, NI, &popMI1, &popMI1Phase);
    //   fprintf(fpInstM1, "%f;%f\n\n", popMI1, popMI1Phase);
    //   fflush(fpInstM1);
    // }

    uNet = 0.0;
    updateNeuronIdx = RandINeuron();
    uNet = netInputVec[updateNeuronIdx] + uExternalI - THRESHOLD_I;
    spinOld = spins[updateNeuronIdx];    
    spins[updateNeuronIdx] = Heavside(uNet);
    // cout << "unet = " << uNet << " spin state:  " << "old= " << spinOld <<  " new= " << spins[updateNeuronIdx] << endl;    
    if(spinOld == 0 && spins[updateNeuronIdx] == 1) {
      unsigned int tmpIdx, cntr;
      cntr = 0;
      tmpIdx = idxVec[updateNeuronIdx];
      while(cntr < nPostNeurons[updateNeuronIdx]) {
	unsigned int kk = sparseConVec[tmpIdx + cntr];
	cntr += 1;
        netInputVec[kk] += JII_K;	
      }
    }
    else if(spinOld == 1 && spins[updateNeuronIdx] == 0) {
      // cout << "flip 1--> 0 " << endl;
      unsigned int tmpIdx, cntr = 0;
      cntr = 0;      
      tmpIdx = idxVec[updateNeuronIdx];
      while(cntr < nPostNeurons[updateNeuronIdx]) {
	unsigned int kk = sparseConVec[tmpIdx + cntr];
	cntr += 1;
        netInputVec[kk] -= JII_K;
      }
    }
    VectorSum(ratesAtInterval, spins);           
    if(i > 0 && (i % intervalLen == 0)){
      // std::string fileName2;     
      // fileName2 = "mean_rates_Inter" + std::to_string((unsigned int)(i / intervalLen)) + "E3_tr" + std::to_string(trialNumber) + ".txt";
      // fpRatesAtInter2 = fopen(fileName2.c_str(), "w");
      VectorDivide(ratesAtInterval, (double)intervalLen);
      M1Component(ratesAtInterval, NI, &popMI1, &popMI1Phase);
      fprintf(fpInstM1, "%f;%f\n\n", popMI1, popMI1Phase);
      fflush(fpInstM1);
      for(unsigned int ijk = 0; ijk < N_NEURONS; ijk++) {
	//	fprintf(fpRatesAtInter2, "%f\n", ratesAtInterval[ijk]);
	ratesAtInterval[ijk] = 0;
      }
      //      fclose(fpRatesAtInter2);
    }
    if(spinOld == 0 && spins[updateNeuronIdx] == 1) {
      spkTimes.push_back(i);
      spkNeuronIdx.push_back(updateNeuronIdx);
    }
    VectorSum(firingRates, spins);
    if(i >= (nSteps - nLastSteps)) {
      VectorSum(frLast, spins);
    }

    spkCntAvg = PopAvg(spins, NE, N_NEURONS);
    fprintf(fpSpkCnt, "%f\n", spkCntAvg);
    
  }
  cout << "\n " << endl;
  //  cout << "\n hello \n" << endl;
  VectorDivide(firingRates, nSteps);
  VectorDivide(frLast, nLastSteps);  
  fclose(fpInstRates);
  fclose(fpInstM1);
  fclose(fpSpkCnt);
  std::string txtFileName = "meanrates_theta0_tr"+ std::to_string(trialNumber) + ".txt";  
  FILE *fpRates = fopen(txtFileName.c_str(), "w");
  for(unsigned int ii = 0; ii < N_NEURONS; ii++) {
    fprintf(fpRates, "%f\n", firingRates[ii]);
  }
  fclose(fpRates);
  txtFileName = "meanrates_theta0_last.txt";  
  fpRates = fopen(txtFileName.c_str(), "w");
  for(unsigned int ii = 0; ii < N_NEURONS; ii++) {
    fprintf(fpRates, "%f\n", frLast[ii]);
  }
  fclose(fpRates);
  txtFileName = "spktimes_tr"+ std::to_string(trialNumber) + ".txt";  
  FILE *fpSpks = fopen(txtFileName.c_str(), "w");
  if(spkTimes.size() > 0) {
    for(unsigned int ii = 0; ii < N_NEURONS; ii++) {
      fprintf(fpSpks, "%f;%f\n", spkNeuronIdx[ii], spkTimes[ii]);
    }
  }
  fclose(fpSpks); 
  ratesAtInterval.clear();  
  spins.clear();
  spkTimes.clear();
  spkNeuronIdx.clear();
  firingRates.clear();
  firingRatesAtT.clear();
  netInputVec.clear();
  printf("\n");
}
  
int main(int argc, char *argv[]) {
  if(argc > 1) {
    m0 = atof(argv[1]);
  }
  if(argc > 2) {
    recModulation = atof(argv[2]); // parameter p
  }
  if(argc > 3) {
    ffModulation = atof(argv[3]); // parameter gamma
  }
  if(argc > 4) {
    trialNumber = atof(argv[4]); // file tag
  }
  if(argc > 5) {
    cout << "ahha!" << endl;
    IF_SAVE_MAT = atoi(argv[5]);
    cout << "but voila!" << endl;    
  }
  tStop = T_STOP;
  if(IF_SAVE_MAT) {
    IF_GEN_MAT = 1;
  }
  else {
    IF_GEN_MAT = IF_GEN_MAT_DEFAULT;
  }
  cout << "NE = " << NE << "NI = " << NI << " K = " << K << " p = " << recModulation << " m0 = " << m0 << endl ;
  cout << "TAU_E = " << TAU_E << " Tau_I = " << TAU_I << endl;
  folderName = "./data/N" + std::to_string(N_NEURONS) + "K" + std::to_string(K) + "m0" + std::to_string((int)(m0 * 1e3)) + "p" + std::to_string((unsigned int)(10 * recModulation)) + "gamma0" + std::to_string((int)(tStop * 1e-3));  
  // string mkdirp = "mkdir -p " ;
  // Create_Dir(folderName);
  nPostNeurons = new unsigned int[N_NEURONS];
  idxVec = new unsigned int[N_NEURONS];
  FILE *fpMatFlag = NULL;
  fpMatFlag = fopen("IF_MatReady.txt", "w");
  fprintf(fpMatFlag, "0\n");
  fclose(fpMatFlag);

  if(IF_GEN_MAT) {
    GenConMat();
  }
  else {
    printf("loading Sparse matrix\n");
    LoadSparseConMat();
  }

  fpMatFlag = fopen("IF_MatReady.txt", "w");
  fprintf(fpMatFlag, "1\n");
  fclose(fpMatFlag);
  
  clock_t timeStart = clock();
  RunSim();
  clock_t timeStop = clock();
  double elapsedTime = (double)(timeStop - timeStart) / CLOCKS_PER_SEC;
  cout << "\n elapsed time = " << elapsedTime << "s, or " << elapsedTime / 60.0 << "min" << endl; 
  FILE *fpet = fopen("elapsedTime.txt", "a");
  fprintf(fpet, "%llu %f\n", nSteps, elapsedTime);
  fclose(fpet);
  delete [] conMat;
  delete [] nPostNeurons;
  delete [] idxVec;
  delete [] sparseConVec;
  return 0; //EXIT_SUCCESS;
}

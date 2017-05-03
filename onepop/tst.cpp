#include <stdio.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <random>
#include <time.h>
#include <string>
#include <cstring>

using namespace::std ;

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

double UniformRand() {
  std::random_device rd;
  std::default_random_engine gen(rd());
  std::uniform_real_distribution<double> UniformRand_Aux(0.0, 1.0);
  return UniformRand_Aux(gen);
  // return (double)rand() / (double)RAND_MAX ;
}



int main(int argc, char *argv[]) {
  // std::random_device rd;
  // std::default_random_engine gen(rd());
  // std::uniform_real_distribution<double> UniformRand(0.0, 1.0);
  for(unsigned int i = 0; i < 10; i++) {
    cout << UniformRand() << endl;
  }
  exit(1);
    
  int n = 10000;
  if(argc > 3) {
    n = atoi(argv[3]);
  }
  double dPhi = M_PI / (double)n;
  double m1 = 0, phase = 0, m1In = 0, phaseIn = 0;
  m1In = atof(argv[1]);
  phaseIn = atof(argv[2]) * M_PI / 180.0;
  vector <double> x(n);  
  for(int i = 0; i < n; i++) {
    x[i] = 0.15 + m1In * cos(2.0 * (i * dPhi - phaseIn));
  }
  M1Component(x, n, &m1, &phase);
  cout << "m1 = " << m1 << " phase in deg = " << phase * 180.0 / M_PI << endl;
  x.clear();
  return 0;
}
  







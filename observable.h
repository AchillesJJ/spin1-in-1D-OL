#ifndef __OBSERVABLE_H
#define __OBSERVABLE_H

#include "itensor/all.h"
#include <typeinfo>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <cmath>
#include <cstring>
#include <sstream>
#include <complex>

using namespace itensor; 
using namespace std;

double particle_num(IQMPS psi, const S1BH_NC sites, string op_name)
{
  int N = sites.N();
  double num = 0.0;
  for (int i = 1; i <= N; i++){
    auto op = sites.op(op_name, i);
    psi.position(i);
    auto C = psi.A(i);
    C *= op;
    C *= dag(prime(psi.A(i), Site));
    num += C.real();
  }
  return num;
}

double correlation(IQMPS psi, const S1BH_NC sites, string op_name_L, string op_name_R, int r)
{
  int N = sites.N();
  int pos_L = int(N/4);
  int pos_R = 3*int(N/4);
  
  if (r>(pos_R-pos_L)){
    cout << "correlation distance must be no more than " << pos_R-pos_L << endl;
    return 0;
  }
  else{
    // calculate average correlation
    double corl = 0.0;
    for (int i = pos_L; i <= pos_R-r; i++){
      psi.position(i);
      auto C = psi.A(i);
      C *= sites.op(op_name_L, i);
      auto iR = commonIndex(psi.A(i), psi.A(i+1), Link);
      C *= dag(prime(prime(psi.A(i), Site), iR));
      for (int j = i+1; j < i+r; j++){
        C *= psi.A(j);
        C *= dag(prime(psi.A(j), Link));
      }
      C *= psi.A(i+r);
      C *= sites.op(op_name_R, i+r);
      auto iL = commonIndex(psi.A(i+r-1), psi.A(i+r), Link);
      C *= dag(prime(prime(psi.A(i+r), Site), iL));
      corl += C.real();
    }
    return corl/(pos_R-pos_L-r+1);
  }
}

#endif






































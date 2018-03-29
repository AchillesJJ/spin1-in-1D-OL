#ifndef __INITIALIZE_H
#define __INITIALIZE_H

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

/*
This routine initializes the MPS by using the results
given by Gutzwiller variational methods
*/

IQMPS initialize(const S1BH_NC sites)
{
  auto N = sites.N();
  
  // preallocate IQIndex set
  std::vector<IQIndex> b;
  for (int i = 1; i < N; i++){
    auto idx = IQIndex(nameint("b", i), 
    Index("m3", 1, Link), QN({-6,1}),
    Index("m2", 2, Link), QN({-4,1}),
    Index("m1", 4, Link), QN({-2,1}),
    Index("z0", 8, Link), QN({-0,1}),
    Index("p1", 4, Link), QN({2,1}),
    Index("p2", 2, Link), QN({4,1}),
    Index("p3", 1, Link), QN({6,1}));
    b.push_back(idx);
  }
  
  // preallocate random IQTensor
  std::vector<IQTensor> wf;
  wf.push_back(randomTensor(QN({0,1}), sites(1), b[0]));
  for (int i = 2; i < N; i++){
    wf.push_back(randomTensor(QN({0,1}), sites(i), b[i-2].dag(), b[i-1]));
  }
  wf.push_back(randomTensor(QN({0,1}), sites(N), b[N-2].dag()));
  
  // initialize MPS 
  IQMPS psi(sites);
  for (int i = 1; i <= N; i++){
    psi.setA(i, wf[i-1]);
  }
  psi.orthogonalize();
  psi = (1.0/normalize(psi))*psi;
  
  return psi;
      
}

#endif


















































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

IQMPS initialize(const S1BH_NC sites)
{
  auto N = 4;
  auto b1 = IQIndex("b1", 
  Index("m3", 1, Link), QN({-6,1}),
  Index("m2", 2, Link), QN({-4,1}),
  Index("m1", 4, Link), QN({-2,1}),
  Index("z0", 8, Link), QN({-0,1}),
  Index("p1", 4, Link), QN({2,1}),
  Index("p2", 2, Link), QN({4,1}),
  Index("p3", 1, Link), QN({6,1}));
  auto b2 = IQIndex("b2", 
  Index("m3", 1, Link), QN({-6,1}),
  Index("m2", 2, Link), QN({-4,1}),
  Index("m1", 4, Link), QN({-2,1}),
  Index("z0", 8, Link), QN({-0,1}),
  Index("p1", 4, Link), QN({2,1}),
  Index("p2", 2, Link), QN({4,1}),
  Index("p3", 1, Link), QN({6,1}));
  auto b3 = IQIndex("b3", 
  Index("m3", 1, Link), QN({-6,1}),
  Index("m2", 2, Link), QN({-4,1}),
  Index("m1", 4, Link), QN({-2,1}),
  Index("z0", 8, Link), QN({-0,1}),
  Index("p1", 4, Link), QN({2,1}),
  Index("p2", 2, Link), QN({4,1}),
  Index("p3", 1, Link), QN({6,1}));
  
  auto wf1 = randomTensor(QN({0,1}), sites(1), b1);
  auto wf2 = randomTensor(QN({0,1}), sites(2), b1.dag(), b2);
  auto wf3 = randomTensor(QN({0,1}), sites(3), b2.dag(), b3);
  auto wf4 = randomTensor(QN({0,1}), sites(4), b3.dag());
  
  IQMPS psi(sites);
  psi.setA(1, wf1);
  psi.setA(2, wf2);
  psi.setA(3, wf3);
  psi.setA(4, wf4);
  psi.orthogonalize();
  psi = (1.0/normalize(psi))*psi;
  
  return psi;
  
}

#endif




































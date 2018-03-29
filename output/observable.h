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

template <typename T1, typename T2>
double particle_num_on_sites(T1 &psi, const T2 &sites, string op_name, int id)
{
  int N = sites.N();
  psi.position(id);
  auto C = psi.A(id);
  C *= sites.op(op_name, id);
  C *= dag(prime(psi.A(id), Site));
  return C.real();
}

template <typename T1, typename T2>
double particle_num(T1 &psi, const T2 &sites, string op_name)
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

template <typename T1, typename T2>
double correlation(T1 &psi, const T2 &sites, string op_name_L, string op_name_R, int r)
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

template <typename T1, typename T2>
double bond_entropy(T1 &psi, const T2 &sites, int bd)
{
  psi.position(bd);
  IQTensor wf = psi.A(bd)*psi.A(bd+1);
  auto U = psi.A(bd);
  IQTensor S, V;
  auto spectrum = svd(wf, U, S, V);
  
  Real SvN = 0.;
  for (auto p : spectrum.eigs())
  {
    if (p>1E-12) SvN += -p*log(p);
  }
  return (double)SvN;
}

template <typename T1, typename T2>
double particle_num_fluc(T1 &psi, const T2 &sites, Args const& args = Args::global())
{
  auto total = args.getBool("total particle number", true);
  int N = sites.N();
  auto ampo1 = AutoMPO(sites);
  auto ampo2 = AutoMPO(sites);
  IQMPO Delta_1, Delta_2;
  if (total){
    for (int i = 1; i <= N; i++){
      for (int j = 1; j <= N; j++){
        ampo1 += 1.0, "Ntot", i, "Ntot", j;
      }
      ampo2 += 1.0, "Ntot", i;
    }
    Delta_1 = IQMPO(ampo1);
    Delta_2 = IQMPO(ampo2);
  }
  else{
    ampo1 += 1.0, "Ntot", int(N/2), "Ntot", int(N/2);
    ampo2 += 1.0, "Ntot", int(N/2);
    Delta_1 = IQMPO(ampo1);
    Delta_2 = IQMPO(ampo2);
  }
  
  double d1 = overlap(psi, Delta_1, psi);
  double d2 = overlap(psi, Delta_2, psi);
  
  return d1-d2*d2;
}

template <typename T1, typename T2>
void central_tensor(const T1 &psi, const T2 &sites, ofstream &ff)
{
  // define variables
  auto N = sites.N();
  auto c = int(N/2);
  auto ic = sites(c);
  auto iL = sites(c-1);
  auto iR = sites(c+1);
  auto pc = psi.A(c);
  auto pL = psi.A(c-1);
  auto pR = psi.A(c+1);
  auto bL = commonIndex(pL, pc, Link);
  auto bR = commonIndex(pc, pR, Link);
  auto dL = commonIndex(pL, psi.A(c-2), Link);
  auto dR = commonIndex(pR, psi.A(c+2), Link);
  
  /* step I : RHS factor SVD */
  IQTensor UR(ic, bL), VR(iR, dR);
  IQTensor wfR = pc*pR;
  factor(wfR, UR, VR, {"Cutoff", 1E-2});
  auto lbd_R = commonIndex(UR, VR, Link);
  
  /* step II : LHS factor SVD */
  IQTensor UL(ic, lbd_R), VL(iL, dL);
  IQTensor wfL = UR*pL;
  factor(wfL, UL, VL, {"Cutoff", 1E-2});
  auto lbd_L = commonIndex(UL, VL, Link);
  
  /* step III : analyze central tensor */
  auto A = UL;
  for (int i1 = 1; i1 <= ic.m(); i1++){
    for (int i2 = 1; i2 <= lbd_L.m(); i2++){
      for (int i3 = 1; i3 <= lbd_R.m(); i3++){
        double res = A.real(ic(i1), lbd_L(i2), lbd_R(i3));
        if (abs(res) > 1E-1){
          ff << i1 << " " << i2 << " " << i3 << " " 
             << res << endl;
        }
      }
    }
  }
  
}

#endif






































// this file consider spin-1 Bose-Hubbard model
#include "itensor/all.h"
#include "observable.h"
#include "initialize.h"
#include <typeinfo>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <cmath>
#include <cstring>
#include <sstream>
#include <random>
#include <algorithm>

using namespace itensor;
using namespace std;

int main(int argc, char* argv[])
{
  // parameters for Rb-87 F=1 atom
  int N = 40;
  double U0 = 1.0;
  double U1 = -0.005;
  double q = 0.0;
  double t = 0.0;
  double mu = 0.0;
  if (argc>0){
    q = (double) atof(argv[1]);
  }
  else{
    cout << "Qudratic Zeeman q is not set. " << endl;
  }
  
  // readme file 
  ofstream readme;
  readme.open("output/readme.txt");
  readme.close();
  readme.open("output/readme.txt", ofstream::app);
  ofstream filling_num;
  filling_num.open("output/filling_num.dat");
  filling_num.close();
  filling_num.open("output/filling_num.dat", ofstream::app);
  ofstream corl;
  corl.open("output/correlation.dat");
  corl.close();
  corl.open("output/correlation.dat", ofstream::app);
  
  // loop over different hopping t
  for (int cnt1 = 0; cnt1 <= 5; cnt1++){
    for (int cnt2 = 1; cnt2 <= 20; cnt2++){
      
      // change parameters
      mu = 0.1*cnt1;
      t = 0.015*cnt2;
      
      auto sites = S1BH_NC(N);
      auto ampo = AutoMPO(sites);
      // hopping term
      for (int i = 1; i < N; i++){
        ampo += -t, "Bpdag", i, "Bp", i+1;
        ampo += -t, "Bpdag", i+1, "Bp", i;
        ampo += -t, "Bzdag", i, "Bz", i+1;
        ampo += -t, "Bzdag", i+1, "Bz", i;
        ampo += -t, "Bmdag", i, "Bm", i+1;
        ampo += -t, "Bmdag", i+1, "Bm", i;
      }
      // on-site term
      for (int i = 1; i <= N; i++){
        // density-density interaction
        ampo += U0/2, "Ntot", i, "Ntot", i;
        ampo += -U0/2, "Ntot", i;
        // spin-spin interaction
        ampo += U1/2.0, "F+", i, "F-", i;
        ampo += U1/2.0, "Fz", i, "Fz", i;
        ampo += -U1/2.0, "Fz", i;
        ampo += -U1, "Ntot", i;
        // quafratic Zeeman term
        ampo += q, "Np", i;
        ampo += q, "Nm", i;
        // chemical potential
        ampo += -mu, "Ntot", i;
      }
      auto Hamil = IQMPO(ampo);
      
      // randomly initialize state with total Sz = 0
      auto psi = IQMPS(sites);
      for (int i = 1; i <= N; i++){
        auto s = sites(i);
        auto wf = IQTensor(s);
        wf.set(s(8), 1.0/sqrt(3.));
        wf.set(s(10), 1.0/sqrt(3.));
        wf.set(s(12), 1.0/sqrt(3.));
        psi.setA(i, wf);
      }
      // auto psi = initialize(sites);
      
      cout << particle_num(psi, sites, "Nz") << endl;
      
      // Stage-1 : DMRG pre-sweep 
      auto sweeps = Sweeps(10);
      sweeps.maxm() = 10,20,30,40,50,60,70,80,90,100; 
      sweeps.cutoff() = 1E-5,1E-5,1E-5,1E-5,1E-6,1E-6,1E-6; 
      sweeps.niter() = 3;
      sweeps.noise() = 1E-8;
      // perform DMRG algorithm
      auto energy = dmrg(psi,Hamil,sweeps,{"Quiet",true}); // ground state
      
      double total_num = 0.0;
      total_num += particle_num(psi, sites, "Np");
      total_num += particle_num(psi, sites, "Nz");
      total_num += particle_num(psi, sites, "Nm");
      filling_num << mu << " " << t << " " << total_num/N << " " 
      << particle_num(psi, sites, "Np")/N << " "
      << particle_num(psi, sites, "Nz")/N << " "
      << particle_num(psi, sites, "Nm")/N << " " << endl;
      
      corl << mu << " " << t << " "
      << correlation(psi, sites, "Bpdag", "Bp", 10) << " "
      << correlation(psi, sites, "Bzdag", "Bz", 10) << " "
      << correlation(psi, sites, "Bmdag", "Bm", 10) << endl;
      
      readme << "cnt1 = " << cnt1 << ", " 
      << "cnt2 = " << cnt2 << " is accomplished " << endl;  
      
    } // loop over hopping t
    
  }// loop over chemical potential mu
  
  // close readme 
  readme.close();
  filling_num.close();
  
} //end main













































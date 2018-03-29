// this file consider spin-1 Bose-Hubbard model
#include "itensor/all.h"
#include "observable.h"
#include "initialize.h"
#include "utility.h"
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
  double t1 = 0.0;
  double t2 = 0.0;
  double t = 0.0;
  double mu = 0.0;
  double params[10];
  int flag = 0;
  IQMPS psi;
  // read command-line arguments
  q = (double) atof(argv[1]);
  t1 = (double) atof(argv[2]);
  t2 = (double) atof(argv[3]);
  mu = (double) atof(argv[4]);
  
  // readme file 
  ofstream readme;
  readme.open("output/readme.txt");
  readme.close();
  readme.open("output/readme.txt", ofstream::app);
  ofstream corl;
  corl.open("output/correlation.dat");
  corl.close();
  corl.open("output/correlation.dat", ofstream::app);
  
  linspace(t1, t2, 10, params);
  
  // loop on parameters: hooping t
  for (auto x : params){
  
    // hopping strength
    t = x; 
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
  
    // the first point using full DMRG sweep procedure
    if (flag == 0){
      // randomly initialize state with total Sz = 0
      psi = IQMPS(sites);
      for (int i = 1; i <= N; i++){
        auto s = sites(i);
        auto wf = IQTensor(s);
        wf.set(s(8), 1.0/sqrt(3.));
        wf.set(s(10), 1.0/sqrt(3.));
        wf.set(s(12), 1.0/sqrt(3.));
        psi.setA(i, wf);
      }
  
      // Stage-1 : DMRG pre-sweep 
      auto sweeps = Sweeps(10);
      sweeps.maxm() = 10,20,30,40,50,60,70,80,90,100; 
      sweeps.cutoff() = 1E-5,1E-5,1E-5,1E-5,1E-6,1E-6,1E-6; 
      sweeps.niter() = 3;
      sweeps.noise() = 1E-8;
      // perform DMRG algorithm
      auto energy = dmrg(psi,Hamil,sweeps,{"Quiet",true});
  
      // Stage-2 : DMRG sweep until convergence
      auto en0_old = energy;
      auto rtol = 100.0;
      auto num_sweep = 0;
      auto num_sweep_max = 100;
  
      do {
        auto sweeps_new = Sweeps(1);
        sweeps_new.cutoff() = 1E-7;
        sweeps_new.niter() = 2;
        if (num_sweep < 10){
          sweeps_new.maxm() = 110+10*num_sweep;        
          sweeps_new.noise() = 1E-10;     
        }
        else if (num_sweep >= 10 && num_sweep < 20){
          sweeps_new.maxm() = 200;        
          sweeps_new.noise() = 1E-10;     
        }
        else if (num_sweep >= 20 && num_sweep <30){
          sweeps_new.maxm() = 300;
          sweeps_new.noise() = 1E-12;
        } 
        else {
          sweeps_new.maxm() = 400;
          sweeps_new.noise() = 1E-12;
        }  
        // perform DMRG algorithm 
        auto en0_new = dmrg(psi,Hamil,sweeps_new,{"Quiet",true});
        num_sweep += 1;
        rtol = abs((en0_old-en0_new)/en0_old);
        println("relative error of energy = ", rtol);
        println("extra step ", num_sweep);
        en0_old = en0_new;
      } while (rtol>1E-8); 
      
      flag += 1;
    }
    // other point use former wave function
    else{
      // Stage-1 : DMRG pre-sweep 
      auto sweeps = Sweeps(1);
      sweeps.maxm() = 400;
      sweeps.cutoff() = 1E-7; 
      sweeps.niter() = 2;
      sweeps.noise() = 1E-12;
      // perform DMRG algorithm
      auto energy = dmrg(psi,Hamil,sweeps,{"Quiet",true});
  
      // Stage-2 : DMRG sweep until convergence
      auto en0_old = energy;
      auto rtol = 100.0;
      auto num_sweep = 0;
      auto num_sweep_max = 100;
  
      do {
        auto sweeps_new = Sweeps(1);
        sweeps_new.cutoff() = 1E-7;
        sweeps_new.niter() = 2;
        sweeps_new.maxm() = 400;
        sweeps_new.noise() = 1E-12;  
        // perform DMRG algorithm 
        auto en0_new = dmrg(psi,Hamil,sweeps_new,{"Quiet",true});
        num_sweep += 1;
        rtol = abs((en0_old-en0_new)/en0_old);
        println("relative error of energy = ", rtol);
        println("extra step ", num_sweep);
        en0_old = en0_new;
      } while (rtol>1E-8);
      
      flag += 1;
    }
    
    // save wave function
    stringstream mid_site, mid_psi;
    mid_site << "output/sites_" << flag;
    mid_psi << "output/psi_" << flag;
    writeToFile(mid_site.str(), sites);
    writeToFile(mid_psi.str(), psi);
    // output correlation
    corl << x << " " 
         << correlation(psi, sites, "Bpdag", "Bp", 10) << " "
         << correlation(psi, sites, "Bzdag", "Bz", 10) << " "
         << correlation(psi, sites, "Bmdag", "Bm", 10) << endl;
  
  }
  
  
  // close readme 
  readme.close();
  corl.close();
  
} //end main













































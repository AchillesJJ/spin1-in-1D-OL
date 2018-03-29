// this file consider spin-1 Bose-Hubbard model
#include "itensor/all.h"
#include "observable.h"
#include <typeinfo>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <cmath>
#include <cstring>
#include <sstream> 

using namespace itensor;
using namespace std;

int main()
{
  S1BH_NC sites;
  readFromFile("sites_1", sites);
  IQMPS psi(sites);
  readFromFile("psi_1", psi);
  auto N = sites.N();
  
  ofstream num_on_sites;
  num_on_sites.open("particle_num_on_sites.dat");
  
  for (int i = 1; i <= N; i++){
    num_on_sites << particle_num_on_sites(psi, sites, "Np", i) << endl;
    num_on_sites << particle_num_on_sites(psi, sites, "Nz", i) << endl;
    num_on_sites << particle_num_on_sites(psi, sites, "Nm", i) << endl;
    cout << "particle number for step " << i << " is done." << endl;
  }
  
  num_on_sites.close();
  
  // calculate bond entropy
  ofstream ben;
  ben.open("bond_entropy.dat");
  
  for (int i = 1; i < N; i++){
    ben << bond_entropy(psi, sites, i) << endl;
    cout << "bond entropy for step " << i << " is done." << endl;
  }
  
  ben.close();
  
  // calculate correlation function 
  ofstream corl;
  corl.open("correlation.dat");
  
  for (int i = 1; i <= int(N/2); i++){
    corl << correlation(psi, sites, "Bpdag", "Bp", i) << endl;
    corl << correlation(psi, sites, "Bmdag", "Bm", i) << endl;
    corl << correlation(psi, sites, "Bzdag", "Bz", i) << endl;
    cout << "correlation for step " << i << " is done." << endl;
  }
  
  corl.close();
  
  // calculate particle number fluctuation
  cout << particle_num_fluc(psi, sites, {"total particle number", false}) << endl;
  
  // analyze central tensor 
  ofstream ctensor;
  ctensor.open("central_tensor.dat");
  
  central_tensor(psi, sites, ctensor);
  
  ctensor.close();
  
  
}  













































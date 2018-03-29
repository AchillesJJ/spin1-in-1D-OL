#ifndef __UTILITY_H
#define __UTILITY_H

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

// function declaration
void linspace(const double r1, const double r2, const int slice, double *ls);


// function definition
void linspace(const double r1, const double r2, const int slice, double *ls)
{
  for (auto i = 0; i < slice; i++){    
    ls[i] = r1 + i*(r2-r1)/(slice-1);
  }  
}

#endif





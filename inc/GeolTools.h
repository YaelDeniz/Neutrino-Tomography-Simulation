#ifndef GEOLTOOLS_H
#define GEOLTOOLS_H


//C++
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

std::vector< std::vector<double> >  NuPATHS3D (std::string PREM_MODEL = "prem_15layers.txt", bool LLVP = false);

#endif



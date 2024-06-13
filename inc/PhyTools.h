#ifndef PHYTOOLS_H
#define PHYTOOLS_H

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

//ROOT

#include "TH2.h"

using namespace std;


/* 

Particle Physics tool for simulations

*/

//Neutrino Cross Sections
double XSec(double E, int ptype);

//Neutrino Flux Generation

class NuFlux
{

    public: 

    //Neutrino Energies
    double Enu_min= 1.0; //GeV
    double Enu_max= 100.0;

    //Detector Locations
    std::string  DetLoc= "South Pole";

    bool DataHeader (const std::string &aline); // Get data into tabular form but ignore Headers from file

    std::vector<std::vector<double>> reshape(const std::vector<double>& vec, int a, int b); //Reshape colum vector into matrix
    
    // By default we use South Pole Flux with Sol Min
    std::vector< std::vector<double> > SetFluxData( std::string nuFluxFile = "../NuFlux/SP_AziAveraged_solmn/spl-nu-20-01-000.d" );
    
    TH2D*GetFluxHist(int flavor, std::vector<std::vector<double>> FluxData); // Create 2D Histogram

    
};




#endif

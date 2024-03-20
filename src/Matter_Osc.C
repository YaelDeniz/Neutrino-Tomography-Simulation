/*
 Changes:
    1.-HondaFlux was chaged to Flux Avaraged for al Azimuthal directions
    2.-Events account for neutrinos and Antineutrinos

*/
#include "NuRateGenerator.h"

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

//Cern Root
#include <math.h>
#include "TF1.h"
#include "TMath.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TFile.h"
//#include "TTree.h"
#include <Math/Interpolator.h>
//#include "TGraphTrueD.h"

//OSCPROB
#ifndef __CINT__
#include "PremModel.h"
#include "PMNS_Fast.h"
bool isCINT = false;
#else
bool isCINT = true;
#endif


//Delete later

// Define the PREM tables path
//#include "../prem_default.hpp"

// Macro to load OscProb library
//#include "LoadOscProb.C"

// Some functions to make nice plots
//#include "SetNiceStyle.C"

// My Math Tools
#include "MathTools.h"
#include "PhyTools.h"

using namespace std;


// Some Constants i need
//# define myPi 3.14159265358979323846  /* pi */
//# define mN   1.67492749804E-27  // Nucleons mass ~ Neutron mass
//# define MTon  1E9  //Metric MegaTon
//# define years2sec 3.154E7 // Years in Seconds

// Make oscillogram for given final flavour and MH
TGraph*  Matter_Osc(int flvi, int flvf, double Emin , double Emax,double rho , double L)
{
   
   

// Initialize your objects
  OscProb::PMNS_Fast myPMNS;

  myPMNS.SetStdPars(); // Set PDG 3-flavor parameters

 
  myPMNS.SetDensity(rho); // Set the matter density in g/cm^3
 
  
  myPMNS.SetLength(L); // Set baseline in km
 
  
  // Say whether you want antineutrinos

  // Default is false, i.e. neutrinos
  // myPMNS.SetIsNuBar(true); 

  n = 100;

  hE = (Emax-Emin)/n;
   
  TGraph *plot = new TGraph(n);

  double Ei = Emin;

  double Pif = 0.0;
   

  for (int i = 0; i < n; ++i)
  {
   
     Ei += i*hE;
    // Calculate probability of 
    // nu_i to nu_f transition
     Pif = myPMNS.Prob(flvi,flvf,Ei);

    plot->SetPoint(i,Ei,Pif) 

  }


  return plot;

}

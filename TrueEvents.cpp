#define _USE_MATH_DEFINES
#include <cmath>

#include <iostream>
#include <vector>
#include <math.h>

//Cern ROOT
#include "TH2.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "Math/Interpolator.h"
//My Libraries
#include "MathTools.h"
#include "PhyTools.h"
#include "GetTrueEvents.h"


# define myPi 3.14159265358979323846  /* pi */
// Some Constants i need
# define mN   1.67492749804E-27  // Nucleons mass ~ Neutron mass
# define MTon  1E9  //Metric MegaTon
# define years2sec 3.154E7 // Years in Seconds

using namespace std;


int main()
{
    // Plot and Histogram2D settings

    std::cout << " Neutrino Oscillation tomography. " << std::endl;
    std::cout << " Generating True/Expected  Mu-like Events. "<< std::endl; 

    // Number of bins
    int Ebins =101; // # of Bins of Energy
    int Tbins =101; // # of Bins of cosThetaZ

    //Range in Theta and E for Event Oscillogram

    //Energy interval (in GeV):
    double Emin = 1.0; 
    double Emax = 25.0;

    //Zenith Angle Interval:
    double Etamin = 33.0;
    double Etamax = 38.0;

    double Phim = 0.0;
    double PhiM = 80.0 ; 

     //DETECTOR PROPERTIES
    double DetMass = 10.0*MTon; //Mass in megaton units
    double Nn      = DetMass/mN; //Number of target nucleons in the detector (Detector Mass / Nucleons Mass)
    double T       = 20.0*years2sec; //Detector Exposure time in sec: One Year

    double NnT = Nn*T; 


    /*  
    double Emin = 4.0; 
    double Emax = 6.0;

    double Etamin = 10.0;
    double Etamax = 30.0;
    
    double Phim    = 0.0;
    double PhiM    = 360.0 ;
    */

    double dAz = PhiM-Phim;


    double E_GeV[] = {Emin,Emax};
    double Eta[] = {Etamin,Etamax};


    // neutrino state options nue (0), numu (1) or nutau (2)
    int flvf = 1;

    std::string prem_llsvp;
    
    prem_llsvp   = "/home/dehy0499/OscProb/PremTables/Prem_LLSVPhigh.txt"; //Specify PREM table from OscProb


    std::string prem_default;
            
    prem_default   = "/home/dehy0499/OscProb/PremTables/prem_default.txt"; //Specify PREM table from OscProb


    TH2D* EventOsc = GetTrueEvents(prem_llsvp, flvf, E_GeV, Eta , dAz , Ebins, Tbins ,NnT);

     std::cout << "Observed data is stored inside './SimulationResults" << std::endl;

    TCanvas *c = new TCanvas();
    
    EventOsc->Draw("SURF1 Z");
    
    gStyle->SetPalette(55);

    gPad->SetTheta(30.0); // default is 30
    gPad->SetPhi(330.0); // default is 30 --Azimuthal
    gPad->Update();
    
   
    
    gPad->SetRightMargin(0.18);
  
    return 0;
    



}
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
#include "GetObservedEvents.h"


# define myPi 3.14159265358979323846  /* pi */

using namespace std;


int main()
{
    // Plot and Histogram2D settings

    std::cout << " Neutrino Oscillation tomography. " << std::endl;
    std::cout << " Generating Observed  Mu-like Events. "<< std::endl; 

    // Number of bins
    int Ebins =101; // # of Bins of Energy
    int Tbins =101; // # of Bins of cosThetaZ


    //Energy interval (in GeV):
    double Emin = 3.0; 
    double Emax = 5.0;

    //Zenith Angle Interval:
    double Etamin = 33.0;
    double Etamax = 38.0;

    double Phim = 0.0;
    double PhiM = 80.0 ; 


    /*
    double Emin = 4.0; 
    double Emax = 6.0;

    double Etamin = 10.0;
    double Etamax = 30.0;
    
    double Phim    = 0.0;
    double PhiM    = 360.0 ;
    */

    

    double dAz = PhiM-Phim;


    double E[] = {Emin,Emax};
    double Eta[] = {Etamin,Etamax};

    int flvf = 1;

    TH2D* EventOsc = ObservedEvents(flvf, E, Eta , dAz , Ebins, Tbins);

    std::cout << "Observed data is stored inside './SimulationResults' " << std::endl;

    TCanvas *c = new TCanvas();
    
    EventOsc->Draw("SURF1 Z");
    
    gStyle->SetPalette(55);

    gPad->SetTheta(30.0); // default is 30
    gPad->SetPhi(330.0); // default is 30 --Azimuthal
    gPad->Update();
    gPad->SetRightMargin(0.18);


    return 0;

}
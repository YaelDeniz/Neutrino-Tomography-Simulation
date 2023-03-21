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
#include "NuRateGenerator.h"


# define myPi 3.14159265358979323846  /* pi */

using namespace std;


int main()
{
    // Plot and Histogram2D settings

    // Number of bins
    int nbinsx =100; // # of Bins of Energy
    int nbinsy =100; // # of Bins of cosThetaZ

    //Range in Theta and E for Event Oscillogram

    //Zenith angle interval:
    //double ctmin = cos(147.0*myPi/180); 
    //double ctmax = cos(142*myPi/180);

    double ctmin = cos(170.0*myPi/180); 
    double ctmax = cos(150*myPi/180);
    
    double dcT = (ctmax - ctmin)/(2.0*nbinsy); //< Bin width/2

    //Energy interval (in GeV):
    //double Emin = 1.0; 
    //double Emax = 20.0;


    //Energy interval (in GeV):
    double Emin = 4.0; 
    double Emax = 6.0;

    double Es[] = {Emin,Emax};
    double cts[] = {ctmin,ctmax};

    int flvf = 1;

    TH2D* EventOsc = GetTrueEvents(flvf, Es, cts , nbinsx, nbinsy);

    TCanvas *c = new TCanvas();
    
    //EventOsc->Draw("colz");    
    EventOsc->Draw("SURF2");
    gStyle->SetPalette(55);

    gPad->SetTheta(45.0); // default is 30
    gPad->SetPhi(120.0); // default is 30 --Azimuthal
    gPad->Update();
    
    //gPad->SetLogx();
    
    gPad->SetRightMargin(0.18);

    c->Print("True_Events.png");
  
    return 0;

}
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

//OSCPROB
#ifndef __CINT__
#include "PremModel.h"
#include "PMNS_Fast.h"
bool isCINT = false;
#else
bool isCINT = true;
#endif

# define myPi 3.14159265358979323846  /* pi */

using namespace std;


int main()
{

    double N_ij = 0; //<  Total # Events for standard Earth

    // Plot and Histogram2D settings

    // Number of bins
    int nbinsx = 100; // # of Bins of Energy
    int nbinsy = 10; // # of Bins of cosThetaZ

    //Range in Theta and E for Event Oscillogram

    //Zenith:
    double ctmin = cos(147.0*myPi/180); 
    double ctmax = cos(142*myPi/180);
    double dcT = (ctmax - ctmin)/(2.0*nbinsy); //< Bin width/2

    //Energy (in GeV):
    double Emin = 1.0; 
    double Emax = 20.0;

    std::cout << "Energy intervals: " <<  Emin << "-" << Emax << std::endl;

    //Create Energy vectors with  Logarithmic spaced points from Emin GeV to Emax GeV
    std::vector<double> xbins = logspace(log10(Emin),log10(Emax),nbinsx + 1);
    double Energies[nbinsx + 1];
    for (int i = 0; i < nbinsx + 1; ++i){ Energies[i] = xbins[i]; }

    //Create 2D histogram for True Events for a Standard Earth
    TH2D* N_std = new TH2D("N_std"," True Events for a Standar Earth",nbinsx, Energies ,nbinsy,ctmin,ctmax); // xbins correspond to energy values and ybins to zenith angle cost
    




    //Compute Neutrino Oscillation Probabilities: 

    // Loop over cos(theta_z) and E
    
    for(int ct=1; ct<= nbinsy ; ct++)
    {   
    
        double cosT = N_std->GetYaxis()->GetBinCenter(ct); //< This will defined a constant L por different values of ct provided Dct is Small

        //std::cout << cosT << std::endl;
        // Loop of  E  
        for(int e=1; e<= nbinsx; e++)
        {   

            // Set E from bin center
            double E  = N_std->GetXaxis()->GetBinCenter(21); // E is the midpoint of Energy bin # "e"

            std::cout << "Bin("<<e <<","<< ct <<")  -> "<< "Emid: " << E << " CosT_mid:" << cosT <<std::endl;
            //std::cout << "Elow: " << Energies[e-1] << " Emax: " << Energies[e] << " Emid: "<< (Energies[e] + Energies[e-1])/2 <<std::endl;
            //std::cout << " " << std::endl;

            double cosT_ij[] = {cosT-dcT , cosT+dcT};
            double E_ij[] = {Energies[e-1], Energies[e]}; 
            N_ij = Nflvf_ij(1, cosT_ij , E_ij);
            std::cout << "     "<< "Nij: " << N_ij <<std::endl;

            //N_ij = e;
            N_std->SetBinContent(e,ct,N_ij); //Add Number of events in bin (i,j) to the histogram

            //--->Monte Carlo Random Sample
        
        }// e loop

        //---->Add to bin
 
    }// cosT loop


    N_std->SetTitle(" #mu-like events (#nu_#mu + #bar{#nu}_#mu );E (GeV);cos#theta_{z};N ");
   
    N_std->SetTitle("Expected events for #nu_{#mu} and #bar{#nu}_#mu ");
    N_std->GetXaxis()->SetTitle("E(GeV)");
    N_std->GetYaxis()->SetTitle("Cos(#theta_{z})");
    N_std->GetZaxis()->SetTitle("N");
    N_std->SetStats(0);

    TCanvas *c = new TCanvas();
    
    N_std->Draw("colz");    
    
    gStyle->SetPalette(55);
    
    gPad->SetLogx();
    
    gPad->SetRightMargin(0.18);

    c->Print("True_Events.png");
  
    return 0;

}
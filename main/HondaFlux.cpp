#define _USE_MATH_DEFINES
#include <cmath>

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

//Cern ROOT
#include "TCanvas.h"
#include "TH2.h"
#include "TGraph.h"
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
# define R_earth 6368.0 //km
# define mN   1.67492749804E-27  // Nucleons mass ~ Neutron mass
# define MTon  1E9  //Metric MegaTon
# define years2sec 3.154E7 // Years in Seconds

using namespace std;

//------------------
int main()
{
 
NuFlux SPflux;

std::vector< std::vector<double> > FluxData = SPflux.SetFluxData();


//Matrix for Histogram & Histogram Draw

//-------------------------------------------------------------------------------------------------------------------------------
TH2D* muflux =  SPflux.GetFluxHist(1,FluxData); //MuFlux
TH2D* mubflux =  SPflux.GetFluxHist(2,FluxData); //MuBarFlux
TH2D* eflux =  SPflux.GetFluxHist(3,FluxData); //EFlux
TH2D* ebflux =  SPflux.GetFluxHist(4,FluxData); //EBarFlux

//---------------------------------------------------------------------------------------------------------------------------------------
muflux->SetStats(0);
muflux-> SetTitle("#nu_{#mu} Flux");
muflux->GetXaxis()->SetTitle("log_{10}(E/GeV)");
muflux->GetYaxis()->SetTitle("Cos(#theta_{z})");
//muflux->GetZaxis()->SetTitle("log_{10}(#Phi_{#nu})");

mubflux->SetStats(0);
mubflux-> SetTitle("#bar{#nu}_{#mu} Flux");
mubflux->GetXaxis()->SetTitle("log_{10}(E/GeV)");
mubflux->GetYaxis()->SetTitle("Cos(#theta_{z})");
//mubflux->GetZaxis()->SetTitle("log_{10}(#Phi_{#nu})");

eflux->SetStats(0);
eflux-> SetTitle("#nu_{e} Flux");
eflux->GetXaxis()->SetTitle("log_{10}(E/GeV)");
eflux->GetYaxis()->SetTitle("Cos(#theta_{z})");
//eflux->GetZaxis()->SetTitle("log_{10}(#Phi_{#nu})");

ebflux->SetStats(0);
ebflux-> SetTitle("#bar{#nu}_{e} Flux");
ebflux->GetXaxis()->SetTitle("log_{10}(E/GeV)");
ebflux->GetYaxis()->SetTitle("Cos(#theta_{z})");
//ebflux->GetZaxis()->SetTitle("log_{10}(#Phi_{#nu})");

TCanvas *c = new TCanvas();
c->Divide(2,2);

c->cd(1);
eflux->Draw("COLZ");


c->cd(2);
ebflux->Draw("COLZ");

c->cd(3);
muflux->Draw("COLZ");

c->cd(4);
mubflux->Draw("COLZ");

//c->SetCanvasSize(500, 500);

c->SetWindowSize(1200, 800);

//gStyle->SetPadLeftMargin(10000); 
//gStyle->SetPadRightMargin(0.1);




//Interpolation of the flux





    double Emax = 100;

    double Emin = 1;

    double pts = 1000; 

    double hE = (100 - 1)/pts;

    double Eo = 1;

    

    TGraph  *numu = new TGraph(pts);
    TGraph  *numub = new TGraph(pts);
    TGraph  *nue = new TGraph(pts);
    TGraph  *nueb = new TGraph(pts);

 std::cout << "----------------" << std::endl;

    for (int i = 0; i < pts; ++i)
     {
          

 

         double Ei =  Eo + i*hE;

         double logEi = log10(Ei);
         
         double logphimu = muflux->Interpolate(logEi,-0.9);
         double logphimub = mubflux->Interpolate(logEi,-0.9);
         double logphie = eflux->Interpolate(logEi,-0.9);
         double logphieb = ebflux->Interpolate(logEi,-0.9);

         

         double phimu = pow(10,logphimu)*Ei*Ei*Ei;
         double phimub = pow(10,logphimub)*Ei*Ei*Ei;
         double phie = pow(10,logphie)*Ei*Ei*Ei;
         double phieb = pow(10,logphieb)*Ei*Ei*Ei;

         numu -> SetPoint (i,  logEi,  phimu);
         numub -> SetPoint (i,  logEi,  phimub);
         nue -> SetPoint (i,  logEi,  phie);
         nueb -> SetPoint (i,  logEi,  phieb);
     } 



    TCanvas *c1 = new TCanvas();
    TMultiGraph * multi = new TMultiGraph();

    numu->SetLineColor(4);
    nue->SetLineColor(2);
    numub->SetLineColor(4);
    nueb->SetLineColor(2);

    numub->SetLineStyle(9);
    nueb->SetLineStyle(9);

    multi->Add(numu);
    multi->Add(numub);
    multi->Add(nue);
    multi->Add(nueb);



    multi->GetXaxis()->SetLimits(0.0,2.0);
    multi->GetXaxis()->SetTitle("log_{10}(E/GeV)");

    multi->GetYaxis()->SetLimits(0.0,0.03);
    multi->GetYaxis()->SetTitle("#Phi_{#nu}*E^{3}");

    multi->SetTitle("Bilinear interpolation of Honda flux");

    gPad->SetGrid(1,1);
    multi->Draw("A");

    c->SaveAs("SimulationResults/NeutrinoFlux/FluxHistSP.jpg");
    c1->SaveAs("SimulationResults/NeutrinoFlux/FluxInterpolSP.jpg");

   return 0;
   
}

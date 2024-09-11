#define _USE_MATH_DEFINES
#include <cmath>

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

//Cern ROOT
#include "TApplication.h"
#include "TRootCanvas.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TH3.h"
#include "TLine.h"
#include "TGraph.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "Math/Interpolator.h"
//My Libraries
#include "MathTools.h"
#include "PhyTools.h"
#include "GetObservedEvents.h"


# define Rcmb 3480.0 // Core-Mantle-Boundary Radius (km)
# define Rocean 6368.0 // Earth radius in (km)
# define Rearth 6371.0 // Earth radius in (km)
# define Ratm 6386.0 // Radius of the Atmosphere (km)
# define mN   1.67492749804E-27  // Nucleons mass (~ Neutron mass) in kg
# define MTon  1E9  //Metric MegaTon in kg
# define years2sec 3.154E7 // Years in Seconds


using namespace std;

int main(int argc, char **argv)
{
    //Detector setting 
    double DetMass = 10.0*MTon; //Mass in megaton units
    double Nn      = DetMass/mN; //Number of target nucleons in the detector (Detector Mass / Nucleons Mass)
    double T       = 10.0*years2sec; //Detector Exposure time in sec: One Year
    double NT = Nn*T; // Exposure [Mton*years]

    double Rdet = Rearth;

    //Set detector location
    double rdet[3] = {0.0,0.0, -1*Rdet };

    //Set Neutrino Flux
    std::string FluxFolder = "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/NuFlux/";
    //std::string FluxFile = "GRN_AziAveraged_solmin/grn-nu-20-01-000.d"; // Sea
    std::string FluxFile = "SP_AziAveraged_solmin/spl-nu-20-01-000.d"; //South Pole
    std::string FluxTable = FluxFolder + FluxFile; //Class Assumes Sout Pole flux as default

    // Earth Models
    std::string PremFolder ="/home/dehy0499/OscProb/PremTables/";
    std::string PremName = "prem_44layers.txt";
    std::string path2prem = PremFolder+PremName;

    //Pile Set Up---------------------------------------------------------------
    double PileHeight = 600; // Height of LLVP in km
    double PileRadius = Rcmb + PileHeight; //km
    double DepthMin = Rcmb-1000;
    double DepthMax = PileRadius+200; // Distance from the center of the Earth
    double PileDensityPct = 2.0; // 2% density contrats for LLVP
    std::string shape = "pancake"; // LLVP shape


    //Simulation Setup

    // Binning------------------------------------------------------------------
    int cthtruebins=20 ; // # Bins in zenith/cos(zenith)
    int atruebins =100; // # Bins in azimuth (optimal bins are 110 or 22)
    int etruebins =20; // bins in energy

      // Binning------------------------------------------------------------------
    int cthrecobins=20 ; // # Bins in zenith/cos(zenith)
    int arecobins =100; // # Bins in azimuth (optimal bins are 110 or 22)
    int erecobins =20; // bins in energy

    //Reconstructed Energy interval (in GeV):
    double Emin = 2.0; 
    double Emax = 10.0;

    //Reconstructed Directions

    //Zenith Angle Interval-----------------------------------------------------
    double thmin = 180-TMath::ASin( DepthMax/Rocean )*(180.0/TMath::Pi()) ; // min 90
    double thmax = 180-TMath::ASin( (DepthMin)/Rocean )*(180.0/TMath::Pi()) ; // max 180

    double cthmin = cos(thmax*TMath::Pi()/180.0);// Cos(thmin) - Minimum possible is -1    
    double cthmax = cos(thmin*TMath::Pi()/180.0);// Cos(thmax) - Maximum possible is 0

    //Azimuthal Interval---------------------------------------------------------
    double phimin = -55.0;
    double phimax =  55.0 ;

    int nuflv = 1; // neutrino  final state: nue (0), numu (1) or nutau (2)
    
   //Standart Earth-------------------------------------------------------------
   AsimovObsSimulation StandardEarth;

   StandardEarth.ThePremTable = path2prem;
   StandardEarth.TheHondaTable = FluxTable;
   StandardEarth.SetDetectorXYZ(rdet);
   StandardEarth.PileInModel = false;

   //Simulation of Standard Earth
   StandardEarth.SetIntervals(thmin,thmax,phimin,phimax,Emin,Emax);
   StandardEarth.SetTrueBinning(cthtruebins,atruebins,etruebins);
   StandardEarth.SetRecoBinning(cthrecobins,arecobins,erecobins);
   StandardEarth.SetExposure(NT);
   StandardEarth.flvf=nuflv;

    //Alternative Earth-------------------------------------------------------------
   AsimovObsSimulation AlternativeEarth;

   AlternativeEarth.ThePremTable = path2prem;
   AlternativeEarth.TheHondaTable = FluxTable;
   AlternativeEarth.SetDetectorXYZ(rdet);
   AlternativeEarth.SetEnergyResolution(0.001,0.0);
   AlternativeEarth.SetAngularResolution(0.001,0.0);
   
   //Set LLVP
   AlternativeEarth.PileInModel = true;
   AlternativeEarth.ShapeOfPile = shape;
   AlternativeEarth.ThePileHight = 1000.0;
   AlternativeEarth.ThePileAperture = 45;
   AlternativeEarth.ThePileDensityContrast = 2.0;
   AlternativeEarth.ThePileChemicalContrast = 0.0;
   
   //Simulation of Alternative Earth
   AlternativeEarth.SetIntervals(thmin,thmax,phimin,phimax,Emin,Emax);
   AlternativeEarth.SetTrueBinning(cthtruebins,atruebins,etruebins);
   AlternativeEarth.SetRecoBinning(cthrecobins,arecobins,erecobins);
   AlternativeEarth.SetExposure(NT);
   AlternativeEarth.flvf=nuflv;


   std::vector<TH2D*>  ObsStd = StandardEarth.GetObsEvents3Dcth();
   std::vector<TH2D*>  ObsAlt = AlternativeEarth.GetObsEvents3Dcth();

   TApplication app("app", &argc, argv);


   //double cthbottom = -0.8376; 
   double cthbottom = TMath::Cos(TMath::Pi()-TMath::ASin( 3480.0/Rdet ));
   double cthmid_bottom = TMath::Cos(TMath::Pi()-TMath::ASin( (3480.0+300)/Rdet ));;
   double cthmid_top = TMath::Cos(TMath::Pi()-TMath::ASin( (3480.0+600)/Rdet ));;
   double cthtop = TMath::Cos(TMath::Pi()-TMath::ASin( (3480.0 + 1000)/Rdet ));;

   TLine * lbottom = new TLine(cthbottom,Emin,cthbottom,Emax);
   TLine * lmid_bottom = new TLine(cthmid_bottom,Emin,cthmid_bottom,Emax);
   TLine * lmid_top = new TLine(cthmid_top,Emin,cthmid_top,Emax);
   TLine * ltop = new TLine(cthtop,Emin,cthtop,Emax);

   lbottom -> SetLineColor(kRed);
   lmid_bottom -> SetLineColor(kRed);
   lmid_top -> SetLineColor(kRed);
   ltop -> SetLineColor(kRed);

   lbottom -> SetLineStyle(2);
   lmid_bottom -> SetLineStyle(2);
   lmid_top -> SetLineStyle(2);
   ltop -> SetLineStyle(2);

   TCanvas *c1 = new TCanvas();

   gStyle->SetPalette(kBlueGreenYellow);

   c1->Divide(5,5);

   for(int i = 0 ; i < 25 ; ++i )

   {

    int nhist = (51 + 2*i)-1;
    c1->cd(i+1);

    TH2D * TrueDiff2D = new TH2D(Form("ObsDiff2D%d",nhist),Form("OscObs%d",nhist), cthrecobins,cthmin,cthmax,erecobins,Emin,Emax); //binning in cth 
    GetDiff2D( ObsStd[nhist] , ObsAlt[nhist], TrueDiff2D );
    TrueDiff2D->Draw("COLZ");
    TrueDiff2D->SetStats(0);
    lbottom->Draw("same");
    lmid_bottom->Draw("same");
    lmid_top->Draw("same");
    ltop->Draw("same");


    //std::string BinLabel = std::to_string(cthrecobins)+std::to_string(arecobins)+std::to_string(erecobins);
    //std::string filename  = "ObsDiff_"+std::to_string(nuflv)+shape+BinLabel+std::to_string(nhist)+".txt";
    //ExportToCSV(TrueDiff2D, filename);

 
   }

    gPad->Update();
    c1->Modified(); c1->Update();
    TRootCanvas *rc = (TRootCanvas *)c1->GetCanvasImp();
    rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
    app.Run();

    return 0;
   


}


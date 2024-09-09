
/* Direction of neutrinos a re given in degrees*/
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
#include "PlotTool.h"
#include "StatsAna.h"
#include "GetTrueEvents.h"


# define myPi 3.14159265358979323846  /* pi */
// Some Constants i need

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
    // Plot and Histogram2D settings

    std::cout << " Neutrino Oscillation tomography. " << std::endl;

    //Detetcor Settings
    
    //Detector size-------------------------------------------------------------
    double DetMass = 10.0*MTon; //Mass in megaton units
    double Nn      = DetMass/mN; //Number of target nucleons in the detector (Detector Mass / Nucleons Mass)
    double T       = 10.0*years2sec; //Detector Exposure time in sec: One Year
    double NnT = Nn*T; // Exposure [Mton*years]

    //Detector location 
    double Rdet = Rearth; //South Pole

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
    double PileHeight = 1000; // Height of LLVP in km
    double PileRadius = Rcmb + PileHeight; //km
    double DepthMin = Rcmb-2500;
    double DepthMax = PileRadius+500; // Distance from the center of the Earth
    double PileDensityPct = 2.0; // 2% density contrats for LLVP

    //SIMULATION SET UP---------------------------------------------------------

    // Binning------------------------------------------------------------------
    int czbins=50 ; // # Bins in zenith/cos(zenith)
    int abins =100; // # Bins in azimuth (optimal bins are 110 or 22)
    int ebins =50; // bins in energy

    //Energy interval (in GeV)--------------------------------------------------
    double EnuMin=1.0 ; 
    double EnuMax=10.0 ;

    //Zenith Angle Interval-----------------------------------------------------
    double zenmin = 180-TMath::ASin( DepthMax/Rdet )*(180.0/TMath::Pi()) ; // min 90
    double zenmax = 180-TMath::ASin( (DepthMin)/Rdet )*(180.0/TMath::Pi()) ; // max 180

    double Czmin = cos(zenmax*TMath::Pi()/180.0);// Cos(zenmin) - Minimum possible is -1    
    double Czmax = cos(zenmin*TMath::Pi()/180.0);// Cos(zenmax) - Maximum possible is 0

    //Azimuthal Interal---------------------------------------------------------
    double phi_a = -50;
    double phi_b =  50;
    double dPhi = (phi_b-phi_a)/(2.0*abins);
    
    // for bins centered
    double phimin = phi_a - dPhi; 
    double phimax = phi_b + dPhi;

    

    

   //SIMULATION-----------------------------------------------------------------
   int nuflv = 1; // neutrino  final state: nue (0), numu (1) or nutau (2)
    
   //Standart Earth-------------------------------------------------------------
   AsimovSimulation StandardEarth;

   StandardEarth.PremTable = path2prem;
   StandardEarth.HondaTable = FluxTable;
   StandardEarth.SetDetectorPosition(rdet);
   StandardEarth.MantleAnomaly = false;
   StandardEarth.SetIntervals(zenmin,zenmax,phimin,phimax,EnuMin,EnuMax);
   StandardEarth.SetBinning(czbins,abins,ebins);
   StandardEarth.SetExposure(NnT);
   StandardEarth.flvf=nuflv;

   //TH2D * TrueStd = StandardEarth.GetTrueEvents2D();

    std::vector< TH2D* >  TrueStd = StandardEarth.GetTrueEvents3D();

    //Alternative Earth--------------------------------------------------------- 

   AsimovSimulation AlternativeEarth;

   AlternativeEarth.PremTable = path2prem;
   AlternativeEarth.HondaTable = FluxTable;
   AlternativeEarth.SetDetectorPosition(rdet);
   AlternativeEarth.PileHeight = 1000;
   AlternativeEarth.aperture=45;
   AlternativeEarth.MantleAnomaly = true;
   AlternativeEarth.AnomalyShape="pancake";
   AlternativeEarth.PileDensityContrast = 2;
   AlternativeEarth.PileChemContrast = 0.0;
   AlternativeEarth.SetIntervals(zenmin,zenmax,phimin,phimax,EnuMin,EnuMax);
   AlternativeEarth.SetBinning(czbins,abins,ebins);
   AlternativeEarth.SetExposure(NnT);
   AlternativeEarth.flvf=nuflv;

   //TH2D * TrueAlt = AlternativeEarth.GetTrueEvents2D();

    std::vector< TH2D* > TrueAlt= AlternativeEarth.GetTrueEvents3D();


  //Data Analysis---------------------------------------------------------------
   
    //TH3D * TrueDiff3D = new TH3D("TrueHist","True Event Histrogram", czbins,Czmin,Czmax,abins,phimin,phimax,ebins,EnuMin,EnuMax); //binning in cth 

/*

    TGraph * chi2plot = new TGraph(7);

    for (int i = 0; i <= 6; ++i)
    {

        AlternativeEarth.PileDensityContrast = -3 + i;
        std::vector< TH2D* > TrueAlt= AlternativeEarth.GetTrueEvents3D();

        double chi2tot = 0;

        for (int n = 0; n < TrueAlt.size(); ++n)
        {
            chi2tot += Get2DChi2( TrueStd[n] , TrueAlt[n]);
        }

        chi2plot->AddPoint(-3+i,chi2total);

        
    }
    
*/

    std::cout << " vector size "  << TrueStd.size() << std::endl;  

    TH2D * TrueDiff2D = new TH2D("TrueHist","True Event Histrogram", czbins,Czmin,Czmax,ebins,EnuMin,EnuMax); //binning in cth 


   TApplication app("app", &argc, argv);


   //double cthbottom = -0.8376; 
   double cthbottom = TMath::Cos(TMath::Pi()-TMath::ASin( 3480.0/Rdet ));
   double cthmid_bottom = TMath::Cos(TMath::Pi()-TMath::ASin( (3480.0+300)/Rdet ));;
   double cthmid_top = TMath::Cos(TMath::Pi()-TMath::ASin( (3480.0+600)/Rdet ));;
   double cthtop = TMath::Cos(TMath::Pi()-TMath::ASin( (3480.0 + 1000)/Rdet ));;

   std::cout << cthtop << std::endl;

   TLine * lbottom = new TLine(cthbottom,EnuMin,cthbottom,EnuMax);
   TLine * lmid_bottom = new TLine(cthmid_bottom,EnuMin,cthmid_bottom,EnuMax);
   TLine * lmid_top = new TLine(cthmid_top,EnuMin,cthmid_top,EnuMax);
   TLine * ltop = new TLine(cthtop,EnuMin,cthtop,EnuMax);

   TCanvas *c = new TCanvas();

   gStyle->SetPalette(kBlueGreenYellow);

    c->Divide(2,3,0.01,0.01);

    c->cd(1);

    GetDiff2D( TrueStd[51] , TrueAlt[51], TrueDiff2D );
    TrueDiff2D->Draw("COLZ");
    TrueDiff2D->SetStats(0);
    lbottom->Draw("same");
    lmid_bottom->Draw("same");
    lmid_top->Draw("same");
    ltop->Draw("same");

    c->cd(2);

    GetDiff2D( TrueStd[61] , TrueAlt[61], TrueDiff2D );
    TrueDiff2D->Draw("COLZ");
    TrueDiff2D->SetStats(0);
    lbottom->Draw("same");
    lmid_bottom->Draw("same");
    lmid_top->Draw("same");
    ltop->Draw("same");
    
    c->cd(3);
    GetDiff2D( TrueStd[71] , TrueAlt[71], TrueDiff2D );
    TrueDiff2D->Draw("COLZ");
    TrueDiff2D->SetStats(0);
    lbottom->Draw("same");
    lmid_bottom->Draw("same");
    lmid_top->Draw("same");
    ltop->Draw("same");

    c->cd(4);
    GetDiff2D( TrueStd[81] , TrueAlt[81], TrueDiff2D );
    TrueDiff2D->Draw("COLZ");
    TrueDiff2D->SetStats(0);
    lbottom->Draw("same");
    lmid_bottom->Draw("same");
    lmid_top->Draw("same");
    ltop->Draw("same");


    c->cd(5);
    GetDiff2D( TrueStd[91] , TrueAlt[91], TrueDiff2D );
    TrueDiff2D->Draw("COLZ");
    TrueDiff2D->SetStats(0);
    lbottom->Draw("same");
    lmid_bottom->Draw("same");
    lmid_top->Draw("same");
    ltop->Draw("same");



    gPad->Update();
    c->Modified(); c->Update();
    TRootCanvas *rc = (TRootCanvas *)c->GetCanvasImp();
    rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
    app.Run();

    return 0;

}


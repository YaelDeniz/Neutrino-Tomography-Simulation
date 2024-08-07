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
void GetDiff2D( TH2D * histstd , TH2D * histalt, TH2D * diff);


int main(int argc, char **argv)
{
    // Plot and Histogram2D settings

    std::cout << " Neutrino Oscillation tomography. " << std::endl;
    std::cout << " Observed Events: track-like only. "<< std::endl; 

    int ErecoBins=40; // # of Bins of Energy
    int CthRecoBins=40; // # of Bins of cosEta

    int EtrueBins=100; // # of Bins of Energy
    int CthTrueBins=10; // # of Bins of cosEta
    int AziTrueBins=100;


    //Target SetUp

    //Pile Set Up---------------------------------------------------------------
    double PileHeight = 600; // Height of LLVP in km
    double PileRadius = Rcmb + PileHeight; //km
    double DepthMin = Rcmb-1000;
    double DepthMax = PileRadius+200; // Distance from the center of the Earth
    double PileDensityPct = 2.0; // 2% density contrats for LLVP

    //Detector SetUp

    //Set detector location
    double xyz[3] = {0.0,0.0, -1*Rocean };

    //Set Neutrino Flux
    std::string FluxFolder = "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/NuFlux/";
    std::string FluxFile = "GRN_AziAveraged_solmin/grn-nu-20-01-000.d";
    std::string FluxTable = FluxFolder + FluxFile; //Class Assumes Sout Pole flux as default


    //Detector size-------------------------------------------------------------
    double DetMass = 10.0*MTon; //Mass in megaton units
    double Nn      = DetMass/mN; //Number of target nucleons in the detector (Detector Mass / Nucleons Mass)
    double T       = 10.0*years2sec; //Detector Exposure time in sec: One Year
    double NnT = Nn*T; // Exposure [Mton*years]

    //Detector Parameters

    //....

    //Reconstructed Energy interval (in GeV):
    double Emin = 2.0; 
    double Emax = 10.0;

    //Reconstructed Directions

    //Zenith Angle Interval-----------------------------------------------------
    double zenmin = 180-TMath::ASin( DepthMax/Rocean )*(180.0/TMath::Pi()) ; // min 90
    double zenmax = 180-TMath::ASin( (DepthMin)/Rocean )*(180.0/TMath::Pi()) ; // max 180

    double Czmin = cos(zenmax*TMath::Pi()/180.0);// Cos(zenmin) - Minimum possible is -1    
    double Czmax = cos(zenmin*TMath::Pi()/180.0);// Cos(zenmax) - Maximum possible is 0

    //Azimuthal Interal---------------------------------------------------------
    double phimin = -55.0;
    double phimax =  55.0 ;

    // Earth Models
    std::string PremFolder ="/home/dehy0499/OscProb/PremTables/";
    std::string PremName = "prem_44layers.txt";
    std::string path2prem = PremFolder+PremName;

    //SIMULATION-----------------------------------------------------------------
   int nuflv = 1; // neutrino  final state: nue (0), numu (1) or nutau (2)
    
   //Standart Earth-------------------------------------------------------------
   AsimovObsSimulation StandardEarth;

   StandardEarth.ThePremTable = path2prem;
   StandardEarth.TheHondaTable = FluxTable;
   StandardEarth.SetDetectorXYZ(xyz);
   StandardEarth.PileInModel = false;
   StandardEarth.SetIntervals(zenmin,zenmax,phimin,phimax,Emin,Emax);
   StandardEarth.SetTrueBinning(CthTrueBins,AziTrueBins,EtrueBins);
   StandardEarth.SetRecoBinning(CthRecoBins,AziTrueBins,ErecoBins);
   StandardEarth.SetExposure(NnT);
   StandardEarth.flvf=nuflv;

      //Standart Earth-------------------------------------------------------------
   AsimovObsSimulation AlternativeEarth;

   AlternativeEarth.ThePremTable = path2prem;
   AlternativeEarth.TheHondaTable = FluxTable;
   AlternativeEarth.SetDetectorXYZ(xyz);
   AlternativeEarth.SetEnergyResolution(0.001,0.0);
   AlternativeEarth.SetAngularResolution(0.001,0.0);
   
   //Set LLVP
   AlternativeEarth.PileInModel = true;
   AlternativeEarth.ShapeOfPile = "cake";
   AlternativeEarth.ThePileHight = 1000.0;
   AlternativeEarth.ThePileAperture = 45;
   AlternativeEarth.ThePileDensityContrast = 2.0;
   AlternativeEarth.ThePileChemicalContrast = 0.0;
   

   AlternativeEarth.SetIntervals(zenmin,zenmax,phimin,phimax,Emin,Emax);
   AlternativeEarth.SetTrueBinning(CthTrueBins,AziTrueBins,EtrueBins);
   AlternativeEarth.SetRecoBinning(CthRecoBins,AziTrueBins,ErecoBins);
   AlternativeEarth.SetExposure(NnT);
   AlternativeEarth.flvf=nuflv;


   TH2D * ObsStd = StandardEarth.GetObsEvents2Dcth();
   TH2D * ObsAlt = AlternativeEarth.GetObsEvents2Dcth();
   TH2D * ObsDiff2D = new TH2D("ObsHist","Obs Event Histrogram", CthRecoBins,Czmin,Czmax,ErecoBins,Emin,Emax); //binning in cth 
   GetDiff2D( ObsStd , ObsAlt, ObsDiff2D );


   TApplication app("app", &argc, argv);

   TCanvas *c = new TCanvas();
   gStyle->SetPalette(kBlueGreenYellow);

   //double cthbottom = -0.8376; 
   double cthbottom = TMath::Cos(TMath::Pi()-TMath::ASin( 3480.0/Rocean ));
   double cthmid_bottom = TMath::Cos(TMath::Pi()-TMath::ASin( (3480.0+300)/Rocean ));;
   double cthmid_top = TMath::Cos(TMath::Pi()-TMath::ASin( (3480.0+600)/Rocean ));;
   double cthtop = TMath::Cos(TMath::Pi()-TMath::ASin( (3480.0 + 1000)/Rocean ));;


   TLine * lbottom = new TLine(cthbottom,1.0,cthbottom,25);
   TLine * lmid_bottom = new TLine(cthmid_bottom,1.0,cthmid_bottom,25);
   TLine * lmid_top = new TLine(cthmid_top,1.0,cthmid_top,25);
   TLine * ltop = new TLine(cthtop,1.0,cthtop,25);



   ObsDiff2D->Draw("COLZ");
   ObsDiff2D->SetStats(0);

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

void GetDiff2D( TH2D * histstd , TH2D * histalt, TH2D * diff)
{


   std::ofstream TrueDiffFile("SimulationResults/TrueEventsResults/EarthSensitivity2D.csv"); 

   double cth, e, Nexp, Nobs, dN;


        for (int i = 1; i <= diff->GetXaxis()->GetNbins(); ++i)
        {
            cth=diff->GetXaxis()->GetBinCenter(i);
         
            for (int k = 1; k <= diff->GetYaxis()->GetNbins(); ++k)
            {

                e=diff->GetYaxis()->GetBinCenter(k);
                Nexp = histstd->GetBinContent(i,k);
                Nobs = histalt->GetBinContent(i,k);
                dN = 100*(Nobs-Nexp)/Nexp;

                diff->SetBinContent(i,k, dN);

                TrueDiffFile << cth  << " , " << e << " , " << dN << std::endl;

            }
        }
        
    TrueDiffFile.close();


}

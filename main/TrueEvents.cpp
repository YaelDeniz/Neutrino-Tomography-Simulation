
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

void GetDiff3D( TH3D * histstd , TH3D * histalt, TH3D * diff);
void GetDiff2D( TH2D * histstd , TH2D * histalt, TH2D * diff); 
double Get2DChi2( TH2D * histstd, TH2D * histalt);
double Get3DChi2( TH3D * histstd, TH3D * histalt);

int main(int argc, char **argv)
{
    // Plot and Histogram2D settings

    std::cout << " Neutrino Oscillation tomography. " << std::endl;

    //Set detector location
    double rdet[3] = {0.0,0.0, -1*Rocean };

    //Detector size-------------------------------------------------------------
    double DetMass = 10.0*MTon; //Mass in megaton units
    double Nn      = DetMass/mN; //Number of target nucleons in the detector (Detector Mass / Nucleons Mass)
    double T       = 10.0*years2sec; //Detector Exposure time in sec: One Year
    double NnT = Nn*T; // Exposure [Mton*years]

    //Set Neutrino Flux
    std::string FluxFolder = "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/NuFlux/";
    std::string FluxFile = "GRN_AziAveraged_solmin/grn-nu-20-01-000.d";
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
    int czbins=100 ; // # Bins in zenith/cos(zenith)
    int abins =1; // # Bins in azimuth
    int ebins =100; // bins in energy

    //Energy interval (in GeV)--------------------------------------------------
    double EnuMin=1.0 ; 
    double EnuMax=25.0 ;

    //Zenith Angle Interval-----------------------------------------------------
    double zenmin = 180-TMath::ASin( DepthMax/Rocean )*(180.0/TMath::Pi()) ; // min 90
    double zenmax = 180-TMath::ASin( (DepthMin)/Rocean )*(180.0/TMath::Pi()) ; // max 180

    double Czmin = cos(zenmax*TMath::Pi()/180.0);// Cos(zenmin) - Minimum possible is -1    
    double Czmax = cos(zenmin*TMath::Pi()/180.0);// Cos(zenmax) - Maximum possible is 0

    //Azimuthal Interal---------------------------------------------------------
    double phimin = -55.0;
    double phimax =  55.0 ;



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

   TH2D * TrueStd = StandardEarth.GetTrueEvents2D();

    //Alternative Earth--------------------------------------------------------- 

   AsimovSimulation AlternativeEarth;

   AlternativeEarth.PremTable = path2prem;
   AlternativeEarth.HondaTable = FluxTable;
   AlternativeEarth.SetDetectorPosition(rdet);
   AlternativeEarth.PileHeight = 1000;
   AlternativeEarth.aperture=45;
   AlternativeEarth.MantleAnomaly = true;
   AlternativeEarth.AnomalyShape="cake";
   AlternativeEarth.PileDensityContrast = 2;
   AlternativeEarth.PileChemContrast = 0.0;
   //AlternativeEarth.ModifyLayer(23,5.0,0.0);
   //AlternativeEarth.AnomalousLayers = layers;
   AlternativeEarth.SetIntervals(zenmin,zenmax,phimin,phimax,EnuMin,EnuMax);
   AlternativeEarth.SetBinning(czbins,abins,ebins);
   AlternativeEarth.SetExposure(NnT);
   AlternativeEarth.flvf=nuflv;

   TH2D * TrueAlt = AlternativeEarth.GetTrueEvents2D();

  //Data Analysis---------------------------------------------------------------
   
    //TH3D * TrueDiff3D = new TH3D("TrueHist","True Event Histrogram", czbins,Czmin,Czmax,abins,phimin,phimax,ebins,EnuMin,EnuMax); //binning in cth 

    TH2D * TrueDiff2D = new TH2D("TrueHist","True Event Histrogram", czbins,Czmin,Czmax,ebins,EnuMin,EnuMax); //binning in cth 


   GetDiff2D( TrueStd , TrueAlt, TrueDiff2D );

   TApplication app("app", &argc, argv);

   TCanvas *c = new TCanvas();
   gStyle->SetPalette(kBlueGreenYellow);

   //double cthbottom = -0.8376; 
   double cthbottom = TMath::Cos(TMath::Pi()-TMath::ASin( 3480.0/Rocean ));
   double cthmid_bottom = TMath::Cos(TMath::Pi()-TMath::ASin( (3480.0+300)/Rocean ));;
   double cthmid_top = TMath::Cos(TMath::Pi()-TMath::ASin( (3480.0+600)/Rocean ));;
   double cthtop = TMath::Cos(TMath::Pi()-TMath::ASin( (3480.0 + 1000)/Rocean ));;

   std::cout << cthtop << std::endl;

   TLine * lbottom = new TLine(cthbottom,1.0,cthbottom,25);
   TLine * lmid_bottom = new TLine(cthmid_bottom,1.0,cthmid_bottom,25);
   TLine * lmid_top = new TLine(cthmid_top,1.0,cthmid_top,25);
   TLine * ltop = new TLine(cthtop,1.0,cthtop,25);



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



void GetDiff3D( TH3D * histstd , TH3D * histalt, TH3D * diff)
{


   double cth, azi, e, Nexp, Nobs, dN;


    std::string filename = "dNTrueAziFixLayered.csv";
    std::ofstream TrueDiffFile("SimulationResults/TrueEventsResults/3DEarth/"+filename); 
        

    for (int j = 1; j <= diff->GetYaxis()->GetNbins(); ++j)
    {
        
        azi=diff->GetYaxis()->GetBinCenter(j);

        for (int i = 1; i <= diff->GetXaxis()->GetNbins(); ++i)
        {
            cth=diff->GetXaxis()->GetBinCenter(i);
         
            for (int k = 1; k <= diff->GetZaxis()->GetNbins(); ++k)
            {

                e=diff->GetZaxis()->GetBinCenter(k);
                Nexp = histstd->GetBinContent(i,j,k);
                Nobs = histalt->GetBinContent(i,j,k);
                dN = 100*(Nobs-Nexp)/Nexp;

                diff->SetBinContent(i,j,k,dN);

                std::cout<< diff->GetXaxis()->GetNbins() << " " << azi << " , " <<  cth  << " , " << e << " , " << dN << std::endl;

                TrueDiffFile << cth << " , " <<  azi  << " , " << e << " , " << dN << std::endl;

            }
        }
        
    }

    TrueDiffFile.close();
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

double Get3DChi2( TH3D * histstd, TH3D * histalt)
{

    double  Nexp, Nobs;

    double chi2 = 0;

    for (int j = 1; j <= histstd->GetYaxis()->GetNbins(); ++j) //loop in azimuth
    {
        
        for (int i = 1; i <= histstd->GetXaxis()->GetNbins(); ++i) //loop in zenith
        {
         
            for (int k = 1; k <= histstd->GetZaxis()->GetNbins(); ++k) //loop in energy
            {

                Nexp = histstd->GetBinContent(i,j,k);
                Nobs = histalt->GetBinContent(i,j,k);
                chi2 +=  2*( Nexp - Nobs + Nobs*TMath::Log(Nobs/Nexp) ); // LLRT

            }
        }
        
    }

    return chi2;

}

double Get2DChi2( TH2D * histstd, TH2D * histalt)
{

    double  Nexp, Nobs;

    double chi2 = 0;

    for (int j = 1; j <= histstd->GetYaxis()->GetNbins(); ++j) //loop in azimuth
    {
        
        for (int i = 1; i <= histstd->GetXaxis()->GetNbins(); ++i) //loop in zenith
        {
        

                Nexp = histstd->GetBinContent(i,j);
                Nobs = histalt->GetBinContent(i,j);
                chi2 +=  2*( Nexp - Nobs + Nobs*TMath::Log(Nobs/Nexp) ); // LLRT

            
        }
        
    }

    return chi2;

}

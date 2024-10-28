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

void ExportToCSV(TH2D* hist, std::string filename);

int main(int argc, char **argv)
{
    //Detector setting 
    double M  = 10.0*MTon; //Mass in megaton units
    double T  = 10.0*years2sec; //Detector Exposure time in sec: One Year
    double MT = M*T; // Exposure [Mton*years]

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
    double PileHeight = 1000; // Height of LLVP in km
    double PileRadius = Rcmb + PileHeight; //km
    double DepthMin = Rcmb-2500;
    double DepthMax = PileRadius+500; // Distance from the center of the Earth
    double PileDensityPct = 2.0; // 2% density contrats for LLVP
    std::string shape = "cake";


    //Simulation Setup

    // Binning------------------------------------------------------------------
    int cthtruebins=100; // # Bins in zenith/cos(zenith)
    int atruebins =100; // # Bins in azimuth (optimal bins are 110 or 22)
    int etruebins =100; // bins in energy

      // Binning------------------------------------------------------------------
    int cthrecobins=40; // # Bins in zenith/cos(zenith)
    int arecobins =100; // # Bins in azimuth (optimal bins are 110 or 22)
    int erecobins =40; // bins in energy

    //Reconstructed Energy interval (in GeV):
    double Emin = 2.0; 
    double Emax = 10.0;

    //Reconstructed Directions

    //Zenith Angle Interval (deg)-----------------------------------------------------
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
   StandardEarth.SetExposure(MT);
   StandardEarth.flvf=nuflv;

    //Alternative Earth-------------------------------------------------------------
   AsimovObsSimulation AlternativeEarth;

   AlternativeEarth.ThePremTable = path2prem;
   AlternativeEarth.TheHondaTable = FluxTable;
   AlternativeEarth.SetDetectorXYZ(rdet);
   //AlternativeEarth.SetEnergyResolution(0.001,0.0);
   //AlternativeEarth.SetAngularResolution(0.001,0.0);
   
   //Set LLVP
   AlternativeEarth.PileInModel = true;
   AlternativeEarth.ShapeOfPile = shape;
   AlternativeEarth.ThePileHight = PileHeight;
   AlternativeEarth.ThePileAperture = 45;
   AlternativeEarth.ThePileDensityContrast = PileDensityPct;
   AlternativeEarth.ThePileChemicalContrast = 0.0;
   
   //Simulation of Alternative Earth
   AlternativeEarth.SetIntervals(thmin,thmax,phimin,phimax,Emin,Emax);
   AlternativeEarth.SetTrueBinning(cthtruebins,atruebins,etruebins);
   AlternativeEarth.SetRecoBinning(cthrecobins,arecobins,erecobins);
   AlternativeEarth.SetExposure(MT);
   AlternativeEarth.flvf=nuflv;


   std::vector<TH2D*>  ObsStd = StandardEarth.GetObsEvents3Dcth();
   std::cout << "Generating Observed events" << std::endl;
   std::vector<TH2D*>  ObsAlt = AlternativeEarth.GetObsEvents3Dcth();

   
   
   //Sensitivity
    std::string NuTomoPath= "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation";
    std::string SenvFolder = "/SimulationResults/PreliminaryResults/ObsChi2/";
    std::string BinLabeltrue = std::to_string(cthtruebins)+std::to_string(atruebins)+std::to_string(etruebins);
    std::string BinLabelreco = std::to_string(cthrecobins)+std::to_string(arecobins)+std::to_string(erecobins);
    //std::string SimLabel = std::to_string(thmin)+std::to_string(thmax)+std::to_string(EnuMin)+std::to_string(EnuMax);
    std::string chi2nameth  = "ObsChi2_cth_LLVP_"+shape+"nu"+std::to_string(nuflv)+"_"+BinLabelreco+"_"+BinLabeltrue+".txt";
    std::string chi2pathth  = NuTomoPath+SenvFolder+chi2nameth;

    
    std::ofstream SenvData(chi2pathth); 


    for (int rpct = -4; rpct <= 4; ++ rpct)
    {

        AlternativeEarth.ThePileDensityContrast = rpct; //Adjust the LLVP contrats density

        std::vector< TH2D* > NewAlt= AlternativeEarth.GetObsEvents3Dcth();

        double chi2totth = 0;

        for (int n = 0; n < ObsAlt.size(); ++n)
        {
            chi2totth += Get2DChi2( ObsStd[n] , NewAlt[n]);
        }

    //std::cout << "Density %: " << rpct << " Chi2: " << chi2tot << std::endl;

        SenvData << rpct << " " << chi2totth << " " << cthrecobins << " " << arecobins << " " << erecobins << std::endl; 

        
    }

    SenvData.close();
   
  

   std::cout << "DATA DISPLAY" << std::endl;


   TApplication app("app", &argc, argv);


   //double cthbottom = -0.8376; 
   double thbottom     = cos( TMath::Pi()-TMath::ASin( 3480.0/Rdet ) );
   double thmid_bottom = cos(TMath::Pi()-TMath::ASin( (3480.0+300)/Rdet ));;
   double thmid_top    = cos(TMath::Pi()-TMath::ASin( (3480.0+600)/Rdet ));;
   double thtop        = cos(TMath::Pi()-TMath::ASin( (3480.0 + 1000)/Rdet ));;

   TLine * lbottomth = new TLine(thbottom,Emin,thbottom,Emax);
   TLine * lmid_bottomth = new TLine(thmid_bottom,Emin,thmid_bottom,Emax);
   TLine * lmid_topth = new TLine(thmid_top,Emin,thmid_top,Emax);
   TLine * ltopth = new TLine(thtop,Emin,thtop,Emax);

   lbottomth -> SetLineColor(kRed);
   lmid_bottomth -> SetLineColor(kRed);
   lmid_topth -> SetLineColor(kRed);
   ltopth -> SetLineColor(kRed);

   lbottomth -> SetLineStyle(2);
   lmid_bottomth -> SetLineStyle(2);
   lmid_topth -> SetLineStyle(2);
   ltopth -> SetLineStyle(2);

   TCanvas *c1 = new TCanvas();

   gStyle->SetPalette(kBlueGreenYellow);

   c1->Divide(5,5);

   for(int i = 0 ; i < 25 ; ++i )

   {

    int nhist = (51 + 2*i)-1;
    c1->cd(i+1);

    TH2D * ObsDiff2D = new TH2D(Form("ObsDiff2D%d",nhist),Form("OscObs%d",nhist), cthrecobins,cthmin,cthmax,erecobins,Emin,Emax); //binning in cth 
    GetDiff2D( ObsStd[nhist] , ObsAlt[nhist], ObsDiff2D );
    ObsDiff2D->Draw("COLZ");
    ObsDiff2D->SetStats(0);
    lbottomth->Draw("same");
    lmid_bottomth->Draw("same");
    lmid_topth->Draw("same");
    ltopth->Draw("same");

    std::string BinLabel = std::to_string(cthrecobins)+std::to_string(arecobins)+std::to_string(erecobins);
    std::string eventsfile = "Obs_cth_LLVP_"+shape+"nu"+std::to_string(nuflv)+"_"+BinLabelreco+"_"+BinLabeltrue+"_"+std::to_string(nhist)+".txt";
    ExportToCSV(ObsDiff2D, eventsfile);

 
   }

    gPad->Update();
    c1->Modified(); c1->Update();
    TRootCanvas *rc = (TRootCanvas *)c1->GetCanvasImp();
    rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
    app.Run();

    return 0;
   


}

void ExportToCSV(TH2D* hist, std::string filename)
{
    std::string NuTomoPath= "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation";
    std::string IntFolder = "/SimulationResults/PreliminaryResults/ObsEvents/";
    std::string path2file  = NuTomoPath+IntFolder+filename;
   
    std::ofstream outfile(path2file);

    // Recorre los bins y guarda el contenido
    for (int i = 1; i <= hist->GetNbinsX(); ++i)
    {
        for (int j = 1; j <= hist->GetNbinsY(); ++j)
        {
            double x = hist->GetXaxis()->GetBinCenter(i);
            double y = hist->GetYaxis()->GetBinCenter(j);
            double content = hist->GetBinContent(i, j);

            outfile << x << "," << y << "," << content << std::endl;
        }
    }

    outfile.close();
}



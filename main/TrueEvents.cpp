
/* Direction of neutrinos a re given in degrees*/
#define _USE_MATH_DEFINES
#include <cmath>

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <sys/stat.h>


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

void ExportToCSV(TH2D* stdhist, TH2D* althist, TH2D* hist, std::string filename);

int main(int argc, char **argv)
{
    // Plot and Histogram2D settings

    std::cout << " Neutrino Oscillation tomography. " << std::endl;

    //Detetcor Settings
    
    //Detector size-------------------------------------------------------------
    int exposue = 10*10;
    double M = 10.0*MTon; //Mass in megaton units
    double T       = 10.0*years2sec; //Detector Exposure time in sec: One Year
    double MT = M*T; // Exposure [Mton*years]

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
    std::string shape = "pancake";

    //SIMULATION SET UP---------------------------------------------------------

    // Binning------------------------------------------------------------------
    int cthbins=50 ; // # Bins in zenith/cos(zenith)
    int abins =100; // # Bins in azimuth (optimal bins are 110 or 22)
    int ebins =50; // bins in energy

    //Energy interval (in GeV)--------------------------------------------------
    double EnuMin=1.0 ; 
    double EnuMax=10.0 ;

    //Zenith Angle Interval-----------------------------------------------------
    double thmin = 180-TMath::ASin( DepthMax/Rdet )*(180.0/TMath::Pi()) ; // min 90
    double thmax = 180-TMath::ASin( (DepthMin)/Rdet )*(180.0/TMath::Pi()) ; // max 180

    double cthmin = cos(thmax*TMath::Pi()/180.0);// Cos(thmin) - Minimum possible is -1    
    double cthmax = cos(thmin*TMath::Pi()/180.0);// Cos(thmax) - Maximum possible is 0

    //Azimuthal Interal---------------------------------------------------------
    double phi_a = -50;
    double phi_b =  50;
    double dPhi = (phi_b-phi_a)/(2.0*abins);
    
    // for bins centered
    double phimin = phi_a - dPhi; 
    double phimax = phi_b + dPhi;

    //Neutrino
    int nuflv = 1; // neutrino  final state: nue (0), numu (1) or nutau (2)
    
   //Standart Earth-------------------------------------------------------------
   AsimovSimulation StandardEarth;

   StandardEarth.PremTable = path2prem;
   StandardEarth.HondaTable = FluxTable;
   StandardEarth.SetDetectorPosition(rdet);
   StandardEarth.MantleAnomaly = false;
   
   //Simulation of Standard Earth
   StandardEarth.SetIntervals(thmin,thmax,phimin,phimax,EnuMin,EnuMax);
   StandardEarth.SetBinning(cthbins,abins,ebins);
   StandardEarth.SetExposure(MT);
   StandardEarth.flvf=nuflv;

   //TH2D * TrueStd = StandardEarth.GetTrueEvents2D();

    std::vector< TH2D* >  TrueStd = StandardEarth.GetTrueEvents3D();
    std::vector< TH2D* >  TrueStdth = StandardEarth.GetTrueEvents3Dth();

    //Alternative Earth--------------------------------------------------------- 

   AsimovSimulation AlternativeEarth;

   AlternativeEarth.PremTable = path2prem;
   AlternativeEarth.HondaTable = FluxTable;
   AlternativeEarth.SetDetectorPosition(rdet);

   // Construct LLVPs
   AlternativeEarth.MantleAnomaly = true;
   AlternativeEarth.PileHeight = PileHeight;
   AlternativeEarth.aperture=45;
   AlternativeEarth.AnomalyShape=shape;
   AlternativeEarth.PileDensityContrast = PileDensityPct;
   AlternativeEarth.PileChemContrast = 0.0;

   //Simulation of Standard Earth
   AlternativeEarth.SetIntervals(thmin,thmax,phimin,phimax,EnuMin,EnuMax);
   AlternativeEarth.SetBinning(cthbins,abins,ebins);
   AlternativeEarth.SetExposure(MT);
   AlternativeEarth.flvf=nuflv;

   //TH2D * TrueAlt = AlternativeEarth.GetTrueEvents2D();
   std::vector< TH2D* > TrueAlt= AlternativeEarth.GetTrueEvents3D(); //cth binning
    std::vector< TH2D* > TrueAltth= AlternativeEarth.GetTrueEvents3Dth(); //th binning



   
   //Sensitivity
    std::string NuTomoPath= "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation";
    std::string SenvFolder = "/SimulationResults/PreliminaryResults/IntChi2/";
    std::string BinLabel = std::to_string(cthbins)+std::to_string(abins)+std::to_string(ebins);
    //std::string SimLabel = std::to_string(thmin)+std::to_string(thmax)+std::to_string(EnuMin)+std::to_string(EnuMax);
    std::string chi2namecth  = "Int_cth_Chi2LLVP"+std::to_string(exposue)+"myrs_"+shape+"nu"+std::to_string(nuflv)+"_"+BinLabel+".txt";
    std::string chi2nameth  = "Int_th_Chi2LLVP"+std::to_string(exposue)+"myrs_"+shape+"nu"+std::to_string(nuflv)+"_"+BinLabel+".txt";
    std::string chi2pathcth  = NuTomoPath+SenvFolder+chi2namecth;
    std::string chi2pathth  = NuTomoPath+SenvFolder+chi2nameth;


    std::ofstream SenvDatacth(chi2pathcth);
    std::ofstream SenvDatath(chi2pathth); 

    
    for (int rpct = -4; rpct <= 4; ++ rpct)
    {

        AlternativeEarth.PileDensityContrast = rpct; //Adjust the LLVP contrats density

        std::vector< TH2D* > NewAltcth= AlternativeEarth.GetTrueEvents3D();
        std::vector< TH2D* > NewAltth= AlternativeEarth.GetTrueEvents3Dth();

        double chi2totcth = 0;
        double chi2totth = 0;

        for (int n = 0; n < TrueAlt.size(); ++n)
        {
            chi2totcth += Get2DChi2( TrueStd[n] , NewAltcth[n]);
            chi2totth += Get2DChi2( TrueStdth[n] , NewAltth[n]);
        }

    //std::cout << "Density %: " << rpct << " Chi2: " << chi2tot << std::endl;

        SenvDatacth << rpct << " " << chi2totcth << " " << cthbins << " " << abins << " " << ebins << std::endl;
        SenvDatath << rpct << " " << chi2totth << " " << cthbins << " " << abins << " " << ebins << std::endl;  

        
    }

    SenvDatacth.close();
    SenvDatath.close();
   

   // Visualization of Events

   TApplication app("app", &argc, argv);

   //cth lines
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

   lbottom -> SetLineColor(kRed);
   lmid_bottom -> SetLineColor(kRed);
   lmid_top -> SetLineColor(kRed);
   ltop -> SetLineColor(kRed);

   lbottom -> SetLineStyle(2);
   lmid_bottom -> SetLineStyle(2);
   lmid_top -> SetLineStyle(2);
   ltop -> SetLineStyle(2);

   //th lines 
   double thbottom     = (180.0/TMath::Pi())*(TMath::Pi()-TMath::ASin( 3480.0/Rdet ));
   double thmid_bottom = (180.0/TMath::Pi())*(TMath::Pi()-TMath::ASin( (3480.0 + 300)/Rdet ));
   double thmid_top    = (180.0/TMath::Pi())*(TMath::Pi()-TMath::ASin( (3480.0 + 600)/Rdet ));
   double thtop        = (180.0/TMath::Pi())*(TMath::Pi()-TMath::ASin( (3480.0 + 1000)/Rdet ));

   TLine * lbottomth = new TLine(thbottom,EnuMin,thbottom,EnuMax);
   TLine * lmid_bottomth = new TLine(thmid_bottom,EnuMin,thmid_bottom,EnuMax);
   TLine * lmid_topth = new TLine(thmid_top,EnuMin,thmid_top,EnuMax);
   TLine * ltopth = new TLine(thtop,EnuMin,thtop,EnuMax);

   lbottomth -> SetLineColor(kBlue);
   lmid_bottomth -> SetLineColor(kBlue);
   lmid_topth -> SetLineColor(kBlue);
   ltopth -> SetLineColor(kBlue);

 
   TCanvas *c1 = new TCanvas();
   gStyle->SetPalette(kBlueGreenYellow);
   
   TCanvas *c2 = new TCanvas();
   gStyle->SetPalette(kBlueGreenYellow);


   c1->Divide(5,5);
   c2->Divide(5,5);

   for(int i = 0 ; i < 25 ; ++i )

   {

    int nhist = (51 + 2*i)-1;
    

    std::string BinLabel = std::to_string(cthbins)+std::to_string(abins)+std::to_string(ebins);
    std::string eventsfile = "IntcthLLVP_"+shape+"nu"+std::to_string(nuflv)+"_"+BinLabel+"_"+std::to_string(nhist)+".txt";
    std::string eventsfileth = "IntthLLVP_"+shape+"nu"+std::to_string(nuflv)+"_"+BinLabel+"_"+std::to_string(nhist)+".txt";


    TH2D * TrueDiff2D = new TH2D(Form("OscDiff2D%d",nhist),Form("Oscillogram%d",nhist), cthbins,cthmin,cthmax,ebins,EnuMin,EnuMax); //binning in cth 
    TH2D * TrueDiff2Dth = new TH2D(Form("OscDiff2Dth%d",nhist),Form("Oscillogramth%d",nhist), cthbins,thmin,thmax,ebins,EnuMin,EnuMax); //binning in th 
    
    GetDiff2D( TrueStd[nhist] , TrueAlt[nhist], TrueDiff2D);
    ExportToCSV(TrueStd[nhist] ,TrueAlt[nhist], TrueDiff2D, eventsfile);

    
    GetDiff2D( TrueStdth[nhist] , TrueAltth[nhist], TrueDiff2Dth);
    ExportToCSV(TrueStdth[nhist] , TrueAltth[nhist], TrueDiff2Dth, eventsfileth);
    
    c1->cd(i+1);
    TrueDiff2D->Draw("COLZ");
    TrueDiff2D->SetStats(0);
    lbottom->Draw("same");
    lmid_bottom->Draw("same");
    lmid_top->Draw("same");
    ltop->Draw("same");
    

    c2->cd(i+1);
    TrueDiff2Dth->Draw("COLZ");
    TrueDiff2Dth->SetStats(0);
    lbottomth->Draw("same");
    lmid_bottomth->Draw("same");
    lmid_topth->Draw("same");
    ltopth->Draw("same");


    //std::string BinLabel = std::to_string(cthbins)+std::to_string(abins)+std::to_string(ebins);
    //std::string eventsfile = "IntLLVP_"+shape+"nu"+std::to_string(nuflv)+"_"+BinLabel+"_"+std::to_string(nhist)+".txt";

 
   }

   gPad->Update();
   c1->Modified(); c1->Update();
   TRootCanvas *rc = (TRootCanvas *)c1->GetCanvasImp();
   rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");

   gPad->Update();
   c2->Modified(); c2->Update();
   TRootCanvas *rc2 = (TRootCanvas *)c2->GetCanvasImp();
   rc2->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");

   app.Run();

    //c1->Print("AzimuthHist.png");

    return 0;

}


void ExportToCSV(TH2D* stdhist, TH2D* althist, TH2D* hist, std::string filename)
{
    std::string NuTomoPath= "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation";
    std::string IntFolder = "/SimulationResults/PreliminaryResults/IntEvents/";
    std::string path2file  = NuTomoPath+IntFolder+filename;
   
    std::ofstream outfile(path2file);

    // Recorre los bins y guarda el contenido


    for (int i = 1; i <= hist->GetNbinsX(); ++i)
    {
        for (int j = 1; j <= hist->GetNbinsY(); ++j)
        {
            double x = hist->GetXaxis()->GetBinCenter(i);
            double y = hist->GetYaxis()->GetBinCenter(j);

            double nexp = stdhist->GetBinContent(i, j);
            double nobs = althist->GetBinContent(i, j); 
            double diff = hist->GetBinContent(i, j);

            outfile << x << "," << y << "," <<  nexp <<" , " << nobs << " , " << diff << std::endl;
        }
    }

    outfile.close();
}


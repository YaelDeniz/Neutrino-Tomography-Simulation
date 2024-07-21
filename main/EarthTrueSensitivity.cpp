
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
#include "TGraph.h"
#include "TMultiGraph.h"
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

//Some import Earth Radius
# define Rcmb 3480.0 // Core-Mantle-Boundary Radius (km)
# define Rocean 6368.0 // Earth radius in (km)
# define Rearth 6371.0 // Earth radius in (km)
# define Ratm 6386.0 // Radius of the Atmosphere (km)

//Some units
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

    //std
    std::string PremFolder ="/home/dehy0499/OscProb/PremTables/";
    std::string PremName = "prem_44layers.txt";
    std::string path2prem = PremFolder+PremName;

    //Modified
    std::string PremFolderAlt ="/home/dehy0499/OscProb/PremTables/Prem44pct10/";
    

    //Interval of integration in Zenith, Azimuth and E for Neutrino Events

    //Energy interval (in GeV)--------------------------------------------
    double EnuMin=1.0 ; 
    double EnuMax=20.0 ;

    //Zenith Angle Interval-----------------------------------------------------------------------------------
    double zenmin = 90.01; // Not Able to use 90 degrees
    double zenmax = 180.0;

    double czmin = cos(zenmax*TMath::Pi()/180.0);
    double czmax = cos(zenmin*TMath::Pi()/180.0);

    //Azimuthal Interal-------------------------------------------------------------------------
    double phimin = 0.0;
    double phimax = 360.0 ;
    
    // Number of bins
    int czbins=50; // # Bins in zenith/cos(zenith)
    int abins=50 ; // # Bins in azimuth
    int ebins=50; // bins in energy

    //double pct = 10;


    /*
    std::vector<double> chi2data;
    std::string ChiFolder = "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/SimulationResults/chi2results/chi2true/";
    std::ChiName = "C2T"+to_string(int(EnuMin))+std::to_string(int(EnuMax))_to_string(int(zenmin))+std::to_string(int(zenmax))_to_string(int(czbins))+std::to_string(int(ebins))+".csv";
    std::string chi2file = ChiFolder + ChiName;
    double chi2 = 0;
  */


    //Event generation-----------------------------------------------------------------------------------------
    int nuflv = 1; // neutrino  final state options nue (0), numu (1) or nutau (2)

    //ofstream EarthChi2(chi2file, std::ofstream::trunc); //Opens a file and rewrite content, if files does not exist it Creates new file


    TMultiGraph *Pnu = new TMultiGraph(); // Oscillations to muon neutrinos
    TMultiGraph *Pnubar = new TMultiGraph(); // Oscillation to electron neutrinos

    // Standard Earth model

    AsimovSimulation StandardEarth;
    StandardEarth.SetDetectorPosition(rdet);
    StandardEarth.PremTable = path2prem;
    StandardEarth.HondaTable = FluxTable;
    //StandardEarth.MantleAnomaly = false;
    StandardEarth.SetIntervals(zenmin,zenmax,phimin,phimax,EnuMin,EnuMax);
    StandardEarth.SetBinning(czbins,abins,ebins);
    StandardEarth.SetExposure(NnT);
    StandardEarth.flvf=nuflv;

    TH2D * TrueStd = StandardEarth.GetTrueEvents2D();

    // Neutrino final flavour
    int nue        = 0;  // electron neutrino  
    int numu       = 1; // muon neutrino
  
    bool  nunubar      = false; //antineutrino
    //double cth  = -0.827;//correspond to layer 25
    double cth  = -0.885;//correspond to layer 18

    //THStack * OscHist = new THStack("OscHist","Oscillograms ");
    TH2D * OscProb_mumu = StandardEarth.GetOscProb2D(numu,numu, nunubar); 
    TH2D * OscProb_emu = StandardEarth.GetOscProb2D(nue,numu, nunubar); 

    TGraph * OscPlot_mumu =  StandardEarth.GetOscProb( numu, numu, nunubar, cth); //To be Deleted
    OscPlot_mumu->SetLineColor(kBlue);
    OscPlot_mumu->SetLineWidth(5);
    Pnu->Add(OscPlot_mumu);

    TGraph * OscPlot_emu =  StandardEarth.GetOscProb( nue, numu, nunubar, cth ); //To be Deleted
    OscPlot_emu->SetLineColor(kRed);
    OscPlot_emu->SetLineWidth(5);
    Pnu->Add(OscPlot_emu);
    //TH2D * TrueStd = StandardEarth.TestTrueEvents2D("/home/dehy0499/OscProb/PremTables/prem_default.txt","/home/dehy0499/OscProb/PremTables/prem_default.txt");

   

    // Alternative Earth Model

    //std::string PremNameAlt = "prem_44layers_18_10pct.txt";
    //std::string path2premAlt = PremFolderAlt+PremNameAlt;

    std::string path2premAlt = PremFolder+PremName; // std prem
    AsimovSimulation AlternativeEarth;

    TGraph * chi2plot = new TGraph(44);

    for (int i = 0; i < 44; ++i)
    {
        
    int layer = i+1;

    AlternativeEarth.SetDetectorPosition(rdet);
    AlternativeEarth.PremTable = path2premAlt;
    AlternativeEarth.HondaTable = FluxTable;
    //AlternativeEarth.MantleAnomaly = false;
    //AlternativeEarth.AnomalyShape="pancake";
    AlternativeEarth.ModifyLayer(layer,10,0.0);
    AlternativeEarth.SetIntervals(zenmin,zenmax,phimin,phimax,EnuMin,EnuMax);
    AlternativeEarth.SetBinning(czbins,abins,ebins);
    AlternativeEarth.SetExposure(NnT);
    AlternativeEarth.flvf=nuflv;

    TH2D * TrueAlt = AlternativeEarth.GetTrueEvents2D();

    double chi2 = Get2DChi2(TrueStd,TrueAlt);

    chi2plot->SetPoint(i-1,chi2,i);


    }


    AlternativeEarth.PremTable = path2premAlt;
    AlternativeEarth.ModifyLayer(18,10,0.0);

    //TH2D * OscProbAlt = AlternativeEarth.GetOscProb2D(nue,numu, nubar); 

    TGraph * OscPlotAlt_mumu =  AlternativeEarth.GetOscProb( numu, numu, nunubar , cth ); //To be Deleted
    OscPlotAlt_mumu->SetLineColor(kBlue);
    OscPlotAlt_mumu->SetLineWidth(5);
    OscPlotAlt_mumu->SetLineStyle(9);
    Pnu->Add(OscPlotAlt_mumu);
    TGraph * OscPlotAlt_emu =  AlternativeEarth.GetOscProb( nue, numu, nunubar , cth ); //To be Deleted
    OscPlotAlt_emu->SetLineColor(kRed);
    OscPlotAlt_emu->SetLineWidth(5);
    OscPlotAlt_emu->SetLineStyle(9);
    Pnu->Add(OscPlotAlt_emu);

    //TH2D * TrueAlt = AlternativeEarth.GetTrueEvents2D();
    //TH2D * TrueAlt= AlternativeEarth.SensitivityTrueEvents2D( i, pct[j] );
    //TH2D * TrueAlt = AlternativeEarth.TestTrueEvents2D("/home/dehy0499/OscProb/PremTables/prem_default.txt","/home/dehy0499/OscProb/PremTables/prem_llsvp.txt")

  // Oscillation Probability
  
    TApplication app("app", &argc, argv);

    TCanvas *c1 = new TCanvas();
    c1->cd();

    Pnu->Draw("apl");
    //Pnu->GetXaxis()->SetNdivisions(5,kFALSE);
    gPad->Update();
    c1->Modified(); 
    c1->Update();
    TRootCanvas *rc1 = (TRootCanvas *)c1->GetCanvasImp();
    rc1->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");


    //Oscillation inside standard Earth

    TCanvas *c2 = new TCanvas();
    c2->Divide(1,2,0.01,0.01);
    c2->cd(1);
  
    gStyle->SetPalette(kBlueGreenYellow);
    OscProb_mumu->Draw("COLZ");
    OscProb_mumu->SetStats(0);
 
    c2->cd(2);
    
    gStyle->SetPalette(kBlueGreenYellow);
    OscProb_emu->Draw("COLZ");
    OscProb_emu->SetStats(0);

    gPad->Update();
    c2->Modified(); 
    c2->Update();
    TRootCanvas *rc2 = (TRootCanvas *)c2->GetCanvasImp();
    rc2->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");



    TCanvas *c3 = new TCanvas();
    c3->cd();
    
    chi2plot->Draw();
    gPad->SetLogx();

    gPad->Update();
    c3->Modified(); 
    c3->Update();
    TRootCanvas *rc3 = (TRootCanvas *)c3->GetCanvasImp();
    rc3->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
    
    


    app.Run();

    std::cout << PremName << std::endl;



  return 0;

}

void GetDiff3D( TH3D * histstd , TH3D * histalt, TH3D * diff)
{


   std::ofstream TrueDiffFile("SimulationResults/TrueEventsResults/dN_true_3D.csv"); 

   double cth, azi, e, Nexp, Nobs, dN;


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

                TrueDiffFile << cth << " , " << azi << " , " << e << " , " << dN << std::endl;

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

    for (int j = 1; j <= histstd->GetYaxis()->GetNbins(); ++j) //loop in energy
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
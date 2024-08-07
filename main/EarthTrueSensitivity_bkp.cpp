
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
#include "TLine.h"
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
    double EnuMax=25.0 ;

    //Zenith Angle Interval-----------------------------------------------------------------------------------
    double zenmin = 90.01; // Not Able to use 90 degrees
    double zenmax = 180.0;

    double czmin = cos(zenmax*TMath::Pi()/180.0);
    double czmax = cos(zenmin*TMath::Pi()/180.0);

    //Azimuthal Interal-------------------------------------------------------------------------
    double phimin = 0.0;
    double phimax = 360.0 ;
    
    // Number of bins
    int czbins=100; // # Bins in zenith/cos(zenith)
    int abins= 50 ; // # Bins in azimuth
    int ebins= 100; // bins in energy

    //double pct = 10;


   


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
    StandardEarth.SetIntervals(zenmin,zenmax,phimin,phimax,EnuMin,EnuMax);
    StandardEarth.SetBinning(czbins,abins,ebins);
    StandardEarth.SetExposure(NnT);
    
    StandardEarth.flvf=1;
    StandardEarth.nunubar = 1;
    TH2D * TrueStd_mu = StandardEarth.GetTrueEvents2D();

    //StandardEarth.flvf=1;
    //StandardEarth.nunubar = -1;
    //TH2D * TrueStd_mub = StandardEarth.GetTrueEvents2D();

    StandardEarth.flvf=0;
    StandardEarth.nunubar = 1;
    TH2D * TrueStd_e = StandardEarth.GetTrueEvents2D();

    //StandardEarth.flvf=0;
    //StandardEarth.nunubar = -1;
    //TH2D * TrueStd_eb = StandardEarth.GetTrueEvents2D();





    // Neutrino final flavour
    int nue        = 0;  // electron neutrino  
    int numu       = 1; // muon neutrino
  
    bool  nunubar      = false; //antineutrino
    //double cth  = -0.827;//correspond to layer 25
    double cth  = -0.885;//correspond to layer 18 in p44, 47 in p64

    //THStack * OscHist = new THStack("OscHist","Oscillograms ");
    
    //TH2D * OscProb_mumu = StandardEarth.GetOscProb2D(numu,numu, nunubar); 
    //TH2D * OscProb_emu = StandardEarth.GetOscProb2D(nue,numu, nunubar); 

    TGraph * OscPlot_mumu =  StandardEarth.GetOscProb( numu, numu, nunubar, cth); //To be Deleted
    OscPlot_mumu->SetLineColor(kBlue);
    OscPlot_mumu->SetLineWidth(5);
    Pnu->Add(OscPlot_mumu);

    TGraph * OscPlot_emu =  StandardEarth.GetOscProb( nue, numu, nunubar, cth ); //To be Deleted
    OscPlot_emu->SetLineColor(kRed);
    OscPlot_emu->SetLineWidth(5);
    Pnu->Add(OscPlot_emu);
    //TH2D * TrueStd = StandardEarth.TestTrueEvents2D("/home/dehy0499/OscProb/PremTables/prem_default.txt","/home/dehy0499/OscProb/PremTables/prem_default.txt");

   

    
    std::vector<double> chi2data;
    std::string ChiFolder = "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/SimulationResults/chi2results/chi2true/";
    std::string ChiName = "chi2test.csv";
    std::string chi2path = ChiFolder + ChiName;
    
    std::ofstream chi2file(chi2path);


     
    //std::string path2premAlt = PremFolder+PremName; // std prem
    AsimovSimulation AlternativeEarth;

    TGraph * chi2plot = new TGraph(44);

    for (int i = 0; i < 44; ++i)
    {
        
    int layer = i+1;

    //PremName = "Prem44pct10/prem_44layers_"+std::to_string(layer) +"_10pct.txt";
    //std::string path2premAlt = PremFolder+PremName; // std prem

    AlternativeEarth.SetDetectorPosition(rdet);
    AlternativeEarth.PremTable = path2prem;
    AlternativeEarth.HondaTable = FluxTable;
    AlternativeEarth.MantleAnomaly = false;
    AlternativeEarth.AnomalyShape="pancake";
    AlternativeEarth.ModifyLayer(layer,10,0.0);
    AlternativeEarth.SetIntervals(zenmin,zenmax,phimin,phimax,EnuMin,EnuMax);
    AlternativeEarth.SetBinning(czbins,abins,ebins);
    AlternativeEarth.SetExposure(NnT);
    
    AlternativeEarth.flvf=1;
    AlternativeEarth.nunubar = 1;
    TH2D * TrueAlt_mu = StandardEarth.GetTrueEvents2D();
    double chi2mu = Get2DChi2(TrueStd_mu,TrueAlt_mu);

    //AlternativeEarth.flvf=1;
    //AlternativeEarth.nunubar = -1;
    //TH2D * TrueAlt_mub = StandardEarth.GetTrueEvents2D();
    //double chi2mub = Get2DChi2(TrueStd_mub,TrueAlt_mub);
    double chi2mub = 0;
    
    AlternativeEarth.flvf=0;
    AlternativeEarth.nunubar = 1;
    TH2D * TrueAlt_e = StandardEarth.GetTrueEvents2D();
    double chi2e = Get2DChi2(TrueStd_e,TrueAlt_e);
    //double chi2e=0;

    //StandardEarth.flvf=0;
    //AlternativeEarth.nunubar = -1;
    //TH2D * TrueAlt_eb = StandardEarth.GetTrueEvents2D();
    //double chi2eb = Get2DChi2(TrueStd_eb,TrueAlt_eb);
    double chi2eb=0;

    double chi2total = chi2e+chi2eb+chi2mu+chi2mub;

    chi2file << std::setprecision(10) << i+1 << " " << chi2e << " " << chi2eb << " " << chi2mu << " " << chi2mub << " " << chi2total << std::endl;

    chi2plot->SetPoint(i,chi2total,i+1);


    }

    chi2file.close();

    


    //AlternativeEarth.PremTable = path2premAlt;
    //AlternativeEarth.ModifyLayer(18,10,0.0);


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



  // Oscillation Probability
    TApplication app("app", &argc, argv);
    TCanvas *c1 = new TCanvas();
    c1->cd();
    Pnu->Draw("apl");
    gPad->Update();
    c1->Modified(); 
    c1->Update();
    TRootCanvas *rc1 = (TRootCanvas *)c1->GetCanvasImp();
    rc1->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");

    //Oscillation inside standard Earth
    TCanvas *c2 = new TCanvas();

    c2->Divide(2,2,0.01,0.01);
    c2->cd(1);
    gStyle->SetPalette(kBlueGreenYellow);
    TrueStd_mu->Draw("COLZ");
    TrueStd_mu->SetStats(0);
    gPad->SetLogy();
    gPad->Update();

    c2->cd(2);
 //   gStyle->SetPalette(kBlueGreenYellow);
 //   TrueStd_mub->Draw("COLZ");
 //   TrueStd_mub->SetStats(0);
 //   gPad->SetLogy();
 //   gPad->Update();

    c2->cd(3);
    gStyle->SetPalette(kBlueGreenYellow);
    TrueStd_e->Draw("COLZ");
    TrueStd_e->SetStats(0);
    gPad->SetLogy();
    gPad->Update();

//    c2->cd(4);
//    gStyle->SetPalette(kBlueGreenYellow);
//    TrueStd_eb->Draw("COLZ");
//    TrueStd_eb->SetStats(0);
//    gPad->SetLogy();
//    gPad->Update();
    
    c2->Modified(); 
    c2->Update();
    TRootCanvas *rc2 = (TRootCanvas *)c2->GetCanvasImp();



/*
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
*/

    //Oscillation inside standard Earth


    


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
                //chi2 +=  2*( Nexp - Nobs + Nobs*TMath::Log(Nobs/Nexp) ); // LLRT

                chi2 +=((Nexp-Nobs)*(Nexp-Nobs))/Nexp; 

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
                //chi2 +=  2*( Nexp - Nobs + Nobs*TMath::Log(Nobs/Nexp) ); // LLRT
                chi2 +=((Nexp-Nobs)*(Nexp-Nobs))/Nexp; 

            
        }
        
    }

    return chi2;

}
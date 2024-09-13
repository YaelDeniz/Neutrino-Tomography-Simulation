
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

//Some important Earth Radius
# define Rcmb 3480.0 // Core-Mantle-Boundary Radius (km)
# define Rocean 6368.0 // Earth radius in (km)
# define Rearth 6371.0 // Earth radius in (km)
# define Ratm 6386.0 // Radius of the Atmosphere (km)

//Some units
# define mN   1.67492749804E-27  // Nucleons mass (~ Neutron mass) in kg
# define MTon  1E9  //Metric MegaTon in kg
# define years2sec 3.154E7 // Years in Seconds

using namespace std;

void ExportToCSV(TH2D* hist, std::string filename); //Export File to CSV format


int main(int argc, char **argv)
{
    //Setup the detector
    double R =  Rocean;                // Underwater detector
    double rdet[3] = {0.0,0.0, -1*R }; // Detector location
    double DetMass = 10.0*MTon;        // Mass in megaton units
    double T       = 20.0*years2sec;   // Detector Exposure time in sec: One Year
    double NT = DetMass*T;             // Exposure [Mton*years]

    //Set the neutrino Neutrino Flux
    std::string FluxFolder = "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/NuFlux/";
    std::string FluxFile = "GRN_AziAveraged_solmin/grn-nu-20-01-000.d"; //Gran Sasso
    std::string FluxTable = FluxFolder + FluxFile; 

    //Set the Earth models
    std::string PremFolder ="/home/dehy0499/OscProb/PremTables/";
    std::string PremName = "prem_44layers";
    std::string PremFile = PremName+".txt";
    std::string path2prem = PremFolder+PremFile;

    //Set Layer Density contrats
    int rhopct = 10;

    //int MaxLayers = 44;
    //std::string PremFolderAlt ="/home/dehy0499/OscProb/PremTables/Prem44pct10/";

    //Simulation Setup

    //Energy interval (in GeV)
    double Emin=1.0 ; 
    double Emax=40.0 ;

    //Zenith Angle Interval
    double thmin = 90.01; // Not Able to use 90 degrees
    double thmax = 180.0;
    double cthmin = cos(thmax*TMath::Pi()/180.0);
    double cthmax = cos(thmin*TMath::Pi()/180.0);

    //Azimuthal angle Interval
    double phimin = 0.0;
    double phimax = 360.0 ; //Whole Earth
    
    // Number of bins
    int cthbins=100; // # Bins in zenith/cos(zenith)
    int abins = 1;
    int ebins= 100; // bins in energy

    //TMultiGraph *Pnu = new TMultiGraph(); // Oscillations to muon neutrinos
    //TMultiGraph *Pnubar = new TMultiGraph(); // Oscillation to electron neutrinos

    // Signal from a Standard Earth model
    AsimovSimulation TrueEvents;
    TrueEvents.HondaTable = FluxTable;
    TrueEvents.SetIntervals(thmin,thmax,phimin,phimax,Emin,Emax);
    TrueEvents.SetBinning(cthbins,abins,ebins);
    TrueEvents.SetDetectorPosition(rdet);
    TrueEvents.SetExposure(NT);
    TrueEvents.PremTable = path2prem;

    //Muon neutrinio like events
    TrueEvents.flvf=1;
    std::vector< TH2D * > TrueStd_mu = TrueEvents.GetTrueEvents2D();

    // Electron netrino like events
    TrueEvents.flvf=0;
    std::vector< TH2D * > TrueStd_e = TrueEvents.GetTrueEvents2D();

    /* Export to CSV */
    //Expected Interacting Events for a Standard Earth.
    std::string ResultFolder = "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/SimulationResults/";
    std::string IntEvntsFolder = "/TrueEventsResults/IntStdEarth/";
    
    for (int nhist = 0; nhist < TrueStd_mu.size(); nhist++)
    {
        std::string  muname= "IntStdEarth"+std::to_string(1)+PremName+std::to_string(rhopct)+"_"+std::to_string(cthbins)+std::to_string(abins)+".txt"
        std::string mufile = ResultFolder + IntEvntsFolder  + muname;
        ExportToCSV(TrueStd_mu[nhist],mufile);

        std::string  ename= "IntStdEarth"+std::to_string(0)+PremName+std::to_string(rhopct)+"_"+std::to_string(cthbins)+std::to_string(abins)+".txt"
        std::string efile = ResultFolder + IntEvntsFolder  + ename;
        ExportToCSV(TrueStd_e[nhist],efile);
    }
    

    //Sentivitity to Earth layers.
    std::string ResultFolder = "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/SimulationResults/";
    std::string SenvFolder = "/Sensitivity/Sensitivity2Layers/";
    std::string chi2name = "EarthSenv"+PremName+std::to_string(rhopct)+"_"+std::to_string(cthbins)+std::to_string(abins)+".txt"
    std::string chi2path = ResultFolder + SenvFolder  + chi2name;
    
    std::ofstream chi2file(chi2path);

    //TGraph * MultiChi     = new TGraph(MaxLayers);
    //TGraph * ChiPoints_e  = new TGraph(MaxLayers);
    //TGraph * ChiPoints_mu = new TGraph(MaxLayers);
    //TGraph * ChiPoints = new TGraph(MaxLayers);

    //Alternative Earth
    for (int layer = 1; layer <= MaxLayers; ++layer)
    {
    //PremName = "Prem44pct10/prem_44layers_"+std::to_string(layer) +"_10pct.txt";
    //path2prem = PremFolder+PremName; // std prem
    //TrueEvents.label = "Alt"+std::to_string(layer);
    
    TrueEvents.PremTable = path2prem;
    TrueEvents.MantleAnomaly = false;
    TrueEvents.ModifyLayer(layer,rhopct,0.0);

    //Muon-neutrino signal
    TrueEvents.flvf=1;
    std::vector <TH2D *> TrueAlt_mu = TrueEvents.GetTrueEvents2D();

    //Electron Neutrino signal
    TrueEvents.flvf=0;
    std::vector <TH2D *> TrueAlt_e = TrueEvents.GetTrueEvents2D();

    //Calculate chi2 values
    double chi2e = Get2DChi2( TrueStd_e[4],  TrueAlt_e[4]);   // Electron neutrino
    double chi2mu = Get2DChi2( TrueStd_mu[4], TrueAlt_mu[4]); // Muon Neutrino
    double chi2total = chi2e + chi2mu;

    chi2file << std::setprecision(10) << chi2e << " " << chi2mu  << " " << chi2total << std::endl;

    //ChiPoints_e->SetPoint(layer-1, chi2e, layer);
    //ChiPoints_mu->SetPoint(layer-1, chi2mu, layer);
    //ChiPoints->SetPoint(layer-1, chi2total, layer);

    
    }

    //MultiChi->Add(ChiPoints);
    //MultiChi->Add(ChiPoints_e);
    //MultiChi->Add(ChiPoints_mu);

    chi2file.close();

    /*
    //Some Neutrino Oscillations---------------------------------------------------------------------------------------------------------------------------------------

    // Neutrino final flavour
    int nue        = 0;  // electron neutrino  
    int numu       = 1; // muon neutrino
  
    bool  nunubar      = false; //antineutrino
    //double cth  = -0.827;//correspond to layer 25
    double cth  = -0.885;//correspond to layer 18 in p44, 47 in p64

    TGraph * OscPlot_mumu =  TrueEvents.GetOscProb( numu, numu, nunubar, cth); //To be Deleted
    OscPlot_mumu->SetLineColor(kBlue);
    OscPlot_mumu->SetLineWidth(5);
    Pnu->Add(OscPlot_mumu);

    TGraph * OscPlot_emu =  TrueEvents.GetOscProb( nue, numu, nunubar, cth ); //To be Deleted
    OscPlot_emu->SetLineColor(kRed);
    OscPlot_emu->SetLineWidth(5);
    Pnu->Add(OscPlot_emu);
    //TH2D * TrueStd = StandardEarth.TestTrueEvents2D("/home/dehy0499/OscProb/PremTables/prem_default.txt","/home/dehy0499/OscProb/PremTables/prem_default.txt");

   


   //Oscillation Probability
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

    c2->Divide(1,2,0.01,0.01);
    c2->cd(1);
    gStyle->SetPalette(kBlueGreenYellow);
    TrueStd_mu->Draw("COLZ");
    TrueStd_mu->SetStats(0);
    gPad->SetLogy();
    gPad->Update();

    c2->cd(2);
    gStyle->SetPalette(kBlueGreenYellow);
    TrueStd_e->Draw("COLZ");
    TrueStd_e->SetStats(0);
    gPad->SetLogy();
    gPad->Update();


    
    c2->Modified(); 
    c2->Update();
    TRootCanvas *rc2 = (TRootCanvas *)c2->GetCanvasImp();



    //Chi2 points
    TCanvas *c3 = new TCanvas();

    c3->Divide(2,2,0.01,0.01);
    c3->cd(1);

    ChiPoints_e->Draw();
    gPad->SetLogx();
    gPad->Update();

    c3->cd(2);
    ChiPoints_mu->Draw();
    gPad->SetLogx();
    gPad->Update();
    
    c3->cd(3);
    ChiPoints->Draw();
    gPad->SetLogx();
    gPad->Update();


    
    c3->Modified(); 
    c3->Update();
    TRootCanvas *rc3 = (TRootCanvas *)c3->GetCanvasImp();


    app.Run();

    std::cout << PremName << std::endl;
    */


  return 0;

}

void ExportToCSV(TH2D* hist, std::string filename)
{
   
    std::ofstream outfile(filename);

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

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
# define Rearth 6371.0 //km
# define Ratm 6368.0 //km
# define mN   1.67492749804E-27  // Nucleons mass ~ Neutron mass
# define MTon  1E9  //Metric MegaTon
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

    // Number of bins
    
    int czbins=50 ; // # Bins in zenith/cos(zenith)
    int abins=50 ; // # Bins in azimuth
    int ebins=50 ; // bins in energy
    

    //Interval of integration in Zenith, Azimuth and E for Neutrino Events

    //Energy interval (in GeV)--------------------------------------------
    double EnuMin=1.0 ; 
    double EnuMax=10.0 ;

    //Zenith Angle Interval-----------------------------------------------

    double zenmin = 91.0; // Not Able to use 90 degrees
    double zenmax = 180.0;
    double czmin = cos(zenmax*TMath::Pi()/180.0);
    double czmax = cos(zenmin*TMath::Pi()/180.0);


    //Azimuthal Interal-------------------------------------------------------------------------
    double phimin = 0.0;
    double phimax = 360.0 ;
   

    //Detector size-----------------------------------------------------------------------------

    double DetMass = 10.0*MTon; //Mass in megaton units
    double Nn      = DetMass/mN; //Number of target nucleons in the detector (Detector Mass / Nucleons Mass)
    double T       = 10.0*years2sec; //Detector Exposure time in sec: One Year

    double NnT = Nn*T; // Exposure Mton*years

    std::vector<double> chi2data;

    std::string chi2directory = "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/SimulationResults/chi2results/chi2true/";

    ofstream EarthChi2(chi2directory+"chi2Earth.csv", std::ofstream::trunc); //Opens a file and rewrite content, if files does not exist it Creates new file

    double chi2 = 0;

    //Event generation-----------------------------------------------------------------------------------------
    int nuflv = 1; // neutrino  final state options nue (0), numu (1) or nutau (2)
    
    std::string PremName = "prem_44layers";
    
    // Standard Earth model

    AsimovSimulation StandardEarth;
    StandardEarth.PremModel = PremName;
    //StandardEarth.MantleAnomaly = false;
    StandardEarth.SetIntervals(zenmin,zenmax,phimin,phimax,EnuMin,EnuMax);
    StandardEarth.SetBinning(czbins,abins,ebins);
    StandardEarth.SetExposure(NnT);
    StandardEarth.flvf=nuflv;

    TH2D * TrueStd = StandardEarth.GetTrueEvents2D();
    //TH2D * TrueStd = StandardEarth.TestTrueEvents2D("/home/dehy0499/OscProb/PremTables/prem_default.txt","/home/dehy0499/OscProb/PremTables/prem_default.txt");

   

    // Alternative Earth Model


    int TotalLayers = 44;
   
    AsimovSimulation AlternativeEarth;

        //std::string PremAltName = PremName + "_" + std::to_string(i);
       //std::cout << PremAltName << std::endl;

    for (int i = 1; i <= TotalLayers; i++)
    {

       AlternativeEarth.PremModel = PremName;
       //AlternativeEarth.MantleAnomaly = false;
       //AlternativeEarth.AnomalyShape="pancake";
       AlternativeEarth.ModifyLayer(i,10.0,0.0);
       AlternativeEarth.SetIntervals(zenmin,zenmax,phimin,phimax,EnuMin,EnuMax);
       AlternativeEarth.SetBinning(czbins,abins,ebins);
       AlternativeEarth.SetExposure(NnT);
       AlternativeEarth.flvf=nuflv;

       TH2D * TrueAlt = AlternativeEarth.GetTrueEvents2D();
       //TH2D * TrueAlt = AlternativeEarth.TestTrueEvents2D("/home/dehy0499/OscProb/PremTables/prem_default.txt","/home/dehy0499/OscProb/PremTables/prem_llsvp.txt");

       chi2 = Get2DChi2( TrueStd, TrueAlt);

       chi2data.push_back(chi2);

       EarthChi2 << i << " , " << chi2 <<  std::endl;

    }
   EarthChi2.close();

  // EVENT VISUALIZATION


   //TH2D * TrueDiff2D = new TH2D("TrueHist","True Event Histrogram", czbins,czmin,czmax,ebins,EnuMin,EnuMax); //binning in cth 

   //GetDiff2D( TrueStd , TrueAlt, TrueDiff2D );

  /*
    TApplication app("app", &argc, argv);

    TCanvas *c = new TCanvas();
    gStyle->SetPalette(kBird);
    TrueDiff2D->Draw("COLZ");
    TrueDiff2D->SetStats(0);
    gPad->Update();

    c->Modified(); c->Update();

    TRootCanvas *rc = (TRootCanvas *)c->GetCanvasImp();
    rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
    app.Run();
*/




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